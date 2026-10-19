"""Device-resident mesh and transport-kernel data (ROADMAP C02 item 2; CUDA env).

``DeviceDataset.load`` puts one audited LiouvilleCOPS run on one GPU in a
structure-of-arrays (SoA) layout. The C06 pressure and Jacobian kernels
consume this layout. The loader works in four steps, and nothing is
allocated on the device until the first two have passed:

1. **Validate identities.** Read the table and ``run_manifest.json``. Check
   ``elements.bin`` completely and compare its SHA-256 with the manifest
   (otherwise it is the wrong mesh). Compute the mesh identity with the
   reference definition (``activity_templates.mesh_identity``) and compare
   it with any identity the caller expects, for example from a template
   bundle or a gravity cache.
2. **Scan.** Read only the 92-byte header of every requested kernel file.
   This checks version, iPoint, record count and exact file size, and
   compares the header geometry with the table row. It yields exact record
   counts and the int64 CSR offsets.
3. **Plan.** Declare every device array with its exact size in a
   ``MemoryPlan`` and admit it through the ``DeviceBudget`` (20% headroom).
   An insufficient device is refused here, before any allocation.
4. **Stream.** Read the records file by file into a reusable host buffer.
   Validate them (face range, duplicates, finite values, shadow flag) and
   hash them. De-interleave them into two alternating pinned staging sets
   and copy them asynchronously to their final device offsets. Disk reading
   overlaps the transfer of the previous chunk. A failure at any point frees
   everything allocated so far, so a partially loaded store never exists.

After the upload, a device-side audit kernel (``kernels/store_audit.cu``)
re-checks the record invariants in device memory. ``verify_roundtrip``
downloads every field and compares its SHA-256 with the digest computed on
the host during the upload. That proves that load/release cycles preserve the
data bit for bit.

Device layout
-------------
Per contribution r, in kernel-file order within a point and in caller row
order across points (this is the reference summation order):

==========  ========  ====================================================
``face``    int32     facet index, checked to lie in [0, nfaces)
``kn``      float64   Kn
``kphi``    float64   KPhi
``kv``      float64   KV
``cos``     float64   cos(SZA) of the facet at this point's epoch
``shadow``  uint8     cast-shadow flag, checked to be 0 or 1
==========  ========  ====================================================

Per point p:

==============  ==============  ===============================================
``offsets``     int64 [npt+1]   CSR row pointer; records of p are
                                ``[offsets[p], offsets[p+1])``
``geometry``    float64 [10,npt] kernel header (``formats.KERNEL_GEOMETRY_NAMES``),
                                one contiguous row per quantity
==============  ==============  ===============================================

Per facet f (only the fields requested, see ``MESH_FIELD_SETS``):
``centroid`` and ``normal`` [3, nfaces]; ``Y`` [25, nfaces]; ``area``,
``AH2O`` and ``ACO2`` [nfaces]. Mesh vectors are stored component-major. In
the C06/C03 kernels each thread handles one facet, so thread f reads
``Y[j, f]`` for fixed j and the threads of a warp read consecutive addresses.

A bound per-facet activity basis (harmonics and templates) is
``enhancement_basis`` [ncols, nfaces], also component-major.

Records occupy 37 bytes on the device instead of 40 on disk. The campaign
cache of 245,515,533 records is 9.08 GB, which fits on one A6000 together
with the mesh.
"""

import hashlib
import time
import traceback
from pathlib import Path

import numpy as np

from . import cuda_sources, formats, precision, provenance, reference
from .errors import BackendError, CacheFormatError, IdentityError
from .memory import DeviceBudget, MemoryPlan, human


DATASET_SCHEMA = "liouville_cops_cuda.device_dataset"
DATASET_SCHEMA_VERSION = 1

# (device field, file field) in the order the arrays are declared/uploaded.
RECORD_FIELDS = (("face", "face"), ("kn", "Kn"), ("kphi", "KPhi"),
                 ("kv", "KV"), ("cos", "cos"), ("shadow", "shadow"))
AUDIT_COUNTERS = ("face_out_of_range", "nonfinite_value", "invalid_shadow",
                  "decreasing_offsets", "bad_offset_endpoints")
# Mesh fields needed by each consumer.  The evaluator (fit_cops.Model) reads
# Y and ACO2 per facet; gravity/slope work reads centroid, normal and area.
MESH_FIELD_SETS = {
    "evaluator": ("Y", "ACO2"),
    "gravity": ("centroid", "normal", "area"),
    "all": ("centroid", "normal", "area", "Y", "AH2O", "ACO2"),
    "none": (),
}
COMPONENT_MAJOR = frozenset(("centroid", "normal", "Y"))
DEFAULT_CHUNK_RECORDS = 1 << 22        # 4.2M records: 168 MB read buffer
AUDIT_THREADS = 256


def _device_shape(name, nfaces):
    """Device shape of one mesh field (component-major for vector fields)."""
    widths = {"centroid": 3, "normal": 3, "Y": formats.ELEMENT_NBASIS}
    return (widths[name], nfaces) if name in widths else (nfaces,)


# ------------------------------------------------------------------ run inputs
class RunInputs:
    """Validated host-side description of one base run (no device memory).

    Holds the parsed table, the run manifest, the validated element fields
    and their identities. Creating it reads ``elements.bin`` completely
    (214 MB for the 788k mesh), but no kernel records.
    """

    def __init__(self, run, expected_mesh_identity=None):
        cops_io = reference.load("cops_io")
        templates = reference.load("activity_templates")
        self.run = Path(run).expanduser().resolve()
        table_path = self.run / "COPS_model.dat"
        try:
            (self.headers, self.names, self.rows,
             self.columns) = cops_io.read_model_table(
                table_path, required_columns=("iPoint", "UTC", "et", "x", "y", "z",
                                              "altitude", "rh_AU", "VdotB_RG"))
        except cops_io.ProductError as exc:
            raise CacheFormatError(str(exc)) from None
        self.point_ids = [int(row[self.columns["iPoint"]]) for row in self.rows]
        self.table_identity = provenance.file_identity(table_path, role="model_table")

        manifest_path = self.run / "run_manifest.json"
        self.manifest = None
        if manifest_path.is_file():
            try:
                self.manifest = cops_io.load_manifest(manifest_path, cops_io.RUN_SCHEMA)
            except cops_io.ProductError as exc:
                raise CacheFormatError(str(exc)) from None
            if self.manifest.get("states", {}).get("process") != "complete":
                raise IdentityError("base run manifest does not mark the process complete: %s"
                                    % manifest_path)
        recorded = (self.manifest or {}).get("outputs", {})
        if recorded.get("model_table") and (
                recorded["model_table"].get("sha256") != self.table_identity["sha256"]):
            raise IdentityError("COPS_model.dat does not match run_manifest.json")

        self.elements, self.elements_identity = formats.load_elements(self.run / "elements.bin")
        if recorded.get("elements") and (
                recorded["elements"].get("sha256") != self.elements_identity["sha256"] or
                recorded["elements"].get("size_bytes") != self.elements_identity["size_bytes"]):
            raise IdentityError("elements.bin (sha256 %s) does not match run_manifest.json "
                                "(sha256 %s): a different or reordered mesh" %
                                (self.elements_identity["sha256"],
                                 recorded["elements"].get("sha256")))
        self.nfaces = self.elements_identity["nfaces"]
        if self.nfaces >= 2 ** 31:
            raise CacheFormatError("mesh has %d faces; int32 face indices cannot address it"
                                   % self.nfaces)
        # The reference's definition of "the same facet order": exact
        # centroid/area bytes.  Template bundles carry this value.
        self.mesh_identity = templates.mesh_identity(self.elements["centroid"],
                                                     self.elements["area"])
        if expected_mesh_identity is not None and expected_mesh_identity != self.mesh_identity:
            raise IdentityError("mesh identity %s of %s differs from the expected %s" %
                                (self.mesh_identity, self.run, expected_mesh_identity))

    def resolve_rows(self, rows=None, point_ids=None):
        """Table row indices to load, in the caller's order.

        Pass either ``rows`` (table row indices, the convention of
        ``fit_cops.Model.rows``) or ``point_ids``, not both. Without either,
        every table row is loaded in table order. Duplicates are refused,
        because a point loaded twice would be counted twice in any sum.
        """
        if rows is not None and point_ids is not None:
            raise ValueError("pass rows or point_ids, not both")
        if point_ids is not None:
            lookup = {point: index for index, point in enumerate(self.point_ids)}
            missing = [p for p in point_ids if int(p) not in lookup]
            if missing:
                raise IdentityError("point ids not in the table: %s" % missing[:5])
            rows = [lookup[int(p)] for p in point_ids]
        if rows is None:
            rows = range(len(self.rows))
        rows = [int(r) for r in rows]
        if any(r < 0 or r >= len(self.rows) for r in rows):
            raise IndexError("row index outside the table of %d rows" % len(self.rows))
        if len(set(rows)) != len(rows):
            raise ValueError("duplicate rows requested")
        return rows

    def is_full_table(self, rows):
        return list(rows) == list(range(len(self.rows)))


# ------------------------------------------------------------------------ scan
class KernelScan:
    """Header-only pass over the requested kernel files (no records read)."""

    def __init__(self, inputs, rows):
        self.rows = list(rows)
        self.point_ids = [inputs.point_ids[r] for r in self.rows]
        self.paths, self.headers = [], []
        counts = np.empty(len(self.rows), dtype=np.int64)
        self.geometry = np.empty((len(self.rows), 10), dtype=np.float64)
        kernel_dir = inputs.run / "kernels"
        names = [name for name, _, _, _ in formats.TABLE_HEADER_COMPARISONS]
        for index, (row, point) in enumerate(zip(self.rows, self.point_ids)):
            path = kernel_dir / formats.kernel_file_name(point)
            header, count, geometry = formats.read_kernel_header(path, point, inputs.nfaces)
            values = {name: inputs.rows[row][inputs.columns[name]] for name in names}
            formats.check_table_geometry(point, geometry, values, path)
            self.paths.append(path)
            self.headers.append(header)
            counts[index] = count
            self.geometry[index] = geometry
        self.counts = counts
        # int64 offsets: the running sum exceeds 2^31 for caches about 9x the
        # campaign's.  np.cumsum on int64 input stays int64.
        self.offsets = np.zeros(len(self.rows) + 1, dtype=np.int64)
        np.cumsum(counts, out=self.offsets[1:])
        self.nrecords = int(self.offsets[-1])


# ------------------------------------------------------------------------ plan
def dataset_plan(nfaces, npt, nrecords, mesh_fields, basis_columns=0):
    """Exact device memory of a dataset (static data only)."""
    plan = MemoryPlan("device dataset")
    for name, _ in RECORD_FIELDS:
        plan.add_array("kernel." + name, "static", (nrecords,), precision.FIELD_DTYPES[name])
    plan.add_array("kernel.offsets", "static", (npt + 1,), precision.FIELD_DTYPES["offsets"])
    plan.add_array("kernel.geometry", "static", (10, npt), precision.FIELD_DTYPES["geometry"])
    for name in mesh_fields:
        plan.add_array("mesh." + name, "static", _device_shape(name, nfaces),
                       precision.FIELD_DTYPES[name])
    if basis_columns:
        plan.add_array("basis.enhancement", "basis", (basis_columns, nfaces),
                       precision.FIELD_DTYPES["enhancement_basis"])
    return plan


# --------------------------------------------------------------------- dataset
class DeviceDataset:
    """One run's mesh and transport kernels resident on one device.

    Use as a context manager, or call ``release()``. After release every
    device array is gone, the reservation is returned and the CuPy pool
    blocks are freed.
    """

    def __init__(self):
        self.arrays = {}
        self.host_digests = {}
        self.reservation = None
        self.budget = None
        self.owns_budget = False
        self.released = False
        self.timing = {}
        self.audit_result = None
        self.basis_identity = None

    # ...................................................................load
    @classmethod
    def load(cls, run, rows=None, point_ids=None, device_id=0, mesh_fields="evaluator",
             budget=None, enhancement_basis=None, expected_mesh_identity=None,
             expected_kernel_digests=None, chunk_records=DEFAULT_CHUNK_RECORDS,
             audit=True):
        """Validate, plan, allocate and stream; see the module docstring.

        Parameters
        ----------
        run : base run directory (``COPS_model.dat``, ``elements.bin``,
            ``kernels/``, ``run_manifest.json``).
        rows / point_ids : the observations to load, in evaluation order.
        mesh_fields : a key of ``MESH_FIELD_SETS`` or an explicit tuple.
        budget : a shared ``DeviceBudget``. Without one, a private budget
            with the default 20% headroom is created and closed on release.
        enhancement_basis : optional dict with keys ``matrix`` (nfaces x k),
            ``mesh_identity``, ``scale`` and ``description``; see
            ``bind_basis``. It is part of the plan from the start.
        expected_mesh_identity : refuse a mesh whose identity differs.
        expected_kernel_digests : optional ``{file_name: sha256}``. Each
            streamed file is compared with it, for example the C00
            ``identities/kernel_files.json``.
        """
        self = cls()
        started = time.perf_counter()
        inputs = RunInputs(run, expected_mesh_identity=expected_mesh_identity)
        rows = inputs.resolve_rows(rows, point_ids)
        if isinstance(mesh_fields, str):
            mesh_fields = MESH_FIELD_SETS[mesh_fields]
        mesh_fields = tuple(mesh_fields)
        scan = KernelScan(inputs, rows)
        self.timing["validate_and_scan_s"] = time.perf_counter() - started
        basis_matrix = None
        if enhancement_basis is not None:
            basis_matrix = self._check_basis(enhancement_basis, inputs)

        self.inputs, self.scan, self.device_id = inputs, scan, int(device_id)
        self.mesh_fields = mesh_fields
        self.plan = dataset_plan(inputs.nfaces, len(rows), scan.nrecords, mesh_fields,
                                 0 if basis_matrix is None else basis_matrix.shape[1])
        import cupy
        self.cupy = cupy
        if budget is None:
            budget, self.owns_budget = DeviceBudget(self.device_id), True
        elif budget.device_id != self.device_id:
            raise ValueError("budget belongs to device %d, dataset to %d" %
                             (budget.device_id, self.device_id))
        self.budget = budget
        try:
            # Admission happens here: CapacityError leaves nothing allocated.
            self.reservation = budget.reserve(self.plan)
            with cupy.cuda.Device(self.device_id), budget.allocation_guard(self.plan):
                t0 = time.perf_counter()
                self._upload_mesh(mesh_fields)
                if basis_matrix is not None:
                    self._upload_basis(basis_matrix, enhancement_basis)
                self.timing["mesh_upload_s"] = time.perf_counter() - t0
                t0 = time.perf_counter()
                self._stream_kernels(chunk_records, expected_kernel_digests)
                self.timing["kernel_stream_s"] = time.perf_counter() - t0
                if audit:
                    t0 = time.perf_counter()
                    self.audit()
                    self.timing["device_audit_s"] = time.perf_counter() - t0
        except BaseException as exc:
            # The traceback keeps the failed frames alive, and with them any
            # device array held in a local variable (e.g. the per-field views
            # in _stream_kernels).  Clear those frames first so the arrays are
            # really freed before the reservation returns the pool blocks.
            self.arrays.clear()
            traceback.clear_frames(exc.__traceback__)
            self.release()
            raise
        self.timing["load_total_s"] = time.perf_counter() - started
        return self

    # .................................................................. mesh
    def _alloc(self, name, shape, dtype):
        """Allocate one planned device array, enforcing the precision policy."""
        precision.require_device_dtype(name.split(".", 1)[1], dtype)
        array = self.cupy.empty(shape, dtype=dtype)
        self.arrays[name] = array
        return array

    def _upload_mesh(self, mesh_fields):
        for name in mesh_fields:
            host = self.inputs.elements[name]
            # Component-major layout: an exact transpose, values unchanged.
            host = np.ascontiguousarray(host.T if name in COMPONENT_MAJOR else host)
            device = self._alloc("mesh." + name, host.shape, host.dtype)
            device.set(host)
            self.host_digests["mesh." + name] = hashlib.sha256(host.tobytes()).hexdigest()

    @staticmethod
    def _check_basis(basis, inputs):
        """Validate a per-facet basis against this mesh before planning it."""
        matrix = np.asarray(basis["matrix"])
        if matrix.dtype != np.float64:
            raise TypeError("enhancement basis must be float64 (precision policy); got %s"
                            % matrix.dtype)
        if matrix.ndim != 2 or matrix.shape[0] != inputs.nfaces or matrix.shape[1] == 0:
            raise IdentityError("enhancement basis has shape %s; this mesh needs (%d, k>0)" %
                                (matrix.shape, inputs.nfaces))
        if basis.get("mesh_identity") != inputs.mesh_identity:
            raise IdentityError("enhancement basis was built for mesh %s, not %s" %
                                (basis.get("mesh_identity"), inputs.mesh_identity))
        if not np.isfinite(matrix).all():
            raise CacheFormatError("enhancement basis contains non-finite values")
        return matrix

    def _upload_basis(self, matrix, basis):
        host = np.ascontiguousarray(matrix.T)
        device = self._alloc("basis.enhancement_basis", host.shape, host.dtype)
        device.set(host)
        digest = provenance.array_sha256(matrix)
        self.host_digests["basis.enhancement_basis"] = hashlib.sha256(host.tobytes()).hexdigest()
        self.basis_identity = {"columns": int(matrix.shape[1]), "array_sha256": digest,
                               "mesh_identity": basis["mesh_identity"],
                               "scale": float(basis.get("scale", 1.0)),
                               "description": basis.get("description")}

    # ............................................................... kernels
    def _stream_kernels(self, chunk_records, expected_digests):
        """Stream every scanned file to the device (see module docstring)."""
        cupy, scan, nfaces = self.cupy, self.scan, self.inputs.nfaces
        import cupyx
        nrec = scan.nrecords
        device = {}
        for name, _ in RECORD_FIELDS:
            device[name] = self._alloc("kernel." + name, (nrec,), precision.FIELD_DTYPES[name])
        offsets = self._alloc("kernel.offsets", scan.offsets.shape, np.int64)
        offsets.set(scan.offsets)
        geometry_host = np.ascontiguousarray(scan.geometry.T)
        geometry = self._alloc("kernel.geometry", geometry_host.shape, np.float64)
        geometry.set(geometry_host)
        self.host_digests["kernel.offsets"] = hashlib.sha256(scan.offsets.tobytes()).hexdigest()
        self.host_digests["kernel.geometry"] = hashlib.sha256(geometry_host.tobytes()).hexdigest()

        chunk_records = int(max(1, min(chunk_records, max(nrec, 1))))
        raw = np.empty(chunk_records, dtype=formats.KERNEL_RECORD_DTYPE)
        raw_bytes = raw.view(np.uint8)
        record_size = formats.KERNEL_RECORD_DTYPE.itemsize
        # Two pinned staging sets: while the GPU copies set k, the host reads
        # and de-interleaves into set 1-k.  An event per set guards reuse.
        staging = [{name: cupyx.empty_pinned((chunk_records,), dtype=precision.FIELD_DTYPES[name])
                    for name, _ in RECORD_FIELDS} for _ in range(2)]
        events = [None, None]
        stream = cupy.cuda.Stream(non_blocking=True)
        digests = {name: hashlib.sha256() for name, _ in RECORD_FIELDS}
        face_seen = np.zeros(nfaces, dtype=bool)
        state = {"fill": 0, "start": 0, "set": 0}
        self.host_staging_bytes = 2 * chunk_records * sum(
            precision.FIELD_DTYPES[name].itemsize for name, _ in RECORD_FIELDS)

        def flush():
            fill = state["fill"]
            if fill == 0:
                return
            k = state["set"]
            if events[k] is not None:
                events[k].synchronize()
            begin = state["start"]
            for name, file_field in RECORD_FIELDS:
                target = staging[k][name][:fill]
                # shadow: int32 -> uint8 is lossless because every value was
                # checked to be 0 or 1; the other fields keep their dtype.
                np.copyto(target, raw[file_field][:fill], casting="unsafe")
                digests[name].update(target)
                device[name][begin:begin + fill].set(target, stream=stream)
            events[k] = stream.record()
            state["start"] = begin + fill
            state["fill"] = 0
            state["set"] = 1 - k

        entries, verified = [], 0
        try:
            for index, path in enumerate(scan.paths):
                count = int(scan.counts[index])
                file_digest = hashlib.sha256()
                with path.open("rb") as handle:
                    header = handle.read(formats.KERNEL_HEADER_DTYPE.itemsize)
                    # The header must be byte-identical to the scanned one:
                    # the plan and offsets were computed from it.
                    if header != scan.headers[index]:
                        raise IdentityError("%s changed between scan and upload" % path)
                    file_digest.update(header)
                    face_seen[:] = False
                    first = 0
                    while first < count:
                        fill = state["fill"]
                        n = min(count - first, chunk_records - fill)
                        view = raw_bytes[fill * record_size:(fill + n) * record_size]
                        got = handle.readinto(view)
                        if got != n * record_size:
                            raise CacheFormatError("%s is shorter than its header declares" % path)
                        file_digest.update(view)
                        formats.validate_records(raw[fill:fill + n], nfaces, path, first, face_seen)
                        state["fill"] = fill + n
                        first += n
                        if state["fill"] == chunk_records:
                            flush()
                    if handle.read(1):
                        raise CacheFormatError("%s has bytes beyond its declared records" % path)
                name = path.name
                digest = file_digest.hexdigest()
                if expected_digests is not None:
                    if name not in expected_digests:
                        raise IdentityError("no expected digest for %s" % name)
                    if expected_digests[name] != digest:
                        raise IdentityError("%s sha256 %s differs from the expected %s" %
                                            (path, digest, expected_digests[name]))
                    verified += 1
                entries.append((name, digest))
            flush()
            stream.synchronize()
        finally:
            # Pinned buffers are host memory owned by CuPy's pinned pool;
            # drop them only after the copies that read them have completed.
            stream.synchronize()
            del staging
            cupy.get_default_pinned_memory_pool().free_all_blocks()
        if state["start"] != nrec:
            raise BackendError("internal: streamed %d records, planned %d" % (state["start"], nrec))
        for name, _ in RECORD_FIELDS:
            self.host_digests["kernel." + name] = digests[name].hexdigest()

        self.kernel_set_sha256 = formats.kernel_set_digest(entries)
        self.file_digests_verified = verified
        self.manifest_set_verified = False
        recorded = (self.inputs.manifest or {}).get("outputs", {}).get("kernels", {})
        if self.inputs.is_full_table(scan.rows) and recorded.get("set_sha256"):
            if recorded["set_sha256"] != self.kernel_set_sha256:
                raise IdentityError("kernel set sha256 %s differs from run_manifest.json %s" %
                                    (self.kernel_set_sha256, recorded["set_sha256"]))
            if recorded.get("records") not in (None, nrec):
                raise IdentityError("record count %d differs from run_manifest.json %s" %
                                    (nrec, recorded.get("records")))
            self.manifest_set_verified = True

    # ................................................................ checks
    def audit(self):
        """Re-check record and offset invariants in device memory.

        Returns the counter dictionary and raises ``CacheFormatError`` if any
        counter is non-zero.
        """
        cupy = self.cupy
        with cupy.cuda.Device(self.device_id):
            counters = cupy.zeros(len(AUDIT_COUNTERS), dtype=cupy.uint64)
            nrec, npt = self.scan.nrecords, len(self.scan.rows)
            a = self.arrays
            device = cupy.cuda.Device(self.device_id)
            # Enough blocks to fill the device; the grid-stride loops cover
            # any size, so the grid does not grow with the data.
            blocks = int(max(1, min(-(-max(nrec, 1) // AUDIT_THREADS),
                                    32 * device.attributes["MultiProcessorCount"])))
            if nrec:
                cuda_sources.function("store_audit", "audit_records", self.device_id)(
                    (blocks,), (AUDIT_THREADS,),
                    (a["kernel.face"], a["kernel.kn"], a["kernel.kphi"], a["kernel.kv"],
                     a["kernel.cos"], a["kernel.shadow"], np.int64(nrec),
                     np.int32(self.inputs.nfaces), counters))
            cuda_sources.function("store_audit", "audit_offsets", self.device_id)(
                (max(1, min(-(-max(npt, 1) // AUDIT_THREADS), 1024)),), (AUDIT_THREADS,),
                (a["kernel.offsets"], np.int64(npt), np.int64(nrec), counters))
            values = [int(v) for v in counters.get()]
            del counters
        result = dict(zip(AUDIT_COUNTERS, values))
        self.audit_result = result
        if any(values):
            raise CacheFormatError("device audit failed: %s" % result)
        return result

    def device_digests(self, chunk_bytes=256 << 20):
        """SHA-256 of every device array's bytes, downloaded in chunks."""
        cupy = self.cupy
        result = {}
        with cupy.cuda.Device(self.device_id):
            for name, array in self.arrays.items():
                flat = array.reshape(-1)
                step = max(1, chunk_bytes // max(flat.itemsize, 1))
                digest = hashlib.sha256()
                for start in range(0, flat.size, step):
                    digest.update(cupy.asnumpy(flat[start:start + step]))
                result[name] = digest.hexdigest()
        return result

    def verify_roundtrip(self):
        """Prove the device bytes equal the validated host bytes, field by field."""
        t0 = time.perf_counter()
        device = self.device_digests()
        mismatched = sorted(name for name in self.host_digests
                            if device.get(name) != self.host_digests[name])
        if mismatched or set(device) != set(self.host_digests):
            raise BackendError("device data differ from the uploaded host data: %s" % mismatched)
        self.timing["roundtrip_verify_s"] = time.perf_counter() - t0
        return device

    # ............................................................... release
    def release(self):
        """Free every device array and return the reservation (idempotent)."""
        if self.released:
            return
        self.released = True
        self.arrays.clear()
        if self.reservation is not None:
            self.reservation.release()
            self.reservation = None
        if self.owns_budget and self.budget is not None:
            self.budget.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.release()
        return False

    # .............................................................. identity
    @property
    def npoints(self):
        return len(self.scan.rows)

    @property
    def nrecords(self):
        return self.scan.nrecords

    def layout(self):
        """Dtype, shape and size of every resident array."""
        return {name: {"dtype": array.dtype.str, "shape": list(array.shape),
                       "bytes": int(array.nbytes),
                       "host_sha256": self.host_digests.get(name)}
                for name, array in sorted(self.arrays.items())}

    def identity(self):
        """Manifest section describing exactly what is resident on the device."""
        inputs, scan = self.inputs, self.scan
        return {
            "schema": DATASET_SCHEMA, "schema_version": DATASET_SCHEMA_VERSION,
            "run": str(inputs.run), "legacy_manifest": inputs.manifest is None,
            "model_table": inputs.table_identity, "elements": inputs.elements_identity,
            "mesh_identity": inputs.mesh_identity,
            "selection": {"rows": len(scan.rows),
                          "full_table_in_order": inputs.is_full_table(scan.rows),
                          "row_indices_sha256": provenance.canonical_sha256(scan.rows),
                          "point_ids_sha256": provenance.canonical_sha256(scan.point_ids),
                          "first_point": scan.point_ids[0] if scan.point_ids else None,
                          "last_point": scan.point_ids[-1] if scan.point_ids else None},
            "kernels": {"version": formats.KERNEL_VERSION, "points": len(scan.rows),
                        "records": scan.nrecords,
                        "loaded_set_sha256": self.kernel_set_sha256,
                        "verified_against_run_manifest": self.manifest_set_verified,
                        "files_verified_against_expected": self.file_digests_verified},
            "mesh_fields": list(self.mesh_fields),
            "enhancement_basis": self.basis_identity,
            "layout": self.layout(),
            "device_audit": self.audit_result,
            "memory_plan": self.plan.as_dict(),
            "host_pinned_staging_bytes": getattr(self, "host_staging_bytes", 0),
            "timing_s": dict(self.timing),
        }

    def summary(self):
        return ("%d points, %d records, %d faces, %s planned on device %d" %
                (self.npoints, self.nrecords, self.inputs.nfaces,
                 human(self.plan.total_bytes), self.device_id))
