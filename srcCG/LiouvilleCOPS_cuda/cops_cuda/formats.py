"""Binary formats of the LiouvilleCOPS caches and their vectorized validation.

Python 3.8 compatible; uses NumPy only.

The CPU reference reads these files in two places: ``cops_io`` validates them
record by record with ``struct``, and ``fit_cops`` loads them with NumPy. This
module restates the two layouts as NumPy dtypes. The test suite checks that
their sizes equal the reference definitions, ``cops_io.KERNEL1_HEADER``,
``cops_io.KERNEL1_RECORD`` and ``fit_cops.REC``, so the two cannot drift
apart unnoticed.

The checks here are those of ``cops_io.validate_elements`` and
``cops_io.validate_kernel``, vectorized so that the 245M-record campaign cache
can be validated while it is streamed to the device. The per-record Python
loop of ``cops_io`` takes hours on that cache.

The layouts, as written by ``LiouvilleCOPS/src/main.cpp:WriteKernels``:

``elements.bin`` (version 1)
    header ``<i4 version, <i4 nfaces, <i4 nbasis`` (nbasis must be 25), then
    ``nfaces`` records of 34 little-endian doubles:
    centroid[3], normal[3], area, Y[25], A_H2O, A_CO2.

``kernels/point=NNNN.bin`` (version 1, FlightTimeRotation=none)
    header ``<i4 version, <i4 iPoint, <f8 geometry[10], <i4 count``, where
    geometry = (et, V.b_RG, x_sc, y_sc, z_sc, sun_x, sun_y, sun_z,
    heliocentric distance [m], altitude [m]); then ``count`` records
    ``<i4 face, <f8 Kn, <f8 KPhi, <f8 KV, <f8 cos(SZA), <i4 shadow``.

The record order inside a kernel file is the order in which the reference
evaluator accumulates its sums (``np.bincount`` follows the input order). It
is therefore part of the reference and is preserved exactly by every loader.
"""

import hashlib
from pathlib import Path

import numpy as np

from .errors import CacheFormatError, IdentityError


AU_M = 149597870700.0

ELEMENT_VERSION = 1
ELEMENT_NBASIS = 25
ELEMENT_DOUBLES = 34
ELEMENT_HEADER_DTYPE = np.dtype([("version", "<i4"), ("nfaces", "<i4"),
                                 ("nbasis", "<i4")])
# Column ranges of one element record, in the order fit_cops.read_elements
# slices them.  A slice maps to a (nfaces, k) block, an int to one column.
ELEMENT_COLUMNS = (("centroid", slice(0, 3)), ("normal", slice(3, 6)),
                   ("area", 6), ("Y", slice(7, 32)), ("AH2O", 32),
                   ("ACO2", 33))

KERNEL_VERSION = 1
KERNEL_HEADER_DTYPE = np.dtype([("version", "<i4"), ("point", "<i4"),
                                ("geometry", "<f8", (10,)), ("count", "<i4")])
KERNEL_GEOMETRY_NAMES = ("et", "VdotB_RG", "sc_x", "sc_y", "sc_z", "sun_x",
                         "sun_y", "sun_z", "rh_m", "altitude_m")
KERNEL_RECORD_DTYPE = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"),
                                ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])

# NumPy structured dtypes are packed (no padding) unless align=True is
# requested, which is what the C++ fwrite sequence produces.  Checked here so
# that a NumPy behaviour change fails at import instead of mis-parsing data.
assert ELEMENT_HEADER_DTYPE.itemsize == 12
assert KERNEL_HEADER_DTYPE.itemsize == 92
assert KERNEL_RECORD_DTYPE.itemsize == 40

# Table-versus-kernel-header consistency checks of cops_io.validate_fit_base:
# (table column, geometry index, relative tolerance, divisor applied to the
# cached value).  The table prints fewer digits than the binary header, hence
# the tolerances; the comparison is a guard against pairing a table row with
# another point's kernel, not a precision check.
TABLE_HEADER_COMPARISONS = (("et", 0, 2e-6, 1.0), ("VdotB_RG", 1, 2e-6, 1.0),
                            ("x", 2, 2e-6, 1.0), ("y", 3, 2e-6, 1.0),
                            ("z", 4, 2e-6, 1.0), ("rh_AU", 8, 2e-8, AU_M),
                            ("altitude", 9, 2e-6, 1.0))


def kernel_file_name(point):
    """File name of one point's kernel, identical to the C++ writer."""
    return "point=%04d.bin" % int(point)


# --------------------------------------------------------------------- elements
def read_element_header(path):
    """Return ``(nfaces, nbasis, size_bytes)`` after checking version and size.

    The size check is exact, so it rejects truncation and appended bytes.
    """
    path = Path(path)
    if not path.is_file():
        raise CacheFormatError("element cache is missing: %s" % path)
    size = path.stat().st_size
    with path.open("rb") as stream:
        raw = stream.read(ELEMENT_HEADER_DTYPE.itemsize)
    if len(raw) != ELEMENT_HEADER_DTYPE.itemsize:
        raise CacheFormatError("truncated element header: %s" % path)
    header = np.frombuffer(raw, dtype=ELEMENT_HEADER_DTYPE)[0]
    version, nfaces, nbasis = (int(header["version"]), int(header["nfaces"]),
                               int(header["nbasis"]))
    if version != ELEMENT_VERSION or nfaces <= 0 or nbasis != ELEMENT_NBASIS:
        raise CacheFormatError(
            "unsupported element header in %s: version=%d faces=%d basis=%d" %
            (path, version, nfaces, nbasis))
    expected = ELEMENT_HEADER_DTYPE.itemsize + nfaces * ELEMENT_DOUBLES * 8
    if size != expected:
        raise CacheFormatError("element cache size mismatch in %s: %d bytes; expected %d" %
                               (path, size, expected))
    return nfaces, nbasis, size


ELEMENT_CHUNK_ROWS = 1 << 20      # faces per validation/hash chunk (272 MB of the file)


def _element_table(path, nfaces):
    """The element rows of ``path`` as a read-only memory map (page cache, not process memory)."""
    return np.memmap(path, dtype="<f8", mode="r", offset=ELEMENT_HEADER_DTYPE.itemsize,
                     shape=(nfaces, ELEMENT_DOUBLES))


def _file_key(path):
    st = Path(path).stat()
    return (st.st_size, st.st_mtime_ns, st.st_ino)


class LazyElementFields(dict):
    """``load_elements(..., lazy=True)`` fields: each one a C-contiguous float64 copy made on first
    access from the memory-mapped file, refused if the file changed since it was hashed."""

    def __init__(self, path, nfaces, key):
        super().__init__()
        self._path, self._nfaces, self._key = Path(path), nfaces, key
        self._columns = dict(ELEMENT_COLUMNS)

    def __missing__(self, name):
        if name not in self._columns:
            raise KeyError(name)
        if _file_key(self._path) != self._key:
            raise IdentityError("element cache %s changed after it was validated" % self._path)
        table = _element_table(self._path, self._nfaces)
        value = np.ascontiguousarray(table[:, self._columns[name]], dtype=np.float64)
        del table
        self[name] = value
        return value

    def __contains__(self, name):
        return name in self._columns

    def __iter__(self):
        return iter(self._columns)

    def __len__(self):
        return len(self._columns)

    def keys(self):
        return self._columns.keys()

    def items(self):
        return [(name, self[name]) for name in self._columns]

    def values(self):
        return [self[name] for name in self._columns]


def load_elements(path, lazy=False):
    """Read, validate and split ``elements.bin``; return fields and identity.

    Returns ``(fields, identity)``. ``fields`` maps the names of
    ``ELEMENT_COLUMNS`` to C-contiguous float64 arrays shaped like the
    arrays of ``fit_cops.read_elements``. ``identity`` holds the SHA-256 of
    the whole file, which is the definition used by ``cops_io`` and recorded
    in ``run_manifest.json``.

    The file is hashed and validated in chunks of rows through a memory map, so the
    process never holds the whole file (5.4 GB at SHAP8 20M) besides the fields.
    ``lazy=True`` makes each field only when it is first used (``LazyElementFields``):
    a device dataset needs centroid, area and its mesh fields, not Y (4 GB at 20M).
    """
    path = Path(path)
    nfaces, nbasis, size = read_element_header(path)
    key = _file_key(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        digest.update(stream.read(ELEMENT_HEADER_DTYPE.itemsize))
    table = _element_table(path, nfaces)
    for f0 in range(0, nfaces, ELEMENT_CHUNK_ROWS):
        chunk = np.asarray(table[f0:f0 + ELEMENT_CHUNK_ROWS])
        digest.update(chunk.tobytes())
        finite = np.isfinite(chunk).all(axis=1)
        if not finite.all():
            raise CacheFormatError("non-finite value in %s face %d" %
                                   (path, f0 + int(np.argmin(finite))))
        area = chunk[:, 6]
        if not (area > 0.0).all():
            raise CacheFormatError("non-positive area in %s face %d" %
                                   (path, f0 + int(np.argmin(area > 0.0))))
    if _file_key(path) != key:
        raise CacheFormatError("element cache %s changed while being read" % path)
    if lazy:
        fields = LazyElementFields(path, nfaces, key)
    else:
        fields = {name: np.ascontiguousarray(table[:, columns], dtype=np.float64)
                  for name, columns in ELEMENT_COLUMNS}
    del table
    identity = {"path": str(path.resolve()), "size_bytes": size,
                "sha256": digest.hexdigest(),
                "version": ELEMENT_VERSION, "nfaces": nfaces, "nbasis": nbasis}
    return fields, identity


def map_elements(path):
    """``fit_cops.read_elements``'s dict as read-only views of a memory map of ``path``.

    For multi-GPU workers, whose model needs Y only to copy its own facets' columns: the values
    are the file's bytes, as in the reference reader, but they stay in the (shared, reclaimable)
    page cache instead of 9.4 GB of process memory per worker at SHAP8 20M (the reference
    reader keeps the whole file plus a contiguous Y). ``Y`` is a strided view here; the
    caller checks it against the driver's Y digest (``provenance.array_sha256``).
    """
    nfaces = read_element_header(path)[0]
    table = _element_table(path, nfaces)
    return {name: table[:, columns] for name, columns in ELEMENT_COLUMNS}


# ---------------------------------------------------------------------- kernels
def read_kernel_header(path, expected_point, nfaces):
    """Read and validate one kernel header without reading its records.

    Returns ``(header_bytes, count, geometry)``. This is the cheap first pass
    of the loader: it establishes the exact record count of every requested
    point, so the device memory plan can be checked before any allocation.
    """
    path = Path(path)
    if not path.is_file():
        raise CacheFormatError("kernel file is missing: %s" % path)
    size = path.stat().st_size
    with path.open("rb") as stream:
        raw = stream.read(KERNEL_HEADER_DTYPE.itemsize)
    return validate_kernel_header(raw, size, path, expected_point, nfaces)


def validate_kernel_header(raw, size, path, expected_point, nfaces):
    """Validate header bytes ``raw`` of a file of ``size`` bytes (see above)."""
    if len(raw) != KERNEL_HEADER_DTYPE.itemsize:
        raise CacheFormatError("truncated kernel header: %s" % path)
    header = np.frombuffer(raw, dtype=KERNEL_HEADER_DTYPE)[0]
    version, point, count = (int(header["version"]), int(header["point"]),
                             int(header["count"]))
    if version != KERNEL_VERSION:
        raise CacheFormatError("%s has kernel version %d; static fitting needs version 1" %
                               (path, version))
    if point != int(expected_point):
        raise IdentityError("%s stores iPoint %d; expected %d" %
                            (path, point, int(expected_point)))
    # A point sees each facet at most once, so a count above nfaces can only
    # be corruption.  The bound also keeps one point's count within int32.
    if count < 0 or count > nfaces:
        raise CacheFormatError("%s has invalid record count %d for %d faces" %
                               (path, count, nfaces))
    expected_size = KERNEL_HEADER_DTYPE.itemsize + count * KERNEL_RECORD_DTYPE.itemsize
    if size != expected_size:
        raise CacheFormatError("kernel size mismatch in %s: %d bytes; expected %d" %
                               (path, size, expected_size))
    geometry = np.array(header["geometry"], dtype=np.float64)
    if not np.isfinite(geometry).all():
        raise CacheFormatError("non-finite value in %s header" % path)
    return bytes(raw), count, geometry


def check_table_geometry(point, geometry, table_values, path):
    """Compare a kernel header with its table row (cops_io tolerances).

    ``table_values`` maps the table column names of
    ``TABLE_HEADER_COMPARISONS`` to floats. A mismatch means that the table
    and the kernel describe different observations.
    """
    for name, index, relative, divisor in TABLE_HEADER_COMPARISONS:
        cached = float(geometry[index]) / divisor
        table_value = float(table_values[name])
        if abs(table_value - cached) > relative * max(1.0, abs(cached)):
            raise IdentityError("point %d %s differs between table and kernel %s (%r vs %r)" %
                                (point, name, path, table_value, cached))


def validate_records(records, nfaces, path, first_record, face_seen):
    """Validate one contiguous block of records belonging to one point.

    ``face_seen`` is a reusable boolean array of length ``nfaces``. It carries
    the faces already seen earlier in the same point, so a duplicate split
    across two streaming chunks is still detected. The caller clears it
    between points. The checks are those of ``cops_io.validate_kernel``.
    """
    if len(records) == 0:
        return
    face = records["face"]
    low, high = int(face.min()), int(face.max())
    if low < 0 or high >= nfaces:
        bad = int(np.argmax((face < 0) | (face >= nfaces)))
        raise CacheFormatError("%s record %d has out-of-range face %d" %
                               (path, first_record + bad, int(face[bad])))
    # Duplicates inside this block: a bincount over the mesh is cheaper than
    # sorting for blocks of ~1e5-1e6 records and keeps the check exact.
    if np.bincount(face, minlength=nfaces).max(initial=0) > 1 or face_seen[face].any():
        raise CacheFormatError("%s contains a duplicate face (records %d..%d)" %
                               (path, first_record, first_record + len(face) - 1))
    face_seen[face] = True
    for name in ("Kn", "KPhi", "KV", "cos"):
        values = records[name]
        if not np.isfinite(values).all():
            bad = int(np.argmin(np.isfinite(values)))
            raise CacheFormatError("non-finite %s in %s record %d" %
                                   (name, path, first_record + bad))
    shadow = records["shadow"]
    if not ((shadow == 0) | (shadow == 1)).all():
        bad = int(np.argmax((shadow != 0) & (shadow != 1)))
        raise CacheFormatError("%s record %d has invalid shadow flag %d" %
                               (path, first_record + bad, int(shadow[bad])))


def read_kernel(path, expected_point, nfaces):
    """Read and fully validate one small kernel file (tests and tools).

    Returns ``(geometry, records, sha256)``. The streaming loader in
    ``device_store`` applies the same checks without holding a whole cache in
    host memory.
    """
    path = Path(path)
    raw = path.read_bytes()
    header, count, geometry = validate_kernel_header(
        raw[:KERNEL_HEADER_DTYPE.itemsize], len(raw), path, expected_point, nfaces)
    records = np.frombuffer(raw, dtype=KERNEL_RECORD_DTYPE,
                            offset=KERNEL_HEADER_DTYPE.itemsize, count=count)
    validate_records(records, nfaces, path, 0, np.zeros(nfaces, dtype=bool))
    return geometry, records, hashlib.sha256(raw).hexdigest()


def kernel_set_digest(entries):
    """Combined kernel-set digest of ``cops_io.validate_fit_base``.

    ``entries`` is an iterable of ``(file_name, sha256_hex)`` in table row
    order. Over the full table this is the ``set_sha256`` stored in
    ``run_manifest.json``; over a subset it identifies the loaded subset.
    """
    digest = hashlib.sha256()
    for name, sha in entries:
        digest.update(name.encode("ascii"))
        digest.update(sha.encode("ascii"))
    return digest.hexdigest()

