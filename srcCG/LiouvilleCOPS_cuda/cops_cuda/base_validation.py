"""Vectorized fitting-base validation, identical in result to the reference.

The reference ``cops_io.validate_fit_base`` (called by ``fit_cops.main``
before any fit) unpacks every kernel record with ``struct`` in Python.  On
SHAP8 12M (4.3e9 records, 173 GB) that takes about 1.8 h per fit, as long as
the GPU fit itself.  ``validate_fit_base`` here runs the *reference function
itself* with only its two per-byte workers, ``cops_io.validate_kernel`` and
``cops_io.validate_elements``, swapped for NumPy equivalents for the duration
of the call.  The table/kernel comparisons, manifest checks and the returned
audit therefore come from the reference code; the workers perform the same
checks in the same order, raise the same first error with the same message,
and hash the same bytes (``tests/test_base_validation.py`` compares both,
including corrupted files).  ``cops_io.py`` and ``fit_cops.py`` are not
modified; ``gpu_fit.py`` installs this in its own process only.
"""

import functools
import hashlib
from pathlib import Path

import numpy as np

POLICY = "numpy_workers_in_reference_validate_fit_base_v1"
KERNEL_RECORD = np.dtype({"names": ["face", "kn", "kphi", "kv", "cos", "shadow"],
                          "formats": ["<i4", "<f8", "<f8", "<f8", "<f8", "<i4"],
                          "offsets": [0, 4, 12, 20, 28, 36], "itemsize": 40})
ELEMENT_FLOATS = 34
ELEMENT_CHUNK = 1 << 20                # faces per read


def _check_structs(io):
    if io.KERNEL1_RECORD.format not in ("<i4di", b"<i4di") or io.KERNEL1_RECORD.size != 40:
        raise RuntimeError("reference kernel record layout changed: %r" % io.KERNEL1_RECORD.format)
    if io.ELEMENT_HEADER.format not in ("<iii", b"<iii"):
        raise RuntimeError("reference element header changed: %r" % io.ELEMENT_HEADER.format)


def validate_elements(io, path):
    """``cops_io.validate_elements`` with the per-face loop vectorized."""
    path = Path(path)
    size = path.stat().st_size if path.is_file() else -1
    with path.open("rb") as stream:
        header = stream.read(io.ELEMENT_HEADER.size)
        if len(header) != io.ELEMENT_HEADER.size:
            raise io.ProductError("truncated elements header: %s" % path)
        version, nfaces, nbasis = io.ELEMENT_HEADER.unpack(header)
        if version != 1 or nfaces <= 0 or nbasis != 25:
            raise io.ProductError("unsupported elements header in %s: version=%d faces=%d basis=%d" %
                                  (path, version, nfaces, nbasis))
        expected = io.ELEMENT_HEADER.size + nfaces * ELEMENT_FLOATS * 8
        if size != expected:
            raise io.ProductError("elements size mismatch in %s: %d bytes; expected %d" %
                                  (path, size, expected))
        digest = hashlib.sha256(header)
        for start in range(0, nfaces, ELEMENT_CHUNK):
            n = min(ELEMENT_CHUNK, nfaces - start)
            raw = stream.read(n * ELEMENT_FLOATS * 8)
            digest.update(raw)
            values = np.frombuffer(raw, dtype="<f8").reshape(n, ELEMENT_FLOATS)
            finite = np.isfinite(values).all(axis=1)
            bad = ~finite | ~(values[:, 6] > 0.0)
            if bad.any():
                i = int(np.argmax(bad))
                if not finite[i]:
                    raise io.ProductError("non-finite value in %s face %d" % (path, start + i))
                raise io.ProductError("non-positive area in %s face %d" % (path, start + i))
    return {"path": str(path.resolve()), "size_bytes": size,
            "sha256": digest.hexdigest(), "version": version,
            "nfaces": nfaces, "nbasis": nbasis}


def _first_duplicates(face):
    """True at every record whose face occurred at an earlier record."""
    order = np.argsort(face, kind="stable")
    s = face[order]
    later = np.empty(len(s), dtype=bool)
    later[:1] = False
    later[1:] = s[1:] == s[:-1]
    dup = np.zeros(len(face), dtype=bool)
    dup[order] = later
    return dup


def validate_kernel(io, path, expected_point, nfaces):
    """``cops_io.validate_kernel`` with the per-record loop vectorized.

    The reference raises at the first record failing any of (range, duplicate,
    finite, shadow), checked in that order within a record; every earlier
    record has then passed, so the earliest flagged record and its first
    failing check give the identical error.
    """
    path = Path(path)
    size = path.stat().st_size if path.is_file() else -1
    with path.open("rb") as stream:
        raw_header = stream.read(io.KERNEL1_HEADER.size)
        if len(raw_header) != io.KERNEL1_HEADER.size:
            raise io.ProductError("truncated kernel header: %s" % path)
        unpacked = io.KERNEL1_HEADER.unpack(raw_header)
        version, point = unpacked[:2]
        values, count = unpacked[2:12], unpacked[12]
        if version != 1:
            raise io.ProductError("%s has kernel version %d; static fitting needs version 1" %
                                  (path, version))
        if point != expected_point:
            raise io.ProductError("%s stores iPoint %d; expected %d" %
                                  (path, point, expected_point))
        if count < 0 or count > nfaces:
            raise io.ProductError("%s has invalid record count %d for %d faces" %
                                  (path, count, nfaces))
        expected_size = io.KERNEL1_HEADER.size + count * io.KERNEL1_RECORD.size
        if size != expected_size:
            raise io.ProductError("kernel size mismatch in %s: %d bytes; expected %d" %
                                  (path, size, expected_size))
        io._finite(values, "%s header" % path)
        raw = stream.read(count * io.KERNEL1_RECORD.size)
    digest = hashlib.sha256(raw_header)
    digest.update(raw)
    rec = np.frombuffer(raw, dtype=KERNEL_RECORD, count=count)
    face = rec["face"]
    out_of_range = (face < 0) | (face >= nfaces)
    floats = np.stack([rec["kn"], rec["kphi"], rec["kv"], rec["cos"]], axis=1)
    nonfinite = ~np.isfinite(floats).all(axis=1)
    bad_shadow = (rec["shadow"] != 0) & (rec["shadow"] != 1)
    counts = np.bincount(face[~out_of_range], minlength=nfaces)
    duplicate = (_first_duplicates(face) & ~out_of_range if (counts > 1).any()
                 else np.zeros(count, dtype=bool))
    bad = out_of_range | duplicate | nonfinite | bad_shadow
    if bad.any():
        r = int(np.argmax(bad))
        if out_of_range[r]:
            raise io.ProductError("%s record %d has out-of-range face %d" % (path, r, int(face[r])))
        if duplicate[r]:
            raise io.ProductError("%s contains duplicate face %d" % (path, int(face[r])))
        if nonfinite[r]:
            raise io.ProductError("non-finite value in %s record %d" % (path, r))
        raise io.ProductError("%s record %d has invalid shadow flag %d" %
                              (path, r, int(rec["shadow"][r])))
    return {"point": point, "records": count, "size_bytes": size,
            "sha256": digest.hexdigest(), "header": values}


class _Prefetch:
    """Kernel validations computed ahead on threads (hashlib and NumPy release the GIL).

    The reference calls ``validate_kernel(path, point, nfaces)`` in table
    order.  Once ``validate_elements`` gives ``nfaces``, every
    ``point=NNNN.bin`` is submitted with the point its name encodes; a call is
    answered from that work only for the identical (path, point, nfaces), and
    computed directly otherwise.  Exceptions are re-raised at the call, so
    the first error is still the one the reference order meets first.
    """

    def __init__(self, io, threads):
        self.io, self.threads = io, int(threads)
        self.pool, self.futures, self.nfaces = None, {}, None

    def elements(self, path):
        info = validate_elements(self.io, path)
        self.nfaces = info["nfaces"]
        if self.threads > 1:
            import concurrent.futures as cf
            self.pool = cf.ThreadPoolExecutor(self.threads, thread_name_prefix="cops-validate")
            for kernel in sorted((Path(path).parent / "kernels").glob("point=*.bin")):
                try:
                    point = int(kernel.name[6:-4])
                except ValueError:
                    continue
                key = (str(kernel), point, self.nfaces)
                self.futures[key] = self.pool.submit(validate_kernel, self.io, kernel, point, self.nfaces)
        return info

    def kernel(self, path, expected_point, nfaces):
        future = self.futures.pop((str(Path(path)), expected_point, nfaces), None)
        if future is None:
            return validate_kernel(self.io, path, expected_point, nfaces)
        return future.result()

    def close(self):
        if self.pool is not None:
            self.pool.shutdown(wait=True, cancel_futures=True)


def validate_fit_base(io, run, threads=8):
    """The reference ``io.validate_fit_base(run)`` with vectorized workers."""
    _check_structs(io)
    prefetch = _Prefetch(io, threads)
    saved = io.validate_kernel, io.validate_elements
    io.validate_kernel, io.validate_elements = prefetch.kernel, prefetch.elements
    try:
        return io.validate_fit_base(run)
    finally:
        io.validate_kernel, io.validate_elements = saved
        prefetch.close()


def install(fit_cops):
    """Make ``fit_cops.main`` use the vectorized validation; returns a record for the manifest."""
    import sys
    io = sys.modules[fit_cops.validate_fit_base.__module__]
    reference = fit_cops.validate_fit_base
    fit_cops.validate_fit_base = functools.partial(validate_fit_base, io)
    return {"policy": POLICY, "reference": "%s.%s" % (reference.__module__, reference.__name__),
            "module": str(Path(__file__).resolve())}
