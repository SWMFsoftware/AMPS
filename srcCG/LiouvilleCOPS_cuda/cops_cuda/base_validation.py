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
import json
import os
from pathlib import Path

import numpy as np

POLICY = "numpy_workers_in_reference_validate_fit_base_v1"
CACHE_POLICY = "validate_once_v1"   # user decision 2026-10-21 (SCALING_PLAN Q4)
CACHE_SCHEMA = "liouville_cops_cuda.base_audit_cache"
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


def _stat_key(path):
    st = Path(path).stat()
    return [int(st.st_size), int(st.st_mtime_ns), int(st.st_ino)]


class AuditCache:
    """Validate once (``CACHE_POLICY``): per-file worker results reused while a file is unchanged.

    The cache maps each kernel file (and ``elements.bin``) to its (size, mtime_ns, inode) and
    the result of the worker that validated it the first time (full read, record checks,
    SHA-256). A later run reuses a result only when the file's stat key and the call's
    arguments (point, face count) are unchanged; anything else is validated again. The
    reference ``validate_fit_base`` still runs in full on these results (table, header
    comparisons, kernel-set digest against ``run_manifest.json``), so the audit it returns
    is the one a full validation returns for unchanged files.
    """

    def __init__(self, path, revalidate=False):
        self.path = Path(path)
        self.entries, self.hits, self.misses = {}, 0, 0
        if self.path.is_file() and not revalidate:
            data = json.loads(self.path.read_text())
            if data.get("schema") == CACHE_SCHEMA and data.get("policy") == CACHE_POLICY:
                self.entries = data["entries"]

    def get(self, path, args):
        entry = self.entries.get(str(Path(path).resolve()))
        if entry is not None and entry["stat"] == _stat_key(path) and entry["args"] == list(args):
            self.hits += 1
            result = dict(entry["result"])
            if "header" in result:
                result["header"] = tuple(result["header"])
            return result
        return None

    def put(self, path, args, stat, result):
        self.misses += 1
        stored = dict(result)
        if "header" in stored:
            stored["header"] = list(stored["header"])
        self.entries[str(Path(path).resolve())] = {"stat": stat, "args": list(args), "result": stored}

    def save(self, run):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".tmp%d" % os.getpid())
        tmp.write_text(json.dumps({"schema": CACHE_SCHEMA, "policy": CACHE_POLICY,
                                   "run": str(Path(run).resolve()), "entries": self.entries}))
        os.replace(tmp, self.path)

    def summary(self):
        return {"policy": CACHE_POLICY, "path": str(self.path), "reused": self.hits,
                "validated": self.misses}


class _Prefetch:
    """Kernel validations computed ahead on threads (hashlib and NumPy release the GIL).

    The reference calls ``validate_kernel(path, point, nfaces)`` in table
    order.  Once ``validate_elements`` gives ``nfaces``, every
    ``point=NNNN.bin`` is submitted with the point its name encodes; a call is
    answered from that work only for the identical (path, point, nfaces), and
    computed directly otherwise.  Exceptions are re-raised at the call, so
    the first error is still the one the reference order meets first.
    With an ``AuditCache``, unchanged files are answered from it without reading them.
    """

    def __init__(self, io, threads, cache=None):
        self.io, self.threads, self.cache = io, int(threads), cache
        self.pool, self.futures, self.nfaces = None, {}, None

    def _cached(self, fn, path, *args):
        if self.cache is None:
            return fn(self.io, path, *args)
        hit = self.cache.get(path, args)
        if hit is not None:
            return hit
        stat = _stat_key(path)
        result = fn(self.io, path, *args)
        if _stat_key(path) == stat:                  # unchanged while it was read
            self.cache.put(path, args, stat, result)
        return result

    def elements(self, path):
        info = self._cached(validate_elements, path)
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
                self.futures[key] = self.pool.submit(self._cached, validate_kernel, kernel, point,
                                                     self.nfaces)
        return info

    def kernel(self, path, expected_point, nfaces):
        future = self.futures.pop((str(Path(path)), expected_point, nfaces), None)
        if future is None:
            return self._cached(validate_kernel, path, expected_point, nfaces)
        return future.result()

    def close(self):
        if self.pool is not None:
            self.pool.shutdown(wait=True, cancel_futures=True)


def validate_fit_base(io, run, threads=8, cache=None):
    """The reference ``io.validate_fit_base(run)`` with vectorized workers.

    ``cache``: an ``AuditCache`` (validate once); saved after a successful validation only.
    """
    _check_structs(io)
    prefetch = _Prefetch(io, threads, cache)
    saved = io.validate_kernel, io.validate_elements
    io.validate_kernel, io.validate_elements = prefetch.kernel, prefetch.elements
    try:
        audit = io.validate_fit_base(run)
    finally:
        io.validate_kernel, io.validate_elements = saved
        prefetch.close()
    if cache is not None:
        cache.save(run)
    return audit


def default_cache_path(run):
    """``$COPS_OUTPUT_ROOT/cuda_cache/base_audit/<base name>-<path digest>.json``, or None."""
    root = os.environ.get("COPS_OUTPUT_ROOT")
    if not root:
        return None
    run = Path(run).resolve()
    tag = hashlib.sha256(str(run).encode()).hexdigest()[:12]
    return Path(root) / "cuda_cache" / "base_audit" / ("%s-%s.json" % (run.name, tag))


def install(fit_cops, cache="auto", revalidate=False):
    """Make ``fit_cops.main`` use the vectorized validation; returns a record for the manifest.

    ``cache``: "auto" (``default_cache_path`` of the base being validated), "off", or a path.
    """
    import sys
    io = sys.modules[fit_cops.validate_fit_base.__module__]
    reference = fit_cops.validate_fit_base
    record = {"policy": POLICY, "reference": "%s.%s" % (reference.__module__, reference.__name__),
              "module": str(Path(__file__).resolve()), "cache": None}

    def run_validation(run):
        path = (default_cache_path(run) if cache == "auto" else None if cache in (None, "off")
                else Path(cache))
        audit_cache = None if path is None else AuditCache(path, revalidate)
        audit = validate_fit_base(io, run, cache=audit_cache)
        if audit_cache is not None:
            record["cache"] = audit_cache.summary()
        return audit

    fit_cops.validate_fit_base = run_validation
    return record
