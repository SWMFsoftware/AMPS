"""Registry of the CUDA C sources and their compile identities (CUDA env only).

Naming: in this package a *CUDA kernel* is a ``__global__`` function in
``cops_cuda/kernels/*.cu``. A *transport kernel* is the per-point Kn/KPhi/KV
data of the CPU application, handled by ``device_store``. The two meanings
never share a module.

Every module is compiled at run time by NVRTC through ``cupy.RawModule`` from
a source file in ``kernels/``. Sources are never assembled from Python
strings at the call site, so the bytes that run on the device are exactly
those a manifest identifies. For every module compiled in the process the
registry records:

* the source file's SHA-256 and size;
* the NVRTC options, checked against ``precision.check_nvrtc_options``;
* the NVRTC version and device compute capability. The same source can
  produce different machine code under a different compiler or for a
  different architecture.

``loaded_identities()`` returns these records for manifests. CuPy caches
compiled binaries on disk, keyed by source text, options and compiler
version, so a source edit always triggers a recompile; the cache is an
optimization, not an identity.
"""

import hashlib
from pathlib import Path

from . import precision


KERNEL_DIRECTORY = Path(__file__).resolve().parent / "kernels"

_modules = {}


def source_path(name):
    """Path of ``kernels/<name>.cu``; only plain module names are accepted."""
    if not name.replace("_", "").isalnum():
        raise ValueError("invalid CUDA source name %r" % name)
    path = KERNEL_DIRECTORY / (name + ".cu")
    if not path.is_file():
        raise FileNotFoundError("CUDA source not found: %s" % path)
    return path


def source_identity(name):
    path = source_path(name)
    data = path.read_bytes()
    return {"name": name, "path": str(path), "size_bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def module(name, device_id=0, options=precision.REFERENCE_NVRTC_OPTIONS):
    """Compile (or reuse) ``kernels/<name>.cu`` for one device.

    Modules are cached per (name, device, options) in this process. A loaded
    module belongs to the device's context, so a module compiled for device 0
    is never reused on device 1.
    """
    import cupy
    from cupy.cuda import nvrtc
    options = precision.check_nvrtc_options(options)
    key = (name, int(device_id), options)
    entry = _modules.get(key)
    if entry is None:
        identity = source_identity(name)
        code = source_path(name).read_bytes().decode("utf-8")
        with cupy.cuda.Device(int(device_id)) as device:
            compiled = cupy.RawModule(code=code, options=options)
            # Compilation is lazy in CuPy; force it now so that compile errors
            # surface here and the identity describes a module that exists.
            compiled.compile()
            identity.update({"device_id": int(device_id), "options": list(options),
                             "nvrtc": "%d.%d" % nvrtc.getVersion(),
                             "compute_capability": device.compute_capability})
        entry = (compiled, identity)
        _modules[key] = entry
    return entry


def function(name, function_name, device_id=0):
    """One ``__global__`` function of a registered module."""
    compiled, _ = module(name, device_id)
    return compiled.get_function(function_name)


def loaded_identities():
    """Identities of every module compiled in this process (for manifests)."""
    return [dict(identity) for _, identity in
            sorted(_modules.values(), key=lambda item: (item[1]["name"], item[1]["device_id"]))]


def all_source_identities():
    """Identities of every source file shipped in ``kernels/``."""
    return [source_identity(path.stem) for path in sorted(KERNEL_DIRECTORY.glob("*.cu"))]
