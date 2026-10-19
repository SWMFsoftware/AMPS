"""CUDA environment discovery for the LiouvilleCOPS GPU backend (CUDA env only).

Everything reported here is observed at run time from the loaded libraries and
devices.  Paths or versions from documentation are never substituted: a found
shared library is not evidence that a kernel compiles or runs, so device
execution is tested separately by ``tools/cuda_doctor.py``.

Thread policy
-------------
The host has 36 logical CPUs and the user prefers 36 threads for CPU work
that is run alone (C++ replay, CPU gravity reference).  GPU workers are
different: each persistent worker process would otherwise start its own
OpenBLAS/OpenMP pool sized to all CPUs, oversubscribing the host when several
workers or CPU jobs run together.  ``worker_thread_environment`` gives the
variables a spawned GPU worker must set before NumPy is imported.
"""

import importlib.metadata
import os
from pathlib import Path


THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS")


def worker_thread_environment(threads=1):
    """Environment for a spawned GPU worker: small, explicit host thread pools."""
    return {name: str(int(threads)) for name in THREAD_VARIABLES}


def cupy_distributions():
    """All installed distributions that provide CuPy; exactly one is valid.

    Installing two CuPy wheels (e.g. cupy-cuda11x and cupy-cuda12x) leaves one
    ``cupy`` package overwriting the other's files, which fails in confusing
    ways; detect it explicitly instead.
    """
    found = []
    for dist in importlib.metadata.distributions():
        name = (dist.metadata.get("Name") or "").lower()
        if name == "cupy" or name.startswith("cupy-") or name.startswith("cupy_"):
            found.append({"name": name, "version": dist.version})
    return sorted(found, key=lambda item: item["name"])


def loaded_cuda_libraries():
    """Resolved paths of CUDA shared libraries mapped into this process."""
    wanted = ("libcudart", "libnvrtc", "libcublas", "libcusparse", "libcusolver",
              "libcurand", "libcufft", "libcuda.so", "libnvJitLink", "libnccl")
    found = set()
    try:
        for line in Path("/proc/self/maps").read_text().splitlines():
            parts = line.split()
            if len(parts) >= 6 and any(key in parts[5] for key in wanted):
                found.add(parts[5])
    except OSError:
        return None
    return sorted(found)


def format_uuid(raw):
    """Format the 16-byte CUDA device UUID like nvidia-smi (GPU-xxxxxxxx-...)."""
    if isinstance(raw, str):
        raw = raw.encode("latin-1")
    hexed = bytes(raw).hex()
    if len(hexed) != 32:
        return None
    return "GPU-%s-%s-%s-%s-%s" % (hexed[0:8], hexed[8:12], hexed[12:16],
                                   hexed[16:20], hexed[20:32])


def device_inventory():
    """Stable identities and capacities of every visible CUDA device."""
    import cupy
    runtime = cupy.cuda.runtime
    devices = []
    for index in range(runtime.getDeviceCount()):
        props = runtime.getDeviceProperties(index)
        with cupy.cuda.Device(index) as device:
            free, total = device.mem_info
        name = props["name"]
        devices.append({
            "index": index,
            "name": name.decode() if isinstance(name, bytes) else str(name),
            "uuid": format_uuid(props.get("uuid", b"")),
            "pci_bus_id": runtime.deviceGetPCIBusId(index),
            "compute_capability": "%d.%d" % (props["major"], props["minor"]),
            "multiprocessors": props["multiProcessorCount"],
            "total_memory_bytes": int(total), "free_memory_bytes": int(free),
            "ecc_enabled": bool(props.get("ECCEnabled", 0)),
            "unified_addressing": bool(props.get("unifiedAddressing", 0)),
        })
    return devices


def runtime_versions():
    """CUDA runtime/driver/NVRTC versions as reported by the loaded libraries."""
    import cupy
    from cupy.cuda import nvrtc

    def version(value):
        return "%d.%d" % (value // 1000, (value % 1000) // 10)

    return {
        "cupy": cupy.__version__,
        "cuda_path": cupy.cuda.get_cuda_path(),
        "nvcc_path": cupy.cuda.get_nvcc_path(),
        "runtime": version(cupy.cuda.runtime.runtimeGetVersion()),
        "driver_api": version(cupy.cuda.runtime.driverGetVersion()),
        "nvrtc": "%d.%d" % nvrtc.getVersion(),
        "cupy_cache_dir": os.environ.get("CUPY_CACHE_DIR",
                                         str(Path.home() / ".cupy" / "kernel_cache")),
        "distributions": cupy_distributions(),
    }
