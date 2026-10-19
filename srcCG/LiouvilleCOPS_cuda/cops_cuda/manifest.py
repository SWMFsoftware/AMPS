"""Backend/software/device identity sections for manifests (ROADMAP C02 item 5).

Every product of the CUDA backend records enough to reproduce or refute it:

* ``software``: content digests of this package's tree and of the bound CPU
  reference tree, plus Python, NumPy, SciPy and CuPy versions. A dirty or
  unversioned tree is still identified exactly, because digests cover the
  file contents and not version-control state.
* ``device``: GPU UUID, PCI bus, name, compute capability, memory and ECC
  state of the device actually used. The device index alone is not stable
  across ``CUDA_VISIBLE_DEVICES`` settings.
* ``runtime``: the CUDA runtime, driver API and NVRTC versions in this
  process.
* ``precision``: the declared FP64 policy (``precision.describe``).
* ``cuda_sources``: every CUDA module compiled in this process, with its
  source digest, NVRTC options and target architecture
  (``cuda_sources.loaded_identities``), and the digest of every shipped
  ``.cu`` file.

``backend_manifest`` composes these with a dataset's ``identity()``. Writers
use ``provenance.write_json_exclusive``, so an existing manifest is never
replaced.
"""

from pathlib import Path

from . import cuda_sources, environment, precision, provenance, reference


MANIFEST_SCHEMA = "liouville_cops_cuda.backend_manifest"
MANIFEST_SCHEMA_VERSION = 1
PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def software_identity(include_reference=True):
    """Digests of this application and the bound reference, plus versions."""
    import numpy
    import scipy
    import sys
    package = provenance.tree_identity(PACKAGE_ROOT, exclude_dirs=("third_party",))
    result = {"application_root": str(PACKAGE_ROOT),
              "application_tree_sha256": package["tree_sha256"],
              "application_file_count": package["file_count"],
              "python": sys.version.split()[0], "numpy": numpy.__version__,
              "scipy": scipy.__version__}
    try:
        import cupy
        result["cupy"] = cupy.__version__
    except ImportError:
        result["cupy"] = None
    if include_reference:
        root = reference.bind()
        tree = provenance.tree_identity(root, exclude_dirs=reference.SOURCE_EXCLUDES)
        result.update({"reference_root": str(root),
                       "reference_tree_sha256": tree["tree_sha256"],
                       "reference_file_count": tree["file_count"]})
    return result


def device_identity(device_id):
    """Inventory entry of one device (stable UUID, not just the index)."""
    for entry in environment.device_inventory():
        if entry["index"] == int(device_id):
            return entry
    raise ValueError("no CUDA device with index %d" % int(device_id))


def backend_manifest(backend_name, device_id, dataset=None, extra=None):
    """Complete identity section for one backend instance on one device."""
    value = {"schema": MANIFEST_SCHEMA, "schema_version": MANIFEST_SCHEMA_VERSION,
             "created_utc": provenance.utc_now(), "backend": backend_name,
             "software": software_identity(), "device": device_identity(device_id),
             "runtime": environment.runtime_versions(), "precision": precision.describe(),
             "cuda_sources": {"compiled": cuda_sources.loaded_identities(),
                              "shipped": cuda_sources.all_source_identities()},
             "process": provenance.process_context()}
    if dataset is not None:
        value["dataset"] = dataset.identity()
        if dataset.budget is not None and not dataset.released:
            value["memory_state"] = dataset.budget.state()
    if extra:
        value.update(extra)
    return value
