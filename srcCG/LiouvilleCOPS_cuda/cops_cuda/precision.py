"""Precision policy of the reference GPU backend (Python 3.8 compatible).

ROADMAP C02 item 4 requires the following until C08 validates anything else:

* mesh geometry, parameters and baseline references stay FP64;
* every transport-kernel contribution field is staged in FP64.

The policy is a declared object, not a convention. It is recorded in every
manifest, and each device field is checked against it when it is
allocated. That way an accidental FP32 copy (``cp.asarray`` of a float32
array, for example) is refused instead of silently lowering the precision of
a comparison.

The integer and Boolean fields are lossless representations of the file
content, not precision choices:

* ``face``: int32, because the file stores int32 and the range is checked.
* ``shadow``: uint8, because the file stores int32 0/1 and the value is
  checked to be 0 or 1 before narrowing.
* CSR ``offsets``: int64. ROADMAP C02 item 2 requires 64-bit offsets, since
  the campaign already has 2.46e8 records, and larger caches can exceed
  2^31.

NVRTC options of reference kernels
----------------------------------
``--fmad=false`` stops the compiler from contracting ``a*b+c`` into one fused
multiply-add. NumPy evaluates the reference products and sums as separately
rounded operations. Without contraction each per-term product is the same
correctly rounded IEEE result as on the CPU, and only the reduction order
differs. Fast-math (``--use_fast_math``) is forbidden: it flushes denormals
and replaces division, square root and transcendental functions with
approximations. A kernel may still call ``fma()`` explicitly where its
documentation justifies that.
"""

import numpy as np


POLICY_NAME = "fp64_reference_v1"
REFERENCE_NVRTC_OPTIONS = ("--std=c++14", "--fmad=false")
FORBIDDEN_NVRTC_OPTIONS = ("--use_fast_math", "-use_fast_math", "--ftz=true",
                           "--prec-div=false", "--prec-sqrt=false", "--fmad=true")

# Device dtype of every field a store may hold.  Anything not listed is
# refused by ``require_device_dtype``.
FIELD_DTYPES = {
    # transport-kernel records (one entry per point-facet contribution)
    "face": np.dtype(np.int32), "kn": np.dtype(np.float64),
    "kphi": np.dtype(np.float64), "kv": np.dtype(np.float64),
    "cos": np.dtype(np.float64), "shadow": np.dtype(np.uint8),
    # per observation point
    "offsets": np.dtype(np.int64), "geometry": np.dtype(np.float64),
    # per mesh facet
    "centroid": np.dtype(np.float64), "normal": np.dtype(np.float64),
    "area": np.dtype(np.float64), "Y": np.dtype(np.float64),
    "AH2O": np.dtype(np.float64), "ACO2": np.dtype(np.float64),
    # per-facet activity basis bound by a fit (harmonics and templates)
    "enhancement_basis": np.dtype(np.float64),
}


def require_device_dtype(name, dtype):
    """Refuse a field whose dtype differs from the declared policy."""
    expected = FIELD_DTYPES.get(name)
    if expected is None:
        raise ValueError("field %r has no declared precision policy" % name)
    if np.dtype(dtype) != expected:
        raise TypeError("field %r must be %s under %s; got %s" %
                        (name, expected, POLICY_NAME, np.dtype(dtype)))
    return expected


def check_nvrtc_options(options):
    """Refuse compiler options that would change reference arithmetic."""
    options = tuple(options)
    bad = [option for option in options if option in FORBIDDEN_NVRTC_OPTIONS]
    if bad:
        raise ValueError("options %s are forbidden by %s" % (", ".join(bad), POLICY_NAME))
    if "--fmad=false" not in options:
        raise ValueError("reference kernels must be compiled with --fmad=false")
    return options


def describe():
    """Machine-readable policy for manifests."""
    return {"name": POLICY_NAME, "arithmetic": "float64",
            "nvrtc_options": list(REFERENCE_NVRTC_OPTIONS),
            "forbidden_nvrtc_options": list(FORBIDDEN_NVRTC_OPTIONS),
            "fields": {name: dtype.str for name, dtype in sorted(FIELD_DTYPES.items())},
            "mixed_precision": "not permitted before C08 validation"}
