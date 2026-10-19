"""Exception classes of the CUDA backend (Python 3.8 compatible).

The categories are separate because they require different responses:

* ``CacheFormatError``: the bytes on disk are malformed (truncated or appended
  data, a wrong version, an out-of-range face index, a duplicate face, a
  non-finite value, an invalid shadow flag).  The cache must be regenerated.
* ``IdentityError``: the bytes are well formed but are not the expected
  product.  Examples: an ``elements.bin`` whose digest disagrees with
  ``run_manifest.json``, a template matrix built for another facet order, or a
  kernel header that disagrees with its table row.  Using such data would
  silently pair a facet or observation with the wrong values.
* ``CapacityError``: the declared memory plan does not fit the device under
  the headroom policy, or an allocation exceeded the plan.  The error message
  carries the full plan so the user can see which array is responsible.
* ``UnsupportedConfiguration``: a backend was asked for a model setting or a
  Jacobian column it does not implement.  It is raised before any evaluation;
  a backend never substitutes a different model.

All derive from ``BackendError``, so callers can catch every backend refusal
at once while tests can still check for the specific category.
"""


class BackendError(Exception):
    """Base class of every refusal raised by the CUDA backend."""


class CacheFormatError(BackendError, ValueError):
    """A cache file is malformed; it is never repaired or partially used."""


class IdentityError(BackendError, ValueError):
    """Well-formed data that is not the product the caller declared."""


class CapacityError(BackendError, MemoryError):
    """The memory plan does not fit, or an allocation exceeded the plan."""


class UnsupportedConfiguration(BackendError, NotImplementedError):
    """A requested model setting or derivative is not implemented."""
