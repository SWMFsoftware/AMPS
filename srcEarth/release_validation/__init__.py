"""Roadmap Step-12 cross-path and release validation package.

The package contains no AMPS, MPI, or SWMF dependency.  It validates immutable
artifacts emitted by those builds, which makes the final decision reproducible on a
review machine while leaving every scientific acceptance threshold unchanged.
"""

from .contract import COMMON_PHYSICS_TAG, ReleaseValidationError

__all__ = ["COMMON_PHYSICS_TAG", "ReleaseValidationError"]
