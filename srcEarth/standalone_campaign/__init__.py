"""Reproducible standalone event campaigns for the SEP-in-geospace model.

Roadmap Step 8 deliberately lives above the Step 2--7 physics kernels.  The
package prepares immutable inputs, launches the existing standalone executable,
and evaluates declared evidence; it never changes a trajectory, quadrature, or
validation tolerance in order to make an event comparison pass.
"""

from .campaign import CampaignError, load_and_validate_manifest

__all__ = ["CampaignError", "load_and_validate_manifest"]
