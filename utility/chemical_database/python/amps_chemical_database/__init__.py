"""Application-independent, offline AMPS chemical catalogue access.

The numerical archive and its pinned runtime export are immutable inputs.
This package queries authenticated source tables; it does not infer reaction
kinetics, product states, illumination, plasma distributions or material laws.
Use the shared C++ production services for supported numerical evaluation.
"""

from .api import Database, DatabaseError, Table

__all__ = ["Database", "DatabaseError", "Table"]
__version__ = "1.0.0"
