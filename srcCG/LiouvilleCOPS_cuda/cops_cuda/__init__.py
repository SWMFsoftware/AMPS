"""CUDA/CuPy backend development package for LiouvilleCOPS.

Modules in this package must stay importable by the CPU reference
interpreter (Python 3.8) unless their docstring states otherwise, because the
C00 reference-capture tools share the provenance and reference-binding code
with the CUDA environment.  CuPy is imported only by modules that need a GPU.
"""
