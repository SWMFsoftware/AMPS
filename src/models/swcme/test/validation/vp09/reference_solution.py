#!/usr/bin/env python3
"""Independent ordered-region classifier for VP09."""
def classify(time: float) -> int:
    """Map shock-relative normalized time to ambient/sheath/ejecta codes."""
    if time < 0.0 or time >= 2.0: return 0
    if time < 1.0: return 1
    return 2
def region_name(code: int) -> str:
    return {0:"ambient",1:"sheath",2:"ejecta"}[code]

