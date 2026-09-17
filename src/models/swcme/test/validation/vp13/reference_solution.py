#!/usr/bin/env python3
"""Scale-aware comparison primitive for VP13."""
def relative_difference(a: float,b: float) -> float:
    return abs(a-b)/max(abs(a),abs(b),1.0e-300)

