#!/usr/bin/env python3
"""Independent finite-width front/intersection oracle for VP08."""
import math
def intersection_time_h(angle_deg: float, half_width_deg: float, apex_time_h: float):
    """Return cosine-flank time, or None when the observer misses the cap."""
    if abs(angle_deg)>half_width_deg: return None
    return apex_time_h/math.cos(math.radians(angle_deg))

