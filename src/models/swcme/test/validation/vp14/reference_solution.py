#!/usr/bin/env python3
"""Independent normalized Gaussian diffusion kernel for VP14."""
import math
def gaussian_spread_ratio(longitude_deg: float, diffusion_variance: float) -> float:
    phi=math.radians(longitude_deg)
    return math.exp(-phi*phi/(4.0*diffusion_variance))

