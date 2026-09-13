#!/usr/bin/env python3
"""Independent vector geometry and DSA algebra for VP07."""
import math
B = (4.0, 3.0, 1.0)
def theta_bn_deg(normal_angle_deg: float) -> float:
    """Return the acute field/normal angle from an explicit dot product."""
    a=math.radians(normal_angle_deg); n=(math.cos(a),math.sin(a),0.0)
    cosine=abs(sum(x*y for x,y in zip(B,n)))/math.sqrt(sum(x*x for x in B))
    return math.degrees(math.acos(min(1.0,max(0.0,cosine))))
def dsa_q(compression: float) -> float:
    """Test-particle phase-space slope q=3r/(r-1)."""
    if compression <= 1.0: raise ValueError("compression must exceed one")
    return 3.0*compression/(compression-1.0)

