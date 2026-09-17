#!/usr/bin/env python3
"""Independent Parker antiderivative and logarithmic-gradient oracle."""
import math
AU_M=149597870700.0
RS_M=695700000.0
# Use the documented nominal rotation rate as an input constant.  The oracle
# remains independent because all integration/differentiation algebra below is
# separate from the production implementation.
OMEGA=2.86533e-6
VSW=4.0e5
SOURCE=2.5*RS_M
def parker_metrics(sin_theta: float, radius_au: float):
    """Return path and focusing lengths in AU without SWCME imports."""
    r=radius_au*AU_M; k=OMEGA*sin_theta/VSW
    if k==0.0: path=abs(r-SOURCE)
    else:
        def primitive(radius):
            x=radius-SOURCE; kx=k*x
            return 0.5*(x*math.sqrt(1+kx*kx)+math.asinh(kx)/k)
        path=abs(primitive(r)-primitive(SOURCE))
    x=r-SOURCE; z=k*x; one=1+z*z
    focus=r*one*math.sqrt(one)/(2*one-k*k*r*x)
    return path/AU_M,focus/AU_M
