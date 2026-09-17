#!/usr/bin/env python3
"""Controlled focused-transport Green-function profile for VP15."""
import math
def profile(time_h: float, energy_mev: float, q: float) -> float:
    if time_h<=0 or energy_mev<=0: raise ValueError("positive time and energy required")
    return time_h**-1.5*math.exp(-0.7/time_h-time_h/24.0)*(energy_mev/20.0)**(-(q-2.0)/2.0)
def profile_log_form(time_h: float, energy_mev: float, q: float) -> float:
    """Log-space evaluation avoids sharing the direct multiplication path."""
    return math.exp(-1.5*math.log(time_h)-0.7/time_h-time_h/24.0-0.5*(q-2.0)*math.log(energy_mev/20.0))

