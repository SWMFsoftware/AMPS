#!/usr/bin/env python3
"""Independent relativistic momentum-spectrum oracle for VP12."""
import math
REST_MEV=938.27208816
def momentum_proxy(energy_mev: float) -> float:
    return math.sqrt(energy_mev*(energy_mev+2.0*REST_MEV))
def intensity_ratio(energy_mev: float, reference_mev: float, q: float) -> float:
    return (momentum_proxy(energy_mev)/momentum_proxy(reference_mev))**(2.0-q)
def production_equivalent_ratio(energy_mev: float, reference_mev: float, q: float) -> float:
    """Evaluate the same physical law through log space as an independent path."""
    return math.exp((2.0-q)*(math.log(momentum_proxy(energy_mev))-math.log(momentum_proxy(reference_mev))))

