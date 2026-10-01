"""Continuous polynomial wind gates with directed interval certificates.

This is an explicit bounded polynomial-interpolant family, not a dense-grid
plausibility heuristic. Every arithmetic enclosure is rounded outward. An
uncertifiable extremum fails rather than being clipped or silently sampled.
"""
from __future__ import annotations
import math
import struct
from .core import require, finite, digest, freeze_record, verify_frozen
from .research import content_hash


def outward(value, direction):
    """One IEEE-754 ULP outward, including Python 3.7/3.8 login nodes.

    ``math.nextafter`` first appeared in Python 3.9. Manipulating the binary64
    representation provides the same directed rounding without a dependency.
    """
    if math.isnan(value) or value == direction: return value
    if value == 0: return math.copysign(struct.unpack('>d', struct.pack('>Q', 1))[0], direction)
    bits = struct.unpack('>Q', struct.pack('>d', value))[0]
    bits += 1 if (direction > value) == (value > 0) else -1
    return struct.unpack('>d', struct.pack('>Q', bits))[0]


def add(a, b):
    return (outward(a[0]+b[0], -math.inf), outward(a[1]+b[1], math.inf))


def multiply(a, b):
    values = [x*y for x in a for y in b]
    return (outward(min(values), -math.inf), outward(max(values), math.inf))


def evaluate(coefficients, interval):
    out = (0.0, 0.0)
    for coefficient in reversed(coefficients):
        out = add(multiply(out, interval), (coefficient, coefficient))
    return out


def derivative(coefficients):
    return [j*x for j,x in enumerate(coefficients)][1:] or [0.0]


def polynomial_product(a, b):
    out = [0.0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j] += x*y
    return out


def certified_extrema(coefficients, support, tolerance):
    finite(coefficients+support+[tolerance])
    require(support[0] < support[1] and tolerance > 0, "invalid extrema support/SI enclosure tolerance")
    slope = derivative(coefficients)
    pending = [(support[0],support[1])]; certified = []; nodes = 0
    while pending:
        lo, hi = pending.pop(); nodes += 1
        require(nodes < 200000, "interval extrema certificate budget exceeded")
        d = evaluate(slope, (lo,hi))
        if d[0] > 0 or d[1] < 0:
            a,b = evaluate(coefficients,(lo,lo)),evaluate(coefficients,(hi,hi))
            certified.append((min(a[0],b[0]),max(a[1],b[1]),lo,hi))
        else:
            value = evaluate(coefficients,(lo,hi))
            if value[1]-value[0] <= tolerance:
                certified.append((value[0],value[1],lo,hi))
            else:
                mid = (lo+hi)/2
                require(lo < mid < hi, "roundoff prevents requested extrema certificate")
                pending.extend([(lo,mid),(mid,hi)])
    minimum = min(certified,key=lambda item:item[0]); maximum = max(certified,key=lambda item:item[1])
    # Certified monotonic cells use endpoint values; critical cells use their
    # full interval enclosure. All cells cover the original support exactly.
    endpoints = [evaluate(coefficients,(x,x)) for cell in certified for x in (cell[2],cell[3])]
    minimum_upper = min(x[1] for x in endpoints); maximum_lower = max(x[0] for x in endpoints)
    require(minimum_upper-minimum[0] <= tolerance and maximum[1]-maximum_lower <= tolerance,
            "global extrema enclosure exceeds its declared SI tolerance")
    return dict(minimum_interval=[minimum[0],minimum_upper],maximum_interval=[maximum_lower,maximum[1]],
        minimum_location=[minimum[2],minimum[3]],maximum_location=[maximum[2],maximum[3]],
        certificate="directed-interval-branch-and-bound-v1",cells=nodes)


def gate_wind(asset, profiles):
    verify_frozen(asset)
    require(asset["schema"] == "sccm-wind-envelope-v6" and asset["data_use_role"] in {"construction","qualification","withheld-validation"}, "wrong wind envelope role/schema")
    require(asset["covariance_sha256"] and asset["inference_provenance_sha256"], "missing envelope uncertainty/provenance")
    content_hash(asset["covariance_sha256"]);content_hash(asset["inference_provenance_sha256"])
    channels = asset["channels"]
    require(len(channels) == len({c["id"] for c in channels}), "duplicate wind channel")
    known = {c["id"]:c for c in channels}
    gates = asset["gates"]
    require(len(gates) == len({g["selector"] for g in gates}) and set(g["selector"] for g in gates) == set(p["selector"] for p in profiles), "missing/ambiguous topology or exact-ID gate")
    rows=[]
    for profile in profiles:
        gate=next(g for g in gates if g["selector"]==profile["selector"])
        require(gate["velocity_channel"] in known and gate["acceleration_channel"] in known, "cross-asset/missing paired channels")
        velocity,acceleration=known[gate["velocity_channel"]],known[gate["acceleration_channel"]]
        require(velocity["quantity"] == "field-aligned-speed" and velocity["units"] == "m/s" and
                acceleration["quantity"] == "quasi-steady-field-aligned-advective-acceleration" and acceleration["units"] == "m/s^2", "incompatible single-quantity channel pair")
        require(profile["quasi_steady"] is True, "full material derivative requires its own time-dependent authority")
        require(velocity["frame"] == acceleration["frame"] == profile["frame"] and
                velocity["support"] == acceleration["support"] == profile["support"], "wind frame/support/channel pair mismatch")
        coefficients=profile["speed_coefficients"]
        finite(coefficients+[profile["unsigned_flux_wb"]]);require(profile["unsigned_flux_wb"]>=0 and
            profile["tube_id"] and profile["segment_id"],"invalid stable wind-support/flux identity")
        advective=polynomial_product(coefficients,derivative(coefficients))
        pair=[]
        for channel,values in ((velocity,coefficients),(acceleration,advective)):
            finite([channel["lower_bound"],channel["upper_bound"],channel["enclosure_tolerance_si"]])
            require(channel["interpolant"] == "polynomial-in-s-SI-v1" and channel["lower_bound"] <= channel["upper_bound"], "unsupported wind interpolant/bounds")
            certificate=certified_extrema(values,profile["support"],channel["enclosure_tolerance_si"])
            low=max(0,channel["lower_bound"]-certificate["minimum_interval"][0]); high=max(0,certificate["maximum_interval"][1]-channel["upper_bound"])
            pair.append(dict(channel_id=channel["id"],certificate=certificate,lower_excursion_si=low,upper_excursion_si=high,passed=low==high==0))
        passed=all(row["passed"] for row in pair)
        rows.append(dict(selector=profile["selector"],stable_tube_id=profile["tube_id"],stable_segment_id=profile["segment_id"],
            channels=pair,passed=passed,rejected_magnetic_flux_wb=0 if passed else profile["unsigned_flux_wb"]))
    return freeze_record(dict(schema="sccm-wind-envelope-result-v6",asset_identity=asset["identity"],rows=rows,passed=all(row["passed"] for row in rows)))
