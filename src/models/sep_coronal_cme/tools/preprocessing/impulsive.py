"""Open-footpoint impulsive attribution with independent physical normalization."""
from __future__ import annotations
import math
from .core import require, finite, freeze_record, verify_frozen
from .research import content_hash


def impulsive_plan(asset, total_number, front_generation, birth_side):
    verify_frozen(asset)
    require(asset["schema"] == "sccm-impulsive-source-v6" and asset["intent"] == "attribution-sensitivity", "unknown impulsive source/intent")
    require(asset["observation_response_sha256"] and asset["mapping_covariance_sha256"] and asset["frame"], "source timing response/mapping uncertainty/frame absent")
    content_hash(asset["observation_response_sha256"]);content_hash(asset["mapping_covariance_sha256"])
    require(asset["random_stream"] == "impulsive-only" and asset["normalization_authority"] == "independent-impulsive", "shock/impulsive calibration contamination")
    require(total_number >= 0 and math.isfinite(total_number), "invalid physical impulsive rate")
    require(birth_side in {"no-front-yet","upstream"} and (birth_side != "upstream" or front_generation > 0), "unsupported downstream/stale front birth")
    dimensions=(asset["footprint"],asset["time_response"],asset["momentum_spectrum"],asset["pitch_distribution"],asset["species_mixture"])
    for records in dimensions:
        finite(r["probability"] for r in records)
        require(records and abs(math.fsum(r["probability"] for r in records)-1) < 1e-12 and all(r["probability"] >= 0 for r in records), "unnormalized source measure")
    require(all(r["open_mapping"] is True and r["supported"] is True for r in asset["footprint"]), "closed/unsupported impulsive mapping")
    require(all(-1 <= r["mu"] <= 1 for r in asset["pitch_distribution"]), "invalid source pitch support")
    require(all(r["momentum_si"] >= 0 for r in asset["momentum_spectrum"]), "invalid momentum spectrum")
    finite(r["momentum_si"] for r in asset["momentum_spectrum"])
    finite(r["mu"] for r in asset["pitch_distribution"])
    finite(r["time_s"] for r in asset["time_response"])
    finite(r["mass_kg"] for r in asset["species_mixture"])
    births=[]
    import itertools
    for index,records in enumerate(itertools.product(*dimensions)):
        footprint,time,momentum,pitch,species=records
        require(species["mass_kg"] > 0 and species["species_id"], "unknown impulsive species")
        number=total_number
        for record in records: number*=record["probability"]
        c=299792458.0; pc=momentum["momentum_si"]*c; rest=species["mass_kg"]*c*c
        kinetic=pc*pc/(math.hypot(rest,pc)+rest)
        births.append(dict(immutable_birth_id=index+1,origin="ImpulsiveCoronalRelease",number=number,kinetic_energy_j=kinetic,
            species_id=species["species_id"],footprint_id=footprint["id"],time_s=time["time_s"],momentum_si=momentum["momentum_si"],mu=pitch["mu"],
            birth_side=birth_side,front_generation=front_generation))
    require(abs(math.fsum(b["number"] for b in births)-total_number) <= 1e-12*max(1,total_number), "physical source weight closure")
    return freeze_record(dict(schema="sccm-impulsive-birth-plan-v6",asset_identity=asset["identity"],births=births,
        physical_number=total_number,shock_first_passage_number=0,shock_calibration_mutated=False))
