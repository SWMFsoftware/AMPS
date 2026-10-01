"""Frozen Stage-14 joint-source, renewal and drift-report products.

This module is host independent. A tabulated conditional kernel is an explicit
input authority, not an inferred replacement for unresolved shock physics.
Every physical branch is followed with its existing momentum and ancestry.
No routine here mutates the first-passage source or samples its original g(p).
"""
from __future__ import annotations
import copy
import math
from .core import require, finite, freeze_record, verify_frozen, digest, load_json, canonical
from .research import content_hash


def kinetic(mass, momentum):
    c=299792458.0; pc=math.sqrt(math.fsum(x*x for x in momentum))*c; rest=mass*c*c
    return pc*(pc/(math.hypot(rest,pc)+rest))


def validate_family(bundle):
    verify_frozen(bundle)
    require(bundle["schema"] == "sep-field-line-family-bundle-v4", "unsupported family bundle major")
    groups={}
    for m in bundle["members"]:
        finite([m[key] for key in ("momentum_lower_si","momentum_upper_si","time_lower_s","time_upper_s","offset_m","joint_measure","integrated_peclet")])
        require(m["momentum_lower_si"] >= 0 and m["momentum_upper_si"] > m["momentum_lower_si"] and
                m["time_upper_s"] > m["time_lower_s"] and m["offset_m"] > 0 and m["joint_measure"] > 0,
                "invalid joint-source support/measure")
        require(m["source_is_separable"] is False and m["geometry_generation"] > 0,
                "separable/stale joint-source member")
        require(set(m["geometry_certificate"]) == {"positive_jacobian","unique_normal_root","reach","overlap","mask","solar_clearance"} and
                all(v is True for v in m["geometry_certificate"].values()), "uncertified family member geometry")
        content_hash(m["joint_measure_derivation_sha256"])
        key=(m["species_id"],m["patch_id"],m["time_lower_s"],m["time_upper_s"])
        groups.setdefault(key,[]).append(m)
    require(groups,"empty family")
    for rows in groups.values():
        ordered=sorted(rows,key=lambda r:r["momentum_lower_si"])
        for a,b in zip(ordered,ordered[1:]):
            require(a["momentum_upper_si"] == b["momentum_lower_si"] and a["geometry_generation"] == b["geometry_generation"],
                    "joint-source momentum gap/overlap/generation mismatch")
    return bundle


def publish_family(bundle, path):
    """Atomic single-product publication; an existing path is never overwritten."""
    from pathlib import Path
    import os
    path=Path(path); validate_family(bundle)
    require(not path.exists(), "family output already exists")
    temp=path.with_name(path.name+".partial")
    with temp.open('xb') as stream: stream.write(canonical(bundle))
    try:
        # link+unlink gives no-overwrite atomic publication on POSIX filesystems.
        os.link(str(temp),str(path))
    finally: temp.unlink()


def read_family(path, supported_major=4):
    require(supported_major == 4,"legacy readers must reject the new joint-measure bundle major")
    return validate_family(load_json(path))


CONDITION_KEYS={"momentum_si","pitch_cosine","surface_position_m","time_s","species_id",
                "frame_id","wave_authority","front_generation","reference_family_identity"}


def renewal_checkpoint(first_passage, kernel, cohort, transactions=()):
    """Restart state binds immutable first passage, kernel and physical cohort."""
    verify_frozen(first_passage); verify_frozen(kernel)
    require(kernel["schema"] == "sccm-conditional-renewal-kernel-v6" and kernel["validated_work_authority"] is True,
            "missing conditional kernel/work qualification")
    content_hash(kernel["work_validation_sha256"])
    require(set(kernel["conditioning"]) == CONDITION_KEYS,"incomplete renewal conditioning")
    finite([kernel["relative_tolerance"],kernel["conditioning"]["pitch_cosine"],kernel["conditioning"]["time_s"]]+
        kernel["conditioning"]["momentum_si"]+kernel["conditioning"]["surface_position_m"])
    require(kernel["relative_tolerance"]>0 and abs(kernel["conditioning"]["pitch_cosine"])<=1,"invalid conditional kernel tolerance/pitch")
    content_hash(kernel["conditioning"]["wave_authority"]);content_hash(kernel["conditioning"]["reference_family_identity"])
    require(kernel["conditioning"]["front_generation"] > 0,"stale renewal front")
    finite([cohort["number"],cohort["mass_kg"],cohort["cycle"],cohort["time_s"]]+cohort["momentum_si"])
    require(cohort["number"] > 0 and cohort["mass_kg"] > 0 and cohort["cycle"] >= 0 and cohort["ancestry_id"],"invalid cohort ancestry/weight")
    return freeze_record(dict(schema="sccm-renewal-checkpoint-v6",first_passage=first_passage,kernel=kernel,
        cohort=cohort,transactions=list(transactions),first_passage_identity=first_passage["identity"]))


def renew(checkpoint):
    """Deterministic weighted transition of one conditional tabulated cohort.

    Splitting the represented cohort changes no probability or energy. A host
    stochastic realization may use the same physical branch measure, but must
    independently qualify its counter-based stream and MPI/restart ownership.
    Absorbed/downstream reservoirs retain four-momentum. Only ReReleased creates
    a continuation, with positive residence time and a *conditioned* new p'.
    """
    verify_frozen(checkpoint)
    require(checkpoint["schema"] == "sccm-renewal-checkpoint-v6","wrong renewal checkpoint")
    kernel=checkpoint["kernel"];verify_frozen(kernel); first=checkpoint["first_passage"];verify_frozen(first)
    # Revalidate the authority on restart, even if a caller has deliberately
    # re-frozen a changed JSON object with a valid new content hash.
    renewal_checkpoint(first,kernel,checkpoint["cohort"],checkpoint["transactions"])
    require(first["identity"] == checkpoint["first_passage_identity"],"first-passage mutation")
    c=checkpoint["cohort"]; condition=kernel["conditioning"]
    require(c["momentum_si"] == condition["momentum_si"] and c["time_s"] == condition["time_s"] and
            c["species_id"] == condition["species_id"] and c["frame_id"] == condition["frame_id"],
            "conditional kernel cannot redraw original source momentum/species/time/frame")
    branches=kernel["branches"]; require(branches,"empty conditional kernel")
    probabilities=[b["probability"] for b in branches];finite(probabilities)
    require(all(x>=0 for x in probabilities) and abs(math.fsum(probabilities)-1)<=1e-12,"conditional probability does not normalize")
    incoming_k=kinetic(c["mass_kg"],c["momentum_si"]); outgoing_k=0.; work=0.; outgoing_p=[0.,0.,0.]; impulse=[0.,0.,0.]
    sinks={"Absorbed":0.,"Downstream":0.,"ReReleased":0.}; continuations=[]
    for index,b in enumerate(branches):
        require(b["outcome"] in sinks,"unknown conditional outcome")
        p=b["momentum_si"];finite(p+[b["residence_time_s"],b["shock_work_j"]]+b["shock_impulse_si"])
        require(len(p)==len(b["shock_impulse_si"])==3 and b["residence_time_s"] >= 0,"invalid four-momentum branch")
        k=kinetic(c["mass_kg"],p); dw=k-incoming_k; dp=[a-z for a,z in zip(p,c["momentum_si"])]
        scale=max(incoming_k,k,abs(dw),1e-300)
        require(abs(dw-b["shock_work_j"])<=kernel["relative_tolerance"]*scale and
                all(abs(a-z)<=kernel["relative_tolerance"]*max(abs(a),abs(z),abs(v),1e-300)
                    for a,z,v in zip(dp,b["shock_impulse_si"],c["momentum_si"])), "conditional shock-work/four-momentum mismatch")
        number=c["number"]*b["probability"];sinks[b["outcome"]]+=number
        outgoing_k+=number*k;work+=number*dw
        outgoing_p=[a+number*z for a,z in zip(outgoing_p,p)];impulse=[a+number*z for a,z in zip(impulse,dp)]
        if b["outcome"] == "ReReleased":
            require(b["residence_time_s"]>0 and p != c["momentum_si"],"zero-delay/unchanged renewal double counts first passage")
            continuations.append(dict(c,number=number,momentum_si=p,cycle=c["cycle"]+1,time_s=c["time_s"]+b["residence_time_s"],
                parent_transaction=len(checkpoint["transactions"]),branch_id=index))
    ledger=dict(cycle=c["cycle"],ancestry_id=c["ancestry_id"],incoming_number=c["number"],outcome_number=sinks,
        incoming_kinetic_energy_j=c["number"]*incoming_k,outgoing_kinetic_energy_j=outgoing_k,shock_work_j=work,
        incoming_momentum_si=[c["number"]*x for x in c["momentum_si"]],outgoing_momentum_si=outgoing_p,shock_impulse_si=impulse)
    require(abs(outgoing_k-c["number"]*incoming_k-work)<=kernel["relative_tolerance"]*max(outgoing_k,c["number"]*incoming_k,1e-300),"cohort energy closure")
    return freeze_record(dict(schema="sccm-renewal-result-v6",checkpoint_identity=checkpoint["identity"],
        first_passage_identity=first["identity"],ledger=ledger,continuations=continuations,
        transactions=checkpoint["transactions"]+[ledger],first_passage_mutated=False))


def drift_exposure(segments, denominators, frame, generation):
    """Report net displacement and path exposure with explicit width conventions.

    Input velocities already come from the same equation-frame coherent
    provider. Widths are measured by the drift-disabled perpendicular-diffusion
    experiment, not inferred from these drift trajectories.
    """
    require(frame and generation>0 and segments,"missing drift report identity")
    net=[0.,0.,0.];exposure=0.;duration=0.
    for segment in segments:
        finite(segment["velocity_m_per_s"]+[segment["dt_s"]])
        require(len(segment["velocity_m_per_s"])==3 and segment["dt_s"]>0 and
                segment["frame"]==frame and segment["generation"]==generation,"drift frame/generation/cadence mismatch")
        delta=[segment["dt_s"]*x for x in segment["velocity_m_per_s"]]
        net=[x+y for x,y in zip(net,delta)];exposure+=math.sqrt(math.fsum(x*x for x in delta));duration+=segment["dt_s"]
    needed={"tube_radius_m","rms_1d_perpendicular_width_m","rms_plane_perpendicular_width_m","diffusion_registration_sha256"}
    require(set(denominators)==needed,"missing/unlabeled drift-disabled diffusion denominator")
    content_hash(denominators["diffusion_registration_sha256"])
    widths={k:v for k,v in denominators.items() if k.endswith('_m')};finite(widths.values());require(all(v>0 for v in widths.values()),"invalid diffusion width")
    magnitude=math.sqrt(math.fsum(x*x for x in net))
    return freeze_record(dict(schema="sccm-drift-exposure-v6",frame=frame,generation=generation,duration_s=duration,
        net_displacement_m=net,net_magnitude_m=magnitude,path_exposure_m=exposure,denominators=denominators,
        net_ratios={k:magnitude/v for k,v in widths.items()},exposure_ratios={k:exposure/v for k,v in widths.items()},
        width_convention="sigma_1d=sqrt(2*kappa_perp*T); sigma_plane=sqrt(4*kappa_perp*T)"))
