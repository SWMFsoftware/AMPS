#!/usr/bin/env python3
"""Finite-time response of test particles to prescribed circular Alfvén fields.

This is an independent Lorentz/orbit reference, not a production mover or a
diffusion model. Every trajectory contains ONE real positive-k circular wave:
B/B0=(a*cos(psi),-h*a*sin(psi),1), E/(c*B0)=-beta_u*ez cross deltaB/B0.
psi=k*(z-u*t)+phi, h=+/-1. A k mixture is an ensemble of separate fields,
not a coherent multimode field. Uniform-phase moments are deterministic
quadrature; their variance/(2T) is deliberately NOT called a diffusivity.

Particles obey the relativistic Lorentz equation. The analytic reference
integrates the first-order force along the unperturbed mean-field gyro-orbit,
independently of any SEP QLT coefficient. Exact E-u*pz is a nonlinear invariant
of the prescribed one-speed field. Collective stability, damping, kinetic
dispersion and the Markov/random-phase closure are outside this experiment.
"""
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import math
from typing import Any, Dict, List, Tuple

import mpmath as mp
import numpy as np


class Unresolved(RuntimeError):
    """A requested numerical orbit/reference cannot support a resolved claim."""


def explicit_number(value: Any, name: str) -> float:
    if not isinstance(value, str):
        raise ValueError(name+" must be an explicit decimal string")
    try:exact=Decimal(value)
    except InvalidOperation as error:raise ValueError(name+" must be an explicit finite decimal") from error
    result=float(value)
    if not math.isfinite(result):
        raise ValueError(name+" must be finite")
    if not exact.is_zero() and result==0:
        raise Unresolved(name+" is nonzero in the raw decimal input but underflows binary64; no zero substituted")
    return result


def exact_mp(value: float) -> mp.mpf:
    """Convert the ACTUAL numerical coordinate, without decimal round-trip loss."""
    numerator,denominator=float(value).as_integer_ratio()
    return mp.mpf(numerator)/denominator


@dataclass(frozen=True)
class Probe:
    mass_kg: float
    charge_c: float
    light_speed_mps: float
    field_t: float
    wave_speed_mps: float
    helicity: int
    momentum: Tuple[float,float,float]  # Actual p/(mc), mean-B coordinates.

    def __post_init__(self) -> None:
        values=(self.mass_kg,self.charge_c,self.light_speed_mps,self.field_t,self.wave_speed_mps,*self.momentum)
        if not all(math.isfinite(x) for x in values) or min(self.mass_kg,self.light_speed_mps,self.field_t)<=0 or self.charge_c==0:
            raise ValueError("probe requires explicit finite positive mass/c/B0 and nonzero signed charge")
        if abs(self.wave_speed_mps)>=self.light_speed_mps or self.helicity not in (-1,1) or self.momentum[0]<=0 or self.momentum[1]!=0:
            raise ValueError("subluminal signed wave speed, circular helicity and initial gyroanglezero/interior pitch required")
        if not math.isfinite(self.momentum_norm) or self.momentum_norm==0:
            raise Unresolved("probe momentum norm is outside resolved representation")

    @property
    def sign(self) -> int:
        return 1 if self.charge_c>0 else -1

    @property
    def gyro_rest_per_s(self) -> float:
        return abs(self.charge_c)*self.field_t/self.mass_kg

    @property
    def gamma(self) -> float:
        return math.sqrt(1+sum(x*x for x in self.momentum))

    @property
    def momentum_norm(self) -> float:
        return math.sqrt(sum(x*x for x in self.momentum))

    @property
    def beta_u(self) -> float:
        return self.wave_speed_mps/self.light_speed_mps

    @property
    def mean_pitch(self) -> float:
        return self.momentum[2]/self.momentum_norm

    @property
    def kinetic_dimensionless(self) -> float:
        # Stable gamma-1, including nonrelativistic probes near rest.
        return self.momentum_norm**2/(self.gamma+1)


def resolve_probe(numbers: Dict[str,Any],scales: Dict[str,Any],control: Dict[str,Any]) -> Probe:
    """Use the existing equilibrium's numbers/scales; no copied physical constants.

    Thermal n/B/T/model data belong to the shared equilibrium helper. They
    determine the explicitly selected wave speed, but this test-particle
    calculation neither samples that thermal bath nor evolves its pressure.
    """
    source=control["alfven_speed_source"]
    if source not in ("alfven_ion_mps","alfven_total_mass_mps"):
        raise ValueError("explicit ion-inertia or total-mass Alfven speed source required")
    species=control["probe"]["species"]
    if species not in ("proton","electron"):
        raise ValueError("controlled probe species must be proton or electron")
    def physical_float(value: Any,name: str) -> float:
        result=float(value)
        if not math.isfinite(result) or (value!=0 and result==0):raise Unresolved(name+" is not representable in the orbit's binary64 coordinates")
        return result
    mass=physical_float(numbers[species+"_mass_kg"],"particle mass")
    magnitude=physical_float(numbers["charge_magnitude_c"],"charge")
    charge=magnitude if species=="proton" else -magnitude
    c=physical_float(numbers["light_speed_mps"],"light speed")
    kinetic=explicit_number(control["probe"]["kinetic_energy_ev"],"probe energy")
    mu=explicit_number(control["probe"]["mean_pitch"],"mean pitch")
    branch=control["probe"]["branch"]
    helicity=control["probe"]["helicity"]
    if kinetic<=0 or not -1<mu<1 or branch not in (-1,1) or helicity not in (-1,1):
        raise ValueError("positive probe energy, interior mean pitch, signed branch/helicity required")
    kappa=kinetic*magnitude/(mass*c*c)
    if not kappa>0 or not math.isfinite(kappa):
        raise Unresolved("positive requested particle kinetic energy cannot be resolved in momentum units")
    p=math.sqrt(kappa*(kappa+2))
    probe=Probe(mass,charge,c,physical_float(numbers["magnetic_field_t"],"mean B"),branch*physical_float(scales[source],"selected wave speed"),helicity,
                (p*math.sqrt((1-mu)*(1+mu)),0.0,p*mu))
    if not all(math.isfinite(x) for x in probe.momentum) or not 0<probe.momentum_norm or not abs(probe.beta_u)<1:
        raise Unresolved("probe momentum or prescribed wave frame is unresolved/nonphysical")
    return probe


def resonance_scale_per_m(probe: Probe) -> float:
    relative=probe.momentum[2]/probe.gamma-probe.beta_u
    if relative==0:
        raise Unresolved("finite k cannot resonate at v_parallel=u; no resonance-gap floor supplied")
    result=probe.gyro_rest_per_s/(probe.gamma*probe.light_speed_mps*abs(relative))
    if not math.isfinite(result) or result<=0:raise Unresolved("resonant coordinate is outside resolved finite k representation")
    return result


def observation_times(probe: Probe,cycles: List[str]) -> np.ndarray:
    values=[explicit_number(value,"observation gyroperiods") for value in cycles]
    if not values or values[0]<=0 or any(a>=b for a,b in zip(values,values[1:])):
        raise ValueError("positive strictly increasing observation gyroperiods required")
    return np.array([value*2*math.pi*probe.gamma for value in values])


def _rk_rhs(state: np.ndarray,kappa: np.ndarray,amplitude: float,probe: Probe) -> np.ndarray:
    # Independent full Lorentz forces in tau=|q|B0*t/m and pbar=p/(mc).
    # Integrating phase directly avoids subtracting large kz and kut values;
    # phase_dot=k*c/Omega0*(v_z/c-u/c) is the exact physical field phase.
    p=state[:,1:]
    gamma=np.sqrt(1+np.sum(p*p,axis=1))
    vx,vy,vz=(p/gamma[:,None]).T
    bx=amplitude*np.cos(state[:,0]);by=-probe.helicity*amplitude*np.sin(state[:,0])
    result=np.empty_like(state)
    result[:,0]=kappa*(vz-probe.beta_u)
    result[:,1]=probe.sign*(probe.beta_u*by+vy-vz*by)
    result[:,2]=probe.sign*(-probe.beta_u*bx+vz*bx-vx)
    result[:,3]=probe.sign*(vx*by-vy*bx)
    return result


def integrate_rk4(probe: Probe,k_per_m: np.ndarray,amplitude: float,phases: int,
                  steps_per_period: int,cycles: List[str],maximum_steps: int) -> Dict[str,Any]:
    if not len(k_per_m) or not np.all(np.isfinite(k_per_m)) or np.any(k_per_m<=0) or not math.isfinite(amplitude) or amplitude<0:
        raise ValueError("finite positive k support and explicit finite nonnegative wave amplitude required")
    if not isinstance(phases,int) or phases<4 or not isinstance(steps_per_period,int) or steps_per_period<=0:
        raise ValueError("explicit phase quadrature>=4 and positive integer gyroperiod step count required")
    times=observation_times(probe,cycles)
    counts=[]
    for value in cycles:
        count=explicit_number(value,"gyroperiods")*steps_per_period
        if count<=0 or count!=round(count):
            raise ValueError("observation times must align exactly with the explicitly selected gyroperiod step grid")
        counts.append(int(count))
    if max(counts)>maximum_steps:
        raise Unresolved("declared per-trajectory orbit resource bound exceeded")
    phases0=2*math.pi*np.arange(phases,dtype=float)/phases
    state=np.empty((len(k_per_m)*phases,4));state[:,0]=np.tile(phases0,len(k_per_m));state[:,1:]=probe.momentum
    khat=np.repeat(k_per_m*probe.light_speed_mps/probe.gyro_rest_per_s,phases)
    dt=2*math.pi*probe.gamma/steps_per_period
    saved=[];next_sample=0
    try:
        with np.errstate(over="raise",invalid="raise",divide="raise"):
            for step in range(max(counts)):
                first=_rk_rhs(state,khat,amplitude,probe)
                second=_rk_rhs(state+dt*first/2,khat,amplitude,probe)
                third=_rk_rhs(state+dt*second/2,khat,amplitude,probe)
                fourth=_rk_rhs(state+dt*third,khat,amplitude,probe)
                state=state+dt*(first+2*second+2*third+fourth)/6
                if not np.all(np.isfinite(state)):
                    raise Unresolved("nonfinite relativistic Lorentz trajectory")
                if step+1==counts[next_sample]:
                    saved.append(state.copy());next_sample+=1
                    if next_sample==len(counts):break
    except FloatingPointError as error:
        raise Unresolved("Lorentz arithmetic failed without a field/model fallback: "+str(error)) from error
    return {"states":saved,"times":times,"phases":phases0,"dt_tau":dt,
            "steps":counts,"trajectory_updates":max(counts)*len(state)}


def linear_pitch(probe: Probe,k_per_m: np.ndarray,amplitude: float,tau: float,phases: np.ndarray) -> np.ndarray:
    p=probe.momentum_norm;gamma=probe.gamma;beta=p/gamma;mu=probe.mean_pitch
    mismatch=k_per_m*probe.light_speed_mps/probe.gyro_rest_per_s*(probe.momentum[2]/gamma-probe.beta_u)-probe.helicity*probe.sign/gamma
    argument=mismatch*tau/2
    # np.sinc(x/pi)=sin(x)/x with its exact zero-argument limit. This is an
    # analytic finite-time force integral, not resonance broadening or a floor.
    pz_amplitude=-probe.helicity*probe.sign*amplitude*(probe.momentum[0]/gamma)*tau*np.sinc(argument/np.pi)
    coefficient=(1-probe.beta_u*mu/beta)/p
    result=coefficient*pz_amplitude[:,None]*np.sin(phases[None,:]+argument[:,None])
    if not np.all(np.isfinite(result)):raise Unresolved("linear finite-time response exceeds numerical representation")
    return result


def moments(probe: Probe,orbit: Dict[str,Any],k_per_m: np.ndarray,weights: np.ndarray,
            amplitude: float) -> List[Dict[str,float]]:
    rows=[];p0=np.array(probe.momentum);p0norm=probe.momentum_norm;gamma0=probe.gamma
    phase_count=len(orbit["phases"]);statistical=np.repeat(weights/phase_count,phase_count)
    for state,tau in zip(orbit["states"],orbit["times"]):
        p=state[:,1:];change=p-p0;norm=np.sqrt(np.sum(p*p,axis=1));gamma=np.sqrt(1+norm*norm)
        if not np.all(np.isfinite(p)) or np.any(norm<=0) or not np.all(np.isfinite(norm)):
            raise Unresolved("pitch is undefined or numerically unresolved at zero/nonfinite plasma momentum")
        # Difference identities avoid subtracting rest-mass dominated energies
        # or almost equal momentum norms. This does not change the force law.
        norm_change=np.sum(change*(p+p0),axis=1)/(norm+p0norm)
        dmu=change[:,2]/norm-p0[2]*norm_change/(norm*p0norm)
        dgamma=np.sum(change*(p+p0),axis=1)/(gamma+gamma0)
        dpz=change[:,2]
        linear=linear_pitch(probe,k_per_m,amplitude,float(tau),orbit["phases"]).reshape(-1)
        mean=float(np.sum(statistical*dmu));variance=float(np.sum(statistical*(dmu-mean)**2))
        invariant=dgamma-probe.beta_u*dpz
        # A response below the encoded momentum's angular resolution cannot
        # be declared physically identified just because two algorithms both
        # round it away. This is a diagnostic representation uncertainty,
        # never a nonzero wave/scattering floor or an edit to the trajectory.
        ulp=np.nextafter(np.abs(p),np.inf)-np.abs(p)
        encoding=np.sqrt(np.sum(ulp*ulp,axis=1))
        if np.any(norm<=encoding):raise Unresolved("pitch is unresolved by encoded momentum precision")
        pitch_resolution=2*encoding/(norm-encoding)
        rows.append({"mean_pitch_change":mean,"pitch_change_variance":variance,
            "linear_pitch_outside_physical_range_fraction":float(np.sum(statistical*((probe.mean_pitch+linear<-1)|(probe.mean_pitch+linear>1)))),
            "rms_pitch_linear_residual":float(np.sqrt(np.sum(statistical*(dmu-linear)**2))),
            "parallel_momentum_change_mean_kg_mps":float(np.sum(statistical*dpz))*probe.mass_kg*probe.light_speed_mps,
            "energy_change_mean_j":float(np.sum(statistical*dgamma))*probe.mass_kg*probe.light_speed_mps**2,
            "energy_change_variance_j2":float(np.sum(statistical*(dgamma-np.sum(statistical*dgamma))**2))*(probe.mass_kg*probe.light_speed_mps**2)**2,
            "momentum_encoding_pitch_resolution_estimate":float(np.max(pitch_resolution)),
            "maximum_sampled_invariant_residual_over_initial_kinetic":float(np.max(np.abs(invariant)))/probe.kinetic_dimensionless,
            "maximum_sampled_invariant_residual_j":float(np.max(np.abs(invariant)))*probe.mass_kg*probe.light_speed_mps**2})
    return rows


def analytic_variance(probe: Probe,k_per_m: np.ndarray,weights: np.ndarray,amplitude: float,
                      tau: float,digits: int,mixture: Tuple[float,float]=None) -> mp.mpf:
    """Independent high-precision phase average; optional continuous UNIFORM k.

    The integral is of a deterministic finite-time sinc-squared response. The
    optional mixture integrates separate single-wave experiments; it neither
    constructs nor qualifies a coherent/random-phase broadband wave field.
    Actual binary64 coordinates are used, with their mapping recorded outside.
    """
    with mp.workdps(digits):
        px,py,pz=[exact_mp(x) for x in probe.momentum]
        norm=mp.sqrt(px*px+py*py+pz*pz);gamma=mp.sqrt(1+norm*norm)
        bu=exact_mp(probe.beta_u);time=exact_mp(tau)
        relative=pz/gamma-bu;factor=(1-bu*(pz/norm)/(norm/gamma))/norm
        coefficient=-probe.helicity*probe.sign*exact_mp(amplitude)*(px/gamma)*time*factor
        def value(k: mp.mpf) -> mp.mpf:
            mismatch=k*exact_mp(probe.light_speed_mps)/exact_mp(probe.gyro_rest_per_s)*relative-probe.helicity*probe.sign/gamma
            argument=mismatch*time/2
            sinc=mp.mpf(1) if argument==0 else mp.sin(argument)/argument
            return coefficient*coefficient*sinc*sinc/2
        if mixture is None:
            return +mp.fsum(exact_mp(float(w))*value(exact_mp(float(k))) for k,w in zip(k_per_m,weights))
        low,high=[exact_mp(x) for x in mixture]
        resonant=(probe.helicity*probe.sign/gamma)/(exact_mp(probe.light_speed_mps)/exact_mp(probe.gyro_rest_per_s)*relative) if relative!=0 else None
        cuts=[low,high]
        if resonant is not None and low<resonant<high:cuts.insert(1,resonant)
        return +mp.quad(value,cuts)/(high-low)


def dop853_comparison(probe: Probe,k_per_m: np.ndarray,amplitude: float,orbit: Dict[str,Any],
                      rtol: float,atol: float,phase_samples: int) -> Dict[str,float]:
    """A second ODE algorithm and separately written vector Lorentz force.

    Agreement is numerical crosschecking, not an independent physical model.
    The analytic force integral, gyro solution and invariant supply governing
    references. Only the explicitly selected k/phase subset is claimed here.
    """
    from scipy.integrate import solve_ivp
    maximum_momentum=0.0;maximum_phase=0.0;maximum_pitch=0.0;selected=0
    phases=len(orbit["phases"])
    indices=sorted(set(int(i*phases/phase_samples) for i in range(min(phase_samples,phases))))
    for mode in sorted(set((0,len(k_per_m)-1))):
        kappa=float(k_per_m[mode])*probe.light_speed_mps/probe.gyro_rest_per_s
        for phase in indices:
            def rhs(_time: float,y: np.ndarray) -> np.ndarray:
                momentum=y[1:];velocity=momentum/math.sqrt(1+float(np.dot(momentum,momentum)))
                magnetic=np.array([amplitude*math.cos(y[0]),-probe.helicity*amplitude*math.sin(y[0]),1.0])
                electric=-probe.beta_u*np.cross(np.array([0.0,0.0,1.0]),magnetic)
                return np.r_[kappa*(velocity[2]-probe.beta_u),probe.sign*(electric+np.cross(velocity,magnetic))]
            initial=np.array([orbit["phases"][phase],*probe.momentum])
            solution=solve_ivp(rhs,(0,float(orbit["times"][-1])),initial,method="DOP853",
                t_eval=orbit["times"],rtol=rtol,atol=atol)
            if not solution.success or not np.all(np.isfinite(solution.y)):
                raise Unresolved("independent DOP853 trajectory did not resolve: "+solution.message)
            for sample,state in enumerate(orbit["states"]):
                rk=state[mode*phases+phase];other=solution.y[:,sample]
                maximum_momentum=max(maximum_momentum,float(np.linalg.norm(rk[1:]-other[1:]))/probe.momentum_norm)
                maximum_phase=max(maximum_phase,abs(float(rk[0]-other[0])))
                maximum_pitch=max(maximum_pitch,abs(rk[3]/np.linalg.norm(rk[1:])-other[3]/np.linalg.norm(other[1:])))
            selected+=1
    return {"selected_trajectories":selected,"maximum_relative_momentum_difference":maximum_momentum,
            "maximum_phase_difference_rad":maximum_phase,"maximum_absolute_pitch_difference":float(maximum_pitch)}
