#include "../pic.h"

void PIC::GYROKINETIC::GetDriftVelocity(double *vDrift,double vParallel,double vNorm,int spec) {
  //==================================================================================
  // Implementation notes (see header for full documentation)
  //----------------------------------------------------------------------------------
  // This routine reproduces the historical drift calculation that lived in
  // PIC::CPLR::GetDriftVelocity(), but with a reduced-state API and with the
  // "gamma bug" removed (non-relativistic gamma=1).
  //
  // Drifts included:
  //   (1) E×B drift       : (E×B)/B^2   (computed as stencil-averaged per-cell drift)
  //   (2) ∇B drift        : (μ/q) (B×∇_⊥B) / B^2
  //   (3) curvature drift : (p_∥^2/(q m)) B×((B·∇)B) / B^4
  //
  // Reduced particle state mapping (non-rel):
  //   p_∥ = m vParallel
  //   μ   = (1/2) m vNorm^2 / B
  //
  // The function assumes the caller already called PIC::CPLR::InitInterpolationStencil()
  // for the position of interest. It does NOT take coordinates as an argument.
  //==================================================================================

  // ---------------------------
  // Initialize output to zero
  // ---------------------------
  vDrift[0]=0.0; vDrift[1]=0.0; vDrift[2]=0.0;

  // ---------------------------
  // Species mass and charge
  // ---------------------------
  // Tables are stored in SI units in AMPS. Convert when the solver input is normalized,
  // matching the Boris mover behavior for ECSIM.
  double ParticleMass   = PIC::MolecularData::MolMass[spec];
  double ParticleCharge = PIC::MolecularData::ElectricChargeTable[spec];

  if (_PIC_FIELD_SOLVER_INPUT_UNIT_==_PIC_FIELD_SOLVER_INPUT_UNIT_NORM_) {
    ParticleMass   = picunits::si2no_m(ParticleMass,   PIC::Units::Factors);
    ParticleCharge = picunits::si2no_q(ParticleCharge, PIC::Units::Factors);
  }

  // Neutrals have no charge => grad-B and curvature drifts are undefined; E×B is also not
  // a "particle drift" in that case. For safety, return zero for all neutral species.
  if (ParticleCharge==0.0) return;

  // ---------------------------
  // Background fields at stencil
  // ---------------------------
  // These return the *interpolated* fields at the stencil point (already set by InitInterpolationStencil()).
  double E[3],B[3];
  PIC::CPLR::GetBackgroundMagneticField(B);
  PIC::CPLR::GetBackgroundElectricField(E);

  const double absB2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
  if (absB2<=0.0) {
    // Guiding-center drift is not defined for B=0; return zero safely.
    return;
  }

  const double absB  = sqrt(absB2);
  const double absB4 = absB2*absB2;

  // ---------------------------
  // (1) E×B drift
  // ---------------------------
  // Historical behavior: compute drift per stencil cell and average the drift velocities directly:
  //     v_EB ≈ Σ w_i (E_i×B_i)/|B_i|^2
  // This is not identical to (⟨E⟩×⟨B⟩)/|⟨B⟩|^2 when B varies strongly across the stencil; however
  // it matches legacy AMPS behavior and is stable in practice.
  // IMPORTANT:
  // Do NOT copy the stencil. Use a reference/pointer to the actual thread-local stencil
  // object that was filled by PIC::CPLR::InitInterpolationStencil(...).
  // This avoids (a) unnecessary copying and (b) accidental divergence from the true
  // stencil used by other interpolation calls in the same thread.
  PIC::InterpolationRoutines::CellCentered::cStencil *StencilPtr;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  StencilPtr = PIC::InterpolationRoutines::CellCentered::StencilTable + omp_get_thread_num();
#else
  StencilPtr = PIC::InterpolationRoutines::CellCentered::StencilTable;
#endif

  PIC::InterpolationRoutines::CellCentered::cStencil &Stencil = *StencilPtr;

  for (int iStencil=0;iStencil<Stencil.Length;iStencil++) {
    double cellB[3],cellE[3];

    // cell center coordinates are stored by the stencil
    PIC::Mesh::cDataCenterNode* Cell = Stencil.cell[iStencil]; 
    double *x = Cell->GetX();
    const double w  = Stencil.Weight[iStencil];

    // Pull cell-centered fields depending on coupler mode
    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      PIC::CPLR::SWMF::GetBackgroundElectricField(cellE,Cell);
      PIC::CPLR::SWMF::GetBackgroundMagneticField(cellB,Cell);
      break;

    case _PIC_COUPLER_MODE__DATAFILE_:
      // This helper does not accept Time as an argument. Passing NAN is treated by the
      // DATAFILE coupler as "use the currently selected time slice / internal time state".
      PIC::CPLR::DATAFILE::GetBackgroundElectricField(cellE,Cell,NAN);
      PIC::CPLR::DATAFILE::GetBackgroundMagneticField(cellB,Cell,NAN);
      break;

    case _PIC_COUPLER_MODE__T96_:
    case _PIC_COUPLER_MODE__T05_:
    case _PIC_COUPLER_MODE__KMAG_:
      for (int i=0;i<3;i++) cellE[i]=0.0;
      PIC::CPLR::DATAFILE::GetBackgroundData(cellB,3,PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset,Cell);
      break;

    case _PIC_COUPLER_MODE__OFF_:
      // For ECSIM (fields from solver), E and B are not provided via CPLR stencil cells.
      // This drift helper is meant for CPLR-based fields. If you need ECSIM-based drifts,
      // implement an ECSIM-aware version that samples solver fields at x.
      exit(__LINE__,__FILE__,"PIC::GYROKINETIC::GetDriftVelocity(): _PIC_COUPLER_MODE__OFF_ not supported by this CPLR/stencil implementation");
      break;

    default:
      exit(__LINE__,__FILE__,"PIC::GYROKINETIC::GetDriftVelocity(): coupler mode not implemented");
    }

    const double cellAbsB2 = cellB[0]*cellB[0] + cellB[1]*cellB[1] + cellB[2]*cellB[2];
    if (cellAbsB2<=0.0) continue;

    double ExB[3];
    Vector3D::CrossProduct(ExB,cellE,cellB);

    // (E×B)/B^2 times weight
    const double c = w / cellAbsB2;
    vDrift[0] += c*ExB[0];
    vDrift[1] += c*ExB[1];
    vDrift[2] += c*ExB[2];
  }

  // ---------------------------
  // Particle moments needed for drifts
  // ---------------------------
  // Non-relativistic:
  //   p_parallel = m vParallel
  //   mu = (1/2) m v_perp^2 / B
  const double gamma = 1.0;
  const double pParallel = gamma*ParticleMass*vParallel;
  const double Mu = (gamma*gamma)*0.5*ParticleMass*vNorm*vNorm/absB;

  // ---------------------------
  // (2) Grad-B drift
  // ---------------------------
  // Compute ∇|B| and keep only perpendicular component.
  double gradAbsB_perp[3];
  // No explicit Time argument: pass NAN (see DATAFILE note above).
  PIC::CPLR::GetAbsBackgroundMagneticFieldGradient(gradAbsB_perp,NAN);
  Vector3D::Orthogonalize(B,gradAbsB_perp);

  // v_gradB = (Mu/(q*gamma*B^2)) (B×∇_⊥|B|)
  double BxGradB[3];
  Vector3D::CrossProduct(BxGradB,B,gradAbsB_perp);

  {
    const double c = Mu / (ParticleCharge*gamma*absB2);
    vDrift[0] += c*BxGradB[0];
    vDrift[1] += c*BxGradB[1];
    vDrift[2] += c*BxGradB[2];
  }

  // ---------------------------
  // (3) Curvature drift
  // ---------------------------
  // Uses (B·∇)B from the background magnetic-field gradient tensor.
  double gradB[9];
  // No explicit Time argument: pass NAN (see DATAFILE note above).
  PIC::CPLR::GetBackgroundMagneticFieldGradient(gradB,NAN);

  // t1 = (B·∇)B / B^4  (component-wise)
  double t1[3];
  t1[0] = (B[0]*gradB[0] + B[1]*gradB[1] + B[2]*gradB[2]) / absB4;
  t1[1] = (B[0]*gradB[3] + B[1]*gradB[4] + B[2]*gradB[5]) / absB4;
  t1[2] = (B[0]*gradB[6] + B[1]*gradB[7] + B[2]*gradB[8]) / absB4;

  // B×t1 = B×((B·∇)B)/B^4
  double BxT1[3];
  Vector3D::CrossProduct(BxT1,B,t1);

  {
    const double c = (pParallel*pParallel) / (ParticleCharge*gamma*ParticleMass);
    vDrift[0] += c*BxT1[0];
    vDrift[1] += c*BxT1[1];
    vDrift[2] += c*BxT1[2];
  }

  // ---------------------------
  // Safety: clamp non-finite outputs
  // ---------------------------
  // Can happen if B is extremely small at one stencil cell, etc.
  if ((isfinite(vDrift[0])==false) || (isfinite(vDrift[1])==false) || (isfinite(vDrift[2])==false)) {
    vDrift[0]=0.0; vDrift[1]=0.0; vDrift[2]=0.0;
  }
}

