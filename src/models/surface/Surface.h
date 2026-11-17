/*
 * Surface.h
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 *
 * Gas-surface interaction (GSI) models.
 *
 * Every model is a function of the incident velocity v, the face frame (n = TriangleCutFace->ExternalNormal, the unit
 * normal that points out of the body into the gas, and the orthonormal tangents e0 = e0Orthogonal, e1 = e1Orthogonal),
 * the species and its parameters.  An incident particle has v.n < 0; the model replaces v with the re-emitted velocity
 * (v.n > 0) and returns _PARTICLE_REJECTED_ON_THE_FACE_, or returns _PARTICLE_DELETED_ON_THE_FACE_ (adsorption).
 *
 * Model selection.  Each species has a model, Surface::Model::Selected[spec].  Its initial value is the compile-time
 * choice _SURFACE__GAS_SURFACE_INTERACTION_MODEL_ (Surface.dfn), so existing applications are unchanged.  It is changed
 * at run time with Surface::Model::Set(), or with the runtime input block 'physics: gas_surface_interaction'
 * (Surface::RuntimeInput::Apply, Surface_runtime_input.cpp; grammar src/INPUT_FILE.md).  The parameters are per species.
 * An unset wall temperature (0 K) or an unset or out-of-range parameter stops the run with a message naming the
 * parameter and the species (parameters are NaN until they are set).
 *
 * Notation: a = sqrt(2 k Tw/m) is the most probable speed at the wall temperature Tw; R, R1, ... are independent
 * uniform random numbers in (0,1].  Velocities in the CLL formulas are in units of a.
 *
 *  specular      u = v - 2 (v.n) n.
 *  diffuse       the flux-weighted half-range Maxwellian at Tw,
 *                  f(u) ~ (u.n) exp(-u^2/a^2), u.n > 0:
 *                  u.n = a sqrt(-ln R1),  u.e0 + i u.e1 = a sqrt(-ln R2) exp(2 pi i R3)
 *                (moments: <u.n> = sqrt(pi)/2 a, <(u.n)^2> = a^2, <(u.e0)^2> = <(u.e1)^2> = a^2/2, <u^2> = 2 a^2,
 *                i.e. the mean energy 2 k Tw of the molecules that a wall at Tw emits).
 *  maxwell       Maxwell (1879), Phil. Trans. R. Soc. 170, 231: diffuse with the probability alpha (the accommodation
 *                coefficient), else specular.
 *  adsorption    the particle is removed (sticking probability 1).
 *  cll           Cercignani & Lampis (1971), Transp. Theory Stat. Phys. 1, 101; in the form of Lord (1991), Phys.
 *                Fluids A 3, 706, and Lord (1991), Rarefied Gas Dynamics 17, VCH, 1427, as used in DSMC by Bird
 *                (1994), Molecular Gas Dynamics and the Direct Simulation of Gas Flows (Oxford), and by Padilla &
 *                Boyd (2009), J. Thermophys. Heat Transfer 23, 96 (MONACO).  The kernel factorizes into the normal part
 *                (w' = |v.n|, w = u.n, in units of a)
 *                  R_n(w'->w) = (2w/alpha_n) I0(2 sqrt(1-alpha_n) w w'/alpha_n) exp(-(w^2 + (1-alpha_n) w'^2)/alpha_n)
 *                and the tangential part (v_t, u_t: the tangential velocity vectors)
 *                  R_t(v_t->u_t) = 1/(pi alpha_t) exp(-|u_t - (1-sigma_t) v_t|^2/alpha_t),  alpha_t = sigma_t (2-sigma_t),
 *                with alpha_n the accommodation coefficient of the normal kinetic energy (CLL::
 *                NormalPartKineticEnergyAccommodationCoefficient, 0..1) and sigma_t the tangential momentum
 *                accommodation coefficient (CLL::TangentialMomentumAccommodationCoefficient, 0..2; sigma_t > 1 is
 *                back-scattering; alpha_t is the tangential energy accommodation coefficient).  Lord's sampling:
 *                  w   = sqrt(r^2 + w_m^2 + 2 r w_m cos th),  r = sqrt(-alpha_n ln R1), th = 2 pi R2, w_m = sqrt(1-alpha_n) w'
 *                  u_t = (1-sigma_t) v_t + r' (cos th' e0 + sin th' e1),  r' = sqrt(-alpha_t ln R3), th' = 2 pi R4
 *                (R_n is the Rice distribution of |w_m + r exp(i th)|; the tangential noise is an isotropic Gaussian, so
 *                both tangential components and their signs are kept, in any tangent frame).  Moments:
 *                  <w^2> = (1-alpha_n) w'^2 + alpha_n,  <u_t> = (1-sigma_t) v_t,  <|u_t|^2> = (1-sigma_t)^2 |v_t|^2 + alpha_t.
 *                alpha_n = sigma_t = 0 is specular reflection, alpha_n = sigma_t = 1 is diffuse reflection.
 *  quasi_specular_cll   the CLL kernel with its own (small) coefficients, QuasiSpecularReflection::CLL, sigma_t <= 1.
 *  nocilla       Nocilla (1963), "The surface re-emission law in free molecule flow", Rarefied Gas Dynamics 3, ed.
 *                J. A. Laurmann, Academic Press, vol. 1, 327; Hurlbut & Sherman (1968), Phys. Fluids 11, 486.  The
 *                re-emitted flux is a half-range drifting Maxwellian
 *                  f(u) ~ (u.n) exp(-|u - U_r d|^2/c_r^2),  u.n > 0,  c_r = sqrt(2 k T_r/m),
 *                  d = cos(theta_r) n + sin(theta_r) t,
 *                with t the unit tangential direction of the incident velocity (the drift is in the plane of incidence,
 *                forward; a random tangent direction for normal incidence).  The parameters U_r (drift speed),
 *                theta_r (drift angle from the normal, 0..90 deg) and T_r are given per species
 *                (QuasiSpecularReflection::Nocilla::DriftSpeed, DriftAngle [rad], Temperature), or from accommodation
 *                coefficients (below).  Sampling: u = c_r (x n + y1 t + y2 (n x t)) + U_r sin(theta_r) t, y1 and y2
 *                normal with the variance 1/2, and x > 0 with the density x exp(-(x-s)^2), s = U_r cos(theta_r)/c_r
 *                (exactly, by the decomposition in SampleFluxNormalSpeedRatio()).  U_r = 0, T_r = Tw is diffuse
 *                reflection.
 *                Accommodation coefficients (Mode[spec] = AccommodationCoefficients; the parameters of the Nocilla
 *                distribution are fixed by the momentum and energy transfer, as Hurlbut & Sherman (1968) apply the
 *                model; the equations are those of Cook & Hoffbauer (1997), NASA TP-3701, Sec. 2-4, Eqs. 1, 2, 10 and
 *                13-22).  With the tangential momentum accommodation coefficient sigma (Eq. 1), the normal one sigma'
 *                (Eq. 2) and the energy accommodation coefficient epsilon (Eq. 10), applied to each incident molecule
 *                (speed v, angle theta_i from the normal), the re-emitted molecules have the mean velocity and the mean
 *                square speed
 *                  <u>.t  = Q = (1-sigma) v sin(theta_i)                                         (Eq. 1)
 *                  <u>.n  = A = (1-sigma') v cos(theta_i) + sigma' sqrt(pi k Tw/(2m))              (Eq. 2)
 *                  <u^2>  = V2 = (1-epsilon) v^2 + epsilon 4 k Tw/m                              (Eq. 10)
 *                (sqrt(pi k Tw/2m) and 4 k Tw/m are the mean normal velocity and the mean square speed of the diffuse
 *                re-emission at Tw).  For the Nocilla distribution with S = U/c_r, Lambda = S cos(theta_u) (the drift
 *                U at the angle theta_u from the normal, Eqs. 18-20):
 *                  Q = U sin(theta_u)                                                           (Eq. 15)
 *                  A = c_r [Lambda e^{-Lambda^2} + sqrt(pi)(1/2 + Lambda^2) Sigma]/D,  Sigma = 1 + erf(Lambda) (Eq. 16)
 *                  V2 = c_r^2 [(2 + S^2) e^{-Lambda^2} + sqrt(pi)(5/2 + S^2) Lambda Sigma]/D,               (Eq. 17)
 *                  D = e^{-Lambda^2} + sqrt(pi) Lambda Sigma.
 *                Eliminating U: (V2 - Q^2)/A^2 = (<x^2>(Lambda) + 1)/<x>(Lambda)^2 =: F(Lambda), with <x>, <x^2> the
 *                moments of x exp(-(x-Lambda)^2) on x > 0.  F decreases monotonically from infinity (Lambda -> -inf)
 *                through 8/pi (Lambda = 0, diffuse) to 1 (Lambda -> +inf, a cold beam), so Lambda is found by bisection
 *                (instead of the iteration of Eqs. 21-22), then c_r = A/<x>(Lambda), U cos(theta_u) = Lambda c_r, and
 *                T_r = m c_r^2/(2k).  theta_u may exceed 90 deg (Lambda < 0: a drift toward the wall, allowed by Eqs.
 *                23-25).  The coefficients are 0..1 (epsilon = sigma = sigma' = 1 is diffuse reflection).  They are
 *                not consistent for every incident velocity: when V2 - Q^2 <= A^2 no distribution with a positive
 *                temperature has these moments (e.g. sigma = sigma' = 0 with epsilon > 0); the molecule is then
 *                re-emitted with the mean velocity A n + Q t (the cold limit, which keeps both momentum transfers) and
 *                counted in Nocilla::nInconsistentReemissions; F - 1 < 1.5e-8 (a nearly cold beam) is also re-emitted
 *                with the mean velocity (Nocilla::nColdLimitReemissions).
 *  schamberg     Schamberg (1959), "A new analytic representation of surface interaction for hyperthermal free
 *                molecule flow with application to neutral-particle drag estimates of satellites", RAND RM-2313;
 *                as summarized by Mostaza Prieto, Graziano & Roberts (2014), Prog. Aerosp. Sci. 64, 56, Sec. 3.2 and
 *                4.1.1, Eqs. (16), (17), (20)-(24) of the accepted manuscript.  Each molecule is re-emitted with one
 *                speed V_r into a conical beam of half-angle phi0 about an axis in the plane of incidence:
 *                  axis angle from the normal:  sin(theta_r) = (sin(theta_i))^nu   (Eq. 17, there with angles measured
 *                                               from the surface: cos(theta_r') = cos(theta_i')^nu)
 *                  speed: alpha = (E_i - E_r)/(E_i - E_w) with E_i = m V_i^2/2, E_r = m V_r^2/2 and E_w = 2 k Tw, the
 *                         mean energy of the molecules re-emitted diffusely by the wall (Eq. 15, 23), so
 *                         V_r^2 = (1-alpha) V_i^2 + alpha 4 k Tw/m   (Eq. 24, V_r = sqrt(1-alpha) V_i, when E_w << E_i)
 *                  beam: the flux per unit solid angle at the angle phi from the axis is ~ cos(pi phi/(2 phi0)),
 *                        phi <= phi0.  Its mean cosine is
 *                          Phi(phi0) = (1-x^2)/(1-4x^2) (sin(2 phi0)/2 - x)/(sin(phi0) - x),  x = 2 phi0/pi,
 *                        Schamberg's beam-shape function (Eq. 16), with Phi(0) = 1 and Phi(pi/2) = 2/3 (the cosine
 *                        law); the mean re-emitted velocity is V_r Phi(phi0) times the axis (Eqs. 20, 21).
 *                nu = 1, phi0 = 0, alpha = 0 is specular reflection; nu = inf, phi0 = 90 deg, alpha = 1 re-emits with
 *                the cosine law at the single speed sqrt(4 k Tw/m) (Schamberg's quasi-diffuse limit).  A direction of
 *                the cone below the surface (possible when theta_r + phi0 > 90 deg) is rejected and the direction is
 *                sampled again, i.e. the beam is truncated by the surface (an assumption of this implementation; the
 *                mean is then not Phi(phi0) times the axis).  Parameters per species:
 *                QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient (0..1), Nu (> 0, may be inf),
 *                BeamHalfAngle [rad] (0..pi/2).
 */

#ifndef _MODELS_SURFACE_SURFACE_H_
#define _MODELS_SURFACE_SURFACE_H_

#include <string>

#include "Surface.dfn"
#include "pic.h"

namespace AMPS { namespace Input { struct Block; } }

//the model of the gas-surface interaction
namespace Surface {

   //the models of the surface temeprature
   namespace Temperature {

     //isothermal: the wall temperature [K]; 0 (the initial value) means that it is not set
     namespace Isothremal {
       extern double Temp;
     }

     //Radiative Equilibrium
     namespace RadiativeEquilibrium {

     }
   }

   //the wall temperature of a face; stops the run when it is not set (<= 0 K)
   double GetSurfaceTemperature(CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

   //the model of each species (codes: SurfaceDefinitions.dfn)
   namespace Model {
     const int Diffuse=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__DIFFUSE_REFLECTION_;
     const int Specular=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__SPECULAR_REFLECTION_;
     const int CLL=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__CLL_;
     const int Maxwell=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__MAXWELL_;
     const int QuasiSpecularCLL=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__QUASI_SPECULAR_CLL_;
     const int Adsorption=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__ADSORPTION_;
     const int QuasiSpecularNocilla=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__QUASI_SPECULAR_NOCILLA_;
     const int QuasiSpecularSchamberg=_SURFACE__GAS_SURFACE_INTERACTION_MODEL__QUASI_SPECULAR_SCHAMBERG_;
     const int nTotalModels=8;

     extern int Selected[PIC::nTotalSpecies];

     //names used in the runtime input: diffuse, specular, cll, maxwell, quasi_specular_cll, adsorption, nocilla, schamberg
     const char* Name(int model);
     int Code(const std::string& name);   //-1: unknown name

     //select the model of a species (spec<0: all species); false (and the reason) for an unknown model code.
     //The parameters are checked when they are set and at every interaction (CheckParameters).
     bool Set(int spec,int model,std::string *error=NULL);

     //does the model of the species use the wall temperature
     bool UsesWallTemperature(int spec);
   }

   //check the parameters of a model (of the selected model) of a species: true, or false and the reason
   bool CheckParameters(int spec,int model,std::string *error=NULL);
   inline bool CheckParameters(int spec,std::string *error=NULL) {return CheckParameters(spec,Model::Selected[spec],error);}

   //diffuse refrelction
   namespace DiffuseReflection {
     //re-emit v from the face with the flux-weighted half-range Maxwellian at the temperature Temp
     void Reemit(double *v,double Temp,int spec,CutCell::cTriangleFace *TriangleCutFace);
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //specular refrelction
   namespace SpecularReflection {
     void Reflect(double *v,CutCell::cTriangleFace *TriangleCutFace);
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //adsorption
   namespace Adsorption {
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //Cercignani-Lampis-Lord model
   namespace CLL {
     extern double TangentialMomentumAccommodationCoefficient[PIC::nTotalSpecies];       //sigma_t, 0..2
     extern double NormalPartKineticEnergyAccommodationCoefficient[PIC::nTotalSpecies];  //alpha_n, 0..1

     //set (and check) the coefficients of a species (spec<0: all species)
     bool SetCoefficients(int spec,double alpha_n,double sigma_t,std::string *error=NULL);

     //the CLL kernel with the given coefficients and wall temperature
     void Reflect(double *v,double alpha_n,double sigma_t,double Temp,int spec,CutCell::cTriangleFace *TriangleCutFace);

     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //Maxwell model
   namespace MaxwellReflection {
     extern double AccommodationCoefficient[PIC::nTotalSpecies];   //0..1

     bool SetAccommodationCoefficient(int spec,double alpha,std::string *error=NULL);

     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //quasi-specular reflection
   namespace QuasiSpecularReflection {
     //quasi-specular CLL: the CLL kernel with small accommodation coefficients (alpha_n 0..1, sigma_t 0..1)
     namespace CLL {
       extern double TangentialMomentumAccommodationCoefficient[PIC::nTotalSpecies];
       extern double NormalPartKineticEnergyAccommodationCoefficient[PIC::nTotalSpecies];

       bool SetCoefficients(int spec,double alpha_n,double sigma_t,std::string *error=NULL);
       int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
     }

     //Nocilla drifting half-range Maxwellian
     namespace Nocilla {
       extern double DriftSpeed[PIC::nTotalSpecies];    //U_r [m/s], >= 0
       extern double DriftAngle[PIC::nTotalSpecies];    //theta_r [rad], 0..pi/2, from the normal
       extern double Temperature[PIC::nTotalSpecies];   //T_r [K], > 0

       bool SetParameters(int spec,double U_r,double theta_r,double T_r,std::string *error=NULL);
       void Reemit(double *v,double U_r,double theta_r,double T_r,int spec,CutCell::cTriangleFace *TriangleCutFace);

       //the parameterization of each species: the direct parameters above, or accommodation coefficients
       const int Direct=0,AccommodationCoefficients=1;
       extern int Mode[PIC::nTotalSpecies];
       extern double NormalMomentumAccommodation[PIC::nTotalSpecies];      //sigma' (Cook & Hoffbauer Eq. 2), 0..1
       extern double TangentialMomentumAccommodation[PIC::nTotalSpecies];  //sigma  (Eq. 1), 0..1
       extern double EnergyAccommodation[PIC::nTotalSpecies];              //epsilon (Eq. 10), 0..1
       extern long int nInconsistentReemissions,nColdLimitReemissions;    //re-emissions with the mean velocity (this process)

       //select the accommodation coefficients for a species (spec<0: all species)
       bool SetAccommodationCoefficients(int spec,double sigma_n,double sigma_t,double epsilon,std::string *error=NULL);

       //the Nocilla parameters for an incident molecule: Lambda = U cos(theta_u)/c_r, c_r = sqrt(2kT_r/m) and
       //Q = U sin(theta_u) along the incident tangential direction; returns 0 (solved), 1 (cold limit: c_r=0,
       //Lambda=+inf) or 2 (inconsistent coefficients: c_r=0, the mean velocity A n + Q t)
       int SolveAccommodation(double vn,double vt,double sigma_n,double sigma_t,double epsilon,double Tw,double mass,
           double *Lambda,double *c_r,double *A,double *Q);

       //the moments <x> and <x^2> of the density ~ x exp(-(x-Lambda)^2), x > 0, for any Lambda
       void FluxNormalMoments(double Lambda,double *mean,double *meanSquare);

       //re-emission with accommodation coefficients
       void ReemitAccommodation(double *v,double sigma_n,double sigma_t,double epsilon,double Tw,int spec,CutCell::cTriangleFace *TriangleCutFace);
       int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

       //a sample x > 0 of the density ~ x exp(-(x-s)^2) (any s)
       double SampleFluxNormalSpeedRatio(double s);
     }

     //Schamberg conical beam
     namespace Schamberg {
       extern double EnergyAccommodationCoefficient[PIC::nTotalSpecies];   //alpha, 0..1
       extern double Nu[PIC::nTotalSpecies];                               //nu > 0 (inf allowed)
       extern double BeamHalfAngle[PIC::nTotalSpecies];                    //phi0 [rad], 0..pi/2

       bool SetParameters(int spec,double alpha,double nu,double phi0,std::string *error=NULL);
       void Reflect(double *v,double alpha,double nu,double phi0,double Temp,int spec,CutCell::cTriangleFace *TriangleCutFace);
       int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

       //Schamberg's beam-shape function Phi(phi0): the mean cosine of the angle between a re-emitted molecule and the axis
       double BeamShapeFunction(double phi0);
     }

     //the generic quasi-specular processor (the compile-time option _QUASI_SPECULAR_REFLECTION_): quasi-specular CLL
     int Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);
   }

   //process gas/surface interaction event with the model of the species
   int ParticleInteractionProcessor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode);

   //the runtime input block 'physics: gas_surface_interaction' (Surface_runtime_input.cpp):
   //
   //   #begin physics: gas_surface_interaction
   //     model = cll                     ! diffuse, specular, maxwell, cll, quasi_specular_cll, nocilla, schamberg,
   //                                     ! adsorption (all species)
   //     wall_temperature_k = 300        ! required by diffuse, maxwell, cll, quasi_specular_cll, schamberg, and nocilla
   //                                     ! with accommodation coefficients
   //     #begin species: H2O             ! the parameters of the selected model, for every species
   //       normal_energy_accommodation = 0.3          ! cll, quasi_specular_cll: alpha_n
   //       tangential_momentum_accommodation = 0.8    ! cll, quasi_specular_cll: sigma_t
   //       ! maxwell:   accommodation
   //       ! nocilla:   drift_speed_m_s, drift_angle_deg, temperature_k, or (accommodation coefficients, Hurlbut &
   //       !            Sherman; requires wall_temperature_k) normal_momentum_accommodation,
   //       !            tangential_momentum_accommodation, energy_accommodation
   //       ! schamberg: energy_accommodation, nu, beam_half_angle_deg
   //     #end species
   //   #end physics
   //
   //('species.H2O.KEY = VALUE' is the same as KEY inside '#begin species: H2O').  Every key of the selected model is
   //required for every species, a key of another model is an error, and nothing is installed unless the whole block
   //is valid.  Apply() also sets the wall temperature (Surface::Temperature::Isothremal::Temp).
   namespace RuntimeInput {
     bool Apply(const AMPS::Input::Block &block,std::string *error);

     //one line per species with the model and its parameters (for run receipts)
     std::string Describe();

     //the parameters of the selected model of every species as "species.SPECIES.KEY=VALUE" lines, with the keys and
     //the units of the runtime input (angles in degrees), for result files that are checked against the input
     std::string ParameterLines();
   }
}



#endif /* _MODELS_SURFACE_SURFACE_H_ */
