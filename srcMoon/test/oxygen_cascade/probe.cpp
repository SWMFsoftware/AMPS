// Thin adapter: all network and recoil operations call maintained production
// services. RNG belongs to this local mathematical test, never to observations.
#include "../../species/OxygenCascade.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
namespace O=Moon::Species::OxygenCascade;
namespace D=AMPS::ChemicalDatabase;
void need(const D::Status& s){if(!s.ok())throw std::runtime_error(s.message);}
void boolean(const char* key,bool value,bool* comma){if(*comma)std::cout<<',';*comma=true;std::cout<<'"'<<key<<"\":"<<(value?"true":"false");}
long double open(std::mt19937_64& r){return (static_cast<long double>(r()>>12)+0.5L)/4503599627370496.0L;}
int main(int argc,char** argv){try{
 if(argc!=8)throw std::runtime_error("expected root/export/exportSHA/release/releaseSHA/profile/profileSHA");
 D::BundleRequest q;q.root=argv[1];q.exportPath=argv[2];q.expectedExportSha256=argv[3];q.expectedReleaseId=argv[4];q.expectedReleaseSha256=argv[5];q.expectedProfileId=argv[6];q.expectedProfileSha256=argv[7];q.requiredCapabilities={"source-access","table-access"};
 D::Catalogue catalogue;need(D::OpenRuntimeBundle(q,&catalogue));O::Network n;need(n.LoadFrozenSelectedChiantI(catalogue));
 std::cout<<std::setprecision(21)<<"{\"scope\":\"source-bound radiative-network mathematics; physical activation unavailable\",\"guards\":{";
 bool comma=false;O::DecayLaw sentinel;sentinel.totalRateS=731;auto s=n.GetDecayLaw(6,&sentinel);boolean("level-domain-transaction",!s.ok()&&sentinel.totalRateS==731,&comma);
 std::size_t ix=99;s=n.SelectDecay(5,1,&ix);boolean("cdf-upper-exclusion",!s.ok()&&ix==99,&comma);
 s=n.SelectDecay(1,0.5,&ix);boolean("stable-no-fabricated-decay",s.code==D::Code::OutsideDomain&&ix==99,&comma);
 long double wait=731;s=n.WaitingTime(5,0,&wait);boolean("wait-open-support",!s.ok()&&wait==731,&comma);
 O::History h;h.finalLevel=731;s=n.AdvanceRadiativeOnly(5,1,{},&h);boolean("missing-draws-transaction",!s.ok()&&h.finalLevel==731,&comma);
 O::Network old=n;D::Catalogue empty;s=old.LoadFrozenSelectedChiantI(empty);boolean("missing-source-preserves-network",!s.ok()&&old.loaded()&&old.links().size()==9,&comma);
 D::Catalogue no;q.expectedExportSha256=std::string(64,'0');s=D::OpenRuntimeBundle(q,&no);boolean("source-hash-refusal",s.code==D::Code::CorruptSource&&!no.loaded(),&comma);
 O::PhotonConstants constants{6.62607015e-34L,299792458.L,"SI-2019-exact-h","SI-exact-c"};O::StateGap gap;gap.upper=731;auto bad=constants;bad.planckSource.clear();s=n.GetStateGap(8,O::EnergyBasis::Observed,bad,&gap);boolean("missing-constant-source",!s.ok()&&gap.upper==731,&comma);
 std::cout<<"},\"laws\":[";
 for(unsigned l=1;l<=5;++l){if(l>1)std::cout<<',';O::DecayLaw law;need(n.GetDecayLaw(l,&law));std::cout<<"{\"level\":"<<l<<",\"rate\":"<<law.totalRateS<<",\"lifetime\":";if(law.lifetimeS)std::cout<<*law.lifetimeS;else std::cout<<"null";std::cout<<'}';}
 std::cout<<"],\"links\":[";
 for(std::size_t k=0;k<n.links().size();++k){if(k)std::cout<<',';const auto& link=n.links()[k];O::StateGap a,b;need(n.GetStateGap(k,O::EnergyBasis::Observed,constants,&a));need(n.GetStateGap(k,O::EnergyBasis::Theoretical,constants,&b));O::DecayLaw law;need(n.GetDecayLaw(link.upper,&law));long double low=0;for(std::size_t j=0;j<law.links.size();++j){if(law.links[j]==k){need(n.SelectDecay(link.upper,low+law.branchProbabilities[j]/2,&ix));break;}low+=law.branchProbabilities[j];}
 std::cout<<"{\"lower\":"<<link.lower<<",\"upper\":"<<link.upper<<",\"rate\":"<<link.rateS<<",\"branch\":"<<link.branchProbability<<",\"vacuum_angstrom\":"<<link.signedVacuumWavelengthAngstrom<<",\"observed_gap_j\":"<<a.releasedStateEnergyJ<<",\"theoretical_gap_j\":"<<b.releasedStateEnergyJ<<",\"midpoint_selected\":"<<ix<<'}';}
 std::cout<<"],\"recoil\":[";
 const long double cases[3][3]={{2,1,3},{4,3,5},{2,1.5,2}};
 for(unsigned k=0;k<3;++k){if(k)std::cout<<',';O::RecoilInput in;in.physicalParentMassKg=cases[k][0];in.physicalDaughterMassKg=cases[k][1];in.constants={1,cases[k][2],"analytical-exact-h","analytical-exact-c"};in.stateGap.planckJS=1;in.stateGap.speedOfLightMS=in.constants.speedOfLightMS;in.stateGap.planckSource=in.constants.planckSource;in.stateGap.lightSource=in.constants.lightSource;in.stateGap.releasedStateEnergyJ=(in.physicalParentMassKg-in.physicalDaughterMassKg)*in.constants.speedOfLightMS*in.constants.speedOfLightMS;in.parentMassSource="analytical-mass";in.daughterMassSource="analytical-mass";in.isotopeIdentity="MATHEMATICAL-NOT-AN-OXYGEN-ISOTOPE";in.frameIdentity="parent-rest";in.photonDirection={1,0,0};in.massEnergyConsistencyJ={0,1e-15L};in.energyJ={0,1e-15L};in.momentumKgMS={0,1e-15L};in.unitDirectionTolerance=0;O::RecoilLedger v;need(O::SolveRestFramePhotonRecoil(in,&v));std::cout<<"{\"mi\":"<<in.physicalParentMassKg<<",\"mf\":"<<in.physicalDaughterMassKg<<",\"c\":"<<in.constants.speedOfLightMS<<",\"gap\":"<<in.stateGap.releasedStateEnergyJ<<",\"photon_j\":"<<v.photonEnergyJ<<",\"p\":"<<v.photonMomentumMagnitudeKgMS<<",\"ke\":"<<v.atomicKineticEnergyJ<<",\"velocity\":"<<v.atomicVelocityMS[0]<<",\"energy_residual\":"<<v.energyResidualJ<<",\"guards\":{";bool g=false;O::RecoilLedger untouched;untouched.photonEnergyJ=731;auto x=in;x.isotopeIdentity.clear();s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("missing-isotope",!s.ok()&&untouched.photonEnergyJ==731,&g);x=in;x.photonDirection={0,0,0};s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("zero-direction",!s.ok()&&untouched.photonEnergyJ==731,&g);x=in;x.stateGap.releasedStateEnergyJ*=2;s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("mass-gap-mismatch",!s.ok()&&untouched.photonEnergyJ==731,&g);x=in;x.constants.speedOfLightMS*=2;s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("constant-mismatch",!s.ok()&&untouched.photonEnergyJ==731,&g);x=in;x.energyJ.absolute=-1;s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("negative-bound",!s.ok()&&untouched.photonEnergyJ==731,&g);if(k==0){x=in;x.photonDirection={0.358026310430055414845557L,0.631744247396772009217267L,0.687543720006900881743395L};x.unitDirectionTolerance=1e-17L;x.energyJ={0,0};x.momentumKgMS={0,0};s=O::SolveRestFramePhotonRecoil(x,&untouched);boolean("non-axis-zero-bound-refusal",s.code==D::Code::NumericalFailure&&untouched.photonEnergyJ==731,&g);}std::cout<<"}}";}
 std::cout<<"],\"families\":[";const long double times[]={1,100,1e5L,1e10L};unsigned family=0;
 for(unsigned initial:{4u,5u})for(long double duration:times){if(family)std::cout<<',';std::mt19937_64 rng(1846053+family);std::array<unsigned,5> populations{};std::array<unsigned,9> photons{};for(unsigned sample=0;sample<20000;++sample){std::vector<O::UniformPair> draws;for(unsigned i=0;i<4;++i)draws.push_back({open(rng),open(rng)});need(n.AdvanceRadiativeOnly(initial,duration,draws,&h));++populations[h.finalLevel-1];for(unsigned j=0;j<9;++j)photons[j]+=h.photonCounts[j];}std::cout<<"{\"initial\":"<<initial<<",\"duration_s\":"<<duration<<",\"seed\":"<<(1846053+family)<<",\"n\":20000,\"population\":[";for(unsigned i=0;i<5;++i){if(i)std::cout<<',';std::cout<<populations[i];}std::cout<<"],\"photons\":[";for(unsigned i=0;i<9;++i){if(i)std::cout<<',';std::cout<<photons[i];}std::cout<<"]}";++family;}
 std::cout<<"]}\n";return 0;
 }catch(const std::exception& e){std::cerr<<"ERROR oxygen cascade probe: "<<e.what()<<'\n';return 2;}}
