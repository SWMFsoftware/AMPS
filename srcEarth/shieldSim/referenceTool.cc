/*
 * shieldSimReferenceTool
 * ----------------------
 * A narrow, machine-readable probe of the Geant4 electromagnetic tables used
 * by Stage 4.  It is intentionally separate from slab transport: NIST
 * PSTAR/ASTAR/XCOM comparisons need to distinguish a physics-table discrepancy
 * from a geometry, scoring, or Monte Carlo discrepancy.
 *
 * The input names every particle/material/energy/quantity query.  No energy is
 * clamped and no material alias is guessed.  Output is in the mass-normalized
 * units used by the NIST tables. The registered tests query exact table nodes,
 * making the Python layer a direct lookup-and-residual calculation with no
 * interpolation or additional unit-conversion model.
 */

#include "StrictJson.hh"

#include <G4Box.hh>
#include <G4EmCalculator.hh>
#include <G4EmParameters.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Event.hh>
#include <G4Geantino.hh>
#include <G4IonQMDPhysics.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleGun.hh>
#include <G4PhysListFactory.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4Version.hh>

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

class ReferenceDetector final : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist=G4NistManager::Instance();
    auto* vacuum=nist->FindOrBuildMaterial("G4_Galactic");
    auto* worldSolid=new G4Box("reference_world",50*cm,50*cm,50*cm);
    auto* worldLogical=new G4LogicalVolume(worldSolid,vacuum,"reference_world");
    auto* world=new G4PVPlacement(nullptr,{},worldLogical,"reference_world",nullptr,false,0,true);

    // Material couples must exist before G4EmCalculator can query CSDA tables.
    // Three disjoint boxes create the validated Stage-4 material domain without
    // introducing a particle source or any scored transport.
    const std::vector<std::string> materials={"G4_Al","G4_Si","G4_WATER"};
    for(std::size_t i=0;i<materials.size();++i){
      auto* material=nist->FindOrBuildMaterial(materials[i]);
      auto* solid=new G4Box((materials[i]+"_solid").c_str(),1*cm,1*cm,1*cm);
      auto* logical=new G4LogicalVolume(solid,material,(materials[i]+"_logical").c_str());
      new G4PVPlacement(nullptr,G4ThreeVector((static_cast<double>(i)-1.0)*5*cm,0,0),
          logical,(materials[i]+"_physical").c_str(),worldLogical,false,
          static_cast<int>(i),true);
    }
    return world;
  }
};

class ReferencePrimary final : public G4VUserPrimaryGeneratorAction {
public:
  ReferencePrimary() : gun_(1) {
    gun_.SetParticleDefinition(G4Geantino::Definition());
    gun_.SetParticleEnergy(1*MeV);
    gun_.SetParticlePosition(G4ThreeVector(0,0,-40*cm));
    gun_.SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  }
  void GeneratePrimaries(G4Event* event) override { gun_.GeneratePrimaryVertex(event); }
private:
  G4ParticleGun gun_;
};

std::string RequiredString(const shield_sim::Json::Object& object,
                           const std::string& name,const std::string& context){
  const auto found=object.find(name);
  if(found==object.end()) throw std::runtime_error(context+" is missing '"+name+"'");
  const std::string value=found->second.AsString(context+"."+name);
  if(value.empty()) throw std::runtime_error(context+"."+name+" is empty");
  return value;
}

double RequiredPositive(const shield_sim::Json::Object& object,
                        const std::string& name,const std::string& context){
  const auto found=object.find(name);
  if(found==object.end()) throw std::runtime_error(context+" is missing '"+name+"'");
  const double value=found->second.AsNumber(context+"."+name);
  if(!std::isfinite(value) || value<=0.0)
    throw std::runtime_error(context+"."+name+" must be finite and > 0");
  return value;
}

double IntegratedMassRange(G4EmCalculator& calculator,double energyMeV,
                           const G4ParticleDefinition* particle,
                           const G4Material* material,double densityGPerCm3){
  // CSDA mass range is the defining integral int dE/(S/rho). Geant4 11.4's
  // reference lists do not populate GetCSDARange for these hadrons in every
  // material couple, so returning that zero would be misleading. Integrating
  // the same unrestricted total stopping model is explicit and independently
  // testable. The lower 1 keV limit matches the lower STAR table boundary; its
  // omitted residual range is negligible for the registered >=10 MeV cases.
  constexpr double lowerMeV=1.0e-3;
  if(energyMeV<=lowerMeV) return 0.0;
  constexpr int intervals=2000; // even, for composite Simpson quadrature
  const double x0=std::log(lowerMeV);
  const double x1=std::log(energyMeV);
  const double step=(x1-x0)/intervals;
  auto transformedIntegrand=[&](double x){
    const double energy=std::exp(x);
    const double stopping=(calculator.ComputeTotalDEDX(
        energy*MeV,particle,material)/(MeV/cm))/densityGPerCm3;
    if(!std::isfinite(stopping) || stopping<=0.0)
      throw std::runtime_error("non-positive stopping power during CSDA quadrature");
    return energy/stopping;
  };
  double sum=transformedIntegrand(x0)+transformedIntegrand(x1);
  for(int i=1;i<intervals;++i)
    sum+=(i%2?4.0:2.0)*transformedIntegrand(x0+i*step);
  return step*sum/3.0;
}

} // namespace

int main(int argc,char** argv){
  try {
    if(argc!=3){
      std::cerr<<"usage: shieldSimReferenceTool INPUT.json OUTPUT.json\n";
      return 2;
    }
    const auto input=shield_sim::Json::ParseFile(argv[1]);
    const auto& root=input.AsObject("reference request");
    shield_sim::RejectUnknownFields(root,{"schema_version","physics_list","em_physics",
                                          "ion_physics","queries"},"reference request");
    if(RequiredString(root,"schema_version","reference request")!=
       "shieldSim.stage4_reference_request.v1")
      throw std::runtime_error("unsupported reference request schema");
    const std::string physicsName=RequiredString(root,"physics_list","reference request");
    const std::string emName=RequiredString(root,"em_physics","reference request");
    const std::string ionName=RequiredString(root,"ion_physics","reference request");
    if(emName!="reference_list_default" && emName!="G4EmStandardPhysics_option4")
      throw std::runtime_error("unsupported em_physics");
    if(ionName!="reference_list_default" && ionName!="G4IonQMDPhysics")
      throw std::runtime_error("unsupported ion_physics");

    auto runManager=std::make_unique<G4RunManager>();
    // Reference physics lists do not build unrestricted CSDA tables unless
    // requested before initialization. Stage 4 compares those tables directly
    // with PSTAR/ASTAR, so enable them explicitly and archive the tool result.
    G4EmParameters::Instance()->SetBuildCSDARange(true);
    runManager->SetUserInitialization(new ReferenceDetector());
    G4PhysListFactory factory;
    auto* physics=factory.GetReferencePhysList(physicsName);
    if(!physics) throw std::runtime_error("Geant4 cannot construct physics_list '"+physicsName+"'");
    if(emName=="G4EmStandardPhysics_option4")
      physics->ReplacePhysics(new G4EmStandardPhysics_option4());
    if(ionName=="G4IonQMDPhysics") physics->ReplacePhysics(new G4IonQMDPhysics());
    physics->SetVerboseLevel(0);
    runManager->SetUserInitialization(physics);
    runManager->SetUserAction(new ReferencePrimary());
    runManager->Initialize();
    // G4EmCalculator's lookup and CSDA tables are finalized during run
    // initialization, not merely kernel initialization. A single non-
    // interacting geantino exercises that lifecycle without introducing an
    // electromagnetic or stochastic contribution to any queried value.
    runManager->BeamOn(1);

    const std::set<std::string> allowedMaterials={"G4_Al","G4_Si","G4_WATER"};
    const std::set<std::string> allowedParticles={"proton","alpha","gamma"};
    const std::set<std::string> allowedQuantities={
      "electronic_stopping_power","nuclear_stopping_power","total_stopping_power",
      "csda_range","mass_attenuation_coefficient"};
    G4EmCalculator calculator;
    shield_sim::Json::Array results;
    const auto& queries=root.at("queries").AsArray("reference request.queries");
    if(queries.empty()) throw std::runtime_error("reference request.queries must not be empty");
    for(std::size_t index=0;index<queries.size();++index){
      const std::string context="reference request.queries["+std::to_string(index)+"]";
      const auto& query=queries[index].AsObject(context);
      shield_sim::RejectUnknownFields(query,{"id","particle","material","energy_MeV","quantities"},context);
      const std::string id=RequiredString(query,"id",context);
      const std::string particleName=RequiredString(query,"particle",context);
      const std::string materialName=RequiredString(query,"material",context);
      const double energyMeV=RequiredPositive(query,"energy_MeV",context);
      if(!allowedMaterials.count(materialName)) throw std::runtime_error(context+" material is outside the Stage-4 validated domain");
      if(!allowedParticles.count(particleName)) throw std::runtime_error(context+" particle is outside the Stage-4 reference domain");
      auto* material=G4NistManager::Instance()->FindOrBuildMaterial(materialName);
      auto* particle=G4ParticleTable::GetParticleTable()->FindParticle(particleName);
      if(!material || !particle) throw std::runtime_error(context+" could not resolve Geant4 particle/material");
      const double density=material->GetDensity()/(g/cm3);
      shield_sim::Json::Object values;
      std::set<std::string> seen;
      for(const auto& entry:query.at("quantities").AsArray(context+".quantities")){
        const std::string quantity=entry.AsString(context+".quantities[]");
        if(!allowedQuantities.count(quantity) || !seen.insert(quantity).second)
          throw std::runtime_error(context+" has unsupported or duplicate quantity '"+quantity+"'");
        if(particleName=="gamma" && quantity!="mass_attenuation_coefficient")
          throw std::runtime_error(context+" gamma supports only mass_attenuation_coefficient");
        if(particleName!="gamma" && quantity=="mass_attenuation_coefficient")
          throw std::runtime_error(context+" mass attenuation requires gamma");
        double value=0.0;
        std::string unit;
        if(quantity=="electronic_stopping_power"){
          value=(calculator.ComputeElectronicDEDX(energyMeV*MeV,particle,material)/(MeV/cm))/density;
          unit="MeV cm2/g";
        } else if(quantity=="nuclear_stopping_power"){
          value=(calculator.ComputeNuclearDEDX(energyMeV*MeV,particle,material)/(MeV/cm))/density;
          unit="MeV cm2/g";
        } else if(quantity=="total_stopping_power"){
          value=(calculator.ComputeTotalDEDX(energyMeV*MeV,particle,material)/(MeV/cm))/density;
          unit="MeV cm2/g";
        } else if(quantity=="csda_range"){
          value=IntegratedMassRange(calculator,energyMeV,particle,material,density);
          unit="g/cm2";
        } else {
          const double length=calculator.ComputeGammaAttenuationLength(energyMeV*MeV,material);
          value=1.0/(length/cm)/density;
          unit="cm2/g";
        }
        if(!std::isfinite(value) || value<0.0)
          throw std::runtime_error(context+" produced a non-finite or negative value");
        values[quantity]=shield_sim::Json(shield_sim::Json::Object{{"unit",unit},{"value",value}});
      }
      results.emplace_back(shield_sim::Json::Object{{"density_g_cm3",density},
        {"energy_MeV",energyMeV},{"id",id},{"material",materialName},
        {"particle",particleName},{"values",shield_sim::Json(std::move(values))}});
    }
    shield_sim::Json output(shield_sim::Json::Object{
      {"csda_range_method","composite_Simpson_quadrature_in_log_energy_of_G4EmCalculator_ComputeTotalDEDX_from_0.001_MeV"},
      {"em_physics",emName},{"geant4_version",G4Version},{"ion_physics",ionName},
      {"physics_list",physicsName},{"results",shield_sim::Json(std::move(results))},
      {"schema_version","shieldSim.stage4_reference_result.v1"}});
    std::ofstream stream(argv[2]);
    if(!stream) throw std::runtime_error(std::string("cannot write output '")+argv[2]+"'");
    stream<<output.DumpCanonical()<<'\n';
    return 0;
  } catch(const std::exception& error){
    std::cerr<<"shieldSimReferenceTool: "<<error.what()<<'\n';
    return 1;
  }
}
