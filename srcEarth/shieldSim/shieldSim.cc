/* ============================================================================
 * shieldSim.cc -- application assembly and independent-case driver
 *
 * Scope
 * -----
 * The versioned Stage-1 path supports spacecraft-anomaly environment-effect
 * analysis and preliminary mission assessment. It is not an instrument-design,
 * device-response, or anatomical-phantom tool. CaseDefinition resolves every
 * unit, dimension, material, source, accuracy request, and transport control
 * before Geant4 is created. A rejected case never reaches transport.
 *
 * A case contains one ordered shield stack and one or more requested absorber
 * slabs. Absorbers are alternatives, not a physical downstream stack. This
 * driver rebuilds the geometry for each stable absorber ID, resets tallies and
 * the random seed, and writes ID-qualified outputs. Consequently an absorber
 * cannot attenuate another absorber and absorber-array permutation is exactly
 * testable. DetectorConstruction builds every shield layer explicitly; the
 * last layer supplies the shield-rear surface used for exit scoring.
 *
 * Materials and geometry
 * ----------------------
 * Versioned cases use data/material_library.v1.json plus immutable case-local
 * records. MaterialFactory instantiates the resolved density, elemental or
 * explicit-isotope composition, state, temperature, and optional mean
 * excitation energy. Templates cannot transport. SlabAbsorberGeometry is the
 * current implementation of the AbsorberGeometry provider boundary; future
 * mesh/CAD/voxel providers must preserve stable sensitive-region IDs, mass
 * properties, material mappings, placement, and canonical geometry hashes.
 *
 * Scoring and accuracy
 * --------------------
 * SteppingAction keeps total and non-ionizing deposition separately and derives
 * ionizing deposition by ledger closure. TID is ionizing energy divided by the
 * finite absorber mass. H100/10 counts exact shield-exit crossing energies for
 * primary protons only; alphas, neutrons, and secondary protons are excluded.
 * EventAction supplies event-level samples to the accuracy controller. BeamOn
 * is issued in bounded batches until all requested statistical conditions are
 * met, a known systematic floor makes the request unachievable, or the maximum
 * history budget is reached. An incomplete combined systematic budget is
 * reported as not_assessable.
 *
 * Output and status
 * -----------------
 * Each absorber gets a resolved-case manifest, a shieldSim.result.v1 JSON, and
 * compatibility spectra/tables. The JSON records hashes, code/build/Geant4
 * provenance, uncertainty state, raw H100/10 counts, deposition ledgers, and
 * an explicit capability status. Published-spectrum validation cases remain
 * labeled as test evidence rather than as physically validated mission cases.
 * DDD, n_eq, and LET are disabled for versioned cases until later roadmap
 * stages. Their historical analytic implementations remain reachable only
 * through the legacy CLI for regression compatibility.
 * ========================================================================== */

#include "CLI.hh"
#include "CaseDefinition.hh"
#include "BuildProvenance.hh"
#include "DetectorConstruction.hh"
#include "ComputedQuantities.hh"
#include "EventAction.hh"
#include "MaterialCatalog.hh"
#include "OutputUtils.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include <G4Exception.hh>
#include <G4PhysListFactory.hh>
#include <G4RunManager.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4VModularPhysicsList.hh>
#include <G4ios.hh>
#include <G4Version.hh>
#include <Randomize.hh>

#include <cmath>
#include <cctype>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string SafeId(std::string value){
  for(char& c:value)
    if(!(std::isalnum(static_cast<unsigned char>(c)) || c=='_' || c=='-')) c='_';
  return value;
}

void WriteResolvedCase(const Options& opts){
  if(!opts.caseDefinition) return;
  shield_sim::Json::Object manifest=opts.caseDefinition->resolvedJson.AsObject();
  const auto& absorber=opts.caseDefinition->absorbers.at(opts.activeAbsorberIndex);
  manifest["active_absorber_id"]=absorber.id;
  manifest["active_absorber_geometry_hash"]=absorber.geometryHash;
  manifest["geant4_version"]=G4Version;
  manifest["code_revision"]=shield_sim_build::GitRevision;
  manifest["source_tree_state"]=shield_sim_build::GitState;
  manifest["compiler"]=shield_sim_build::Compiler;
  manifest["build_type"]=shield_sim_build::BuildType;
  manifest["compiler_flags"]=shield_sim_build::CompilerFlags;
  manifest["platform"]=shield_sim_build::Platform;
  shield_sim::Json::Object dataLibraries;
  for(const char* name:{"G4LEDATA","G4NEUTRONHPDATA","G4PARTICLEXSDATA","G4SAIDXSDATA",
                        "G4ABLADATA","G4ENSDFSTATEDATA","G4INCLDATA","G4PIIDATA",
                        "G4REALSURFACEDATA"}){
    const char* value=std::getenv(name);
    dataLibraries[name]=value ? value : "not_set";
  }
  manifest["geant4_data_libraries"]=shield_sim::Json(std::move(dataLibraries));
  manifest["output_prefix"]=opts.outputPrefix;
  manifest["run_status"]=opts.caseDefinition->sourceFormat=="ams02_hapi_component"
      ? "stage2_published_spectrum_transport_test_not_physical_validation"
      : "stage1_core_not_yet_quantity_validated";
  const std::string path=opts.outputPrefix+"_resolved_case.json";
  std::ofstream output(path);
  if(!output) G4Exception("WriteResolvedCase","Output",FatalException,
                          ("Cannot write resolved case provenance '"+path+"'").c_str());
  output<<shield_sim::Json(std::move(manifest)).DumpCanonical()<<'\n';
}

int RunSimulation(Options opts){

  const std::string caseBasePrefix=opts.outputPrefix;
  if(opts.caseDefinition){
    const auto& absorber=opts.caseDefinition->absorbers.at(opts.activeAbsorberIndex);
    opts.outputPrefix=caseBasePrefix+"__"+SafeId(absorber.id);
    opts.dumpRunSummaryFile=opts.outputPrefix+"_run_summary.dat";
    opts.scoringMaterials.clear();
    opts.scoringMaterials.push_back({absorber.materialId,absorber.thicknessMm*mm});
  }

  // Optional deterministic random seed.  This is especially important for
  // automated tests: source-sampling tests can then compare angular moments
  // against fixed tolerances, and geometry/scoring tests can be reproduced
  // exactly if a regression changes.
  if(opts.useRandomSeed){
    G4Random::setTheSeed(opts.randomSeed);
  }

  // Collect scoring names and thicknesses once.  RunAction uses these for
  // reports and output variable names, while DetectorConstruction owns the
  // actual logical volumes.  The names are kept as the user provided them for
  // compact Tecplot variable names; configuration printouts use the catalog
  // description so it is clear which material definition was selected.
  std::vector<std::string> scNames;
  std::vector<G4double>    scThick;
  for(const auto& s:opts.scoringMaterials){
    scNames.push_back(s.first);
    scThick.push_back(s.second);
  }

  G4cout<<"============ shieldSim configuration ============"<<G4endl;
  G4cout<<"Physics  : "<<opts.physicsList<<G4endl;
  G4cout<<"Spectrum  : "
        <<(opts.spectrumFile.empty()?"GCR Badhwar-O'Neill phi=550MV"
                                    :opts.spectrumFile)<<G4endl;
  G4cout<<"Source    : "<<opts.sourceMode<<G4endl;
  G4cout<<"E proton  : ["<<opts.eMinProton<<", "<<opts.eMaxProton<<"] MeV"<<G4endl;
  G4cout<<"E alpha   : ["<<opts.eMinAlpha <<", "<<opts.eMaxAlpha <<"] MeV"<<G4endl;
  G4cout<<"Events    : "<<opts.nEvents<<" per run"<<G4endl;
  G4cout<<"Output    : prefix="<<opts.outputPrefix<<G4endl;
  if(opts.useRandomSeed) G4cout<<"Random    : seed="<<opts.randomSeed<<G4endl;
  if(!opts.dumpSourceSamplesFile.empty())
    G4cout<<"Diag src  : "<<opts.dumpSourceSamplesFile<<G4endl;
  if(!opts.dumpExitParticlesFile.empty())
    G4cout<<"Diag exit : "<<opts.dumpExitParticlesFile<<G4endl;
  if(!opts.dumpRunSummaryFile.empty())
    G4cout<<"Diag sum  : "<<opts.dumpRunSummaryFile<<G4endl;
  if(opts.productionCut>0.0)
    G4cout<<"Prod cut  : "<<opts.productionCut/mm<<" mm"<<G4endl;
  if(opts.maxStepLength>0.0)
    G4cout<<"Max step  : "<<opts.maxStepLength/mm<<" mm"<<G4endl;
  ComputedQuantities::Selection qsel;
  qsel.tid      = opts.calcTID;
  qsel.ddd      = opts.calcDDD;
  qsel.neq      = opts.calcNEq;
  qsel.let      = opts.calcLET;
  qsel.hardness = opts.calcHardness;
  G4cout<<"Quantities: "<<ComputedQuantities::SelectionSummary(qsel)<<G4endl;
  G4cout<<"Scoring   :";
  for(const auto& s:opts.scoringMaterials)
    G4cout<<"  "<<DescribeDetectorMaterial(s.first)<<":"<<s.second/mm<<"mm";
  G4cout<<G4endl;
  if(opts.doSweep){
    G4cout<<"Mode      : SWEEP"<<G4endl;
    G4cout<<"Material  : "<<DescribeShieldMaterial(opts.sweepMaterial)<<G4endl;
    G4cout<<"Thickness : "<<opts.sweepTmin<<" - "<<opts.sweepTmax<<" mm  ("
          <<opts.sweepN<<" points, "<<(opts.sweepLog?"log":"linear")<<")"<<G4endl;
  } else {
    G4cout<<"Mode      : single run"<<G4endl;
    G4cout<<"Shield    : "<<DescribeShieldMaterial(opts.shieldMaterial)
          <<", "<<opts.shieldThickness/mm<<" mm"<<G4endl;
  }
  G4cout<<"================================================="<<G4endl;

  // Build the requested thickness list.  Single-run mode uses the same loop
  // with exactly one element so that the control flow remains nearly identical.
  std::vector<G4double> thicknesses;
  if(opts.doSweep){
    if(opts.sweepN==1){
      thicknesses.push_back(opts.sweepTmin);
    } else {
      for(G4int i=0;i<opts.sweepN;++i){
        G4double frac=G4double(i)/(opts.sweepN-1);
        G4double t;
        if(opts.sweepLog)
          t=opts.sweepTmin*std::pow(opts.sweepTmax/opts.sweepTmin,frac);
        else
          t=opts.sweepTmin+(opts.sweepTmax-opts.sweepTmin)*frac;
        thicknesses.push_back(t);
      }
    }
  } else {
    thicknesses.push_back(opts.shieldThickness/mm);
  }

  // Set the initial geometry for sweep mode before DetectorConstruction is
  // created.  Later points are handled by SetShieldThickness + ReinitializeGeometry.
  if(opts.doSweep){
    opts.shieldMaterial  = opts.sweepMaterial;
    opts.shieldThickness = thicknesses[0]*mm;
  }

  auto* runManager = new G4RunManager();

  auto* detector = new DetectorConstruction(opts);
  runManager->SetUserInitialization(detector);

  // Create the requested Geant4 reference physics list.  The CLI parser
  // validates the allowed names, but we still check the factory return value
  // here so that build/runtime differences in Geant4 installations fail
  // cleanly with an explicit message.
  G4PhysListFactory factory;
  auto* physList=factory.GetReferencePhysList(opts.physicsList);
  if(!physList){
    const G4String msg = G4String("Cannot create requested physics list: ") + opts.physicsList;
    G4Exception("main","PhysList",FatalException,msg.c_str());
  }
  // Optional production range cut and step limiter for numerical convergence
  // testing.  The range cut controls when Geant4 creates explicit secondaries
  // instead of depositing their energy locally.  The step limiter is registered
  // only when --max-step is requested; DetectorConstruction then attaches the
  // corresponding G4UserLimits to the shield and scoring slabs.
  if(opts.productionCut>0.0)
    physList->SetDefaultCutValue(opts.productionCut);
  if(opts.maxStepLength>0.0)
    physList->RegisterPhysics(new G4StepLimiterPhysics());

  physList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physList);

  auto* runAction=new RunAction(opts,detector,scNames,scThick);
  auto* primaryAction=new PrimaryGeneratorAction(opts,runAction,detector);
  runAction->SetSourceNormalization(primaryAction->GetSourceNormNoAngular(),
                                    primaryAction->GetSourceAngularFactor(),
                                    primaryAction->GetSourceMode());

  runManager->SetUserAction(primaryAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(new EventAction(runAction));
  runManager->SetUserAction(new SteppingAction(runAction));

  runManager->Initialize();

  if(opts.caseDefinition){
    const auto& absorbers=opts.caseDefinition->absorbers;
    for(std::size_t absorberIndex=0;absorberIndex<absorbers.size();++absorberIndex){
      const auto& absorber=absorbers[absorberIndex];
      const std::string prefix=caseBasePrefix+"__"+SafeId(absorber.id);
      if(absorberIndex>0){
        detector->SetActiveAbsorber(absorberIndex);
        runAction->ConfigureIndependentCase(absorberIndex,prefix,absorber.materialId,
                                            absorber.thicknessMm*mm);
        runManager->ReinitializeGeometry();
      }
      Options manifestOptions=opts;
      manifestOptions.activeAbsorberIndex=absorberIndex;
      manifestOptions.outputPrefix=prefix;
      WriteResolvedCase(manifestOptions);
      // Reset to the same seed for every absorber. With a common source this
      // makes results exactly invariant to absorber array ordering, while each
      // absorber still has an independent Geant4 geometry and tally state.
      G4Random::setTheSeed(opts.randomSeed);
      const long maximum=opts.caseDefinition->accuracy.maxHistories;
      const long batch=opts.caseDefinition->accuracy.batchSize;
      while(runAction->GetCumulativeEvents()<maximum){
        const long remaining=maximum-runAction->GetCumulativeEvents();
        runManager->BeamOn(static_cast<G4int>(std::min(batch,remaining)));
        if(runAction->AccuracyStatisticallySatisfied()) break;
      }
    }
    delete runManager;
    G4cout<<"\nAll independent absorber cases completed."<<G4endl;
    return 0;
  }

  for(G4int step=0;step<(G4int)thicknesses.size();++step){
    G4double tMM=thicknesses[step];

    if(step>0){
      detector->SetShieldThickness(tMM*mm);
      runManager->ReinitializeGeometry();
    }

    if(opts.doSweep){
      runAction->SetSweepMode(true,opts.sweepMaterial,tMM);
      G4cout<<"\n>>> Sweep point "<<step+1<<"/"<<thicknesses.size()
            <<"  t = "<<std::fixed<<std::setprecision(2)<<tMM<<" mm"<<G4endl;
    }

    runManager->BeamOn(opts.nEvents);

    if(opts.doSweep){
      auto* mat=FindOrBuildShieldMaterial(opts.sweepMaterial);
      G4double rho_gcc = mat ? mat->GetDensity()/(g/cm3) : 0.;
      G4double areal   = rho_gcc * tMM * 0.1;  // g/cm^2, since 1 mm = 0.1 cm

      SweepPoint sp;
      sp.thickMM       = tMM;
      sp.arealDensity  = areal;
      sp.dose_Gy       = runAction->GetLastDose();
      sp.doseRate_Gy_s = runAction->GetLastDoseRate();
      sp.ddd_MeV_g_perPrimary = runAction->GetLastDDD();
      sp.dddRate_MeV_g_s      = runAction->GetLastDDDRate();
      sp.neq_cm2_perPrimary   = runAction->GetLastNEq();
      sp.neqRate_cm2_s        = runAction->GetLastNEqRate();
      sp.hardness_H100_10     = runAction->GetLastHardness();
      runAction->AppendSweepPoint(sp);
    }
  }

  if(opts.doSweep)
    WriteDoseSweepTecplot(runAction->GetSweepData(),scNames,
                          opts.sweepMaterial,opts);

  delete runManager;
  G4cout<<"\nSimulation completed."<<G4endl;
  return 0;
}

} // namespace

int main(int argc,char** argv){
  Options opts=ParseArguments(argc,argv);
  if(opts.showHelp){ PrintHelp(); return 0; }
  if(opts.listMaterials){
    G4cout<<"Available shield-material keys and aliases:\n"
          <<ShieldMaterialCatalogText()<<G4endl;
    return 0;
  }
  if(opts.listTargetMaterials){
    G4cout<<"Available detector/absorber/target-material keys and aliases:\n"
          <<DetectorMaterialCatalogText()<<G4endl;
    return 0;
  }
  if(opts.listQuantities){
    G4cout<<ComputedQuantities::QuantityCatalogText()<<G4endl;
    return 0;
  }

  if(!opts.caseDefinition) return RunSimulation(opts);

  if(opts.validateCaseOnly){
    const auto& absorbers=opts.caseDefinition->absorbers;
    for(std::size_t i=0;i<absorbers.size();++i){
      Options independent=opts;
      independent.activeAbsorberIndex=i;
      independent.outputPrefix=opts.outputPrefix+"__"+SafeId(absorbers[i].id);
      WriteResolvedCase(independent);
    }
    G4cout<<"Validated case '"<<opts.caseDefinition->caseId<<"' with "
          <<absorbers.size()<<" independent absorber case(s)."<<G4endl;
    return 0;
  }
  return RunSimulation(opts);
}
