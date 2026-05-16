/* ============================================================================
 * shieldSim.cc
 *
 * Geant4 energetic-particle shielding prototype
 * ------------------------------------------------
 *
 * Purpose
 * -------
 * shieldSim transports energetic protons and alpha particles through a planar
 * shielding slab and estimates quantities behind the shield:
 *
 *   1. Differential transmitted spectra of protons, alpha particles, and
 *      secondary neutrons immediately after the downstream shield face.
 *
 *   2. Energy-deposition dose in one or more downstream scoring slabs, such as
 *      water and silicon detector layers.
 *
 * Geometry
 * --------
 * The detector is a one-dimensional slab stack embedded in vacuum:
 *
 *   upstream source plane -> shield slab -> scoring slab(s) -> downstream gap
 *
 * The default shield is 2 mm Al.  The default scoring stack is 1 mm water plus
 * 1 mm silicon.  Materials are resolved using Geant4 NIST material names.
 *
 * Source and normalization
 * ------------------------
 * The generator samples one primary particle per event from either a built-in
 * approximate GCR-like spectrum or a user-provided three-column spectrum file:
 *
 *   E[MeV]   protonSpectrum   alphaSpectrum
 *
 * E is total kinetic energy per particle.  Alpha energy is total alpha kinetic
 * energy, not MeV/nucleon.  Sampling weights are proportional to J(E)dE rather
 * than J(E) alone, and the integral of J(E)dE is stored as a source-normalization
 * factor.  RunAction uses this factor to convert raw Monte Carlo spectra in
 * counts/(MeV primary) to source-normalized spectra.
 *
 * Source modes
 * ------------
 *   beam:
 *     Normal-incidence pencil beam, x=y=0, direction +z.
 *
 *   isotropic:
 *     Uniform finite upstream source plane with cosine-law inward-hemisphere
 *     directions, p(mu)=2mu.  If the input spectrum is particles/(cm2 s sr MeV),
 *     the code applies the pi angular factor to get plane-crossing flux.
 *
 * Algorithm
 * ---------
 *   1. Parse command-line options into Options.
 *   2. Create the Geant4 run manager, detector, physics list, primary generator,
 *      run action, event action, and stepping action.
 *   3. For each requested shield thickness:
 *        a. Build or reinitialize the geometry.
 *        b. Run BeamOn(nEvents).
 *        c. During tracking, SteppingAction scores energy deposition in the
 *           scoring slabs and transmitted particles crossing the shield rear
 *           face in shield-local coordinates.
 *        d. At end of run, RunAction computes dose per primary, source-normalized
 *           dose rate, and writes spectra.
 *   4. In sweep mode, write the dose-vs-thickness Tecplot table.
 * ========================================================================== */

#include "CLI.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "OutputUtils.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include <G4Exception.hh>
#include <G4NistManager.hh>
#include <G4PhysListFactory.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VModularPhysicsList.hh>
#include <G4ios.hh>

#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

int main(int argc,char** argv){
  Options opts=ParseArguments(argc,argv);
  if(opts.showHelp){ PrintHelp(); return 0; }

  // Collect scoring names and thicknesses once.  RunAction uses these for
  // reports and output variable names, while DetectorConstruction owns the
  // actual logical volumes.
  std::vector<std::string> scNames;
  std::vector<G4double>    scThick;
  for(const auto& s:opts.scoringMaterials){
    scNames.push_back(s.first);
    scThick.push_back(s.second);
  }

  G4cout<<"============ shieldSim configuration ============"<<G4endl;
  G4cout<<"Spectrum  : "
        <<(opts.spectrumFile.empty()?"GCR Badhwar-O'Neill phi=550MV"
                                    :opts.spectrumFile)<<G4endl;
  G4cout<<"Source    : "<<opts.sourceMode<<G4endl;
  G4cout<<"E proton  : ["<<opts.eMinProton<<", "<<opts.eMaxProton<<"] MeV"<<G4endl;
  G4cout<<"E alpha   : ["<<opts.eMinAlpha <<", "<<opts.eMaxAlpha <<"] MeV"<<G4endl;
  G4cout<<"Events    : "<<opts.nEvents<<" per run"<<G4endl;
  G4cout<<"Scoring   :";
  for(const auto& s:opts.scoringMaterials)
    G4cout<<"  "<<s.first<<":"<<s.second/mm<<"mm";
  G4cout<<G4endl;
  if(opts.doSweep){
    G4cout<<"Mode      : SWEEP"<<G4endl;
    G4cout<<"Material  : "<<opts.sweepMaterial<<G4endl;
    G4cout<<"Thickness : "<<opts.sweepTmin<<" - "<<opts.sweepTmax<<" mm  ("
          <<opts.sweepN<<" points, "<<(opts.sweepLog?"log":"linear")<<")"<<G4endl;
  } else {
    G4cout<<"Mode      : single run"<<G4endl;
    G4cout<<"Shield    : "<<opts.shieldMaterial
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

  G4PhysListFactory factory;
  auto* physList=factory.GetReferencePhysList("FTFP_BERT");
  if(!physList) G4Exception("main","PhysList",FatalException,
                            "Cannot create FTFP_BERT.");
  physList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physList);

  auto* runAction=new RunAction(opts,detector,scNames,scThick);
  auto* primaryAction=new PrimaryGeneratorAction(opts,runAction,detector);
  runAction->SetSourceNormalization(primaryAction->GetSourceNormNoAngular(),
                                    primaryAction->GetSourceAngularFactor(),
                                    primaryAction->GetSourceMode());

  runManager->SetUserAction(primaryAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new SteppingAction(runAction));

  runManager->Initialize();

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
      auto* nist=G4NistManager::Instance();
      auto* mat=nist->FindOrBuildMaterial(opts.sweepMaterial);
      G4double rho_gcc = mat ? mat->GetDensity()/(g/cm3) : 0.;
      G4double areal   = rho_gcc * tMM * 0.1;  // g/cm^2, since 1 mm = 0.1 cm

      SweepPoint sp;
      sp.thickMM       = tMM;
      sp.arealDensity  = areal;
      sp.dose_Gy       = runAction->GetLastDose();
      sp.doseRate_Gy_s = runAction->GetLastDoseRate();
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
