#include "OutputUtils.hh"

#include "MaterialCatalog.hh"

#include <G4SystemOfUnits.hh>
#include <G4ios.hh>

#include <fstream>
#include <iomanip>

void WriteDoseSweepTecplot(const std::vector<SweepPoint>& data,
                           const std::vector<std::string>& matNames,
                           const std::string& shieldMat,
                           const Options& opts)
{
  const std::string fname="shieldSim_dose_sweep.dat";
  std::ofstream out(fname);
  if(!out){ G4cerr<<"Cannot write "<<fname<<G4endl; return; }

  out<<"TITLE = \"Dose vs Shield Thickness - Geant4/"<<opts.physicsList<<"\"\n";
  out<<"# Physics list: "<<opts.physicsList<<"\n";
  out<<"# Dose_perPrimary columns are Gy/primary.\n";
  out<<"# DoseRate columns are Gy/s if the input source spectrum units are physical.\n";
  out<<"# Beam mode assumes the source spectrum columns are differential beam rates [particles/s/MeV].\n";
  out<<"# Isotropic mode assumes the source spectrum columns are differential intensities [particles/(cm2 s sr MeV)].\n";
  out<<"# In isotropic mode the code applies the pi angular factor and multiplies by the finite source-plane area.\n";
  out<<"# Shield material: "<<DescribeShieldMaterial(shieldMat)<<"\n";
  out<<"# Source mode: "<<opts.sourceMode<<"\n";
  out<<"# Events per point: "<<opts.nEvents<<"\n";
  out<<"# Proton energy range: ["<<opts.eMinProton<<", "<<opts.eMaxProton<<"] MeV total kinetic energy\n";
  out<<"# Alpha  energy range: ["<<opts.eMinAlpha <<", "<<opts.eMaxAlpha <<"] MeV total kinetic energy per alpha\n";
  out<<"# Scoring volumes:";
  for(const auto& n:matNames) out<<" "<<n;
  out<<"\n";
  out<<"# Scoring material descriptions:";
  for(const auto& n:matNames) out<<" "<<DescribeDetectorMaterial(n);
  out<<"\n";

  out<<"VARIABLES = \"Thickness [mm]\" \"Areal_Density [g/cm2]\"";
  for(const auto& n:matNames)
    out<<" \"Dose_"<<n<<"_perPrimary [Gy]\"";
  for(const auto& n:matNames)
    out<<" \"DoseRate_"<<n<<" [Gy/s]\"";
  out<<"\n";

  out<<"ZONE T=\""<<DescribeShieldMaterial(shieldMat)<<" shield | source-mode="<<opts.sourceMode<<" | "
     <<opts.nEvents<<" evt/pt\", I="<<data.size()<<", DATAPACKING=POINT\n";

  out<<std::scientific<<std::setprecision(6);
  for(const auto& pt:data){
    out<<std::setw(14)<<pt.thickMM
       <<std::setw(14)<<pt.arealDensity;
    for(G4double d:pt.dose_Gy)
      out<<std::setw(14)<<d/gray;
    for(G4double r:pt.doseRate_Gy_s)
      out<<std::setw(14)<<r/gray;
    out<<"\n";
  }
  out.close();
  G4cout<<"\nDose sweep written: "<<fname
        <<"  ("<<data.size()<<" thickness points)"<<G4endl;
}
