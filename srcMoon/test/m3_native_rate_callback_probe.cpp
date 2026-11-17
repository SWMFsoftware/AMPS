// Native callback wiring only. Real generated AMPS species/record layout and
// the normal Moon particle callback are used; there is no duplicate rate law,
// particle mover, reaction executor, synthetic observation or product event.
#include "pic.h"
#include "Moon.h"
#include "MoonInput.h"

#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#if _PIC_ELECTRON_IMPACT_IONIZATION_REACTION_MODE_ != _PIC_MODE_OFF_
#error M3 rate probe requires the unqualified particle executor to remain OFF
#endif

int main(int argc,char** argv) {
  PIC::InitMPI();
  SetExitErrorCode(__LINE__,1);
  int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  if (argc!=3) {MPI_Finalize();return 2;}
  Moon::Runtime::Configuration configuration;
  std::string error;
  if (!Moon::Runtime::ParseApplicationInput(argv[1],&configuration,&error) ||
      !Moon::Runtime::InstallConfiguration(configuration,&error)) {
    std::cerr<<error<<'\n';MPI_Finalize();return 3;
  }
  if (!Moon::Runtime::GetElectronImpactTable()) {MPI_Finalize();return 4;}
  struct Pair {const char* name;int neutral;int ion;};
  const Pair species[]={{"He",_HE_SPEC_,_HE_PLUS_SPEC_},
      {"Na",_NA_SPEC_,_NA_PLUS_SPEC_},{"Ne",_NE_SPEC_,_NE_PLUS_SPEC_},
      {"Ar",_AR_SPEC_,_AR_PLUS_SPEC_}};
  std::ofstream output(std::string(argv[2])+".rank"+std::to_string(rank)+".csv");
  output<<std::setprecision(17)
      <<"rank,atom,neutral_index,expected_ion_index,actual_ion_index,callback_rate_s,scalar_rate_s,neutral_mass_kg,ion_mass_kg,neutral_charge_c,ion_charge_c\n";
  int failed=0;
  for (const Pair& atom : species) {
    if (atom.neutral<0 || atom.ion<0) {failed=1;break;}
    // The callback only reads the species field. A zeroed native record is a
    // deterministic input fixture, not an allocated/advanced model particle.
    std::vector<PIC::ParticleBuffer::byte> record(PIC::ParticleBuffer::ParticleDataLength,0);
    PIC::ParticleBuffer::SetI(atom.neutral,record.data());
    int product=-1,scalarProduct=-1;
    const double rate=Moon::ElectronImpactIonizationRate(record.data(),product,nullptr);
    const double scalar=Moon::ElectronImpactIonizationRateForSpecies(atom.neutral,scalarProduct);
    const double neutralMass=PIC::MolecularData::GetMass(atom.neutral);
    const double ionMass=PIC::MolecularData::GetMass(atom.ion);
    const double neutralCharge=PIC::MolecularData::GetElectricCharge(atom.neutral);
    const double ionCharge=PIC::MolecularData::GetElectricCharge(atom.ion);
    // A compiled species ID alone cannot establish a physically valid table.
    // In particular a misspelled input mass key used to leave He+ at zero.
    // Read the actual generated parameters, rather than assuming the intended
    // source declaration survived native preprocessing.
    if (!std::isfinite(neutralMass) || !(neutralMass>0) ||
        !std::isfinite(ionMass) || !(ionMass>0) || neutralCharge!=0.0 ||
        !std::isfinite(ionCharge) || !(ionCharge>0)) failed=1;
    if (product!=atom.ion || scalarProduct!=atom.ion || rate!=scalar) failed=1;
    output<<rank<<','<<atom.name<<','<<atom.neutral<<','<<atom.ion<<','<<product
        <<','<<rate<<','<<scalar<<','<<neutralMass<<','<<ionMass
        <<','<<neutralCharge<<','<<ionCharge<<'\n';
    // The configured registry contains neutral first ionization only. An ion
    // must not accidentally receive another neutral's rate or product index.
    PIC::ParticleBuffer::SetI(atom.ion,record.data());
    product=999;
    if (Moon::ElectronImpactIonizationRate(record.data(),product,nullptr)!=0.0 || product!=-1) failed=1;
  }
  output.flush();if (!output.good()) failed=1;
  int globalFailed=0;MPI_Allreduce(&failed,&globalFailed,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  std::cout<<"M3 native rate callback rank="<<rank<<" status="<<(failed?"FAIL":"PASS")<<'\n';
  MPI_Finalize();return globalFailed;
}
