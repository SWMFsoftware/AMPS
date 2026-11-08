//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>


#include "pic.h"
#include "constants.h"
#include "Exosphere.h"
#include "Dust.h"
#include "Dust.dfn"

using namespace ElectricallyChargedDust::Charging;
using namespace std;

int main(int argc,char **argv) {
  //initialize MPI
  PIC::InitMPI();

  rnd_seed();
  PIC::Init_BeforeParser();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  //set distanc to Sun
  Exosphere::xSun_SO[0]=2.0*_AU_;

  // The nightly adapter supplies a private artifact path so this executable
  // never overwrites a user's PT/plots result.  Direct legacy runs retain the
  // historical OutputDataFileDirectory destination.  This environment
  // variable changes only the diagnostic receipt, never the charging physics.
  const char* TestOutputOverride=getenv("AMPS_TEST_MODEL_DUST_OUTPUT");
  string TestFileName=(TestOutputOverride!=NULL && TestOutputOverride[0]!='\0') ?
      string(TestOutputOverride) :
      string(PIC::OutputDataFileDirectory)+"/test_model-dust.dat";

  ofstream OutFile;
  if(PIC::ThisThread==0)
    OutFile.open(TestFileName.c_str());

  bool TestPassed=true;
  if (PIC::ThisThread==0) {
    OutFile<<setprecision(17);
    cout<<"[CASE] RUN  CORE-DUST-CHARGING-DYNAMICS-001 "
          "configured equilibrium and transient grain charging"<<endl;
    if (!OutFile.is_open()) TestPassed=false;
  }

  //grain parameters
  double GrainElectricChargeInit = 0.0;
  double GrainElectricCharge, EquilibriumCharge;
  const int nSize=4;
  double GrainRadii[nSize] = {1E-7, 1E-6, 1E-5, 1E-4};
  //plasma parameters
  double PlasmaVelocity[3] = {4E5,0,0};
  double PlasmaTemperature = 1E6;
  double PlasmaNumberDensity = 1E7;
  double PlasmaPressure = 1.5*Kbol*PlasmaNumberDensity*PlasmaTemperature;
  //time step
  double dt=1E+0;
  //particle velocity
  double ParticleVelocity[3] = {0};
  double ReferenceChargePerRadius=numeric_limits<double>::quiet_NaN();

  // Independent published-equation fixture for every current enabled by the
  // configured dust profile.  The numbers below were evaluated separately
  // from Horanyi (1996), ARA&A 34, equations 4 and 10--15, using the physical
  // constants in constants.h and this test's declared plasma state:
  // a=1 um, ni=ne=1e7 m^-3, Ti=1e6 K, Te=1.5e6 K, ion drift=4e5 m/s,
  // heliocentric distance=2 AU, photoelectron energy=3 eV, and the retained
  // 250 eV/unit-yield secondary-emission closure.  They are immutable numeric
  // references, not values emitted by GetGrainCurrent or the equilibrium
  // solver.  Comparing signed component currents at negative, zero, and
  // positive potential detects sign/branch mistakes that a net-current root
  // or q/a dimensional scaling alone cannot expose.
  struct PublishedCurrentReference {
    double potential_v;
    double electron_a;
    double ion_a;
    double photoelectron_a;
    double secondary_a;
  };
  const PublishedCurrentReference PublishedCurrents[] = {
    {-5.0,-3.68446179158128898E-17,2.12927679452796571E-18,
      3.14586632897358053E-17,3.24913339871086425E-17},
    { 0.0,-3.82977536529403464E-17,2.11722342997076852E-18,
      3.14586632897358053E-17,3.37727780957568249E-17},
    { 5.0,-3.97791754485218468E-17,2.10517282823571167E-18,
      5.94177399331272733E-18,4.40558203279479803E-16}
  };
  const double ReferenceRadius=1.0E-6;
  const double ElectronTemperature=
      PlasmaPressure/(Kbol*PlasmaNumberDensity);
  for (size_t i=0;
       i<sizeof(PublishedCurrents)/sizeof(PublishedCurrents[0]);++i) {
    const double charge=4.0*Pi*VacuumPermittivity*ReferenceRadius*
        PublishedCurrents[i].potential_v;
    double electronCurrent,electronDerivative;
    double ionCurrent,ionDerivative;
    double photoCurrent,photoDerivative;
    double secondaryCurrent,secondaryDerivative;
    GetGrainCurrent(ReferenceRadius,charge,ParticleVelocity,2.0,
        PlasmaNumberDensity,PlasmaTemperature,PlasmaVelocity,
        PlasmaNumberDensity,ElectronTemperature,
        electronCurrent,electronDerivative,ionCurrent,ionDerivative,
        photoCurrent,photoDerivative,secondaryCurrent,secondaryDerivative);

    const double calculated[4]={electronCurrent,ionCurrent,
        photoCurrent,secondaryCurrent};
    const double expected[4]={PublishedCurrents[i].electron_a,
        PublishedCurrents[i].ion_a,PublishedCurrents[i].photoelectron_a,
        PublishedCurrents[i].secondary_a};
    for (int component=0;component<4;++component) {
      const double scale=std::max(std::fabs(expected[component]),1.0E-30);
      if (!std::isfinite(calculated[component]) ||
          std::fabs(calculated[component]-expected[component])>
              2.0E-12*scale) {
        TestPassed=false;
        if (PIC::ThisThread==0)
          cerr<<setprecision(17)<<"published dust-current mismatch: phi="
              <<PublishedCurrents[i].potential_v<<" component="<<component
              <<" calculated="<<calculated[component]
              <<" expected="<<expected[component]<<endl;
      }
    }
    if(PIC::ThisThread==0)
      OutFile<<"Published-current fixture phi="
          <<PublishedCurrents[i].potential_v<<" V: Je="<<electronCurrent
          <<" A, Ji="<<ionCurrent<<" A, Jpe="<<photoCurrent
          <<" A, Jse="<<secondaryCurrent<<" A\n";
  }

  // Independent bisection of the same published current equations gives the
  // stable positive floating potential 13.64710731500525 V for this profile.
  // The production Newton iteration declares a 1e-5 relative charge-increment
  // stopping scale, so the acceptance tolerance is that solver contract, not
  // roundoff equality with the high-accuracy reference root.
  const double PublishedEquilibriumPotentialV=13.64710731500525;

  for(int iSize=0; iSize<nSize; iSize++){
    double GrainRadius = GrainRadii[iSize];
    
    //equlibrium charging
    EquilibriumCharge = GrainElectricChargeInit;
    UpdateGrainCharge(dt, ParticleVelocity,
		      PlasmaVelocity, PlasmaPressure,
		      PlasmaTemperature, PlasmaNumberDensity, 
		      GrainRadius, EquilibriumCharge,
		      CHARGE_INTEGRATION_MODE__EQUILIBRIUM_POTENTIAL);

    // For identical plasma state the retained collection/emission currents
    // scale with area (a^2) and depend on charge only through the isolated-
    // sphere potential q/(4*pi*eps0*a).  The equilibrium therefore has q/a
    // invariant across radius.  This is an independent dimensional oracle,
    // not a second evaluation of the Newton residual.  The comparison uses
    // the production equilibrium solver's declared 1e-5 relative increment
    // stopping scale; demanding roundoff equality would overstate the
    // accuracy of that nonlinear solve.  The observed pre-qualification
    // discrepancy (3.61e-7) is retained in the first diagnostic receipt.
    if (!isfinite(EquilibriumCharge) || !(EquilibriumCharge>0.0))
      TestPassed=false;
    const double EquilibriumPotential=
        GetDustGrainPotential(GrainRadius,EquilibriumCharge);
    if (!isfinite(EquilibriumPotential) ||
        fabs(EquilibriumPotential/PublishedEquilibriumPotentialV-1.0)>
            1.0E-5) {
      TestPassed=false;
      if (PIC::ThisThread==0)
        cerr<<setprecision(17)<<"published floating-potential mismatch: radius="
            <<GrainRadius<<" potential_v="<<EquilibriumPotential
            <<" expected_v="<<PublishedEquilibriumPotentialV<<endl;
    }
    const double ChargePerRadius=EquilibriumCharge/GrainRadius;
    if (iSize==0) ReferenceChargePerRadius=ChargePerRadius;
    else if (!isfinite(ChargePerRadius) ||
             fabs(ChargePerRadius/ReferenceChargePerRadius-1.0)>1.0E-5) {
      TestPassed=false;
      if (PIC::ThisThread==0)
        cerr<<setprecision(17)<<"dust q/a invariant failure: radius="
            <<GrainRadius<<" q/a="<<ChargePerRadius<<" reference="
            <<ReferenceChargePerRadius<<endl;
    }
    
    if(PIC::ThisThread==0){
      OutFile<<"\nCharging of a grain with radius "<<GrainRadius<<" m"<<
	"-------------------------------\n"<<
      "Equilibrium charge:   "<<EquilibriumCharge<<endl;
    }

    //time dependent charging
    GrainElectricCharge=GrainElectricChargeInit;
    int iIter=0;
    const int MaxChargingIterations=200000;
    double RelativeEquilibriumError;
    do{
      UpdateGrainCharge(dt,ParticleVelocity, 
			PlasmaVelocity, PlasmaPressure, 
			PlasmaTemperature, PlasmaNumberDensity, 
			GrainRadius, GrainElectricCharge, 
			CHARGE_INTEGRATION_MODE__TIME_DEPENDENT);
      
      ++iIter;
      if(PIC::ThisThread==0)
	OutFile<<"Charging iteration "<<iIter<<": " <<
	  GrainElectricCharge<<endl;
      RelativeEquilibriumError=fabs(
          GrainElectricCharge/EquilibriumCharge-1.0);
    }while (isfinite(RelativeEquilibriumError) &&
             RelativeEquilibriumError>1E-5 &&
             iIter<MaxChargingIterations);

    // The transient implicit update must converge to the independently
    // computed steady-current root.  A nonfinite iterate, failure to meet the
    // declared relative tolerance, or iteration-cap exhaustion is a test
    // failure rather than an unbounded loop.
    if (!isfinite(GrainElectricCharge) ||
        !isfinite(RelativeEquilibriumError) ||
        RelativeEquilibriumError>1.0E-5 ||
        iIter>=MaxChargingIterations) {
      TestPassed=false;
      if (PIC::ThisThread==0)
        cerr<<setprecision(17)<<"dust transient convergence failure: radius="
            <<GrainRadius<<" equilibrium="<<EquilibriumCharge<<" final="
            <<GrainElectricCharge<<" relative_error="
            <<RelativeEquilibriumError<<" iterations="<<iIter<<endl;
    }
    
  }

  if(PIC::ThisThread==0)
    OutFile.close();

  int LocalFailure=TestPassed ? 0 : 1;
  int GlobalFailure=0;
  MPI_Allreduce(&LocalFailure,&GlobalFailure,1,MPI_INT,MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread==0)
    cout<<"[CORE-DUST-CHARGING-DYNAMICS-001] "
        <<(GlobalFailure==0 ? "PASS" : "FAIL")<<endl;
  
  //finish the run
  MPI_Finalize();
  //  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return GlobalFailure==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}





