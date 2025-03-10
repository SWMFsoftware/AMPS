
//$Id$
//functions that are used for interpolation of the dataset obtainedwith other models into AMPS' grid

/*
 * pic_datafile.cpp
 *
 *  Created on: May 5, 2015
 *      Author: vtenishe
 */

#include "pic.h"
//#include <algorithm>

//number of file to be loaded
int PIC::CPLR::DATAFILE::MULTIFILE::nFile=0;

//schedule for loading multiple data files
vector<PIC::CPLR::DATAFILE::MULTIFILE::cScheduleItem> PIC::CPLR::DATAFILE::MULTIFILE::Schedule;

//name of file with table defining schedule
char PIC::CPLR::DATAFILE::MULTIFILE::FileTable[_MAX_STRING_LENGTH_PIC_]="Schedule";

//variable to track whether to break simulation at the last datafile
bool PIC::CPLR::DATAFILE::MULTIFILE::BreakAtLastFile  = true;

//variable to track whether the last datafile has been reached
bool PIC::CPLR::DATAFILE::MULTIFILE::ReachedLastFile  = false;

//offset to the current datafiles relative to beginning of data buffer
int _TARGET_DEVICE_ _CUDA_MANAGED_ PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset = -1;
int PIC::CPLR::DATAFILE::MULTIFILE::NextDataFileOffset = -1;

// next file to load
int PIC::CPLR::DATAFILE::MULTIFILE::iFileLoadNext = -1;

//path to the location of the datafiles
char PIC::CPLR::DATAFILE::path[_MAX_STRING_LENGTH_PIC_]=".";

//the offset from the cell->AssociatedData()
int _TARGET_DEVICE_ _CUDA_MANAGED_ PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin=-1;
int PIC::CPLR::DATAFILE::nTotalBackgroundVariables=0;
bool PIC::CPLR::DATAFILE::Offset::InitFlag=false;

//the unmber of ion fluids used in PIC::CPLR::DATAFILE 
int PIC::CPLR::DATAFILE::nIonFluids=1;

//Physical quantaties offsets that could be read and srored
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity={false,false,1,"Plasma number density",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity={false,false,3,"\"vPlasmaX\", \"vPlasmaY\", \"vPlasmaZ\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaTemperature={false,false,1,"Plasma temperature",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure={false,false,1,"Plasma pressure",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaElectronPressure={false,false,1,"\"Plasma electron pressure\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement _TARGET_DEVICE_ _CUDA_MANAGED_ PIC::CPLR::DATAFILE::Offset::MagneticField={false,false,3,"\"Bx\", \"By\", \"Bz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement _TARGET_DEVICE_ _CUDA_MANAGED_ PIC::CPLR::DATAFILE::Offset::ElectricField={false,false,3,"\"Ex\", \"Ey\", \"Ez\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::PlasmaDivU={false,false,1,"\"Plasma DivU\"",-1};


PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::Current={false,false,3,"\"Jx\", \"Jy\", \"Jz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::b_dot_grad_b={false,false,3,"\"b_dot_grad_bx\", \"b_dot_grad_by\", \"b_dot_grad_bz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b={false,false,3,"\"vE_dot_grad_bx\", \"vE_dot_grad_by\", \"vE_dot_grad_bz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE={false,false,3,"\"b_dot_grad_vEx\", \"b_dot_grad_vEy\", \"b_dot_grad_vEz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE={false,false,3,"\"vE_dot_grad_vEx\", \"vE_dot_grad_vEy\", \"vE_dot_grad_vEz\"",-1};
PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::grad_kappaB={false,false,3,"\"grad_kappaBx\", \"grad_kappaBy\", \"grad_kappaBz\"",-1};

PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient={false,false,9,"\"dBx/dx\", \"dBx/dy\", \"dBx/dz\", \"dBy/dx\", \"dBy/dy\", \"dBy/dz\", \"dBz/dx\", \"dBz/dy\", \"dBz/dz\"",-1};

PIC::CPLR::DATAFILE::cOffsetElement PIC::CPLR::DATAFILE::Offset::MagneticFluxFunction={false,false,1,"\"FluxFunction\"",-1};


//==============================================================================
void PIC::CPLR::DATAFILE::MULTIFILE::Init(bool BreakAtLastFileIn,int  FileNumberFirst) {
  //load schedule from file
  GetSchedule();
  iFileLoadNext = FileNumberFirst;

  CurrDataFileOffset = 0;

  if (_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_) {
    NextDataFileOffset = nTotalBackgroundVariables*sizeof(double);
  } else {
    NextDataFileOffset = 0;
  }

  BreakAtLastFile= BreakAtLastFileIn;
  // set the global time counter value
  PIC::SimulationTime::SetInitialValue(Schedule[iFileLoadNext].Time);
  // load the first file
  UpdateDataFile();
  CopyCurrDataFile2NextDataFile();

  if (_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_) {
    // load the second file
    UpdateDataFile();
  }
}


//=============================================================================
void PIC::CPLR::DATAFILE::MULTIFILE::GetSchedule() {
  char fullname[_MAX_STRING_LENGTH_PIC_];

  //compose full name of the file table file with path
  sprintf(fullname,"%s/%s", PIC::CPLR::DATAFILE::path,PIC::CPLR::DATAFILE::MULTIFILE::FileTable);

  //read the schedule; the format is the following:
  //---------------------------------------------------------------------------
  // Comment lines             | Comment lines           
  // #NFILE                    | #NFILE                  
  // nFile                     | nFile                   
  // #FILELIST                OR #FILESCHEDULE               
  // <name of the file 1>      | <time of file1> <name of the file 1>    
  // ...                       | ...                     
  // <name of the file nFile>  | <time of file nFile> <name of the file nFile>
  // Ignored part of the file  | Ignored part of the file
  //---------------------------------------------------------------------------
  // indicating times is optional provided that reader is able to extract them
  // from the data files themselves, otherwise using #FILESCHEDULE is mandatory
  //---------------------------------------------------------------------------
  // read schedule using class CiFileOperations (see src/general/ifileopr.h)
  CiFileOperations fin;
  // containers for a line from file's contents
  char str[_MAX_STRING_LENGTH_PIC_], str1[_MAX_STRING_LENGTH_PIC_];
  fin.openfile(fullname);

  // read the file's contents; first find number of files
  //---------------------------------------------------------------------------
  while (fin.eof()==false) {
    // get a line from the file
    fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);

    // check if found keyword "#NFILE"
    fin.CutInputStr(str1, str);

    if (strcmp("#NFILE",str1)==0) {
      fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
      nFile = strtol(str, NULL,10);
      break;
    }
  }

  // check if successfully found number of input files
  if (nFile == -1) exit(__LINE__,__FILE__,"Can't read number of input data files from schedule file");

  // check correctness of the found value
  if (nFile < 1) exit(__LINE__,__FILE__,"Number of input data files in schedule file is invalid");
  //---------------------------------------------------------------------------

  // now read the actual table
  //---------------------------------------------------------------------------
  //first, clear the schedule
  Schedule.clear();

  // whether times are provided in the table
  bool IsSchedule;

  //find the beginning of the table
  while (fin.eof()==false) {
    // get a line from the file
    fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_);
    // check if found keyword "#NFILE"
    fin.CutInputStr(str1, str);
    if(strcmp("#FILELIST",    str1)==0) {IsSchedule = false; break;}
    if(strcmp("#FILESCHEDULE",str1)==0) {IsSchedule = true;  break;}
  }

  //check if successfully found the beginning of the table
  if (fin.eof()==true) exit(__LINE__,__FILE__,"Can't locate the beginning of the actual file table in the schedule file");

  //read file names
  cScheduleItem Item;

  for (int iFile=0;iFile<nFile;iFile++) {
    fin.GetInputStr(str,_MAX_STRING_LENGTH_PIC_, false);

    if (IsSchedule) {
      //time is provided in the table
      fin.CutInputStr(str1, str);
      Item.Time=strtod(str1, NULL);
      sprintf(Item.FileName,"%s",str);
    }
    else {
      //time has to be extracted from the file itself
      sprintf(Item.FileName,"%s",str);
      Item.Time=GetFileTime(str);
    }

    Schedule.push_back(Item);
  }

  //---------------------------------------------------------------------------
  // the rest of the file is ignored; close file
  fin.closefile();

  // as a final step, need to sort schedule by time
  sort(Schedule.begin(), Schedule.end(), _compare);
}

//=============================================================================
double PIC::CPLR::DATAFILE::MULTIFILE::GetFileTime(const char* FileName) {
  double res=-1.0;

  //the particular  reader
  switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
  case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
    exit(__LINE__,__FILE__,"TECPLOT reader mode is not able to extract time from the input data file");
    break;

  case _PIC_COUPLER_DATAFILE_READER_MODE__ARMS_:
    res=ARMS::GetFileTime(FileName);
    break;

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  return res; 
}

//=============================================================================
void PIC::CPLR::DATAFILE::MULTIFILE::UpdateDataFile() {
  if (MULTIFILE::ReachedLastFile==false) {
    if (_PIC_DATAFILE__TIME_INTERPOLATION_MODE_ == _PIC_MODE_ON_) {
      //swap data offsets
      if (CurrDataFileOffset==0) {
        CurrDataFileOffset=NextDataFileOffset;
        NextDataFileOffset=0;
      }
      else {
        NextDataFileOffset=CurrDataFileOffset;
        CurrDataFileOffset=0;
      }
    }

    //load the new data file
    if (iFileLoadNext<nFile) {
      PIC::CPLR::DATAFILE::ImportData(Schedule[iFileLoadNext].FileName);

      if (PIC::ThisThread==0) std::cout << "Data file " << Schedule[iFileLoadNext].FileName << " has been loaded" << std::endl;
      iFileLoadNext++;
    }

    if (iFileLoadNext>=nFile) ReachedLastFile=true;  //the last file has been reached
  }
}


//=============================================================================
//load new data file
//IMPORTANT! The list of the data that are loaded has to be indicated before PIC::Init_BeforeParser
void PIC::CPLR::DATAFILE::ImportData(const char *fname) {

  //the particular  reader
  switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
  case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
    TECPLOT::ImportData(fname);
    break;
  case _PIC_COUPLER_DATAFILE_READER_MODE__ARMS_:
    ARMS::LoadDataFile(fname);
    break;
  case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
    BATSRUS::LoadDataFile(fname);
    break;

#ifndef _NO_KAMELEON_CALLS_
  case _PIC_COUPLER_DATAFILE_READER_MODE__KAMELEON_:
    KAMELEON::LoadDataFile(fname);
    break;
#endif

  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  //may need to generate additional data
  if (Offset::MagneticFieldGradient.allocate) {
    if (_PIC_COUPLER__INTERPOLATION_MODE_==_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_)
      exit(__LINE__,__FILE__,"ERROR: magnetic field gradient can't be computed with 0th order interpolation method");

    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
      GenerateMagneticFieldGradient(node);
    }
  }

  if ( _PIC_MOVER_INTEGRATOR_MODE_ == _PIC_MOVER_INTEGRATOR_MODE__RELATIVISTIC_GCA_){
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
      GenerateVarForRelativisticGCA(node);
    }   
  }

}

//initialize the data reading namespace
void PIC::CPLR::DATAFILE::Init() {

  //initialize offset to the current datafile storage
  MULTIFILE::CurrDataFileOffset = 0;


  /*  // put a dummy item in the schedule
  MULTIFILE::nFile = 1; MULTIFILE::iFileLoadNext=0;
  MULTIFILE::cScheduleItem ItemDummy;
  ItemDummy.Time=NAN;
  sprintf(ItemDummy.FileName,"%s","dummy_name");
  MULTIFILE::Schedule.push_back(ItemDummy);*/

  //set the initialization flag and call the init procedures of the particular file reader
  Offset::InitFlag=true;

  //initialize common offsets

  PIC::CPLR::DATAFILE::Offset::ElectricField.allocate=true;
  PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;

  if ( _PIC_MOVER_INTEGRATOR_MODE_ == _PIC_MOVER_INTEGRATOR_MODE__RELATIVISTIC_GCA_){
    PIC::CPLR::DATAFILE::Offset::Current.allocate=true;

    PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.allocate=true;
    PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b.allocate=true;   
    PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE.allocate=true;   
    PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE.allocate=true;   
    PIC::CPLR::DATAFILE::Offset::grad_kappaB.allocate=true; 
  }

  PIC::CPLR::DATAFILE::Offset::PlasmaIonPressure.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaNumberDensity.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaTemperature.allocate=true;
  PIC::CPLR::DATAFILE::Offset::PlasmaBulkVelocity.allocate=true;

  //offset for magnetic field gradient (if needed)
  if ((_PIC_MOVER_INTEGRATOR_MODE_==_PIC_MOVER_INTEGRATOR_MODE__GUIDING_CENTER_)||(_PIC_OUTPUT__DRIFT_VELOCITY__MODE_==_PIC_MODE_ON_)) { 
    PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate=true;
  }


  if (_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__DATAFILE_) {
    switch (_PIC_COUPLER_DATAFILE_READER_MODE_) {
    case _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_:
      TECPLOT::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__ARMS_:
      ARMS::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__ICES_:
      ICES::Init();
      break;
    case _PIC_COUPLER_DATAFILE_READER_MODE__KAMELEON_:
      KAMELEON::Init();
      break;
#if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_
    case _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_:
      BATSRUS::Init();
      break;
#endif //  _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__BATSRUS_

    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  }
  else if ((_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__T96_)||(_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__T05_)||(_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__KMAG_)) {
    PIC::CPLR::DATAFILE::Offset::MagneticField.allocate=true;
  }
  else exit(__LINE__,__FILE__,"Error: wrong option");

  //init the offset table and request memory
  int RelativeOffset=0;
  CenterNodeAssociatedDataOffsetBegin=PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength;

  //init the data offsets
  if (Offset::PlasmaNumberDensity.allocate==true) {
    Offset::PlasmaNumberDensity.active=true;
    Offset::PlasmaNumberDensity.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=nIonFluids*Offset::PlasmaNumberDensity.nVars*sizeof(double);
    RelativeOffset+=nIonFluids*Offset::PlasmaNumberDensity.nVars*sizeof(double);
    nTotalBackgroundVariables+=nIonFluids*Offset::PlasmaNumberDensity.nVars;
  }


  if (Offset::PlasmaBulkVelocity.allocate==true) {
    if (Offset::PlasmaBulkVelocity.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaBulkVelocity");

    Offset::PlasmaBulkVelocity.active=true;
    Offset::PlasmaBulkVelocity.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=nIonFluids*Offset::PlasmaBulkVelocity.nVars*sizeof(double);
    RelativeOffset+=nIonFluids*Offset::PlasmaBulkVelocity.nVars*sizeof(double);
    nTotalBackgroundVariables+=nIonFluids*Offset::PlasmaBulkVelocity.nVars;
  }


  if (Offset::PlasmaTemperature.allocate==true) {
    if (Offset::PlasmaTemperature.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaTemperature");

    Offset::PlasmaTemperature.active=true;
    Offset::PlasmaTemperature.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=nIonFluids*Offset::PlasmaTemperature.nVars*sizeof(double);
    RelativeOffset+=nIonFluids*Offset::PlasmaTemperature.nVars*sizeof(double);
    nTotalBackgroundVariables+=nIonFluids*Offset::PlasmaTemperature.nVars;
  }


  if (Offset::PlasmaIonPressure.allocate==true) {
    if (Offset::PlasmaIonPressure.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaIonPressure");

    Offset::PlasmaIonPressure.active=true;
    Offset::PlasmaIonPressure.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=nIonFluids*Offset::PlasmaIonPressure.nVars*sizeof(double);
    RelativeOffset+=nIonFluids*Offset::PlasmaIonPressure.nVars*sizeof(double);
    nTotalBackgroundVariables+=nIonFluids*Offset::PlasmaIonPressure.nVars;
  }

  if (Offset::PlasmaDivU.allocate==true) {
    if (Offset::PlasmaDivU.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaIonPressure");

    Offset::PlasmaDivU.active=true;
    Offset::PlasmaDivU.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::PlasmaDivU.nVars*sizeof(double);
    RelativeOffset+=sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaDivU.nVars;
  }

  if (Offset::PlasmaElectronPressure.allocate==true) {
    if (Offset::PlasmaElectronPressure.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaElectronPressure");

    Offset::PlasmaElectronPressure.active=true;
    Offset::PlasmaElectronPressure.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::PlasmaElectronPressure.nVars*sizeof(double);
    RelativeOffset+=Offset::PlasmaElectronPressure.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::PlasmaElectronPressure.nVars;
  }

  if (Offset::MagneticField.allocate==true) {
    if (Offset::MagneticField.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::MagneticField");

    Offset::MagneticField.active=true;
    Offset::MagneticField.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::MagneticField.nVars*sizeof(double);
    RelativeOffset+=Offset::MagneticField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::MagneticField.nVars;
  }

  if (Offset::ElectricField.allocate==true) {
    if (Offset::ElectricField.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::ElectricField");

    Offset::ElectricField.active=true;
    Offset::ElectricField.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::ElectricField.nVars*sizeof(double);
    RelativeOffset+=Offset::ElectricField.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::ElectricField.nVars;
  }


  if (Offset::Current.allocate==true) {
    if (Offset::Current.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::Current");

    Offset::Current.active=true;
    Offset::Current.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::Current.nVars*sizeof(double);
    RelativeOffset+=Offset::Current.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::Current.nVars;
  }

  if (Offset::b_dot_grad_b.allocate==true) {
    if (Offset::b_dot_grad_b.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::b_dot_grad_b");

    Offset::b_dot_grad_b.active=true;
    Offset::b_dot_grad_b.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::b_dot_grad_b.nVars*sizeof(double);
    RelativeOffset+=Offset::b_dot_grad_b.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::b_dot_grad_b.nVars;
  }

  if (Offset::vE_dot_grad_b.allocate==true) {
    if (Offset::vE_dot_grad_b.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::vE_dot_grad_b");

    Offset::vE_dot_grad_b.active=true;
    Offset::vE_dot_grad_b.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::vE_dot_grad_b.nVars*sizeof(double);
    RelativeOffset+=Offset::vE_dot_grad_b.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::vE_dot_grad_b.nVars;
  }

  if (Offset::b_dot_grad_vE.allocate==true) {
    if (Offset::b_dot_grad_vE.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::b_dot_grad_vE");

    Offset::b_dot_grad_vE.active=true;
    Offset::b_dot_grad_vE.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::b_dot_grad_vE.nVars*sizeof(double);
    RelativeOffset+=Offset::b_dot_grad_vE.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::b_dot_grad_vE.nVars;
  }

  if (Offset::vE_dot_grad_vE.allocate==true) {
    if (Offset::vE_dot_grad_vE.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::PlasmaElectronPressure");

    Offset::vE_dot_grad_vE.active=true;
    Offset::vE_dot_grad_vE.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::vE_dot_grad_vE.nVars*sizeof(double);
    RelativeOffset+=Offset::vE_dot_grad_vE.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::vE_dot_grad_vE.nVars;
  }

  if (Offset::grad_kappaB.allocate==true) {
    if (Offset::grad_kappaB.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::grad_kappaB");

    Offset::grad_kappaB.active=true;
    Offset::grad_kappaB.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::grad_kappaB.nVars*sizeof(double);
    RelativeOffset+=Offset::grad_kappaB.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::grad_kappaB.nVars;
  }


  if (Offset::MagneticFieldGradient.allocate==true) {
    if (Offset::MagneticFieldGradient.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::MagneticFieldGradient");

    Offset::MagneticFieldGradient.active=true;
    Offset::MagneticFieldGradient.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::MagneticFieldGradient.nVars*sizeof(double);
    RelativeOffset+=Offset::MagneticFieldGradient.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::MagneticFieldGradient.nVars;
  }

  if (Offset::MagneticFluxFunction.allocate==true) {
    if (Offset::MagneticFluxFunction.RelativeOffset!=-1) exit(__LINE__,__FILE__,"Error: reinitialization of Offset::MagneticFluxFunction");

    Offset::MagneticFluxFunction.active=true;
    Offset::MagneticFluxFunction.RelativeOffset=RelativeOffset;

    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=Offset::MagneticFluxFunction.nVars*sizeof(double);
    RelativeOffset+=Offset::MagneticFluxFunction.nVars*sizeof(double);
    nTotalBackgroundVariables+=Offset::MagneticFluxFunction.nVars;
  }


  if (_PIC_DATAFILE__TIME_INTERPOLATION_MODE_==_PIC_MODE_ON_) {
    // double the reserved memory for time inteprolation mode
    PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength+=nTotalBackgroundVariables*sizeof(double);
  }

  if (nTotalBackgroundVariables==0) {
    exit(__LINE__,__FILE__,"Error: no background variables will be loaded. The background variables to load has to be indicated before call of the PIC::Init_BeforeParser()");
  }

  //print out of the output file
  PIC::Mesh::PrintVariableListCenterNode.push_back(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
}


//save and load the binary data saved in the AMPS data structure
bool PIC::CPLR::DATAFILE::BinaryFileExists(const char *fNameBase) {
  FILE *fData=NULL;
  char fname[400];

  sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",path,PIC::Mesh::mesh->getMeshSignature(),fNameBase);
  fData=fopen(fname,"r");

  if (fData!=NULL) {
    fclose(fData);
    return true;
  }

  return false;
}


void PIC::CPLR::DATAFILE::SaveBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static CMPI_channel pipe;
  static FILE *fout=NULL;

  if (startNode==PIC::Mesh::mesh->rootTree) {
    char fname[400];
    sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",path,PIC::Mesh::mesh->getMeshSignature(),fNameBase);

    pipe.init(1000000);

    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.openRecvAll();
      fout=fopen(fname,"w");
    }
    else pipe.openSend(0);
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char SendCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((PIC::ThisThread==0)||(startNode->Thread==PIC::ThisThread)) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      SendCellFlag=true;

      //determine whether the cell data neede to be saved
      if (startNode->Thread==PIC::ThisThread) {
        //locate the cell
        if (startNode->block==NULL) SendCellFlag=false,offset=NULL;

        nd=_getCenterNodeLocalNumber(i,j,k);
        if (SendCellFlag==true) {
          if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) offset=CenterNode->GetAssociatedDataBufferPointer();
          else SendCellFlag=false,offset=NULL;
        }

        if (startNode->Thread!=0) pipe.send(SendCellFlag);
      }
      else {
        pipe.recv(SendCellFlag,startNode->Thread);
      }

      //save the cell data saving flag
      if (PIC::ThisThread==0) fwrite(&SendCellFlag,sizeof(char),1,fout);

      //save the cell data
      if (SendCellFlag==true) {
        if (startNode->Thread==PIC::ThisThread) {
          if (startNode->Thread==0) {
            fwrite(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset,nTotalBackgroundVariables*sizeof(double),1,fout);

            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)             {
              PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset),nTotalBackgroundVariables,__LINE__,__FILE__);
            }
          }
          else {
            pipe.send((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset),nTotalBackgroundVariables);
          }
        }
        else {
          double data[nTotalBackgroundVariables];

          pipe.recv(data,nTotalBackgroundVariables,startNode->Thread);
          fwrite(data,nTotalBackgroundVariables*sizeof(double),1,fout);

          if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)           {
            PIC::Debugger::CatchOutLimitValue(data,nTotalBackgroundVariables,__LINE__,__FILE__);
          }
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SaveBinaryFile(NULL,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh->rootTree) {
    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.closeRecvAll();
      fclose(fout);
    }
    else pipe.closeSend();

    pipe.remove();
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }

}


void PIC::CPLR::DATAFILE::LoadBinaryFile(const char *fNameBase,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static FILE *fData=NULL;

  if (startNode==PIC::Mesh::mesh->rootTree) {
    char fname[400];
    sprintf(fname,"%s/amr.sig=0x%lx.f=%s.CenterNodeBackgroundData.bin",path,PIC::Mesh::mesh->getMeshSignature(),fNameBase);

    fData=fopen(fname,"r");
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  char savedLoadCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->block==NULL) {
      //the block belongs to a other processor -> skip all data
      for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++)  {
        if (fread(&savedLoadCellFlag,sizeof(char),1,fData)!=1) exit(__LINE__,__FILE__,"Error: fread failed"); 

        if (savedLoadCellFlag==true) {
          //the cell data is saved -> skip it
          fseek(fData,nTotalBackgroundVariables*sizeof(double),SEEK_CUR);
        }
      }

    }
    else for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
      if (fread(&savedLoadCellFlag,sizeof(char),1,fData)!=1) exit(__LINE__,__FILE__,"Error: fread failed");  

      if (savedLoadCellFlag==true) {
        //determine whether the cell data needed to be read
        //locate the cell
        nd=_getCenterNodeLocalNumber(i,j,k);

        if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) {
          offset=CenterNode->GetAssociatedDataBufferPointer();

          //read the data
          if (fread(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset,nTotalBackgroundVariables*sizeof(double),1,fData)!=1) exit(__LINE__,__FILE__,"Error: fread failed"); 

          if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)           {
            PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset),nTotalBackgroundVariables,__LINE__,__FILE__);
          }
        }
        else {
          //the cell data is saved -> skip it
          fseek(fData,nTotalBackgroundVariables*sizeof(double),SEEK_CUR);
        }
      }

    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) LoadBinaryFile(NULL,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh->rootTree) {
    fclose(fData);

    //initialize derived data
    if (PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.allocate==true) {
      if (_PIC_COUPLER__INTERPOLATION_MODE_==_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_) {
        exit(__LINE__,__FILE__,"ERROR: magnetic field gradient can't be computed with 0th order interpolation method");
      }

      for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        PIC::CPLR::DATAFILE::GenerateMagneticFieldGradient(node);
      }

      //Exchange derived data betwenn the boundary nodes
      PIC::Mesh::mesh->ParallelBlockDataExchange();
    }

    if ( _PIC_MOVER_INTEGRATOR_MODE_ == _PIC_MOVER_INTEGRATOR_MODE__RELATIVISTIC_GCA_){
      for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        PIC::CPLR::DATAFILE::GenerateVarForRelativisticGCA(node);
      }      

      PIC::Mesh::mesh->ParallelBlockDataExchange();
    }
  }
}

//====================================================
//output background parameters loaded with ICES from SWMF output
void PIC::CPLR::DATAFILE::PrintVariableList(FILE* fout,int DataSetNumber) {
  if (Offset::PlasmaNumberDensity.active) {
    if (nIonFluids!=1) for (int i=0;i<nIonFluids;i++) fprintf(fout,", \"%s[%i]\"",Offset::PlasmaNumberDensity.VarList,i);
    else fprintf(fout,", \"%s\"",Offset::PlasmaNumberDensity.VarList);
  }

  if (Offset::PlasmaBulkVelocity.active) {
    if (nIonFluids!=1) for (int i=0;i<nIonFluids;i++) fprintf(fout,", \"vPlasmaX[%i]\", \"vPlasmaY[%i]\", \"vPlasmaZ[%i]\"",i,i,i);
    else fprintf(fout,", \"vPlasmaX\", \"vPlasmaY\", \"vPlasmaZ\"");
  }

  if (Offset::PlasmaTemperature.active) {
    if (nIonFluids!=1) for (int i=0;i<nIonFluids;i++) fprintf(fout,", \"%s[%i]\"",Offset::PlasmaTemperature.VarList,i);
    else fprintf(fout,", \"%s\"",Offset::PlasmaTemperature.VarList);
  }

  if (Offset::PlasmaIonPressure.active) {
    if (nIonFluids!=1) for (int i=0;i<nIonFluids;i++) fprintf(fout,", \"%s[%i]\"",Offset::PlasmaIonPressure.VarList,i);
    else fprintf(fout,", \"%s\"",Offset::PlasmaIonPressure.VarList);
  }

  if (Offset::PlasmaDivU.active) fprintf(fout,", %s",Offset::PlasmaDivU.VarList);  

  if (Offset::PlasmaElectronPressure.active) fprintf(fout,", %s",Offset::PlasmaElectronPressure.VarList);

  if (Offset::MagneticField.active) fprintf(fout,", %s",Offset::MagneticField.VarList);
  if (Offset::ElectricField.active) fprintf(fout,", %s",Offset::ElectricField.VarList);
  if (Offset::Current.active) fprintf(fout,", %s",Offset::Current.VarList);

  if (Offset::b_dot_grad_b.active) fprintf(fout,", %s",Offset::b_dot_grad_b.VarList);
  if (Offset::vE_dot_grad_b.active) fprintf(fout,", %s",Offset::vE_dot_grad_b.VarList);
  if (Offset::b_dot_grad_vE.active) fprintf(fout,", %s",Offset::b_dot_grad_vE.VarList);
  if (Offset::vE_dot_grad_vE.active) fprintf(fout,", %s",Offset::vE_dot_grad_vE.VarList);
  if (Offset::grad_kappaB.active) fprintf(fout,", %s",Offset::grad_kappaB.VarList);

  if (Offset::MagneticFieldGradient.active) fprintf(fout,", %s",Offset::MagneticFieldGradient.VarList);
  if (Offset::MagneticFluxFunction.active) fprintf(fout,", %s",Offset::MagneticFluxFunction.VarList);
}

void PIC::CPLR::DATAFILE::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) {
  int iData,iCoefficient;
  double c,data[nTotalBackgroundVariables];
  double *offset;

  //init the interpolation buffer
  for (iData=0;iData<nTotalBackgroundVariables;iData++) data[iData]=0.0;

  //interpolation loop
  for (iCoefficient=0;iCoefficient<nInterpolationCoeficients;iCoefficient++) {
    c=InterpolationCoeficients[iCoefficient];
    offset=(double*)(MULTIFILE::CurrDataFileOffset+CenterNodeAssociatedDataOffsetBegin+InterpolationList[iCoefficient]->GetAssociatedDataBufferPointer());

    for (iData=0;iData<nTotalBackgroundVariables;iData++) data[iData]+=c*offset[iData];
  }

  //copy the interpolated data
  memcpy(MULTIFILE::CurrDataFileOffset+CenterNodeAssociatedDataOffsetBegin+CenterNode->GetAssociatedDataBufferPointer(),data,nTotalBackgroundVariables*sizeof(double));
}

void PIC::CPLR::DATAFILE::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double data[nTotalBackgroundVariables];
  int i;

  bool gather_output_data=false;

  if (pipe==NULL) gather_output_data=true;
  else if (pipe->ThisThread==CenterNodeThread) gather_output_data=true;

  //get the interpolated data
  if (gather_output_data==true) { // (pipe->ThisThread==CenterNodeThread) {
    memcpy(data,MULTIFILE::CurrDataFileOffset+CenterNodeAssociatedDataOffsetBegin+CenterNode->GetAssociatedDataBufferPointer(),nTotalBackgroundVariables*sizeof(double));
  }

  //send the data to the root processor if needed and print them into a file
  if ((PIC::ThisThread==0)||(pipe==NULL)) {
    if ((CenterNodeThread!=0)&&(pipe!=NULL)) pipe->recv(data,nTotalBackgroundVariables,CenterNodeThread);

    for (i=0;i<nTotalBackgroundVariables;i++) fprintf(fout,"%e ",data[i]);
  }
  else pipe->send(data,nTotalBackgroundVariables);
}

//====================================================
//print the ion flux at a sphere
void PIC::CPLR::DATAFILE::PrintSphereSurfaceIonFlux(char const* fname,double SphereRadius) {
  FILE *fout=NULL;
  int nZenithAngle,nPolarAngle,i,j,k;
  double ZenithAngle,PolarAngle,x[3],n,v[3],flux;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->rootTree;
  CMPI_channel pipe(1000000);

  const int nTotalZenithPoints=100;
  const double dZenithAngle=Pi/(nTotalZenithPoints-1);
  const int nTotalPolarPoints=2*nTotalZenithPoints;
  const double dPolarAngle=2.0*Pi/(nTotalPolarPoints-1);

  if (PIC::ThisThread==0) {
    pipe.openRecvAll();

    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"Lon\", \"Lat\", \"Ion Flux\"\n");
    fprintf(fout,"ZONE I=%i, J=%i, DATAPACKING=POINT\n",nTotalPolarPoints,nTotalZenithPoints);
  }
  else pipe.openSend(0);

  //calculate the flux distribution
  for (nZenithAngle=0,ZenithAngle=-Pi/2.0;nZenithAngle<nTotalZenithPoints;nZenithAngle++,ZenithAngle+=dZenithAngle) {
    for (nPolarAngle=0,PolarAngle=-Pi;nPolarAngle<nTotalPolarPoints;nPolarAngle++,PolarAngle+=dPolarAngle) {
      x[0]=SphereRadius*cos(PolarAngle)*cos(ZenithAngle);
      x[1]=SphereRadius*sin(PolarAngle)*cos(ZenithAngle);
      x[2]=SphereRadius*sin(ZenithAngle);

      node=PIC::Mesh::mesh->findTreeNode(x,node);
      if (PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find cell");

      if (node->Thread==PIC::ThisThread) {
        PIC::CPLR::InitInterpolationStencil(x,node);
        n=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
        PIC::CPLR::GetBackgroundPlasmaVelocity(v);

        flux=-n*(v[0]*x[0]+v[1]*x[1]+v[2]*x[2]);

        if (PIC::ThisThread!=0) pipe.send(flux);
      }

      if (PIC::ThisThread==0) {
        if (node->Thread!=0) flux=pipe.recv<double>(node->Thread);

        fprintf(fout,"%e  %e  %e\n",PolarAngle/Pi*180.0,ZenithAngle/Pi*180.0,flux);
      }

    }
  }


  if (PIC::ThisThread==0) {
    pipe.closeRecvAll();
    fclose(fout);
  }
  else pipe.closeSend();

}


//====================================================
//evaluate the ion flux at the surface
void PIC::CPLR::DATAFILE::EvaluateSurfaceIonFlux(double ShiftFactor) {
  int nTotalElements,el;
  list<cInternalBoundaryConditionsDescriptor>::iterator d;
  cInternalSphericalData* Sphere;
  double x[3],norm[3],l,n,v[3],flux;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->rootTree;
  int i,j,k;


#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_
  for (d=PIC::Mesh::mesh->InternalBoundaryList.begin();d!=PIC::Mesh::mesh->InternalBoundaryList.end();d++) {
    if (d->BondaryType==_INTERNAL_BOUNDARY_TYPE_SPHERE_) {
      Sphere=(cInternalSphericalData*)d->BoundaryElement;
      nTotalElements=Sphere->GetTotalSurfaceElementsNumber();

      //assemble the table of ion fluxes
      double FluxTable[nTotalElements],TotalSurfaceFlux;
      CMPI_channel pipe(1000000);

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else pipe.openSend(0);

      for (el=0;el<nTotalElements;el++) {
        Sphere->GetSurfaceElementNormal(norm,el);
        Sphere->GetSurfaceElementMiddlePoint(x,el);

        l=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
        x[0]+=max(ShiftFactor-1.0,0.0)*l*norm[0],x[1]+=max(ShiftFactor-1.0,0.0)*l*norm[1],x[2]+=max(ShiftFactor-1.0,0.0)*l*norm[2];

        node=PIC::Mesh::mesh->findTreeNode(x,node);
        if (PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find cell");

        if (node->Thread==PIC::ThisThread) {
          PIC::CPLR::InitInterpolationStencil(x,node);
          n=PIC::CPLR::GetBackgroundPlasmaNumberDensity();
          PIC::CPLR::GetBackgroundPlasmaVelocity(v);

          flux=-n*(v[0]*x[0]+v[1]*x[1]+v[2]*x[2]);

          if (PIC::ThisThread!=0) pipe.send(flux);
        }

        if (PIC::ThisThread==0) {
          if (node->Thread!=0) flux=pipe.recv<double>(node->Thread);

          FluxTable[el]=flux;
        }
      }

      if (PIC::ThisThread==0) {
        pipe.closeRecvAll();
      }
      else pipe.closeSend();

      //broundcast the table across all processors
      MPI_Bcast(FluxTable,nTotalElements,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

      //save the table on the sphere
      for (el=0,TotalSurfaceFlux=0.0;el<nTotalElements;el++) {
        Sphere->SolarWindSurfaceFlux[el]=FluxTable[el];
        if (FluxTable[el]>0.0) TotalSurfaceFlux+=FluxTable[el]*Sphere->SurfaceElementArea[el];
      }

      Sphere->TotalSolarWindSurfaceFlux=TotalSurfaceFlux;

    }
    else {
      exit(__LINE__,__FILE__,"Error: unknown internal boundary type");
    }
  }
#endif
}

//====================================================
//compare the extracted data with the referece data set (extracted previously)
void PIC::CPLR::DATAFILE::SaveTestReferenceData(const char *fName,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  static CMPI_channel pipe;
  static FILE *fout=NULL;

  if (startNode==PIC::Mesh::mesh->rootTree) {
    char fname[400];
    sprintf(fname,"%s",fName);

    pipe.init(1000000);

    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.openRecvAll();
      fout=fopen(fname,"w");
    }
    else pipe.openSend(0);
  }

  //loop through all points
  //create the list of the points
  //perform the interpolation loop
  int i,j,k,nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char *offset;

  char SendCellFlag;

  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((PIC::ThisThread==0)||(startNode->Thread==PIC::ThisThread)) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
      SendCellFlag=true;

      //determine whether the cell data neede to be saved
      if (startNode->Thread==PIC::ThisThread) {
        //locate the cell
        if (startNode->block==NULL) SendCellFlag=false,offset=NULL;

        nd=_getCenterNodeLocalNumber(i,j,k);
        if (SendCellFlag==true) {
          if ((CenterNode=startNode->block->GetCenterNode(nd))!=NULL) offset=CenterNode->GetAssociatedDataBufferPointer();
          else SendCellFlag=false,offset=NULL;
        }

        if (startNode->Thread!=0) pipe.send(SendCellFlag);
      }
      else {
        pipe.recv(SendCellFlag,startNode->Thread);
      }

      //save the cell data
      if (SendCellFlag==true) {
        if (startNode->Thread==PIC::ThisThread) {
          if (startNode->Thread==0) {
            for (int iVar=0;iVar<nTotalBackgroundVariables;iVar++) fprintf(fout,"%e ",((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset))[iVar]);
            fprintf(fout,"\n");

            if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)            {
              PIC::Debugger::CatchOutLimitValue((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset),nTotalBackgroundVariables,__LINE__,__FILE__);
            }
          }
          else {
            pipe.send((double*)(offset+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset),nTotalBackgroundVariables);
          }
        }
        else {
          double data[nTotalBackgroundVariables];

          pipe.recv(data,nTotalBackgroundVariables,startNode->Thread);
          for (int iVar=0;iVar<nTotalBackgroundVariables;iVar++) fprintf(fout,"%e ",data[iVar]);
          fprintf(fout,"\n");

          if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)           {
            PIC::Debugger::CatchOutLimitValue(data,nTotalBackgroundVariables,__LINE__,__FILE__);
          }
        }
      }
    }
  }
  else {
    for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) SaveTestReferenceData(NULL,startNode->downNode[nDownNode]);
  }

  if (startNode==PIC::Mesh::mesh->rootTree) {
    if (PIC::Mesh::mesh->ThisThread==0) {
      pipe.closeRecvAll();
      fclose(fout);
    }
    else pipe.closeSend();

    pipe.remove();
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
}



//====================================================================
// the function to produce additional data based on the imported parameters

//generate magnetic field gradient in the given block
void PIC::CPLR::DATAFILE::GenerateMagneticFieldGradient(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {

  if (DIM < 3) exit(__LINE__,__FILE__,"This function is tested for 3D case, may require further testing and development for lower dimensional case!");

  // variables to access data storage
  int nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char* offset;

  //block's mesh parameters
  double  *xNodeMin=node->xmin;
  double  *xNodeMax=node->xmax;
  double  dXCell[3]= {(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_,(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_,(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_};

  //  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  //  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  //  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  const int iMin=0,iMax=_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=_BLOCK_CELLS_Z_-1;


  double xCenter[3]={0.0,0.0,0.0},x[3]={0.0,0.0,0.0};   //locations
  double Bplus[3]={0.0,0.0,0.0}, Bminus[3]={0.0,0.0,0.0}; //values of the field

  for (int k=kMin;k<=kMax;k++) for (int j=jMin;j<=jMax;j++) for (int i=iMin;i<=iMax;i++) {
    //cell center
    xCenter[0]=xNodeMin[0]+dXCell[0]*(0.5+i);
    xCenter[1]=xNodeMin[1]+dXCell[1]*(0.5+j);
    xCenter[2]=xNodeMin[2]+dXCell[2]*(0.5+k);

    //get data storage location
    nd=_getCenterNodeLocalNumber(i,j,k);
    if ((CenterNode=node->block->GetCenterNode(nd))==NULL) continue;

    offset=CenterNode->GetAssociatedDataBufferPointer()+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;

    // compute components of the gradient
    for (int idim=0; idim < DIM; idim++) {
      x[0]=xCenter[0],x[1]=xCenter[1],x[2]=xCenter[2];

      // value of B on face in -idim direction
      x[idim] -= dXCell[idim]/2.0;
      PIC::CPLR::InitInterpolationStencil(x,node);
      PIC::CPLR::GetBackgroundMagneticField(Bminus);

      // value of B on face in +idim direction
      x[idim] += dXCell[idim];
      PIC::CPLR::InitInterpolationStencil(x,node);
      PIC::CPLR::GetBackgroundMagneticField(Bplus);

      //compute and write gradient's components
      *(0+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[0]-Bminus[0]) / dXCell[idim];
      *(3+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[1]-Bminus[1]) / dXCell[idim];
      *(6+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[2]-Bminus[2]) / dXCell[idim];
    }

    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
      PIC::Debugger::CatchOutLimitValue((double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset),9,__LINE__,__FILE__);
    }
  }
}


void PIC::CPLR::DATAFILE::GenerateVarForRelativisticGCA(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {

  // variables to access data storage
  int nd;
  PIC::Mesh::cDataCenterNode *CenterNode;
  char* offset;

  //block's mesh parameters
  double  *xNodeMin=node->xmin;
  double  *xNodeMax=node->xmax;
  double  dXCell[3]= {(xNodeMax[0]-xNodeMin[0])/_BLOCK_CELLS_X_,(xNodeMax[1]-xNodeMin[1])/_BLOCK_CELLS_Y_,(xNodeMax[2]-xNodeMin[2])/_BLOCK_CELLS_Z_};

  //  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  //  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  //  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  const int iMin=0,iMax=_BLOCK_CELLS_X_-1;
  const int jMin=0,jMax=_BLOCK_CELLS_Y_-1;
  const int kMin=0,kMax=_BLOCK_CELLS_Z_-1;


  double xCenter[3]={0.0,0.0,0.0},x[3]={0.0,0.0,0.0};   //locations
  double Bplus[3]={0.0,0.0,0.0}, Bminus[3]={0.0,0.0,0.0}; //values of the field
  double Eplus[3]={0.0,0.0,0.0}, Eminus[3]={0.0,0.0,0.0}; //values of the field
  for (int k=kMin;k<=kMax;k++) for (int j=jMin;j<=jMax;j++) for (int i=iMin;i<=iMax;i++) {
    double grad_b[9], grad_vE[9], grad_kappaB[3];
    //cell center
    xCenter[0]=xNodeMin[0]+dXCell[0]*(0.5+i);
    xCenter[1]=xNodeMin[1]+dXCell[1]*(0.5+j);
    xCenter[2]=xNodeMin[2]+dXCell[2]*(0.5+k);

    //get data storage location
    nd=_getCenterNodeLocalNumber(i,j,k);
    if ((CenterNode=node->block->GetCenterNode(nd))==NULL) continue;

    offset=CenterNode->GetAssociatedDataBufferPointer()+CenterNodeAssociatedDataOffsetBegin+MULTIFILE::CurrDataFileOffset;

    // compute components of the gradient
    for (int idim=0; idim < DIM; idim++) {
      x[0]=xCenter[0],x[1]=xCenter[1],x[2]=xCenter[2];

      // value of B on face in -idim direction
      x[idim] -= dXCell[idim]*0.5;
      PIC::CPLR::InitInterpolationStencil(x,node);
      PIC::CPLR::GetBackgroundMagneticField(Bminus);
      PIC::CPLR::GetBackgroundElectricField(Eminus);

      // value of B on face in +idim direction
      x[idim] += dXCell[idim]*1.0;
      PIC::CPLR::InitInterpolationStencil(x,node);
      PIC::CPLR::GetBackgroundMagneticField(Bplus);
      PIC::CPLR::GetBackgroundElectricField(Eplus);

      double Bminus_norm = Vector3D::Length(Bminus);
      double Bplus_norm = Vector3D::Length(Bplus);
      if (Bminus_norm==0) Bminus_norm +=1e-15;
      if (Bplus_norm==0) Bplus_norm +=1e-15;
      double vE_plus[3], vE_minus[3];
      Vector3D::CrossProduct(vE_plus, Eplus, Bplus);
      Vector3D::CrossProduct(vE_minus, Eminus, Bminus);
      double B2plus = Bplus_norm *Bplus_norm;
      double B2minus= Bminus_norm*Bminus_norm;

      for (int ii=0; ii<3; ii++){
        vE_plus[ii] /= (B2plus+1e-15);
        vE_minus[ii] /= (B2minus+1e-15);
      }
      double vE_plus_norm = Vector3D::Length(vE_plus);
      double vE_minus_norm =Vector3D::Length(vE_minus);
      double c2 =  SpeedOfLight* SpeedOfLight;
      double kappa_plus = sqrt(1- vE_plus_norm*vE_plus_norm/c2 );
      double kappa_minus = sqrt(1- vE_minus_norm*vE_minus_norm/c2 );

      //compute and write gradient's components
      /*
       *(0+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[0]-Bminus[0]) / dXCell[idim];
       *(3+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[1]-Bminus[1]) / dXCell[idim];
       *(6+idim+(double*)(offset+PIC::CPLR::DATAFILE::Offset::MagneticFieldGradient.RelativeOffset))=(Bplus[2]-Bminus[2]) / dXCell[idim];
       */
      grad_b[0+idim] = (Bplus[0]/Bplus_norm-Bminus[0]/Bminus_norm) / dXCell[idim];
      grad_b[3+idim] = (Bplus[1]/Bplus_norm-Bminus[1]/Bminus_norm) / dXCell[idim];
      grad_b[6+idim] = (Bplus[2]/Bplus_norm-Bminus[2]/Bminus_norm) / dXCell[idim];

      grad_vE[0+idim] =  (vE_plus[0]-vE_minus[0])/ dXCell[idim];
      grad_vE[3+idim] =  (vE_plus[1]-vE_minus[1])/ dXCell[idim];
      grad_vE[6+idim] =  (vE_plus[2]-vE_minus[2])/ dXCell[idim];

      grad_kappaB[idim] = (kappa_plus*Bplus_norm - kappa_minus*Bminus_norm) / dXCell[idim];

    }

    double B[3],E[3];
    PIC::CPLR::InitInterpolationStencil(xCenter,node);
    PIC::CPLR::GetBackgroundMagneticField(B);
    PIC::CPLR::GetBackgroundElectricField(E);

    double b_norm = Vector3D::Length(B);
    double vE[3];
    Vector3D::CrossProduct(vE,E,B);


    for (int ii=0;ii<3; ii++) {
      B[ii] /= b_norm+1e-15; //get unit vector of b
      vE[ii] /= b_norm*b_norm+1e-15;
    }
    //(bx*d/dx  + by* d/dy+ bz*b/dz) bx
    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.RelativeOffset)) = B[0]*grad_b[0]+B[1]*grad_b[1]+B[2]*grad_b[2];
    //(bx*d/dx  + by* d/dy+ bz*b/dz) by    
    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.RelativeOffset)) = B[0]*grad_b[3]+B[1]*grad_b[4]+B[2]*grad_b[5];
    //(bx*d/dx  + by* d/dy+ bz*b/dz) bz
    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.RelativeOffset)) = B[0]*grad_b[6]+B[1]*grad_b[7]+B[2]*grad_b[8];
    /*
	if (isnan(*(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.RelativeOffset)))) {
	  printf("b_dot_grad_b is nan\n");
	  printf("xCenter:%e,%e,%e\n", xCenter[0],xCenter[1],xCenter[2]);
	  printf("B:%e,%e,%e, grad_b:%e,%e,%e,%e,%e,%e,%e,%e,%e\n",B[0], B[1],B[2], grad_b[0],grad_b[1],grad_b[2],
		 grad_b[3],grad_b[4],grad_b[5],grad_b[6],grad_b[7],grad_b[8] );
	}
     */
    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) bx
    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b.RelativeOffset)) = vE[0]*grad_b[0]+vE[1]*grad_b[1]+vE[2]*grad_b[2];
    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) by  
    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b.RelativeOffset)) = vE[0]*grad_b[3]+vE[1]*grad_b[4]+vE[2]*grad_b[5];
    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) bz
    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b.RelativeOffset)) = vE[0]*grad_b[6]+vE[1]*grad_b[7]+vE[2]*grad_b[8];

    /*
	if (isnan(*(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_b.RelativeOffset)))) {
	  printf("vE_dot_grad_b is nan\n");
	  printf("xCenter:%e,%e,%e\n", xCenter[0],xCenter[1],xCenter[2]);
	}
     */


    //(bx*d/dx  + by* d/dy+ bz*b/dz) vEx
    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE.RelativeOffset)) = B[0]*grad_vE[0]+B[1]*grad_vE[1]+B[2]*grad_vE[2];
    //(bx*d/dx  + by* d/dy+ bz*b/dz) vEy    
    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE.RelativeOffset)) = B[0]*grad_vE[3]+B[1]*grad_vE[4]+B[2]*grad_vE[5];
    //(bx*d/dx  + by* d/dy+ bz*b/dz) vEz
    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE.RelativeOffset)) = B[0]*grad_vE[6]+B[1]*grad_vE[7]+B[2]*grad_vE[8];

    /*
	if (isnan(*(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_vE.RelativeOffset)))) {
	  printf("b_dot_grad_vE is nan\n");
	  printf("xCenter:%e,%e,%e\n", xCenter[0],xCenter[1],xCenter[2]);
	}
     */


    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) vEx
    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE.RelativeOffset)) = vE[0]*grad_vE[0]+vE[1]*grad_vE[1]+vE[2]*grad_vE[2];
    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) vEy    
    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE.RelativeOffset)) = vE[0]*grad_vE[3]+vE[1]*grad_vE[4]+vE[2]*grad_vE[5];
    //(vEx*d/dx  + vEy* d/dy+ vEz*b/dz) vEz
    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE.RelativeOffset)) = vE[0]*grad_vE[6]+vE[1]*grad_vE[7]+vE[2]*grad_vE[8];

    /*
	if (isnan(*(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::vE_dot_grad_vE.RelativeOffset)))) {
	  printf("vE_dot_grad_vE is nan\n");
	  printf("xCenter:%e,%e,%e\n", xCenter[0],xCenter[1],xCenter[2]);
	}
     */


    //grad_kappaB
    *(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::grad_kappaB.RelativeOffset)) = grad_kappaB[0];
    *(1+(double*)(offset+PIC::CPLR::DATAFILE::Offset::grad_kappaB.RelativeOffset)) = grad_kappaB[1];
    *(2+(double*)(offset+PIC::CPLR::DATAFILE::Offset::grad_kappaB.RelativeOffset)) = grad_kappaB[2];

    /*
	if (isnan(*(0+(double*)(offset+PIC::CPLR::DATAFILE::Offset::grad_kappaB.RelativeOffset)))) {
	  printf("grad_kappaB is nan\n");
	  printf("xCenter:%e,%e,%e\n", xCenter[0],xCenter[1],xCenter[2]);
	}
     */

    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
      PIC::Debugger::CatchOutLimitValue((double*)(offset+PIC::CPLR::DATAFILE::Offset::b_dot_grad_b.RelativeOffset),15,__LINE__,__FILE__);
    }
  }
}


//=====================================================================================================
//calculate drift velocity
//calculate a particle drift velocity (Elkington-2002-JASTP)
void PIC::CPLR::GetDriftVelocity(double *vDrift,double *ParticleVelocity,double ParticleMass,double ParticleCharge,double Time) {
  double E[3],B[3],absB2,absB,absB4,t[3],c;
  double M,gamma,gradAbsB_perp[3],ParticleMomentum[3],ParticleMomentum_normB[3],pParallel;
  int idim;

  //get the particle momentum
  gamma=1.0/sqrt(1.0-sqrt(ParticleVelocity[0]*ParticleVelocity[0]+ParticleVelocity[1]*ParticleVelocity[1]+ParticleVelocity[2]*ParticleVelocity[2])/SpeedOfLight);
  for (idim=0;idim<3;idim++) ParticleMomentum[idim]=gamma*ParticleMass*ParticleVelocity[idim],vDrift[idim]=0.0;

  //get the background fields
  PIC::CPLR::GetBackgroundMagneticField(B);
  PIC::CPLR::GetBackgroundElectricField(E);
  absB2=B[0]*B[0]+B[1]*B[1]+B[2]*B[2];
  absB=sqrt(absB2);
  absB4=absB2*absB2;

  //E cross B drift (avarage the drift velocities directly)
  PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
  int iStencil;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable+omp_get_thread_num(),sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
#else
  memcpy(&Stencil,PIC::InterpolationRoutines::CellCentered::StencilTable,sizeof(PIC::InterpolationRoutines::CellCentered::cStencil));
#endif

  for (iStencil=0;iStencil<Stencil.Length;iStencil++) {
    double cellB[3],cellE[3];

    //loop through all cells of the stencil
    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      SWMF::GetBackgroundElectricField(cellE,Stencil.cell[iStencil]);
      SWMF::GetBackgroundMagneticField(cellB,Stencil.cell[iStencil]);

      break;
    case _PIC_COUPLER_MODE__T96_: case _PIC_COUPLER_MODE__T05_: case _PIC_COUPLER_MODE__KMAG_: 
      for (int i=0;i<3;i++) cellE[i]=0.0;
      DATAFILE::GetBackgroundData(cellB,3,DATAFILE::Offset::MagneticField.RelativeOffset,Stencil.cell[iStencil]);

      break;
    case _PIC_COUPLER_MODE__DATAFILE_:
      DATAFILE::GetBackgroundElectricField(cellE,Stencil.cell[iStencil],Time);
      DATAFILE::GetBackgroundMagneticField(cellB,Stencil.cell[iStencil],Time);

      break;
    default:
      exit(__LINE__,__FILE__,"not implemented");
    }

    Vector3D::CrossProduct(t,cellE,cellB);
    c=Stencil.Weight[iStencil]/(cellB[0]*cellB[0]+cellB[1]*cellB[1]+cellB[2]*cellB[2]);

    for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];
  }

  //next: Grad B drift
  memcpy(ParticleMomentum_normB,ParticleMomentum,3*sizeof(double));
  Vector3D::Orthogonalize(B,ParticleMomentum_normB);
  M=pow(Vector3D::Length(ParticleMomentum_normB),2)/(2.0*ParticleMass*absB);

  GetAbsBackgroundMagneticFieldGradient(gradAbsB_perp,Time);
  Vector3D::Orthogonalize(B,gradAbsB_perp);
  Vector3D::CrossProduct(t,B,gradAbsB_perp);

  c=M/(ParticleCharge*gamma*absB2);
  for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];

  //next drift coeffecient: curvature drift
  double gradB[9],t1[3];

  // structure of gradB is the following
  //   gradB[0:2] = {d/dx, d/dy, d/dz} B_x
  //   gradB[3:5] = {d/dx, d/dy, d/dz} B_y
  //   gradB[6:8] = {d/dx, d/dy, d/dz} B_z
  GetBackgroundMagneticFieldGradient(gradB,Time);

  //t1=(b\dot)b
  t1[0]=(B[0]*gradB[0]+B[1]*gradB[1]+B[2]*gradB[2])/absB4;
  t1[1]=(B[0]*gradB[3]+B[1]*gradB[4]+B[2]*gradB[5])/absB4;
  t1[2]=(B[0]*gradB[6]+B[1]*gradB[7]+B[2]*gradB[8])/absB4;

  Vector3D::CrossProduct(t,B,t1);
  pParallel=Vector3D::ParallelComponentLength(ParticleMomentum,B);

  c=pow(pParallel,2)/(ParticleCharge*gamma*ParticleMass);
  for (idim=0;idim<3;idim++) vDrift[idim]+=c*t[idim];

  //if any of the velocity components are not finite -> set the entire vector to zero (can be caused by magnetic field be zero)
  if ((isfinite(vDrift[0])==false)||(isfinite(vDrift[1])==false)||(isfinite(vDrift[2])==false)) {
    for (idim=0;idim<3;idim++) vDrift[idim]=0.0;
  }

}

//=========================================================================================================================
//COPY dataset 'CurrDataFileOffset' to the 'NextDataFileOffset
void PIC::CPLR::DATAFILE::MULTIFILE::CopyCurrDataFile2NextDataFile(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  if (startNode==NULL) startNode=PIC::Mesh::mesh->rootTree;

  const int iMin=-_GHOST_CELLS_X_,iMax=_GHOST_CELLS_X_+_BLOCK_CELLS_X_-1;
  const int jMin=-_GHOST_CELLS_Y_,jMax=_GHOST_CELLS_Y_+_BLOCK_CELLS_Y_-1;
  const int kMin=-_GHOST_CELLS_Z_,kMax=_GHOST_CELLS_Z_+_BLOCK_CELLS_Z_-1;

  if ((CurrDataFileOffset!=-1)&&(NextDataFileOffset!=-1)&&(CurrDataFileOffset!=NextDataFileOffset)) {
    if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      int i,j,k,nd;
      char *offset;
      PIC::Mesh::cDataCenterNode *cell;

      if (startNode->block!=NULL) for (k=kMin;k<=kMax;k++) for (j=jMin;j<=jMax;j++) for (i=iMin;i<=iMax;i++) {
        nd=_getCenterNodeLocalNumber(i,j,k);
        cell=startNode->block->GetCenterNode(nd);

        if (cell!=NULL) {
          offset=cell->GetAssociatedDataBufferPointer();

          if (offset!=NULL) {
            memcpy(offset+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+NextDataFileOffset,offset+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+CurrDataFileOffset,nTotalBackgroundVariables*sizeof(double));
          }
        }
      }
    }
    else for (int nDownNode=0;nDownNode<(1<<3);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) CopyCurrDataFile2NextDataFile(startNode->downNode[nDownNode]);
  }
}




