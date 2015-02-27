//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf

#include "pic.h"

//int PIC::nTotalSpecies=0;
int PIC::ThisThread=0,PIC::nTotalThreads=1;

//the list containing the functions used to exchange the run time execution statistics
vector<PIC::fExchangeExecutionStatistics> PIC::ExchangeExecutionStatisticsFunctions;

//the list contains the functions used for user defined sampling procedures
vector<PIC::IndividualModelSampling::fRequestSamplingData> PIC::IndividualModelSampling::RequestSamplingData;
vector<PIC::IndividualModelSampling::fSamplingProcedure> PIC::IndividualModelSampling::SamplingProcedure;
vector<PIC::IndividualModelSampling::fPrintVariableList> PIC::IndividualModelSampling::PrintVariableList;
vector<PIC::IndividualModelSampling::fInterpolateCenterNodeData> PIC::IndividualModelSampling::InterpolateCenterNodeData;
vector<PIC::IndividualModelSampling::fPrintSampledData> PIC::IndividualModelSampling::PrintSampledData;
vector<PIC::IndividualModelSampling::fRequestStaticCellData> PIC::IndividualModelSampling::RequestStaticCellData;

//generic particle transformation
//PIC::ChemicalReactions::GenericParticleTranformation::fTransformationIndicator *PIC::ChemicalReactions::GenericParticleTranformation::TransformationIndicator=NULL;
//PIC::ChemicalReactions::GenericParticleTranformation::fTransformationProcessor *PIC::ChemicalReactions::GenericParticleTranformation::TransformationProcessor=NULL;

//execution alarm
bool PIC::Alarm::AlarmInitialized=false,PIC::Alarm::WallTimeExeedsLimit=false;
double PIC::Alarm::StartTime=0.0,PIC::Alarm::RequestedExecutionWallTime=0.0;

//the file descriptor and the prefix for output of the diagnostic
FILE* PIC::DiagnospticMessageStream=stdout;
char PIC::DiagnospticMessageStreamName[_MAX_STRING_LENGTH_PIC_]="stdout";

//the directory for the output and input files
char PIC::OutputDataFileDirectory[_MAX_STRING_LENGTH_PIC_]=".";
char PIC::InputDataFileDirectory[_MAX_STRING_LENGTH_PIC_]=".";

//define the test-run parameters
bool PIC::ModelTestRun::mode=false;
int PIC::ModelTestRun::nTotalIteraction=-1;

//the name of the resdtart files
char PIC::Restart::SavedSamplingDataRestartFileName[_MAX_STRING_LENGTH_PIC_]="";
