//Versioned run-manifest support and strict post-write cache validation.
//
//This module records byte identities and execution state only.  It does not
//interpret a successful process as scientific validation, and it contains no
//physical constants or model decisions.
#ifndef LIOUVILLE_COPS_PROVENANCE_H
#define LIOUVILLE_COPS_PROVENANCE_H

#include <map>
#include <string>
#include <vector>

#include "Input.h"

struct cProvenanceFile {
  std::string Role,Path;
};

struct cRunOutputAudit {
  bool CacheVerified=false;
  std::string Error;
  int KernelVersion=0,nPoints=0,nFaces=0;
  long long KernelRecords=0,KernelBytes=0;
  std::string KernelSetSHA256;
};

//UTC timestamp for state transitions in a manifest.
std::string ProvenanceUTCNow();

//Validate the table/cache correspondence by reading the completed files back.
//Geometry-only runs have no binary cache and are reported as not applicable.
cRunOutputAudit AuditRunOutputs(const std::string& outputDirectory,
    const std::vector<int>& expectedPoints,bool writeKernels,int kernelVersion);

//Write through a temporary sibling and rename, so a reader never observes a
//partially serialized JSON document.  Required input/output files are SHA-256
//hashed with bounded memory.
void WriteRunManifest(const std::string& filename,const std::string& createdUTC,
    const std::string& updatedUTC,const std::string& processState,
    const std::string& runKind,const std::string& inputFile,
    const std::string& executable,const std::string& cspiceVersion,
    const cInput& input,
    const std::vector<cProvenanceFile>& dependencies,
    const std::string& outputDirectory,const cRunOutputAudit& audit,
    const std::string& error="");

#endif
