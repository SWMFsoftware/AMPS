#ifndef _GENERIC_COMET_RUNTIME_INPUT_H_
#define _GENERIC_COMET_RUNTIME_INPUT_H_

#include <string>

namespace Comet {
namespace RuntimeInput {

// A compile-time comet configuration sets Required=true when the executable
// contains capabilities/layout only and must receive its physical run model
// through ``-i FILE``. Legacy Halley and Hartley inputs leave this false and
// retain their generated constants until they are migrated deliberately.
extern bool Required;

// Parse the generic AMPS runtime-input grammar, select a named case through
// the shared AMPS::Input::ResolveCase contract, validate the comet-owned keys,
// and install mutable run state before mesh/block/particle allocation. This
// routine is collective after PIC::Init_BeforeParser(): every rank verifies
// that it resolved byte-identical canonical input before applying it.
bool LoadAndApply(const std::string &path,const std::string &caseName,
    std::string *error=NULL);

const std::string& LoadedPath();
const std::string& SelectedCase();

}
}

#endif
