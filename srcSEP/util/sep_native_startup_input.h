#ifndef SEP_NATIVE_STARTUP_INPUT_H
#define SEP_NATIVE_STARTUP_INPUT_H

namespace SEP { namespace Testing { namespace NativeStartupInput {

// A descriptor carrying this explicit capability requests the existing full
// schema-6 Initialization::LoadFile/Install startup pipeline for its
// --test-input. Other tests use that option for their own reviewed manifests;
// their files must never be guessed to be scientific initialization records.
// The standalone driver additionally requires a single selected test and an
// explicit path. This grants no mutable initialization override, changes no
// production CLI grammar, and does not qualify the resulting physical model.
inline const char* DescriptorTag() {
  return "native-startup-input-schema6";
}

} } }

#endif
