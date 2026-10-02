// Model registration is independent of the AMPS mesh/publication adapter.
#ifndef SEP3D_RUNTIME_BACKGROUND_FACTORY_H
#define SEP3D_RUNTIME_BACKGROUND_FACTORY_H
#include "run_configuration.h"
#include "../background/bg_provider.h"
#include <functional>
#include <memory>
namespace SEP3D { namespace RuntimeModel {
using BackgroundFactory=std::function<Core::Status(const RunConfiguration3D&,
    std::shared_ptr<Background::BackgroundProvider>*)>;
// The callback receives the immutable application configuration and creates a
// provider in its output. Model parameters and derivative/heating policies
// belong in the provider's resolved manifest; registration does not define a
// parser or an external-process protocol for model-specific text parameters.
// Register once on EACH rank, before initialization. Duplicate names (including
// built-ins) are rejected. Factories must freeze/validate their configuration;
// their provider owns Prepare/evaluation, never PIC buffers or MPI calls.
Core::Status RegisterBackgroundModel(const std::string& id,BackgroundFactory factory);
// Dispatch analytic-parker/swcme built-ins or the registered runtime-model ID.
// Validate the candidate before assigning output. On any error, *output is
// unchanged. SWMF acquisition is host-owned and does not use this factory.
Core::Status CreateBackgroundProvider(const RunConfiguration3D& configuration,
    std::shared_ptr<Background::BackgroundProvider>* output);
// Compare provenance tags, not provider names or filenames. A runtime-model
// extension must report ProviderKind::RuntimeModel even if its adapter reuses
// another model internally; a SWMF import must remain explicitly identifiable.
bool BackgroundAuthorityMatches(BackgroundAuthority authority,
                                Background::ProviderKind kind);
} }
#endif
