// Default table used by direct source-tree builds. The AMPS input parser
// replaces this file in the generated build tree.
#ifndef _COMET_CHEMISTRY_GENERATED_H_
#define _COMET_CHEMISTRY_GENERATED_H_

namespace Comet { namespace Chemistry { namespace Generated {
static const std::size_t TermCount=0;
static const std::size_t ReactionCount=0;
// Standard C++ does not permit zero-length arrays.  The two sentinels below
// provide storage for direct source-tree compilation; the zero counts ensure
// that no network operation observes them as chemistry input.
static const ChemicalSolver::cStoichiometricTerm Terms[]={{0,1.0}};
static const ChemicalSolver::cReaction Reactions[]={
    {"none",0,0,0,0,ChemicalSolver::RateLaw::Constant,
        0.0,0.0,0.0,1.0,0.0}};
// Chemistry variants (srcHartley2/CHEMISTRY_VARIANTS.md): no groups or rate
// tables; the sentinels below are never read because the counts are zero.
static const std::size_t GroupCount=0;
static const Comet::Chemistry::cCompiledGroup Groups[]={{"none",0,0,0}};
static const std::size_t RateTableCount=0;
static const Comet::Chemistry::cCompiledRateTable RateTables[]={
    {"none","none",0.0,0.0}};
static const Comet::Chemistry::cCompiledReactionVariant ReactionVariants[]={
    {-1,-1,-1,0}};
}}}

#endif
