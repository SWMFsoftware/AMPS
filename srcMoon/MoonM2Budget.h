#ifndef MOON_M2_BUDGET_H
#define MOON_M2_BUDGET_H
#include <filesystem>
namespace Moon { namespace M2Budget {
// Supporting I19 diagnostics. These functions observe production events;
// they do not advance physics or establish formal I19 qualification.
void Begin(const std::filesystem::path& csv);
void Snapshot(long long completedSteps);
bool Valid();
void End();
} }
#endif
