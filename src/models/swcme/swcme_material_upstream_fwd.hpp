#pragma once

// Keep the typed construction boundary in the SWCME source owner. Consumers
// may share an immutable event without importing its implementation or
// opening a canonical SWCME namespace in an application source file.
namespace swcme { namespace material1d {
class Event;
} }
