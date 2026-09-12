#include "test_framework.hpp"

#include "swcme_output.hpp"

#include <cstddef>
#include <cstring>
#include <iostream>
#include <string>

namespace {

// A minimal allocation-owning backend keeps OUT08's runtime assertion focused
// on record semantics rather than the filesystem.  The strict build compiles
// this same translation unit twice, so callback signatures, size conversions,
// and format arguments are checked in both debug and optimized configurations.
struct LiteralSink {
  std::string contents;
  bool opened=false;
  bool closed=false;
};

void* open_sink(void* user_data,const char*,const char*) {
  LiteralSink& sink=*static_cast<LiteralSink*>(user_data);
  sink.opened=true;
  return &sink;
}

std::size_t write_sink(void*,void* handle,const char* data,
                       std::size_t count) {
  LiteralSink& sink=*static_cast<LiteralSink*>(handle);
  sink.contents.append(data,count);
  return count;
}

int flush_sink(void*,void*) { return 0; }
int error_sink(void*,void*) { return 0; }

int close_sink(void*,void* handle) {
  LiteralSink& sink=*static_cast<LiteralSink*>(handle);
  sink.closed=true;
  return 0;
}

swcme::output::FileOperations sink_operations(LiteralSink& sink) {
  swcme::output::FileOperations operations;
  operations.user_data=&sink;
  operations.open=open_sink;
  operations.write=write_sink;
  operations.flush=flush_sink;
  operations.error=error_sink;
  operations.close=close_sink;
  return operations;
}

}  // namespace

// OUT08 combines a compile-time and runtime contract.  The Makefile's strict
// target supplies the compile-time half (-Werror plus format/conversion/shadow
// diagnostics); this registered assertion proves that the literal-only path
// preserves percent tokens verbatim while the annotated formatted path still
// renders correctly through the shared checked lifecycle.
void test_out08(swcme_test::Context& context) {
  std::cout << "OUT08 strict warning writer build\n";

  LiteralSink sink;
  const swcme::output::FileOperations operations=sink_operations(sink);
  swcme::output::CheckedTextFile output(operations);
  context.expect_true(output.open("unused-out08-path"),
                      "strict-writer probe opens its injected backend");

  const bool literal_ok=output.write_literal(
      "OUT08 literal record",swcme::ModelStatus::npos,
      "literal percent tokens: 100% %s %zu\n");
  const std::size_t item_count=9U;
  const bool formatted_ok=output.print(
      "OUT08 formatted record",swcme::ModelStatus::npos,
      "formatted values: %d %zu\n",7,item_count);
  const swcme::ModelStatus status=output.finish(
      "OUT08 flush","OUT08 stream error","OUT08 close");

  context.expect_true(literal_ok && formatted_ok && status.ok() &&
                          sink.opened && sink.closed,
                      "literal and formatted records complete successfully");
  context.expect_true(
      sink.contents==
          "literal percent tokens: 100% %s %zu\n"
          "formatted values: 7 9\n",
      "literal percent tokens are never interpreted as a format string");
}
