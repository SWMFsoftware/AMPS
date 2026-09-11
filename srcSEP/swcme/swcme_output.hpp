#pragma once

// ============================================================================
// swcme_output.hpp
// ----------------------------------------------------------------------------
// Shared checked text-output support for the 1-D and 3-D Tecplot writers.
//
// stdio is commonly buffered, so fprintf/fwrite can appear to succeed even
// when the destination will reject the data during fflush or fclose (the
// canonical example is /dev/full).  A science-output writer must therefore
// check the complete lifecycle: open, every formatted/raw write, flush, stream
// error state, and close.  CheckedTextFile records the first failure and never
// lets a later cleanup error overwrite its more useful byte/row/zone context.
// ============================================================================

#include "swcme_status.hpp"

#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <vector>

namespace swcme {
namespace output {

// FileOperations is a deliberately small C-compatible backend.  Production
// uses the stdio implementation below.  Validation can supply an in-memory
// backend that fails at an exact byte, flush, stream-error query, or close,
// making otherwise rare storage failures deterministic and reproducible.
struct FileOperations {
  void* user_data = nullptr;
  void* (*open)(void* user_data, const char* path, const char* mode) = nullptr;
  std::size_t (*write)(void* user_data, void* handle,
                       const char* bytes, std::size_t count) = nullptr;
  int (*flush)(void* user_data, void* handle) = nullptr;
  int (*error)(void* user_data, void* handle) = nullptr;
  int (*close)(void* user_data, void* handle) = nullptr;
};

namespace detail {

inline void* stdio_open(void*, const char* path, const char* mode) {
  return std::fopen(path,mode);
}

inline std::size_t stdio_write(void*, void* handle, const char* bytes,
                               std::size_t count) {
  return std::fwrite(bytes,1,count,static_cast<std::FILE*>(handle));
}

inline int stdio_flush(void*, void* handle) {
  return std::fflush(static_cast<std::FILE*>(handle));
}

inline int stdio_error(void*, void* handle) {
  return std::ferror(static_cast<std::FILE*>(handle));
}

inline int stdio_close(void*, void* handle) {
  return std::fclose(static_cast<std::FILE*>(handle));
}

}  // namespace detail

inline const FileOperations& stdio_file_operations() noexcept {
  // Function-local static initialization is thread-safe in C++11 and later;
  // the immutable table can therefore be shared by independent output calls.
  static const FileOperations operations{
      nullptr,detail::stdio_open,detail::stdio_write,detail::stdio_flush,
      detail::stdio_error,detail::stdio_close};
  return operations;
}

class CheckedTextFile {
 public:
  explicit CheckedTextFile(const FileOperations& operations) noexcept
      : operations_(operations) {}

  CheckedTextFile(const CheckedTextFile&)=delete;
  CheckedTextFile& operator=(const CheckedTextFile&)=delete;

  ~CheckedTextFile() {
    // Normal callers use finish() so close errors can be reported.  This path
    // exists only for exceptions/early returns and prevents descriptor leaks;
    // a destructor cannot safely surface a new status.
    if (handle_ && operations_.close)
      (void)operations_.close(operations_.user_data,handle_);
  }

  bool open(const char* path) noexcept {
    if (!path || !operations_.open || !operations_.write ||
        !operations_.flush || !operations_.error || !operations_.close)
      return false;
    handle_=operations_.open(operations_.user_data,path,"w");
    return handle_!=nullptr;
  }

  bool good() const noexcept { return failure_.ok(); }

  // Format one record into an owned buffer before passing bytes to the backend.
  // This makes partial writes observable at an exact byte offset, unlike
  // fprintf's single negative return, and avoids forwarding nonliteral format
  // strings to fprintf (which also removes -Wformat-security warnings).
  bool print(const char* failure_context, std::size_t item_index,
             const char* format, ...) noexcept {
    if (!failure_.ok()) return false;
    if (!handle_ || !format) {
      record_failure(failure_context,item_index,bytes_written_);
      return false;
    }

    va_list arguments;
    va_start(arguments,format);
    bool arguments_active=true;
    va_list measure_arguments;
    va_copy(measure_arguments,arguments);
    const int required=std::vsnprintf(nullptr,0,format,measure_arguments);
    va_end(measure_arguments);
    if (required<0) {
      va_end(arguments);
      record_failure(failure_context,item_index,bytes_written_);
      return false;
    }

    try {
      std::vector<char> buffer(static_cast<std::size_t>(required)+1U);
      const int formatted=std::vsnprintf(
          buffer.data(),buffer.size(),format,arguments);
      va_end(arguments);
      arguments_active=false;
      if (formatted!=required) {
        record_failure(failure_context,item_index,bytes_written_);
        return false;
      }

      const std::size_t count=static_cast<std::size_t>(required);
      const std::size_t written=operations_.write(
          operations_.user_data,handle_,buffer.data(),count);
      if (written!=count) {
        // Clamp a malformed backend result before forming the diagnostic; a
        // backend may never report more bytes than it was asked to consume.
        const std::size_t accepted=written<count ? written : 0U;
        record_failure(failure_context,item_index,bytes_written_+accepted);
        bytes_written_+=accepted;
        return false;
      }
      bytes_written_+=written;
      return true;
    } catch (...) {
      if (arguments_active) va_end(arguments);
      record_failure(failure_context,item_index,bytes_written_);
      return false;
    }
  }

  // Complete the buffered-stream lifecycle.  Flush is checked first, then the
  // persistent stream error flag, then close.  Close is always attempted, but
  // record_failure() is first-error-wins so cleanup cannot mask the operation
  // that actually truncated the output.
  swcme::ModelStatus finish(const char* flush_context,
                            const char* stream_context,
                            const char* close_context) noexcept {
    if (!handle_) return failure_.ok()
        ? swcme::ModelStatus::make(
              swcme::StatusCode::FileOpenFailure,"output stream not open")
        : failure_;

    if (failure_.ok() &&
        operations_.flush(operations_.user_data,handle_)!=0)
      record_failure(flush_context,swcme::ModelStatus::npos,bytes_written_);
    if (failure_.ok() &&
        operations_.error(operations_.user_data,handle_)!=0)
      record_failure(stream_context,swcme::ModelStatus::npos,bytes_written_);

    const int close_result=
        operations_.close(operations_.user_data,handle_);
    handle_=nullptr;
    if (failure_.ok() && close_result!=0)
      record_failure(close_context,swcme::ModelStatus::npos,bytes_written_);
    return failure_;
  }

 private:
  void record_failure(const char* context, std::size_t item_index,
                      std::size_t byte_offset) noexcept {
    if (!failure_.ok()) return;
    failure_=swcme::ModelStatus::file_write_failure(
        context,byte_offset,item_index);
  }

  // Copy the small immutable callback table instead of retaining a reference
  // to a caller-owned local object.  The user_data target must remain valid
  // only for the synchronous writer call, while the operation pointers and
  // dispatch metadata themselves cannot dangle inside this object.
  FileOperations operations_;
  void* handle_=nullptr;
  std::size_t bytes_written_=0;
  swcme::ModelStatus failure_=swcme::ModelStatus::success();
};

}  // namespace output
}  // namespace swcme
