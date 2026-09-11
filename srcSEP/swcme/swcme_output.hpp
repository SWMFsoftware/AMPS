#pragma once

// ============================================================================
// swcme_output.hpp
// ----------------------------------------------------------------------------
// Shared checked and transactional text-output support for the 1-D and 3-D
// Tecplot writers.
//
// stdio is commonly buffered, so fprintf/fwrite can appear to succeed even
// when the destination will reject the data during fflush or fclose (the
// canonical example is /dev/full).  A science-output writer must therefore
// check the complete lifecycle: open, every formatted/raw write, flush, stream
// error state, and close.  OUT03 additionally stages regular-file output in
// the destination directory and publishes it only through a final atomic
// rename.  CheckedTextFile records the first failure and never lets a later
// cleanup error overwrite its more useful byte/row/zone context.
// ============================================================================

#include "swcme_status.hpp"

#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

namespace swcme {
namespace output {

// FileOperations is a deliberately small C-compatible backend.  Production
// uses the stdio implementation below.  Validation can supply an in-memory
// backend that fails at an exact byte, flush, stream-error query, close,
// commit, or removal, making otherwise rare storage failures deterministic
// and reproducible.
struct FileOperations {
  void* user_data = nullptr;
  void* (*open)(void* user_data, const char* path, const char* mode) = nullptr;
  std::size_t (*write)(void* user_data, void* handle,
                       const char* bytes, std::size_t count) = nullptr;
  int (*flush)(void* user_data, void* handle) = nullptr;
  int (*error)(void* user_data, void* handle) = nullptr;
  int (*close)(void* user_data, void* handle) = nullptr;
  // commit must atomically replace destination with the already closed
  // temporary file when the host platform supports regular-file replacement.
  int (*commit)(void* user_data, const char* temporary_path,
                const char* destination_path) = nullptr;
  // remove discards an uncommitted temporary file after any earlier failure.
  int (*remove)(void* user_data, const char* temporary_path) = nullptr;
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

inline int stdio_commit(void*,const char* temporary_path,
                        const char* destination_path) {
  // Replacing a directory, device, FIFO, or symbolic link is surprising and
  // can be destructive (for example, replacing /dev/full while validating
  // OUT02 as a privileged user).  Only an absent destination or an existing
  // regular file is eligible for the transactional rename.
  std::error_code status_error;
  const std::filesystem::file_status destination_status=
      std::filesystem::symlink_status(destination_path,status_error);
  if (status_error &&
      status_error!=std::errc::no_such_file_or_directory)
    return -1;
  if (!status_error && std::filesystem::exists(destination_status) &&
      !std::filesystem::is_regular_file(destination_status))
    return -1;

  // The temporary name is constructed by appending to the destination path,
  // so both names are in the same directory.  On the POSIX platforms targeted
  // by SWCME, rename therefore performs one atomic namespace replacement.
  return std::rename(temporary_path,destination_path);
}

inline int stdio_remove(void*,const char* temporary_path) {
  return std::remove(temporary_path);
}

}  // namespace detail

inline const FileOperations& stdio_file_operations() noexcept {
  // Function-local static initialization is thread-safe in C++11 and later;
  // the immutable table can therefore be shared by independent output calls.
  static const FileOperations operations{
      nullptr,detail::stdio_open,detail::stdio_write,detail::stdio_flush,
      detail::stdio_error,detail::stdio_close,detail::stdio_commit,
      detail::stdio_remove};
  return operations;
}

inline std::uint64_t next_temporary_sequence() noexcept {
  // Exclusive creation, rather than this counter alone, supplies collision
  // safety across processes.  The atomic sequence prevents needless retries
  // among concurrent writers in one process and never contributes to output
  // physics or reproducibility.
  static std::atomic<std::uint64_t> sequence{1};
  return sequence.fetch_add(1,std::memory_order_relaxed);
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
    cancel();
  }

  bool open(const char* path) noexcept {
    if (!path || !operations_.open || !operations_.write ||
        !operations_.flush || !operations_.error || !operations_.close)
      return false;
    handle_=operations_.open(operations_.user_data,path,"w");
    return handle_!=nullptr;
  }

  // Open a same-directory staging file with exclusive creation.  Appending a
  // private suffix to the complete destination path guarantees that rename
  // never crosses a filesystem boundary, while "wx" prevents one writer from
  // truncating another writer's staging file after a name collision.
  bool open_transactional(const char* destination_path) noexcept {
    if (!destination_path || !operations_.open || !operations_.write ||
        !operations_.flush || !operations_.error || !operations_.close ||
        !operations_.commit || !operations_.remove)
      return false;

    try {
      destination_path_=destination_path;
      constexpr unsigned max_name_attempts=64;
      for (unsigned attempt=0; attempt<max_name_attempts; ++attempt) {
        temporary_path_=destination_path_+".swcme-tmp-"+
            std::to_string(next_temporary_sequence());
        handle_=operations_.open(
            operations_.user_data,temporary_path_.c_str(),"wx");
        if (handle_) {
          transactional_=true;
          return true;
        }
      }
    } catch (...) {
      // Allocation/path-construction failure is indistinguishable from an
      // inability to create the staging file at the public writer boundary.
    }
    temporary_path_.clear();
    destination_path_.clear();
    return false;
  }

  bool good() const noexcept { return failure_.ok(); }

  // Abandon an output whose surrounding physics operation failed before a
  // complete product existed.  Unlike finish(), cancel never flushes or
  // commits: it closes the handle and removes a transactional staging file so
  // an incomplete header cannot replace a valid previous destination.
  void cancel() noexcept {
    if (handle_ && operations_.close)
      (void)operations_.close(operations_.user_data,handle_);
    handle_=nullptr;
    if (transactional_ && !temporary_path_.empty() && operations_.remove)
      (void)operations_.remove(operations_.user_data,temporary_path_.c_str());
    transactional_=false;
    temporary_path_.clear();
    destination_path_.clear();
  }

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
  // that actually truncated the output.  For a transactional stream, commit
  // follows successful close; any earlier failure instead removes staging.
  swcme::ModelStatus finish(const char* flush_context,
                            const char* stream_context,
                            const char* close_context,
                            const char* commit_context="output commit") noexcept {
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

    if (transactional_) {
      if (failure_.ok() && operations_.commit(
              operations_.user_data,temporary_path_.c_str(),
              destination_path_.c_str())!=0) {
        // A rename failure is not a byte-stream failure: the staged file was
        // complete, but the public destination was never changed.
        failure_=swcme::ModelStatus::file_commit_failure(
            commit_context,bytes_written_);
      }

      if (!failure_.ok()) {
        // Cleanup follows first-error-wins semantics.  A failed removal may
        // leave a private staging file for diagnosis, but cannot replace the
        // earlier write/close/commit status returned to the caller.
        (void)operations_.remove(
            operations_.user_data,temporary_path_.c_str());
      }
      transactional_=false;
      temporary_path_.clear();
      destination_path_.clear();
    }
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
  bool transactional_=false;
  std::string temporary_path_;
  std::string destination_path_;
  swcme::ModelStatus failure_=swcme::ModelStatus::success();
};

}  // namespace output
}  // namespace swcme
