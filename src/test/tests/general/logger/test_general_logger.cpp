#include "test_harness.h"

#include "logger.h"

#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <sys/types.h>
#include <sys/wait.h>

int ThisThread = 0;
int TotalThreadsNumber = 1;
int deviceThisThread = 0;
int deviceTotalThreadsNumber = 1;
int cudaThreadLimitMallocHeapSize = 0;
unsigned int ExitErrorCode = 0;
MPI_Comm MPI_GLOBAL_COMMUNICATOR = MPI_COMM_WORLD;

long int nint(double value) {
  const long int lower = static_cast<long int>(value);
  return value - lower > 0.5 ? lower + 1 : lower;
}

void PrintErrorLog(const char*) {}
void PrintErrorLog(long int, const char*, const char*) {}
void exit(long int line, const char* file, const char* message) {
  std::string text = AMPS::Test::Location(file, static_cast<int>(line));
  text += message == NULL ? "AMPS trap exit" : message;
  throw AMPS::Test::Failure(text);
}

namespace {

std::string Read(const std::string& path) {
  std::ifstream stream(path.c_str(), std::ios::binary);
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return buffer.str();
}

void Clean(const std::string& path) {
  std::remove(path.c_str());
  for (int generation = 1; generation <= 4; ++generation)
    std::remove((path + "." + std::to_string(generation)).c_str());
}

int CountLines(const std::string& value) {
  int count = 0;
  for (std::size_t index = 0; index < value.size(); ++index)
    if (value[index] == '\n') ++count;
  return count;
}

}  // namespace

AMPS_TEST_CASE("GENERAL-LOGGER-FORMAT-LONG-001",
               "UTC level formatting preserves one complete long message") {
  const std::string path = "amps-test-logger-format.log";
  Clean(path);
  const std::string payload(8192, 'L');
  {
    AMPS::Diagnostics::cFileLogger logger(path);
    logger.Log("INFO", "long-marker=" + payload);
    logger.Flush();
    const std::string observed = Read(path);
    AMPS_EXPECT_EQ(CountLines(observed), 1);
    AMPS_EXPECT_TRUE(observed.find("Z [INFO] long-marker=") != std::string::npos);
    AMPS_EXPECT_TRUE(observed.find(payload) != std::string::npos);
    AMPS_EXPECT_TRUE(!observed.empty() && observed[observed.size() - 1] == '\n');
  }
}

AMPS_TEST_CASE("GENERAL-LOGGER-CONCURRENT-001",
               "concurrent writers produce complete unique lines") {
  const std::string path = "amps-test-logger-concurrent.log";
  Clean(path);
  const int thread_count = 6;
  const int lines_per_thread = 75;
  AMPS::Diagnostics::cFileLogger logger(path);
  std::vector<std::thread> threads;
  for (int thread = 0; thread < thread_count; ++thread) {
    threads.push_back(std::thread([thread, lines_per_thread, &logger]() {
      for (int sequence = 0; sequence < lines_per_thread; ++sequence)
        logger.Log("TRACE", "thread=" + std::to_string(thread) +
                              " sequence=" + std::to_string(sequence));
    }));
  }
  for (std::size_t index = 0; index < threads.size(); ++index)
    threads[index].join();
  logger.Flush();

  const std::string observed = Read(path);
  AMPS_EXPECT_EQ(CountLines(observed), thread_count * lines_per_thread);
  for (int thread = 0; thread < thread_count; ++thread)
    for (int sequence = 0; sequence < lines_per_thread; ++sequence) {
      const std::string token = "thread=" + std::to_string(thread) +
                                " sequence=" + std::to_string(sequence) + "\n";
      AMPS_EXPECT_TRUE(observed.find(token) != std::string::npos);
    }
}

AMPS_TEST_CASE("GENERAL-LOGGER-ROTATE-FLUSH-001",
               "byte-bounded generations rotate and flush before destruction") {
  const std::string path = "amps-test-logger-rotate.log";
  Clean(path);
  {
    AMPS::Diagnostics::cFileLogger logger(path, 220, 2);
    for (int sequence = 0; sequence < 18; ++sequence)
      logger.Log("DEBUG", "rotation-sequence=" + std::to_string(sequence) +
                            " payload=abcdefghijklmnopqrstuvwxyz");
    logger.Flush();
    const std::string active = Read(path);
    const std::string first = Read(path + ".1");
    const std::string second = Read(path + ".2");
    AMPS_EXPECT_TRUE(!active.empty());
    AMPS_EXPECT_TRUE(!first.empty());
    AMPS_EXPECT_TRUE(!second.empty());
    AMPS_EXPECT_TRUE(active.find("rotation-sequence=17") != std::string::npos);
    AMPS_EXPECT_TRUE(active[active.size() - 1] == '\n');
    AMPS_EXPECT_TRUE(first[first.size() - 1] == '\n');
    AMPS_EXPECT_TRUE(second[second.size() - 1] == '\n');
  }
}

AMPS_TEST_CASE("GENERAL-LOGGER-FATAL-001",
               "fatal diagnostics are durable and terminate a subprocess nonzero") {
  const std::string path = "amps-test-logger-fatal.log";
  Clean(path);
  const pid_t child = fork();
  AMPS_EXPECT_TRUE(child >= 0);
  if (child == 0)
    AMPS::Diagnostics::cFileLogger::Fatal(path, "intentional-fatal-marker", 23);

  int status = 0;
  AMPS_EXPECT_EQ(waitpid(child, &status, 0), child);
  AMPS_EXPECT_TRUE(WIFEXITED(status));
  AMPS_EXPECT_EQ(WEXITSTATUS(status), 23);
  const std::string observed = Read(path);
  AMPS_EXPECT_TRUE(observed.find("[FATAL] intentional-fatal-marker") !=
                   std::string::npos);
  AMPS_EXPECT_TRUE(!observed.empty() && observed[observed.size() - 1] == '\n');
}
