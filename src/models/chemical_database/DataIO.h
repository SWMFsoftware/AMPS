#ifndef AMPS_CHEMICAL_DATABASE_DATA_IO_H
#define AMPS_CHEMICAL_DATABASE_DATA_IO_H
#include <cstddef>
#include <string>
#include <vector>
namespace AMPS { namespace ChemicalDatabase { namespace ReactionData {
using Row=std::vector<std::string>;
using Csv=std::vector<Row>;
// Original production source-authentication service. Hashing authenticates
// bytes, never scientific permission/approval. These bool APIs preserve output
// on failure and remain source-compatible with the historical Moon aliases.
bool FileSha256(const std::string& path,std::string* digest,std::string* error);
bool VerifyFileSha256(const std::string& path,const std::string& expected,std::string* error);
bool DigestSyntax(const std::string& digest);
bool BytesSha256(const std::string& bytes,std::string* digest,std::string* error);
// Strict quoted CSV: exact ordered header, field counts, quotes/newlines;
// no trimming, locale conversion, gap filling or interpretation of science.
bool ParseCsv(const std::string& bytes,const Row& header,Csv* rows,std::string* error);
// Retained PHIDRATES path hashes the SAME fixed-size bytes that it parses.
bool ReadPinnedCsv(const std::string& path,const std::string& expected,std::size_t bytes,
    const Row& header,Csv* rows,std::string* error);
// General runtime read uses the caller's explicit ceiling (0 = no artificial
// ceiling), rejects unrepresentable lengths and never truncates source bytes.
bool ReadFileBytes(const std::string& path,std::size_t maximumBytes,
    std::string* bytes,std::string* error);
bool ValidUtf8(const std::string& bytes);
} } }
#endif
