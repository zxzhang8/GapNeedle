#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <vector>

namespace gapneedle {

using FastaMap = std::unordered_map<std::string, std::string>;

class FastaIndexedReader {
 public:
  explicit FastaIndexedReader(std::string fastaPath);
  ~FastaIndexedReader();
  FastaIndexedReader(FastaIndexedReader&&) noexcept;
  FastaIndexedReader& operator=(FastaIndexedReader&&) noexcept;
  FastaIndexedReader(const FastaIndexedReader&) = delete;
  FastaIndexedReader& operator=(const FastaIndexedReader&) = delete;

  const std::string& fastaPath() const;
  std::vector<std::string> listNames() const;
  int length(const std::string& seqName) const;
  std::string fetch(const std::string& seqName, int start, int end) const;  // [start, end)

 private:
  struct Impl;
  std::string fastaPath_;
  std::unique_ptr<Impl> impl_;
};

FastaMap readFasta(const std::string& path);
FastaMap readFastaSelected(const std::string& path, const std::vector<std::string>& names);
std::vector<std::string> readFastaNames(const std::string& path);
std::vector<std::string> readFastaNamesIndexed(const std::string& path);
std::string readFastaSliceIndexed(const std::string& path,
                                  const std::string& seqName,
                                  int start,
                                  int end);  // [start, end)
void writeFasta(const std::string& path, const FastaMap& records);
std::string reverseComplement(const std::string& seq);

}  // namespace gapneedle
