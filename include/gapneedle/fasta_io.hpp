#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace gapneedle {

using FastaMap = std::unordered_map<std::string, std::string>;

FastaMap readFasta(const std::string& path);
FastaMap readFastaSelected(const std::string& path, const std::vector<std::string>& names);
void writeFasta(const std::string& path, const FastaMap& records);
std::string reverseComplement(const std::string& seq);

}  // namespace gapneedle
