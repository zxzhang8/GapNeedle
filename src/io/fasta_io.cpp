#include "gapneedle/fasta_io.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace gapneedle {

namespace {

std::string normalizeName(const std::string& raw) {
  std::string name = raw;
  auto sp = name.find(' ');
  if (sp != std::string::npos) {
    name = name.substr(0, sp);
  }
  return name;
}

std::string trim(const std::string& s) {
  auto begin = s.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return {};
  }
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(begin, end - begin + 1);
}

}  // namespace

FastaMap readFasta(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open FASTA: " + path);
  }

  FastaMap out;
  std::string line;
  std::string current;
  std::ostringstream seq;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == '>') {
      if (!current.empty()) {
        out[current] = seq.str();
        seq.str("");
        seq.clear();
      }
      current = normalizeName(trim(line.substr(1)));
      continue;
    }
    if (current.empty()) {
      continue;
    }
    for (char ch : trim(line)) {
      if (!std::isspace(static_cast<unsigned char>(ch))) {
        seq << static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
      }
    }
  }
  if (!current.empty()) {
    out[current] = seq.str();
  }
  return out;
}

FastaMap readFastaSelected(const std::string& path, const std::vector<std::string>& names) {
  if (names.empty()) {
    return {};
  }
  FastaMap all = readFasta(path);
  FastaMap out;
  for (const auto& name : names) {
    auto it = all.find(name);
    if (it != all.end()) {
      out.emplace(it->first, it->second);
    }
  }
  return out;
}

void writeFasta(const std::string& path, const FastaMap& records) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Failed to write FASTA: " + path);
  }
  for (const auto& [name, seq] : records) {
    out << '>' << name << '\n';
    for (std::size_t i = 0; i < seq.size(); i += 80) {
      out << seq.substr(i, std::min<std::size_t>(80, seq.size() - i)) << '\n';
    }
  }
}

std::string reverseComplement(const std::string& seq) {
  std::string out;
  out.resize(seq.size());
  for (std::size_t i = 0; i < seq.size(); ++i) {
    char c = static_cast<char>(std::toupper(static_cast<unsigned char>(seq[seq.size() - 1 - i])));
    switch (c) {
      case 'A': out[i] = 'T'; break;
      case 'T': out[i] = 'A'; break;
      case 'C': out[i] = 'G'; break;
      case 'G': out[i] = 'C'; break;
      default: out[i] = 'N'; break;
    }
  }
  return out;
}

}  // namespace gapneedle
