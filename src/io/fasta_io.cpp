#include "gapneedle/fasta_io.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
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

struct FaiEntry {
  std::string name;
  long long length{0};
  long long offset{0};
  int lineBases{0};
  int lineWidth{0};
};

std::string faiPathOf(const std::string& fastaPath) {
  return fastaPath + ".fai";
}

std::unordered_map<std::string, FaiEntry> parseFai(const std::string& fastaPath, std::vector<std::string>* namesOut) {
  std::ifstream in(faiPathOf(fastaPath));
  if (!in) {
    throw std::runtime_error("Failed to open FASTA index (.fai): " + faiPathOf(fastaPath));
  }

  std::unordered_map<std::string, FaiEntry> out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string name, len, off, lbase, lwidth;
    if (!std::getline(ss, name, '\t') ||
        !std::getline(ss, len, '\t') ||
        !std::getline(ss, off, '\t') ||
        !std::getline(ss, lbase, '\t') ||
        !std::getline(ss, lwidth, '\t')) {
      continue;
    }
    FaiEntry e;
    e.name = name;
    e.length = std::stoll(len);
    e.offset = std::stoll(off);
    e.lineBases = std::stoi(lbase);
    e.lineWidth = std::stoi(lwidth);
    if (e.name.empty() || e.length < 0 || e.offset < 0 || e.lineBases <= 0 || e.lineWidth <= 0) {
      continue;
    }
    if (out.emplace(e.name, e).second && namesOut) {
      namesOut->push_back(e.name);
    }
  }
  if (out.empty()) {
    throw std::runtime_error("No valid entries in FASTA index (.fai): " + faiPathOf(fastaPath));
  }
  return out;
}

void buildFai(const std::string& fastaPath) {
  std::ifstream in(fastaPath, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open FASTA: " + fastaPath);
  }
  std::ofstream out(faiPathOf(fastaPath), std::ios::trunc);
  if (!out) {
    throw std::runtime_error("Failed to write FASTA index (.fai): " + faiPathOf(fastaPath));
  }

  std::string line;
  std::string currName;
  long long currLen = 0;
  long long seqOffset = -1;
  int lineBases = 0;
  int lineWidth = 0;

  std::streampos lineStart = in.tellg();
  while (std::getline(in, line)) {
    std::streampos nextPos = in.tellg();
    if (!line.empty() && line.back() == '\r') line.pop_back();

    if (!line.empty() && line[0] == '>') {
      if (!currName.empty()) {
        out << currName << '\t' << currLen << '\t' << seqOffset << '\t'
            << lineBases << '\t' << lineWidth << '\n';
      }
      currName = normalizeName(trim(line.substr(1)));
      currLen = 0;
      seqOffset = -1;
      lineBases = 0;
      lineWidth = 0;
    } else if (!currName.empty()) {
      int bases = 0;
      for (char c : line) {
        if (!std::isspace(static_cast<unsigned char>(c))) ++bases;
      }
      if (bases > 0) {
        if (seqOffset < 0) {
          seqOffset = static_cast<long long>(lineStart);
          lineBases = bases;
          if (nextPos != std::streampos(-1)) {
            lineWidth = static_cast<int>(nextPos - lineStart);
          } else {
            lineWidth = bases;
          }
        }
        currLen += bases;
      }
    }

    if (nextPos == std::streampos(-1)) break;
    lineStart = nextPos;
  }
  if (!currName.empty()) {
    out << currName << '\t' << currLen << '\t' << seqOffset << '\t'
        << lineBases << '\t' << lineWidth << '\n';
  }
}

std::unordered_map<std::string, FaiEntry> loadOrBuildFai(const std::string& fastaPath, std::vector<std::string>* namesOut) {
  try {
    return parseFai(fastaPath, namesOut);
  } catch (...) {
    buildFai(fastaPath);
    return parseFai(fastaPath, namesOut);
  }
}

}  // namespace

struct FastaIndexedReader::Impl {
  std::unordered_map<std::string, FaiEntry> entries;
  std::vector<std::string> names;
};

FastaIndexedReader::FastaIndexedReader(std::string fastaPath) : fastaPath_(std::move(fastaPath)), impl_(std::make_unique<Impl>()) {
  impl_->entries = loadOrBuildFai(fastaPath_, &impl_->names);
}

FastaIndexedReader::~FastaIndexedReader() = default;

FastaIndexedReader::FastaIndexedReader(FastaIndexedReader&& other) noexcept
    : fastaPath_(std::move(other.fastaPath_)), impl_(std::move(other.impl_)) {}

FastaIndexedReader& FastaIndexedReader::operator=(FastaIndexedReader&& other) noexcept {
  if (this == &other) return *this;
  fastaPath_ = std::move(other.fastaPath_);
  impl_ = std::move(other.impl_);
  return *this;
}

const std::string& FastaIndexedReader::fastaPath() const { return fastaPath_; }

std::vector<std::string> FastaIndexedReader::listNames() const {
  return impl_->names;
}

int FastaIndexedReader::length(const std::string& seqName) const {
  auto it = impl_->entries.find(seqName);
  return it == impl_->entries.end() ? -1 : static_cast<int>(it->second.length);
}

std::string FastaIndexedReader::fetch(const std::string& seqName, int start, int end) const {
  auto it = impl_->entries.find(seqName);
  if (it == impl_->entries.end()) {
    throw std::runtime_error("Sequence not found in FASTA index: " + seqName);
  }
  const FaiEntry& e0 = it->second;
  const int seqLen = static_cast<int>(e0.length);
  const int s = std::max(0, start);
  const int e = std::min(end, seqLen);
  if (e <= s) return {};

  std::ifstream in(fastaPath_, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open FASTA: " + fastaPath_);
  }
  const long long baseOffset =
      e0.offset + static_cast<long long>(s / e0.lineBases) * e0.lineWidth + (s % e0.lineBases);
  in.seekg(baseOffset, std::ios::beg);
  if (!in) {
    throw std::runtime_error("Failed to seek FASTA for sequence slice: " + seqName);
  }

  const int want = e - s;
  std::string out;
  out.reserve(static_cast<std::size_t>(want));
  char ch = 0;
  while (static_cast<int>(out.size()) < want && in.get(ch)) {
    if (ch == '\n' || ch == '\r') continue;
    out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
  }
  if (static_cast<int>(out.size()) != want) {
    throw std::runtime_error("Failed to fetch full sequence slice: " + seqName);
  }
  return out;
}

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

std::vector<std::string> readFastaNames(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Failed to open FASTA: " + path);
  }

  std::vector<std::string> names;
  std::unordered_set<std::string> seen;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] != '>') {
      continue;
    }
    std::string name = normalizeName(trim(line.substr(1)));
    if (name.empty()) {
      continue;
    }
    if (seen.insert(name).second) {
      names.push_back(std::move(name));
    }
  }
  return names;
}

std::vector<std::string> readFastaNamesIndexed(const std::string& path) {
  FastaIndexedReader reader(path);
  return reader.listNames();
}

std::string readFastaSliceIndexed(const std::string& path,
                                  const std::string& seqName,
                                  int start,
                                  int end) {
  FastaIndexedReader reader(path);
  return reader.fetch(seqName, start, end);
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
