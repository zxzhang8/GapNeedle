#include "gapneedle/aligner.hpp"

#include "gapneedle/fasta_io.hpp"
#include "minimap2_bridge.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace gapneedle {

namespace {

std::string safePart(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char ch : s) {
    if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || (ch >= '0' && ch <= '9') ||
        ch == '.' || ch == '_' || ch == '-') {
      out.push_back(ch);
    } else {
      out.push_back('_');
    }
  }
  return out;
}

std::string defaultPafPath(const AlignmentRequest& req) {
  namespace fs = std::filesystem;
  const fs::path target(req.targetFasta);
  const fs::path query(req.queryFasta);
  const std::string base = safePart(query.stem().string()) + "." + safePart(req.querySeq) +
                           "_vs_" + safePart(target.stem().string()) + "." + safePart(req.targetSeq);
  const fs::path dir = fs::path("resources") / base;
  fs::create_directories(dir);
  return (dir / (base + "." + safePart(req.preset) + ".paf")).string();
}

std::string defaultIndexCacheDir() {
  return (std::filesystem::path("resources") / "mm2_index").string();
}

void ensureParentDirectoryExists(const std::string& outputPath) {
  const std::filesystem::path pafPath(outputPath);
  if (!pafPath.has_parent_path()) {
    return;
  }

  std::error_code ec;
  std::filesystem::create_directories(pafPath.parent_path(), ec);
  if (ec) {
    throw std::runtime_error(std::string("failed to create output directory: ") + ec.message());
  }
}

}  // namespace

AlignmentResult Minimap2Aligner::align(const AlignmentRequest& request) const {
  const std::string pafPath = request.outputPafPath.empty() ? defaultPafPath(request) : request.outputPafPath;
  ensureParentDirectoryExists(pafPath);

  if (request.reuseExisting && std::filesystem::exists(pafPath)) {
    std::error_code ec;
    const auto sz = std::filesystem::file_size(pafPath, ec);
    if (!ec && sz > 0) {
      return AlignmentResult{pafPath, true, {"reused existing paf"}};
    }
  }

  gn_mm2_request cReq{};
  cReq.target_fasta = request.targetFasta.c_str();
  cReq.query_fasta = request.queryFasta.c_str();
  cReq.target_seq = request.targetSeq.c_str();
  cReq.query_seq = request.querySeq.c_str();
  cReq.reverse_target = request.reverseTarget ? 1 : 0;
  cReq.reverse_query = request.reverseQuery ? 1 : 0;
  cReq.preset = request.preset.c_str();
  cReq.threads = request.threads;
  cReq.allow_index_cache = request.useIndexCache ? 1 : 0;
  const std::string indexCacheDir = request.indexCacheDir.empty() ? defaultIndexCacheDir() : request.indexCacheDir;
  cReq.index_cache_dir = indexCacheDir.c_str();
  cReq.output_paf = pafPath.c_str();

  char errBuf[512] = {0};
  char traceBuf[4096] = {0};
  const int rc = gn_mm2_align_to_paf(&cReq, errBuf, 512, traceBuf, 4096);
  if (rc != 0) {
    throw std::runtime_error(std::string("minimap2 alignment failed: ") + errBuf);
  }

  AlignmentResult out{pafPath, false, {}};
  if (traceBuf[0] != '\0') {
    std::istringstream in(traceBuf);
    for (std::string line; std::getline(in, line);) {
      if (!line.empty()) {
        out.warnings.push_back(line);
      }
    }
  }
  return out;
}

}  // namespace gapneedle
