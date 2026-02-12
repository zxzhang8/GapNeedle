#include "minimap2_bridge.h"

#include <stdio.h>
#include <string.h>

#include "minimap.h"

int gn_mm2_align_to_paf(const gn_mm2_request* req, char* err_buf, int err_buf_len) {
  if (!req || !req->target_fasta || !req->query_fasta || !req->output_paf) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len, "invalid minimap2 request");
    }
    return -1;
  }

  // TODO: Complete direct minimap2 C API wiring for single-sequence filtering.
  // Current fallback returns not implemented so behavior is explicit during migration.
  if (err_buf && err_buf_len > 0) {
    snprintf(err_buf, (size_t)err_buf_len,
             "minimap2 bridge wiring pending: use stub or complete bridge implementation");
  }
  return -2;
}
