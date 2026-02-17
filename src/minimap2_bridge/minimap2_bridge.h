#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  const char* target_fasta;
  const char* query_fasta;
  const char* target_seq;
  const char* query_seq;
  int reverse_target;
  int reverse_query;
  const char* preset;
  int threads;
  const char* output_paf;
} gn_mm2_request;

// Returns 0 on success, non-zero on failure.
int gn_mm2_align_to_paf(const gn_mm2_request* req, char* err_buf, int err_buf_len);

#ifdef __cplusplus
}
#endif
