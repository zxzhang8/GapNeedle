#include "minimap2_bridge.h"

#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <io.h>
#define GN_DUP _dup
#define GN_DUP2 _dup2
#define GN_FILENO _fileno
#define GN_CLOSE _close
#else
#include <unistd.h>
#define GN_DUP dup
#define GN_DUP2 dup2
#define GN_FILENO fileno
#define GN_CLOSE close
#endif

#include "minimap.h"

static void gn_set_err(char* err_buf, int err_buf_len, const char* msg) {
  if (err_buf && err_buf_len > 0) {
    snprintf(err_buf, (size_t)err_buf_len, "%s", msg ? msg : "unknown minimap2 error");
  }
}

typedef struct {
  char* s;
  size_t n;
  size_t m;
} gn_str_t;

static void gn_str_free(gn_str_t* s) {
  if (!s) return;
  free(s->s);
  s->s = 0;
  s->n = 0;
  s->m = 0;
}

static int gn_str_reserve(gn_str_t* s, size_t cap) {
  if (cap <= s->m) return 0;
  size_t new_cap = s->m > 0 ? s->m : 256;
  while (new_cap < cap) new_cap <<= 1;
  char* p = (char*)realloc(s->s, new_cap);
  if (!p) return -1;
  s->s = p;
  s->m = new_cap;
  return 0;
}

static int gn_str_set(gn_str_t* s, const char* src, size_t len) {
  if (gn_str_reserve(s, len + 1) != 0) return -1;
  if (len > 0) memcpy(s->s, src, len);
  s->s[len] = '\0';
  s->n = len;
  return 0;
}

static int gn_str_append_char(gn_str_t* s, char ch) {
  if (gn_str_reserve(s, s->n + 2) != 0) return -1;
  s->s[s->n++] = ch;
  s->s[s->n] = '\0';
  return 0;
}

static size_t gn_header_name(const char* line, char* out, size_t out_cap) {
  size_t i = 0;
  if (out_cap == 0) return 0;
  while (line[i] == ' ' || line[i] == '\t') ++i;
  size_t j = 0;
  while (line[i] != '\0' && line[i] != '\r' && line[i] != '\n' &&
         line[i] != ' ' && line[i] != '\t') {
    if (j + 1 < out_cap) out[j] = line[i];
    ++j;
    ++i;
  }
  if (j < out_cap) out[j] = '\0';
  else out[out_cap - 1] = '\0';
  return j;
}

static int gn_load_fasta_seq(const char* path, const char* wanted_name, gn_str_t* out_name, gn_str_t* out_seq) {
  FILE* fp = fopen(path, "rb");
  if (!fp) return -1;

  char line[8192];
  int collecting = 0;
  int seen_header = 0;
  int found = 0;
  while (fgets(line, (int)sizeof(line), fp)) {
    if (line[0] == '>') {
      char name[4096];
      (void)gn_header_name(line + 1, name, sizeof(name));
      if (collecting) {
        found = 1;
        break;
      }
      seen_header = 1;
      if (!wanted_name || wanted_name[0] == '\0') {
        collecting = 1;
      } else if (strcmp(name, wanted_name) == 0) {
        collecting = 1;
      } else {
        collecting = 0;
      }
      if (collecting && gn_str_set(out_name, name, strlen(name)) != 0) {
        fclose(fp);
        return -2;
      }
      continue;
    }
    if (!collecting) continue;
    for (size_t i = 0; line[i] != '\0'; ++i) {
      unsigned char c = (unsigned char)line[i];
      if (c == '\r' || c == '\n' || c == ' ' || c == '\t') continue;
      if (gn_str_append_char(out_seq, (char)toupper(c)) != 0) {
        fclose(fp);
        return -2;
      }
    }
  }
  if (!found && collecting && out_name->n > 0) found = 1;
  fclose(fp);
  if (!seen_header || !found || out_seq->n == 0) return -3;
  return 0;
}

static void gn_write_cigar(FILE* out, const mm_reg1_t* r) {
  if (!r->p || r->p->n_cigar == 0) return;
  fputs("\tcg:Z:", out);
  for (uint32_t i = 0; i < r->p->n_cigar; ++i) {
    uint32_t c = r->p->cigar[i];
    uint32_t len = c >> 4;
    uint32_t op = c & 0xfU;
    char op_ch = MM_CIGAR_STR[op];
    fprintf(out, "%u%c", len, op_ch);
  }
}

static void gn_free_regs(mm_reg1_t* regs, int n_regs) {
  if (!regs) return;
  for (int i = 0; i < n_regs; ++i) free(regs[i].p);
  free(regs);
}

static int gn_write_single_fasta(const char* path, const char* name, const char* seq) {
  FILE* fp = fopen(path, "wb");
  if (!fp) return -1;
  fprintf(fp, ">%s\n", name ? name : "query");
  size_t n = strlen(seq);
  for (size_t i = 0; i < n; i += 80) {
    size_t len = n - i;
    if (len > 80) len = 80;
    fwrite(seq + i, 1, len, fp);
    fputc('\n', fp);
  }
  fclose(fp);
  return 0;
}

static int gn_filter_paf_target(const char* in_path, const char* out_path, const char* target, int* n_written) {
  FILE* in = fopen(in_path, "rb");
  FILE* out = fopen(out_path, "wb");
  if (!in || !out) {
    if (in) fclose(in);
    if (out) fclose(out);
    return -1;
  }
  char line[1 << 16];
  int written = 0;
  while (fgets(line, (int)sizeof(line), in)) {
    char* cols[7] = {0};
    int c = 0;
    cols[c++] = line;
    for (char* p = line; *p && c < 7; ++p) {
      if (*p == '\t') cols[c++] = p + 1;
    }
    if (c < 6) continue;
    char* tname = cols[5];
    char* end = tname;
    while (*end && *end != '\t' && *end != '\n' && *end != '\r') ++end;
    *end = '\0';
    if (!target || target[0] == '\0' || strcmp(tname, target) == 0) {
      *end = '\t';
      fputs(line, out);
      ++written;
    }
  }
  fclose(in);
  fclose(out);
  *n_written = written;
  return 0;
}

int gn_mm2_align_to_paf(const gn_mm2_request* req, char* err_buf, int err_buf_len) {
  if (!req || !req->target_fasta || !req->query_fasta || !req->output_paf) {
    gn_set_err(err_buf, err_buf_len, "invalid minimap2 request");
    return -1;
  }

  mm_idxopt_t idx_opt;
  mm_mapopt_t map_opt;
  const char* preset = (req->preset && req->preset[0]) ? req->preset : "asm10";
  const int n_threads = req->threads > 0 ? req->threads : 1;
  mm_idx_reader_t* rdr = 0;
  mm_idx_t* mi = 0;
  int n_written = 0;
  int saved_stdout_fd = -1;
  gn_str_t q_name = {0, 0, 0};
  gn_str_t q_seq = {0, 0, 0};
  gn_str_t tmp_query = {0, 0, 0};
  gn_str_t tmp_raw_paf = {0, 0, 0};
  int rc = -2;

  mm_verbose = 2;
  mm_set_opt(0, &idx_opt, &map_opt);
  if (mm_set_opt(preset, &idx_opt, &map_opt) < 0) {
    gn_set_err(err_buf, err_buf_len, "unknown minimap2 preset");
    return -3;
  }
  map_opt.flag |= MM_F_CIGAR | MM_F_OUT_CG;
  map_opt.flag |= MM_F_NO_PRINT_2ND;

  rc = gn_load_fasta_seq(req->query_fasta, req->query_seq, &q_name, &q_seq);
  if (rc != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to load selected query sequence from FASTA");
    rc = -8;
    goto cleanup;
  }

  if (gn_str_set(&tmp_query, req->output_paf, strlen(req->output_paf)) != 0 ||
      gn_str_append_char(&tmp_query, '.') != 0 ||
      gn_str_set(&tmp_raw_paf, req->output_paf, strlen(req->output_paf)) != 0 ||
      gn_str_append_char(&tmp_raw_paf, '.') != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory for temporary paths");
    rc = -11;
    goto cleanup;
  }
  if (gn_str_set(&tmp_query, req->output_paf, strlen(req->output_paf)) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory for temporary query path");
    rc = -11;
    goto cleanup;
  }
  if (gn_str_set(&tmp_raw_paf, req->output_paf, strlen(req->output_paf)) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory for temporary paf path");
    rc = -11;
    goto cleanup;
  }
  if (gn_str_reserve(&tmp_query, tmp_query.n + 16) != 0 ||
      gn_str_reserve(&tmp_raw_paf, tmp_raw_paf.n + 16) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory for temporary suffix");
    rc = -11;
    goto cleanup;
  }
  strcat(tmp_query.s, ".query.tmp.fa");
  tmp_query.n = strlen(tmp_query.s);
  strcat(tmp_raw_paf.s, ".raw.tmp.paf");
  tmp_raw_paf.n = strlen(tmp_raw_paf.s);

  if (gn_write_single_fasta(tmp_query.s, q_name.s, q_seq.s) != 0) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len, "failed to write temp query FASTA: %s", strerror(errno));
    }
    rc = -12;
    goto cleanup;
  }

  FILE* raw_out = fopen(tmp_raw_paf.s, "wb");
  if (!raw_out) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len, "failed to open temp output PAF: %s", strerror(errno));
    }
    rc = -4;
    goto cleanup;
  }
  fflush(stdout);
  saved_stdout_fd = GN_DUP(GN_FILENO(stdout));
  if (saved_stdout_fd < 0 || GN_DUP2(GN_FILENO(raw_out), GN_FILENO(stdout)) < 0) {
    fclose(raw_out);
    gn_set_err(err_buf, err_buf_len, "failed to redirect minimap2 output");
    rc = -13;
    goto cleanup;
  }
  fclose(raw_out);

  rdr = mm_idx_reader_open(req->target_fasta, &idx_opt, 0);
  if (!rdr) {
    gn_set_err(err_buf, err_buf_len, "failed to open target fasta/index");
    rc = -5;
    goto cleanup;
  }

  int target_found_any_part = 0;
  while ((mi = mm_idx_reader_read(rdr, n_threads)) != 0) {
    int target_id = -1;
    mm_mapopt_update(&map_opt, mi);
    if (req->target_seq && req->target_seq[0] != '\0') {
      mm_idx_index_name(mi);
      target_id = mm_idx_name2id(mi, req->target_seq);
      if (target_id < 0) {
        mm_idx_destroy(mi);
        mi = 0;
        continue;  // target sequence is not in this index part; skip costly mapping
      }
      target_found_any_part = 1;
    }

    if (mm_map_file(mi, tmp_query.s, &map_opt, n_threads) < 0) {
      gn_set_err(err_buf, err_buf_len, "minimap2 mapping failed");
      rc = -6;
      mm_idx_destroy(mi);
      mi = 0;
      goto cleanup;
    }
    mm_idx_destroy(mi);
    mi = 0;
  }

  fflush(stdout);
  GN_DUP2(saved_stdout_fd, GN_FILENO(stdout));
  GN_CLOSE(saved_stdout_fd);
  saved_stdout_fd = -1;

  if (gn_filter_paf_target(tmp_raw_paf.s, req->output_paf, req->target_seq, &n_written) != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to finalize filtered PAF");
    rc = -14;
    goto cleanup;
  }

  if (req->target_seq && req->target_seq[0] != '\0' && !target_found_any_part) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len,
               "target sequence '%s' not found in target FASTA/index", req->target_seq);
    }
    rc = -10;
  } else if (n_written == 0) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len,
               "no alignments produced for query '%s' against target '%s'",
               req->query_seq ? req->query_seq : "*",
               req->target_seq ? req->target_seq : "*");
    }
    rc = -9;
  } else {
    rc = 0;
  }

cleanup:
  if (mi) {
    mm_idx_destroy(mi);
  }
  if (saved_stdout_fd >= 0) {
    fflush(stdout);
    GN_DUP2(saved_stdout_fd, GN_FILENO(stdout));
    GN_CLOSE(saved_stdout_fd);
  }
  if (rdr) {
    mm_idx_reader_close(rdr);
  }
  gn_str_free(&q_name);
  gn_str_free(&q_seq);
  if (tmp_query.s) remove(tmp_query.s);
  if (tmp_raw_paf.s) remove(tmp_raw_paf.s);
  gn_str_free(&tmp_query);
  gn_str_free(&tmp_raw_paf);
  if (rc != 0 && (!err_buf || err_buf_len <= 0 || err_buf[0] == '\0')) {
    gn_set_err(err_buf, err_buf_len, "minimap2 bridge failed");
  }
  return rc;
}
