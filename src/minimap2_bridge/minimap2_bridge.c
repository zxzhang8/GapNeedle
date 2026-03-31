#include "minimap2_bridge.h"

#include <errno.h>
#include <ctype.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#ifdef _WIN32
#include <direct.h>
#include <io.h>
#include <process.h>
#include <windows.h>
#define GN_OPEN _open
#define GN_CLOSE _close
#define GN_WRITE _write
#define GN_MKDIR(path) _mkdir(path)
#define GN_STAT _stat64
#define GN_STAT_STRUCT struct _stat64
#define GN_PID _getpid
#else
#include <unistd.h>
#define GN_OPEN open
#define GN_CLOSE close
#define GN_WRITE write
#define GN_MKDIR(path) mkdir(path, 0777)
#define GN_STAT stat
#define GN_STAT_STRUCT struct stat
#define GN_PID getpid
#endif

#include "minimap.h"
#include "bseq.h"
#include "kseq.h"
#include "mmpriv.h"

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
  if (len > 0 && src) memcpy(s->s, src, len);
  s->s[len] = '\0';
  s->n = len;
  return 0;
}

static int gn_str_set_cstr(gn_str_t* s, const char* src) {
  return gn_str_set(s, src ? src : "", src ? strlen(src) : 0);
}

static int gn_str_append(gn_str_t* s, const char* src) {
  const size_t len = src ? strlen(src) : 0;
  if (gn_str_reserve(s, s->n + len + 1) != 0) return -1;
  if (len > 0) memcpy(s->s + s->n, src, len);
  s->n += len;
  s->s[s->n] = '\0';
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

static char gn_comp_base(char c) {
  switch ((unsigned char)toupper((unsigned char)c)) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'U': return 'A';
    case 'R': return 'Y';
    case 'Y': return 'R';
    case 'S': return 'S';
    case 'W': return 'W';
    case 'K': return 'M';
    case 'M': return 'K';
    case 'B': return 'V';
    case 'D': return 'H';
    case 'H': return 'D';
    case 'V': return 'B';
    default: return 'N';
  }
}

static void gn_reverse_complement_inplace(gn_str_t* seq) {
  if (!seq || !seq->s || seq->n == 0) return;
  size_t i = 0;
  size_t j = seq->n - 1;
  while (i < j) {
    const char li = seq->s[i];
    const char rj = seq->s[j];
    seq->s[i] = gn_comp_base(rj);
    seq->s[j] = gn_comp_base(li);
    ++i;
    --j;
  }
  if (i == j) seq->s[i] = gn_comp_base(seq->s[i]);
}

static void gn_free_regs(mm_reg1_t* regs, int n_regs) {
  if (!regs) return;
  for (int i = 0; i < n_regs; ++i) free(regs[i].p);
  free(regs);
}

static int gn_sleep_ms(int ms) {
#ifdef _WIN32
  Sleep((DWORD)ms);
  return 0;
#else
  struct timespec ts;
  ts.tv_sec = ms / 1000;
  ts.tv_nsec = (long)(ms % 1000) * 1000000L;
  while (nanosleep(&ts, &ts) != 0) {
    if (errno != EINTR) return -1;
  }
  return 0;
#endif
}

static double gn_now_ms(void) {
#ifdef _WIN32
  return (double)GetTickCount64();
#else
  struct timespec ts;
#ifdef CLOCK_MONOTONIC
  clock_gettime(CLOCK_MONOTONIC, &ts);
#else
  clock_gettime(CLOCK_REALTIME, &ts);
#endif
  return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1000000.0;
#endif
}

static void gn_trace_appendf(char* trace_buf, int trace_buf_len, const char* fmt, ...) {
  if (!trace_buf || trace_buf_len <= 1 || !fmt) return;
  size_t used = strlen(trace_buf);
  if (used >= (size_t)trace_buf_len - 1) return;
  va_list args;
  va_start(args, fmt);
  const int n = vsnprintf(trace_buf + used, (size_t)trace_buf_len - used, fmt, args);
  va_end(args);
  if (n <= 0) return;
  used = strlen(trace_buf);
  if (used < (size_t)trace_buf_len - 1 && trace_buf[used - 1] != '\n') {
    trace_buf[used++] = '\n';
    trace_buf[used] = '\0';
  }
}

static int gn_file_stat(const char* path, long long* size_out, long long* mtime_out) {
  GN_STAT_STRUCT st;
  if (GN_STAT(path, &st) != 0) return -1;
  if (size_out) *size_out = (long long)st.st_size;
  if (mtime_out) *mtime_out = (long long)st.st_mtime;
  return 0;
}

static int gn_file_exists(const char* path) {
  GN_STAT_STRUCT st;
  return path && GN_STAT(path, &st) == 0;
}

static int gn_remove_if_exists(const char* path) {
  if (!path || !path[0]) return 0;
  if (!gn_file_exists(path)) return 0;
  return remove(path);
}

static uint64_t gn_hash_bytes(uint64_t h, const void* data, size_t n) {
  const unsigned char* p = (const unsigned char*)data;
  for (size_t i = 0; i < n; ++i) {
    h ^= (uint64_t)p[i];
    h *= 1099511628211ULL;
  }
  return h;
}

static uint64_t gn_hash_cstr(uint64_t h, const char* s) {
  if (!s) return gn_hash_bytes(h, "-", 1);
  return gn_hash_bytes(h, s, strlen(s));
}

static int gn_append_u64_hex(gn_str_t* out, uint64_t v) {
  char buf[17];
  snprintf(buf, sizeof(buf), "%016llx", (unsigned long long)v);
  return gn_str_append(out, buf);
}

static int gn_make_dirs(const char* path) {
  if (!path || !path[0]) return 0;
  gn_str_t cur = {0, 0, 0};
  const size_t len = strlen(path);
  size_t i = 0;
#ifdef _WIN32
  if (len >= 2 && path[1] == ':') {
    if (gn_str_set(&cur, path, 2) != 0) {
      gn_str_free(&cur);
      return -1;
    }
    i = 2;
  }
#endif
  if (i < len && (path[i] == '/' || path[i] == '\\')) {
    if (gn_str_append_char(&cur, path[i]) != 0) {
      gn_str_free(&cur);
      return -1;
    }
    ++i;
  }
  for (; i < len; ++i) {
    const char ch = path[i];
    if (gn_str_append_char(&cur, ch) != 0) {
      gn_str_free(&cur);
      return -1;
    }
    if (ch == '/' || ch == '\\') {
      if (cur.n > 0 && GN_MKDIR(cur.s) != 0 && errno != EEXIST) {
        gn_str_free(&cur);
        return -1;
      }
    }
  }
  if (cur.n > 0 && GN_MKDIR(cur.s) != 0 && errno != EEXIST) {
    gn_str_free(&cur);
    return -1;
  }
  gn_str_free(&cur);
  return 0;
}

static int gn_path_join(gn_str_t* out, const char* dir, const char* base) {
  if (gn_str_set_cstr(out, dir) != 0) return -1;
  if (out->n > 0 && out->s[out->n - 1] != '/' && out->s[out->n - 1] != '\\') {
    if (gn_str_append_char(out, '/') != 0) return -1;
  }
  return gn_str_append(out, base);
}

static int gn_atomic_write_text(const char* path, const char* text) {
  gn_str_t tmp = {0, 0, 0};
  int rc = -1;
  FILE* fp = 0;
  if (gn_str_set_cstr(&tmp, path) != 0 || gn_str_append(&tmp, ".tmp") != 0) goto cleanup;
  fp = fopen(tmp.s, "wb");
  if (!fp) goto cleanup;
  if (text && text[0] != '\0') fwrite(text, 1, strlen(text), fp);
  if (fclose(fp) != 0) {
    fp = 0;
    goto cleanup;
  }
  fp = 0;
  if (rename(tmp.s, path) != 0) goto cleanup;
  rc = 0;
cleanup:
  if (fp) fclose(fp);
  if (rc != 0 && tmp.s) remove(tmp.s);
  gn_str_free(&tmp);
  return rc;
}

static int gn_try_acquire_lock(const char* lock_path) {
  int fd = GN_OPEN(lock_path, O_WRONLY | O_CREAT | O_EXCL, 0666);
  if (fd < 0) return -1;
  char buf[64];
  const int len = snprintf(buf, sizeof(buf), "%ld\n", (long)GN_PID());
  if (len > 0) (void)GN_WRITE(fd, buf, (unsigned int)len);
  GN_CLOSE(fd);
  return 0;
}

static void gn_release_lock(const char* lock_path) {
  if (lock_path && lock_path[0]) remove(lock_path);
}

static int gn_build_index_file(const char* fasta_path,
                               const mm_idxopt_t* idx_opt,
                               int n_threads,
                               const char* out_index,
                               char* err_buf,
                               int err_buf_len) {
  mm_idx_reader_t* rdr = mm_idx_reader_open(fasta_path, idx_opt, out_index);
  mm_idx_t* mi = 0;
  int n_parts = 0;
  if (!rdr) {
    gn_set_err(err_buf, err_buf_len, "failed to open source FASTA for index build");
    return -1;
  }
  while ((mi = mm_idx_reader_read(rdr, n_threads)) != 0) {
    ++n_parts;
    mm_idx_destroy(mi);
    mi = 0;
  }
  mm_idx_reader_close(rdr);
  if (n_parts == 0) {
    gn_set_err(err_buf, err_buf_len, "minimap2 index build produced no parts");
    return -1;
  }
  return 0;
}

static int gn_write_meta_file(const char* meta_path,
                              const char* target_fasta,
                              const char* preset,
                              const mm_idxopt_t* idx_opt,
                              int reverse_target,
                              const char* target_seq,
                              const char* index_path) {
  char buf[2048];
  snprintf(buf,
           sizeof(buf),
           "{\n"
           "  \"version\": \"gn_mm2_index_v1\",\n"
           "  \"target_fasta\": \"%s\",\n"
           "  \"preset\": \"%s\",\n"
           "  \"reverse_target\": %s,\n"
           "  \"target_seq\": \"%s\",\n"
           "  \"index\": \"%s\",\n"
           "  \"idxopt\": {\"k\": %d, \"w\": %d, \"flag\": %d, \"bucket_bits\": %d, \"mini_batch_size\": %lld, \"batch_size\": %llu}\n"
           "}\n",
           target_fasta ? target_fasta : "",
           preset ? preset : "",
           reverse_target ? "true" : "false",
           target_seq ? target_seq : "",
           index_path ? index_path : "",
           idx_opt->k,
           idx_opt->w,
           idx_opt->flag,
           idx_opt->bucket_bits,
           (long long)idx_opt->mini_batch_size,
           (unsigned long long)idx_opt->batch_size);
  return gn_atomic_write_text(meta_path, buf);
}

static int gn_prepare_index_cache(const gn_mm2_request* req,
                                  const mm_idxopt_t* idx_opt,
                                  const char* preset,
                                  int n_threads,
                                  const char* source_fasta,
                                  const char* reverse_target_name,
                                  const char* reverse_target_seq,
                                  gn_str_t* out_index_path,
                                  gn_str_t* out_cached_target_fasta,
                                  char* trace_buf,
                                  int trace_buf_len,
                                  char* err_buf,
                                  int err_buf_len) {
  gn_str_t hash = {0, 0, 0};
  gn_str_t index_path = {0, 0, 0};
  gn_str_t meta_path = {0, 0, 0};
  gn_str_t target_fasta_path = {0, 0, 0};
  gn_str_t lock_path = {0, 0, 0};
  gn_str_t tmp_index = {0, 0, 0};
  long long target_size = 0;
  long long target_mtime = 0;
  uint64_t h = 1469598103934665603ULL;
  int have_lock = 0;
  int waited_for_lock = 0;
  int rc = -1;
  const double t0 = gn_now_ms();

  if (gn_make_dirs(req->index_cache_dir) != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to create minimap2 index cache directory");
    goto cleanup;
  }
  if (gn_file_stat(req->target_fasta, &target_size, &target_mtime) != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to stat target FASTA for index cache");
    goto cleanup;
  }

  h = gn_hash_cstr(h, "gn_mm2_index_v1|");
  h = gn_hash_cstr(h, req->target_fasta);
  h = gn_hash_bytes(h, &target_size, sizeof(target_size));
  h = gn_hash_bytes(h, &target_mtime, sizeof(target_mtime));
  h = gn_hash_cstr(h, preset ? preset : "");
  h = gn_hash_bytes(h, idx_opt, sizeof(*idx_opt));
  h = gn_hash_bytes(h, &req->reverse_target, sizeof(req->reverse_target));
  if (req->reverse_target) {
    h = gn_hash_cstr(h, reverse_target_name ? reverse_target_name : req->target_seq);
    h = gn_hash_cstr(h, reverse_target_seq ? reverse_target_seq : "");
  }
  if (gn_append_u64_hex(&hash, h) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory while preparing index cache");
    goto cleanup;
  }

  if (gn_path_join(&index_path, req->index_cache_dir, hash.s) != 0 ||
      gn_str_append(&index_path, ".mmi") != 0 ||
      gn_path_join(&meta_path, req->index_cache_dir, hash.s) != 0 ||
      gn_str_append(&meta_path, ".meta.json") != 0 ||
      gn_path_join(&lock_path, req->index_cache_dir, hash.s) != 0 ||
      gn_str_append(&lock_path, ".lock") != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory while building cache paths");
    goto cleanup;
  }

  if (req->reverse_target) {
    if (gn_path_join(&target_fasta_path, req->index_cache_dir, hash.s) != 0 ||
        gn_str_append(&target_fasta_path, ".target.fa") != 0) {
      gn_set_err(err_buf, err_buf_len, "out of memory while preparing cached target FASTA path");
      goto cleanup;
    }
  }

  for (;;) {
    if (gn_file_exists(index_path.s)) {
      if (gn_str_set_cstr(out_index_path, index_path.s) != 0) {
        gn_set_err(err_buf, err_buf_len, "out of memory while returning cached index path");
        goto cleanup;
      }
      if (req->reverse_target && gn_str_set_cstr(out_cached_target_fasta, target_fasta_path.s) != 0) {
        gn_set_err(err_buf, err_buf_len, "out of memory while returning cached target FASTA path");
        goto cleanup;
      }
      {
        char msg[512];
        snprintf(msg, sizeof(msg), "cache hit (%s)", index_path.s);
        gn_trace_appendf(trace_buf, trace_buf_len, "Index: cache hit: %s", index_path.s);
      }
      rc = 0;
      goto cleanup;
    }
    if (gn_try_acquire_lock(lock_path.s) == 0) {
      have_lock = 1;
      break;
    }
    if (errno != EEXIST) {
      gn_set_err(err_buf, err_buf_len, "failed to create minimap2 index cache lock");
      goto cleanup;
    }
    waited_for_lock = 1;
    if (gn_sleep_ms(200) != 0) {
      gn_set_err(err_buf, err_buf_len, "interrupted while waiting for index cache lock");
      goto cleanup;
    }
  }

  if (req->reverse_target && !gn_file_exists(target_fasta_path.s)) {
    gn_str_t tmp_target = {0, 0, 0};
    if (gn_str_set_cstr(&tmp_target, target_fasta_path.s) != 0 || gn_str_append(&tmp_target, ".tmp") != 0) {
      gn_str_free(&tmp_target);
      gn_set_err(err_buf, err_buf_len, "out of memory while staging cached target FASTA");
      goto cleanup;
    }
    if (gn_write_single_fasta(tmp_target.s, reverse_target_name, reverse_target_seq) != 0) {
      gn_str_free(&tmp_target);
      if (err_buf && err_buf_len > 0) {
        snprintf(err_buf, (size_t)err_buf_len, "failed to write cached reverse target FASTA: %s", strerror(errno));
      }
      goto cleanup;
    }
    if (rename(tmp_target.s, target_fasta_path.s) != 0) {
      remove(tmp_target.s);
      gn_str_free(&tmp_target);
      gn_set_err(err_buf, err_buf_len, "failed to finalize cached reverse target FASTA");
      goto cleanup;
    }
    gn_str_free(&tmp_target);
  }

  if (gn_str_set_cstr(&tmp_index, index_path.s) != 0 || gn_str_append(&tmp_index, ".building") != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory while staging minimap2 index file");
    goto cleanup;
  }
  gn_remove_if_exists(tmp_index.s);
  if (gn_build_index_file(req->reverse_target ? target_fasta_path.s : source_fasta,
                          idx_opt,
                          n_threads,
                          tmp_index.s,
                          err_buf,
                          err_buf_len) != 0) {
    goto cleanup;
  }
  if (rename(tmp_index.s, index_path.s) != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to finalize cached minimap2 index");
    goto cleanup;
  }
  if (gn_write_meta_file(meta_path.s,
                         req->target_fasta,
                         preset,
                         idx_opt,
                         req->reverse_target,
                         req->reverse_target ? reverse_target_name : req->target_seq,
                         index_path.s) != 0) {
    gn_set_err(err_buf, err_buf_len, "failed to write minimap2 index cache metadata");
    goto cleanup;
  }
  if (gn_str_set_cstr(out_index_path, index_path.s) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory while returning built index path");
    goto cleanup;
  }
  if (req->reverse_target && gn_str_set_cstr(out_cached_target_fasta, target_fasta_path.s) != 0) {
    gn_set_err(err_buf, err_buf_len, "out of memory while returning built target FASTA path");
    goto cleanup;
  }
  {
    gn_trace_appendf(trace_buf,
                     trace_buf_len,
                     "Index: cache built: %s (%.1f ms)%s",
                     index_path.s,
                     gn_now_ms() - t0,
                     waited_for_lock ? ", waited for existing builder" : "");
  }
  rc = 0;

cleanup:
  if (rc == 0 && waited_for_lock) {
    gn_trace_appendf(trace_buf, trace_buf_len, "Index: cache lock wait completed");
  }
  if (have_lock) gn_release_lock(lock_path.s);
  if (rc != 0 && tmp_index.s) remove(tmp_index.s);
  gn_str_free(&hash);
  gn_str_free(&index_path);
  gn_str_free(&meta_path);
  gn_str_free(&target_fasta_path);
  gn_str_free(&lock_path);
  gn_str_free(&tmp_index);
  return rc;
}

static int gn_run_single_query_mapping(const mm_idx_t* mi,
                                       const mm_mapopt_t* map_opt,
                                       const char* q_name,
                                       const char* q_seq,
                                       const char* filter_target_name,
                                       FILE* out_fp,
                                       int* n_written,
                                       char* err_buf,
                                       int err_buf_len) {
  mm_tbuf_t* tbuf = 0;
  mm_reg1_t* regs = 0;
  mm_bseq1_t query = {0, 0, 0, 0, 0};
  kstring_t paf = {0, 0, 0};
  int n_regs = 0;
  int rc = -1;

  query.l_seq = (int)strlen(q_seq);
  query.name = (char*)(q_name ? q_name : "query");
  query.seq = (char*)q_seq;

  tbuf = mm_tbuf_init();
  if (!tbuf) {
    gn_set_err(err_buf, err_buf_len, "failed to allocate minimap2 thread buffer");
    goto cleanup;
  }
  regs = mm_map(mi, query.l_seq, query.seq, &n_regs, tbuf, map_opt, query.name);
  if (n_regs > 0 && !regs) {
    gn_set_err(err_buf, err_buf_len, "minimap2 mapping returned null records");
    goto cleanup;
  }
  for (int i = 0; i < n_regs; ++i) {
    const mm_reg1_t* r = &regs[i];
    if ((map_opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent) continue;
    if (filter_target_name && filter_target_name[0] != '\0') {
      const char* t_name = mi->seq[r->rid].name;
      if (!t_name || strcmp(t_name, filter_target_name) != 0) continue;
    }
    mm_write_paf(&paf, mi, &query, r, mm_tbuf_get_km(tbuf), map_opt->flag);
    if (paf.l > 0) {
      fwrite(paf.s, 1, (size_t)paf.l, out_fp);
      fputc('\n', out_fp);
      ++(*n_written);
    }
  }
  rc = 0;

cleanup:
  gn_free_regs(regs, n_regs);
  free(paf.s);
  if (tbuf) mm_tbuf_destroy(tbuf);
  return rc;
}

int gn_mm2_align_to_paf(const gn_mm2_request* req,
                        char* err_buf,
                        int err_buf_len,
                        char* trace_buf,
                        int trace_buf_len) {
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
  FILE* out_fp = 0;
  gn_str_t q_name = {0, 0, 0};
  gn_str_t q_seq = {0, 0, 0};
  gn_str_t t_name = {0, 0, 0};
  gn_str_t t_seq = {0, 0, 0};
  gn_str_t tmp_target = {0, 0, 0};
  gn_str_t cached_index = {0, 0, 0};
  gn_str_t cached_target = {0, 0, 0};
  const char* mapping_source = req->target_fasta;
  const char* filter_target_name = req->target_seq;
  int target_found_any_part = 0;
  int n_written = 0;
  int n_parts_mapped = 0;
  int rc = -2;
  const double total_t0 = gn_now_ms();
  double stage_t0 = total_t0;

  if (trace_buf && trace_buf_len > 0) trace_buf[0] = '\0';
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
  gn_trace_appendf(trace_buf,
                   trace_buf_len,
                   "Timing: load query FASTA '%s' seq '%s': %.1f ms",
                   req->query_fasta,
                   req->query_seq ? req->query_seq : "",
                   gn_now_ms() - stage_t0);
  if (req->reverse_query) {
    const double rc_t0 = gn_now_ms();
    gn_reverse_complement_inplace(&q_seq);
    gn_trace_appendf(trace_buf, trace_buf_len, "Timing: reverse-complement query: %.1f ms", gn_now_ms() - rc_t0);
  }

  if (req->reverse_target) {
    stage_t0 = gn_now_ms();
    if (!req->target_seq || req->target_seq[0] == '\0') {
      gn_set_err(err_buf, err_buf_len, "reverse_target requires target_seq");
      rc = -15;
      goto cleanup;
    }
    rc = gn_load_fasta_seq(req->target_fasta, req->target_seq, &t_name, &t_seq);
    if (rc != 0) {
      gn_set_err(err_buf, err_buf_len, "failed to load selected target sequence from FASTA");
      rc = -16;
      goto cleanup;
    }
    gn_reverse_complement_inplace(&t_seq);
    filter_target_name = t_name.s;
    gn_trace_appendf(trace_buf,
                     trace_buf_len,
                     "Timing: load+reverse target FASTA '%s' seq '%s': %.1f ms",
                     req->target_fasta,
                     req->target_seq ? req->target_seq : "",
                     gn_now_ms() - stage_t0);
  }

  stage_t0 = gn_now_ms();
  if (req->allow_index_cache && req->index_cache_dir && req->index_cache_dir[0] != '\0') {
    if (gn_prepare_index_cache(req,
                               &idx_opt,
                               preset,
                               n_threads,
                               req->target_fasta,
                               t_name.s,
                               t_seq.s,
                               &cached_index,
                               &cached_target,
                               trace_buf,
                               trace_buf_len,
                               err_buf,
                               err_buf_len) != 0) {
      rc = -18;
      goto cleanup;
    }
    mapping_source = cached_index.s;
    gn_trace_appendf(trace_buf, trace_buf_len, "Timing: prepare index cache: %.1f ms", gn_now_ms() - stage_t0);
  } else if (req->reverse_target) {
    if (gn_str_set_cstr(&tmp_target, req->output_paf) != 0 || gn_str_append(&tmp_target, ".target.tmp.fa") != 0) {
      gn_set_err(err_buf, err_buf_len, "out of memory for temporary target path");
      rc = -11;
      goto cleanup;
    }
    if (gn_write_single_fasta(tmp_target.s, t_name.s, t_seq.s) != 0) {
      if (err_buf && err_buf_len > 0) {
        snprintf(err_buf, (size_t)err_buf_len, "failed to write temp target FASTA: %s", strerror(errno));
      }
      rc = -17;
      goto cleanup;
    }
    mapping_source = tmp_target.s;
    gn_trace_appendf(trace_buf, trace_buf_len, "Timing: write temporary reverse target FASTA: %.1f ms", gn_now_ms() - stage_t0);
  }

  stage_t0 = gn_now_ms();
  out_fp = fopen(req->output_paf, "wb");
  if (!out_fp) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len, "failed to open output PAF: %s", strerror(errno));
    }
    rc = -4;
    goto cleanup;
  }
  gn_trace_appendf(trace_buf, trace_buf_len, "Timing: open output PAF '%s': %.1f ms", req->output_paf, gn_now_ms() - stage_t0);

  stage_t0 = gn_now_ms();
  rdr = mm_idx_reader_open(mapping_source, &idx_opt, 0);
  if (!rdr) {
    gn_set_err(err_buf, err_buf_len, "failed to open target fasta/index");
    rc = -5;
    goto cleanup;
  }
  gn_trace_appendf(trace_buf, trace_buf_len, "Timing: open target FASTA/index '%s': %.1f ms", mapping_source, gn_now_ms() - stage_t0);

  while ((mi = mm_idx_reader_read(rdr, n_threads)) != 0) {
    const double part_t0 = gn_now_ms();
    mm_mapopt_update(&map_opt, mi);
    if (filter_target_name && filter_target_name[0] != '\0') {
      mm_idx_index_name(mi);
      if (mm_idx_name2id(mi, filter_target_name) < 0) {
        mm_idx_destroy(mi);
        mi = 0;
        continue;
      }
      target_found_any_part = 1;
    }
    if (gn_run_single_query_mapping(mi,
                                    &map_opt,
                                    q_name.s,
                                    q_seq.s,
                                    filter_target_name,
                                    out_fp,
                                    &n_written,
                                    err_buf,
                                    err_buf_len) != 0) {
      rc = -6;
      mm_idx_destroy(mi);
      mi = 0;
      goto cleanup;
    }
    ++n_parts_mapped;
    gn_trace_appendf(trace_buf,
                     trace_buf_len,
                     "Timing: map index part %d (%u target seqs): %.1f ms",
                     n_parts_mapped,
                     mi->n_seq,
                     gn_now_ms() - part_t0);
    mm_idx_destroy(mi);
    mi = 0;
  }

  stage_t0 = gn_now_ms();
  fflush(out_fp);
  fclose(out_fp);
  out_fp = 0;
  gn_trace_appendf(trace_buf, trace_buf_len, "Timing: flush/finalize output PAF: %.1f ms", gn_now_ms() - stage_t0);

  if (filter_target_name && filter_target_name[0] != '\0' && !target_found_any_part) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len,
               "target sequence '%s' not found in target FASTA/index", filter_target_name);
    }
    rc = -10;
  } else if (n_written == 0) {
    if (err_buf && err_buf_len > 0) {
      snprintf(err_buf, (size_t)err_buf_len,
               "no alignments produced for query '%s' against target '%s'",
               req->query_seq ? req->query_seq : "*",
               filter_target_name ? filter_target_name : "*");
    }
    rc = -9;
  } else {
    rc = 0;
  }
  if (rc == 0) {
    gn_trace_appendf(trace_buf,
                     trace_buf_len,
                     "Timing: total minimap2 pipeline: %.1f ms, parts=%d, records=%d",
                     gn_now_ms() - total_t0,
                     n_parts_mapped,
                     n_written);
  }

cleanup:
  if (mi) mm_idx_destroy(mi);
  if (rdr) mm_idx_reader_close(rdr);
  if (out_fp) fclose(out_fp);
  if (rc != 0) gn_remove_if_exists(req->output_paf);
  gn_str_free(&q_name);
  gn_str_free(&q_seq);
  gn_str_free(&t_name);
  gn_str_free(&t_seq);
  if (tmp_target.s) remove(tmp_target.s);
  gn_str_free(&tmp_target);
  gn_str_free(&cached_index);
  gn_str_free(&cached_target);
  if (rc != 0 && (!err_buf || err_buf_len <= 0 || err_buf[0] == '\0')) {
    gn_set_err(err_buf, err_buf_len, "minimap2 bridge failed");
  }
  return rc;
}
