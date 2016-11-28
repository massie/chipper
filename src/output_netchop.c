#include "output.h"
#include "project.h"
#include "util.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct chipper_netchop_output {
  struct chipper_output output;
  FILE *file;
};

static char cleavage_code(char *seq, double *predictions, int offset,
                          double cutoff, int *cleavage_sites) {
  double prediction = predictions[offset];
  char code = '?';
  int is_cleavage_site;

  if (!isnan(prediction)) {
    is_cleavage_site = prediction >= cutoff;
    if (is_cleavage_site) {
      code = 'S';
      (*cleavage_sites)++;
    } else {
      code = '.';
    }
  }
  return code;
}

static void write_netchop_footer(FILE *file, int cleavage_sites, int seq_len,
                                 const char *name) {
  fprintf(file, "--------------------------------------\n\n");
  fprintf(file, "Number of cleavage sites %d. Number of amino acids %d. "
                "Protein name %s\n\n",
          cleavage_sites, seq_len, name);
  fprintf(file, "--------------------------------------\n");
}

static void netchop_short_write_record(struct chipper_output *co, char *name,
                                       int name_len, char *comment,
                                       int comment_len, char *seq, int seq_len,
                                       double *predictions) {
  struct chipper_netchop_output *nc =
      container_of(co, struct chipper_netchop_output, output);
  const int output_width = 80;
  int i, j, cleavage_sites = 0, remaining_positions;

  fprintf(nc->file, "%d %s\n", seq_len, name);
  for (i = 1; i < seq_len; i++) {
    fprintf(nc->file, "%c", seq[i - 1]);
    if (i % output_width == 0) {
      fprintf(nc->file, "\n");
      for (j = i - output_width; j < i; j++) {
        fprintf(nc->file, "%c", cleavage_code(seq, predictions, j, co->cutoff,
                                              &cleavage_sites));
      }
      fprintf(nc->file, "\n");
    }
  }
  remaining_positions = seq_len % output_width;
  if (remaining_positions) {
    fprintf(nc->file, "\n");
    for (i = seq_len - remaining_positions; i < seq_len - 1; i++) {
      fprintf(nc->file, "%c",
              cleavage_code(seq, predictions, i, co->cutoff, &cleavage_sites));
    }
    fprintf(nc->file, "\n");
  }
  write_netchop_footer(nc->file, cleavage_sites, seq_len, name);
}

static void netchop_write_record(struct chipper_output *co, char *name,
                                 int name_len, char *comment, int comment_len,
                                 char *seq, int seq_len, double *predictions) {
  struct chipper_netchop_output *nc =
      container_of(co, struct chipper_netchop_output, output);
  int i, cleavage_sites = 0;

  fprintf(nc->file, "--------------------------------------\n");
  fprintf(nc->file, " pos  AA  C      score      Ident\n");
  fprintf(nc->file, "--------------------------------------\n");
  for (i = 0; i < seq_len; i++) {
    fprintf(nc->file, "%4d   %c  %c    %f %s\n", i + 1, seq[i],
            cleavage_code(seq, predictions, i, co->cutoff, &cleavage_sites),
            predictions[i], name);
  }
  write_netchop_footer(nc->file, cleavage_sites, seq_len, name);
}

static void netchop_close(struct chipper_output *co) {
  struct chipper_netchop_output *nc =
      container_of(co, struct chipper_netchop_output, output);
  fclose(nc->file);
}

struct chipper_output *new_netchop_output(const char *output_path,
                                          double cutoff, int short_format) {
  struct chipper_netchop_output *nc =
      malloc(sizeof(struct chipper_netchop_output));
  if (nc == NULL) {
    fprintf(stderr, "Can't allocate memory for netchop output\n");
    exit(EXIT_FAILURE);
  }
  if (output_path) {
    nc->file = fopen(output_path, "w");
  } else {
    nc->file = stdout;
  }
  if (nc->file == NULL) {
    fprintf(stderr, "Unable to open file '%s' for writing. Exiting.\n",
            output_path);
    exit(EXIT_FAILURE);
  }
  /* Write the netchop header */
  fprintf(nc->file,
          "chipper %s predictions using version C-term. Threshold %f\n\n",
          CHIPPER_VERSION, cutoff);

  struct chipper_output *out = &nc->output;

  if (short_format) {
    out->format = NETCHOP_SHORT_OUTPUT;
    out->write_record_cb = netchop_short_write_record;
  } else {
    out->format = NETCHOP_OUTPUT;
    out->write_record_cb = netchop_write_record;
  }
  out->close_cb = netchop_close;
  out->cutoff = cutoff;
  return out;
}
