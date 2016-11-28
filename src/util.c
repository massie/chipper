/** Generic utility functions */
#include "util.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

gzFile gzFile_create(const char *filename, char *mode) {
  FILE *stream = NULL;
  if (filename == NULL) {
    if (strcmp(mode, "r") == 0) {
      stream = stdin;
    } else if (strcmp(mode, "w") == 0) {
      stream = stdout;
    } else {
      fprintf(stderr, "Unknown file mode='%s'\n", mode);
      exit(EXIT_FAILURE);
    }
  } else {
    stream = fopen(filename, mode);
    if (!stream) {
      fprintf(stderr, "file:%s, error: %s", filename, strerror(errno));
      exit(EXIT_FAILURE);
    }
  }
  return gzdopen(fileno(stream), mode);
}
