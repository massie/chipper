#ifndef CHIPPER_UTIL_H_DEFINED
#define CHIPPER_UTIL_H_DEFINED

#define container_of(ptr_, type_, member_)                                     \
  ((type_ *)((char *)ptr_ - (size_t) & ((type_ *)0)->member_))

#define CHECK(x)                                                               \
  do {                                                                         \
    int retval = (x);                                                          \
    if (retval <= 0) {                                                         \
      fprintf(stderr, "Runtime error: %s returned %d at %s:%d", #x, retval,    \
              __FILE__, __LINE__);                                             \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

#define ARRAY_LEN(x) (sizeof(x) / sizeof(x[0]))

#include <zlib.h>
gzFile gzFile_create(const char *filename, char *mode);

#define ARGUMENT_IS_SET(arg) (arg->count == 1)

#endif
