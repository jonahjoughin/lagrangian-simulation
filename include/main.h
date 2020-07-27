#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

#define CHECK_NC_NOERR(x)                                                                          \
  do {                                                                                             \
    int retval = (x);                                                                              \
    if (retval != NC_NOERR) {                                                                      \
      fprintf(stderr, "Runtime error: %s returned %d at %s:%d\n", #x, retval, __FILE__, __LINE__); \
      fprintf(stderr, "%s\n", nc_strerror(retval));                                                \
      return retval;                                                                               \
    }                                                                                              \
  } while (0)

#endif