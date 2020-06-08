int displayVars(const char * path);
int varShape(int ncid, int varid);

#define CHECK_NC_NOERR(x) do {\
  int retval = (x); \
  if (retval != NC_NOERR) { \
    fprintf(stderr, "Runtime error: %s returned %d at %s:%d", #x, retval, __FILE__, __LINE__); \
    return retval;\
  }\
} while (0)