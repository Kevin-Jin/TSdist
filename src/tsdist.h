#define ELT(array, rowLength, i, j) ((array)[(i) * (rowLength) + (j)])
#define PATH_X_GAP 1
#define PATH_Y_GAP -1
#define PATH_NO_GAP 0

#ifndef MAX
#define MAX(a,b) ({ \
  __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; \
})
#endif

#ifndef MIN
#define MIN(a,b) ({ \
  __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a < _b ? _a : _b; \
})
#endif

#ifndef ABS
#define ABS(x) ({ \
  __typeof__ (x) _x = (x); \
  _x > 0 ? _x : -_x; \
})
#endif






