#ifndef __PYX_HAVE__bem__pytriangle
#define __PYX_HAVE__bem__pytriangle


#ifndef __PYX_HAVE_API__bem__pytriangle

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(int) triunsuitable(double *, double *, double *, double);

#endif /* !__PYX_HAVE_API__bem__pytriangle */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initpytriangle(void);
#else
PyMODINIT_FUNC PyInit_pytriangle(void);
#endif

#endif /* !__PYX_HAVE__bem__pytriangle */
