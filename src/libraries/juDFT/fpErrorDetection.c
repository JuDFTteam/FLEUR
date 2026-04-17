
#include <fenv.h>

#if defined(__APPLE__)
/* macOS does not provide feenableexcept/fedisableexcept (Linux/glibc extensions). */
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>

int startFPErrorDetection()
{
   _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
      ~(_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW));
   return 0;
}

int stopFPErrorDetection()
{
   _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() |
      (_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW));
   return 0;
}

#else
/* ARM64 macOS: fenv exception trapping not reliably supported */
int startFPErrorDetection() { return -1; }
int stopFPErrorDetection()  { return -1; }
#endif

#else
/* Linux and other platforms with glibc */
int startFPErrorDetection()
{
   return feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

int stopFPErrorDetection()
{
   return fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#endif

