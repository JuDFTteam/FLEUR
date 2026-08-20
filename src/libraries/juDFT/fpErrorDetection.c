#if defined(__APPLE__) && (defined(__x86_64__) || defined(__i386__))
/* macOS does not provide feenableexcept/fedisableexcept (Linux/glibc extensions). */
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
#elif defined(__GLIBC__)
/* Linux and other platforms with glibc */
#include <fenv.h>
int startFPErrorDetection()
{
   return feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

int stopFPErrorDetection()
{
   return fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#else
/* ARM64 macOS: fenv exception trapping not reliably supported */
int startFPErrorDetection() { return -1; }
int stopFPErrorDetection()  { return -1; }
#endif
