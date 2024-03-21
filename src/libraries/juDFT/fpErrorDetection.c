//#ifdef CPP_DEBUG // This is commented out because at the moment this symbol is not defined for the compilation of c files.

#include <fenv.h>

int startFPErrorDetection()
{
   return feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

int stopFPErrorDetection()
{
   return fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

//#endif
