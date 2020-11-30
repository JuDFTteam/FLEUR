/* Code from Agner Fog to improve performance of INTEL compiler/libraries on AMD */



/***********************  intel_cpu_feature_patch.c  **************************
* Author:           Agner Fog
* Date created:     2014-07-30
* Last modified:    2019-12-29
* Source URL:       https://www.agner.org/optimize/intel_dispatch_patch.zip
* Language:         C or C++
*
* Description:
* Patch for Intel compiler version 13.0 and later, including the general
* libraries, LIBM and SVML, but not MKL and VML.
*
* Example of how to patch Intel's CPU feature dispatcher in order to improve
* compatibility of generated code with non-Intel processors.
* In Windows: Use the static link libraries (*.lib), not the dynamic link
* librarise (*.DLL).
* In Linux and Mac: use static linking (*.a) or dynamic linking (*.so).
*
* Include this code in your C or C++ program and call intel_cpu_patch();
* before any call to the library functions.
*
* Copyright (c) 2014-2019. BSD License 2.0
******************************************************************************/
#include <stdint.h>

#ifdef __cplusplus  // use C-style linking
extern "C" {
#endif

// link to Intel libraries
extern int64_t __intel_cpu_feature_indicator;    // CPU feature bits
extern int64_t __intel_cpu_feature_indicator_x;  // CPU feature bits
void __intel_cpu_features_init();                // unfair dispatcher: checks CPU features for Intel CPU's only
void __intel_cpu_features_init_x();              // fair dispatcher: checks CPU features without discriminating by CPU brand

#ifdef __cplusplus
}  // end of extern "C"
#endif

void intel_cpu_patch() {
    // force a re-evaluation of the CPU features without discriminating by CPU brand
    __intel_cpu_feature_indicator = 0;
    __intel_cpu_feature_indicator_x = 0;
    __intel_cpu_features_init_x();
    __intel_cpu_feature_indicator = __intel_cpu_feature_indicator_x;
}

/***********************  intel_mkl_cpuid_patch.c  **************************
* Author:           Agner Fog
* Date created:     2019-12-29
* Source URL:       https://www.agner.org/optimize/intel_dispatch_patch.zip
* Language:         C or C++
*
* Description:
* Patch for Intel Math Kernel Library (MKL) version 14.0 and later, except
* the Vector Math Library (VML).
*
* Example of how to override Intel's CPU feature dispatcher in order to improve
* compatibility of Intel function libraries with non-Intel processors.
*
* Include this code in your C or C++ program and make sure it is linked before
* any Intel libraries. You may need to include intel_mkl_feature_patch.c as well.
*
* Copyright (c) 2019. BSD License 2.0
******************************************************************************/
//#include <stdint.h>
#ifdef __cplusplus  // use C-style linking
extern "C" {
#endif
/*
    // detect if Intel CPU
    int mkl_serv_intel_cpu() {
        return 1;
    }

    // detect if Intel CPU
    int mkl_serv_intel_cpu_true() {
        return 1;
    }


    int mkl_serv_cpuhaspnr_true() {
        return 0;
    }

    int mkl_serv_cpuhaspnr() {
        return 0;
    }

    int mkl_serv_cpuhasnhm() {
        return 0;
    }

    int mkl_serv_cpuisbulldozer() {
        return 0;
    }

    int mkl_serv_cpuiszen() {
        return 1;
    }

    int mkl_serv_cpuisatomsse4_2() {
        return 0;
    }

    int mkl_serv_cpuisatomssse3() {
        return 0;
    }

    int mkl_serv_cpuisitbarcelona() {
        return 0;
    }

    int mkl_serv_cpuisskl() {
        return 0;
    }

    int mkl_serv_cpuisknm() {
        return 0;
    }

    int mkl_serv_cpuisclx() {
        return 0;
    }

    int mkl_serv_get_microarchitecture() {
        // I don't know what this number means
        return 33;
    }
   int mkl_serv_cpuisclx() {
        return 1;
    }
    int mkl_serv_cpuiszen() {
        return 1;
    }
*/
    int mkl_serv_intel_cpu() {
        return 1;
    }

    int mkl_serv_intel_cpu_true() {
        return 1;
    }

#ifdef __cplusplus
}  // end of extern "C"
#endif
/***********************  intel_mkl_feature_patch.c  **************************
* Author:           Agner Fog
* Date created:     2014-07-30
* Last modified:    2019-12-29
* Source URL:       https://www.agner.org/optimize/intel_dispatch_patch.zip
* Language:         C or C++
*
* Description:
* Patch for Intel Math Kernel Library (MKL) version 14.0 and later, except
* the Vector Math Library (VML).
*
* Example of how to patch Intel's CPU feature dispatcher in order to improve
* compatibility of Intel function libraries with non-Intel processors.
* In Windows: Use the static link libraries (*.lib), not the dynamic link
* librarise (*.DLL).
* In Linux and Mac: use static linking (*.a) or dynamic linking (*.so).
*
* Include this code in your C or C++ program and call intel_mkl_patch();
* before any call to the MKL functions. You may need to include
* intel_mkl_cpuid_patch.c as well.
*
* Copyright (c) 2014-2019. BSD License 2.0
******************************************************************************/
//#include <stdint.h>

#ifdef __cplusplus  // use C-style linking
extern "C" {
#endif

// link to MKL libraries
extern int64_t __intel_mkl_feature_indicator;       // CPU feature bits
extern int64_t __intel_mkl_feature_indicator_x;     // CPU feature bits
void __intel_mkl_features_init();                   // unfair dispatcher: checks CPU features for Intel CPU's only
void __intel_mkl_features_init_x();                 // fair dispatcher: checks CPU features without discriminating by CPU brand

#ifdef __cplusplus
}  // end of extern "C"
#endif

void intel_mkl_patch() {
    // force a re-evaluation of the CPU features without discriminating by CPU brand
    __intel_mkl_feature_indicator = 0;
    __intel_mkl_feature_indicator_x = 0;
    __intel_mkl_features_init_x();
    __intel_mkl_feature_indicator = __intel_mkl_feature_indicator_x;
}
