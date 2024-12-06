#ifndef MAGMA_HEADER_INCLUDED
#define MAGMA_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define MAGMA_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define MAGMA_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define MAGMA_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define MAGMA_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_##name##_

#endif
