PROGRAM test
  USE xc_f90_lib_m
  TYPE(xc_f90_func_t) :: xc_func
  CALL xc_f90_func_init(xc_func, 2, XC_UNPOLARIZED)	  
END PROGRAM test
