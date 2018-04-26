PROGRAM test
  USE xc_f03_lib_m
  TYPE(xc_f03_func_t) :: xc_func
  CALL xc_f03_func_init(xc_func, 2, XC_UNPOLARIZED)	  
END PROGRAM test
