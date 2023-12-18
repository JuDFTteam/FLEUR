program test
  use MKL_DFT_TYPE
  use MKL_test
  type(dfti_descriptor),pointer :: dfti_handle
  integer :: dfti_status
  dfti_status = DftiCommitDescriptor(dfti_handle)
end program test
