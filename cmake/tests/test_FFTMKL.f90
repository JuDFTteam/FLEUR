MODULE  MKL_DFT_TYPE
  TYPE, PUBLIC :: DFTI_DESCRIPTOR
     PRIVATE
     INTEGER :: dontuse
     ! Structure of this type is not used in Fortran code
     ! the pointer to this type is used only
  END TYPE DFTI_DESCRIPTOR
END MODULE MKL_DFT_TYPE
MODULE MKL_TEST
INTERFACE DftiCommitDescriptor
           
   FUNCTION dfti_commit_descriptor_external(desc)
     USE MKL_DFT_TYPE
     !DEC$ ATTRIBUTES C :: dfti_commit_descriptor_external
     !DEC$ ATTRIBUTES REFERENCE :: dfti_commit_descriptor_external
     INTEGER dfti_commit_descriptor_external
     TYPE(DFTI_DESCRIPTOR), POINTER :: desc
   END FUNCTION dfti_commit_descriptor_external
   
END INTERFACE
END MODULE

program test
  use MKL_DFT_TYPE
  use MKL_test
  type(dfti_descriptor),pointer :: dfti_handle
  integer :: dfti_status
  dfti_status = DftiCommitDescriptor(dfti_handle)
end program test
