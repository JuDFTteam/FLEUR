module m_types_mpbasis
   implicit none

   type t_mpbasis
      INTEGER, ALLOCATABLE   ::  gptm(:,:) ! (3, num_gpts)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER                ::  gptmd
   end type t_mpbasis
end module m_types_mpbasis
