module m_types_mpbasis
   implicit none

   type t_mpbasis
      INTEGER, ALLOCATABLE   ::  gptm(:,:) ! (3, num_gpts)
      INTEGER, ALLOCATABLE   ::  ngptm(:)
      INTEGER, ALLOCATABLE   ::  gptm_ptr(:,:)
      REAL                   ::  g_cutoff
   CONTAINS
      procedure :: num_gpts => mpbasis_num_gpts
   end type t_mpbasis
contains
   function mpbasis_num_gpts(mpbasis)
      implicit NONE
      class(t_mpbasis), intent(in) :: mpbasis
      integer    :: mpbasis_num_gpts

      mpbasis_num_gpts = size(mpbasis%gptm,dim=2)
   end function mpbasis_num_gpts
end module m_types_mpbasis
