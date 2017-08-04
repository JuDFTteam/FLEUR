!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_corespec

  implicit none

! PARAMETERS

  complex, parameter :: cone = cmplx(1.d0,0.d0)
  complex, parameter :: cimu = cmplx(0.d0,1.d0)
  real, parameter :: alpha = 7.29735257d-3
  real, parameter :: mec2 = 0.51099891d6
  real, parameter :: ecoredeep = 0.5d0

  integer, parameter :: edgel(11) = (/0,1,1,2,2,3,3,4,4,5,5/)
  integer, parameter :: edgej(11) = (/1,1,3,3,5,5,7,7,9,9,11/)
  integer, parameter :: sign(4) = (/1,-1,0,0/)

  character(len=2), parameter :: ssep="> "
  character(len=3), parameter :: fsos1="t55"
  character(len=3), parameter :: fsos2="t90"
  character(len=15), parameter :: fsb="(a,('"//ssep//"'),a,"//fsos1
  character(len=9), parameter :: fse=fsos2//",2x,a)"
  character(len=11), parameter :: csmsgerr=" ... STOP !"
  character(len=14), parameter :: csmsgwar=" ... WARNING !"
  character(len=32), parameter :: csmsgs = fsb//")"
  character(len=64), parameter :: csmsgsss = fsb//",a,a5,"//fse
  character(len=64), parameter :: csmsgsis = fsb//",a,i5,"//fse
  character(len=64), parameter :: csmsgsisis = fsb//",a,i5,a,i5,"//fse
  character(len=64), parameter :: csmsgsfs = fsb//",a,f8.3,"//fse
  character(len=64), parameter :: csmsgses = fsb//",a,es12.3,"//fse

! INPUT PARAMETERS

  type csitype
     sequence
     integer :: verb  ! output verbosity
     integer :: type  ! atomic type used for calculation of core spectra
     character(len=1) :: edge  ! edge character (K,L,M,N,O,P)
     integer :: edgeidx(11)  ! l-j edges
     integer :: lx  ! maximum lmax considered in spectra calculation
     real :: ek0  ! kinetic energy of incoming electrons
     real :: emn  ! energy spectrum lower bound
     real :: emx  ! energy spectrum upper bound
     real :: ein  ! energy spectrum increment
  end type csitype

! VARIABLES

  logical :: l_cs

  integer :: l1,l2,la1,la2,li
  integer :: m1,m2,mu1,mu2,mi
  integer :: lx,ln,lax,lan,lix,lin

  character(len=32) :: smeno

  type (csitype) :: csi

  type csvtype
     sequence
     integer :: nc  ! main quantum no. of the core level of atomic type
     integer :: nljc  ! number of l-j edge lines
     integer, allocatable :: lc(:)  ! edge angular quantum nos.; nljc elements
     real, allocatable :: eedge(:)  ! lc-dep. edge energy; nljc elements
     real, allocatable :: occ(:)  ! lc-dep. occupation; nljc elements
     integer :: nex  ! no. of energy sampling points
     real, allocatable :: egrid(:)  ! energy grid; 0:nex elements
     real, allocatable :: eos(:)  ! energy grid / sigma
     real, allocatable :: eloss(:,:)  ! efermi-eedge+egrid
     integer :: nen  ! minimum index for which egrid >=0
     integer :: nqv  ! no. of q vectors
     real :: qv0  ! |q| of incoming electrons
     real, allocatable :: qv1(:,:,:)  ! |q| of outgoing electrons
     real, allocatable :: qv(:,:,:,:)  ! delta q vectors
     real :: gamma  ! gamma = 1+ek0/mc2
     real :: beta  ! beta = v/c = 1/sqrt(1-1/gamma^2)
     real, allocatable :: gaunt(:,:,:,:,:,:)  ! gaunt coefficients
     real, allocatable :: fc(:,:,:,:)  ! core radial function
     real, allocatable :: fv(:,:,:,:)  ! valence radial function
     real, allocatable :: fb(:,:,:,:,:)  ! bessel function
     real, allocatable :: rmeA(:,:,:,:,:,:,:)  ! matrix elements
     real, allocatable :: rmeB(:,:,:,:,:,:,:)  ! matrix elements
     real, allocatable :: rmeC(:,:,:,:,:,:,:)  ! matrix elements
     real, allocatable :: dose(:,:,:,:,:)  ! dos (bands)
     real, allocatable :: dosb(:,:,:,:,:)  ! dos (bands)
     complex, allocatable :: ddscs(:,:,:,:,:)  ! dos (bands)
     
  end type csvtype

  type (csvtype) :: csv

end module m_corespec
