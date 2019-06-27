!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_stars
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_stars
  !The stars
  TYPE:: t_stars
      !max-length of star
      REAL :: gmax
      !no of 3d-stars
      INTEGER :: ng3
      !no of 2d-stars
      INTEGER :: ng2
      !dim of box
      INTEGER ::mx1
      INTEGER ::mx2
      INTEGER ::mx3
      !No of elements in FFT
      INTEGER ::kimax
      !No of elements in 2D-FFT
      INTEGER ::kimax2

      !Box for FFT in pwden
      INTEGER :: kq1_fft
      INTEGER :: kq2_fft
      INTEGER :: kq3_fft
      INTEGER :: kmxq_fft !no of g-vectors in sphere
      INTEGER, ALLOCATABLE :: igq_fft(:)
      INTEGER, ALLOCATABLE :: igq2_fft(:)

      !fft box for xc-pot
      INTEGER :: kxc1_fft
      INTEGER :: kxc2_fft
      INTEGER :: kxc3_fft

      INTEGER :: ng3_fft
      INTEGER :: kmxxc_fft !<number of g-vectors forming the nxc3_fft stars in the charge density or xc-density sphere

      INTEGER :: nxc3_fft !< number of stars in the  charge density  fft-box

      !rep. g-vector of star
      INTEGER, ALLOCATABLE ::kv3(:, :)
      !length of star
      REAL, ALLOCATABLE    ::sk3(:)
      !mapping of g-vectors to stars
      INTEGER, ALLOCATABLE ::ig(:, :, :)
      !No of g-vectors in star
      INTEGER, ALLOCATABLE ::nstr(:)
      !rep. g-vector of 2D-star
      INTEGER, ALLOCATABLE ::kv2(:, :)
      !length of 2D-star
      REAL, ALLOCATABLE    ::sk2(:)
      !No of g-vecs in 2D-star
      INTEGER, ALLOCATABLE ::nstr2(:)
      !mapping of
      INTEGER, ALLOCATABLE ::ig2(:)
      !
      REAL, ALLOCATABLE:: phi2(:) !<(n2d)
      !phase phactor of g-vector
      COMPLEX, ALLOCATABLE    ::rgphs(:, :, :)
      !mapping of stars to FFT-box
      INTEGER, ALLOCATABLE :: igfft(:, :)
      !same for 2D
      INTEGER, ALLOCATABLE :: igfft2(:, :)
      !phasefactors for mapping
      COMPLEX, ALLOCATABLE  :: pgfft(:)
      !same of 2D
      COMPLEX, ALLOCATABLE  :: pgfft2(:)
      !
      REAL, ALLOCATABLE     :: ft2_gfx(:), ft2_gfy(:)
      COMPLEX, ALLOCATABLE :: ustep(:)
      REAL, ALLOCATABLE    :: ufft(:)
    contains
      procedure :: init
   END TYPE t_stars

   subroutine init(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot,oneD,mpi)
     class(t_stars),intent(INOUT) :: stars
     type(t_sym),intent(in)::sym
     type(t_atoms),intent(in)::atoms
     type(t_vacuum),intent(in)::vacuum
     type(t_sphhar),intent(in)::sphhar
     type(t_input),intent(in)::input
     type(t_cell),intent(in)::cell
     type(t_xcpot),intent(in)::xcpot
     type(t_oneD),intent(in)::oneD
     type(t_mpi),intent(in)::mpi
     ! Generate stars
  
  CALL timestart("strgn") 
  IF (input%film) THEN
     CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
     IF (oneD%odd%d1) THEN
        CALL od_strgn1(xcpot,cell,sym,oneD)
     END IF
  ELSE
     CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
  END IF
  
  ALLOCATE (stars%igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))
  ALLOCATE (stars%igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1))
  
  ! Set up pointer for backtransformation from g-vector in positive 
  ! domain of carge density fftibox into stars
  CALL prp_qfft_map(stars,sym,input,stars%igq2_fft,stars%igq_fft)
  CALL prp_qfft(stars,cell,noco,input)
  
  CALL timestop("strgn") 

  CALL timestart("stepf") 
  CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,mpi)
  CALL timestop("stepf") 
  

 END MODULE m_types_stars
