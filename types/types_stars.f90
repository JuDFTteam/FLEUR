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
     REAL :: gmax=0.0
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
     INTEGER :: kmxxc_fft !number of g-vectors forming the nxc3_fft stars in the charge density or xc-density sphere

     INTEGER :: nxc3_fft ! number of stars in the  charge density  fft-box

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
     REAL, ALLOCATABLE:: phi2(:) !(n2d)
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
   CONTAINS
     PROCEDURE :: mpi_bc=>mpi_bc_stars
  END TYPE t_stars
CONTAINS
  SUBROUTINE mpi_bc_stars(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_stars),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank,myrank,ierr,n
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF

    CALL mpi_bc(this%gmax,rank,mpi_comm)
    CALL mpi_bc(this%ng3,rank,mpi_comm)
    CALL mpi_bc(this%ng2,rank,mpi_comm)
    CALL mpi_bc(this%mx1,rank,mpi_comm)
    CALL mpi_bc(this%mx2,rank,mpi_comm)
    CALL mpi_bc(this%mx3,rank,mpi_comm)
    CALL mpi_bc(this%kimax,rank,mpi_comm)
    CALL mpi_bc(this%kimax2,rank,mpi_comm)
    CALL mpi_bc(this%kq1_fft,rank,mpi_comm)
    CALL mpi_bc(this%kq2_fft,rank,mpi_comm)
    CALL mpi_bc(this%kq3_fft,rank,mpi_comm)
    CALL mpi_bc(this%kmxq_fft ,rank,mpi_comm)
    CALL mpi_bc(this%igq_fft,rank,mpi_comm)
    CALL mpi_bc(this%igq2_fft,rank,mpi_comm)
    CALL mpi_bc(this%kxc1_fft,rank,mpi_comm)
    CALL mpi_bc(this%kxc2_fft,rank,mpi_comm)
    CALL mpi_bc(this%kxc3_fft,rank,mpi_comm)
    CALL mpi_bc(this%ng3_fft,rank,mpi_comm)
    CALL mpi_bc(this%kmxxc_fft,rank,mpi_comm)
    CALL mpi_bc(this%nxc3_fft ,rank,mpi_comm)
    CALL mpi_bc(this%kv3,rank,mpi_comm)
    CALL mpi_bc(this%sk3,rank,mpi_comm)
    CALL mpi_bc(this%ig,rank,mpi_comm)
    CALL mpi_bc(this%nstr,rank,mpi_comm)
    CALL mpi_bc(this%kv2,rank,mpi_comm)
    CALL mpi_bc(this%sk2,rank,mpi_comm)
    CALL mpi_bc(this%nstr2,rank,mpi_comm)
    CALL mpi_bc(this%ig2,rank,mpi_comm)
    CALL mpi_bc(this%phi2,rank,mpi_comm)
    CALL mpi_bc(this%rgphs,rank,mpi_comm)
    CALL mpi_bc(this%igfft,rank,mpi_comm)
    CALL mpi_bc(this%igfft2,rank,mpi_comm)
    CALL mpi_bc(this%pgfft,rank,mpi_comm)
    CALL mpi_bc(this%pgfft2,rank,mpi_comm)
    CALL mpi_bc(this%ft2_gfx,rank,mpi_comm)
    CALL mpi_bc(this%ustep,rank,mpi_comm)
    CALL mpi_bc(this%ufft,rank,mpi_comm)


  END SUBROUTINE mpi_bc_stars
END MODULE m_types_stars
