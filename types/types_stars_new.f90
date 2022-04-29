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
     !phase factor of g-vector
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

  subroutine init_stars()

    gmax2=stars%gmax**2
    stars%ig=0
    associate(k=>stars%ng3)
    !Generate 3D stars
    x_dim: DO k1 = stars%mx1,-stars%mx1,-1
      kv(1) = k1
      y_dim: DO k2 = stars%mx2,-stars%mx2,-1
        kv(2) = k2
        z_dim: DO k3 = stars%mx3,-stars%mx3,-1
          IF ( stars%ig(k1,k2,k3) .NE. 0 ) CYCLE z_dim  ! belongs to another star
          kv(3) = k3

          g=matmul(kv,cell%bmat)
          s=dot_product(g,g)
          if (s>gmax2) cycle z_dim !not in sphere
          k=k+1
          stars%kv3(:,k)=kv
          stars%sk3(k)=sqrt(s)
          ! secondary key for equal length stars
          gsk3(k) = (stars%mx1+stars%kv3(1,k)) +&
            &           (stars%mx2+stars%kv3(2,k))*(2*stars%mx1+1) +&
            &           (stars%mx3+stars%kv3(3,k))*(2*stars%mx1+1)*(2*stars%mx2+1)
          !Now generate all equivalent g-vectors
          CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),&
                             kr,phas)
          DO n = 1,nop
            stars%ig(kr(1,n),kr(2,n),kr(2,n))=k
          ENDDO
        ENDDO z_dim
      ENDDO y_dim  
    ENDDO x_dim
    end associate
    !sort for increasing length sk3
    CALL sort(index,stars%sk3,gsk3)
    stars%kv3=stars%kv3(:,index)
    stars%sk3=stars%sk3(:,index)
    ! set up the pointers and phases for 3d stars
    DO  k = 1,nq3
      CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),&
      kr,phas)
      symloop: DO n = 1,nop    
        DO n1 = 1,n-1
          if (all(kr(:,n)==kr(:,n1))) cycle symloop !This g-vector we have already
        ENDDO
        stars%ig(kr(1,n),kr(2,n),kr(2,n))=k
        stars%rgph(kr(1,n),kr(2,n),kr(2,n))=stars%rgph(kr(1,n),kr(2,n),kr(2,n))+phas(n)

END MODULE m_types_stars
