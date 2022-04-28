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
     
   
     INTEGER :: ng3_fft !number of stars in fft-box of size 2*rkmax
   
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
     INTEGER, ALLOCATABLE ::i2g(:,:) !map g_x,g_y to 2D-star
     INTEGER, ALLOCATABLE ::ig2(:)   !map 3D-star to 2D-star
     INTEGER, ALLOCATABLE :: igvac(:,:)  ! map 2D-star and g_z to 3D-star

  !
     REAL, ALLOCATABLE:: phi2(:) !for oneD code, currently not in use
     !phase factor of g-vector
     COMPLEX, ALLOCATABLE    ::rgphs(:, :, :)
     !phase factor of 2d-g vectors
     COMPLEX, ALLOCATABLE    ::r2gphs(:, :)

     COMPLEX, ALLOCATABLE :: ustep(:)
     REAL, ALLOCATABLE    :: ufft(:)
   CONTAINS
     PROCEDURE :: mpi_bc=>mpi_bc_stars
     PROCEDURE :: init=>init_stars
     PROCEDURE :: dim=>dim_stars
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
    !CALL mpi_bc(this%kimax,rank,mpi_comm)
    !CALL mpi_bc(this%kimax2,rank,mpi_comm)
    !CALL mpi_bc(this%kq1_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kq2_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kq3_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kmxq_fft ,rank,mpi_comm)
    !CALL mpi_bc(this%igq_fft,rank,mpi_comm)
    !CALL mpi_bc(this%igq2_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kxc1_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kxc2_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kxc3_fft,rank,mpi_comm)
    CALL mpi_bc(this%ng3_fft,rank,mpi_comm)
    !CALL mpi_bc(this%kmxxc_fft,rank,mpi_comm)
    !CALL mpi_bc(this%nxc3_fft ,rank,mpi_comm)
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
    CALL mpi_bc(this%r2gphs,rank,mpi_comm)

    CALL mpi_bc(this%ustep,rank,mpi_comm)
    CALL mpi_bc(this%ufft,rank,mpi_comm)


  END SUBROUTINE mpi_bc_stars

  subroutine init_stars(stars,cell,sym,film,rkmax)
    USE m_spgrot
    USE m_types_cell
    USE m_types_sym
    USE m_sort
    CLASS(t_stars),INTENT(INOUT)  :: stars 
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_sym),INTENT(IN)        :: sym 
    LOGICAL,INTENT(IN)            :: film
    REAL,INTENT(IN)               :: rkmax

    INTEGER :: k1,k2,k3,n,n1,k
    REAL    :: s,g(3),gmax2
    INTEGER :: kr(3,sym%nop),kv(3)
    COMPLEX :: phas(sym%nop)
    INTEGER,ALLOCATABLE :: index(:)
    REAL,ALLOCATABLE :: gsk3(:)

    gmax2=stars%gmax**2
    allocate(gsk3(stars%ng3),index(stars%ng3))

    ALLOCATE(stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
    ALLOCATE(stars%kv3(3,stars%ng3),stars%sk3(stars%ng3),stars%nstr(stars%ng3))
    ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
    
    stars%rgphs=0.0
    stars%ig=0
    stars%nstr=0
    
    k=0
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
                             kr)
          DO n = 1,sym%nop
            stars%ig(kr(1,n),kr(2,n),kr(3,n))=k
          ENDDO
        ENDDO z_dim
      ENDDO y_dim  
    ENDDO x_dim
    if (k.ne.stars%ng3) call judft_error("BUG inconsistency in star setup")

    !sort for increasing length sk3
    CALL sort(index,stars%sk3,gsk3)
    stars%kv3=stars%kv3(:,index)
    stars%sk3=stars%sk3(index)
    ! set up the pointers and phases for 3d stars
    DO  k = 1,stars%ng3
      CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),&
      kr,phas)
      symloop: DO n = 1,sym%nop    
        stars%ig(kr(1,n),kr(2,n),kr(3,n))=k
        stars%rgphs(kr(1,n),kr(2,n),kr(3,n))=stars%rgphs(kr(1,n),kr(2,n),kr(3,n))+phas(n)
        DO n1 = 1,n-1
          if (all(kr(:,n)==kr(:,n1))) cycle symloop !This g-vector we have already
        ENDDO
        stars%nstr(k)=stars%nstr(k)+1
      ENDDO symloop
    ENDDO    
    !Adjust phases
    if (sym%symor) THEN
      stars%rgphs=1.0
    ELSE
      DO k1 = stars%mx1,-stars%mx1,-1
        DO k2 = stars%mx2,-stars%mx2,-1
          DO k3 = stars%mx3,-stars%mx3,-1
            IF ( stars%ig(k1,k2,k3)==0 ) CYCLE 
            stars%rgphs(k1,k2,k3)=stars%rgphs(k1,k2,k3)*stars%nstr(stars%ig(k1,k2,k3))/sym%nop
          enddo
        ENDDO
      ENDDO
    ENDIF    
    !count number of stars in 2*rkmax (stars are ordered)
    associate(i=>stars%ng3_fft)
      DO i=stars%ng3,1,-1
        IF ( stars%sk3(i).LE.2.0*rkmax ) EXIT
      ENDDO
    end associate 

    if (.not.film) return
    !
    ! Now the same for the 2D stars
    !
    ALLOCATE(stars%r2gphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2))
    ALLOCATE(stars%kv2(2,stars%ng2),stars%sk2(stars%ng2),stars%nstr2(stars%ng2))
    ALLOCATE(stars%i2g(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2))
    ALLOCATE(stars%ig2(stars%ng3))
    ALLOCATE(stars%igvac(stars%ng2,-stars%mx3:stars%mx3))
    stars%r2gphs=0.0 
    stars%i2g=0
    stars%nstr2=0
    kv(3)=0
    k=0
    !Generate 2D stars
    x_dim2: DO k1 = stars%mx1,-stars%mx1,-1
      kv(1) = k1
      y_dim2: DO k2 = stars%mx2,-stars%mx2,-1
        kv(2) = k2
        IF ( stars%i2g(k1,k2) .NE. 0 ) CYCLE y_dim2  ! belongs to another star
        g(:2)=matmul(kv(:2),cell%bmat(:2,:2))
        s=dot_product(g(:2),g(:2))
        if (s>gmax2) cycle y_dim2 !not in sphere
        k=k+1
        stars%kv2(:,k)=kv(:2)
        stars%sk2(k)=sqrt(s)
        ! secondary key for equal length stars
        gsk3(k) = (stars%mx1+stars%kv3(1,k)) +(stars%mx2+stars%kv3(2,k))*(2*stars%mx1+1) 
        !Now generate all equivalent g-vectors
        CALL spgrot(sym%nop2,sym%symor,sym%mrot(:2,:2,:),sym%tau(:2,:),sym%invtab,stars%kv2(:2,k),kr(:2,:))
        DO n = 1,sym%nop2
            stars%i2g(kr(1,n),kr(2,n))=k
        ENDDO
     
      ENDDO y_dim2  
    ENDDO x_dim2
    if (k.ne.stars%ng2) call judft_error("BUG in init_stars: inconsistency in ng2")
    !sort for increasing length sk2
    CALL sort(index(:stars%ng2),stars%sk2,gsk3(:stars%ng2))
    stars%kv2(:,:)=stars%kv2(:,index(:stars%ng2))
    stars%sk2=stars%sk2(index(:stars%ng2))
    ! set up the pointers and phases for 2d stars
    DO  k = 1,stars%ng2
      DO k3= stars%mx3,-stars%mx3,-1
         stars%igvac(k,k3) = stars%ig(stars%kv2(1,k),stars%kv2(2,k),k3)
      ENDDO  
      CALL spgrot(sym%nop2,sym%symor,sym%mrot(:2,:2,:),sym%tau(:2,:),sym%invtab,stars%kv2(:2,k),kr(:2,:),phas)
      symloop2: DO n = 1,sym%nop2    
        stars%i2g(kr(1,n),kr(2,n))=k
        stars%r2gphs(kr(1,n),kr(2,n))=stars%r2gphs(kr(1,n),kr(2,n))+phas(n)
        DO n1 = 1,n-1
          if (all(kr(:2,n)==kr(:2,n1))) cycle symloop2 !This g-vector we have already
        ENDDO
        stars%nstr2(k)=stars%nstr2(k)+1
      ENDDO symloop2
    ENDDO
    DO k=1,stars%ng3
      stars%ig2(k)=stars%i2g(stars%kv3(1,k),stars%kv3(2,k))
    ENDDO    
    !Adjust phases
    IF (sym%symor) THEN
      stars%r2gphs=1.0
    ELSE
      DO k1 = stars%mx1,-stars%mx1,-1
        DO k2 = stars%mx2,-stars%mx2,-1
          IF ( stars%i2g(k1,k2)==0 ) CYCLE 
          stars%r2gphs(k1,k2)=stars%r2gphs(k1,k2)*stars%nstr2(stars%i2g(k1,k2))/sym%nop2
        ENDDO
      ENDDO
    ENDIF     
  END SUBROUTINE init_stars

  subroutine dim_stars(stars,sym,cell,film)
    !! determine the key dimensions of the stars:
    !! mx1,mx2,mx3
    !! ng3, ng2
    USE m_spgrot
    USE m_types_cell
    USE m_types_sym
    USE m_boxdim
    CLASS(t_stars),INTENT(INOUT)  :: stars 
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_sym),INTENT(IN)        :: sym 
    LOGICAL,INTENT(IN)            :: film
    
    INTEGER :: k1,k2,k3,n,kv(3),kr(3,sym%nop)
    REAL    :: s,g(3),gmax2
    REAL    :: arltv1,arltv2,arltv3
   
    gmax2=stars%gmax**2
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    stars%mx1 = int(stars%gmax/arltv1) + 1
    stars%mx2 = int(stars%gmax/arltv2) + 1
    stars%mx3 = int(stars%gmax/arltv3) + 1

    ALLOCATE(stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
    ALLOCATE(stars%i2g(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2))
    
   
    
    stars%i2g=0
    stars%ig=0
    stars%ng2=0
    stars%ng3=0
    x_dim: DO k1 = stars%mx1,-stars%mx1,-1
      kv(1) = k1
      y_dim: DO k2 = stars%mx2,-stars%mx2,-1
        kv(2) = k2
        kv(3) = 0
        !Check 2d-star
        IF (stars%i2g(k1,k2)==0) THEN
            g=matmul(kv,cell%bmat)
            s=dot_product(g,g)
            IF (.not.s>gmax2) THEN !in sphere
              stars%ng2=stars%ng2+1
              CALL spgrot(sym%nop2,sym%symor,sym%mrot,sym%tau,sym%invtab,kv,kr)
              DO n = 1,sym%nop2
                if (kr(3,n).ne.0) cycle
                stars%i2g(kr(1,n),kr(2,n))=stars%ng2
              ENDDO
            ENDIF
        ENDIF      
        z_dim: DO k3 = stars%mx3,-stars%mx3,-1
          IF ( stars%ig(k1,k2,k3) .NE. 0 ) CYCLE z_dim  ! belongs to another star
          kv(3) = k3
          g=matmul(kv,cell%bmat)
          s=dot_product(g,g)
          if (s>gmax2) cycle z_dim !not in sphere
          stars%ng3=stars%ng3+1
          CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,kv(:),kr)
          DO n = 1,sym%nop
            stars%ig(kr(1,n),kr(2,n),kr(3,n))=stars%ng3
          ENDDO
        ENDDO z_dim
      ENDDO y_dim
    ENDDO x_dim

    if (.not.film) stars%ng2=0
    DEALLOCATE(stars%ig)
    DEALLOCATE(stars%i2g)
    
   
  END SUBROUTINE dim_stars
END MODULE m_types_stars
