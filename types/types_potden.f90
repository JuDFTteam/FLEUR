!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_potden

  !> Data type for the density or the potential
   TYPE t_potden
     INTEGER             :: iter  
     INTEGER             :: potdenType
     COMPLEX,ALLOCATABLE :: pw(:,:),pw_w(:,:)
     REAL,ALLOCATABLE    :: mt(:,:,:,:)
     REAL,ALLOCATABLE    :: vacz(:,:,:)
     COMPLEX,ALLOCATABLE :: vacxy(:,:,:,:)
     !For angles of density/potential in noco case
     REAL,ALLOCATABLE  :: theta_pw(:)
     REAL,ALLOCATABLE  :: phi_pw(:)
     REAL,ALLOCATABLE  :: theta_vacz(:,:)
     REAL,ALLOCATABLE  :: phi_vacz(:,:)
     REAL,ALLOCATABLE  :: theta_vacxy(:,:,:)
     REAL,ALLOCATABLE  :: phi_vacxy(:,:,:)
     REAL,ALLOCATABLE  :: theta_mt(:,:)
     REAL,ALLOCATABLE  :: phi_mt(:,:)
     

     ! For density matrix and associated potential matrix
     COMPLEX, ALLOCATABLE :: mmpMat(:,:,:,:)

     !this type contains two init routines that should be used to allocate
     !memory. You can either specify the datatypes or give the dimensions as integers
     !See implementation below!
   CONTAINS
     PROCEDURE :: init_potden_types
     PROCEDURE :: init_potden_simple
     PROCEDURE :: resetpotden
     GENERIC   :: init=>init_potden_types,init_potden_simple
     PROCEDURE :: copy_both_spin
     PROCEDURE :: sum_both_spin
     procedure :: SpinsToChargeAndMagnetisation
     procedure :: ChargeAndMagnetisationToSpins
     procedure :: addPotDen
     procedure :: subPotDen
     procedure :: distribute
     procedure :: collect
  END TYPE t_potden

CONTAINS
  subroutine collect(this,mpi_comm)
    use m_mpi_bc_tool
    implicit none
    class(t_potden),INTENT(INOUT) :: this
    integer :: mpi_comm
#ifdef CPP_MPI
    include 'mpif.h'
    real,ALLOCATABLE::rtmp(:)
    complex,ALLOCATABLE::ctmp(:)
    !pw
    ALLOCATE(ctmp(size(this%pw)))
    CALL MPI_REDUCE(this%pw,ctmp,size(this%pw),MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm,ierr)
    if (irank==0) this%pw=reshape(ctmp,shape(this%pw))
    deallocate(ctmp)
    !mt
    ALLOCATE(rtmp(size(this%mt)))
    CALL MPI_REDUCE(this%mt,rtmp,size(this%mt),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,ierr)
    if (irank==0) this%mt=reshape(rtmp,shape(this%mt))
    deallocate(rtmp)
    !vac
    if (allocated(this%vacz)) THEN
       ALLOCATE(rtmp(size(this%vacz)))
       CALL MPI_REDUCE(this%vacz,rtmp,size(this%vacz),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,ierr)
       if (irank==0) this%vacz=reshape(rtmp,shape(this%vacz))
       deallocate(rtmp)
       ALLOCATE(ctmp(size(this%vacxy)))
       CALL MPI_REDUCE(this%vacxy,ctmp,size(this%vacxy),MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm,ierr)
       if (irank==0) this%vacxy=reshape(ctmp,shape(this%vacxy))
       deallocate(ctmp)
    endif
    !density matrix
    if (allocated(this%mmpMat)) then
       ALLOCATE(ctmp(size(this%mmpMat)))
       CALL MPI_REDUCE(this%mmpMat,ctmp,size(this%mmpMat),MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm,ierr)
       if (irank==0) this%mmpMat=reshape(ctmp,shape(this%mmpMat))
       deallocate(ctmp)
    endif
#endif
  end subroutine collect
  
  subroutine distribute(this,mpi_comm)
    use m_mpi_bc_tool
    implicit none
    class(t_potden),INTENT(INOUT) :: this
    integer :: mpi_comm
#ifdef CPP_MPI
    include 'mpif.h'
    call mpi_bc(this%iter,0,mpi_comm)
    call mpi_bc(this%potdentype,0,mpi_comm)
    call mpi_bc(this%pw,0,mpi_comm)
    call mpi_bc(this%pw_w ,0,mpi_comm)
    call mpi_bc(this%mt ,0,mpi_comm)
    call mpi_bc(this%vacz,0,mpi_comm)
    call mpi_bc(this%vacxy,0,mpi_comm)
    call mpi_bc(this%theta_pw,0,mpi_comm)
    call mpi_bc(this%phi_pw,0,mpi_comm)
    call mpi_bc(this%theta_vacz,0,mpi_comm)
    call mpi_bc(this%phi_vacz,0,mpi_comm)
    call mpi_bc(this%theta_vacxy,0,mpi_comm)
    call mpi_bc(this%phi_vacxy,0,mpi_comm)
    call mpi_bc(this%theta_mt,0,mpi_comm)
    call mpi_bc(this%phi_mt,0,mpi_comm)
#endif
  end subroutine distribute
  
  SUBROUTINE sum_both_spin(this,that)
    IMPLICIT NONE
    CLASS(t_potden),INTENT(INOUT)   :: this
    TYPE(t_potden),INTENT(INOUT),OPTIONAL :: that

    IF (PRESENT(that)) THEN
       IF (SIZE(this%pw,2)>1) THEN
          that%mt(:,0:,:,1)=this%mt(:,0:,:,1)+this%mt(:,0:,:,2)
          that%pw(:,1)=this%pw(:,1)+this%pw(:,2)
          that%vacz(:,:,1)=this%vacz(:,:,1)+this%vacz(:,:,2)
          that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)+this%vacxy(:,:,:,2)
          IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)+this%pw_w(:,2)
       ELSE
          that%mt(:,0:,:,1)=this%mt(:,0:,:,1)
          that%pw(:,1)=this%pw(:,1)
          that%vacz(:,:,1)=this%vacz(:,:,1)
          that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)
          IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)
       ENDIF
    ELSE
       IF (SIZE(this%pw,2)>1) THEN
          this%mt(:,0:,:,1)=this%mt(:,0:,:,1)+this%mt(:,0:,:,2)
          this%pw(:,1)=this%pw(:,1)+this%pw(:,2)
          this%vacz(:,:,1)=this%vacz(:,:,1)+this%vacz(:,:,2)
          this%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)+this%vacxy(:,:,:,2)
          IF (ALLOCATED(this%pw_w)) this%pw_w(:,1)=this%pw_w(:,1)+this%pw_w(:,2)
       ENDIF
    END IF
  END SUBROUTINE sum_both_spin
    
  SUBROUTINE copy_both_spin(this,that)
    IMPLICIT NONE
    CLASS(t_potden),INTENT(IN)   :: this
    TYPE(t_potden),INTENT(INOUT) :: that

    that%mt(:,0:,:,1)=this%mt(:,0:,:,1)
    that%pw(:,1)=this%pw(:,1)
    that%vacz(:,:,1)=this%vacz(:,:,1)
    that%vacxy(:,:,:,1)=this%vacxy(:,:,:,1)
    IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,1)=this%pw_w(:,1)
    
    IF (SIZE(that%mt,4)>1) THEN
       that%mt(:,0:,:,2)=this%mt(:,0:,:,1)
       that%pw(:,2)=this%pw(:,1)
       that%vacz(:,:,2)=this%vacz(:,:,1)
       that%vacxy(:,:,:,2)=this%vacxy(:,:,:,1)
       IF (ALLOCATED(that%pw_w).AND.ALLOCATED(this%pw_w)) that%pw_w(:,2)=this%pw_w(:,1)
    END IF
  END SUBROUTINE copy_both_spin

  subroutine SpinsToChargeAndMagnetisation( den )
    implicit none
    class(t_potden), intent(inout)    :: den
    !type(t_potden),  intent(inout) :: charge_magn

    type(t_potden) :: copy

    copy = den

    den%mt(:,0:,:,  1) = copy%mt(:,0:,:,  1) + copy%mt(:,0:,:,  2)
    den%mt(:,0:,:,  2) = copy%mt(:,0:,:,  1) - copy%mt(:,0:,:,  2)
    den%pw(:,       1) = copy%pw(:,       1) + copy%pw(:,       2)
    den%pw(:,       2) = copy%pw(:,       1) - copy%pw(:,       2)
    den%vacz(:,:,   1) = copy%vacz(:,:,   1) + copy%vacz(:,:,   2)
    den%vacz(:,:,   2) = copy%vacz(:,:,   1) - copy%vacz(:,:,   2)
    den%vacxy(:,:,:,1) = copy%vacxy(:,:,:,1) + copy%vacxy(:,:,:,2)
    den%vacxy(:,:,:,2) = copy%vacxy(:,:,:,1) - copy%vacxy(:,:,:,2)

  end subroutine

  subroutine ChargeAndMagnetisationToSpins( den )
    implicit none
    class(t_potden), intent(inout)    :: den
    !type(t_potden),  intent(inout) :: spins

    type(t_potden) :: copy

    copy = den

    den%mt(:,0:,:,  1) = ( copy%mt(:,0:,:,  1) + copy%mt(:,0:,:,  2) ) / 2
    den%mt(:,0:,:,  2) = ( copy%mt(:,0:,:,  1) - copy%mt(:,0:,:,  2) ) / 2
    den%pw(:,       1) = ( copy%pw(:,       1) + copy%pw(:,       2) ) / 2
    den%pw(:,       2) = ( copy%pw(:,       1) - copy%pw(:,       2) ) / 2
    den%vacz(:,:,   1) = ( copy%vacz(:,:,   1) + copy%vacz(:,:,   2) ) / 2
    den%vacz(:,:,   2) = ( copy%vacz(:,:,   1) - copy%vacz(:,:,   2) ) / 2
    den%vacxy(:,:,:,1) = ( copy%vacxy(:,:,:,1) + copy%vacxy(:,:,:,2) ) / 2
    den%vacxy(:,:,:,2) = ( copy%vacxy(:,:,:,1) - copy%vacxy(:,:,:,2) ) / 2

  end subroutine

  subroutine addPotDen( PotDen3, PotDen1, PotDen2 )
    implicit none
    class(t_potden), intent(in)    :: PotDen1
    class(t_potden), intent(in)    :: PotDen2
    class(t_potden), intent(inout) :: PotDen3

    PotDen3%iter       = PotDen1%iter
    PotDen3%potdenType = PotDen1%potdenType
    PotDen3%mt         = PotDen1%mt + PotDen2%mt
    PotDen3%pw         = PotDen1%pw + PotDen2%pw
    PotDen3%vacz       = PotDen1%vacz + PotDen2%vacz
    PotDen3%vacxy      = PotDen1%vacxy + PotDen2%vacxy
    if( allocated( PotDen1%pw_w ) .and. allocated( PotDen2%pw_w ) .and. allocated( PotDen3%pw_w ) ) then
      PotDen3%pw_w = PotDen1%pw_w + PotDen2%pw_w
    end if
  
  end subroutine

  subroutine subPotDen( PotDen3, PotDen1, PotDen2 )
    implicit none
    class(t_potden), intent(in)    :: PotDen1
    class(t_potden), intent(in)    :: PotDen2
    class(t_potden), intent(inout) :: PotDen3
 
    PotDen3%iter       = PotDen1%iter
    PotDen3%potdenType = PotDen1%potdenType
    PotDen3%mt         = PotDen1%mt - PotDen2%mt
    PotDen3%pw         = PotDen1%pw - PotDen2%pw
    PotDen3%vacz       = PotDen1%vacz - PotDen2%vacz
    PotDen3%vacxy      = PotDen1%vacxy - PotDen2%vacxy
    if( allocated( PotDen1%pw_w ) .and. allocated( PotDen2%pw_w ) .and. allocated( PotDen3%pw_w ) ) then
      PotDen3%pw_w = PotDen1%pw_w - PotDen2%pw_w
    end if
 
  end subroutine

  SUBROUTINE init_potden_types(pd,stars,atoms,sphhar,vacuum,noco,jspins,potden_type)
    USE m_judft
    USE m_types_setup
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT):: pd 
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_stars),INTENT(IN) :: stars
    TYPE(t_sphhar),INTENT(IN):: sphhar
    TYPE(t_vacuum),INTENT(IN):: vacuum
    TYPE(t_noco),INTENT(IN)  :: noco
    INTEGER,INTENT(IN)       :: jspins, potden_type
 
    CALL init_potden_simple(pd,stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,&
         atoms%n_u,jspins,noco%l_noco,noco%l_mtnocopot,potden_type,&
         vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
  END SUBROUTINE init_potden_types

  SUBROUTINE init_potden_simple(pd,ng3,jmtd,nlhd,ntype,n_u,jspins,nocoExtraDim,nocoExtraMTDim,potden_type,nmzd,nmzxyd,n2d)
    USE m_constants
    USE m_judft
    IMPLICIT NONE
    CLASS(t_potden),INTENT(OUT) :: pd
    INTEGER,INTENT(IN)          :: ng3,jmtd,nlhd,ntype,n_u,jspins,potden_type
    LOGICAL,INTENT(IN)          :: nocoExtraDim,nocoExtraMTDim
    INTEGER,INTENT(IN)          :: nmzd,nmzxyd,n2d

    INTEGER:: err(4)

    err=0
    pd%iter=0
    pd%potdenType=potden_type
    IF(ALLOCATED(pd%pw)) DEALLOCATE (pd%pw)
    IF(ALLOCATED(pd%mt)) DEALLOCATE (pd%mt)
    IF(ALLOCATED(pd%vacz)) DEALLOCATE (pd%vacz)
    IF(ALLOCATED(pd%vacxy)) DEALLOCATE (pd%vacxy)
    IF(ALLOCATED(pd%mmpMat)) DEALLOCATE (pd%mmpMat)
    ALLOCATE (pd%pw(ng3,MERGE(3,jspins,nocoExtraDim)),stat=err(1))
    ALLOCATE (pd%mt(jmtd,0:nlhd,ntype,MERGE(4,jspins,nocoExtraMTDim)),stat=err(2))
    ALLOCATE (pd%vacz(nmzd,2,MERGE(4,jspins,nocoExtraDim)),stat=err(3))
    ALLOCATE (pd%vacxy(nmzxyd,n2d-1,2,MERGE(3,jspins,nocoExtraDim)),stat=err(4))

    ALLOCATE (pd%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,n_u),jspins))

    IF (ANY(err>0)) CALL judft_error("Not enough memory allocating potential or density")
    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
  END SUBROUTINE init_potden_simple
!!$#CPP_TODO_copy !code from brysh1,brysh2... 
!!$  SUBROUTINE get_combined_vector(input,stars,atoms,sphhar,noco,vacuum,sym,oneD,&
!!$                    den,nmap,nmaph,mapmt,mapvac2,sout) 
!!$    !This was brysh1 before
!!$    USE m_types
!!$    IMPLICIT NONE
!!$    TYPE(t_oneD),INTENT(IN)    :: oneD
!!$    TYPE(t_input),INTENT(IN)   :: input
!!$    TYPE(t_vacuum),INTENT(IN)  :: vacuum
!!$    TYPE(t_noco),INTENT(IN)    :: noco
!!$    TYPE(t_sym),INTENT(IN)     :: sym
!!$    TYPE(t_stars),INTENT(IN)   :: stars
!!$    TYPE(t_sphhar),INTENT(IN)  :: sphhar
!!$    TYPE(t_atoms),INTENT(IN)   :: atoms
!!$    TYPE(t_potden),INTENT(IN)  :: den
!!$
!!$    ! Scalar Arguments
!!$    INTEGER, INTENT (OUT) :: mapmt,mapvac2,nmap,nmaph
!!$
!!$    ! Array Arguments
!!$    REAL,ALLOCATABLE,INTENT (OUT) :: sout(:)
!!$
!!$    ! Local Scalars
!!$    INTEGER i,iv,j,js,k,l,n,na,nvaccoeff,nvaccoeff2,mapmtd
!!$
!!$    !Calculation of size
!!$    i=SIZE(den%mt)+MERGE(SIZE(den%pw),2*SIZE(den%pw),sym%invs)+SIZE(den%vacxz)+MERGE(SIZE(den%vacz)*2,SIZE(den%vacz),sym%invs)
!!$    IF (noco%l_mtnocopot.AND.sym%invs) i=i+
!!$    
!!$
!!$    
!!$    !--->  put input into arrays sout 
!!$    !      in the spin polarized case the arrays consist of 
!!$    !      spin up and spin down densities
!!$    
!!$    j=0
!!$    DO  js = 1,input%jspins
!!$       DO i = 1,stars%ng3
!!$          j = j + 1
!!$          sout(j) = REAL(den%pw(i,js))
!!$       END DO
!!$       IF (.NOT.sym%invs) THEN
!!$          DO i = 1,stars%ng3
!!$             j = j + 1
!!$             sout(j) = AIMAG(den%pw(i,js))
!!$          END DO
!!$       ENDIF
!!$       mapmt=0
!!$       na = 1
!!$       DO n = 1,atoms%ntype
!!$          DO l = 0,sphhar%nlh(atoms%ntypsy(na))
!!$             DO i = 1,atoms%jri(n)
!!$                mapmt = mapmt +1
!!$                j = j + 1
!!$                sout(j) = den%mt(i,l,n,js)
!!$             END DO
!!$          END DO
!!$          na = na + atoms%neq(n)
!!$       END DO
!!$       IF (input%film) THEN
!!$          DO iv = 1,vacuum%nvac
!!$             DO k = 1,vacuum%nmz
!!$                j = j + 1
!!$                sout(j) = den%vacz(k,iv,js)
!!$             END DO
!!$             DO k = 1,stars%ng2-1
!!$                DO i = 1,vacuum%nmzxy
!!$                   j = j + 1
!!$                   sout(j) =  REAL(den%vacxy(i,k,iv,js))
!!$                END DO
!!$             END DO
!!$             IF (.NOT.sym%invs2) THEN
!!$                DO k = 1,stars%ng2-1
!!$                   DO i = 1,vacuum%nmzxy
!!$                      j = j + 1
!!$                      sout(j) =  AIMAG(den%vacxy(i,k,iv,js))
!!$                   END DO
!!$                END DO
!!$             END IF
!!$          END DO
!!$       END IF
!!$       IF (js .EQ. 1) nmaph = j
!!$    ENDDO
!!$
!!$    mapvac2=0
!!$    IF (noco%l_noco) THEN
!!$       !--->    off-diagonal part of the density matrix
!!$       DO i = 1,stars%ng3
!!$          j = j + 1
!!$          sout(j) = REAL(den%pw(i,3))
!!$       END DO
!!$       DO i = 1,stars%ng3
!!$          j = j + 1
!!$          sout(j) = AIMAG(den%pw(i,3))
!!$       END DO
!!$       IF (input%film) THEN
!!$          DO iv = 1,vacuum%nvac
!!$             DO k = 1,vacuum%nmz
!!$                mapvac2 = mapvac2 + 1
!!$                j = j + 1
!!$                sout(j) = den%vacz(k,iv,3)
!!$             END DO
!!$             DO k = 1,stars%ng2-1
!!$                DO i = 1,vacuum%nmzxy
!!$                   mapvac2 = mapvac2 + 1
!!$                   j = j + 1
!!$                   sout(j) =  REAL(den%vacxy(i,k,iv,3))
!!$                END DO
!!$             END DO
!!$          END DO
!!$          DO iv = 1,vacuum%nvac
!!$             DO k = 1,vacuum%nmz
!!$                mapvac2 = mapvac2 + 1
!!$                j = j + 1
!!$                sout(j) = den%vacz(k,iv,4)
!!$             END DO
!!$             DO k = 1,stars%ng2-1
!!$                DO i = 1,vacuum%nmzxy
!!$                   mapvac2 = mapvac2 + 1
!!$                   j = j + 1
!!$                   sout(j) =  AIMAG(den%vacxy(i,k,iv,3))
!!$                END DO
!!$             END DO
!!$          END DO
!!$          nvaccoeff2 = 2*vacuum%nmzxy*(stars%ng2-1)*vacuum%nvac + 2*vacuum%nmz*vacuum%nvac
!!$          IF (mapvac2 .NE. nvaccoeff2) THEN
!!$             WRITE (6,*)'The number of vaccum coefficients off the'
!!$             WRITE (6,*)'off-diagonal part of the density matrix is'
!!$             WRITE (6,*)'inconsitent:'
!!$             WRITE (6,8000) mapvac2,nvaccoeff2
!!$8000         FORMAT ('mapvac2= ',i12,'nvaccoeff2= ',i12)
!!$             CALL juDFT_error("brysh1:# of vacuum coeff. inconsistent" ,calledby ="brysh1")
!!$          ENDIF
!!$       END IF
!!$    ENDIF ! noco
!!$
!!$    IF (atoms%n_u > 0 ) THEN     ! lda+U
!!$       DO js = 1,input%jspins
!!$          DO n = 1, atoms%n_u
!!$             DO k = -3, 3
!!$                DO i = -3, 3
!!$                   j = j + 1 
!!$                   sout(j) = REAL(den%mmpMat(i,k,n,js))
!!$                   j = j + 1 
!!$                   sout(j) = AIMAG(den%mmpMat(i,k,n,js))
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDIF
!!$
!!$    mapmtd = atoms%ntype*(sphhar%nlhd+1)*atoms%jmtd
!!$    IF (mapmt .GT. mapmtd) THEN
!!$       WRITE(6,*)'The number of mt coefficients is larger than the'
!!$       WRITE(6,*)'dimensions:'
!!$       WRITE (6,8040) mapmt,mapmtd
!!$8040   FORMAT ('mapmt= ',i12,' > mapmtd= ',i12)
!!$       CALL juDFT_error("brysh1: mapmt > mapmtd (dimensions)",calledby ="brysh1")
!!$    ENDIF
!!$
!!$    nmap = j
!!$    IF (nmap.GT.SIZE(sout)) THEN 
!!$       WRITE(6,*)'The total number of charge density coefficients is'
!!$       WRITE(6,*)'larger than the dimensions:'
!!$       WRITE (6,8030) nmap,SIZE(sout)
!!$8030   FORMAT ('nmap= ',i12,' > size(sout)= ',i12)
!!$       CALL juDFT_error("brysh1: nmap > mmap (dimensions)",calledby ="brysh1")
!!$    ENDIF
!!$
!!$  END SUBROUTINE get_combined_vector
!!$#endif
    
    
  
  SUBROUTINE resetPotDen(pd)

    IMPLICIT NONE

    CLASS(t_potden),INTENT(INOUT) :: pd

    pd%pw=CMPLX(0.0,0.0)
    pd%mt=0.0
    pd%vacz=0.0
    pd%vacxy=CMPLX(0.0,0.0)
    pd%mmpMat = CMPLX(0.0,0.0)
    IF (ALLOCATED(pd%pw_w)) DEALLOCATE(pd%pw_w)
  END SUBROUTINE resetPotDen

END MODULE m_types_potden
