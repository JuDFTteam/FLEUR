!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mixvector
  !TODO!!!
  ! LDA+U
  ! Noco (third spin)
  
  use m_types
  implicit none
#ifdef CPP_MPI
      include 'mpif.h'
#endif
  PRIVATE
  !Here we store the pointers used for metric
  TYPE(t_stars),POINTER  :: stars
  TYPE(t_cell),POINTER   :: cell
  TYPE(t_sphhar),POINTER :: sphhar
  TYPE(t_atoms),POINTER  :: atoms  =>null()
  TYPE(t_sym),POINTER    :: sym  =>null()
  INTEGER                :: jspins,nvac
  LOGICAL                :: l_noco,l_mtnocopot
  INTEGER                :: pw_length !The shape of the local arrays
  INTEGER                :: pw_start(3)=0,pw_stop(3) !First and last index for spin
  INTEGER                :: mt_length,mt_length_g
  INTEGER                :: mt_start(3)=0,mt_stop(3) !First and last index for spin
  INTEGER                :: vac_length,vac_length_g
  INTEGER                :: vac_start(3)=0,vac_stop(3) !First and last index for spin
  INTEGER                :: misc_length=0,misc_length_g
  INTEGER                :: misc_start(3)=0,misc_stop(3) !First and last index for spin
  INTEGER                :: mix_mpi_comm !Communicator for all PEs doing mixing
  LOGICAL                :: spin_here(3)=.TRUE.
  LOGICAL                :: pw_here=.TRUE.
  LOGICAL                :: mt_here=.TRUE.
  LOGICAL                :: vac_here=.TRUE.
  LOGICAL                :: misc_here=.TRUE.
  INTEGER                :: mt_rank=0
  INTEGER                :: mt_size=1
  LOGICAL                :: l_pot=.FALSE. !Is this a potential?
  REAL,ALLOCATABLE       :: g_mt(:),g_vac(:),g_misc(:)
  
  TYPE,PUBLIC:: t_mixvector
     REAL,ALLOCATABLE       :: vec_pw(:)
     REAL,ALLOCATABLE       :: vec_mt(:)
     REAL,ALLOCATABLE       :: vec_vac(:)
     REAL,ALLOCATABLE       :: vec_misc(:)
   CONTAINS
     procedure :: alloc=>mixvector_alloc
     PROCEDURE :: from_density=>mixvector_from_density
     PROCEDURE :: to_density=>mixvector_to_density
     PROCEDURE :: apply_metric=>mixvector_metric
     PROCEDURE :: multiply_dot_mask
     PROCEDURE :: read_unformatted
     PROCEDURE :: write_unformatted
     GENERIC :: READ(UNFORMATTED) =>read_unformatted
     GENERIC :: WRITE(UNFORMATTED) =>write_unformatted
  END TYPE t_mixvector

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE multiply_scalar
     MODULE PROCEDURE multiply_scalar_spin
  END INTERFACE OPERATOR (*)
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_vectors
  END INTERFACE OPERATOR (+)
  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_vectors
  END INTERFACE OPERATOR (-)
  INTERFACE OPERATOR (.dot.)
     MODULE PROCEDURE multiply_dot
  END INTERFACE OPERATOR (.dot.)

  PUBLIC :: OPERATOR(+),OPERATOR(-),OPERATOR(*),OPERATOR(.dot.)
  PUBLIC :: mixvector_init,mixvector_reset

CONTAINS

  SUBROUTINE READ_unformatted(this,unit,iostat,iomsg)
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(INOUT)::this
    INTEGER,INTENT(IN)::unit
    INTEGER,INTENT(OUT)::iostat
    CHARACTER(len=*),INTENT(INOUT)::iomsg

    CALL this%alloc()
    IF (pw_here) READ(unit) this%vec_pw
    IF (mt_here) READ(unit) this%vec_mt
    IF (vac_here) READ(unit) this%vec_vac
    IF (misc_here) READ(unit) this%vec_misc
  END SUBROUTINE READ_unformatted

  SUBROUTINE write_unformatted(this,unit,iostat,iomsg)
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)::this
    INTEGER,INTENT(IN)::unit
    INTEGER,INTENT(OUT)::iostat
    CHARACTER(len=*),INTENT(INOUT)::iomsg
    IF (pw_here) WRITE(unit) this%vec_pw
    IF (mt_here) WRITE(unit) this%vec_mt
    IF (vac_here) WRITE(unit) this%vec_vac
    IF (misc_here) WRITE(unit) this%vec_misc
  END SUBROUTINE write_unformatted

  


  SUBROUTINE mixvector_reset()
    IMPLICIT NONE
    atoms=>NULL()
    sym=>NULL()
    IF (ALLOCATED(g_mt)) DEALLOCATE(g_mt)
    IF (ALLOCATED(g_vac)) DEALLOCATE(g_vac)
    IF (ALLOCATED(g_misc)) DEALLOCATE(g_misc)
  END SUBROUTINE mixvector_reset

  
  SUBROUTINE mixvector_from_density(vec,den,swapspin)
    USE m_types
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(INOUT)    :: vec
    TYPE(t_potden),    INTENT(inout)    :: Den
    LOGICAL,INTENT(IN),OPTIONAL         :: swapspin
    INTEGER:: js,ii,n,l,iv,j
    call den%distribute(mix_mpi_comm)
    DO js=1,MERGE(jspins,3,.NOT.l_noco)
       j=js
       IF (PRESENT(swapspin)) THEN
          IF (swapspin.AND.js<3) j=MERGE(1,2,js==2)
       ENDIF
       IF (spin_here(js)) THEN
          !PW part
          IF (pw_here) THEN
             vec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1)=REAL(den%pw(:,j))
             IF ((.NOT.sym%invs).or.(js==3)) THEN
                vec%vec_pw(pw_start(js)+stars%ng3:pw_start(js)+2*stars%ng3-1)=AIMAG(den%pw(:,j))
             ENDIF
          ENDIF
          IF (vac_here) THEN
             !This PE stores vac-data
             ii=vac_start(js)-1
             DO iv=1,nvac
                vec%vec_vac(ii+1:ii+SIZE(den%vacz,1))=den%vacz(:,iv,j)
                ii=ii+SIZE(den%vacz,1)
                vec%vec_vac(ii+1:ii+SIZE(den%vacxy(:,:,iv,js)))=RESHAPE(REAL(den%vacxy(:,:,iv,j)),(/SIZE(den%vacxy(:,:,iv,j))/))
                ii=ii+SIZE(den%vacxy(:,:,iv,j))
                IF ((.NOT.sym%invs2).or.(js==3))THEN
                   vec%vec_vac(ii+1:ii+SIZE(den%vacxy(:,:,iv,j)))=RESHAPE(AIMAG(den%vacxy(:,:,iv,j)),(/SIZE(den%vacxy(:,:,iv,j))/))
                   ii=ii+SIZE(den%vacxy(:,:,iv,j))
                ENDIF
                IF (js>2)THEN
                   vec%vec_vac(ii+1:ii+SIZE(den%vacz,1))=den%vacz(:,iv,4)
                   ii=ii+SIZE(den%vacz,1)
                ENDIF
             ENDDO
          ENDIF
          IF (js>2.AND..NOT.l_mtnocopot) RETURN
          IF (mt_here) THEN
             !This PE stores some(or all) MT data
             ii=mt_start(js)-1
             DO n=mt_rank+1,atoms%ntype,mt_size
                DO l=0,sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))
                   vec%vec_mt(ii+1:ii+atoms%jri(n))=den%mt(:atoms%jri(n),l,n,j)
                   ii=ii+atoms%jri(n)
                ENDDO
             ENDDO
             IF (js==3) THEN !Imaginary part 
                DO n=mt_rank+1,atoms%ntype,mt_size
                   DO l=0,sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))
                      vec%vec_mt(ii+1:ii+atoms%jri(n))=den%mt(:atoms%jri(n),l,n,4)
                      ii=ii+atoms%jri(n)
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
          IF (js>2) RETURN
          IF (misc_here) THEN
             vec%vec_misc(misc_start(js):misc_start(js)+SIZE(den%mmpMat(:,:,:,j))-1)=RESHAPE(REAL(den%mmpMat(:,:,:,j)),(/SIZE(den%mmpMat(:,:,:,j))/))
             vec%vec_misc(misc_start(js)+SIZE(den%mmpMat(:,:,:,j)):misc_start(js)+2*SIZE(den%mmpMat(:,:,:,j))-1)=RESHAPE(AIMAG(den%mmpMat(:,:,:,j)),(/SIZE(den%mmpMat(:,:,:,j))/))
          END IF
       END IF
    END DO

  END SUBROUTINE mixvector_from_density

  SUBROUTINE mixvector_to_density(vec,den)
    USE m_types
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_potden),    INTENT(INOUT) :: Den
     INTEGER:: js,ii,n,l,iv

     DO js=1,MERGE(jspins,3,.NOT.l_noco)
        IF (spin_here(js)) THEN
           !PW part
           IF (pw_here) THEN
              IF (sym%invs.and.js<3) THEN
                 den%pw(:,js)=vec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1)
              ELSE
                 den%pw(:,js)=CMPLX(vec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1),vec%vec_pw(pw_start(js)+stars%ng3:pw_start(js)+2*stars%ng3-1))
              ENDIF
           ENDIF
           IF (mt_here.AND.(js<3.OR.l_mtnocopot)) THEN
              !This PE stores some(or all) MT data
              ii=mt_start(js)
              DO n=mt_rank+1,atoms%ntype,mt_size
                 DO l=0,sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))
                    den%mt(:atoms%jri(n),l,n,js)=vec%vec_mt(ii:ii+atoms%jri(n)-1)
                    ii=ii+atoms%jri(n)
                 ENDDO
              ENDDO
              IF (js==3) THEN !Imaginary part comes as 4th spin
                 DO n=mt_rank+1,atoms%ntype,mt_size
                    DO l=0,sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))
                       den%mt(:atoms%jri(n),l,n,4)=vec%vec_mt(ii:ii+atoms%jri(n)-1)
                       ii=ii+atoms%jri(n)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
           IF (vac_here) THEN
              !This PE stores vac-data
              ii=vac_start(js)-1
              DO iv=1,nvac
                 den%vacz(:,iv,js)=vec%vec_vac(ii+1:ii+SIZE(den%vacz,1))
                 ii=ii+SIZE(den%vacz,1)
                 IF (sym%invs2.and.js<3)THEN
                    den%vacxy(:,:,iv,js)=RESHAPE(vec%vec_vac(ii+1:ii+SIZE(den%vacxy(:,:,iv,js))),SHAPE(den%vacxy(:,:,iv,js)))
                    ii=ii+SIZE(den%vacxy(:,:,iv,js))
                 ELSE
                    den%vacxy(:,:,iv,js)=RESHAPE(CMPLX(vec%vec_vac(ii+1:ii+SIZE(den%vacxy(:,:,iv,js))),&
                         vec%vec_vac(ii+SIZE(den%vacxy(:,:,iv,js))+1:ii+2*SIZE(den%vacxy(:,:,iv,js)))),&
                         SHAPE(den%vacxy(:,:,iv,js)))
                    ii=ii+2*SIZE(den%vacxy(:,:,iv,js))
                 ENDIF
                 IF (js>2) THEN
                    den%vacz(:,iv,4)=vec%vec_vac(ii+1:ii+SIZE(den%vacz,1))
                    ii=ii+SIZE(den%vacz,1)
              ENDIF
              ENDDO
           ENDIF
           IF (misc_here.AND.js<3) THEN
              den%mmpMat(:,:,:,js)=RESHAPE(CMPLX(vec%vec_misc(misc_start(js):misc_start(js)+SIZE(den%mmpMat(:,:,:,js))-1),vec%vec_misc(misc_start(js)+SIZE(den%mmpMat(:,:,:,js)):misc_start(js)+2*SIZE(den%mmpMat(:,:,:,js))-1)),SHAPE(den%mmpMat(:,:,:,js)))
           END IF
        END IF
     ENDDO
     call den%collect(mix_mpi_comm)
    
  END SUBROUTINE mixvector_to_density


  FUNCTION mixvector_metric(vec)RESULT(mvec)
    USE m_types
    USE m_convol
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_mixvector)                :: mvec

    INTEGER:: js,ii,n,l,iv
    COMPLEX,ALLOCATABLE::pw(:),pw_w(:)
    mvec=vec
    IF (pw_here) ALLOCATE(pw(stars%ng3),pw_w(stars%ng3))
    
    DO js=1,MERGE(jspins,3,.NOT.l_noco)
       IF (spin_here(js)) THEN
          !PW part
          IF (pw_here) THEN
             !Put back on g-grid and use convol
             IF (sym%invs.and.js<3) THEN
                pw(:)=vec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1)
             ELSE
                pw(:)=CMPLX(vec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1),vec%vec_pw(pw_start(js)+stars%ng3:pw_start(js)+2*stars%ng3-1))
             ENDIF
             CALL convol(stars,pw_w,pw,stars%ufft)
             pw_w=pw_w*cell%omtil
             mvec%vec_pw(pw_start(js):pw_start(js)+stars%ng3-1)=REAL(pw_w)
             IF ((.NOT.sym%invs).or.(js==3)) THEN 
                mvec%vec_pw(pw_start(js)+stars%ng3:pw_start(js)+2*stars%ng3-1)=AIMAG(pw_w)
             ENDIF
          ENDIF
          IF (mt_here.AND.(js<3.OR.l_mtnocopot)) THEN
             !This PE stores some(or all) MT data
             mvec%vec_mt(mt_start(js):mt_start(js)+SIZE(g_mt)-1)=g_mt*vec%vec_mt(mt_start(js):mt_start(js)+SIZE(g_mt)-1)
             IF (js==3) THEN    
                !Here we have a the imaginary part as well
                mvec%vec_mt(mt_start(js)+SIZE(g_mt):mt_stop(js))=g_mt*vec%vec_mt(mt_start(js)+SIZE(g_mt):mt_stop(js))
             ENDIF
          ENDIF
          IF (vac_here) THEN
             mvec%vec_vac(vac_start(js):vac_start(js)+SIZE(g_vac)-1)=g_vac*vec%vec_vac(vac_start(js):vac_start(js)+SIZE(g_vac)-1)
             IF (js==3) THEN !We have some extra data that corresponds to first part of metric
                mvec%vec_vac(vac_start(js)+SIZE(g_vac):vac_stop(js))=g_vac(:vac_stop(js)-vac_start(js)-SIZE(g_vac)+1)*vec%vec_vac(vac_start(js)+SIZE(g_vac):vac_stop(js))
             ENDIF
          ENDIF
          IF (misc_here.AND.(js<3)) THEN
             mvec%vec_misc(misc_start(js):misc_stop(js))=g_misc*vec%vec_misc(misc_start(js):misc_stop(js))
          END IF
       ENDIF
    END DO
   
  END FUNCTION mixvector_metric


  SUBROUTINE init_metric(oneD,vacuum)
    USE m_metrz0
    IMPLICIT NONE
    TYPE(t_oned),INTENT(in)::oneD
    TYPE(t_vacuum),INTENT(in)::vacuum

    
    INTEGER:: i,n,l,j,ivac,iz,iv2c,k2,iv2
    REAL:: dxn,dxn2,dxn4,dvol,volnstr2
    REAL,allocatable:: wght(:)
    
    IF (mt_here) THEN
       !This PE stores some(or all) MT data
       ALLOCATE(g_mt(mt_length_g)) 
       i=0
       DO n =mt_rank+1,atoms%ntype,mt_size
          dxn = atoms%neq(n) * atoms%dx(n) / 3.0
          dxn2 =2.0 * dxn
          dxn4 =4.0 * dxn
          DO l = 0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))
             i = i + 1
             g_mt(i) = dxn / atoms%rmsh(1,n)
             IF (.NOT.l_pot) THEN
                DO j = 2, atoms%jri(n) - 1, 2
                   i = i + 2
                   g_mt(i-1) = dxn4 / atoms%rmsh(j,n) 
                   g_mt(i) = dxn2 / atoms%rmsh(j+1,n) 
                END DO
                ! CHANGE JR 96/12/01
                ! take care when jri(n) is even
                i = i + 1 - MOD(atoms%jri(n),2)
                g_mt(i) = dxn / atoms%rmsh(atoms%jri(n),n)
             ELSE
                ! for the potential multiply by r^4
                DO j = 2, atoms%jri(n) - 1, 2
                   i = i + 2
                   g_mt(i-1) = dxn4 * atoms%rmsh(j,n)**3 
                   g_mt(i) = dxn2 * atoms%rmsh(j+1,n)**3
                END DO
                i = i + 1 - MOD(atoms%jri(n),2)
                g_mt(i) = dxn * atoms%rmsh(atoms%jri(n),n)**3
             END IF
          END DO
       END DO
    ENDIF
    i=0
    IF (vac_here) THEN
       iv2 = 2
       IF (sym%invs2) iv2 = 1

       ALLOCATE(g_vac(vac_length_g),wght(vacuum%nmzd))
       dvol = cell%area*vacuum%delz
       ! nvac=1 if (zrfs.or.invs)
       IF (vacuum%nvac.EQ.1) dvol = dvol + dvol
       IF (oneD%odi%d1) dvol = cell%area*vacuum%delz
       DO ivac = 1, vacuum%nvac
          ! G||=0 components
          !
          ! use 7-point simpson integration in accordance to intgz0.f
          ! calculate weights for integration
          CALL metr_z0(vacuum%nmz,wght)
          DO iz = 1, vacuum%nmz
             i = i + 1
             IF (oneD%odi%d1) THEN
                g_vac(i) = wght(iz) * dvol * (cell%z1+(iz-1)*vacuum%delz)
             ELSE
                g_vac(i) = wght(iz) * dvol
             END IF
          END DO
          ! G||.ne.0 components
          !
          ! calculate weights for integration
          CALL metr_z0(vacuum%nmzxy,wght)
          DO iv2c = 1, iv2
             DO k2 = 1, oneD%odi%nq2 - 1
                IF (oneD%odi%d1) THEN
                   DO iz = 1,vacuum%nmzxy
                      i = i + 1
                      g_vac(i) = wght(iz) * oneD%odi%nst2(k2) * dvol * (cell%z1+(iz-1)*vacuum%delz)
                   END DO
                ELSE
                   volnstr2 = dvol * stars%nstr2(k2)
                   DO iz = 1, vacuum%nmzxy
                      i = i + 1
                      g_vac(i) = wght(iz) * volnstr2
                   END DO
                END IF
             END DO
          END DO
       END DO
    END IF
    IF (misc_here) THEN
       ALLOCATE(g_misc(misc_length_g))
       g_misc=1.0
    END IF
    
    
  END SUBROUTINE init_metric
    
  
  
  SUBROUTINE init_storage_mpi(mpi_comm)
    IMPLICIT NONE
    INTEGER,INTENT(in):: mpi_comm
    INTEGER      :: irank,isize,err,js,new_comm
    mix_mpi_comm=mpi_comm
#ifdef CPP_MPI

    CALL mpi_comm_rank(mpi_comm,irank,err)
    CALL mpi_comm_size(mpi_comm,isize,err)

    IF (isize==1) RETURN !No parallelization
    js=MERGE(jspins,3,.NOT.l_noco)!distribute spins
    js=MIN(js,isize)
    CALL MPI_COMM_SPLIT(mpi_comm,MOD(irank,js),irank,new_comm,err)
    spin_here=(/MOD(irank,js)==0,MOD(irank,js)==1,(isize==2.AND.irank==0).or.MOD(irank,js)==2/)

    CALL mpi_comm_rank(new_comm,irank,err)
    CALL mpi_comm_size(new_comm,isize,err)
    CALL mpi_comm_free(new_comm,err)

    !Now distribute data   
    IF(isize==1) return !No further parallelism
    !Split off the pw-part
    pw_here=(irank==0)
    mt_here=(irank>0)
    vac_here=vac_here.AND.(irank>0)
    misc_here=misc_here.AND.(irank>0)
    isize=isize-1
    irank=irank-1
    mt_rank=irank
    mt_size=isize
    IF(isize==1.OR.irank<0) RETURN !No further parallelism
    IF (vac_here.OR.misc_here) THEN !split off-vacuum&misc part
       vac_here=vac_here.AND.(irank==0)
       misc_here=misc_here.AND.(irank==0)
       mt_here=(irank>0)
       isize=isize-1
       irank=irank-1
    ENDIF
    mt_rank=irank
    mt_size=isize
#endif
  END SUBROUTINE init_storage_mpi
      

  
  SUBROUTINE mixvector_init(mpi_comm,l_densitymatrix,oneD,input,vacuum,noco,stars_i,cell_i,sphhar_i,atoms_i,sym_i)
    USE m_types
    IMPLICIT NONE
    INTEGER,INTENT(IN)               :: mpi_comm
    LOGICAL,INTENT(IN)               :: l_densitymatrix
    TYPE(t_oneD),INTENT(IN)          :: oneD
    TYPE(t_input),INTENT(IN)         :: input
    TYPE(t_vacuum),INTENT(IN),TARGET :: vacuum
    TYPE(t_noco),INTENT(IN)          :: noco
    TYPE(t_stars),INTENT(IN),TARGET  :: stars_i
    TYPE(t_cell),INTENT(IN),TARGET   :: cell_i
    TYPE(t_sphhar),INTENT(IN),TARGET :: sphhar_i
    TYPE(t_atoms),INTENT(IN),TARGET  :: atoms_i
    TYPE(t_sym),INTENT(IN),TARGET    :: sym_i

    INTEGER::js,n,len
    

    !Store pointers to data-types
    if (associated(atoms)) return !was done before...
    jspins=input%jspins
    nvac=vacuum%nvac
    l_noco=noco%l_noco
    l_mtnocopot=noco%l_mtnocopot
    stars=>stars_i;cell=>cell_i;sphhar=>sphhar_i;atoms=>atoms_i;sym=>sym_i
    
    vac_here=input%film
    misc_here=l_densitymatrix
    CALL init_storage_mpi(mpi_comm)
    
    pw_length=0;mt_length=0;vac_length=0;misc_length=0
    mt_length_g=0;vac_length_g=0;misc_length_g=0
    DO js=1,MERGE(jspins,3,.NOT.l_noco)
       IF (spin_here(js)) THEN
          !Now calculate the length of the vectors
          IF (pw_here) THEN
             pw_start(js)=pw_length+1
             IF (sym%invs.and.js<3) THEN
                pw_length=pw_length+stars%ng3
             ELSE
                pw_length=pw_length+2*stars%ng3
             ENDIF
          ENDIF
          pw_stop(js)=pw_length
          IF (mt_here) THEN
             IF (js<3.OR.noco%l_mtnocopot) mt_start(js)=mt_length+1
             len=0
             !This PE stores some(or all) MT data
             DO n=mt_rank+1,atoms%ntype,mt_size
                len=len+(sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))+1)*atoms%jri(n)
             ENDDO
             mt_length_g=MAX(len,mt_length_g)
             IF (js==3) THEN
                !need to store imaginary part as well...
                DO n=mt_rank+1,atoms%ntype,mt_size
                   len=len+(sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n-1))+1))+1)*atoms%jri(n)
                ENDDO
             ENDIF
             IF (js<3.OR.noco%l_mtnocopot) mt_length=mt_length+len
             mt_stop(js)=mt_length
          END IF
          IF (vac_here) THEN
             !This PE stores vac-data
             vac_start(js)=vac_length+1
             len=0
             IF (sym%invs2.and.js<3) THEN
                len=len+vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + vacuum%nmzd * vacuum%nvac
             ELSE
                len=len+2*vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + vacuum%nmzd * vacuum%nvac
             ENDIF
             vac_length_g=MAX(vac_length_g,len)
             IF (js==3) len=len+vacuum%nmzd * vacuum%nvac !Offdiagnal potential is complex
             vac_length=vac_length+len
             vac_stop(js)=vac_length
          ENDIF
          IF (misc_here.AND.(js<3)) THEN
             len = 7*7*2*atoms%n_u
             misc_start(js)=misc_length+1
             misc_length = misc_length + len
             misc_stop(js)=misc_length
             misc_length_g = MAX(len,misc_length_g)
          END IF
       END IF
    END DO
    CALL init_metric(oneD,vacuum)
  END SUBROUTINE mixvector_init
  SUBROUTINE mixvector_alloc(vec)
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(OUT)    :: vec
    ALLOCATE( vec%vec_pw(pw_length) )
    ALLOCATE( vec%vec_mt(mt_length) )
    ALLOCATE( vec%vec_vac(vac_length) )
    ALLOCATE( vec%vec_misc(misc_length) )   
  END SUBROUTINE mixvector_alloc


    FUNCTION multiply_scalar(scalar,vec)RESULT(vecout)
      TYPE(t_mixvector),INTENT(IN)::vec
      REAL,INTENT(IN)             ::scalar
      TYPE(t_mixvector)           ::vecout

      vecout=vec
      vecout%vec_pw=vecout%vec_pw*scalar
      vecout%vec_mt=vecout%vec_mt*scalar
      vecout%vec_vac=vecout%vec_vac*scalar
      vecout%vec_misc=vecout%vec_misc*scalar
    END FUNCTION multiply_scalar


    FUNCTION multiply_scalar_spin(scalar,vec)RESULT(vecout)
      TYPE(t_mixvector),INTENT(IN)::vec
      REAL,INTENT(IN)             ::scalar(:)
      TYPE(t_mixvector)           ::vecout

      INTEGER:: js
      REAL:: fac
      
      vecout=vec
      DO js=1,MERGE(jspins,3,.NOT.l_noco)
         IF (SIZE(scalar)<js) THEN
            fac=0.0
         ELSE
            fac=scalar(js)
         ENDIF
         IF (pw_start(js)>0) vecout%vec_pw(pw_start(js):pw_stop(js))=vecout%vec_pw(pw_start(js):pw_stop(js))*fac
         IF (mt_start(js)>0)vecout%vec_mt(mt_start(js):mt_stop(js))=vecout%vec_mt(mt_start(js):mt_stop(js))*fac
         IF (vac_start(js)>0)vecout%vec_vac(vac_start(js):vac_stop(js))=vecout%vec_vac(vac_start(js):vac_stop(js))*fac
         IF (misc_start(js)>0)vecout%vec_misc(misc_start(js):misc_stop(js))=vecout%vec_misc(misc_start(js):misc_stop(js))*fac
      END DO
    END FUNCTION multiply_scalar_spin

    FUNCTION add_vectors(vec1,vec2)RESULT(vecout)
      TYPE(t_mixvector),INTENT(IN)::vec1,vec2
      TYPE(t_mixvector)           ::vecout
      
      vecout=vec1
      vecout%vec_pw=vecout%vec_pw+vec2%vec_pw
      vecout%vec_mt=vecout%vec_mt+vec2%vec_mt
      vecout%vec_vac=vecout%vec_vac+vec2%vec_vac
      vecout%vec_misc=vecout%vec_misc+vec2%vec_misc
    END FUNCTION add_vectors
    
    FUNCTION subtract_vectors(vec1,vec2)RESULT(vecout)
      TYPE(t_mixvector),INTENT(IN)::vec1,vec2
      TYPE(t_mixvector)           ::vecout
      
      vecout=vec1
      vecout%vec_pw=vecout%vec_pw-vec2%vec_pw
      vecout%vec_mt=vecout%vec_mt-vec2%vec_mt
      vecout%vec_vac=vecout%vec_vac-vec2%vec_vac
      vecout%vec_misc=vecout%vec_misc-vec2%vec_misc
    END FUNCTION subtract_vectors
    
    FUNCTION multiply_dot(vec1,vec2)RESULT(dprod)
      TYPE(t_mixvector),INTENT(IN)::vec1,vec2
      REAL                        ::dprod,dprod_tmp
      integer                     ::ierr
      dprod=dot_PRODUCT(vec1%vec_pw,vec2%vec_pw)
      dprod=dprod+dot_PRODUCT(vec1%vec_mt,vec2%vec_mt)
      dprod=dprod+dot_PRODUCT(vec1%vec_vac,vec2%vec_vac)
      dprod=dprod+dot_PRODUCT(vec1%vec_misc,vec2%vec_misc)
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(dprod,dprod_tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mix_mpi_comm,ierr)
      dprod=dprod_tmp
#endif      
    END FUNCTION multiply_dot

    FUNCTION multiply_dot_mask(vec1,vec2,mask,spin)RESULT(dprod)
      CLASS(t_mixvector),INTENT(IN)::vec1
      TYPE(t_mixvector),INTENT(IN)::vec2
      LOGICAL,INTENT(IN)          ::mask(4)
      INTEGER,INTENT(IN)          ::spin
      REAL                        ::dprod,dprod_tmp

      INTEGER:: js,ierr

      dprod=0.0

      DO js=1,3
         IF (mask(1).and.(spin==js.or.spin==0).and.pw_start(js)>0) &
                 dprod=dprod+dot_PRODUCT(vec1%vec_pw(pw_start(js):pw_stop(js)),&
                 vec2%vec_pw(pw_start(js):pw_stop(js)))
         IF (mask(2).and.(spin==js.or.spin==0).and.mt_start(js)>0) &
                 dprod=dprod+dot_PRODUCT(vec1%vec_mt(mt_start(js):mt_stop(js)),&
                 vec2%vec_mt(mt_start(js):mt_stop(js)))
         IF (mask(3).and.(spin==js.or.spin==0).and.vac_start(js)>0) &
                 dprod=dprod+dot_PRODUCT(vec1%vec_vac(vac_start(js):vac_stop(js)),&
                 vec2%vec_vac(vac_start(js):vac_stop(js)))
         IF (mask(4).and.(spin==js.or.spin==0).and.misc_start(js)>0) &
                 dprod=dprod+dot_PRODUCT(vec1%vec_misc(misc_start(js):misc_stop(js)),&
                 vec2%vec_misc(misc_start(js):misc_stop(js)))
      enddo

#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(dprod,dprod_tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM,mix_mpi_comm,ierr)
      dprod=dprod_tmp
#endif
    END FUNCTION multiply_dot_mask
  end MODULE m_types_mixvector
