!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_ssomat
  USE m_judft
  USE m_types_nococonv
  USE m_types_mat
  USE m_types_mpi
  USE m_types_enpara
  USE m_types_potden
  USE m_types_misc
  USE m_types_kpts
  USE m_types_lapw
  USE m_types_abc
  USE m_types_radfun
  USE m_constants
  USE m_eig66_io
  USE m_fermifct
  USE m_types_rsoc
#ifdef CPP_MPI
  USE mpi
#endif
  USE m_types_atoms
  USE m_types_cell
  USE m_types_input
  USE m_types_noco
  USE m_types_sym
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ssomat, ssomatel
CONTAINS
  SUBROUTINE ssomat(seigvso,h_so,theta,phi,eig_id,atoms,kpts,sym,&
       cell,noco,nococonv, input,fmpi,  enpara,v,results,ef )
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)         :: fmpi
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_nococonv),INTENT(IN)    :: nococonv
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_enpara),INTENT(IN)      :: enpara
    TYPE(t_potden),INTENT(IN)      :: v
    TYPE(t_results),INTENT(IN)     :: results
    INTEGER,INTENT(IN)             :: eig_id
    REAL,INTENT(in)                :: theta(:),phi(:) ! more than a single angle at once...
    REAL,INTENT(IN)                :: ef(:) !Multiple Fermi energies (bandfillings)
    REAL,INTENT(OUT)               :: seigvso(:,0:)
    REAL,INTENT(OUT)               :: h_so(0:,:,:)
    !     ..
    !     .. Locals ..
#ifdef CPP_MPI
    INTEGER:: ierr
#endif
    INTEGER :: neigf=1  !not full-matrix
    INTEGER :: jsloc,nk,n,nr,ne,nef

    COMPLEX, ALLOCATABLE :: matel(:,:,:)
    REAL,    ALLOCATABLE :: eig_shift(:,:,:,:)
    Real,    allocatable :: w_iks(:)
    COMPLEX, ALLOCATABLE :: soangl(:,:,:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: chi(:,:,:)

    TYPE(t_rsoc) :: rsoc
    TYPE(t_mat)  :: zmat
    TYPE(t_lapw) :: lapw
    TYPE(t_radfun), ALLOCATABLE :: radfun(:)
    TYPE(t_abc),    ALLOCATABLE :: abc(:,:)

    ! With one atom per type every atom is its own representative, so the
    ! coefficients need no rotation to the global frame (rot_to_unrotated).
    IF (ANY(atoms%neq/=1)) CALL judft_error('(spin spiral + soc) does not work'//&
         ' properly for more than one atom per type!',calledby="ssomat")
    IF (input%jspins/=2) CALL judft_error('ssomat requires a spin-polarized calculation',calledby="ssomat")

    seigvso=0.0
    ALLOCATE(eig_shift(input%neig,0:atoms%ntype,kpts%nkpt,SIZE(theta)));eig_shift=0.0
    ALLOCATE(matel(neigf,input%neig,0:atoms%ntype))

    ! Muffin-tin radial functions for both local spins
    ALLOCATE(radfun(atoms%ntype),abc(2,atoms%ntype))
    DO n=1,atoms%ntype
       CALL radfun(n)%generate_radial_functions(atoms,input,enpara,fmpi,v,n)
    ENDDO

    ! Spin rotation local -> global frame: chi(js,jsloc,n)
    ALLOCATE(chi(2,2,atoms%ntype))
    DO n=1,atoms%ntype
       chi(:,:,n)=nococonv%umat(n)
    ENDDO

    !Calculate the radial spin-orbit matrix elements rsoc%rso. The angular
    !matrix elements are set up below for the requested SOC directions.
    CALL rsoc%init(atoms)
    CALL rsoc%rad_matrix(atoms,noco,nococonv,input,fmpi,enpara,v)

    ALLOCATE(soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,&
         atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,SIZE(theta)))
    soangl=0.0
    DO nr=1,SIZE(theta)
       CALL rsoc%angles(atoms,fmpi,theta(nr),phi(nr))
       soangl(:,:,:,:,:,:,nr) = rsoc%soangl
    ENDDO

    DO nk=fmpi%irank+1,kpts%nkpt,fmpi%isize
       CALL lapw%init(input,noco,nococonv, kpts,atoms,sym,nk,cell)
       zMat%matsize1=lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot
       zmat%matsize2=input%neig
       zmat%l_real=.FALSE.
       IF (ALLOCATED(zmat%data_c)) DEALLOCATE(zmat%data_c)
       ALLOCATE(zmat%data_c(zMat%matsize1,zmat%matsize2))
       CALL read_eig(eig_id,nk,1,neig=ne,eig=eig_shift(:,0,nk,1),zmat=zmat)
       eig_shift(:,0,nk,1)=0.0 !not needed
       ! Coefficients in the local spin frame of each atom
       DO n=1,atoms%ntype
          DO jsloc=1,2
             CALL abc(jsloc,n)%init(input,atoms,input%neig,n)
             CALL abc(jsloc,n)%calc_abc(input,atoms,sym,cell,lapw,ne,radfun(n),noco,nococonv,jsloc,n,zMat)
          ENDDO
       ENDDO
       DO nr=1,size(theta) !loop over angles
          ! matrix elements within k
          CALL ssomatel(neigf,input,atoms,soangl(:,:,:,:,:,:,nr),rsoc%rso,&
               .TRUE.,abc,abc,chi,matel)
          eig_shift(:,0:,nk,nr)=matel(1,:,0:)
       ENDDO
    ENDDO

    !Collect data from distributed k-loop
#ifdef CPP_MPI
    IF (fmpi%irank==0) THEN
       CALL MPI_REDUCE(MPI_IN_PLACE,eig_shift,SIZE(eig_shift),MPI_DOUBLE_PRECISION,MPI_SUM,0,fmpi%mpi_comm,ierr)
    ELSE
       CALL MPI_REDUCE(eig_shift,eig_shift,SIZE(eig_shift),MPI_DOUBLE_PRECISION,MPI_SUM,0,fmpi%mpi_comm,ierr)
    ENDIF
#endif
    h_so=0.0
    IF (fmpi%irank==0) THEN
       !Sum all shift using weights
       DO nr=1,SIZE(theta)
          DO nk=1,kpts%nkpt
            DO nef=1,size(ef)
              w_iks=kpts%wtkpt(nk)*fermifct(results%eig(:,nk,1),ef(nef),input%tkb)
              !for first angle, also add unmodified eigenvalue sum
              if (nr==1) seigvso(nef,0)=seigvso(nef,nr)+dot_PRODUCT(w_iks,eig_shift(:,0,nk,nr)+results%eig(:,nk,1))
              seigvso(nef,nr)=seigvso(nef,nr)+dot_PRODUCT(w_iks,eig_shift(:,0,nk,nr)+results%eig(:,nk,1))
              DO n=0,atoms%ntype
                H_so(n,nef,nr)=H_so(n,nef,nr)+dot_PRODUCT(w_iks,eig_shift(:,n,nk,nr))
              enddo
            enddo
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE ssomat

  ! ==================================================================== !

  !> First-order SOC matrix elements in the global spin frame
  !>   matel(bandf,band2,n) = sum c1(icof;l,m1,js1)*rso(icof,jcof,n,l,jsloc1,jsloc2)
  !>                          *conjg(c2(jcof;l,m2,js2))*soangl(l,m2,js2,l,m1,js1)
  !> with the global-frame coefficients c(js) = chi(js,jsloc,n)*abc(jsloc,n)%cof.
  !> diag: only js1==js2 (both spins); else js1=2,js2=1. neigf==1: band1==band2.
  SUBROUTINE ssomatel(neigf,input,atoms,soangl,rso,diag,abc1,abc2,chi,matel)
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    LOGICAL, INTENT(IN)        :: diag
    INTEGER, INTENT(IN)        :: neigf
    ! rso(icof,jcof,ntype,0:lmaxd,2,2), slots as in t_abc%cof
    REAL,    INTENT(IN)        :: rso(:,:,:,0:,:,:)
    COMPLEX, INTENT(IN)        :: soangl(:,-atoms%lmaxd:,:,:,-atoms%lmaxd:,:)
    TYPE(t_abc), INTENT(IN)    :: abc1(:,:),abc2(:,:)  !(jsloc,ntype)
    COMPLEX, INTENT(IN)        :: chi(:,:,:)           !(js,jsloc,ntype)
    COMPLEX, INTENT(OUT)       :: matel(:,:,0:)        !(neigf,neig,0:ntype)

    INTEGER :: band1,band2,bandf,n,na,l,m1,m2,lm1,lm2,nr,icof,jcof
    INTEGER :: jsloc1,jsloc2,js1,js2,jsnumber
    COMPLEX, ALLOCATABLE :: s(:,:,:),r(:,:,:)

    ! cost per atom and spin: lmd*neig*(lmd+neig)
    matel(:,:,:)= CMPLX(0.,0.)
    nr= atoms%max_radial_functions
    ALLOCATE(s(nr,2,0:atoms%lmaxd*(atoms%lmaxd+2)),r(nr,2,0:atoms%lmaxd*(atoms%lmaxd+2)))

    IF (diag) THEN
       jsnumber= 2
    ELSE
       jsnumber= 1
    ENDIF
    DO js2= 1,jsnumber
       IF (diag) THEN
          js1= js2
       ELSE
          js1= 2
       ENDIF
       DO n= 1,atoms%ntype
          DO na= 1,atoms%neq(n)
             DO band2= 1,input%neig
                ! s = sum_m2 conjg(c2)*soangl
                DO l= 1,atoms%lmax(n)
                   DO m1= -l,l
                      lm1= l*(l+1) + m1
                      DO jsloc2= 1,2
                         DO jcof= 1,abc2(jsloc2,n)%n_r(l)
                            s(jcof,jsloc2,lm1)= CMPLX(0.,0.)
                            DO m2= -l,l
                               lm2= l*(l+1) + m2
                               s(jcof,jsloc2,lm1)= s(jcof,jsloc2,lm1) + &
                                    CONJG(chi(js2,jsloc2,n)*abc2(jsloc2,n)%cof(band2,lm2,jcof,na))&
                                    *soangl(l,m2,js2,l,m1,js1)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
                ! r = rso*s
                DO l= 1,atoms%lmax(n)
                   DO m1= -l,l
                      lm1= l*(l+1) + m1
                      DO jsloc1= 1,2
                         DO icof= 1,abc1(jsloc1,n)%n_r(l)
                            r(icof,jsloc1,lm1)= CMPLX(0.,0.)
                            DO jsloc2= 1,2
                               DO jcof= 1,abc2(jsloc2,n)%n_r(l)
                                  r(icof,jsloc1,lm1)= r(icof,jsloc1,lm1) + &
                                       rso(icof,jcof,n,l,jsloc1,jsloc2)*s(jcof,jsloc2,lm1)
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
                ! matel += c1*r
                DO bandf= 1,neigf
                   IF (neigf==input%neig) THEN
                      band1= bandf
                   ELSE
                      band1= band2
                   ENDIF
                   DO l= 1,atoms%lmax(n)
                      DO m1= -l,l
                         lm1= l*(l+1) + m1
                         DO jsloc1= 1,2
                            DO icof= 1,abc1(jsloc1,n)%n_r(l)
                               matel(bandf,band2,n)= matel(bandf,band2,n) + &
                                    chi(js1,jsloc1,n)*abc1(jsloc1,n)%cof(band1,lm1,icof,na)*r(icof,jsloc1,lm1)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO     ! band2
          ENDDO       ! na
       ENDDO         ! n
    ENDDO           ! js2,js1

    DO n= 1,atoms%ntype
       matel(:,:,0)= matel(:,:,0) + matel(:,:,n)
    ENDDO

    IF (diag) THEN
       DO n= 1,atoms%ntype
          DO band2= 1,input%neig
             IF (neigf==input%neig) THEN
                bandf= band2
             ELSE
                bandf= 1
             ENDIF
             IF (ABS(AIMAG(matel(bandf,band2,n)))>1.e-12) THEN
                PRINT *,bandf,band2,n,AIMAG(matel(bandf,band2,n))
                CALL judft_error('Stop in ssomatel:  diagonal matrix element not real')
             ENDIF
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE ssomatel
END MODULE m_ssomat
