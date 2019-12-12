!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_lapw
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  !These dimensions should be set once per call of FLEUR
  !They can be queried by the functions lapw%dim_nvd,...
  !You probably should avoid using the variables directly
  integer,save :: lapw_dim_nvd
  integer,save :: lapw_dim_nv2d
  integer,save :: lapw_dim_nbasfcn

  TYPE t_lapw
     INTEGER :: nv(2)
     INTEGER :: num_local_cols(2)
     INTEGER :: nv_tot
     INTEGER :: nmat
     INTEGER :: nlotot
     INTEGER,ALLOCATABLE:: k1(:,:)
     INTEGER,ALLOCATABLE:: k2(:,:)
     INTEGER,ALLOCATABLE:: k3(:,:)
     INTEGER,ALLOCATABLE CPP_MANAGED:: gvec(:,:,:) !replaces k1,k2,k3
     INTEGER,ALLOCATABLE:: kp(:,:)
     REAL,ALLOCATABLE::rk(:,:)
     REAL,ALLOCATABLE CPP_MANAGED::gk(:,:,:)
     REAL,ALLOCATABLE::vk(:,:,:)
     INTEGER,ALLOCATABLE::matind(:,:)
     INTEGER,ALLOCATABLE::index_lo(:,:)
     INTEGER,ALLOCATABLE::kvec(:,:,:)
     INTEGER,ALLOCATABLE::nkvec(:,:)
     REAL   :: bkpt(3)
   CONTAINS
     PROCEDURE,PASS :: init =>lapw_init
     PROCEDURE,PASS :: alloc =>lapw_alloc
     PROCEDURE,PASS :: phase_factors =>lapw_phase_factors
     PROCEDURE,NOPASS:: dim_nvd
     PROCEDURE,NOPASS:: dim_nv2d
     PROCEDURE,NOPASS:: dim_nbasfcn
     PROCEDURE,NOPASS:: init_dim=>lapw_init_dim
  END TYPE t_lapw
  PUBLIC :: t_lapw,lapw_dim_nbasfcn,lapw_dim_nvd,lapw_dim_nv2d


CONTAINS

  subroutine lapw_init_dim(nvd_in,nv2d_in,nbasfcn_in)
    IMPLICIT NONE
    INTEGER,INTENT(IN)      :: nvd_in,nv2d_in,nbasfcn_in
    lapw_dim_nvd=nvd_in
    lapw_dim_nv2d=nv2d_in
    lapw_dim_nbasfcn=nbasfcn_in
  end subroutine

  INTEGER function dim_nvd()
    dim_nvd=lapw_dim_nvd
  end function
  INTEGER function dim_nv2d()
    dim_nv2d=lapw_dim_nv2d
  end function
  INTEGER function dim_nbasfcn()
    dim_nbasfcn=lapw_dim_nbasfcn
  end function

  SUBROUTINE lapw_alloc(lapw,cell,input,noco)
    !
    !*********************************************************************
    !     determines dimensions of the lapw basis set with |k+G|<rkmax.
    !     bkpt is the k-point given in internal units
    !*********************************************************************
    USE m_boxdim
    USE m_types_fleurinput

    IMPLICIT NONE
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    CLASS(t_lapw),INTENT(INOUT)  :: lapw


    INTEGER j1,j2,j3,mk1,mk2,mk3,nv,addX,addY,addZ
    INTEGER ispin,nvh(2)

    REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3)
    ! ..
    !
    !------->          ABBREVIATIONS
    !

    !   rkmax       : cut-off for |g+k|
    !   arltv(i)    : length of reciprical lattice vector along
    !                 direction (i)
    !
    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int(input%rkmax/arltv1)+2
    mk2 = int(input%rkmax/arltv2)+2
    mk3 = int(input%rkmax/arltv3)+2

    rkm = input%rkmax
    rk2 = rkm*rkm
    !---> obtain vectors
    !---> in a spin-spiral calculation different basis sets are used for
    !---> the two spin directions, because the cutoff radius is defined
    !---> by |G + k +/- qss/2| < rkmax.
    nvh(2)=0
    DO ispin = 1,MERGE(2,1,noco%l_ss)
       addX = abs(NINT((lapw%bkpt(1)+ (2*ispin - 3)/2.0*noco%qss(1))/arltv1))
       addY = abs(NINT((lapw%bkpt(2)+ (2*ispin - 3)/2.0*noco%qss(2))/arltv2))
       addZ = abs(NINT((lapw%bkpt(3)+ (2*ispin - 3)/2.0*noco%qss(3))/arltv3))
       nv = 0
       DO  j1 = -mk1-addX,mk1+addX
          DO  j2 = -mk2-addY,mk2+addY
             DO  j3 = -mk3-addZ,mk3+addZ
                s = lapw%bkpt + (/j1,j2,j3/) + (2*ispin - 3)/2.0*noco%qss
                r2 = dot_PRODUCT(MATMUL(s,cell%bbmat),s)
                IF (r2.LE.rk2)  nv = nv + 1
             END DO
          END DO
       END DO
       nvh(ispin)  = nv
    END DO
    nv  = MAX(nvh(1),nvh(2))

    IF (ALLOCATED(lapw%rk)) THEN
       IF (SIZE(lapw%rk)==nv) THEN
          RETURN !
       ELSE
          DEALLOCATE(lapw%rk,lapw%gvec,lapw%vk,lapw%gk,lapw%matind)
          DEALLOCATE(lapw%k1,lapw%k2,lapw%k3)
       ENDIF
    ENDIF
    ALLOCATE(lapw%rk(nv,input%jspins) )
    ALLOCATE(lapw%gvec(3,nv,input%jspins))
    ALLOCATE(lapw%vk(3,nv,input%jspins))
    ALLOCATE(lapw%gk(3,nv,input%jspins))
    ALLOCATE(lapw%k1(nv,input%jspins)) !shpuld be removed
    ALLOCATE(lapw%k2(nv,input%jspins)) !
    ALLOCATE(lapw%k3(nv,input%jspins)) !
    ALLOCATE(lapw%matind(nv,2))

    lapw%rk = 0 ; lapw%gvec = 0 ;lapw%nv=0
  END SUBROUTINE lapw_alloc


  SUBROUTINE lapw_init(lapw,input,noco,kpts,atoms,sym,&
       nk,cell,l_zref,mpi)
    USE m_types_mpi
    USE m_sort
    USE m_boxdim
    USE m_types_fleurinput
    USE m_types_kpts
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_mpi),INTENT(IN),OPTIONAL:: mpi
    CLASS(t_lapw),INTENT(INOUT)    :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT  (IN) :: nk
    LOGICAL, INTENT (IN)  :: l_zref
    !     ..
    !     .. Array Arguments ..
    !     ..
    !     .. Local Scalars ..
    REAL arltv1,arltv2,arltv3,r2,rk2,rkm,r2q,gla,eps,t
    INTEGER i,j,j1,j2,j3,k,l ,mk1,mk2,mk3,n,ispin,gmi,m,nred,n_inner,n_bound,itt(3),addX,addY,addZ
    !     ..
    !     .. Local Arrays ..
    REAL                :: s(3),sq(3)
    REAL,ALLOCATABLE    :: rk(:),rkq(:),rkqq(:)
    INTEGER,ALLOCATABLE :: gvec(:,:),index3(:)


    !     ..
    !---> in a spin-spiral calculation different basis sets are used for
    !---> the two spin directions, because the cutoff radius is defined
    !---> by |G + k +/- qss/2| < rkmax.

    IF (nk>kpts%nkpt) THEN
       lapw%bkpt(:)=kpts%bkf(:,nk)
    ELSE
       lapw%bkpt(:) = kpts%bk(:,nk)
    ENDIF


    CALL lapw%alloc(cell,input,noco)

    ALLOCATE(gvec(3,SIZE(lapw%gvec,2)))
    ALLOCATE(rk(SIZE(lapw%gvec,2)),rkq(SIZE(lapw%gvec,2)),rkqq(SIZE(lapw%gvec,2)))
    ALLOCATE(index3(SIZE(lapw%gvec,2)))


    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !     arltv(i) length of reciprical lattice vector along direction (i)
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int( input%rkmax/arltv1 )+4
    mk2 = int( input%rkmax/arltv2 )+4
    mk3 = int( input%rkmax/arltv3 )+4

    rk2 = input%rkmax*input%rkmax
    !---> if too many basis functions, reduce rkmax
    spinloop:DO ispin = 1,input%jspins
       addX = abs(NINT((lapw%bkpt(1)+ (2*ispin - 3)/2.0*noco%qss(1))/arltv1))
       addY = abs(NINT((lapw%bkpt(2)+ (2*ispin - 3)/2.0*noco%qss(2))/arltv2))
       addZ = abs(NINT((lapw%bkpt(3)+ (2*ispin - 3)/2.0*noco%qss(3))/arltv3))
       !--->    obtain vectors
       n = 0
       DO  j1 = -mk1-addX,mk1+addX
          DO  j2 = -mk2-addY,mk2+addY
             DO  j3 = -mk3-addZ,mk3+addZ
                s=lapw%bkpt+(/j1,j2,j3/)+(2*ispin - 3)/2.0*noco%qss
                sq = lapw%bkpt+ (/j1,j2,j3/)
                r2 = dot_PRODUCT(s,MATMUL(s,cell%bbmat))
                r2q = dot_PRODUCT(sq,MATMUL(sq,cell%bbmat))
                IF (r2.LE.rk2) THEN
                   n = n + 1
                   gvec(:,n) = (/j1,j2,j3/)
                   rk(n) = SQRT(r2)
                   rkq(n) = SQRT(r2q)
                END IF
             ENDDO
          ENDDO
       ENDDO
       lapw%nv(ispin) = n

       !Sort according to k+g, first construct secondary sort key
       DO  k = 1,lapw%nv(ispin)
          rkqq(k) = (mk1+gvec(1,k)) + (mk2+gvec(2,k))*(2*mk1+1) +&
               (mk3+gvec(3,k))*(2*mk1+1)*(2*mk2+1)
       ENDDO
       CALL sort(index3(:lapw%nv(ispin)),rkq,rkqq)
       DO n=1,lapw%nv(ispin)
          lapw%gvec(:,n,ispin) = gvec(:,index3(n))
          lapw%rk(n,ispin)     = rk(index3(n))
       ENDDO
       !--->    determine pairs of K-vectors, where K_z = K'_-z to use
       !--->    z-reflection
       IF (l_zref) THEN
          n=0
          DO i=1,lapw%nv(ispin)
             DO j=1,i
                IF (ALL(lapw%gvec(1:2,i,ispin).EQ.lapw%gvec(1:2,j,ispin)).AND.&
                     (lapw%gvec(3,i,ispin).EQ.-lapw%gvec(3,j,ispin))) THEN
                   n=n+1
                   lapw%matind(n,1)=i
                   lapw%matind(n,2)=j
                ENDIF
             ENDDO
          ENDDO
          nred=n
          IF (PRESENT(mpi)) THEN
             IF (mpi%n_size.GT.1) THEN
                !
                !--->     order K's in sequence K_1,...K_n | K_0,... | K_-1,....K_-n
                !
                n_inner = lapw%nv(ispin) - nred
                IF (MOD(nred,mpi%n_size).EQ.0) THEN
                   n_bound = nred
                ELSE
                   n_bound = (1+INT( nred/mpi%n_size ))*mpi%n_size
                ENDIF
                IF (lapw%nv(ispin) - nred + n_bound.GT.SIZE(lapw%gvec,2)) THEN
                   CALL juDFT_error("BUG:z-ref & ev || : dimension too small!" ,calledby ="types_lapw")
                ENDIF

                i = 1
                j = 1
                DO n = 1, nred
                   IF (lapw%matind(n,1).EQ.lapw%matind(n,2)) THEN
                      index3(lapw%matind(n,1)) = n_inner + i
                      i = i + 1
                   ELSE
                      index3(lapw%matind(n,1)) = j
                      index3(lapw%matind(n,2)) = j + n_bound
                      j = j + 1
                   ENDIF
                ENDDO
                !--->          resort the rk,k1,k2,k3 and lapw%matind arrays:
                DO n = 1, lapw%nv(ispin)
                   rk(n)  = lapw%rk(n,ispin)
                   gvec(:,n) = lapw%gvec(:,n,ispin)
                ENDDO
                DO n = lapw%nv(ispin), 1, -1
                   lapw%rk(index3(n),ispin) = rk(n)
                   lapw%gvec(:,index3(n),ispin) = gvec(:,n)
                ENDDO
                DO n = nred + 1, n_bound
                   lapw%rk(n,ispin) = lapw%rk(lapw%nv(ispin),ispin)
                   lapw%gvec(:,n,ispin) = lapw%gvec(:,lapw%nv(ispin),ispin)
                ENDDO
                lapw%nv(ispin) = lapw%nv(ispin) - nred + n_bound
             ENDIF
          ENDIF
       ENDIF
       DO k=1,lapw%nv(ispin)
          lapw%vk(:,k,ispin)=lapw%bkpt+lapw%gvec(:,k,ispin)+(ispin-1.5)*noco%qss
          lapw%gk(:,k,ispin)=MATMUL(TRANSPOSE(cell%bmat),lapw%vk(:,k,ispin))/MAX (lapw%rk(k,ispin),1.0e-30)
       ENDDO


       IF (.NOT.noco%l_ss.AND.input%jspins==2) THEN
          !Second spin is the same
          lapw%nv(2)=lapw%nv(1)
          lapw%gvec(:,:,2)=lapw%gvec(:,:,1)
          lapw%rk(:,2)=lapw%rk(:,1)
          lapw%vk(:,:,2)=lapw%vk(:,:,1)
          lapw%gk(:,:,2)=lapw%gk(:,:,1)
          EXIT spinloop
       END IF

    ENDDO spinloop
    !should be removed later...
    lapw%k1=lapw%gvec(1,:,:)
    lapw%k2=lapw%gvec(2,:,:)
    lapw%k3=lapw%gvec(3,:,:)

    !Count No of lapw distributed to this PE
    lapw%num_local_cols=0
    DO ispin=1,input%jspins
       IF (PRESENT(mpi)) THEN
          DO k=mpi%n_rank+1,lapw%nv(ispin),mpi%n_size
             lapw%num_local_cols(ispin)=lapw%num_local_cols(ispin)+1
          END DO
       ELSE
          lapw%num_local_cols(ispin) = lapw%nv(ispin)
       END IF
    END DO

    IF (ANY(atoms%nlo>0)) CALL priv_lo_basis_setup(lapw,atoms,sym,noco,cell)

    lapw%nv_tot=lapw%nv(1)
    lapw%nmat=lapw%nv(1)+atoms%nlotot
    IF (noco%l_noco) lapw%nv_tot=lapw%nv_tot+lapw%nv(2)
    IF (noco%l_noco) lapw%nmat=lapw%nv_tot+2*atoms%nlotot

  CONTAINS

    SUBROUTINE priv_lo_basis_setup(lapw,atoms,sym,noco,cell)
      USE m_types_fleurinput

      IMPLICIT NONE
      TYPE(t_lapw),INTENT(INOUT):: lapw
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_noco),INTENT(IN)   :: noco


      INTEGER:: n,na,nn,np,lo,nkvec_sv,nkvec(atoms%nlod,2),iindex
      IF (.NOT.ALLOCATED(lapw%kvec)) THEN
         ALLOCATE(lapw%kvec(2*(2*atoms%llod+1),atoms%nlod,atoms%nat))
         ALLOCATE(lapw%nkvec(atoms%nlod,atoms%nat))
         ALLOCATE(lapw%index_lo(atoms%nlod,atoms%nat))
      ENDIF
      iindex=0
      na=0
      nkvec_sv=0
      DO n=1,atoms%ntype
         DO nn=1,atoms%neq(n)
            na=na+1
            if (sym%invsat(na)>1) cycle
            !np = MERGE(oneD%ods%ngopr(na),sym%invtab(sym%ngopr(na)),oneD%odi%d1)
            np=sym%invtab(sym%ngopr(na))
            CALL priv_vec_for_lo(atoms,sym,na,n,np,noco,lapw,cell)
            DO lo = 1,atoms%nlo(n)
               lapw%index_lo(lo,na)=iindex
               iindex=iindex+lapw%nkvec(lo,na)
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE priv_lo_basis_setup

  END SUBROUTINE lapw_init


  SUBROUTINE lapw_phase_factors(lapw,iintsp,tau,qss,cph)
    USE m_constants
    USE m_types_fleurinput
    IMPLICIT NONE
    CLASS(t_lapw),INTENT(in):: lapw
    INTEGER,INTENT(IN)     :: iintsp
    REAL,INTENT(in)        :: tau(3),qss(3)
    COMPLEX,INTENT(out)    :: cph(:)

    INTEGER:: k
    REAL:: th
    DO k = 1,lapw%nv(iintsp)
       th= DOT_PRODUCT(lapw%gvec(:,k,iintsp)+(iintsp-1.5)*qss,tau)
       cph(k) = CMPLX(COS(tpi_const*th),-SIN(tpi_const*th))
    END DO
  END SUBROUTINE lapw_phase_factors


  SUBROUTINE priv_vec_for_lo(atoms,sym,na,&
       n,np,noco, lapw,cell)
    USE m_constants,ONLY: tpi_const,fpi_const
    USE m_orthoglo
    USE m_ylm
    USE m_types_fleurinput
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(INOUT):: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,n,np
    !     ..
    !     .. Array Arguments ..
    !     ..
    !     .. Local Scalars ..
    COMPLEX term1
    REAL th,con1
    INTEGER l,lo ,mind,ll1,lm,iintsp,k,nkmin,ntyp,lmp,m,nintsp
    LOGICAL linind,enough,l_lo1,l_real
    !     ..
    !     .. Local Arrays ..
    INTEGER :: nkvec(atoms%nlod,2)
    REAL qssbti(3),bmrot(3,3),v(3),vmult(3)
    REAL :: gkrot(3,SIZE(lapw%gk,2),2)
    REAL :: rph(SIZE(lapw%gk,2),2)
    REAL :: cph(SIZE(lapw%gk,2),2)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    COMPLEX cwork(-2*atoms%llod:2*atoms%llod+1,2*(2*atoms%llod+1),atoms%nlod ,2)
    !     ..
    !     .. Data statements ..
    REAL, PARAMETER :: eps = 1.0E-8
    REAL, PARAMETER :: linindq = 1.0e-4

    l_real = sym%invs.and..not.noco%l_noco.and..not.(noco%l_soc.and.atoms%n_hia+atoms%n_u>0)

    con1=fpi_const/SQRT(cell%omtil)
    ntyp = n
    nintsp=MERGE(2,1,noco%l_ss)
    DO iintsp = 1,nintsp
       IF (iintsp.EQ.1) THEN
          qssbti = - noco%qss/2
       ELSE
          qssbti = + noco%qss/2
       ENDIF

       !--->    set up phase factors
       DO k = 1,lapw%nv(iintsp)
          th= tpi_const*DOT_PRODUCT((/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+qssbti,atoms%taual(:,na))
          rph(k,iintsp) = COS(th)
          cph(k,iintsp) = -SIN(th)
       END DO

       IF (np.EQ.1) THEN
          gkrot(:,:,iintsp)=lapw%gk(:,:,iintsp)
       ELSE
          bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
          DO k = 1,lapw%nv(iintsp)
             !-->           apply the rotation that brings this atom into the
             !-->           representative (this is the definition of ngopr(na))
             !-->           and transform to cartesian coordinates
             v(:) = lapw%vk(:,k,iintsp)
             gkrot(:,k,iintsp) = MATMUL(v,bmrot)
          END DO
       END IF
       !--->   end loop over interstitial spin
    ENDDO

    nkvec(:,:) = 0
    cwork(:,:,:,:) = CMPLX(0.0,0.0)
    enough=.FALSE.
    DO k = 1,MIN(lapw%nv(1),lapw%nv(nintsp))
       IF (ANY(lapw%rk(k,:nintsp).LT.eps)) CYCLE
       IF (.NOT.enough) THEN
          DO iintsp = 1,nintsp

             !-->        generate spherical harmonics
             vmult(:) =  gkrot(:,k,iintsp)
             CALL ylm4(atoms%lnonsph(ntyp),vmult, ylm)
                enough = .TRUE.
                term1 = con1* ((atoms%rmt(ntyp)**2)/2)* CMPLX(rph(k,iintsp),cph(k,iintsp))
                DO lo = 1,atoms%nlo(ntyp)
                   IF (sym%invsat(na).EQ.0) THEN
                      IF ((nkvec(lo,iintsp)).LT. (2*atoms%llo(lo,ntyp)+1)) THEN
                         enough = .FALSE.
                         nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                         l = atoms%llo(lo,ntyp)
                         ll1 = l*(l+1) + 1
                         DO m = -l,l
                            lm = ll1 + m
                            cwork(m,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                         END DO
                         CALL orthoglo(l_real,atoms,nkvec(lo,iintsp),lo,l,linindq,.FALSE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                         IF (linind) THEN
                            lapw%kvec(nkvec(lo,iintsp),lo,na) = k
                         ELSE
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                         ENDIF
                      ENDIF
                   ELSE
                      IF ((sym%invsat(na).EQ.1) .OR. (sym%invsat(na).EQ.2)) THEN
                         IF (nkvec(lo,iintsp).LT.2*(2*atoms%llo(lo,ntyp)+1)) THEN
                            enough = .FALSE.
                            nkvec(lo,iintsp) = nkvec(lo,iintsp) + 1
                            l = atoms%llo(lo,ntyp)
                            ll1 = l*(l+1) + 1
                            DO m = -l,l
                               lm = ll1 + m
                               mind = -l + m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = term1*ylm(lm)
                               mind = l + 1 + m
                               lmp = ll1 - m
                               cwork(mind,nkvec(lo,iintsp),lo,iintsp) = ((-1)** (l+m))*CONJG(term1*ylm(lmp))
                            END DO
                            CALL orthoglo(l_real,atoms,nkvec(lo,iintsp),lo,l,linindq,.TRUE., cwork(-2*atoms%llod,1,1,iintsp),linind)
                            IF (linind) THEN
                               lapw%kvec(nkvec(lo,iintsp),lo,na) = k
                               !                          write(*,*) nkvec(lo,iintsp),k,' <- '
                            ELSE
                               nkvec(lo,iintsp) = nkvec(lo,iintsp) - 1
                            END IF
                         END IF
                      END IF
                   END IF
                END DO
                IF ((k.EQ.lapw%nv(iintsp)) .AND. (.NOT.enough)) THEN
                   WRITE (6,FMT=*) 'vec_for_lo did not find enough linearly independent'
                   WRITE (6,FMT=*) 'clo coefficient-vectors. the linear independence'
                   WRITE (6,FMT=*) 'quality, linindq, is set: ',linindq
                   WRITE (6,FMT=*) 'this value might be to large.'
                   WRITE(*,*) na,k,lapw%nv
                   CALL juDFT_error("not enough lin. indep. clo-vectors" ,calledby ="vec_for_lo")
                END IF
             ! -- >        end of abccoflo-part
          ENDDO
       ENDIF

       ! -->    check whether we have already enough k-vecs
       enough=.TRUE.
       DO lo = 1,atoms%nlo(ntyp)
          IF (nkvec(lo,1).EQ.nkvec(lo,nintsp)) THEN   ! k-vec accepted by both spin channels
             IF (sym%invsat(na).EQ.0) THEN
                IF ( nkvec(lo,1).LT.(2*atoms%llo(lo,ntyp)+1) ) THEN
                   enough=.FALSE.
                ENDIF
             ELSE
                IF ( nkvec(lo,1).LT.(2*(2*atoms%llo(lo,ntyp)+1))  ) THEN
                   enough=.FALSE.
                ENDIF
             ENDIF
          ELSE
             nkmin = MIN(nkvec(lo,1),nkvec(lo,nintsp)) ! try another k-vec
             nkvec(lo,1) = nkmin ; nkvec(lo,nintsp) = nkmin
             enough=.FALSE.
          ENDIF
       ENDDO
       IF ( enough ) THEN
          lapw%nkvec(:atoms%nlo(ntyp),na)=nkvec(:atoms%nlo(ntyp),1)
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE priv_vec_for_lo
END MODULE m_types_lapw
