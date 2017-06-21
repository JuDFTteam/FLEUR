!
!     Calculates the Coulomb matrix
!
!     v      =  < M    | v | M    >
!      k,IJ        k,I        k,J
!
!     with the mixed-basis functions M (indices I and J).
!
!     Note that
!                 *
!     v      =  v     .
!      k,JI      k,IJ
!
!     In the code: coulomb(IJ,k) = v     where only the upper triangle (I<=J) is stored.
!                                   k,IJ 
!
!     The Coulomb matrix v(IJ,k) diverges at the Gamma-point. Here, we apply the decomposition
!
!              (0)        (1)   *        2-l              (0)*   (0)    (1)*        m  (1)
!     v     = v    + SUM v   * Y  (k) / k        with    v    = v   ,  v      = (-1)  v
!      k,IJ    IJ     lm  IJ    lm                        JI     IJ     JI,lm          IJ,l,-m
!
!     where a = atom index, R  = position vector, T  = Wigner-Seitz radius (scalar).
!                            a                     0
!                                    (0)
!     In the code: coulomb(IJ,1)  = v    where only the upper triangle (I<=J) is stored,
!                                    IJ
!                                    (1)
!                  coulfac(IJ,lm) = v
!                                    IJ,lm
!
!     For the PW contribution we have to construct plane waves within the MT spheres with the help
!     of spherical Bessel functions. The value lexp (LEXP in gwinp) is the corresponding cutoff.
!
      MODULE m_coulomb

      CONTAINS


      SUBROUTINE coulombmatrix(mpi,atoms,kpts,&
           cell, sym, hybrid, xcpot,l_restart)

      USE m_constants    , ONLY : pi_const
      USE m_olap         , ONLY : olap_pw,gptnorm
      USE m_trafo        , ONLY : symmetrize,bramat_trafo
      USE m_util         , ONLY : sphbessel,intgrf,intgrf_init,&
     &                            harmonicsr,primitivef
      USE m_hsefunctional, ONLY : change_coulombmatrix
#ifdef CPP_MPI
      USE m_mpi_work_dist, ONLY : work_dist_coulomb
#endif
      USE m_wrapper
      USE m_icorrkeys

    USE m_types

      IMPLICIT NONE

      TYPE(t_xcpot),INTENT(IN)     :: xcpot
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_hybrid),INTENT(INOUT) :: hybrid
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_kpts),INTENT(IN)      :: kpts
      TYPE(t_atoms),INTENT(IN)     :: atoms

      ! - scalars -
      LOGICAL    , INTENT(IN)    :: l_restart
     

      ! - local scalars -
      INTEGER                    :: maxfac
!       INTEGER                   ::  irecl_coulomb
      INTEGER                    :: inviop
      INTEGER                    :: nqnrm,iqnrm,iqnrm1,iqnrm2,&
     &                              iqnrmstart,iqnrmstep
      INTEGER                    :: itype,l ,ix,iy,iy0,i,j,lm,l1,l2,m1,&
     &                              m2,ineq,idum,ikpt,ikpt0,ikpt1
      INTEGER                    :: lm1,lm2,itype1,itype2,ineq1,ineq2,n,&
     &                              n1,n2,ng
      INTEGER                    :: ic,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8
      INTEGER                    :: igpt,igpt1,igpt2,igptp,igptp1,igptp2
      INTEGER                    :: isym,isym1,isym2,igpt0
      INTEGER                    :: maxlcut,ok
      INTEGER                    :: nbasp,m
      INTEGER                    :: ikptmin,ikptmax,nkminmax
      INTEGER                    :: idisp,icnt
#ifdef CPP_INVERSION
      INTEGER    , PARAMETER     :: coul_size =  8
#else
      INTEGER    , PARAMETER     :: coul_size = 16
#endif

      LOGICAL                    :: lsym

      REAL                       :: rdum,rdum1,rdum2
      REAL                       :: svol,qnorm,qnorm1,qnorm2,gnorm
      REAL                       :: fcoulfac
      REAL                       :: time1,time2

      COMPLEX                    :: cdum,cdum1,cexp,csum
      COMPLEX    , PARAMETER     :: img = (0d0,1d0)

      ! - local arrays -
      INTEGER                    :: g(3)
      INTEGER                    :: nbasm1(kpts%nkptf)
      INTEGER    , ALLOCATABLE   :: pqnrm(:,:)
      INTEGER                    :: rrot(3,3,sym%nsym),invrrot(3,3,sym%nsym)
      INTEGER    , ALLOCATABLE   :: iarr(:),pointer(:,:,:,:)!,pointer(:,:,:)
      INTEGER                    :: igptmin(kpts%nkpt),igptmax(kpts%nkpt)
      INTEGER    , ALLOCATABLE   :: nsym_gpt(:,:),sym_sym1(:,:,:), sym_gpt(:,:,:)
      INTEGER                    :: nsym1(kpts%nkpt+1),ngptm1(kpts%nkpt+1), sym1(sym%nsym,kpts%nkpt+1)
   
      INTEGER                    :: comm(kpts%nkpt)

      LOGICAL                    :: calc_mt(kpts%nkpt)

      REAL                       :: q(3),q1(3),q2(3)
      REAL                       :: integrand(atoms%jmtd),primf1(atoms%jmtd),&
     &                              primf2(atoms%jmtd)
      REAL                       :: mat(hybrid%maxindxm1*(hybrid%maxindxm1+1)/2)
      REAL                       :: moment(hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype),&
     &                              moment2(hybrid%maxindxm1,atoms%ntype)
      REAL                       :: sphbes(atoms%jmtd,0:hybrid%maxlcutm1)
      REAL                       :: sphbesmoment1(atoms%jmtd,0:hybrid%maxlcutm1)
      REAL                       :: rarr(0:hybrid%lexp+1),rarr1(0:hybrid%maxlcutm1)
      REAL       , ALLOCATABLE   :: fac(:),sfac(:),facfac(:)
      REAL       , ALLOCATABLE   :: gmat(:,:),qnrm(:),eig(:)
      REAL       , ALLOCATABLE   :: sphbesmoment(:,:,:)
      REAL       , ALLOCATABLE   :: sphbes0(:,:,:)   
      REAL       , ALLOCATABLE   :: olap(:,:,:,:),integral(:,:,:,:)
      REAL       , ALLOCATABLE   :: gridf(:,:)

      COMPLEX                    :: cexp1(atoms%ntype),csumf(9)
      COMPLEX                    :: structconst((2*hybrid%lexp+1)**2 ,atoms%nat,atoms%nat,&
     &                                          kpts%nkpt)             ! nw = 1
      COMPLEX                    :: y((hybrid%lexp+1)**2),y1((hybrid%lexp+1)**2),&
     &                              y2((hybrid%lexp+1)**2)
      COMPLEX                    :: imgl(0:hybrid%lexp)
      COMPLEX                    :: dwgn(-hybrid%maxlcutm1:hybrid%maxlcutm1,&
     &                                   -hybrid%maxlcutm1:hybrid%maxlcutm1,&
     &                                           0:hybrid%maxlcutm1,sym%nsym)
      COMPLEX    , ALLOCATABLE   :: smat(:,:)
      COMPLEX    , ALLOCATABLE   :: coulmat(:,:)
      COMPLEX    , ALLOCATABLE   :: carr1(:),carr2(:,:),carr2a(:,:),&
     &                              carr2b(:,:)
      COMPLEX    , ALLOCATABLE   :: structconst1(:,:)

      REAL       , ALLOCATABLE   :: eval(:)
      REAL       , ALLOCATABLE   :: coulomb_mt1(:,:,:,:,:)
#if ( defined(CPP_INVERSION) )
      REAL       , ALLOCATABLE   :: coulomb(:,:),coulhlp(:,:),&
     &                              coulhlp1(:,:)
      REAL       , ALLOCATABLE   :: involapm(:,:),olapm(:,:),evec(:,:),&
     &                              invevec(:,:)
      REAL       , ALLOCATABLE   :: coulomb_mt2(:,:,:,:,:),&
     &                              coulomb_mt3(:,:,:,:)
      REAL       , ALLOCATABLE   :: coulomb_mtir(:,:,:),&
     &                              coulombp_mtir(:,:)
#else
      COMPLEX   , ALLOCATABLE   :: coulomb(:,:),coulhlp(:,:),&
     &                              coulhlp1(:,:)
      COMPLEX   , ALLOCATABLE   :: olapm(:,:),evec(:,:),invecec(:,:)
      COMPLEX   , ALLOCATABLE   :: coulomb_mt2(:,:,:,:,:),&
     &                              coulomb_mt3(:,:,:,:)
      COMPLEX   , ALLOCATABLE   :: coulomb_mtir(:,:,:),&
     &                              coulombp_mtir(:,:)
#endif
      INTEGER                    :: ishift,ishift1
      INTEGER                    :: iatom,iatom1
      INTEGER                    :: indx1,indx2,indx3,indx4
      LOGICAL                    :: l_found,l_warn,l_warned,&
     &                              l_plot = .false.!.true.!.false.

#ifdef CPP_NOSPMVEC
      CHARACTER*8 , PARAMETER    :: fname = 'coulomb'
#else
      CHARACTER*8 , PARAMETER    :: fname = 'coulomb1'
#endif

#ifndef CPP_MPI
      INTEGER                    :: irecl_coulomb
#else
      INCLUDE 'mpif.h'
      INTEGER :: fh,ierr,reqd,length,irank2,isize2
      INTEGER, ALLOCATABLE :: req(:)
      INTEGER(KIND=MPI_OFFSET_KIND)  :: filesize,offset0,offset,&
     &                                  irecl_coulomb
      INTEGER(KIND=MPI_ADDRESS_KIND) :: iaddr
#ifdef CPP_INVERSION
      INTEGER, PARAMETER :: datatype = MPI_REAL8
#else
      INTEGER, PARAMETER :: datatype = MPI_COMPLEX16
#endif
      CHARACTER(MPI_MAX_ERROR_STRING) :: errmsg
#endif

      CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,gridf)

      !hybrid%maxgptm = maxval(hybrid%ngptm(:))
      nbasm1   = 0 
      nbasp   = 0
      DO itype=1,atoms%ntype
        DO i=1,atoms%neq(itype)
          DO l=0,hybrid%lcutm1(itype)
            DO M=-l,l
              DO j=1,hybrid%nindxm1(l,itype)
                nbasp = nbasp + 1
              END DO
            END DO
          END DO
        END DO
      END DO
      hybrid%maxbasm1  = nbasp  + hybrid%maxgptm
      nbasm1    = nbasp  + hybrid%ngptm(:)

      svol     = sqrt(cell%vol)
      fcoulfac = 4*pi_const/cell%vol
      maxlcut  = maxval( atoms%lmax ) 
      maxfac   = max(2*maxlcut+hybrid%maxlcutm1+1,4*max(hybrid%maxlcutm1,hybrid%lexp)+1)

      ALLOCATE ( fac( 0:maxfac),sfac( 0:maxfac),facfac(-1:maxfac) )
      fac(0)       = 1                    !
      sfac(0)      = 1                    ! Define:
      facfac(-1:0) = 1                    ! fac(i)    = i!
      DO i=1,maxfac                       ! sfac(i)   = sqrt(i!)
        fac(i)    = fac(i-1)*i            ! facfac(i) = (2i+1)!!
        sfac(i)   = sfac(i-1)*sqrt(i*1d0) !
        facfac(i) = facfac(i-1)*(2*i+1)   !
      END DO

      ! Define imgl(l) = img**l
      imgl(0) = 1
      DO i=1,hybrid%lexp
        imgl(i) = imgl(i-1) * img
      END DO

#     ifdef CPP_NOSPMVEC
        irecl_coulomb = hybrid%maxbasm1 * (hybrid%maxbasm1+1) * coul_size / 2

      ! if the sparse matrix technique is used, several entries of the
      ! matrix vanish so that the size of each entry is smaller
#     else
        irecl_coulomb = ( atoms%ntype*(hybrid%maxlcutm1+1)*(hybrid%maxindxm1-1)**2&
     &                +   atoms%nat *(hybrid%maxlcutm1+2)*(2*hybrid%maxlcutm1+1)*(hybrid%maxindxm1-1)&
     &                +   (hybrid%maxindxm1-1)*atoms%nat**2&
     &                +   ((hybrid%maxlcutm1+1)**2*atoms%nat+hybrid%maxgptm)&
     &                   *((hybrid%maxlcutm1+1)**2*atoms%nat+hybrid%maxgptm+1)/2 )*coul_size
#     endif



!     Calculate the structure constant
      CALL structureconstant(structconst,cell,hybrid,&
     &                       atoms,kpts,&
     &                       mpi)


      ! Don't calculate new Coulomb matrix if old one is found
      IF ( l_restart ) THEN
        INQUIRE( FILE=fname, EXIST=l_found )
        IF ( l_found ) RETURN
      END IF

      ! Open file for output
#     ifdef CPP_MPI
        CALL MPI_FILE_OPEN(mpi,fname,&
     &            MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
!        ! reserve memory for the file
!        filesize = irecl_coulomb * kpts%nkpt(nw)
!        CALL MPI_FILE_PREALLOCATE(fh,filesize,ierr)
!        IF ( ierr /= 0 ) THEN
!          CALL MPI_ERROR_STRING(ierr,errmsg,length,ierr)
!          WRITE(*,*) errmsg(1:length)
!          STOP 'coulombmatrix failed preallocating coulomb matrix'
!        END IF
#     endif


      IF ( mpi%irank == 0 )&
     &  WRITE(6,'(//A)') '### subroutine: coulombmatrix ###'

!
!     Matrix allocation
!

      IF(allocated(coulomb)) DEALLOCATE (coulomb)

      ALLOCATE ( coulomb(hybrid%maxbasm1*(hybrid%maxbasm1+1)/2,kpts%nkpt) , stat = ok )
      IF( ok .ne. 0 )&
     &     STOP 'coulombmatrix: failure allocation coulomb matrix'
      coulomb = 0

      IF ( mpi%irank == 0 ) WRITE(6,'(/A,F6.1," MB")')&
     &              'Size of coulomb matrix:',16d0/1048576*size(coulomb)

!
!     Generate Symmetry:
!     Reduce list of g-Points so that only one of each symm-equivalent is calculated
!
#     ifndef CPP_NOCOULSYM

        IF ( mpi%irank == 0 )&
     &    WRITE(6,'(/A)',advance='no') 'Setup for symmetry...'
        CALL cpu_time(time1)
        ! calculate rotations in reciprocal space
        DO isym = 1,sym%nsym
          IF( isym .le. sym%nop ) THEN
            inviop         = sym%invtab(isym)
            rrot(:,:,isym) = transpose(sym%mrot(:,:,inviop))
            DO l = 0,hybrid%maxlcutm1
              dwgn(:,:,l,isym) = transpose( &
     &             sym%d_wgn(-hybrid%maxlcutm1:hybrid%maxlcutm1,-hybrid%maxlcutm1:hybrid%maxlcutm1,l,isym) )
            END DO
          ELSE
            inviop           = isym - sym%nop
            rrot(:,:,isym)   = -rrot(:,:,inviop)
            dwgn(:,:,:,isym) = dwgn(:,:,:,inviop)
            DO l = 0,hybrid%maxlcutm1
              DO m1 = -l,l
                DO m2 = -l,-1
                  cdum                = dwgn(m1, m2,l,isym)
                  dwgn(m1, m2,l,isym) = dwgn(m1,-m2,l,isym) * (-1)**m2
                  dwgn(m1,-m2,l,isym) = cdum                * (-1)**m2
                END DO
              END DO
            END DO
          END IF
        END DO
        invrrot(:,:,:sym%nop)   = rrot(:,:,sym%invtab)
        IF (sym%nsym > sym%nop) THEN
          invrrot(:,:,sym%nop+1:) = rrot(:,:,sym%invtab+sym%nop)
        END IF

        ! Get symmetry operations that leave bk(:,ikpt) invariant -> sym1
        nsym1 = 0
        DO ikpt = 1,kpts%nkpt
          isym1 = 0
          DO isym = 1,sym%nsym
            ! temporary fix until bramat_trafo is correct
            ! for systems with symmetries including translations
            IF ( isym > sym%nop ) THEN
              isym2 = isym-sym%nop
            ELSE
              isym2 = isym
            END IF
            IF ( any(sym%tau(:,isym2) /= 0) ) CYCLE

            IF(all(abs(matmul(rrot(:,:,isym),kpts%bk(:,ikpt))&
     &                                   -kpts%bk(:,ikpt)).lt.1d-12)) THEN
              isym1            = isym1 + 1
              sym1(isym1,ikpt) = isym
            END IF
          END DO
          nsym1(ikpt) = isym1
        END DO
        ! Define reduced lists of G points -> pgptm1(:,ikpt), ikpt=1,..,nkpt
        ALLOCATE ( hybrid%pgptm1(hybrid%maxgptm,kpts%nkpt),iarr(hybrid%maxgptm),&
     &             pointer(kpts%nkpt,&
     &                     minval(hybrid%gptm(1,:))-1:maxval(hybrid%gptm(1,:))+1,&
     &                     minval(hybrid%gptm(2,:))-1:maxval(hybrid%gptm(2,:))+1,&
     &                     minval(hybrid%gptm(3,:))-1:maxval(hybrid%gptm(3,:))+1))
        hybrid%pgptm1 = 0 ; iarr = 0 ; pointer = 0
        DO ikpt = 1,kpts%nkpt
          DO igpt = 1,hybrid%ngptm(ikpt)
            g = hybrid%gptm(:,hybrid%pgptm(igpt,ikpt))
            pointer(ikpt,g(1),g(2),g(3)) = igpt
          END DO
          iarr = 0
          j    = 0
          DO igpt = hybrid%ngptm(ikpt),1,-1
            IF (iarr(igpt).eq.0) THEN
              j              = j + 1
              hybrid%pgptm1(j,ikpt) = igpt
              DO isym1 = 1,nsym1(ikpt)
                g = matmul ( rrot(:,:,sym1(isym1,ikpt)) ,&
     &                       hybrid%gptm(:,hybrid%pgptm(igpt,ikpt)) )
                i = pointer(ikpt,g(1),g(2),g(3))
                IF(i.eq.0) STOP 'coulombmatrix: zero pointer (bug?)'
                iarr(i) = 1
              END DO
            END IF
          END DO
          hybrid%ngptm1(ikpt) = j
        END DO
        DEALLOCATE ( iarr )

        IF ( mpi%irank == 0 ) WRITE(6,'(12X,A)',advance='no') 'done'
        CALL cpu_time(time2)
        IF ( mpi%irank == 0 )&
     &    WRITE(6,'(2X,A,F8.2,A)') '( Timing:', time2-time1, ' )'

      ! no symmetry used
#     else 

      ALLOCATE ( hybrid%pgptm1(hybrid%maxgptm,kpts%nkpt) )
      DO ikpt = 1,kpts%nkpt
        hybrid%pgptm1(:,ikpt)    = (/ (igpt0, igpt0 = 1,hybrid%maxgptm) /)
        hybrid%ngptm1(ikpt)      = hybrid%ngptm(ikpt)
      END DO

#     endif

      ! Distribute the work as equally as possible over the processes
#     ifdef CPP_MPI
        CALL work_dist_coulomb(&
     &           mpi,kpts%nkpt,hybrid(:kpts%nkpt),&
     &           ikptmin,ikptmax,igptmin,igptmax,comm)
        calc_mt  = ( igptmax == hybrid%ngptm1(:kpts%nkpt) )
        nkminmax = ikptmax - ikptmin + 1
#     else
        ikptmin  = 1
        ikptmax  = kpts%nkpt
        igptmin  = 1
        igptmax  = hybrid%ngptm1(:kpts%nkpt)
        calc_mt  = .true.
        nkminmax = kpts%nkpt
#     endif

      IF ( mpi%irank == 0 ) WRITE(6,'(A)',advance='no') 'Preparations...'
      CALL cpu_time(time1)

      ! Define gmat (symmetric)
      i = (hybrid%lexp+1)**2
      ALLOCATE ( gmat(i,i) )
      gmat = 0
      lm1 = 0
      DO l1=0,hybrid%lexp
        DO m1=-l1,l1
          lm1 = lm1 + 1
          lm2 = 0
      lp1:DO l2=0,l1
            DO m2=-l2,l2
              lm2 = lm2 + 1
              IF(lm2.gt.lm1) EXIT lp1 ! Don't cross the diagonal!
              gmat(lm1,lm2) =&
     &          sfac(l1+l2+m2-m1)*sfac(l1+l2+m1-m2)/&
     &          ( sfac(l1+m1)*sfac(l1-m1)*sfac(l2+m2)*sfac(l2-m2) ) /&
     &          sqrt(1d0*(2*l1+1)*(2*l2+1)*(2*(l1+l2)+1))*(4*pi_const)**1.5d0
              gmat(lm2,lm1) = gmat(lm1,lm2)
            END DO
          END DO LP1
        END DO
      END DO
      ! Calculate moments of MT functions
      DO itype=1,atoms%ntype
        DO l=0,hybrid%lcutm1(itype)
          DO i=1,hybrid%nindxm1(l,itype)
            ! note that hybrid%basm1 already contains the factor rgrid
            moment(i,l,itype) =&
     &        intgrf(atoms%rmsh(:,itype)**(l+1)*hybrid%basm1(:,i,l,itype),&
     &               atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
          END DO
        END DO
        DO i =1,hybrid%nindxm1(0,itype)
          moment2(i,itype) = intgrf(atoms%rmsh(:,itype)**3*hybrid%basm1(:,i,0,itype),&
     &                          atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
        END DO
      END DO
      ! Look for different qnorm = |k+G|, definition of qnrm and pqnrm.
      CALL getnorm(kpts,hybrid%gptm,hybrid%ngptm,hybrid%pgptm,&
           &             qnrm,nqnrm,pqnrm,cell)
      ALLOCATE ( sphbesmoment(0:hybrid%lexp,atoms%ntype,nqnrm),&
     &           olap(hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype,nqnrm),&
     &           integral(hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype,nqnrm) )
      sphbes        = 0
      sphbesmoment  = 0
      sphbesmoment1 = 0
      olap          = 0
      integral      = 0

      ! Calculate moments of spherical Bessel functions (for (2) and (3))              (->sphbesmoment)
      ! Calculate overlap of spherical Bessel functions with basis functions (for (2)) (->olap)
      ! Calculate overlap of sphbesmoment1(r,l)         with basis functions (for (2)) (->integral)
      ! We use           sphbes(r,l) = j_l(qr)
      ! and       sphbesmoment1(r,l) = 1/r**(l-1) * INT(0..r) r'**(l+2) * j_l(qr') dr'
      !                                + r**(l+2) * INT(r..S) r'**(1-l) * j_l(qr') dr' .

      iqnrmstart = mpi%irank + 1
      iqnrmstep  = mpi%isize

      DO iqnrm = iqnrmstart,nqnrm,iqnrmstep
        qnorm = qnrm(iqnrm)
        DO itype = 1,atoms%ntype
          ng            = atoms%jri(itype)
          rdum          = atoms%rmt(itype)
          sphbes        = 0
          sphbesmoment1 = 0 
          IF(qnorm.eq.0) THEN
            sphbesmoment(0,itype,iqnrm) = rdum**3 / 3
            DO i = 1,ng
              sphbes(i,0)        = 1
              sphbesmoment1(i,0) = atoms%rmsh(i,itype)**2 / 3&
     &                           + ( rdum**2 - atoms%rmsh(i,itype)**2 ) / 2
            END DO
          ELSE
            CALL sphbessel(rarr,qnorm*rdum,hybrid%lexp+1)
            DO l = 0,hybrid%lexp
              sphbesmoment(l,itype,iqnrm)&
     &                                 = rdum**(l+2) * rarr(l+1) / qnorm
            END DO
            DO i = ng,1,-1
              rdum = atoms%rmsh(i,itype)
              CALL sphbessel(rarr,qnorm*rdum,hybrid%lcutm1(itype)+1)
              DO l = 0,hybrid%lcutm1(itype)
                sphbes(i,l) = rarr(l)
                IF(l.ne.0) THEN ; rdum1 = -rdum**(1-l) * rarr(l-1)
                ELSE            ; rdum1 = -cos(qnorm*rdum) / qnorm
                ENDIF
                IF(i.eq.ng)  rarr1(l) = rdum1
                sphbesmoment1(i,l) = (rdum**(l+2)*rarr(l+1)/rdum**(l+1)&
     &                        + ( rarr1(l) - rdum1 ) * rdum**l ) / qnorm
              END DO
            END DO
          END IF
          DO l = 0,hybrid%lcutm1(itype)
            DO n = 1,hybrid%nindxm1(l,itype)
              ! note that hybrid%basm1 already contains one factor rgrid
              olap(n,l,itype,iqnrm)     =  &
     &        intgrf(atoms%rmsh(:,itype)*hybrid%basm1(:,n,l,itype)*sphbes(:,l),&
     &               atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

              integral(n,l,itype,iqnrm) = &
     &        intgrf(atoms%rmsh(:,itype)*hybrid%basm1(:,n,l,itype)*sphbesmoment1(:,l),&
     &               atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

            END DO
          END DO
        END DO
      END DO
#     ifdef CPP_MPI
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,sphbesmoment,size(sphbesmoment),&
     &                            MPI_REAL,MPI_SUM,mpi,ierr)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,        olap,size(    olap    ),&
     &                            MPI_REAL,MPI_SUM,mpi,ierr)
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,    integral,size(  integral  ),&
     &                            MPI_REAL,MPI_SUM,mpi,ierr)
#     endif

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(18X,A)',advance='no') 'done'
        CALL cpu_time(time2) 
        WRITE(6,'(2X,A,F8.2,A)',advance='no')&
     &       '( Timing:',time2-time1,' )'
        WRITE(6,*)
      END IF

!
!     (1) Case < MT | v | MT >
!

      IF( mpi%irank == 0 )&
     &  WRITE(6,'(A)',advance='no') '< MT | v | MT > contribution...'

      CALL cpu_time(time1)

      IF ( ANY( calc_mt ) ) THEN

!       (1a) r,r' in same MT

        ix  = 0
        iy  = 0
        iy0 = 0
        DO itype=1,atoms%ntype
          DO ineq=1,atoms%neq(itype) 
            ! Here the diagonal block matrices do not depend on ineq. In (1b) they do depend on ineq, though,
            DO l=0,hybrid%lcutm1(itype)
              DO n2=1,hybrid%nindxm1(l,itype)
                ! note that hybrid%basm1 already contains the factor rgrid
                CALL primitivef(primf1,hybrid%basm1(:,n2,l,itype)&
     &               *atoms%rmsh(:,itype)**(l+1),atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,itype,atoms%ntype)
                ! -itype is to enforce inward integration
                CALL primitivef(primf2,hybrid%basm1(:,n2,l,itype)&
     &               /atoms%rmsh(:,itype)**l,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,-itype,atoms%ntype)

                primf1 = primf1 / atoms%rmsh(:,itype)**l
                primf2 = primf2 * atoms%rmsh(:,itype)**(l+1)

                DO n1=1,n2
                  integrand = hybrid%basm1(:,n1,l,itype) * (primf1 + primf2)
!                 call intgr0( (4*pimach())/(2*l+1)*integrand,rmsh(1,itype),dx(itype),jri(itype),mat(n2*(n2-1)/2+n1) )
                  mat(n2*(n2-1)/2+n1) = (4*pi_const)/(2*l+1)&
     &                         * intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,&
     &                                  atoms%ntype,itype,gridf)
                END DO
              END DO

              ! distribute mat for m=-l,l on coulomb in block-matrix form
              DO M=-l,l
                DO n2=1,hybrid%nindxm1(l,itype)
                  ix = ix + 1
                  iy = iy0
                  DO n1=1,n2
                    iy                   = iy + 1
                    i                    = ix*(ix-1)/2+iy
                    j                    = n2*(n2-1)/2+n1
                    coulomb(i,kpts%nkpt) = mat(j)
                  END DO
                END DO
                iy0 = ix
              END DO

            END DO
          END DO
        END DO

!       (1b) r,r' in different MT

        ALLOCATE( coulmat(nbasp,nbasp), stat=ok)
        IF( ok .ne. 0 ) STOP 'coulombmatrix: failure allocation coulmat'
        coulmat = 0

      END IF

      DO ikpt=ikptmin,ikptmax

        ! only the first rank handles the MT-MT part
        IF ( calc_mt(ikpt) ) THEN

          ix  = 0
          ic2 = 0
          DO itype2=1,atoms%ntype
            DO ineq2=1,atoms%neq(itype2)
              ic2 = ic2 + 1
              lm2 = 0
              DO l2=0,hybrid%lcutm1(itype2)
                DO m2=-l2,l2
                  lm2 = lm2 + 1
                  DO n2=1,hybrid%nindxm1(l2,itype2)
                    ix  = ix + 1

                    iy  = 0
                    ic1 = 0
               lp2: DO itype1=1,itype2
                      DO ineq1=1,atoms%neq(itype1)
                        ic1 = ic1 + 1
                        lm1 = 0
                        DO l1=0,hybrid%lcutm1(itype1)
                          DO m1=-l1,l1
                            lm1 = lm1 + 1
                            DO n1=1,hybrid%nindxm1(l1,itype1)
                              iy = iy + 1
                              IF(iy.gt.ix) EXIT lp2 ! Don't cross the diagonal!
                              rdum            = (-1)**(l2+m2)*&
     & moment(n1,l1,itype1)*moment(n2,l2,itype2)*gmat(lm1,lm2)
                              l               = l1 + l2
                              lm              = l**2 + l + m1 - m2 + 1
                              idum            = ix*(ix-1)/2+iy
                              coulmat(iy,ix)  = coulomb(idum,kpts%nkpt)&
     & + exp(img* 2*pi_const * dot_product(kpts%bk(:,ikpt),&
     &                   atoms%taual(:,ic2)-atoms%taual(:,ic1)))&
     & *rdum * structconst(lm,ic1,ic2,ikpt)

                              coulmat(ix,iy)  = conjg(coulmat(iy,ix))
                            END DO
                          END DO
                        END DO
                      END DO
                    END DO lp2

                  END DO
                END DO
              END DO
            END DO
          END DO

#if ( defined(CPP_INVERSION) )
          !symmetrize makes the Coulomb matrix real symmetric
          CALL symmetrize(coulmat,nbasp,nbasp,3,.false.,&
     &                    atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                    hybrid%nindxm1,sym)
#endif

          coulomb(:nbasp*(nbasp+1)/2,ikpt) = packmat(coulmat)

        END IF

      END DO
      IF ( ANY( calc_mt ) ) DEALLOCATE( coulmat )

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(2X,A)',advance='no') 'done'
        CALL cpu_time(time2) 
        WRITE(6,'(2X,A,F8.2,A)',advance='no')&
     &       '( Timing:',time2-time1,' )'
        WRITE(6,*)
      END IF

      IF(hybrid%maxgptm.eq.0) GOTO 1 ! skip calculation of plane-wave contribution if mixed basis does not contain plane waves

!
!     (2) Case < MT | v | PW >
!

      IF( mpi%irank == 0 )&
     &  WRITE(6,'(A)',advance='no') '< MT | v | PW > contribution...'

      CALL cpu_time(time1)

!     (2a) r in MT, r' everywhere
!     (2b) r,r' in same MT
!     (2c) r,r' in different MT

      ALLOCATE( coulmat(nbasp,hybrid%maxgptm), stat=ok )
      IF( ok .ne. 0 ) STOP 'coulombmatrix: failure allocation coulmat'
      coulmat = 0

      DO ikpt = ikptmin,ikptmax !1,kpts%nkpt

        coulmat = 0

        ! start to loop over interstitial plane waves
        DO igpt0 = igptmin(ikpt),igptmax(ikpt) !1,hybrid%ngptm1(ikpt)
          igpt  = hybrid%pgptm1(igpt0,ikpt)
          igptp = hybrid%pgptm(igpt,ikpt)
          ix    = nbasp + igpt
          q     = matmul ( kpts%bk(:,ikpt) + hybrid%gptm(:,igptp), cell%bmat )
          qnorm = sqrt(sum(q**2))
          iqnrm = pqnrm(igpt,ikpt)
          IF(abs(qnrm(iqnrm)-qnorm).gt.1d-12) STOP 'coulombmatrix: qnorm does not equal corresponding & element in qnrm (bug?)' ! We shouldn't stop here!

          CALL harmonicsr(y1,matmul(kpts%bk(:,kpts%nkpt),cell%bmat),2)
          CALL harmonicsr(y2,matmul(hybrid%gptm(:,igptp)     ,cell%bmat),2)
          CALL harmonicsr(y,q,hybrid%lexp)
          y1 = conjg(y1) ; y2 = conjg(y2) ; y = conjg(y)

          iy = 0
          ic = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              ic = ic + 1
              lm = 0
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  lm = lm + 1

                  ! calculate sum over lm and centers for (2c) -> csum, csumf
                  csum  = 0
                  csumf = 0
                  ic1   = 0
                  DO itype1=1,atoms%ntype
                    DO ineq1=1,atoms%neq(itype1)
                      ic1  = ic1 + 1
                      cexp = 4*pi_const * exp( img * 2*pi_const&
                           * ( dot_product( kpts%bk(:,ikpt)+hybrid%gptm(:,igptp),atoms%taual(:,ic1) )&
                           - dot_product( kpts%bk(:,ikpt),atoms%taual(:,ic)  ) ) )

                      lm1 = 0
                      DO l1=0,hybrid%lexp
                        l2   = l + l1 ! for structconst
                        idum = 1
                        cdum = sphbesmoment(l1,itype1,iqnrm) * imgl(l1) * cexp
                        DO m1=-l1,l1
                          lm1  = lm1 + 1
                          m2   = M - m1              ! for structconst
                          lm2  = l2**2 + l2 + m2 + 1 !
                          csum = csum - idum * gmat(lm1,lm) * y(lm1) * cdum * structconst(lm2,ic,ic1,ikpt)
                          idum = -idum ! factor (-1)*(l1+m1)
                        END DO
                      END DO

                      ! add contribution of (2c) to csum and csumf coming from linear and quadratic orders of Y_lm*(G) / G * j_(l+1)(GS)
                      IF(ikpt.eq.1.and.l.le.2) THEN
                        cexp      = exp(img*2*pi_const * dot_product(hybrid%gptm(:,igptp),atoms%taual(:,ic1)))&
                             * gmat(lm,1) * 4*pi_const/cell%vol
                        csumf(lm) = csumf(lm) - cexp * sqrt(4*pi_const) *&
                             img**l * sphbesmoment(0,itype1,iqnrm) / facfac(l-1)
                        IF(l.eq.0) THEN
                          IF(igpt.ne.1) THEN
                            csum = csum - cexp * ( sphbesmoment(0,itype1,iqnrm)*atoms%rmt(itype1)**2 -&
                                 sphbesmoment(2,itype1,iqnrm)*2d0/3 ) / 10
                          ELSE
                            csum = csum - cexp * atoms%rmt(itype1)**5/30
                          END IF
                        ELSE IF(l.eq.1) THEN
                          csum = csum + cexp * img * sqrt(4*pi_const)&
                               * sphbesmoment(1,itype1,iqnrm) * y(lm) / 3
                        END IF
                      END IF

                    END DO
                  END DO

                  ! add contribution of (2a) to csumf
                  IF(ikpt.eq.1.and.igpt.eq.1.and.l.le.2) THEN
                    csumf(lm)=csumf(lm) + (4*pi_const)**2 * img**l / facfac(l)
                  END IF

                  ! finally define coulomb
                  idum = ix*(ix-1)/2
                  cdum = (4*pi_const)**2 * imgl(l) * y(lm) * exp(img * 2*pi_const&
                       * dot_product(hybrid%gptm(:,igptp),atoms%taual(:,ic)))
                  DO n=1,hybrid%nindxm1(l,itype)
                    iy = iy + 1

                    IF(ikpt.eq.1.and.igpt.eq.1) THEN
                       IF(l.eq.0) coulmat(iy,ix-nbasp) =&
                            - cdum * moment2(n,itype) / 6 / svol         ! (2a)
                       coulmat(iy,ix-nbasp)   = coulmat(iy,ix-nbasp)&
                            + ( - cdum / (2*l+1) * integral(n,l,itype,iqnrm)& ! (2b)&
                            + csum * moment(n,l,itype) ) / svol          ! (2c)
                    ELSE
                       coulmat(iy,ix-nbasp)   = &
                            (   cdum * olap(n,l,itype,iqnrm) / qnorm**2  &  ! (2a)&
                            - cdum / (2*l+1) * integral(n,l,itype,iqnrm)& ! (2b)&
                            + csum * moment(n,l,itype) ) / svol          ! (2c)
                       
                    END IF

                  END DO

                END DO
              END DO
            END DO

          END DO
        END DO

#if ( defined(CPP_INVERSION) )
        CALL symmetrize(coulmat,nbasp,hybrid(ikpt),1,.false.,&
     &                  atoms,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                  hybrid%nindxm1,sym)
#endif

        M = nbasp*(nbasp+1)/2
        DO i=1,hybrid%ngptm(ikpt)
          DO j=1,nbasp+i
            M = M + 1
            IF(j.le. nbasp) coulomb(M,ikpt) = coulmat(j,i)
          END DO
        END DO

      END DO 

      DEALLOCATE( coulmat,olap,integral )

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(2X,A)',advance='no') 'done'
        CALL cpu_time(time2) 
        WRITE(6,'(2X,A,F8.2,A)') '( Timing:',time2-time1,' )'
      END IF

!
!     (3) Case < PW | v | PW >
!

      IF( mpi%irank == 0 )&
     &  WRITE(6,'(A)',advance='no') '< PW | v | PW > contribution...'

      CALL cpu_time(time1)

!     (3a) r,r' everywhere; r everywhere, r' in MT; r in MT, r' everywhere

      CALL cpu_time(time1)
      ! Calculate the hermitian matrix smat(i,j) = sum(a) integral(MT(a)) exp[i(Gj-Gi)r] dr
      ALLOCATE ( smat(hybrid%gptmd,hybrid%gptmd) )
      smat = 0
      DO igpt2=1,hybrid%gptmd
        DO igpt1=1,igpt2
          g     = hybrid%gptm(:,igpt2)-hybrid%gptm(:,igpt1)
          gnorm = gptnorm(g,cell%bmat)
          IF(gnorm.eq.0) THEN
            DO itype=1,atoms%ntype
              smat(igpt1,igpt2) = smat(igpt1,igpt2) + atoms%neq(itype) * 4*pi_const*atoms%rmt(itype)**3/3
            END DO
          ELSE
            ic = 0
            DO itype=1,atoms%ntype
              rdum = atoms%rmt(itype) * gnorm
              rdum = 4*pi_const * ( sin(rdum) - rdum * cos(rdum) ) / gnorm**3
              DO ineq=1,atoms%neq(itype)
                ic                = ic + 1
                smat(igpt1,igpt2) = smat(igpt1,igpt2)&
                     + rdum * exp( img * 2*pi_const * dot_product(atoms%taual(:,ic),g) )
              END DO
            END DO
          END IF
          smat(igpt2,igpt1) = conjg(smat(igpt1,igpt2))
        END DO
      END DO

      ! Coulomb matrix, contribution (3a)
      DO ikpt=ikptmin,ikptmax

        DO igpt0=igptmin(ikpt),igptmax(ikpt)
          igpt2  = hybrid%pgptm1(igpt0,ikpt)
          igptp2 = hybrid%pgptm(igpt2,ikpt)
          ix     = nbasp + igpt2
          iy     = nbasp
          q2     = matmul ( kpts%bk(:,ikpt) + hybrid%gptm(:,igptp2) , cell%bmat )
          rdum2  = sum(q2**2)
          IF( rdum2 .ne. 0 ) rdum2 = 4*pi_const/rdum2

          DO igpt1=1,igpt2
            igptp1 = hybrid%pgptm(igpt1,ikpt)
            iy     = iy + 1
            q1     = matmul ( kpts%bk(:,ikpt) + hybrid%gptm(:,igptp1) , cell%bmat )
            idum   = ix*(ix-1)/2+iy
            rdum1  = sum(q1**2)
            IF( rdum1 .ne. 0 ) rdum1 = 4*pi_const/rdum1

            IF(ikpt.eq.1) THEN
              IF(igpt1.ne.1) THEN
                coulomb(idum,1) = - smat(igptp1,igptp2) * rdum1 / cell%vol
              END IF
              IF(igpt2.ne.1) then
                coulomb(idum,1) = coulomb(idum,1) - smat(igptp1,igptp2) * rdum2 / cell%vol
              END IF
            ELSE
              coulomb(idum,ikpt) = - smat(igptp1,igptp2) * ( rdum1 + rdum2 ) / cell%vol
            END IF
          END DO
          IF(ikpt.ne.1.or.igpt2.ne.1) THEN                  !
            coulomb(idum,ikpt) = coulomb(idum,ikpt) + rdum2 ! diagonal term
          END IF                                            !
        END DO

      END DO

!     (3b) r,r' in different MT

      DO ikpt=ikptmin,ikptmax!1,kpts%nkpt

        ! group together quantities which depend only on l,m and igpt -> carr2a
        ALLOCATE( carr2a((hybrid%lexp+1)**2,hybrid%maxgptm),carr2b(atoms%nat,hybrid%maxgptm) )
        carr2a = 0 ; carr2b = 0
        DO igpt=1,hybrid%ngptm(ikpt)
          igptp = hybrid%pgptm(igpt,ikpt)
          iqnrm = pqnrm(igpt,ikpt)
          q     = matmul ( kpts%bk(:,ikpt) + hybrid%gptm(:,igptp),cell%bmat)
          CALL harmonicsr(y,q,hybrid%lexp)
          y     = conjg(y)
          lm = 0
          DO l=0,hybrid%lexp
            DO M=-l,l
              lm              = lm + 1
              carr2a(lm,igpt) = 4*pi_const * imgl(l) * y(lm)
            END DO
          END DO
          DO ic = 1,atoms%nat
            carr2b(ic,igpt) = exp ( -img * 2*pi_const * &
                 dot_product(kpts%bk(:,ikpt)+hybrid%gptm(:,igptp),atoms%taual(:,ic)) )
          END DO
        END DO

        !finally we can loop over the plane waves (G: igpt1,igpt2)
        ALLOCATE ( carr2(atoms%nat,(hybrid%lexp+1)**2),&
             structconst1(atoms%nat,(2*hybrid%lexp+1)**2) )
        carr2 = 0 ; structconst1 = 0
        DO igpt0=igptmin(ikpt),igptmax(ikpt)!1,hybrid%ngptm1(ikpt)
          igpt2  = hybrid%pgptm1(igpt0,ikpt)
          ix     = nbasp + igpt2
          igptp2 = hybrid%pgptm(igpt2,ikpt)
          iqnrm2 = pqnrm(igpt2,ikpt)
          ic2    = 0
          carr2  = 0
          DO itype2 = 1,atoms%ntype
            DO ineq2 = 1,atoms%neq(itype2)
              ic2   = ic2 + 1
              cexp  = conjg ( carr2b(ic2,igpt2) )
              lm2   = 0
              DO ic1 = 1,atoms%nat
                structconst1(ic1,:) = structconst(:,ic1,ic2,ikpt)
              END DO
              DO l2 = 0,hybrid%lexp
                idum = 1
                DO m2 = -l2,l2
                  lm2  = lm2 + 1
                  cdum = idum * sphbesmoment(l2,itype2,iqnrm2) * cexp * carr2a(lm2,igpt2)
                  IF( cdum .ne. 0 ) THEN
                    lm1 = 0
                    DO l1 = 0,hybrid%lexp
                      l  =  l1 + l2
                      M  = -l1 - m2 !first loop of m1
                      lm = l**2 + l + M
                      DO m1 = -l1,l1
                        lm1  = lm1 + 1
                        lm   = lm  + 1
                        cdum1= cdum * gmat(lm1,lm2)
                        DO ic1 = 1,atoms%nat
                          carr2(ic1,lm1) = carr2(ic1,lm1) + cdum1 * structconst1(ic1,lm)
                        END DO
                      END DO
                    END DO
                  END IF
                  idum = -idum !factor (-1)**(l+M)
                END DO
              END DO
            END DO
          END DO

          iy = nbasp
          DO igpt1=1,igpt2
            iy      = iy + 1
            igptp1  = hybrid%pgptm(igpt1,ikpt)
            iqnrm1  = pqnrm(igpt1,ikpt)
            csum    = 0
            ic      = 0
            DO itype=1,atoms%ntype
              DO ineq=1,atoms%neq(itype)
                ic   = ic  + 1
                cexp = carr2b(ic,igpt1)
                lm   = 0
                DO l=0,hybrid%lexp
                  cdum = cexp * sphbesmoment(l,itype,iqnrm1)
                  DO M=-l,l
                    lm   = lm + 1
                    csum = csum + cdum * carr2(ic,lm) * conjg ( carr2a(lm,igpt1) ) ! for coulomb
                  END DO
                END DO
              END DO
            END DO
            idum = ix*(ix-1)/2+iy
            coulomb(idum,ikpt) = coulomb(idum,ikpt) + csum / cell%vol
          END DO
        END DO
        DEALLOCATE( carr2,carr2a,carr2b,structconst1 )
      END DO !ikpt

!     Add corrections from higher orders in (3b) to coulomb(:,1)
      ! (1) igpt1 > 1 , igpt2 > 1  (finite G vectors)
      rdum = (4*pi_const)**(1.5d0)/cell%vol**2 * gmat(1,1)
      DO igpt0 = 1,hybrid%ngptm1(1)
        igpt2  = hybrid%pgptm1(igpt0,1) ; IF ( igpt2 == 1 ) CYCLE
        ix     = nbasp + igpt2
        iqnrm2 = pqnrm(igpt2,1)
        igptp2 = hybrid%pgptm(igpt2,1)
        q2     = matmul(hybrid%gptm(:,igptp2),cell%bmat)
        qnorm2 = sqrt(sum(q2**2))
        iy     = nbasp + 1
        DO igpt1 = 2,igpt2
          iy     = iy + 1
          idum   = ix*(ix-1)/2+iy
          iqnrm1 = pqnrm(igpt1,1)
          igptp1 = hybrid%pgptm(igpt1,1)
          q1     = matmul(hybrid%gptm(:,igptp1),cell%bmat)
          qnorm1 = sqrt(sum(q1**2))
          rdum1  = dot_product(q1,q2) / (qnorm1*qnorm2)
          ic1    = 0
          DO itype1 = 1,atoms%ntype
            DO ineq1 = 1,atoms%neq(itype1)
              ic1 = ic1 + 1
              ic2 = 0
              DO itype2 = 1,atoms%ntype
                DO ineq2 = 1,atoms%neq(itype2)
                  ic2  = ic2 + 1
                  cdum = exp ( img * 2*pi_const *&
     &                 ( - dot_product(hybrid%gptm(:,igptp1),atoms%taual(:,ic1))&
     &                   + dot_product(hybrid%gptm(:,igptp2),atoms%taual(:,ic2)) ) )
                  coulomb(idum,1) = coulomb(idum,1) + rdum * cdum * (&
     &              - sphbesmoment(1,itype1,iqnrm1) &
     &                * sphbesmoment(1,itype2,iqnrm2) * rdum1  / 3&
     &              - sphbesmoment(0,itype1,iqnrm1)&
     &                * sphbesmoment(2,itype2,iqnrm2)          / 6&
     &              - sphbesmoment(2,itype1,iqnrm1)&
     &                * sphbesmoment(0,itype2,iqnrm2)          / 6 &
     &              + sphbesmoment(0,itype1,iqnrm1)&
     &                * sphbesmoment(1,itype2,iqnrm2) / qnorm2 / 2&
     &              + sphbesmoment(1,itype1,iqnrm1)&
     &                * sphbesmoment(0,itype2,iqnrm2) / qnorm1 / 2 )
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

      ! (2) igpt1 = 1 , igpt2 > 1  (first G vector vanishes, second finite)
      iy = nbasp + 1
      DO igpt0 = 1,hybrid%ngptm1(1)
        igpt2  = hybrid%pgptm1(igpt0,1) ; IF ( igpt2 == 1 ) CYCLE
        ix     = nbasp + igpt2
        iqnrm2 = pqnrm(igpt2,1)
        igptp2 = hybrid%pgptm(igpt2,1)
        qnorm2 = qnrm(iqnrm2)
        idum   = ix*(ix-1)/2+iy
        DO itype1 = 1,atoms%ntype
          DO ineq1 = 1,atoms%neq(itype1)
            ic2 = 0
            DO itype2 = 1,atoms%ntype
              DO ineq2 = 1,atoms%neq(itype2)
                ic2  = ic2 + 1
                cdum = exp ( img * 2*pi_const&
     &               * dot_product(hybrid%gptm(:,igptp2),atoms%taual(:,ic2)) )
                coulomb(idum,1) = coulomb(idum,1)&
     &          + rdum * cdum * atoms%rmt(itype1)**3 * (&
     &            + sphbesmoment(0,itype2,iqnrm2) / 30 * atoms%rmt(itype1)**2&
     &            - sphbesmoment(2,itype2,iqnrm2) / 18 &
     &            + sphbesmoment(1,itype2,iqnrm2) /  6 / qnorm2 )
              END DO
            END DO
          END DO
        END DO
      END DO

      ! (2) igpt1 = 1 , igpt2 = 1  (vanishing G vectors)
      iy   = nbasp + 1
      ix   = nbasp + 1
      idum = ix*(ix-1)/2+iy
      DO itype1 = 1,atoms%ntype
        DO ineq1 = 1,atoms%neq(itype1)
          DO itype2 = 1,atoms%ntype
            DO ineq2 = 1,atoms%neq(itype2)
              coulomb(idum,1) = coulomb(idum,1)&
     &          + rdum * atoms%rmt(itype1)**3 * atoms%rmt(itype2)**3 *&
     &          ( atoms%rmt(itype1)**2 + atoms%rmt(itype2)**2 ) / 90
            END DO
          END DO
        END DO
      END DO


!     (3c) r,r' in same MT

      ! Calculate sphbesintegral
      ALLOCATE ( sphbes0(-1:hybrid%lexp+2,atoms%ntype,nqnrm),&
     &           carr2((hybrid%lexp+1)**2,hybrid%maxgptm) )
      sphbes0 = 0 ; carr2 = 0
      DO iqnrm = 1,nqnrm
        DO itype = 1,atoms%ntype
          rdum = qnrm(iqnrm) * atoms%rmt(itype)
          CALL sphbessel(sphbes0(0,itype,iqnrm),rdum,hybrid%lexp+2)
          IF( rdum.ne.0 ) sphbes0(-1,itype,iqnrm) = cos(rdum)/rdum
        END DO
      END DO

      l_warn = ( mpi%irank == 0 )
      DO ikpt=ikptmin,ikptmax!1,nkpt

        DO igpt = 1,hybrid%ngptm(ikpt)
          igptp = hybrid%pgptm(igpt,ikpt)
          q     = matmul ( kpts%bk(:,ikpt) + hybrid%gptm(:,igptp), cell%bmat )
          CALL harmonicsr(carr2(:,igpt),q,hybrid%lexp)
        END DO

        DO igpt0=igptmin(ikpt),igptmax(ikpt)!1,hybrid%ngptm1(ikpt)
          igpt2  = hybrid%pgptm1(igpt0,ikpt)
          ix     = nbasp + igpt2
          igptp2 = hybrid%pgptm(igpt2,ikpt)
          iqnrm2 = pqnrm(igpt2,ikpt)
          q2     = matmul (kpts%bk(:,ikpt) + hybrid%gptm(:,igptp2),cell%bmat)
          y2     = conjg ( carr2(:,igpt2) )
          iy     = nbasp
          DO igpt1=1,igpt2
            iy     = iy + 1
            igptp1 = hybrid%pgptm(igpt1,ikpt)
            iqnrm1 = pqnrm(igpt1,ikpt)
            q1     = matmul (kpts%bk(:,ikpt) + hybrid%gptm(:,igptp1),cell%bmat)
            y1     = carr2(:,igpt1)
            cexp1  = 0
            ic     = 0
            DO itype=1,atoms%ntype
              DO ineq=1,atoms%neq(itype)
                ic           = ic + 1
                cexp1(itype) = cexp1(itype) +&
     &                         exp(img * 2*pi_const * dot_product(&
     &                    (hybrid%gptm(:,igptp2)-hybrid%gptm(:,igptp1)),atoms%taual(:,ic)) )
              ENDDO
            ENDDO
            lm   = 0
            cdum = 0
            DO l=0,hybrid%lexp
              cdum1 = 0
              DO itype=1,atoms%ntype
                cdum1  = cdum1 + cexp1(itype)*sphbessel_integral(&
     &                                       atoms,itype,qnrm,nqnrm,&
     &                                       iqnrm1,iqnrm2,l,hybrid,&
     &                                       sphbes0,l_warn,l_warned)&
     &                         / (2*l+1)
                l_warn = l_warn .and. .not. l_warned ! only warn once
              END DO
              DO M=-l,l
                lm   = lm + 1
                cdum = cdum + cdum1 * y1(lm) * y2(lm)
              ENDDO
            ENDDO
            idum               = ix*(ix-1)/2+iy
            coulomb(idum,ikpt) = coulomb(idum,ikpt)+(4*pi_const)**3*cdum / cell%vol
          END DO
        END DO

      END DO

      DEALLOCATE( carr2 )

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(2X,A)',advance='no') 'done'
        CALL cpu_time(time2)
        WRITE(6,'(2X,A,F8.2,A)') '( Timing:',time2-time1,' )'
      END IF

!
!     Symmetry-equivalent G vectors
!
#     ifndef CPP_NOCOULSYM

        IF ( mpi%irank == 0 )&
     &    WRITE(6,'(A)',advance='no') 'Symm.-equiv. matrix elements...'
        CALL cpu_time(time1)
        ! All elements are needed so send all data to all processes treating the
        ! respective k-points
#       ifdef CPP_MPI
          idisp = nbasp*(nbasp+1)/2
          DO ikpt = ikptmin,ikptmax
            icnt = nbasm1(ikpt)*(nbasm1(ikpt)+1)/2 - idisp
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,coulomb(idisp+1:,ikpt),icnt,&
     &                         datatype,MPI_SUM,comm(ikpt),ierr)
          END DO
#       endif

        ALLOCATE ( carr2(hybrid%maxbasm1,2),iarr(hybrid%maxgptm) )
        ALLOCATE ( nsym_gpt(hybrid%gptmd,kpts%nkpt),&
     &             sym_gpt(maxval(nsym1),hybrid%gptmd,kpts%nkpt) )
        nsym_gpt = 0 ; sym_gpt = 0
        DO ikpt = ikptmin,ikptmax
          carr2 = 0 ; iarr = 0
          iarr(hybrid%pgptm1(:hybrid%ngptm1(ikpt),ikpt)) = 1
          DO igpt0 = 1,hybrid%ngptm1(ikpt) !igptmin(ikpt),igptmax(ikpt)
            lsym         = ( ( igptmin(ikpt) <= igpt0 ) .AND.&
     &                       ( igptmax(ikpt) >= igpt0 ) )
            igpt2        = hybrid%pgptm1(igpt0,ikpt)
            j            = (nbasp+igpt2-1) * (nbasp+igpt2) / 2
            i            = nbasp+igpt2
            carr2(1:i,2) = coulomb(j+1:j+i,ikpt)
            j            = j + i
            DO i = nbasp+igpt2+1,nbasm1(ikpt)
              j          = j + i - 1
#             ifdef CPP_INVERSION
                    carr2(i,2) =  coulomb(j,ikpt)
#             else
                    carr2(i,2) = conjg( coulomb(j,ikpt) )
#             endif
            END DO
            IF ( lsym ) THEN
              ic                     = 1
              sym_gpt(ic,igpt0,ikpt) = igpt2
            END IF
            DO isym1 = 2,nsym1(ikpt)
              isym  = sym1(isym1,ikpt)
              CALL bramat_trafo(&
                   carr2(:,1),igpt1,&
                   carr2(:,2),igpt2,ikpt,isym,.false.,pointer(ikpt,:,:,:),&
                   sym,rrot(:,:,isym),invrrot(:,:,isym),hybrid,&
                   kpts,hybrid%maxlcutm1,atoms,hybrid%lcutm1,&
                   hybrid%nindxm1,hybrid%maxindxm1,dwgn(:,:,:,isym),&
                   nbasp,nbasm1)
              IF(iarr(igpt1).eq.0) THEN
                CALL bramat_trafo(&
     &            carr2(:,1),igpt1,&
     &            carr2(:,2),igpt2,ikpt,isym,.true.,pointer(ikpt,:,:,:),&
     &            sym,rrot(:,:,isym),invrrot(:,:,isym),hybrid,&
     &            kpts,hybrid%maxlcutm1,atoms,hybrid%lcutm1,&
     &            hybrid%nindxm1,hybrid%maxindxm1,&
     &            dwgn(:,:,:,isym),nbasp,nbasm1)
                l = (nbasp+igpt1-1) * (nbasp+igpt1) / 2
                coulomb(l+1:l+nbasp+igpt1,ikpt) = carr2(:nbasp+igpt1,1)
                iarr(igpt1) = 1
                IF ( lsym ) THEN
                  ic                     = ic + 1
                  sym_gpt(ic,igpt0,ikpt) = igpt1
                END IF
              END IF
            END DO
            nsym_gpt(igpt0,ikpt) = ic
          END DO ! igpt0
        END DO ! ikpt
        DEALLOCATE ( carr2,iarr,hybrid%pgptm1 )

        IF ( mpi%irank == 0 ) THEN
          WRITE(6,'(2X,A)',advance='no') 'done'
          CALL cpu_time(time2)
          WRITE(6,'(2X,A,F8.2,A)') '( Timing:',time2-time1,' )'
        END IF

      ! no symmetry used
#     else

        ALLOCATE( nsym_gpt(hybrid%gptmd,kpts%nkpt),sym_gpt(1,hybrid%gptmd,kpts%nkpt) )
        nsym_gpt = 1
        DO ikpt = 1,kpts%nkpt
          sym_gpt(1,:,ikpt) = (/ (igpt0, igpt0=1,hybrid%gptmd) /)
        END DO

#     endif

  1   DEALLOCATE (qnrm,pqnrm)

      CALL cpu_time(time1)
      IF ( xcpot%icorr == icorr_hse .or. xcpot%icorr == icorr_vhse ) THEN
        !
        ! The HSE functional is realized subtracting erf/r from
        ! the normal Coulomb matrix
        !
#ifdef CPP_NOSPMVEC
        CALL change_coulombmatrix(&
     &     atoms,kpts,kpts,kpts%nkpt,&
     &     cell,cell,hybrid%lcutm1,hybrid%maxlcutm1,&
     &     hybrid%nindxm1,hybrid%maxindxm1,hybrid,&
     &     hybrid%basm1,hybrid%maxbasm1,nbasm1,sym,mpi,&
     &     coulomb)
#endif
      ELSE
        IF ( ikptmin == 1 ) CALL subtract_sphaverage()
      END IF

      ! transform Coulomb matrix to the biorthogonal set
      IF ( mpi%irank == 0 ) WRITE(6,'(A)',advance='no') 'Transform to biorthogonal set...'
      CALL cpu_time(time1)
      DO ikpt = ikptmin,ikptmax

        !calculate IR overlap-matrix
        ALLOCATE( olapm(hybrid%ngptm(ikpt),hybrid%ngptm(ikpt)) )
        olapm = 0

        CALL olap_pw(olapm,hybrid%gptm(:hybrid%ngptm(ikpt),ikpt),hybrid%ngptm(ikpt), atoms,cell )

!         !calculate eigenvalues of olapm
!         ALLOCATE( eval(ngptm(ikpt)),evec(ngptm(ikpt),ngptm(ikpt)) )
!         CALL diagonalize(evec,eval,olapm)
!
!         !
!         ! small eigenvalues lead to inaccuries in the inversion
!         ! however it seems that these do not play an important role
!         ! for the HF exchange
!         ! thus we do not do a SingularValueDecomposition
!         !
!
!         IF( any(eval .le. 1E-06 ) ) THEN
! !           WRITE(*,*) count( eval .le. 1E-06 )
! !           WRITE(*,*) 'eval .le. 1E-06'
! !           ALLOCATE( involapm(ngptm(ikpt),ngptm(ikpt)) )
!           olapm = 0
!           DO i = 1,ngptm(ikpt)
!             IF( eval(i) .le. 1E-06) CYCLE
!             olapm(i,i) = 1/eval(i)
!           END DO
!
!           ALLOCATE( invevec(ngptm(ikpt),ngptm(ikpt)) )
! #ifdef CPP_INVERSION
!           invevec =         transpose(evec)
! #else
!           invevec = conjg( transpose(evec) )
! #endif
!
!           olapm   = matmul(evec,matmul(olapm,invevec) )
!
!           DEALLOCATE(invevec)!,involapm)
!         ELSE
          !calculate inverse overlap-matrix
        CALL inverse(olapm)
!         END IF

        !unpack matrix coulomb
        ALLOCATE( coulhlp(nbasm1(ikpt),nbasm1(ikpt)) )
        coulhlp=unpackmat(coulomb(:nbasm1(ikpt)*(nbasm1(ikpt)+1)/2,ikpt))

        !multiply with inverse olap from right hand side
        coulhlp(:,nbasp+1:) = matmul(coulhlp(:,nbasp+1:),olapm)

        !multiply with inverse olap from left side
        coulhlp(nbasp+1:,:) = matmul(olapm,coulhlp(nbasp+1:,:))

        coulomb(:nbasm1(ikpt)*(nbasm1(ikpt)+1)/2,ikpt)&
     &    = packmatcoul(coulhlp)
        DEALLOCATE ( olapm,coulhlp )!,eval,evec )

      END DO

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(1X,A)',advance='no') 'done'
        CALL cpu_time(time2) 
        WRITE(6,'(2X,A,F8.2,A)') '( Timing:',time2-time1,' )'
      END IF

      !
      ! plot matrix
      !
      IF( l_plot ) THEN
        IF ( mpi%isize /= 1 )&
     &    STOP 'coulombmatrix: l_plot only works with one process'
        DO ikpt = 1,kpts%nkpt

          ic = 0
          DO i = 1,nbasm1(ikpt)
            DO j = 1,i
              ic = ic + 1
              IF( abs(coulomb(ic,ikpt)) .gt. 1E-8 ) THEN 
                WRITE(600+ikpt,'(2i6)') i,j
                WRITE(600+ikpt,'(2i6)') j,i
              END IF
            END DO
          END DO

          ALLOCATE( coulhlp(nbasm1(ikpt),nbasm1(ikpt)) )      
          coulhlp = unpackmat( &
     &                    coulomb(:nbasm1(ikpt)*(nbasm1(ikpt)+1)/2,ikpt) )

          ic = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  WRITE(700+ikpt,*) l,M
                  DO n = 1,hybrid%nindxm1(l,itype)
                    WRITE(700+ikpt,'(16f8.4)')&
     &                             coulhlp(ic+n,ic+1:ic+hybrid%nindxm1(l,itype))
                  END DO
                  ic = ic + hybrid%nindxm1(l,itype)
                ENDDO
              END DO
            END DO
          END DO
          ALLOCATE( coulhlp1(nbasm1(ikpt),nbasm1(ikpt)) )
          coulhlp1 = 0

          ic2 = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  DO n = 1,hybrid%nindxm1(l,itype) - 1
                    ic2 = ic2 + 1
                  END DO
                END DO
              END DO
            END DO
          END DO

          ic  = 0
          ic1 = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  DO n = 1,hybrid%nindxm1(l,itype) - 1
                    coulhlp1(ic+n,ic+1:ic+hybrid%nindxm1(l,itype)-1)&
     &                      = coulhlp(ic1+n,ic1+1:ic1+hybrid%nindxm1(l,itype)-1)
                  END DO

                  coulhlp1(ic2+1,ic+1:ic+hybrid%nindxm1(l,itype)-1)&
     &        = coulhlp(ic1+hybrid%nindxm1(l,itype),ic1+1:ic1+hybrid%nindxm1(l,itype)-1)

                  ic  = ic  + hybrid%nindxm1(l,itype)-1
                  ic1 = ic1 + hybrid%nindxm1(l,itype)
                  ic2 = ic2 + 1

                END DO
              END DO
            END DO
          END DO

          IF( ikpt .eq. 1 ) THEN
            ! add contributions from s channel and G=0 component

            ic = 0; ic1 = 0
            DO itype = 1,atoms%ntype
              DO ineq = 1,atoms%neq(itype)
                WRITE(*,*) ic+1,ic+hybrid%nindxm1(0,itype)-1,ic1+1,&
     &                     ic1+hybrid%nindxm1(0,itype)-1
                coulhlp1(ic+1:ic+hybrid%nindxm1(0,itype)-1,nbasp+1)&
     &            = coulhlp(ic1+1:ic1+hybrid%nindxm1(0,itype)-1,nbasp+1)
                coulhlp1(nbasp+1,ic+1:ic+hybrid%nindxm1(0,itype)-1)&
     &            = coulhlp(nbasp+1,ic1+1:ic1+hybrid%nindxm1(0,itype)-1)
                ic = ic  + sum( (/ ( (2*l+1)*(hybrid%nindxm1(l,itype)-1),&
     &                          l=0,hybrid%lcutm1(itype) ) /) )
                ic1= ic1 + sum( (/ ( (2*l+1)*hybrid%nindxm1(l,itype),&
     &                          l=0,hybrid%lcutm1(itype) ) /) )
              END DO
            END DO

            ic2 = 0
            DO itype = 1,atoms%ntype
              DO ineq = 1,atoms%neq(itype)
                DO l = 0,hybrid%lcutm1(itype)
                  DO M = -l,l
                    DO n = 1,hybrid%nindxm1(l,itype) - 1
                      ic2 = ic2 + 1
                    END DO
                  END DO
                END DO
              END DO
            END DO

            ic1 = 0
            DO itype = 1,atoms%ntype
              DO ineq = 1,atoms%neq(itype)
                DO l = 0,hybrid%lcutm1(itype)
                  DO M = -l,l
                    ic2 = ic2 + 1

                    ic1 = ic1 + hybrid%nindxm1(l,itype)

                    IF( l .ne. 0 ) CYCLE
                    WRITE(900,*) ic2,ic1,itype,ineq

                    ic3 = 0
                    ic4 = 0
                    DO itype1 = 1,atoms%ntype
                      ishift = sum( (/ ( (2*l2+1)* hybrid%nindxm1(l2,itype1)   ,&
     &                              l2 = 0,hybrid%lcutm1(itype1) ) /) )
                      ishift1= sum( (/ ( (2*l2+1)*(hybrid%nindxm1(l2,itype1)-1),&
     &                              l2 = 0,hybrid%lcutm1(itype1) ) /) )
                      DO ineq1 = 1,atoms%neq(itype1)
                        ic5 = ic3 + (ineq1 - 1)*ishift + 1
                        ic6 = ic5 + hybrid%nindxm1(0,itype1) - 2

                        ic7 = ic4 + (ineq1 - 1)*ishift1 + 1
                        ic8 = ic7 + hybrid%nindxm1(0,itype1) - 2
                        WRITE(901,*)ic2,ic7,ic8,ic1,ic5,ic6,itype,itype1


                        coulhlp1(ic2,ic7:ic8) = coulhlp(ic1,ic5:ic6)
#ifdef CPP_INVERSION  
                        coulhlp1(ic7:ic8,ic2) = coulhlp1(ic2,ic7:ic8)
#else
                        coulhlp1(ic7:ic8,ic2) =&
     &                                   conjg(coulhlp1(ic2,ic7:ic8))
#endif
                      END DO
                      ic3 = ic3 + ishift *atoms%neq(itype1)
                      ic4 = ic4 + ishift1*atoms%neq(itype1)
                    END DO

                  END DO
                END DO
              END DO
            END DO

          END IF

          ic2 = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  DO n = 1,hybrid%nindxm1(l,itype) - 1
                    ic2 = ic2 + 1
                  END DO
                END DO
              END DO
            END DO
          END DO

          ic1 = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l
                  ic2 = ic2 + 1

                  ic1 = ic1 + hybrid%nindxm1(l,itype)

                  ic3 = 0
                  ic4 = ic2
                  DO itype1 = 1,atoms%ntype
                    DO ineq1 = 1,atoms%neq(itype1)
                      DO l1 = 0,hybrid%lcutm1(itype1)
                        DO m1 = -l1,l1
                          ic3 = ic3 + hybrid%nindxm1(l1,itype1)
                          IF( ic3 .lt. ic1 ) CYCLE
                          WRITE(300,'(4i6,2f15.10)')&
     &                                  ic2,ic4,ic1,ic3,coulhlp(ic1,ic3)
                          coulhlp1(ic2,ic4) = coulhlp(ic1,ic3)
#ifdef CPP_INVERSION  
                          coulhlp1(ic4,ic2) = coulhlp1(ic2,ic4)
#else
                          coulhlp1(ic4,ic2) = conjg(coulhlp1(ic2,ic4))
#endif
                          ic4 = ic4 + 1
                        END DO
                      END DO
                    END DO
                  END DO

                  DO igpt = 1,hybrid%ngptm(ikpt)
                    coulhlp1(ic2,ic4) = coulhlp(ic1,nbasp+igpt)
#ifdef CPP_INVERSION
                    coulhlp1(ic4,ic2) = coulhlp1(ic2,ic4)
#else
                    coulhlp1(ic4,ic2) = conjg(coulhlp1(ic2,ic4))
#endif
                    ic4 = ic4 + 1
                  END DO

                END DO
              END DO
            END DO
          END DO

          coulhlp1(nbasp+1:,nbasp+1:) = coulhlp(nbasp+1:,nbasp+1:)

          ic = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)

                DO M = -l,l
                  WRITE(800+ikpt,*) l,M
                  DO n = 1,hybrid%nindxm1(l,itype)-1
                    WRITE(800+ikpt,'(16f8.4)')&
     &                          coulhlp1(ic+n,ic+1:ic+hybrid%nindxm1(l,itype)-1)
                  END DO
                  ic = ic + hybrid%nindxm1(l,itype)-1
                ENDDO
              END DO
            END DO
          END DO

          ic = 0
          DO i = 1,nbasm1(ikpt)
            DO j = 1,i
              IF( abs(coulhlp1(i,j)) .gt. 1E-8 ) THEN 
                WRITE(850+ikpt,'(2i6)') i,j
                WRITE(850+ikpt,'(2i6)') j,i
              END IF
            END DO
          END DO
          DEALLOCATE( coulhlp,coulhlp1 )
        END DO
        STOP
      END IF ! lplot

      IF ( mpi%irank == 0 )&
     &  WRITE(6,'(A)',advance='no') 'Writing of data to file...'
      CALL cpu_time(time1)
#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
      !
      ! rearrange coulomb matrix
      !

#     ifndef CPP_MPI
        OPEN(unit=778,file=fname,form='unformatted',access='direct',&
     &       recl=irecl_coulomb)
        ALLOCATE( coulomb_mt1(hybrid%maxindxm1-1,hybrid%maxindxm1-1,0:hybrid%maxlcutm1,atoms%ntype,1))
        ALLOCATE( coulomb_mt2(hybrid%maxindxm1-1,-hybrid%maxlcutm1:hybrid%maxlcutm1,&
     &                        0:hybrid%maxlcutm1+1,atoms%nat,1) )
        ALLOCATE( coulomb_mt3(hybrid%maxindxm1-1,atoms%nat,atoms%nat,1) )
#     else
        ALLOCATE( coulomb_mt1(hybrid%maxindxm1-1,hybrid%maxindxm1-1,0:hybrid%maxlcutm1,atoms%ntype,&
     &                        ikptmin:ikptmax) )
        ALLOCATE( coulomb_mt2(hybrid%maxindxm1-1,-hybrid%maxlcutm1:hybrid%maxlcutm1,&
     &                        0:hybrid%maxlcutm1+1,atoms%nat,ikptmin:ikptmax) )
        ALLOCATE( coulomb_mt3(hybrid%maxindxm1-1,atoms%nat,nat,ikptmin:ikptmax) )
#     endif
      ic = (hybrid%maxlcutm1+1)**2*atoms%nat
#     ifdef CPP_IRCOULOMBAPPROX
#       ifndef CPP_MPI
          ALLOCATE( coulomb_mtir(ic,ic+hybrid%maxgptm,1) )
#       else
          ALLOCATE( coulomb_mtir(ic,ic+hybrid%maxgptm,kpts%nkpt )
#       endif
        ALLOCATE( coulombp_mtir(0,0) )
#     else
        ALLOCATE( coulomb_mtir(ic+hybrid%maxgptm,ic+hybrid%maxgptm,1) )
        idum = ic + hybrid%maxgptm
        idum = ( idum*(idum+1) ) / 2
#       ifndef CPP_MPI
          ALLOCATE( coulombp_mtir(idum,1) )
#       else
          ALLOCATE( coulombp_mtir(idum,ikptmin:ikptmax) )
#       endif
#     endif

#     ifdef CPP_MPI
        ! initialize arrays
        coulomb_mt1 = 0 ; coulomb_mt2 = 0 ; coulomb_mt3 = 0
#       ifdef CPP_IRCOULOMBAPPROX
          coulomb_mtir  = 0
#       else
          coulombp_mtir = 0
#       endif

        ! Allocate memory to store all request handles
        ALLOCATE( req( (sum(nsym_gpt)+4)*kpts%nkpt ) )
        reqd = 0
#     endif

      DO ikpt = ikptmin,ikptmax
#       ifdef CPP_MPI
          ikpt0 = ikpt
#         ifdef CPP_IRCOULOMBAPPROX
            ikpt1 = ikpt
#         else
            ikpt1 = 1
            coulomb_mtir = 0
#         endif
#       else
          ikpt0 = 1
          ikpt1 = 1
          ! initialize arrays
          coulomb_mt1 = 0 ; coulomb_mt2   = 0 
          coulomb_mt3 = 0 ; coulombp_mtir = 0
#       endif

        ! unpack coulomb into coulhlp
        ALLOCATE( coulhlp(nbasm1(ikpt),nbasm1(ikpt)) )
        coulhlp=unpackmat( coulomb(:nbasm1(ikpt)*(nbasm1(ikpt)+1)/2,ikpt))

        ! only one processor per k-point calculates MT convolution
        IF ( calc_mt(ikpt) ) THEN

          !
          ! store m-independent part of Coulomb matrix in MT spheres
          ! in coulomb_mt1(:hybrid%nindxm1(l,itype)-1,:hybrid%nindxm1(l,itype)-1,l,itype)
          !
          indx1 = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              DO l = 0,hybrid%lcutm1(itype)

                IF ( ineq == 1 ) THEN
                  DO n = 1,hybrid%nindxm1(l,itype) - 1
                    coulomb_mt1(n,1:hybrid%nindxm1(l,itype)-1,l,itype,ikpt0)&
     &              = coulhlp(indx1+n,indx1+1:indx1+hybrid%nindxm1(l,itype)-1)
                  END DO
                END IF

                indx1 = indx1 + (2*l+1)*hybrid%nindxm1(l,itype)
              END DO
            END DO
          END DO

          !
          ! store m-dependent and atom-dependent part of Coulomb matrix in MT spheres
          ! in coulomb_mt2(:hybrid%nindxm1(l,itype)-1,-l:l,l,iatom)
          !
          indx1 = 0
          iatom = 0
          DO itype = 1,atoms%ntype
            DO ineq = 1,atoms%neq(itype)
              iatom = iatom + 1
              DO l = 0,hybrid%lcutm1(itype)
                DO M = -l,l

                  coulomb_mt2(:hybrid%nindxm1(l,itype)-1,M,l,iatom,ikpt0)&
     & = coulhlp(indx1+1:indx1+hybrid%nindxm1(l,itype)-1,indx1+hybrid%nindxm1(l,itype))

                  indx1 = indx1 + hybrid%nindxm1(l,itype)

                END DO
              END DO
            END DO
          END DO

          !
          ! due to the subtraction of the divergent part at the Gamma point
          ! additional contributions occur
          !
          IF ( ikpt .eq. 1 ) THEN
            !
            ! store the contribution of the G=0 plane wave with the MT l=0 functions in
            ! coulomb_mt2(:hybrid%nindxm1(l=0,itype),0,hybrid%maxlcutm1+1,iatom)
            !
            ic    = 0
            iatom = 0
            DO itype = 1,atoms%ntype
              DO ineq = 1,atoms%neq(itype)
                iatom = iatom + 1
                DO n = 1,hybrid%nindxm1(0,itype)-1
                  coulomb_mt2(n,0,hybrid%maxlcutm1+1,iatom,ikpt0)&
     &            = coulhlp(ic+n,nbasp+1)
                END DO
                ic = ic + sum( (/ ( (2*l+1)*hybrid%nindxm1(l,itype),l=0,&
     &                         hybrid%lcutm1(itype) ) /) )
              END DO
            END DO

            !
            ! store the contributions between the MT s-like functions at atom1 and
            ! and the constant function at a different atom2
            !
            iatom = 0
            ic    = 0
            DO itype = 1,atoms%ntype
              ishift = sum( (/ ( (2*l+1)*hybrid%nindxm1(l,itype),l = 0,&
     &                      hybrid%lcutm1(itype) ) /) )
              DO ineq = 1,atoms%neq(itype)
                iatom = iatom + 1
                ic1   = ic    + hybrid%nindxm1(0,itype)

                iatom1 = 0
                ic2    = 0
                DO itype1 = 1,atoms%ntype
                  ishift1 = sum( (/ ( (2*l1+1)*hybrid%nindxm1(l1,itype1),l1 = 0,&
     &                           hybrid%lcutm1(itype1) ) /) )
                  DO ineq1 = 1,atoms%neq(itype1)
                    iatom1 = iatom1 + 1
                    ic3 = ic2 +  1
                    ic4 = ic3 + hybrid%nindxm1(0,itype1)   - 2
#ifdef CPP_INVERSION
                    coulomb_mt3(:hybrid%nindxm1(0,itype1)-1,iatom,iatom1,ikpt0)&
     &                = coulhlp(ic1,ic3:ic4)
#else
                    coulomb_mt3(:hybrid%nindxm1(0,itype1)-1,iatom,iatom1,ikpt0)&
     &                = conjg(coulhlp(ic1,ic3:ic4))
#endif
                    ic2 = ic2 + ishift1
                  END DO
                END DO

                ic = ic + ishift
              END DO
            END DO

            !test
            iatom = 0
            DO itype = 1,atoms%ntype
              DO ineq = 1,atoms%neq(itype)
                iatom = iatom + 1
                IF( maxval( abs(coulomb_mt2(:hybrid%nindxm1(0,itype)-1,0,0,&
     &                                      iatom,ikpt0)&
     &                         -coulomb_mt3(:hybrid%nindxm1(0,itype)-1,iatom,&
     &                                      iatom,ikpt0)) )&
     &             .gt. 1E-08 ) THEN
                  STOP 'coulombmatrix: &&
     &                  coulomb_mt2 and coulomb_mt3 are inconsistent'
                END IF
              END DO
            END DO
          END IF

        END IF ! calc_mt

        !
        ! add the residual MT contributions, i.e. those functions with an moment,
        ! to the matrix coulomb_mtir, which is fully occupied
        !
        ic = 0
        DO itype = 1,atoms%ntype
          DO ineq = 1,atoms%neq(itype)
            DO l = 0,hybrid%lcutm1(itype)
              DO M = -l,l
                ic = ic + 1
              END DO
            END DO
          END DO
        END DO

        indx1 = 0; indx2 = 0; indx3 = 0; indx4 = 0
        DO itype = 1,atoms%ntype
          DO ineq = 1,atoms%neq(itype)
            DO l = 0,hybrid%lcutm1(itype)
              DO M = -l,l
                indx1 = indx1 + 1
                indx3 = indx3 + hybrid%nindxm1(l,itype)

                indx2 = 0
                indx4 = 0
                DO itype1 = 1,atoms%ntype
                  DO ineq1 = 1,atoms%neq(itype1)
                    DO l1 = 0,hybrid%lcutm1(itype1)
                      DO m1 = -l1,l1
                        indx2 = indx2 + 1
                        indx4 = indx4 + hybrid%nindxm1(l1,itype1)
                        IF( indx4 .lt. indx3 ) CYCLE
                        IF ( calc_mt(ikpt) ) THEN
                          coulomb_mtir(indx1,indx2,ikpt1)&
     &                    = coulhlp(indx3,indx4)
#ifdef CPP_INVERSION
                          coulomb_mtir(indx2,indx1,ikpt1)&
     &                    = coulomb_mtir(indx1,indx2,ikpt1)
#else
                          coulomb_mtir(indx2,indx1,ikpt1)&
     &                    = conjg(coulomb_mtir(indx1,indx2,ikpt1))
#endif
                        END IF
                      END DO
                    END DO
                  END DO
                END DO

                DO igpt = 1,hybrid%ngptm(ikpt)
                  indx2 = indx2 + 1
                  coulomb_mtir(indx1,indx2,ikpt1)&
     &            = coulhlp(indx3,nbasp+igpt)
#if !defined CPP_IRCOULOMBAPPROX

#ifdef CPP_INVERSION
                  coulomb_mtir(indx2,indx1,ikpt1)&
     &            = coulomb_mtir(indx1,indx2,ikpt1)
#else
                  coulomb_mtir(indx2,indx1,ikpt1)&
     &            = conjg(coulomb_mtir(indx1,indx2,ikpt1))
#endif

#endif
                END DO

              END DO
            END DO
          END DO
        END DO

        IF( indx1 .ne. ic ) STOP 'coulombmatrix: error index counting'

        !
        ! add ir part to the matrix coulomb_mtir
        !
#       ifndef CPP_IRCOULOMBAPPROX
          coulomb_mtir(ic+1:ic+hybrid%ngptm(ikpt),ic+1:ic+hybrid%ngptm(ikpt),ikpt1)&
     &    = coulhlp(nbasp+1:nbasm1(ikpt),nbasp+1:nbasm1(ikpt))
          ic2 = indx1+hybrid%ngptm(ikpt)
          coulombp_mtir(:ic2*(ic2+1)/2,ikpt0)&
     &    = packmat(coulomb_mtir(:ic2,:ic2,ikpt1))
#       endif

#       ifndef CPP_MPI

#         ifdef CPP_IRCOULOMBAPPROX
            WRITE(778,rec=ikpt) coulomb_mt1,coulomb_mt2,coulomb_mt3,&
     &                          coulomb_mtir
#         else
            WRITE(778,rec=ikpt) coulomb_mt1,coulomb_mt2,coulomb_mt3,&
     &                          coulombp_mtir
#         endif

#       else

          ! determine the offset of this k-point
          offset0 = irecl_coulomb * (ikpt - 1)
          ! write MT-MT part
          IF ( calc_mt(ikpt) ) THEN
            ! write different parts of coulomb matrix and change offset accordingly
            icnt = size( coulomb_mt1 ) / nkminmax
            reqd = reqd + 1
            CALL MPI_FILE_IWRITE_AT(fh,offset0,&
     &         coulomb_mt1(:,:,:,:,ikpt0),icnt,MPI_REAL8,req(reqd),ierr)
            ! coulomb_mt1 is always real, so use 8 bytes
            offset0 = offset0 + icnt * 8
            icnt = size( coulomb_mt2 ) / nkminmax
            reqd = reqd + 1
            CALL MPI_FILE_IWRITE_AT(fh,offset0,&
     &          coulomb_mt2(:,:,:,:,ikpt0),icnt,datatype,req(reqd),ierr)
            offset0 = offset0 + icnt * coul_size
            icnt = size( coulomb_mt3 ) / nkminmax
            reqd = reqd + 1
            CALL MPI_FILE_IWRITE_AT(fh,offset0,coulomb_mt3(:,:,:,ikpt0),&
     &                                     icnt,datatype,req(reqd),ierr)
            offset0 = offset0 + icnt * coul_size
            reqd = reqd + 1
#           ifdef CPP_IRCOULOMBAPPROX
              icnt = ic * ic
              CALL MPI_FILE_IWRITE_AT(fh,offset0,&
     &             coulomb_mtir(:,:,ikpt1),icnt,datatype,req(reqd),ierr)
#           else
              icnt = ic * (ic+1) / 2
              CALL MPI_FILE_IWRITE_AT(fh,offset0,coulombp_mtir(:,ikpt0),&
     &                                     icnt,datatype,req(reqd),ierr)
#           endif
          ELSE
            offset0 = offset0 + ( 8 * size(coulomb_mt1) + coul_size *&
     &            ( size(coulomb_mt2) + size(coulomb_mt3) ) ) / nkminmax
          END IF

          ! loop over all g-points handled by this process
          DO igpt0 = igptmin(ikpt),igptmax(ikpt)

              ! loop over g-points which are symmetry equivalent
              DO isym = 1,nsym_gpt(igpt0,ikpt)

                ! Get corresponding index of the field
                igpt = sym_gpt(isym,igpt0,ikpt)

#               ifdef CPP_IRCOULOMBAPPROX
                  idisp  = ic + igpt
                  offset = offset0 + (idisp-1) * ic * coul_size
                  ! write coulomb matrix and store the request handle
                  reqd   = reqd + 1
                  CALL MPI_FILE_IWRITE_AT(fh,offset,&
     &                                    coulomb_mtir(:,idisp,ikpt1),&
     &                                    ic,datatype,req(reqd),ierr)
#               else
                  icnt   = ic + igpt
                  idisp  = (icnt-1) * icnt / 2
                  offset = offset0 + idisp * coul_size
                  ! write coulomb matrix and store the request handle
                  reqd   = reqd + 1
                  CALL MPI_FILE_IWRITE_AT(fh,offset,&
     &                          coulombp_mtir(idisp+1:idisp+icnt,ikpt0),&
     &                          icnt,datatype,req(reqd),ierr)
#               endif

              END DO

          END DO ! igpt0

#       endif

        DEALLOCATE( coulhlp )

      END DO ! ikpt

#     ifndef CPP_MPI

        CLOSE(778)

#     else

        IF ( mpi%irank == 0 .AND. hybrid%ngptm(kpts%nkpt) /= hybrid%maxgptm ) THEN
          offset = irecl_coulomb * kpts%nkpt - 1
          reqd = reqd + 1
          CALL MPI_FILE_IWRITE_AT(fh,offset,0,1,MPI_BYTE,req(reqd),ierr)
        END IF
        CALL MPI_WAITALL(reqd,req,MPI_STATUSES_IGNORE,ierr)
        ! Check if an error occured
        IF ( ierr /= 0 ) THEN
          CALL MPI_ERROR_STRING(ierr,errmsg,length,ierr)
          WRITE(*,*) errmsg(1:length)
          STOP 'coulombmatrix failed with mpi error'
        END IF
        DEALLOCATE( req )

        CALL MPI_GET_ADDRESS(coulomb_mt1, iaddr,ierr)
        CALL MPI_GET_ADDRESS(coulomb_mt2, iaddr,ierr)
        CALL MPI_GET_ADDRESS(coulomb_mt3, iaddr,ierr)
        CALL MPI_GET_ADDRESS(coulomb_mtir,iaddr,ierr)

#     endif

      DEALLOCATE ( coulomb_mt1,  coulomb_mt2, coulomb_mt3,&
     &             coulomb_mtir, coulombp_mtir )

#else

      ! write the coulomb matrix parallel with all processes
#     ifdef CPP_MPI

        ! Allocate memory to store all request handles
        ALLOCATE( req( sum(nsym_gpt)+kpts%nkpt ) )
        reqd = 0
        ! loop over all k-points handled on this process
        DO ikpt = ikptmin,ikptmax

          ! determine the offset of this k-point
          offset0 = irecl_coulomb * (ikpt - 1)

          ! write MT-MT part
          IF ( calc_mt(ikpt) ) THEN
            icnt = nbasp*(nbasp+1)/2
            ! write coulomb matrix and store the request handle
            reqd = reqd + 1
            CALL MPI_FILE_IWRITE_AT(fh,offset0,coulomb(1:icnt,ikpt),&
     &                              icnt,datatype,req(reqd),ierr)
          END IF

          ! loop over all g-points
          DO igpt0 = igptmin(ikpt),igptmax(ikpt)

            ! loop over g-points which are symmetry equivalent
            DO isym = 1,nsym_gpt(igpt0,ikpt)

              ! Determine position of first element
              igpt   = sym_gpt(isym,igpt0,ikpt)
              icnt   = nbasp + igpt
              idisp  = (icnt-1) * icnt / 2
              offset = offset0 + idisp * coul_size

              ! write coulomb matrix and store the request handle
              reqd = reqd + 1
              CALL MPI_FILE_IWRITE_AT(fh,offset,&
     &                                coulomb(idisp+1:idisp+icnt,ikpt),&
     &                                icnt,datatype,req(reqd),ierr)

            END DO

          END DO ! igpt0

        END DO ! ikpt

        ! wait until file processing is done
        CALL MPI_WAITALL(reqd,req,MPI_STATUSES_IGNORE,ierr)
        ! Check if no error occured
        IF ( ierr /= 0 ) THEN
          CALL MPI_ERROR_STRING(ierr,errmsg,length,ierr)
          WRITE(*,*) errmsg(1:length)
          STOP 'coulombmatrix failed with mpi error'
        END IF
        DEALLOCATE( req )

#     else

        !write coulomb matrix to direct access file coulomb
        OPEN(unit=777,file=fname,form='unformatted',access='direct',&
     &       recl=irecl_coulomb)
        DO i=1,kpts%nkpt
          WRITE(777,rec=i) coulomb(:,i)
        END DO
        CLOSE(777)

#     endif

#endif

      ! close the file after the data was written
#     ifdef CPP_MPI
        CALL MPI_FILE_CLOSE(fh,ierr)
#     endif

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(7X,A)',advance='no') 'done'
        CALL cpu_time(time2) 
        WRITE(6,'(2X,A,F8.2,A)') '( Timing:',time2-time1,' )'
      END IF

      CONTAINS

!     Calculate body of Coulomb matrix at Gamma point: v_IJ = SUM(G) c^*_IG c_JG 4*pi/G**2 .
!     For this we must subtract from coulomb(:,1) the spherical average of a term that comes
!     from the fact that MT functions have k-dependent Fourier coefficients (see script).
      SUBROUTINE subtract_sphaverage()

      USE m_wrapper
      USE m_trafo
      USE m_util
      IMPLICIT NONE

      ! - local scalars -
      INTEGER               :: l ,i,j,n,nn,itype,ieq

      ! - local arrays -
#ifdef CPP_INVERSION
      REAL                  :: olap(hybrid%ngptm(1),hybrid%ngptm(1))
      REAL,    ALLOCATABLE  :: constfunc(:)
#else
      COMPLEX               :: olap(hybrid%ngptm(1),hybrid%ngptm(1))
      COMPLEX , ALLOCATABLE :: constfunc(:)
#endif
      COMPLEX               :: coeff(nbasm1(1)),cderiv(nbasm1(1),-1:1),&
     &                         claplace(nbasm1(1))

      n  = nbasm1(1)
      nn = n*(n+1)/2

      CALL olap_pw ( olap,hybrid%pgptm(:hybrid%ngptm(1),1),hybrid%ngptm(1),atoms,cell )

      ! Define coefficients (coeff) and their derivatives (cderiv,claplace)
      coeff    = 0
      cderiv   = 0
      claplace = 0
      j        = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          DO l = 0,hybrid%lcutm1(itype)
            DO M = -l,l
              DO i = 1,hybrid%nindxm1(l,itype)
                j = j + 1
                IF(l.eq.0) THEN
                  coeff(j)    =  sqrt(4*pi_const)&
     &            * intgrf( atoms%rmsh(:,itype)* hybrid%basm1(:,i,0,itype),&
     &                      atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
     &            / sqrt(cell%vol)

                  claplace(j) = -sqrt(4*pi_const)&
     &            * intgrf(atoms%rmsh(:,itype)**3*hybrid%basm1(:,i,0,itype),&
     &                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
     &            / sqrt(cell%vol)

                ELSE IF(l.eq.1) THEN
                  cderiv(j,M) = -sqrt(4*pi_const/3)*img&
     &            * intgrf(atoms%rmsh(:,itype)**2*hybrid%basm1(:,i,1,itype),&
     &                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
     &            / sqrt(cell%vol)
                END IF
              END DO
            END DO
          END DO
        END DO
      END DO
      coeff(nbasp+1:n) = olap(1,1:n-nbasp)
# ifdef CPP_INVERSION
      CALL symmetrize(coeff,       1,nbasm1(1),2,.false.,&
     &                  atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                  hybrid%nindxm1,sym)
      CALL symmetrize(claplace,    1,nbasm1(1),2,.false.,&
     &                atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                hybrid%nindxm1,sym)
      CALL symmetrize(cderiv(:,-1),1,nbasm1(1),2,.false.,&
     &                atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                hybrid%nindxm1,sym)
      CALL symmetrize(cderiv(:, 0),1,nbasm1(1),2,.false.,&
     &                atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                hybrid%nindxm1,sym)
      CALL symmetrize(cderiv(:, 1),1,nbasm1(1),2,.false.,&
     &                atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                hybrid%nindxm1,sym)
# endif
      ALLOCATE ( constfunc(n) )
      ! Subtract head contributions from coulomb(:nn,1) to obtain the body
      l = 0
      DO j = 1,n
        DO i = 1,j
          l            = l + 1
          coulomb(l,1) = coulomb(l,1) - 4*pi_const/3&
     &                 * ( dot_product(cderiv(i,:),cderiv(j,:))&
     &                 + ( conjg(coeff(i)) * claplace(j)&
     &                   + conjg(claplace(i)) * coeff(j) ) / 2 )
        END DO
      END DO
      coeff(nbasp+1)  = 1d0
      coeff(nbasp+2:) = 0d0
# ifdef CPP_INVERSION
      CALL desymmetrize(coeff,1,nbasm1(1),2,&
     &                  atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                  hybrid%nindxm1,sym)
      CALL   symmetrize(coeff,nbasm1(1),1,1,.false.,&
     &                  atoms,ntype,hybrid%lcutm1,hybrid%maxlcutm1,&
     &                  hybrid%nindxm1,sym)
# endif
      ! Explicit normalization here in order to prevent failure of the diagonalization in diagonalize_coulomb
      ! due to inaccuracies in the overlap matrix (which can make it singular).
      constfunc = coeff / sqrt ( ( sum(abs(coeff(:nbasp))**2)&
     &   + dotprod ( coeff(nbasp+1:), matmul(olap,coeff(nbasp+1:)) ) ) )

      END SUBROUTINE subtract_sphaverage


      END SUBROUTINE coulombmatrix

!     -----------------------------------------------------------------------------------------------


!     Calculates the structure constant
!                                                        1               *      ^
!     structconst(lm,ic1,ic2,k) = SUM exp(ikT) -----------------------  Y  ( T + R(ic) )
!                                  T           | T + R(ic1) - R(ic2) |   lm
!
!     with T = lattice vectors
!
!     An Ewald summation method devised by O.K. Andersen is used for l<5 
!     (see e.g. H.L. Skriver, "The LMTO method", Springer 1984).
!     (The real-space function G can be calculated with gfunction.f)
!
 
! Convergence parameter
#define CONVPARAM 1d-18
! Do some additional shells ( real-space and Fourier-space sum )
#define ADDSHELL1 40
#define ADDSHELL2 0

#define HLP9 (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/9*

      SUBROUTINE structureconstant(&
     &               structconst,cell,hybrid,atoms,&
     &               kpts,mpi )

      USE m_constants, ONLY: pi_const
      USE m_util     , ONLY: harmonicsr,rorderp,rorderpf
    USE m_types
      IMPLICIT NONE

      TYPE(t_mpi),INTENT(IN)   :: mpi

      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid

      TYPE(t_cell),INTENT(IN)   :: cell

      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_kpts),INTENT(IN)    :: kpts
      ! - scalars -

      ! - arrays -
      COMPLEX , INTENT(INOUT)   ::  structconst((2*hybrid%lexp+1)**2 ,atoms%nat,atoms%nat, kpts%nkpt)

      ! - local scalars -
      INTEGER                   ::  i,ic1,ic2,lm,ikpt,l ,ilastsh, ishell,nshell,i1,i2
      INTEGER                   ::  ix,iy,iz,n,nn,m
      INTEGER                   ::  nptsh,maxl,maxlm
      
      REAL                      ::  rad,rrad ,rdum,rvol
      REAL                      ::  a,a1,aa
      REAL                      ::  pref,rexp
      REAL                      ::  time1,time2
      REAL                      ::  scale

      COMPLEX                   ::  cdum,cexp
      COMPLEX , PARAMETER       ::  img = (0d0,1d0)

      LOGICAL,    SAVE          ::  first = .true.
      ! - local arrays -
      INTEGER                   ::  conv(0:2*hybrid%lexp) 
      INTEGER , ALLOCATABLE     ::  shell(:)
      INTEGER , ALLOCATABLE     ::  pnt(:),ptsh(:,:)

      REAL                      ::  r(3),rc(3),ra(3),k(3),ki(3),ka(3)
      REAL                      ::  convpar(0:2*hybrid%lexp),g(0:2*hybrid%lexp)
      REAL    , ALLOCATABLE     ::  radsh(:)
      REAL    , ALLOCATABLE     ::  structdum2(:,:)

      COMPLEX                   ::  y((2*hybrid%lexp+1)**2)
      COMPLEX                   ::  shlp((2*hybrid%lexp+1)**2,kpts%nkpt)
      COMPLEX , ALLOCATABLE     ::  structdum(:,:,:),structhlp(:,:), chlp(:),structdum1(:)

      IF ( mpi%irank /= 0 ) first = .false.

      rdum   = cell%vol**(1d0/3) ! define "average lattice parameter"

      ! ewaldlambda = ewaldscale
      scale  = hybrid%ewaldlambda  / rdum

!       lambda = ewaldlambda / rdum

      pref = 4*pi_const / (scale**3 * cell%vol)

      DO l = 0,2*hybrid%lexp
        convpar(l) = CONVPARAM / scale**(l+1)
      END DO

      IF( first ) THEN
        WRITE(6,'(//A)') '### subroutine: structureconstant ###'
        WRITE(6,'(/A)')  'Real-space sum:'
      END IF

!
!     Determine cutoff radii for real-space and Fourier-space summation
      ! (1) real space
      a    = 1
 1    rexp = exp(-a)
      g(0) = rexp / a    * (1+a*11/16*(1+a*3/11*(1+a/9)))
      g(1) = rexp / a**2 * (1+a*(1+a/2*(1+a*7/24*(1+a/7))))
      g(2) = rexp / a**3 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a*3/16&
     &*(1+a/9))))))
      g(3) = rexp / a**4 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6&
     &*(1+a/8)))))))
      g(4) = rexp / a**5 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6&
     &*(1+a/7*(1+a/8*(1+a/10)))))))))
      g(5) = rexp / a**6 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6&
     &*(1+a/7*(1+a/8*(1+a/9*(1+a/10))))))))))
      g(6) = rexp / a**7 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6&
     &*(1+a/7*(1+a/8*(1+a/9*(1+a/10*(1+a/11*(1+a/12))))))))))))
      g(7) = rexp / a**8 * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6&
     &*(1+a/7*(1+a/8*(1+a/9*(1+a/10*(1+a/11*(1+a/12*(1+a/13)))))))))))))
      DO l = 8,2*hybrid%lexp
        g(l) = a**(-l-1)
      END DO
      IF(any(g.gt.convpar/10)) THEN ! one digit more accuracy for real-space sum
        a = a + 1
        GOTO 1
      END IF
      rad = a / scale

      ! (2) Fourier space
      a    = 1
 2    aa   = (1+a**2)**(-1)
      g(0) = pref * aa**4 / a**2
      g(1) = pref * aa**4 / a
      g(2) = pref * aa**5        / 3
      g(3) = pref * aa**5 * a    / 15
      g(4) = pref * aa**6 * a**2 / 105
      g(5) = pref * aa**6 * a**3 / 945
      g(6) = pref * aa**7 * a**4 / 10395
      g(7) = pref * aa**7 * a**5 / 135135
      IF(any(g.gt.convpar)) THEN
        a = a + 1
        GOTO 2
      END IF
      rrad = a*scale

      IF(first) THEN
        WRITE(6,'(/A,2F10.5)') 'Cutoff radii: ',rad,rrad
        WRITE(6,'(/A)')        'Real-space sum'
      END IF

!
!     Determine atomic shells
      CALL getshells(ptsh,nptsh,radsh,nshell,rad,cell%amat,first)

      ALLOCATE( pnt(nptsh) )
      structconst = 0

!
!     Real-space sum
!
      CALL cpu_time(time1)
      
      DO ic2 = 1,atoms%nat
        DO ic1 = 1,atoms%nat
          IF(ic2.ne.1.and.ic1.eq.ic2) CYCLE
          rc = matmul(cell%amat,(atoms%taual(:,ic2)-atoms%taual(:,ic1)))
          DO i=1,nptsh
            ra       = matmul(cell%amat,ptsh(:,i)) + rc 
            a        = sqrt(sum(ra**2))
            radsh(i) = a
          END DO
          CALL rorderpf(pnt,radsh,nptsh,&
     &                  max(0,int(log(nptsh*0.001d0)/log(2d0))))
          ptsh   = ptsh (:,pnt)
          radsh  = radsh(  pnt)
          maxl   = 2*hybrid%lexp
          a1     = huge(a1)  ! stupid initial value
          ishell = 1
          conv   = huge(i)
          shlp   = 0
          DO i = 1,nptsh
            IF(all(conv.ne.huge(i))) EXIT
            IF(i.ne.1) THEN 
              IF(abs(radsh(i)-radsh(i-1)).gt.1d-10) ishell = ishell + 1 
            ENDIF
            ra = matmul(cell%amat,ptsh(:,i)) + rc
            a  = scale * sqrt(sum(ra**2))
            IF(a.eq.0)  THEN
              CYCLE
            ELSE IF(abs(a-a1).gt.1d-10) THEN
              a1   = a
              rexp = exp(-a)
              IF(ishell.le.conv(0)) g(0) = rexp / a    &
     & * (1+a*11/16*(1+a*3/11*(1+a/9)))
              IF(ishell.le.conv(1)) g(1) = rexp / a**2&
     & * (1+a*(1+a/2*(1+a*7/24*(1+a/7))))
              IF(ishell.le.conv(2)) g(2) = rexp / a**3&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a*3/16*(1+a/9))))))
              IF(ishell.le.conv(3)) g(3) = rexp / a**4&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/8)))))))
              IF(ishell.le.conv(4)) g(4) = rexp / a**5&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8&
     & *(1+a/10)))))))))
              IF(ishell.le.conv(5)) g(5) = rexp / a**6&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/9&
     & *(1+a/10))))))))))
              IF(ishell.le.conv(6)) g(6) = rexp / a**7&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/9&
     & *(1+a/10*(1+a/11*(1+a/12))))))))))))
              IF(ishell.le.conv(7)) g(7) = rexp / a**8&
     & * (1+a*(1+a/2*(1+a/3*(1+a/4*(1+a/5*(1+a/6*(1+a/7*(1+a/8*(1+a/9&
     & *(1+a/10*(1+a/11*(1+a/12*(1+a/13)))))))))))))
              DO l = 8,maxl
                IF(ishell.le.conv(l)) g(l) = a**(-l-1)
              END DO
              DO l = 0,maxl
                IF(conv(l).eq.huge(i).and.g(l).lt.convpar(l)/10)&
     &            conv(l) = ishell + ADDSHELL1
              END DO
            END IF
            IF(ishell.gt.conv(maxl).and.maxl.ne.0) maxl = maxl - 1
            CALL harmonicsr(y,ra,maxl)
            y = conjg(y)
            DO ikpt = 1,kpts%nkpt
              rdum = kpts%bk(1,ikpt)*ptsh(1,i) + kpts%bk(2,ikpt)*ptsh(2,i) &
     &             + kpts%bk(3,ikpt)*ptsh(3,i)
              cexp = exp(img*2*pi_const*rdum)
              lm   = 0
              DO l = 0,maxl
                IF(ishell.le.conv(l)) THEN
                  cdum = cexp * g(l)
                  DO M = -l,l
                    lm            = lm + 1
                    shlp(lm,ikpt) = shlp(lm,ikpt) + cdum * y(lm)
                  END DO
                ELSE
                  lm = lm + 2*l + 1
                END IF
              END DO
            END DO
          END DO
          structconst(:,ic1,ic2,:) = shlp
        END DO
      END DO

      DEALLOCATE( ptsh,radsh )

      CALL cpu_time(time2)
      IF(first) WRITE(6,'(A,F7.2)') '  Timing: ',time2-time1
      CALL cpu_time(time1)

      if(first) WRITE(6,'(/A)')  'Fourier-space sum'

!
!     Determine reciprocal shells
!
      call getshells(ptsh,nptsh,radsh,nshell,rrad,cell%bmat,first)
      ! minimum nonzero reciprocal-shell radius (needed in routines concerning the non-local hartree-fock exchange)
      hybrid%radshmin = radsh(2) 
!
!     Fourier-space sum
!
      DO ikpt=1,kpts%nkpt
        k         = kpts%bk(:,ikpt)
        maxl      = min(7,hybrid%lexp*2)
        ishell    = 1
        conv      = huge(i)
        DO i = 1,nptsh
          IF(i.gt.1) THEN 
            IF(abs(radsh(i)-radsh(i-1)).gt.1d-10) ishell = ishell + 1
          ENDIF
          ki = ptsh(:,i) + k - nint(k) ! -nint(...) transforms to Wigner-Seitz cell ( i.e. -0.5 <= x,y,z < 0.5 )
          ka = matmul(ki,cell%bmat)
          a  = sqrt(sum(ka**2)) / scale
          aa = (1+a**2)**(-1)
          IF(abs(a-a1).gt.1d-10) THEN
            a1 = a
            IF(a.eq.0) THEN
              g(0) = pref * (-4) 
              g(1) = 0
            ELSE
              IF(ishell.le.conv(0)) g(0) = pref * aa**4 / a**2
              IF(ishell.le.conv(1)) g(1) = pref * aa**4 / a
            END IF
            IF(ishell.le.conv(2))   g(2) = pref * aa**5        / 3
            IF(ishell.le.conv(3))   g(3) = pref * aa**5 * a    / 15
            IF(ishell.le.conv(4))   g(4) = pref * aa**6 * a**2 / 105
            IF(ishell.le.conv(5))   g(5) = pref * aa**6 * a**3 / 945
            IF(ishell.le.conv(6))   g(6) = pref * aa**7 * a**4 / 10395
            IF(ishell.le.conv(7))   g(7) = pref * aa**7 * a**5 / 135135
            IF(ishell.gt.1) THEN
              DO l = 0,7
                IF(conv(l).eq.huge(i).and.g(l).lt.convpar(l))&
     &            conv(l) = ishell + ADDSHELL2
              END DO
            END IF
          END IF

          IF(ishell.gt.conv(maxl).and.maxl.ne.0) maxl = maxl - 1
          CALL harmonicsr(y,ka,maxl)
          cdum = 1d0
          lm   = 0
          DO l = 0,maxl
            IF(ishell.le.conv(l)) THEN
              DO M = -l,l
                lm    = lm + 1
                y(lm) = conjg(y(lm)) * cdum * g(l)
              END DO
            ELSE
              y(lm+1:lm+2*l+1) = 0
              lm               = lm + 2*l + 1
            END IF
            cdum = cdum * img
          END DO
          DO ic2 = 1,atoms%nat
            DO ic1 = 1,atoms%nat
              IF(ic2.ne.1.and.ic1.eq.ic2) CYCLE
              cexp = exp(img*2*pi_const*&
     &          dot_product(ki,atoms%taual(:,ic1)-atoms%taual(:,ic2)))
              DO lm = 1,(maxl+1)**2
                structconst(lm,ic1,ic2,ikpt)&
     &          = structconst(lm,ic1,ic2,ikpt) + cexp * y(lm)
              END DO
            END DO
          END DO
        END DO
      END DO

      CALL cpu_time(time2) 
      IF(first) WRITE(6,'(A,F7.2)') '  Timing: ',time2-time1

!
!     Add contribution for l=0 to diagonal elements and rescale structure constants
!
      structconst(1,1,1,:) = structconst(1,1,1,:) - 5d0/16/sqrt(4*pi_const)
      DO i = 2,atoms%nat
        structconst(:,i,i,:) = structconst(:,1,1,:)
      END DO
      DO l=0,2*hybrid%lexp
        structconst(l**2+1:(l+1)**2,:,:,:)&
     &  = structconst(l**2+1:(l+1)**2,:,:,:) * scale**(l+1)
      END DO

      rad = (cell%vol*3/4/pi_const)**(1d0/3) ! Wigner-Seitz radius (rad is recycled)

!     Calculate accuracy of Gamma-decomposition 
      IF(all(kpts%bk.eq.0)) GOTO 4
      a = 1d30 ! ikpt = index of shortest non-zero k-point
      DO i = 2,kpts%nkpt 
        rdum = sum(matmul(kpts%bk(:,i),cell%bmat)**2)
        IF(rdum.lt.a) THEN 
          ikpt = i 
          a    = rdum 
        END IF
      END DO
      rdum = sqrt(sum(matmul(kpts%bk(:,ikpt),cell%bmat)**2))
      a    = 0
      DO ic2=1,atoms%nat
        DO ic1=1,max(1,ic2-1)
          a = a + abs( structconst(1,ic1,ic2,ikpt) -&
     &             ( structconst(1,ic1,ic2,1) + sqrt(4*pi_const)/cell%vol/rdum**2 *&
     &               exp(-img*2*pi_const*dot_product(&
     &                    kpts%bk(:,ikpt),atoms%taual(:,ic2)-atoms%taual(:,ic1)))))**2
        END DO
      END DO
      a  = sqrt(a/atoms%nat**2)
      aa = sqrt(sum(abs(structconst(1,:,:,ikpt))**2)/atoms%nat**2)
      IF(first)&
     &  WRITE(6,'(/A,F8.5,A,F8.5,A)') 'Accuracy of Gamma-decomposition&&
     & (structureconstant):',a,' (abs)',a/aa,' (rel)'

 4    DEALLOCATE (ptsh,radsh)

      first = .false.

      END SUBROUTINE structureconstant


!     -----------------

!     Determines all shells of the crystal defined by lat and vol with radii smaller than rad.
!     The lattice points (number = nptsh) are stored in ptsh, their corresponding lengths (shell radii) in radsh.

      SUBROUTINE getshells(ptsh,nptsh,radsh,nshell,rad,lat,lwrite)
      USE m_util     , ONLY : rorderpf
      IMPLICIT NONE
      LOGICAL , INTENT(IN)    :: lwrite
      INTEGER , INTENT(OUT)   :: nptsh,nshell
      INTEGER , ALLOCATABLE   :: ptsh(:,:)
      REAL    , ALLOCATABLE   :: radsh(:)
      REAL    , INTENT(IN)    :: rad,lat(3,3)
      REAL                    :: r(3),rdum
      INTEGER , ALLOCATABLE   :: pnt(:)
      INTEGER                 :: n,i,ix,iy,iz,ok
      LOGICAL                 :: found
      INTEGER , ALLOCATABLE   :: ihelp(:,:)
      REAL    , ALLOCATABLE   :: rhelp(:)

      ALLOCATE ( ptsh(3,100000),radsh(100000), stat=ok )
      IF( ok .ne. 0 ) STOP 'getshells: failure allocation ptsh/radsh'

      ptsh  = 0
      radsh = 0

      n = 0
      i = 0
      DO
        found = .false.
        DO ix = -n,n
          DO iy = -(n-abs(ix)),n-abs(ix)
            iz   = n-abs(ix)-abs(iy)
 1          r    = ix*lat(:,1) + iy*lat(:,2) + iz*lat(:,3)
            rdum = sum(r**2)
            IF(rdum.lt.rad**2) THEN
              found     = .true.
              i         = i + 1
              IF( i .gt. size(radsh) ) THEN 
                ALLOCATE( rhelp(size(radsh)), ihelp(3,size(ptsh,2)),&
     &                    stat = ok )
                IF( ok .ne. 0 )&
     &            STOP 'getshells: failure allocation rhelp/ihelp'
                rhelp = radsh
                ihelp = ptsh
                DEALLOCATE( radsh,ptsh )
                ALLOCATE( radsh(size(rhelp)+100000),&
     &                    ptsh(3,size(ihelp,2)+100000), stat=ok )
                IF( ok .ne. 0 )&
     &            STOP 'getshells: failure re-allocation ptsh/radsh'
                radsh(  1:size(rhelp)  ) = rhelp
                ptsh (:,1:size(ihelp,2)) = ihelp
                DEALLOCATE( rhelp,ihelp )
              END IF
              ptsh(:,i) = (/ ix,iy,iz /)
              radsh(i) = sqrt(rdum)
            END IF
            IF(iz.gt.0) THEN
              iz = -iz
              GOTO 1
            END IF
          END DO
        END DO
        IF(.not.found) EXIT
        n = n + 1
      END DO
      nptsh = i

      ALLOCATE( pnt(nptsh) )

      !reallocate radsh ptsh
      ALLOCATE( rhelp(nptsh), ihelp(3,nptsh) )
      rhelp = radsh(  1:nptsh)
      ihelp = ptsh (:,1:nptsh)
      DEALLOCATE( radsh,ptsh )
      ALLOCATE( radsh(nptsh),ptsh(3,nptsh) )
      radsh = rhelp
      ptsh  = ihelp
      DEALLOCATE( rhelp, ihelp )

      CALL rorderpf(pnt,radsh,nptsh,&
     &              max(0,int(log(nptsh*0.001d0)/log(2d0))))
      radsh = radsh(pnt)
      ptsh  = ptsh(:,pnt)
      nshell = 1
      DO i=2,nptsh
        IF(radsh(i)-radsh(i-1).gt.1d-10) nshell = nshell + 1
      END DO

      IF(lwrite)&
     &  WRITE(6,'(A,F10.5,A,I7,A,I5,A)')&
     &    '  Sphere of radius',rad,' contains',&
     &    nptsh,' lattice points and',nshell,' shells.'

      END SUBROUTINE getshells


!     ---------

!     Returns a list of (k+G) vector lengths in qnrm(1:nqnrm) and the corresponding pointer pqnrm(1:ngpt(ikpt),ikpt)
      SUBROUTINE getnorm(kpts,gpt,ngpt,pgpt, qnrm,nqnrm,pqnrm,cell)
        USE m_types
      IMPLICIT NONE
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_kpts),INTENT(IN)   :: kpts

      INTEGER , INTENT(IN)  :: ngpt(kpts%nkpt),gpt(3,*),pgpt(:,:)!(dim,kpts%nkpt)
      REAL    , ALLOCATABLE :: qnrm(:),help(:)
      INTEGER , INTENT(OUT) :: nqnrm
      INTEGER , ALLOCATABLE :: pqnrm(:,:)
      INTEGER               :: i,j,ikpt,igpt,igptp
      REAL                  :: q(3),qnorm

      ALLOCATE ( qnrm(maxval(ngpt)*kpts%nkpt),pqnrm(maxval(ngpt),kpts%nkpt) )
      i = 0
      DO ikpt = 1,kpts%nkpt
        DO igpt = 1,ngpt(ikpt)
          igptp = pgpt(igpt,ikpt)
          IF(igptp.eq.0) STOP 'getnorm: zero pointer (bug?)'
          q     = matmul (kpts%bk(:,ikpt) + gpt(:,igptp),cell%bmat)
          qnorm = sqrt(sum(q**2))
          DO j=1,i
            IF(abs(qnrm(j)-qnorm).lt.1d-12) THEN
              pqnrm(igpt,ikpt) = j
              GOTO 1
            END IF
          END DO
          i                = i + 1
          qnrm(i)          = qnorm
          pqnrm(igpt,ikpt) = i
 1      END DO
      END DO
      nqnrm = i

      ALLOCATE(help(nqnrm))
      help(1:nqnrm)=qnrm(1:nqnrm)
      DEALLOCATE(qnrm)
      ALLOCATE(qnrm(1:nqnrm))
      qnrm=help

      END SUBROUTINE getnorm


      FUNCTION sphbessel_integral(atoms,itype,qnrm,nqnrm,&
     &                  iqnrm1,iqnrm2,l,hybrid,sphbes0,l_warnin,l_warnout)

    USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_atoms),INTENT(IN)   :: atoms

      INTEGER , INTENT(IN)  :: itype ,nqnrm,iqnrm1,iqnrm2,l 
      REAL    , INTENT(IN)  :: qnrm(nqnrm), sphbes0(-1:hybrid%lexp+2,atoms%ntype,nqnrm)
      LOGICAL , INTENT(IN)  , OPTIONAL  ::  l_warnin
      LOGICAL , INTENT(OUT) , OPTIONAL  ::  l_warnout
      REAL                  :: sphbessel_integral
      REAL                  :: q1,q2,dq,s,sb01,sb11,sb21,sb31,sb02,sb12,&
     &                        sb22,sb32,a1,a2,da,b1,b2,db,c1,c2,dc,r1,r2
      LOGICAL               :: l_warn, l_warned

      IF ( present(l_warnin) ) THEN
        l_warn = l_warnin
      ELSE
        l_warn = .true.
      END IF
      l_warned = .false.

      q1 = qnrm(iqnrm1)
      q2 = qnrm(iqnrm2)
      s  = atoms%rmt(itype)
      IF(q1.eq.0.and.q2.eq.0) THEN
        IF(l.gt.0) THEN
          sphbessel_integral = 0
        ELSE
          sphbessel_integral = 2*s**5/15
        ENDIF
      ELSE IF(q1.eq.0.or.q2.eq.0) THEN
        IF(l.gt.0)       THEN
          sphbessel_integral = 0
        ELSE IF(q1.eq.0) THEN 
          sphbessel_integral = s**3/(3*q2**2)&
     &    * ( q2*s * sphbes0(1,itype,iqnrm2) + sphbes0(2,itype,iqnrm2) )
        ELSE
          sphbessel_integral = s**3/(3*q1**2)&
     &    * ( q1*s * sphbes0(1,itype,iqnrm1) + sphbes0(2,itype,iqnrm1) )
        ENDIF
      ELSE IF(q1.eq.q2) THEN
        sphbessel_integral = s**3/(2*q1**2)&
     &  * ( (2*l+3) * sphbes0(l+1,itype,iqnrm1)**2 -&
     &   (2*l+1) * sphbes0(l,itype,iqnrm1) * sphbes0(l+2,itype,iqnrm1) )
      ELSE ! We use either if two fromulas that are stable for high and small q1/q2 respectively
        sb01 = sphbes0(l-1,itype,iqnrm1)
        sb11 = sphbes0(l  ,itype,iqnrm1)
        sb21 = sphbes0(l+1,itype,iqnrm1)
        sb31 = sphbes0(l+2,itype,iqnrm1)
        sb02 = sphbes0(l-1,itype,iqnrm2)
        sb12 = sphbes0(l  ,itype,iqnrm2)
        sb22 = sphbes0(l+1,itype,iqnrm2)
        sb32 = sphbes0(l+2,itype,iqnrm2)
        dq   = q1**2 - q2**2
        a1   = q2/q1 * sb21 * sb02
        a2   = q1/q2 * sb22 * sb01
        da   = a1 - a2
        b1   = sb31 * sb12
        b2   = sb32 * sb11
        db   = b1 - b2
        c1   = sb21 * sb22 / ( q1*q2 )
        c2   = db / dq * (2*l+1)/(2*l+3)
        dc   = c1 + c2
        r1   = abs(da/a1)
        r2   = min ( abs(db/b1) , abs(dc/c1) )
        ! Ensure numerical stability. If both formulas are not sufficiently stable, the program stops.
        IF(r1.gt.r2) THEN
          IF(r1.lt.1d-6 .and. l_warn) THEN
            WRITE(6,'(A,E10.5,A,E10.5,A)')&
     &         'sphbessel_integral: Warning! Formula One possibly'//&
     &         ' unstable. Ratios:',r1,'(',r2,')'
            WRITE(6,'(A,2F15.10,I4)')&
     &         '                    Current qnorms and atom type:',&
     &         q1,q2,itype
            l_warned = .true.
          END IF
          sphbessel_integral = s**3 / dq * da 
        ELSE
          IF(r2.lt.1d-6 .and. l_warn) THEN
            WRITE(6,'(A,E10.5,A,E10.5,A)')&
     &         'sphbessel_integral: Warning! Formula Two possibly'//&
     &         ' unstable. Ratios:',r2,'(',r1,')'
            WRITE(6,'(A,2F15.10,I4)')&
     &         '                    Current qnorms and atom type:',&
     &         q1,q2,itype
            l_warned = .true.
          END IF
          sphbessel_integral = s**3 * dc
        END IF
      END IF

      IF ( present(l_warnout) ) l_warnout = l_warned

      END FUNCTION sphbessel_integral


!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      END MODULE m_coulomb
