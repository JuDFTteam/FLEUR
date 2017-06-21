!     Calculates the HF exchange term 
!
!                                          s          s*          s            s*
!                                       phi    (r) phi     (r) phi     (r') phi    (r')
!                         occ.             n_1k       n'k+q       n'k+q        n_2k
!     exchange(n,q)  =  - SUM  INT INT  ------------------------------------------- dr dr'
!                         k,n'                           | r - r' |
!
!                         occ                  s          s    ~        ~       s         s
!                    =  - SUM  SUM  v     < phi      | phi     M    > < M    phi     | phi      >
!                         k,n' I,J   k,IJ      n'k+q      n_1k  q,I      q,J    n_2k      n'k+q
!
!     for the different combinations of n_1 and n_2 and where n' runs only over the valence states.     
!     ( n_1,n_2:  valence-valence, core-core,core-valence )
!
!
!     At the Gamma point (k=0) v diverges. After diagonalization of v at k=0 the divergence is
!     restricted to the head element I=1. Furthermore, we expand <...> with kp perturbation theory.
!     As a result, the total I=1 element is given by a sum of a divergent 1/k**2-term and an
!     angular dependent term. The former is separated from the numerical k-summation and treated
!     analytically while the latter is spherically averaged and added to the k=0 contribution of
!     the numerical k-summation. (A better knowledge of the integrand's behavior at the BZ edges
!     might further improve the integration.)
!
!     The divergence at the Gamma point is integrated with one of the following algorithms:
! (1) Switching-Off Function
!     In a sphere of radius k0=radshmin/2 a switching-off function g(k)=1-(k/k0)**n*(n+1-n*k/k0)
!     (n=npot) is defined. The 1/k**2 divergence is subtracted from the BZ integral in the form
!     g(k)/k**2 and integrated analytically. The non-divergent rest is integrated numerically.
! (2) Periodic Function (similar to the one used by Massidda PRB 48, 5058)
!     The function  F(k) = SUM(G) exp(-expo*|k+G|**3) / |k+G|**2  is subtracted from the BZ integral
!     and integrated analytically. The non-divergent rest is integrated numerically.
!     The parameter expo is chosen such that exp(-expo*q**3)=1/2
!     with q = radius of sphere with same volume as BZ.
! (3) Periodic Function (same as Massidda's) with expo->0
!     The function  F(k) = lim(expo->0) SUM(G) exp(-expo*|k+G|**2) / |k+G|**2  is subtracted from
!     the BZ integral and integrated analytically. The contribution to the BZ integral including
!     the "tail" is
!     vol/(8*pi**3) INT F(k) d^3k - P SUM(k) F(k)  ( P = principal value ) .
!     For expo->0 the two terms diverge. Therefore a cutoff radius q0 is introduced and related to
!     expo by exp(-expo*q0**2)=delta  ( delta = small value, e.g., delta = 1d-10 ) .
!     The resulting formula
!     vol/(4*pi**1.5*sqrt(expo)) * erf(sqrt(a)*q0) - sum(q,0<q<q0) exp(-expo*q**2)/q**2
!     converges well with q0. (Should be the default.)
      MODULE m_exchange_valence_hf


#define ALGORITHM 3
#define zero_order .false.
#define ibs_corr   .false.
#define maxmem      600

      CONTAINS

      SUBROUTINE exchange_valence_hf(&
                    nk,kpts,nkpti,nkpt_EIBZ,&
                    sym,atoms,hybrid,&
                    cell,&
                    dimension,input,jsp,&
                    basm,bas1,bas2,bas1_MT,&
                    drbas1_MT,maxlcutm,lcutm,nindxm,maxindxm,nbasp,&
                    nbasm,maxbasm,maxindxp,nindxp,&
                    prod,prodm,mnobd,&
                    nobd,nbands,ne_eig,lapw,&
                    eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,&
                    kveclo,gauntarr,it,xcpot,&
                    noco,nsest,indx_sest,&
                    mpi,irank2,isize2,comm,&
                    div_vv,mat_ex)


      USE m_wrapper
      USE m_constants   
      USE m_trafo
      USE m_util          ,ONLY: cerf
      USE m_wavefproducts
      USE m_olap
      USE m_spmvec
      USE m_hsefunctional ,ONLY: dynamic_hse_adjustment
#ifdef CPP_MPI
      USE m_mpi_work_dist
      USE m_mpi_tags
#endif
      USE m_icorrkeys
      USE m_kp_perturbation
      USE m_types
      IMPLICIT NONE
      TYPE(t_results),INTENT(IN)   :: results
      TYPE(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_noco),INTENT(IN)   :: noco
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)   :: lapw

!     - scalars -
      INTEGER,INTENT(IN)      :: it  ,irank2 ,isize2,comm
      INTEGER,INTENT(IN)      :: jsp
      INTEGER,INTENT(IN)      ::  maxbasm,maxlcutm,maxindxm
      INTEGER,INTENT(IN)      ::  maxindxp
      INTEGER,INTENT(IN)      ::  nk ,nkpti ,nkpt_EIBZ
      INTEGER,INTENT(IN)      ::  nbasp  ,ne_eig   
      INTEGER,INTENT(IN)      :: mnobd,nbands



!     - arrays -
      INTEGER,INTENT(IN)      :: lcutm(atoms%ntype)
      INTEGER,INTENT(IN)      ::  nobd(kpts%nkptf),nbasm(kpts%nkptf)
      INTEGER,INTENT(IN)      ::  n_q(nkpt_EIBZ)
      INTEGER,INTENT(IN)      ::  nindxm(0:maxlcutm,atoms%ntype),&
                                  nindxp(0:maxlcutm,atoms%ntype)
      INTEGER,INTENT(IN)      ::  parent(kpts%nkptf)
      INTEGER,INTENT(IN)      ::  pointer_EIBZ(nkpt_EIBZ)
      INTEGER,INTENT(IN)      ::  kveclo(atoms%nlotot,nkpti)
      INTEGER,INTENT(IN)      ::  nsest(nbands),indx_sest(nbands,nbands)

      REAL   ,INTENT(IN)      ::  basm(atoms%jmtd,maxindxm,0:maxlcutm,atoms%ntype)
      REAL   ,INTENT(IN)      ::  bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
                                  bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
      REAL   ,INTENT(IN)      ::    bas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
                                  drbas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)   
      REAL   ,INTENT(IN)      ::  eig_irr(dimension%neigd,nkpti)
      REAL   ,INTENT(IN)      ::  gauntarr(2,0:atoms%lmaxd,0:atoms%lmaxd,0:maxlcutm,&
                                        -atoms%lmaxd:atoms%lmaxd,-maxlcutm:maxlcutm)
      REAL   ,INTENT(IN)      ::  prodm(maxindxm,maxindxp,&
                                        0:maxlcutm,atoms%ntype)
      REAL   ,INTENT(IN)      ::  wl_iks(dimension%neigd,kpts%nkptf)
      REAL   ,INTENT(OUT)     ::  div_vv(nbands)

#ifdef CPP_INVERSION
      REAL   ,INTENT(OUT)     ::  mat_ex(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#else
      COMPLEX,INTENT(OUT)     ::  mat_ex(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#endif
      TYPE(PRODTYPE)          ::  prod(maxindxp,0:maxlcutm,atoms%ntype)

!     - local scalars -
      INTEGER                 ::  iband,iband1,ibando,ikpt,ikpt0
      INTEGER                 ::  i,ic,ix,iy,iz
      INTEGER                 ::  irecl_coulomb,irecl_coulomb1
      INTEGER                 ::  j
      INTEGER                 :: m1,m2
      INTEGER                 ::  n,n1,n2,nn,nn2
      INTEGER                 ::  nkqpt
      INTEGER                 ::  npot
      INTEGER                 ::  ok
      INTEGER                 ::  psize
      INTEGER                 ::  iqptmin,iqptmax
#ifdef CPP_INVERSION
      INTEGER                 ::  bytes = 8
#else
      INTEGER                 ::  bytes = 16
#endif
      REAL                    :: svol
      REAL                    ::  rws,rrad,rdum
      REAL                    ::  k0,knorm,knorm2
      REAL                    ::  time1,time2,time3
      REAL                    ::  constant1,constant2,constant3
      REAL                    ::  time_mt,time_ir,vol
      REAL                    ::  expo
      REAL , SAVE             ::  divergence

      COMPLEX                 ::  cdum,cdum1,cdum2 
      COMPLEX                 ::  exch0

      LOGICAL                 ::  found
      LOGICAL, SAVE           ::  initialize = .true.

!     - local arrays -
      INTEGER                 ::  kcorner(3,8) = reshape((/ 0,0,0, 1,0,0, 0,1,0, 0,0,1,&
                                             1,1,0, 1,0,1, 0,1,1, 1,1,1 /), (/3,8/) )
      REAL                    ::  k(3),kv1(3),kv2(3),kv3(3),kvec(3)

      COMPLEX,ALLOCATABLE     ::  phase_vv(:,:)
      COMPLEX                 ::  exchcorrect(kpts%nkptf)
      COMPLEX                 ::  dcprod(nbands,nbands,3) 

      COMPLEX(8)              ::  exch_vv(nbands,nbands)
#ifdef CPP_MPI
      COMPLEX(8)              ::  buf_vv(nbands,nbands)
#endif
      COMPLEX                 ::  hessian(3,3)
      COMPLEX                 ::  proj_ibsc(3,mnobd,nbands)
      COMPLEX                 ::  olap_ibsc(3,3,mnobd,mnobd)
#if ( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
      REAL                    ::  coulomb_mt1(maxindxm-1,maxindxm-1, 0:maxlcutm,atoms%ntype)       
#ifdef CPP_INVERSION
      REAL                    ::  coulomb_mt2(maxindxm-1, -maxlcutm:maxlcutm, 0:maxlcutm+1,atoms%nat)
      REAL                    ::  coulomb_mt3(maxindxm-1,atoms%nat,atoms%nat)
#else
      COMPLEX                 ::  coulomb_mt2(maxindxm-1, -maxlcutm:maxlcutm, 0:maxlcutm+1,atoms%nat)
      COMPLEX                 ::  coulomb_mt3(maxindxm-1,atoms%nat,atoms%nat)
#endif

#else

#ifdef CPP_INVERSION
      REAL                    ::  coulomb(maxbasm*(maxbasm+1)/2)
#else
      COMPLEX                 ::  coulomb(maxbasm*(maxbasm+1)/2)
#endif 

#endif

#if ( defined(CPP_INVERSION) )
      REAL   ,ALLOCATABLE     ::  cprod_vv(:,:,:),cprod_cv(:,:,:), carr3_vv(:,:,:),carr3_cv(:,:,:)
      REAL                    ::  carr1_v(maxbasm),carr1_c(maxbasm)
#ifdef CPP_IRCOULOMBAPPROX
      REAL                    ::  coulomb_mtir((maxlcutm+1)**2* , (maxlcutm+1)**2* +maxval(hybrid%ngptm) )
#else
      REAL                    ::  coulomb_mtir(((maxlcutm+1)**2* +maxval(hybrid%ngptm))* ((maxlcutm+1)**2* +maxval(hybrid%ngptm)+1)/2 )
#endif

#else
      COMPLEX,ALLOCATABLE     ::  cprod_vv(:,:,:),cprod_cv(:,:,:), carr3_vv(:,:,:),carr3_cv(:,:,:)
      COMPLEX                 ::  carr1_v(maxbasm),carr1_c(maxbasm)

#ifdef CPP_IRCOULOMBAPPROX
      COMPLEX                 ::  coulomb_mtir((maxlcutm+1)**2* , (maxlcutm+1)**2* +maxval(hybrid%ngptm) )
#else
      COMPLEX                 ::  coulomb_mtir(((maxlcutm+1)**2* +maxval(hybrid%ngptm))* ((maxlcutm+1)**2* +maxval(hybrid%ngptm)+1)/2 )
#endif

#endif
      LOGICAL                 ::  occup(dimension%neigd)
#ifdef CPP_MPI
      INCLUDE "mpif.h"
      INTEGER                 :: ierr,ierr2,length,rank
      CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: errmsg
#endif
      time_mt = 0
      time_ir = 0

      vol  = cell%omtil
      svol = sqrt(cell%omtil)

      rws  = (3*cell%omtil/fpi_const)**(1d0/3)  ! Wigner-Seitz radius
#if   ALGORITHM == 1
      npot = 3                        ! for switching-off function 
      k0   = hybrid%radshmin / 2             ! radius of largest sphere that fits inside the BZ
#elif ALGORITHM == 3
      IF( initialize ) THEN !it .eq. 1 .and. nk .eq. 1) THEN
!         CALL cpu_time(time1)
        expo       = 5d-3
        rrad       = sqrt(-log(5d-3)/expo)
        cdum       = sqrt(expo)*rrad
        divergence = vol / (tpi_const**2) * sqrt(pi_const/expo) * cerf(cdum)
        rrad       = rrad**2
        kv1        = cell%bmat(1,:)/kpts%nkpt3(1)
        kv2        = cell%bmat(2,:)/kpts%nkpt3(2)
        kv3        = cell%bmat(3,:)/kpts%nkpt3(3)
        n          = 1
        found      = .true.
        DO WHILE(found)
          found = .false.
          DO ix = -n,n
            DO iy = -(n-abs(ix)),n-abs(ix)
              iz     = n - abs(ix) - abs(iy)
 1            k(1)   = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
              k(2)   = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
              k(3)   = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
              knorm2 = k(1)**2   + k(2)**2   + k(3)**2
              IF(knorm2.lt.rrad) THEN
                found      = .true.
                divergence = divergence &
                           - exp(-expo*knorm2)/knorm2 / kpts%nkptf
              END IF
              IF(iz.gt.0) THEN
                iz = -iz
                GOTO 1
              END IF
            END DO
          END DO
          n = n + 1
        END DO
!         CALL cpu_time(time2)
!         WRITE(*,*) 'time for calculating periodic function',time2-time1
        initialize = .false.
      END IF
#endif

#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )

#ifdef CPP_INVERSION
      irecl_coulomb1 = ( atoms%ntype*(maxlcutm+1)*(maxindxm-1)**2&
                    +    atoms%nat *(maxlcutm+2)*(2*maxlcutm+1)*(maxindxm-1)&
                    +    (maxindxm-1)*atoms%nat**2&
                    +    ((maxlcutm+1)**2*atoms%nat+maxval(hybrid%ngptm))&
                    *    ((maxlcutm+1)**2*atoms%nat+maxval(hybrid%ngptm)+1)/2) *8
#else
      irecl_coulomb1 = ( atoms%ntype*(maxlcutm+1)*(maxindxm-1)**2&
                     +   atoms%nat *(maxlcutm+2)*(2*maxlcutm+1)*(maxindxm-1)&
                     +   (maxindxm-1)*atoms%nat**2&
                     +   ((maxlcutm+1)**2*atoms%nat+maxval(hybrid%ngptm))&
                     *   ((maxlcutm+1)**2*atoms%nat+maxval(hybrid%ngptm)+1)/2) *16
#endif
      OPEN(unit=676,file='coulomb1',form='unformatted',access='direct',&
           recl=irecl_coulomb1)

#else

      !open direct acces file coulomb/cprod
#ifdef CPP_INVERSION
      irecl_coulomb =  maxbasm*(maxbasm+1)*4 !(maxbasm*maxbasm)* 8+maxbasm*8 + 8
#else
      irecl_coulomb =  maxbasm*(maxbasm+1)*8!(maxbasm*maxbasm)*16+maxbasm*8 + 8
#endif

      OPEN(unit=677,file='coulomb',form='unformatted',access='direct',&
           recl=irecl_coulomb)

#endif

      ! calculate valence-valence-valence-valence, core-valence-valence-valence
      ! and core-valence-valence-core exchange at current k-point
      ! the sum over the inner occupied valence states is restricted to the EIBZ(k)
      ! the contribution of the Gamma-point is treated separately (see below)


      ! determine package size loop over the occupied bands
      rdum  = maxbasm*nbands*bytes/1048576.
      psize = 1
      DO iband = mnobd,1,-1
        ! ensure that the packages have equal size
        IF( modulo(mnobd,iband) .eq. 0 ) THEN
          ! choose packet size such that cprod is smaller than memory threshold
          IF( rdum*iband .le. maxmem ) THEN
            psize = iband
            EXIT
          END IF
        END IF
      END DO

      IF( psize .ne. mnobd ) THEN
        WRITE(6,'(A,A,i3,A,f7.2,A)') ' Divide the loop over the occupied bands in packages', ' of the size',psize,' (cprod=',rdum*psize,'MB)'
      END IF

      ALLOCATE( cprod_vv(maxbasm,psize,nbands),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation cprod'
      ALLOCATE( carr3_vv(maxbasm,psize,nbands),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation carr3'
      ALLOCATE( phase_vv(psize,nbands),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation phase'
      cprod_vv = 0 ; carr3_vv = 0 ; phase_vv = 0

      CALL cpu_time(time1) 
      exch_vv = 0

#     ifndef CPP_MPI
        iqptmin = 1
        iqptmax = nkpt_EIBZ
#     else
        ! read the limits for current k-point
        CALL work_dist_nqpt_limits(nk,iqptmax,iqptmin)
#     endif
      DO ikpt = iqptmin,iqptmax
        ikpt0 = pointer_EIBZ(ikpt)

        n  = nbasp + hybrid%ngptm(ikpt0)
        IF( nbasm(ikpt0) .ne. n ) STOP 'error nbasm'
        nn = n*(n+1)/2

        ! read in coulomb matrix from direct access file coulomb
#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
        READ(676,rec=kpts%bkp(ikpt0)) coulomb_mt1,coulomb_mt2,&
                                 coulomb_mt3,coulomb_mtir
#else
        READ(677,rec=kpts%bkp(ikpt0)) coulomb
#endif

        IF( kpts%bkp(ikpt0) .ne. ikpt0 ) THEN
#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )

#ifndef CPP_INVERSION
          IF( kpts%bksym(ikpt0) .gt. sym%nop ) THEN
!             coulomb_mt1 = conjg(coulomb_mt1)
            coulomb_mt2 = conjg(coulomb_mt2)
            coulomb_mtir= conjg(coulomb_mtir)
          END IF
#endif

#else

#ifndef CPP_INVERSION
          IF( kpts%bksym(ikpt0) .gt. sym%nop ) coulomb = conjg(coulomb)
#endif

#endif
        END IF

        DO ibando = 1,mnobd,psize
#ifdef CPP_INVERSION

#ifdef CPP_IRAPPROX
          CALL wavefproducts_inv(&
                         1,nbands,dimension,jsp,atoms,&
                         lapw,obsolete,nkpti,kpts,kpts,nkpt_EIBZ,&
                         nk,ikpt0,nobd,mnobd,hybrid,maxbasm,&
                         parent,cell,&
                         lcutm,maxlcutm,maxindxp,&
                         maxindxm,nindxp,nindxm,&
                         prodm,prod,gauntarr,nbasp,sym,&
                         time_mt,time_ir,nkqpt,cprod_vv)
#else
          CALL wavefproducts_inv5(&
                         1,nbands,ibando,ibando+psize-1,&
                         dimension,input,jsp,atoms,&
                         lapw,obsolete,nkpti,kpts,kpts,nkpt_EIBZ,&
                         nk,ikpt0,nobd,mnobd,hybrid,maxbasm,&
                         parent,cell,&
                         lcutm,maxlcutm,maxindxp,&
                         maxindxm,nindxp,nindxm,&
                         prodm,prod,gauntarr,nbasp,sym,&
                         noco,noco,&
                         time_mt,time_ir,nkqpt,cprod_vv)
#endif

#else
#ifdef CPP_IRAPPROX
          CALL wavefproducts_noinv(&
                         1,nbands,nk,ikpt0,dimension,jsp,&
                         cell,atoms,hybrid,nindxm,&
                         lcutm,maxlcutm,maxindxp,nindxp,gauntarr,&
                         kpts,maxindxm,&
                         maxbasm,prod,prodm,mnobd,&
                         lapw,sym,&
                         nobd,nbasp,nkqpt,&
                         cprod_vv)
#else
          CALL wavefproducts_noinv5(&
                         1,nbands,ibando,ibando+psize-1,&
                         nk,ikpt0,dimension,input,jsp, &!jsp,&
                         cell,atoms,hybrid,nindxm,&
                         lcutm,maxlcutm,maxindxp,nindxp,gauntarr,&
                         kpts,maxindxm,&
                         maxbasm,prod,prodm,mnobd,&
                         lapw,sym,&
                         nobd,nbasp,&
                         noco,&
                         nkqpt,cprod_vv)
#endif

#endif

          ! The sparse matrix technique is not feasible for the HSE
          ! functional. Thus, a dynamic adjustment is implemented
          ! The mixed basis functions and the potential difference
          ! are Fourier transformed, so that the exchange can be calculated
          ! in Fourier space
#ifndef CPP_NOSPMVEC
          IF ( xcpot%icorr == icorr_hse .OR. xcpot%icorr == icorr_vhse ) THEN
            iband1  = nobd(nkqpt)
            exch_vv = exch_vv + dynamic_hse_adjustment(&
                       atoms%rmsh,atoms%rmt,atoms%dx,atoms%jri,atoms%jmtd,kpts%bk(:,ikpt0),ikpt0,kpts%nkptf,&
                       cell%bmat,vol,atoms%ntype,atoms%neq,atoms%nat,atoms%taual,lcutm,maxlcutm,&
                       nindxm,maxindxm,hybrid%gptm,hybrid%ngptm(ikpt0),hybrid%pgptm(:,ikpt0),&
                       hybrid%gptmd,basm,nbasm(ikpt0),iband1,nbands,nsest,&
                       ibando,psize,indx_sest,atoms%invsat,sym%invsatnr,mpi%irank,&
                       cprod_vv(:nbasm(ikpt0),:,:),&
                       wl_iks(:iband1,nkqpt),n_q(ikpt))
          END IF
#endif

          ! the Coulomb matrix is only evaluated at the irrecuible k-points
          ! bra_trafo transforms cprod instead of rotating the Coulomb matrix
          ! from IBZ to current k-point
          IF( kpts%bkp(ikpt0) .ne. ikpt0 ) THEN
             STOP "INTERFACE to bra_trafo2"
          !  CALL bra_trafo2(&
          !      carr3_vv(:nbasm(ikpt0),:,:),cprod_vv(:nbasm(ikpt0),:,:),&
          !      nbasm(ikpt0),psize,nbands,&
          !      kpts(ikpt0),ikpt0(ikpt0),sym,&
          !      hybrid,cell,maxlcutm,atoms,&
          !      lcutm,nindxm,maxindxm,nw,obsolete,&
          !      nbasp,&
          !      phase_vv)

            cprod_vv(:nbasm(ikpt0),:,:) = carr3_vv(:nbasm(ikpt0),:,:)
          ELSE
            phase_vv(:,:) = (1d0,0d0)
          END IF

          ! calculate exchange matrix at ikpt0

          DO n1=1,nbands
            DO iband = 1,psize
              IF( ibando + iband - 1 .gt. nobd(nkqpt) ) CYCLE
              cdum  = wl_iks(ibando+iband-1,nkqpt)&
                    * conjg(phase_vv(iband,n1))/n_q(ikpt)

#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
              carr1_v(:n) = 0 
              CALL spmvec(atoms,lcutm,maxlcutm,nindxm,maxindxm,&
                          nbasm(ikpt0),nbasp,hybrid,ikpt0,kpts,&
                          cell,&
                          coulomb_mt1,coulomb_mt2,coulomb_mt3,&
                          coulomb_mtir,cprod_vv(:n,iband,n1),&
                          carr1_v(:n))
#else
              carr1_v(:n) = matvec( coulomb(:nn),cprod_vv(:n,iband,n1) )
#endif

              DO n2=1,nsest(n1)!n1
                nn2 = indx_sest(n2,n1)
                exch_vv(nn2,n1) = exch_vv(nn2,n1) &
                                + cdum*phase_vv(iband,nn2)&
                                *dotprod( carr1_v(:n), &
                                          cprod_vv(:n,iband,nn2) )

              END DO !n2
            END DO
          END DO  !n1
        END DO !ibando
      END DO  !ikpt

      ! close direct access file coulomb
#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
      CLOSE(676)
#else
      CLOSE(677)
#endif

!
!     Send data from all processes to master of the subgroup
!
#ifdef CPP_MPI
      IF ( irank2 == 0 ) THEN
        ierr = 0
        DO rank = 1, isize2-1
          buf_vv = 0
          CALL MPI_RECV(buf_vv,nbands*nbands,MPI_COMPLEX16,rank,&
                        TAG_SNDRCV_EXCH_VV,comm,MPI_STATUS_IGNORE,ierr)
          exch_vv = exch_vv + buf_vv
        END DO
      ELSE
        CALL MPI_BSEND(exch_vv,nbands*nbands,MPI_COMPLEX16,0,&
                       TAG_SNDRCV_EXCH_VV,comm,ierr)
      END IF
      IF ( ierr /= 0 ) THEN
        CALL MPI_ERROR_STRING( ierr, errmsg, length, ierr2 )
        WRITE(*,*) errmsg
        STOP
      END IF
      CALL cpu_time(time2)
      IF ( irank2 /= 0 ) RETURN
#endif

      CALL cpu_time(time3)

      !CALL outtime('             time for calculating cprod(IR):', time_ir,mpi,mpi)
      !CALL outtime('             time for calculating cprod(MT):', time_mt,mpi,mpi)
      !CALL outtime('          time for calculating cprod:', time_ir+time_mt,mpi,mpi)

      !CALL outtime('          time for all k-points except G point:', time3-time1-time_ir-time_mt,mpi,mpi)


      !
      ! add contribution of the gamma point to the different cases (exch_vv,exch_cv,exch_cc)
      !
      constant1 =  fpi_const/vol
      constant2 =  fpi_const/sqrt(vol)
      constant3 =  fpi_const/3

      ! valence-valence-valence-valence exchange

      IF ( xcpot%icorr .NE. icorr_hse .AND. xcpot%icorr   .NE. icorr_vhse ) THEN ! no gamma point correction needed for HSE functional
        IF( zero_order .and. .not. ibs_corr ) THEN
          WRITE(6,'(A)') ' Take zero order terms into account.'
        ELSE IF( zero_order .and.  ibs_corr ) THEN
          WRITE(6,'(A)') ' Take zero order terms and ibs-correction into account.'
        END IF
        IF( zero_order ) THEN
          CALL dwavefproducts(  &
                            dcprod,nk,1,nbands,1,nbands,.false., atoms,hybrid,&
                            cell, nbands,kpts,nkpti,lapw,&
                            bas1,bas2,dimension,jsp,&
                            eig_irr,ne_eig)

          ! make dcprod hermitian
          DO n1 = 1,nbands
            DO n2 = 1,n1
              dcprod(n1,n2,:) = (dcprod(n1,n2,:) &
                              - conjg(dcprod(n2,n1,:)))/2   
              dcprod(n2,n1,:) = -conjg(dcprod(n1,n2,:))
            END DO
          END DO

          IF( ibs_corr ) THEN
            CALL ibs_correction(&
                        nk,atoms,&
                        dimension,input,jsp,&
                        bas1,bas2,bas1_MT,drbas1_MT,hybrid,&
                        lapw,kpts,nkpti,&
                        nbands,cell,mnobd,&
                        sym,kveclo,&
                        proj_ibsc,olap_ibsc)
          END IF

        END IF


        occup = .false.
        DO i=1,ne_eig
          IF ( results%ef  .ge. eig_irr(i,nk) ) THEN
            occup(i) = .true.
          ELSE IF ( eig_irr(i,nk) - results%ef .le. 1E-06) THEN
             occup(i) = .true.
          END IF
        END DO


        DO n1 = 1,nbands
          DO n2 = 1,nsest(n1)!n1
            nn2 = indx_sest(n2,n1)
            exchcorrect = 0
            exch0       = 0

            ! if zero_order = .true. add averaged k-dependent term to the numerical
            ! integration at Gamma-point contribution
            !
            ! if we start with a system with a small DFT band gap (like GaAs), the contribution
            ! of the highest occupied and lowest unoccupied state in Hessian is typically
            ! large; a correct numerical integration requires a dense k-point mesh, so
            ! we don't add the contribution exchcorrect for such materials 

            IF( zero_order ) THEN
              hessian = 0
              IF( occup(n1) .and. occup(nn2) ) THEN
                DO i = 1,3
                  j = i

                  DO iband = 1,nbands
                    IF( occup(iband) ) THEN
                      hessian(i,j) = hessian(i,j) + conjg(dcprod(iband,n1,i)) *dcprod(iband,nn2,j)
                    END IF
                    hessian(i,j) = hessian(i,j) - dcprod(iband,nn2,i) * conjg(dcprod(iband,n1,j))
                  END DO

                  ! ibs correction
                  IF( ibs_corr ) THEN 
                    hessian(i,j) = hessian(i,j) - olap_ibsc(i,j,n1,nn2)/vol
                    DO iband = 1,nbands
                      hessian(i,j) = hessian(i,j) + conjg(proj_ibsc(i,nn2,iband)) * proj_ibsc(j,n1,iband)/vol
                    END DO
                  END IF

                END DO
              ELSE

                DO i = 1,3
                  j = i 
                  DO iband = 1,nbands
                    IF( occup(iband) ) THEN
                      hessian(i,j) = hessian(i,j) + conjg(dcprod(iband,n1,i)) * dcprod(iband,nn2,j)
                    END IF
                  END DO
                END DO

              END IF

              exchcorrect(1) = constant3 * (hessian(1,1)+hessian(2,2)+hessian(3,3))
              exch0          = exchcorrect(1)/kpts%nkptf
            END IF


            ! tail correction/contribution from all other k-points (it  goes into exchcorrect )
#if ALGORITHM == 1
            DO ikpt = 2,kpts%nkptf

              ! Calculate distances from the eight reciprocal unit-cell corners
              knorm = k0
              DO i = 1,8
                rdum=sqrt(sum(matmul(kpts%bk(:,ikpt)-kcorner(:,i),cell%bmat)**2))
                IF(rdum.lt.k0) THEN
                  knorm = rdum
                  kvec  = ( kpts%bk(:,ikpt) - kcorner(:,i) ) / knorm
                END IF
              END DO

              ! The tail of the divergent term goes into exchcorrect.
              IF(knorm.lt.k0) THEN
                rdum = 1 - (knorm/k0)**npot * (npot+1-npot*knorm/k0)

                IF ( (n1 .eq. nn2) .and.  occup(n1)  ) THEN 
                  exchcorrect(ikpt) = - rdum * constant1 / knorm**2 
                END IF

                ! the contribution in the case unoccupied/occupied band proportional to 1/k
                ! vanishes after summing up over all k-points

              END IF

            END DO
#endif

            ! Analytic contribution

            cdum2 = 0
#if   ALGORITHM == 1
            IF ( (n1 .eq. nn2) .and. occup(n1) ) THEN
              cdum2 = 2*k0/pi_const * npot/(npot+2) + sum(exchcorrect(2:))/kpts%nkptf ! the tail is subtracted here
            END IF
#elif ALGORITHM == 3
            !multiply divergent contribution with occupation number;
            !this only affects metals 
            IF ( n1 .eq. nn2 ) THEN
               cdum2 = fpi_const/vol * divergence * wl_iks(n1,nk)*kpts%nkptf
            END IF
#endif

            ! due to the symmetrization afterwards the factor 1/n_q(1) must be added

            IF( n1 .eq. nn2 ) div_vv(n1) = real(cdum2) 

            exch_vv(nn2,n1)  = exch_vv(nn2,n1) + (exch0 + cdum2)/n_q(1)

          END DO !n2
        END DO !n1
      ELSE
        div_vv = 0.
      END IF ! xcpot%icorr .ne. icorr_hse


#ifdef CPP_INVERSION
      IF(any( abs(aimag(exch_vv)) .gt. 1E-08)) STOP 'exchange: unusally large imaginary part of exch_vv'
#endif

      ! write exch_vv in mat_ex
      ic = 0
      DO n1=1,nbands
        DO n2=1,n1
          ic = ic + 1
          mat_ex(ic) = mat_ex(ic) + exch_vv(n2,n1)
        END DO
      END DO

      END SUBROUTINE exchange_valence_hf


      END MODULE m_exchange_valence_hf
