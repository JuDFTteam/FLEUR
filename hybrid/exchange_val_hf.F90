!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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

   LOGICAL,PARAMETER:: zero_order=.false.,ibs_corr=.false.
   INTEGER,PARAMETER:: maxmem=600

CONTAINS

SUBROUTINE exchange_valence_hf(nk,kpts,nkpt_EIBZ,sym,atoms,hybrid,cell,dimension,input,jsp,hybdat,mnobd,lapw,&
                               eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,it,xcpot, noco,nsest,indx_sest,&
                               mpi,mat_ex)

   USE m_types
   USE m_wrapper
   USE m_constants   
   USE m_trafo
   USE m_wavefproducts
   USE m_olap
   USE m_spmvec
   USE m_hsefunctional ,ONLY: dynamic_hse_adjustment
#if defined(CPP_MPI)&&defined(CPP_NEVER)
   USE m_mpi_work_dist
   USE m_mpi_tags
#endif
   USE m_io_hybrid
   USE m_kp_perturbation

   IMPLICIT NONE

   TYPE(t_results),       INTENT(IN)    :: results
   TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
   TYPE(t_mpi),           INTENT(IN)    :: mpi
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_hybrid),        INTENT(INOUT) :: hybrid
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_noco),          INTENT(IN)    :: noco
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_lapw),          INTENT(IN)    :: lapw
   TYPE(t_mat),           INTENT(INOUT) :: mat_ex
   TYPE(t_hybdat),        INTENT(INOUT) :: hybdat

   ! scalars
   INTEGER,               INTENT(IN)    :: it
   INTEGER,               INTENT(IN)    :: jsp
   INTEGER,               INTENT(IN)    :: nk,nkpt_EIBZ
   INTEGER,               INTENT(IN)    :: mnobd 

   ! arrays
   INTEGER,               INTENT(IN)    ::  n_q(nkpt_EIBZ)

   INTEGER,               INTENT(IN)    ::  parent(kpts%nkptf)
   INTEGER,               INTENT(IN)    ::  pointer_EIBZ(nkpt_EIBZ)
   INTEGER,               INTENT(IN)    ::  nsest(hybrid%nbands(nk))
   INTEGER,               INTENT(IN)    ::  indx_sest(hybrid%nbands(nk),hybrid%nbands(nk))

 
   REAL,                  INTENT(IN)    ::  eig_irr(dimension%neigd,kpts%nkpt)
   REAL,                  INTENT(IN)    ::  wl_iks(dimension%neigd,kpts%nkptf)    

   ! local scalars
   INTEGER                 ::  iband,iband1,ibando,ikpt,ikpt0
   INTEGER                 ::  i,ic,ix,iy,iz
   INTEGER                 ::  irecl_coulomb,irecl_coulomb1
   INTEGER                 ::  j
   INTEGER                 ::  m1,m2
   INTEGER                 ::  n,n1,n2,nn,nn2
   INTEGER                 ::  nkqpt
   INTEGER                 ::  npot
   INTEGER                 ::  ok
   INTEGER                 ::  psize
   REAL                    ::  rdum
   REAL                    ::  k0
     
   REAL , SAVE             ::  divergence

   COMPLEX                 ::  cdum,cdum1,cdum2 
   COMPLEX                 ::  exch0
 
   LOGICAL, SAVE           ::  initialize = .true.

   ! local arrays
   INTEGER              :: kcorner(3,8) = reshape((/ 0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1, 0,1,1, 1,1,1 /), (/3,8/) )
   COMPLEX              :: exchcorrect(kpts%nkptf)
   COMPLEX              :: dcprod(hybrid%nbands(nk),hybrid%nbands(nk),3) 
   COMPLEX              :: exch_vv(hybrid%nbands(nk),hybrid%nbands(nk))
   COMPLEX              :: hessian(3,3)
   COMPLEX              :: proj_ibsc(3,mnobd,hybrid%nbands(nk))
   COMPLEX              :: olap_ibsc(3,3,mnobd,mnobd)
   REAL                 :: carr1_v_r(hybrid%maxbasm1),carr1_c_r(hybrid%maxbasm1)
   COMPLEX              :: carr1_v_c(hybrid%maxbasm1),carr1_c_c(hybrid%maxbasm1)
   COMPLEX, ALLOCATABLE :: phase_vv(:,:)
   REAL,    ALLOCATABLE :: cprod_vv_r(:,:,:),cprod_cv_r(:,:,:), carr3_vv_r(:,:,:),carr3_cv_r(:,:,:)
   COMPLEX, ALLOCATABLE :: cprod_vv_c(:,:,:),cprod_cv_c(:,:,:), carr3_vv_c(:,:,:),carr3_cv_c(:,:,:)


#if defined(CPP_MPI)&&defined(CPP_NEVER)
   COMPLEX             :: buf_vv(hybrid%nbands(nk),nbands(nk))
#endif

#if ( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
   REAL                 :: coulomb_mt1(hybrid%maxindxm1-1,hybrid%maxindxm1-1, 0:hybrid%maxlcutm1,atoms%ntype)       
   REAL                 :: coulomb_mt2_r(hybrid%maxindxm1-1,-hybrid%maxlcutm1:hybrid%maxlcutm1,0:hybrid%maxlcutm1+1,atoms%nat)
   REAL                 :: coulomb_mt3_r(hybrid%maxindxm1-1,atoms%nat,atoms%nat)
   COMPLEX              :: coulomb_mt2_c(hybrid%maxindxm1-1,-hybrid%maxlcutm1:hybrid%maxlcutm1,0:hybrid%maxlcutm1+1,atoms%nat)
   COMPLEX              :: coulomb_mt3_c(hybrid%maxindxm1-1,atoms%nat,atoms%nat)
#else
   REAL                 :: coulomb_r(hybrid%maxbasm1*(hybrid%maxbasm1+1)/2)
   COMPLEX              :: coulomb_c(hybrid%maxbasm1*(hybrid%maxbasm1+1)/2)
#endif

#ifdef CPP_IRCOULOMBAPPROX
   REAL                 :: coulomb_mtir_r((hybrid%maxlcutm1+1)**2*atoms%nat,&
                                          (hybrid%maxlcutm1+1)**2*atoms%nat+maxval(hybrid%ngptm))
#else
   REAL                 :: coulomb_mtir_r(((hybrid%maxlcutm1+1)**2*atoms%nat +maxval(hybrid%ngptm)) *&
                                          ((hybrid%maxlcutm1+1)**2*atoms%nat +maxval(hybrid%ngptm)+1)/2)
#endif

#ifdef CPP_IRCOULOMBAPPROX
   COMPLEX              :: coulomb_mtir_c((hybrid%maxlcutm1+1)**2*atoms%nat,&
                                          (hybrid%maxlcutm1+1)**2*atoms%nat+maxval(hybrid%ngptm))
#else
   COMPLEX              :: coulomb_mtir_c(((hybrid%maxlcutm1+1)**2*atoms%nat +maxval(hybrid%ngptm)) *&
                                          ((hybrid%maxlcutm1+1)**2*atoms%nat +maxval(hybrid%ngptm)+1)/2)
#endif

   LOGICAL              :: occup(dimension%neigd)
#if defined(CPP_MPI)&&defined(CPP_NEVER)
   INCLUDE "mpif.h"
   INTEGER              :: ierr,ierr2,length,rank
   CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: errmsg
#endif
   CALL timestart("valence exchange calculation")
     
   IF(initialize) THEN !it .eq. 1 .and. nk .eq. 1) THEN
      call calc_divergence(cell,kpts,divergence)
      PRINT *,"Divergence:",divergence
      initialize = .false.
   END IF
   
   ! calculate valence-valence-valence-valence, core-valence-valence-valence
   ! and core-valence-valence-core exchange at current k-point
   ! the sum over the inner occupied valence states is restricted to the EIBZ(k)
   ! the contribution of the Gamma-point is treated separately (see below)

   ! determine package size loop over the occupied bands
   rdum  = hybrid%maxbasm1*hybrid%nbands(nk)*4/1048576.
   psize = 1
   DO iband = mnobd,1,-1
      ! ensure that the packages have equal size
      IF(modulo(mnobd,iband).eq.0) THEN
         ! choose packet size such that cprod is smaller than memory threshold
         IF(rdum*iband.le.maxmem) THEN
            psize = iband
            EXIT
         END IF
      END IF
   END DO

   IF(psize.ne.mnobd) THEN
      WRITE(6,'(A,A,i3,A,f7.2,A)') ' Divide the loop over the occupied hybrid%bands in packages',&
                                   ' of the size',psize,' (cprod=',rdum*psize,'MB)'
   END IF
   ALLOCATE( phase_vv(psize,hybrid%nbands(nk)),stat=ok )
   IF(ok.ne.0) STOP 'exchange_val_hf: error allocation phase'
   phase_vv=0
   IF(ok.ne.0) STOP 'exchange_val_hf: error allocation phase'

   if (mat_ex%l_real) THEN
      ALLOCATE( cprod_vv_c(hybrid%maxbasm1,0,0), carr3_vv_c(hybrid%maxbasm1,0,0))
      ALLOCATE( cprod_vv_r(hybrid%maxbasm1,psize,hybrid%nbands(nk)),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation cprod'
      ALLOCATE( carr3_vv_r(hybrid%maxbasm1,psize,hybrid%nbands(nk)),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation carr3'
      cprod_vv_r = 0 ; carr3_vv_r = 0 
   ELSE
      ALLOCATE( cprod_vv_r(hybrid%maxbasm1,0,0), carr3_vv_r(hybrid%maxbasm1,0,0))
      ALLOCATE( cprod_vv_c(hybrid%maxbasm1,psize,hybrid%nbands(nk)),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation cprod'
      ALLOCATE( carr3_vv_c(hybrid%maxbasm1,psize,hybrid%nbands(nk)),stat=ok )
      IF( ok .ne. 0 ) STOP 'exchange_val_hf: error allocation carr3'
      cprod_vv_c = 0 ; carr3_vv_c = 0
   END IF
         
   exch_vv = 0

   DO ikpt = 1,nkpt_EIBZ

      ikpt0 = pointer_EIBZ(ikpt)

      n  = hybrid%nbasp + hybrid%ngptm(ikpt0)
      IF( hybrid%nbasm(ikpt0).ne.n) STOP 'error hybrid%nbasm'
      nn = n*(n+1)/2

      ! read in coulomb matrix from direct access file coulomb
      IF (mat_ex%l_real) THEN
         CALL read_coulomb_spm_r(kpts%bkp(ikpt0),coulomb_mt1,coulomb_mt2_r,coulomb_mt3_r,coulomb_mtir_r)
      ELSE
         CALL read_coulomb_spm_c(kpts%bkp(ikpt0),coulomb_mt1,coulomb_mt2_c,coulomb_mt3_c,coulomb_mtir_c)
      END IF

      IF(kpts%bkp(ikpt0).ne.ikpt0) THEN
#if( !defined CPP_NOSPMVEC && !defined CPP_IRAPPROX )
         IF((kpts%bksym(ikpt0).gt.sym%nop).and.(.not.mat_ex%l_real)) THEN
            coulomb_mt2_c = conjg(coulomb_mt2_c)
            coulomb_mtir_c = conjg(coulomb_mtir_c)
         END IF
#else
         if (.not.mat_ex%l_real) THEN
            IF( kpts%bksym(ikpt0) .gt. sym%nop ) coulomb = conjg(coulomb)
         endif
#endif
      END IF

      DO ibando = 1, mnobd, psize

         IF (mat_ex%l_real) THEN
#ifdef CPP_IRAPPROX
            CALL wavefproducts_inv(1,hybdat,dimension,input,jsp,atoms,lapw,obsolete,kpts,nk,ikpt0,&
                                   mnobd,hybrid,parent,cell,sym,noco,nkqpt,cprod_vv)
#else
            CALL wavefproducts_inv5(1,hybrid%nbands(nk),ibando,ibando+psize-1,dimension,input,jsp,atoms,&
                                    lapw,kpts,nk,ikpt0,hybdat,mnobd,hybrid,parent,cell,hybrid%nbasp,sym,&
                                    noco,nkqpt,cprod_vv_r)
#endif
         ELSE
#ifdef CPP_IRAPPROX
            CALL wavefproducts_noinv(1,hybdat,nk,ikpt0,dimension,input,jsp,cell,atoms,hybrid, 
                                     kpts,mnobd,lapw,sym,noco,nkqpt,cprod_vv)
#else
            CALL wavefproducts_noinv5(1,hybrid%nbands(nk),ibando,ibando+psize-1,nk,ikpt0,dimension,input,jsp,&!jsp,&
                                      cell,atoms,hybrid,hybdat,kpts,mnobd,lapw,sym,hybrid%nbasp,noco,nkqpt,cprod_vv_c)
#endif
         END IF

         ! The sparse matrix technique is not feasible for the HSE
         ! functional. Thus, a dynamic adjustment is implemented
         ! The mixed basis functions and the potential difference
         ! are Fourier transformed, so that the exchange can be calculated
         ! in Fourier space
#ifndef CPP_NOSPMVEC
         IF (xcpot%is_name("hse").OR.xcpot%is_name("vhse")) THEN
            iband1  = hybrid%nobd(nkqpt)

            exch_vv = exch_vv +&
                      dynamic_hse_adjustment(atoms%rmsh,atoms%rmt,atoms%dx,atoms%jri,atoms%jmtd,kpts%bkf(:,ikpt0),ikpt0,&
                                             kpts%nkptf,cell%bmat,cell%omtil,atoms%ntype,atoms%neq,atoms%nat,atoms%taual,&
                                             hybrid%lcutm1,hybrid%maxlcutm1,hybrid%nindxm1,hybrid%maxindxm1,hybrid%gptm,&
                                             hybrid%ngptm(ikpt0),hybrid%pgptm(:,ikpt0),hybrid%gptmd,hybrid%basm1,&
                                             hybrid%nbasm(ikpt0),iband1,hybrid%nbands(nk),nsest,ibando,psize,indx_sest,&
                                             atoms%invsat,sym%invsatnr,mpi%irank,cprod_vv_r(:hybrid%nbasm(ikpt0),:,:),&
                                             cprod_vv_c(:hybrid%nbasm(ikpt0),:,:),mat_ex%l_real,wl_iks(:iband1,nkqpt),n_q(ikpt))
         END IF
#endif

         ! the Coulomb matrix is only evaluated at the irrecuible k-points
         ! bra_trafo transforms cprod instead of rotating the Coulomb matrix
         ! from IBZ to current k-point
         IF( kpts%bkp(ikpt0) .ne. ikpt0 ) THEN
            CALL bra_trafo2(mat_ex%l_real,carr3_vv_r(:hybrid%nbasm(ikpt0),:,:),cprod_vv_r(:hybrid%nbasm(ikpt0),:,:),&
                            carr3_vv_c(:hybrid%nbasm(ikpt0),:,:),cprod_vv_c(:hybrid%nbasm(ikpt0),:,:),&
                            hybrid%nbasm(ikpt0),psize,hybrid%nbands(nk),kpts%bkp(ikpt0),ikpt0,kpts%bksym(ikpt0),sym,&
                            hybrid,kpts,cell,atoms,phase_vv)
            IF (mat_ex%l_real) THEN
               cprod_vv_r(:hybrid%nbasm(ikpt0),:,:) = carr3_vv_r(:hybrid%nbasm(ikpt0),:,:)
            ELSE
               cprod_vv_c(:hybrid%nbasm(ikpt0),:,:) = carr3_vv_c(:hybrid%nbasm(ikpt0),:,:)
            ENDIF
         ELSE
            phase_vv(:,:) = (1d0,0d0)
         END IF

         ! calculate exchange matrix at ikpt0
   
         call timestart("exchange matrix")
         DO n1=1,hybrid%nbands(nk)
            DO iband = 1,psize
               IF((ibando+iband-1).gt.hybrid%nobd(nkqpt)) CYCLE

               cdum  = wl_iks(ibando+iband-1,nkqpt) * conjg(phase_vv(iband,n1))/n_q(ikpt)

               IF (mat_ex%l_real) THEN
                  carr1_v_r(:n) = 0 
                  CALL spmvec_invs(atoms,hybrid,hybdat,ikpt0,kpts,cell,coulomb_mt1,coulomb_mt2_r,coulomb_mt3_r,&
                                   coulomb_mtir_r,cprod_vv_r(:n,iband,n1),carr1_v_r(:n))
               ELSE
                  carr1_v_c(:n) = 0 
                  CALL spmvec_noinvs(atoms,hybrid,hybdat,ikpt0,kpts,cell,coulomb_mt1,coulomb_mt2_c,coulomb_mt3_c,&
                                     coulomb_mtir_c,cprod_vv_c(:n,iband,n1),carr1_v_c(:n))
               END IF

               IF (mat_ex%l_real) THEN
                  DO n2=1,nsest(n1)!n1
                     nn2 = indx_sest(n2,n1)
                     exch_vv(nn2,n1) = exch_vv(nn2,n1) + cdum*phase_vv(iband,nn2) *&
                                                         dotprod(carr1_v_r(:n),cprod_vv_r(:n,iband,nn2))
                  END DO !n2
               ELSE
                  DO n2=1,nsest(n1)!n1
                     nn2 = indx_sest(n2,n1)
                     exch_vv(nn2,n1) = exch_vv(nn2,n1) + cdum*phase_vv(iband,nn2) *&
                                                         dotprod(carr1_v_c(:n),cprod_vv_c(:n,iband,nn2))
                  END DO !n2
               END IF
            END DO
         END DO  !n1
         call timestop("exchange matrix")
      END DO !ibando
   END DO  !ikpt

!   WRITE(7001,'(a,i7)') 'nk: ', nk
!   DO n1=1,hybrid%nbands(nk)
!      DO n2=1,n1
!         WRITE(7001,'(2i7,2f15.8)') n2, n1, exch_vv(n2,n1)
!     END DO
!   END DO

   ! add contribution of the gamma point to the different cases (exch_vv,exch_cv,exch_cc)

   ! valence-valence-valence-valence exchange

   IF ((.not.xcpot%is_name("hse")).AND.(.not.xcpot%is_name("vhse"))) THEN ! no gamma point correction needed for HSE functional
      IF( zero_order .and. .not. ibs_corr ) THEN
         WRITE(6,'(A)') ' Take zero order terms into account.'
      ELSE IF( zero_order .and.  ibs_corr ) THEN
         WRITE(6,'(A)') ' Take zero order terms and ibs-correction into account.'
      END IF

      IF(zero_order) THEN
         CALL dwavefproducts(dcprod,nk,1,hybrid%nbands(nk),1,hybrid%nbands(nk),.false.,atoms,hybrid,&
                             cell,hybdat,kpts,kpts%nkpt,lapw,dimension,jsp,eig_irr)

         ! make dcprod hermitian
         DO n1 = 1, hybrid%nbands(nk)
            DO n2 = 1,n1
               dcprod(n1,n2,:) = (dcprod(n1,n2,:) - conjg(dcprod(n2,n1,:)))/2   
               dcprod(n2,n1,:) = -conjg(dcprod(n1,n2,:))
            END DO
         END DO

         IF(ibs_corr) THEN
            CALL ibs_correction(nk,atoms,dimension,input,jsp,hybdat,hybrid,lapw,kpts,kpts%nkpt,cell,mnobd,&
                                sym,proj_ibsc,olap_ibsc)
         END IF
      END IF
        
      !This should be done with w_iks I guess!TODO
      occup = .false.
      DO i=1,hybrid%ne_eig(nk)
         IF (results%ef.ge.eig_irr(i,nk)) THEN
            occup(i) = .true.
         ELSE IF ((eig_irr(i,nk)-results%ef).le.1E-06) THEN
            occup(i) = .true.
         END IF
      END DO

      DO n1 = 1, hybrid%nbands(nk)
         DO n2 = 1, nsest(n1)!n1
            nn2 = indx_sest(n2,n1)
            exchcorrect = 0
            exch0 = 0

            ! if zero_order = .true. add averaged k-dependent term to the numerical integration at Gamma-point contribution

            ! if we start with a system with a small DFT band gap (like GaAs), the contribution
            ! of the highest occupied and lowest unoccupied state in Hessian is typically
            ! large; a correct numerical integration requires a dense k-point mesh, so
            ! we don't add the contribution exchcorrect for such materials 

            IF(zero_order) THEN
               hessian = 0
               IF(occup(n1).and.occup(nn2)) THEN
                  DO i = 1,3
                     j = i
                     DO iband = 1, hybrid%nbands(nk)
                        IF(occup(iband)) THEN
                           hessian(i,j) = hessian(i,j) + conjg(dcprod(iband,n1,i)) *dcprod(iband,nn2,j)
                        END IF
                        hessian(i,j) = hessian(i,j) - dcprod(iband,nn2,i) * conjg(dcprod(iband,n1,j))
                     END DO

                     ! ibs correction
                     IF(ibs_corr) THEN 
                        hessian(i,j) = hessian(i,j) - olap_ibsc(i,j,n1,nn2)/cell%omtil
                        DO iband = 1,hybrid%nbands(nk)
                           hessian(i,j) = hessian(i,j) + conjg(proj_ibsc(i,nn2,iband)) * proj_ibsc(j,n1,iband)/cell%omtil
                        END DO
                     END IF
                  END DO
               ELSE
                  DO i = 1,3
                     j = i 
                     DO iband = 1, hybrid%nbands(nk)
                        IF(occup(iband)) THEN
                           hessian(i,j) = hessian(i,j) + conjg(dcprod(iband,n1,i)) * dcprod(iband,nn2,j)
                        END IF
                     END DO
                  END DO
               END IF
 
               exchcorrect(1) = fpi_const/3 * (hessian(1,1)+hessian(2,2)+hessian(3,3))
               exch0 = exchcorrect(1)/kpts%nkptf
            END IF

            ! tail correction/contribution from all other k-points (it  goes into exchcorrect )

            ! Analytic contribution

            cdum2 = 0
            !multiply divergent contribution with occupation number;
            !this only affects metals 
            IF (n1.eq.nn2) THEN
               cdum2 = fpi_const/cell%omtil * divergence * wl_iks(n1,nk)*kpts%nkptf
            END IF

            ! due to the symmetrization afterwards the factor 1/n_q(1) must be added

            IF(n1.EQ.nn2) hybrid%div_vv(n1,nk,jsp) = REAL(cdum2) 
            exch_vv(nn2,n1)  = exch_vv(nn2,n1) + (exch0 + cdum2)/n_q(1)

         END DO !n2
      END DO !n1
   END IF ! xcpot%icorr .ne. icorr_hse


   IF (mat_ex%l_real) THEN
      IF(any(abs(aimag(exch_vv)).gt.1E-08)) CALL judft_warn('unusally large imaginary part of exch_vv',&
                                                            calledby='exchange_val_hf.F90')
   END IF

!   WRITE(7000,'(a,i7)') 'nk: ', nk
!   DO n1=1,hybrid%nbands(nk)
!      DO n2=1,n1
!         WRITE(7000,'(2i7,2f15.8)') n2, n1, exch_vv(n2,n1)
!      END DO
!   END DO

   ! write exch_vv in mat_ex
   CALL mat_ex%alloc(matsize1=hybrid%nbands(nk))
   IF (mat_ex%l_real) THEN
      mat_ex%data_r=exch_vv
   ELSE
      mat_ex%data_c=exch_vv
   END IF
   CALL timestop("valence exchange calculation")
     
END SUBROUTINE exchange_valence_hf




SUBROUTINE calc_divergence(cell,kpts,divergence)

   USE m_util, ONLY: cerf
   USE m_types
   USE m_constants

   IMPLICIT NONE

   TYPE(t_cell), INTENT(IN)  :: cell
   TYPE(t_kpts), INTENT(IN)  :: kpts
   REAL,         INTENT(OUT) :: divergence
        
   INTEGER :: ix,iy,iz,sign,n
   logical :: found
   REAL    :: expo,rrad,k(3),kv1(3),kv2(3),kv3(3),knorm2
   COMPLEX :: cdum
        
   expo       = 5d-3
   rrad       = sqrt(-log(5d-3)/expo)
   cdum       = sqrt(expo)*rrad
   divergence = cell%omtil / (tpi_const**2) * sqrt(pi_const/expo) * cerf(cdum)
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
            DO sign=-1,1,2
               iz=sign*iz
               k(1) = ix*kv1(1) + iy*kv2(1) + iz*kv3(1)
               k(2) = ix*kv1(2) + iy*kv2(2) + iz*kv3(2)
               k(3) = ix*kv1(3) + iy*kv2(3) + iz*kv3(3)
               knorm2 = k(1)**2 + k(2)**2 + k(3)**2
               IF(knorm2.lt.rrad) THEN
                  found = .true.
                  divergence = divergence - exp(-expo*knorm2)/knorm2 / kpts%nkptf
               END IF
               IF(iz==0) exit
            END DO 
         END DO
      END DO
      n = n + 1
   END DO

END SUBROUTINE calc_divergence

END MODULE m_exchange_valence_hf
