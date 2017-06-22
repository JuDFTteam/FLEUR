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
!     for the different combinations of n_1 and n_2 and  n' runs over the core states.
!     ( n_1,n_2:  valence-valence, core-core,core-valence )
!
!     It is done directly without employing the mixed basis set.


      MODULE m_exchange_core

      CONTAINS

      SUBROUTINE exchange_vccv(&
     &                   nk,bkpt,kpts,nkpti,atoms,&
     &                   hybrid,hybdat,&
     &                   dimension,jsp,lapw,&
     &                   maxbands,mnobd,mpi,irank2,&
     &                   degenerat,symequivalent,results,&
     &                   ex_vv)


      USE m_constants   
      USE m_util
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      TYPE(t_hybdat),INTENT(IN)   :: hybdat
      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)   :: lapw

!     -scalars -
      INTEGER,INTENT(IN)      :: jsp 
      INTEGER,INTENT(IN)      ::nk  ,maxbands, mnobd
      INTEGER,INTENT(IN)      :: nkpti ,irank2
!     - arays -
      INTEGER,INTENT(IN)      ::  degenerat(hybdat%ne_eig(nk))
      REAL,INTENT(IN)         ::  bkpt(3)

#ifdef CPP_INVERSION
      REAL    ,INTENT(INOUT)  ::  ex_vv(maxbands,mnobd,nkpti)
#else
      COMPLEX ,INTENT(INOUT)  ::  ex_vv(maxbands,mnobd,nkpti)
#endif
      LOGICAL                 ::  symequivalent(count(degenerat .ge. 1),&
     &                                          count(degenerat .ge. 1))

!     - local scalars -
      INTEGER                 ::  iatom,ieq,itype,ic,l,l1,l2,&
     &                            ll,lm ,m1,m2,p1,p2,n,n1,n2,i,j
      INTEGER                 ::  iband1,iband2,ndb1,ndb2,ic1,ic2
      INTEGER                 ::  irecl_cmt,m

      REAL                    ::  time1,time2
      REAL                    ::  rdum
      REAL                    ::  sum_offdia

      COMPLEX                 ::  cdum

!     - local arrays -
      INTEGER,ALLOCATABLE     ::  larr(:),larr2(:)
      INTEGER,ALLOCATABLE     ::  parr(:),parr2(:)

      REAL                    ::  integrand(atoms%jmtd)
      REAL                    ::  primf1(atoms%jmtd),primf2(atoms%jmtd)
      REAL,ALLOCATABLE        ::  fprod(:,:),fprod2(:,:)
      REAL,ALLOCATABLE        ::  integral(:,:)

      COMPLEX                 ::  cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat)
      COMPLEX                 ::  exchange(hybdat%nbands(nk),hybdat%nbands(nk))
      COMPLEX,ALLOCATABLE     ::  carr(:,:),carr2(:,:),carr3(:,:)

      LOGICAL                 ::  ldum(hybdat%nbands(nk),hybdat%nbands(nk))

      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A)') new_line('n') // new_line('n') // '### valence-core-core-valence exchange ###'
        WRITE(6,'(A)') new_line('n') // '        k-point       band    exchange (core contribution)'
      END IF

      ! read in mt wavefunction coefficients from file cmt
      irecl_cmt = dimension%neigd*hybrid%maxlmindx*atoms%nat*16
      OPEN(unit=777,file='cmt',form='unformatted',access='direct', recl=irecl_cmt)
      READ(777,rec=nk) cmt(:,:,:)
      CLOSE(777)

      ALLOCATE ( fprod(atoms%jmtd,5),larr(5),parr(5) )

       ! generate ldum(nbands(nk),nbands(nk)), which is true if the corresponding matrix entry is non-zero
      ic1  = 0
      ldum = .false.
      DO iband1 = 1,hybdat%nbands(nk)
        ndb1 = degenerat(iband1)
        IF( ndb1 .ge. 1 ) THEN
          ic1 = ic1 + 1
          ic2 = 0
          DO iband2 = 1,hybdat%nbands(nk)
            ndb2 = degenerat(iband2)
            IF( ndb2 .ge. 1 ) THEN
              ic2 = ic2 + 1
              IF( symequivalent(ic2,ic1) ) THEN
                IF( ndb1 .ne. ndb2 ) STOP 'exchange: failure symequivalent'
                DO i = 0,ndb1-1
                  DO j = 0,ndb2 - 1
                    ldum(iband1+i,iband2+j) = .true.
                  END DO
                END DO

              END IF
            END IF
          END DO
        END IF
      END DO

      exchange = 0
      iatom = 0
      rdum  = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l1 = 0,hybdat%lmaxc(itype)
            DO p1 = 1,hybdat%nindxc(l1,itype)

              DO l = 0,hybrid%lcutm1(itype)

              ! Define core-valence product functions

                n = 0
                DO l2 = 0,atoms%lmax(itype)
                  IF(l.lt.abs(l1-l2).or.l.gt.l1+l2) CYCLE

                  DO p2 = 1,hybrid%nindx(l2,itype)
                    n = n + 1
                    M = size(fprod,2)
                    IF(n.gt.M) THEN
                      ALLOCATE ( fprod2(atoms%jmtd,M),larr2(M),parr2(M) )
                      fprod2 = fprod ; larr2 = larr ; parr2 = parr
                      DEALLOCATE ( fprod,larr,parr )
                      ALLOCATE ( fprod(atoms%jmtd,M+5),larr(M+5),parr(M+5) )
                      fprod(:,:M) = fprod2
                      larr(:M)    = larr2
                      parr(:M)    = parr2
                      DEALLOCATE ( fprod2,larr2,parr2 )
                    END IF
                    fprod(:,n) = ( hybdat%core1(:,p1,l1,itype) *hybdat%bas1 (:,p2,l2,itype)&
     &                            +hybdat%core2(:,p1,l1,itype) *hybdat%bas2 (:,p2,l2,itype) )/ atoms%rmsh(:,itype)
                    larr(n)    = l2
                    parr(n)    = p2
                  END DO
                END DO

                ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)

                ALLOCATE ( integral(n,n),carr(n,hybdat%nbands(nk)), carr2(n,lapw%nv(jsp)),carr3(n,lapw%nv(jsp)) )

                DO i = 1,n
                  CALL primitivef(primf1,fprod(:,i)*atoms%rmsh(:,itype)**(l+1) ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd, itype,atoms%ntype)
                  CALL primitivef(primf2,fprod(:,i)/atoms%rmsh(:,itype)**l ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,-itype,atoms%ntype)  ! -itype is to enforce inward integration

                  primf1 = primf1 / atoms%rmsh(:,itype)**l
                  primf2 = primf2 * atoms%rmsh(:,itype)**(l+1)
                  DO j = 1,n
                    integrand     = fprod(:,j) * (primf1 + primf2)
                    integral(i,j) = fpi_const/(2*l+1) * intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,&
     &                                     atoms%dx,atoms%ntype,itype,hybdat%gridf)
                  END DO
                END DO

                ! Add everything up

                DO m1 = -l1,l1
                  DO M = -l,l
                    m2 = m1 + M

                    carr = 0
                    ic   = 0
                    DO n1=1,hybdat%nbands(nk)

                      DO i = 1,n
                        ll      = larr(i)
                        IF(abs(m2).gt.ll) CYCLE

                        lm = sum((/ ((2*l2+1)*hybrid%nindx(l2,itype),l2=0,ll-1) /)) + (m2+ll)*hybrid%nindx(ll,itype) + parr(i)

                        carr(i,n1) = cmt(n1,lm,iatom) * gaunt(l1,ll,l,m1,m2,M,hybdat%maxfac,hybdat%fac,hybdat%sfac)

                      END DO
                      DO n2=1,n1
                        IF( ldum(n2,n1) ) THEN
                          ic =ic + 1
                          exchange(n2,n1) = exchange(n2,n1) + dot_product( carr(:,n1), matmul(integral,carr(:,n2)) )
                        END IF
                      END DO
                    END DO
                  END DO
                END DO

                DEALLOCATE ( integral,carr,carr2,carr3 )

              END DO
            END DO
          END DO
        END DO
      END DO



#ifdef CPP_INVERSION
      IF( any(abs(aimag(exchange)).gt.1d-10) ) THEN
        IF ( mpi%irank == 0 ) WRITE(6,'(A)') 'exchangeCore: Warning! Unusually large imaginary component.'
        WRITE(*,*) maxval(abs(aimag(exchange)))
        STOP 'exchangeCore: Unusually large imaginary component.'
      END IF
#endif

      ! add the core-valence contribution to the exchange matrix ex_vv

!      ic         = 0
      sum_offdia = 0
      DO n1=1,hybdat%nobd(nk)
        DO n2=1,hybdat%nbands(nk)
          ex_vv(n2,n1,nk) = ex_vv(n2,n1,nk) - exchange(n1,n2)
          IF( n1 /= n2) sum_offdia = sum_offdia + 2*abs(exchange(n1,n2))

        END DO
      END DO


      DO n1=1,hybdat%nobd(nk)
        results%te_hfex%core = results%te_hfex%core - results%w_iks(n1,nk,jsp)*exchange(n1,n1)
      END DO

      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A,F20.15)') 'sum of the absolut real part of the non diagonal elements',sum_offdia
      END IF


      END SUBROUTINE exchange_vccv

      SUBROUTINE exchange_vccv1(nk,kpts,nkpti,atoms,&
     &                          hybrid,hybdat,&
     &                          dimension,jsp,&
     &                          lapw,&
     &                          nsymop,nsest,indx_sest,mpi,&
     &                          a_ex,results,&
     &                          mat_ex)
     
      USE m_constants  
      USE m_util
      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      TYPE(t_hybdat),INTENT(IN)   :: hybdat
      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)   :: lapw
      
!     -scalars -
      INTEGER,INTENT(IN)      :: jsp 
      INTEGER,INTENT(IN)      :: nk  
      INTEGER,INTENT(IN)      :: nkpti 
      INTEGER,INTENT(IN)      ::  nsymop
      REAL,INTENT(IN)         ::  a_ex
!     - arays -
      INTEGER,INTENT(IN)      ::  nsest(hybdat%nbands(nk)),indx_sest(hybdat%nbands(nk),hybdat%nbands(nk))
      

#ifdef CPP_INVERSION
      REAL    ,INTENT(INOUT)  ::  mat_ex(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#else
      COMPLEX ,INTENT(INOUT)  ::  mat_ex(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#endif
!     - local scalars - 
      INTEGER                 ::  iatom,ieq,itype,ic,l,l1,l2,&
     &                            ll,lm ,m1,m2,p1,p2,n,n1,n2,nn2,i,j
      INTEGER                 ::  iband1,iband2,ndb1,ndb2,ic1,ic2
      INTEGER                 ::  irecl_cmt,m
      
      REAL                    ::  time1,time2
      REAL                    ::  rdum
      REAL                    ::  sum_offdia

      COMPLEX                 ::  cdum
!     - local arrays -
      INTEGER,ALLOCATABLE     ::  larr(:),larr2(:)
      INTEGER,ALLOCATABLE     ::  parr(:),parr2(:)
      
      REAL                    ::  integrand(atoms%jmtd)
      REAL                    ::  primf1(atoms%jmtd),primf2(atoms%jmtd)
      REAL,ALLOCATABLE        ::  fprod(:,:),fprod2(:,:)
      REAL,ALLOCATABLE        ::  integral(:,:)
     
      COMPLEX                 ::  cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat)
      COMPLEX                 ::  exchange(hybdat%nbands(nk),hybdat%nbands(nk))
      COMPLEX,ALLOCATABLE     ::  carr(:,:),carr2(:,:),carr3(:,:)
      
      LOGICAL                 ::  ldum(hybdat%nbands(nk),hybdat%nbands(nk))
      
      

      ! read in mt wavefunction coefficients from file cmt
      irecl_cmt = dimension%neigd*hybrid%maxlmindx*atoms%nat*16
      OPEN(unit=777,file='cmt',form='unformatted',access='direct', recl=irecl_cmt)
      READ(777,rec=nk) cmt(:,:,:)
      CLOSE(777)

      ALLOCATE ( fprod(atoms%jmtd,5),larr(5),parr(5) )
      
      exchange = 0
      iatom    = 0
      rdum     = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l1 = 0,hybdat%lmaxc(itype)
            DO p1 = 1,hybdat%nindxc(l1,itype)

              DO l = 0,hybrid%lcutm1(itype) 
              
              ! Define core-valence product functions

                n = 0
                DO l2 = 0,atoms%lmax(itype)
                  IF(l.lt.abs(l1-l2).or.l.gt.l1+l2) CYCLE

                  DO p2 = 1,hybrid%nindx(l2,itype)
                    n = n + 1
                    M = size(fprod,2)
                    IF(n.gt.M) THEN
                      ALLOCATE ( fprod2(atoms%jmtd,M),larr2(M),parr2(M) )
                      fprod2 = fprod ; larr2 = larr ; parr2 = parr
                      DEALLOCATE ( fprod,larr,parr )
                      ALLOCATE ( fprod(atoms%jmtd,M+5),larr(M+5),parr(M+5) )
                      fprod(:,:M) = fprod2
                      larr(:M)    = larr2
                      parr(:M)    = parr2
                      DEALLOCATE ( fprod2,larr2,parr2 )
                    END IF
                    fprod(:,n) = ( hybdat%core1(:,p1,l1,itype) *hybdat%bas1 (:,p2,l2,itype) &
     &                            +hybdat%core2(:,p1,l1,itype) *hybdat%bas2 (:,p2,l2,itype) )/ atoms%rmsh(:,itype)
                    larr(n)    = l2
                    parr(n)    = p2
                  END DO
                END DO

                ! Evaluate radial integrals (special part of Coulomb matrix : contribution from single MT)

                ALLOCATE ( integral(n,n),carr(n,hybdat%nbands(nk)), carr2(n,lapw%nv(jsp)),carr3(n,lapw%nv(jsp)) )

                DO i = 1,n
                   CALL primitivef(primf1,fprod(:,i)*atoms%rmsh(:,itype)**(l+1) ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd, itype,atoms%ntype)
                    CALL primitivef(primf2,fprod(:,i)/atoms%rmsh(:,itype)**l ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,-itype,atoms%ntype)  ! -itype is to enforce inward integration
                  
                  primf1 = primf1 / atoms%rmsh(:,itype)**l
                  primf2 = primf2 * atoms%rmsh(:,itype)**(l+1)
                  DO j = 1,n
                    integrand     = fprod(:,j) * (primf1 + primf2)
                    integral(i,j) = fpi_const/(2*l+1) * intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,&
     &                                     atoms%dx,atoms%ntype,itype,hybdat%gridf)
                  END DO
                END DO

                ! Add everything up

                DO m1 = -l1,l1
                  DO M = -l,l
                    m2 = m1 + M

                    carr = 0
                    DO n1=1,hybdat%nbands(nk)

                      DO i = 1,n
                        ll      = larr(i)
                        IF(abs(m2).gt.ll) CYCLE
                        
                        lm = sum((/ ((2*l2+1)*hybrid%nindx(l2,itype),l2=0,ll-1) /))&
                             + (m2+ll)*hybrid%nindx(ll,itype) + parr(i)
                        
                        carr(i,n1) = cmt(n1,lm,iatom) * gaunt(l1,ll,l,m1,m2,M,hybdat%maxfac,hybdat%fac,hybdat%sfac)
                      
                      END DO
                      DO n2=1,nsest(n1)!n1
                        nn2 = indx_sest(n2,n1)
                        exchange(nn2,n1) = exchange(nn2,n1) + dot_product( carr(:,n1), matmul(integral,carr(:,nn2)) )

                      END DO
                    END DO
                  END DO
                END DO

                DEALLOCATE ( integral,carr,carr2,carr3 )

              END DO
            END DO
          END DO
        END DO
      END DO




#ifdef CPP_INVERSION
      IF( any(abs(aimag(exchange)).gt.1d-10) ) THEN
        IF (mpi%irank == 0) WRITE(6,'(A)') 'exchangeCore: Warning! Unusually large imaginary component.'
        WRITE(*,*) maxval(abs(aimag(exchange))) STOP 'exchangeCore: Unusually large imaginary component.'
      END IF
#endif
      
      DO n1=1,hybdat%nobd(nk)
        results%te_hfex%core = results%te_hfex%Core - a_ex * results%w_iks(n1,nk,jsp)*exchange(n1,n1)
      END DO

      ! add the core-valence contribution to the exchange matrix mat_ex
      ! factor 1/nsymop is needed due to the symmetrization in symmetrizeh 
      
      ic         = 0
      sum_offdia = 0
      DO n1=1,hybdat%nbands(nk)
        DO n2=1,n1
          ic = ic + 1
          mat_ex(ic) = mat_ex(ic) + conjg(exchange(n2,n1))/nsymop
        END DO
      END DO


      END SUBROUTINE exchange_vccv1

                  
      SUBROUTINE exchange_cccc(nk,nkpti,atoms,hybdat, ncstd,&
           sym,kpts,a_ex,mpi, results)
     
     
      USE m_constants   
      USE m_util
      USE m_wrapper
      USE m_gaunt
      USE m_trafo
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybdat),INTENT(IN)   :: hybdat
      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      
      ! - scalars - 
      INTEGER,INTENT(IN)    ::  nk,nkpti   ,ncstd

      REAL   ,INTENT(IN)    ::  a_ex

      ! - arays -
                            
      ! - local scalars - 
      INTEGER               ::  itype,ieq,icst,icst1,icst2,iatom,iatom0
      INTEGER               ::  l1,l2,l,ll,llmax
      INTEGER               ::  m1,m2 ,mm,m
      INTEGER               ::  n1,n2,n
      
      REAL                  ::  rdum,rdum1
      ! - local arrays -
      INTEGER               ::  point(hybdat%maxindxc,-hybdat%lmaxcd:hybdat%lmaxcd, 0:hybdat%lmaxcd,atoms%nat)
      REAL                  ::  rprod(atoms%jmtd),primf1(atoms%jmtd),primf2(atoms%jmtd), integrand(atoms%jmtd)
      COMPLEX               ::  exch(ncstd,ncstd)

!       IF ( irank == 0 ) THEN
!         WRITE(6,'(//A)') '### core-core-core-core exchange ###'
!         WRITE(6,'(/A)') '        k-point       band    exchange'
!       END IF

            
      ! set up point
      icst = 0
      iatom= 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l = 0,hybdat%lmaxc(itype)
            DO M = -l,l
              DO n = 1,hybdat%nindxc(l,itype)
                icst = icst + 1
                point(n,M,l,iatom) = icst
              END DO
            END DO
          END DO
        END DO
      END DO
      
      llmax  = 2*hybdat%lmaxcd 
      exch   = 0
      iatom0 = 0
      DO itype = 1,atoms%ntype
                
        DO l1 = 0,hybdat%lmaxc(itype)  ! left core state
          DO l2 = 0,hybdat%lmaxc(itype)  ! right core state
            DO l = 0,hybdat%lmaxc(itype)   ! occupied core state
            
              DO ll = abs(l1-l),l1+l
                IF( ll .lt. abs(l-l2) .or. ll .gt. l+l2 ) CYCLE
                IF( mod(l+l1+ll,2) .ne. 0 ) CYCLE
                IF( mod(l+l2+ll,2) .ne. 0 ) CYCLE
                
                DO m1 = -l1,l1
                  m2 = m1
                  IF( abs(m2) .gt. l2 ) CYCLE
                  DO M = -l,l
                    mm = M - m1
                    IF( abs(mm) .gt. ll ) CYCLE
                    rdum = fpi_const/(2*ll+1)*gaunt1(l,ll,l1,M,mm,m1,llmax) *gaunt1(l,ll,l2,M,mm,m2,llmax)
                      
                    DO n = 1,hybdat%nindxc(l,itype)
                      DO n2 = 1,hybdat%nindxc(l2,itype)
                        rprod(:) = ( hybdat%core1(:,n,l,itype)*hybdat%core1(:,n2,l2,itype)&
                             +hybdat%core2(:,n,l,itype)*hybdat%core2(:,n2,l2,itype) ) / atoms%rmsh(:,itype)

                        CALL primitivef(primf1,rprod(:)*atoms%rmsh(:,itype)**(ll+1) ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd, itype,atoms%ntype)
                        CALL primitivef(primf2,rprod(:)/atoms%rmsh(:,itype)**ll ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,-itype,atoms%ntype)  ! -itype is to enforce inward integration
                  

                        primf1 = primf1 / atoms%rmsh(:,itype)**ll
                        primf2 = primf2 * atoms%rmsh(:,itype)**(ll+1)
                
                        DO n1 = 1,hybdat%nindxc(l1,itype)
                            
                          rprod(:) = ( hybdat%core1(:,n,l,itype)*hybdat%core1(:,n1,l1,itype)&
                               +hybdat%core2(:,n,l,itype)*hybdat%core2(:,n1,l1,itype) ) / atoms%rmsh(:,itype)
      
                          integrand     = rprod * (primf1 + primf2)
                          
                          rdum1 = rdum*intgrf(integrand,atoms%jri,atoms%jmtd,&
                                        atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)
                          
                          iatom = iatom0
                          DO ieq = 1,atoms%neq(itype)
                            iatom = iatom + 1
                            icst1    = point(n1,m1,l1,iatom)
                            icst2    = point(n2,m2,l2,iatom)
                            exch(icst1,icst2) = exch(icst1,icst2)+rdum1
                          END DO
                        END DO  !n1
                        
                      END DO  !n2
                    END DO  !n
                    
                  END DO  !M
                END DO  !m1
                
              END DO  !ll
              
            END DO  !l
          END DO  !l2
        END DO  !l1
        iatom0 = iatom0 + atoms%neq(itype)
      END DO  !itype      
       
#ifdef CPP_INVERSION
      CALL symmetrize(exch,ncstd,ncstd,3,.false.,&
                     atoms,ntype,hybdat,&
                     sym)
       IF( any( abs(aimag(exch)) .gt. 1E-6 ) ) STOP 'exchange_cccc: exch possesses significant imaginary part'
#endif
!       DO icst = 1,ncstd
!         IF ( irank == 0 )
!      &    WRITE(6,'(    ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,1X,F12.5)')bkpt,icst,REAL(exch(icst,icst))*(-27.211608)
!       END DO


      ! add core exchange contributions to the te_hfex
      
      DO icst1 = 1,ncstd
        results%te_hfex%core = results%te_hfex%core - a_ex*kpts%wtkpt(nk)*exch(icst1,icst1)
      END DO
      
      END SUBROUTINE exchange_cccc
      
      SUBROUTINE exchange_cccv( &
     &        nk,nkpti,atoms,hybdat,&
     &        hybrid,dimension,maxbands,ncstd,&
     &        bkpt,sym,mpi,&
     &        exch_cv )
     
      USE m_constants   
      USE m_util
      USE m_wrapper
      USE m_gaunt
      USE m_trafo
      
    USE m_types
      IMPLICIT NONE
      TYPE(t_hybdat),INTENT(IN)   :: hybdat
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_atoms),INTENT(IN)   :: atoms
      ! - scalars - 
      INTEGER,INTENT(IN)    ::  nk,nkpti  ,ncstd
      INTEGER,INTENT(IN)    :: maxbands

      ! - arays -
      REAL   ,INTENT(IN)    ::  bkpt(3)
#ifdef CPP_INVERSION
      REAL   ,INTENT(INOUT) ::  exch_cv(maxbands,ncstd,nkpti)
#else
      COMPLEX,INTENT(INOUT) ::  exch_cv(maxbands,ncstd,nkpti)
#endif
      ! - local scalars - 
      INTEGER               ::  itype,ieq,icst,icst1,icst2,iatom,iatom0,&
     &                          iatom1,iband
      INTEGER               ::  l1,l2,l,ll,llmax
      INTEGER               ::  lm2,lmp2
      INTEGER               ::  m1,m2 ,mm,m
      INTEGER               ::  n1,n2,n,nn
      INTEGER               ::  irecl_cmt
      
      REAL                  ::  rdum0,rdum1,rdum2,rdum3,rdum4
      COMPLEX               ::  cdum
      COMPLEX,PARAMETER     ::  img=(0d0,1d0)
      ! - local arrays -
      INTEGER               ::  point(hybdat%maxindxc,-hybdat%lmaxcd:hybdat%lmaxcd, 0:hybdat%lmaxcd,atoms%nat)
      INTEGER               ::  lmstart(0:atoms%lmaxd,atoms%ntype)
      REAL                  ::  rprod(atoms%jmtd),primf1(atoms%jmtd),primf2(atoms%jmtd),&
     &                          integrand(atoms%jmtd)
      COMPLEX               ::  cexp(atoms%nat)
      COMPLEX               ::  exch(hybdat%nbands(nk),ncstd)
      COMPLEX               ::  cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat),carr(hybdat%nbands(nk))

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(//A)') '### core-core-core-valence exchange  ###'
        WRITE(6,'(/A)') '        k-point       band    exchange'
      END IF
            
      ! set up point
      icst = 0
      iatom= 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l = 0,hybdat%lmaxc(itype)
            DO M = -l,l
              DO n = 1,hybdat%nindxc(l,itype)
                icst = icst + 1
                point(n,M,l,iatom) = icst
              END DO
            END DO
          END DO
        END DO
      END DO
      
      ! lmstart = lm start index for each l-quantum number and atom type (for cmt-coefficients)
      DO itype = 1,atoms%ntype
        DO l = 0,atoms%lmax(itype)
          lmstart(l,itype) = sum( (/ (hybrid%nindx(ll,itype)*(2*ll+1),ll=0,l-1) /) )
        END DO
      END DO

      ! read in cmt coefficient at k-point nk
      irecl_cmt = dimension%neigd*hybrid%maxlmindx*atoms%nat*16
      OPEN(unit=777,file='cmt',form='unformatted',access='direct', recl=irecl_cmt)
      READ(777,rec=nk)  cmt(:,:,:)      
      CLOSE(777)
  
      iatom = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom       = iatom + 1
          cexp(iatom) = exp ( img*tpi_const* dot_product(bkpt(:),atoms%taual(:,iatom)) )
        END DO
      END DO   

      cmt = conjg(cmt)

      llmax = max(2*hybdat%lmaxcd,atoms%lmaxd)

      exch   = 0
      iatom0 = 0
      DO itype = 1,atoms%ntype
                
        DO l1 = 0,hybdat%lmaxc(itype)  ! left core state
          DO l2 = 0,atoms%lmax(itype)  ! right valence state
            DO l = 0,hybdat%lmaxc(itype)   ! occupied core state
            
              DO ll = abs(l1-l),l1+l
                IF( ll .lt. abs(l-l2) .or. ll .gt. l+l2 ) CYCLE
                IF( mod(l+l1+ll,2) .ne. 0 ) CYCLE
                IF( mod(l+l2+ll,2) .ne. 0 ) CYCLE
                
!                 WRITE(*,*) 'l1,l2,l,ll',l1,l2,l,ll
                rdum0 = fpi_const/(2*ll+1)
                
                DO m1 = -l1,l1
                  m2 = m1
                  IF( abs(m2) .GT. l2 ) CYCLE
                  lm2 = lmstart(l2,itype) +(m2+l2)*hybrid%nindx(l2,itype)
                  
                  DO M = -l,l
                    mm = M - m1
                    IF( abs(M-m1) .gt. ll ) CYCLE
                    
                    rdum1 = gaunt1(l,ll,l1,M,mm,m1,llmax)&
                         *gaunt1(l,ll,l2,M,mm,m1,llmax)*rdum0
                     
                    DO n = 1,hybdat%nindxc(l,itype)
                      DO n2 = 1,hybrid%nindx(l2,itype)
                        lmp2 = lm2 + n2
                        
                        rprod(:) = ( hybdat%core1(:,n,l,itype)*hybdat%bas1(:,n2,l2,itype)&
                              +hybdat%core2(:,n,l,itype)*hybdat%bas2(:,n2,l2,itype) ) / atoms%rmsh(:,itype)

                        CALL primitivef(primf1,rprod(:)*atoms%rmsh(:,itype)**(ll+1) ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd, itype,atoms%ntype)
                        CALL primitivef(primf2,rprod(:)/atoms%rmsh(:,itype)**ll ,atoms%rmsh,atoms%dx,atoms%jri,atoms%jmtd,-itype,atoms%ntype)  ! -itype is to enforce inward integration
                  

                       
                        
                        primf1 = primf1 / atoms%rmsh(:,itype)**ll
                        primf2 = primf2 * atoms%rmsh(:,itype)**(ll+1)
                
                        DO n1 = 1,hybdat%nindxc(l1,itype)
                          
                          rprod(:) = ( hybdat%core1(:,n,l,itype)*hybdat%core1(:,n1,l1,itype)&
                               +hybdat%core2(:,n,l,itype)*hybdat%core2(:,n1,l1,itype) ) / atoms%rmsh(:,itype)
      
                          integrand     = rprod * (primf1 + primf2)
                          
                          rdum2 = rdum1*intgrf(integrand,atoms%jri,atoms%jmtd, atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)
                          
                          iatom = iatom0
                          DO ieq = 1,atoms%neq(itype)
                            iatom = iatom + 1
                            icst1 = point(n1,m1,l1,iatom)
                            cdum  = rdum2*cexp(iatom)
                            DO iband = 1,hybdat%nbands(nk)

                              exch(iband,icst1) = exch(iband,icst1) + cdum*cmt(iband,lmp2,iatom)

                            END DO
                          END DO

                        END DO  !n1
                      
                      END DO  !n2
                    END DO  !n
                    
                  END DO  !M
                END DO  !m1
                
              END DO  !ll
              
            END DO  !l
          END DO  !l2
        END DO  !l1
        iatom0 = iatom0 + atoms%neq(itype)
      END DO  !itype      
      

#ifdef CPP_INVERSION
      !symmetrize core-wavefunctions such that phi(-r) = phi(r)*
      CALL symmetrize( exch,hybdat,ncstd,2,.false.,&
     &                 atoms,ntype,&
     &                 sym)
     
      IF( any( abs(aimag(exch)) .gt. 1E-6 ) ) STOP 'exchange_cccv: exch possesses significant imaginary part'
#endif


      DO icst = 1,ncstd
        DO iband = 1,hybdat%nbands(nk)
          exch_cv(iband,icst,nk) = exch_cv(iband,icst,nk) - exch(iband,icst)
        END DO
      END DO

      END SUBROUTINE exchange_cccv
      
      
      END MODULE m_exchange_core
