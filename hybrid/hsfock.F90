MODULE m_hsfock
  USE m_judft
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
!     This module is the driver routine for the calculation of the Hartree    c
!     Fock exchange term by using the mixed basis set.                        c
!                                                                             c
!     hsfock                                                                  c
!         |                                                                   c
!         |- symm.F:                                                          c
!         |  calculates the irreducible representation                        c
!         |                                                                   c
!         |- wavefproducts.F:                 s      s*                       c
!         |  computes the repsentation of phi    phi       in the mixed basis c
!         |                                  n,k    n',k+q                    c
!         |                                                                   c
!         |- exchange.F:                                                      c
!         |  calculates valence-valence part of the exchange matrix (mat_ex), c
!         |                                                                   c
!         |- exchange_core.F                                                  c
!         |  calculate valence-core contribution                              c
!                                                                             c
!     variables:                                                              c
!         nkptf   :=   number of kpoints                                      c
!         nkpti   :=   number of irreducible kpoints                          c
!         nbands  :=   number of bands for which the exchange matrix (mat_ex) c
!                      in the space of the wavefunctions is calculated        c
!         te_hfex :=   hf exchange contribution to the total energy           c
!         mnobd   :=   maximum number of occupied bands                       c
!         parent  :=   parent(ikpt) points to the symmetry equivalent point   c
!                      under the little group of kpoint nk                    c
!         symop   :=   symop(ikpt) points to the symmetry operation, which    c
!                      maps parent(ikpt) on ikpt                              c
!                                                                             c
!                                                                             c
!                                               M.Betzinger (09/07)           c
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
      CONTAINS

      SUBROUTINE hsfock(&
     &             nk,atoms,lcutm,obsolete,lapw,dimension,&
     &             kpts,nkpti,jsp,input,&
     &             hybrid,maxbasm,maxindxp,maxlcutm,&
     &             maxindxm,nindxm,&
     &             basm,bas1,bas2,bas1_MT,drbas1_MT,&
     &             ne_eig,eig_irr,n_size,sym,&
     &             cell,noco,oneD,&
     &             nbasp,nbasm,&
     &             results,&
     &             it,nbands,maxbands,nobd,mnobd,&
     &             xcpot,&
     &             core1,core2,nindxc,maxindxc,lmaxc,lmaxcd,&
     &             kveclo_eig,&
     &             maxfac,fac,sfac,gauntarr,&
     &             nindxp,prod,prodm,gwc,&
     &             mpi,irank2,isize2,comm,&
     &             a)

      USE m_symm_hf       ,ONLY: symm_hf
      USE m_util          ,ONLY: intgrf,intgrf_init
      USE m_exchange_valence_hf
      USE m_exchange_core
      USE m_symmetrizeh
      USE m_wrapper
      USE m_hsefunctional ,ONLY: exchange_vccvHSE,exchange_ccccHSE
      USE m_hybridmix
      USE m_icorrkeys
      USE m_types
      IMPLICIT NONE

      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_xcpot),INTENT(IN)   :: xcpot
      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_obsolete),INTENT(IN)   :: obsolete
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_noco),INTENT(IN)   :: noco
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)   :: lapw

!     - scalars -
      INTEGER,INTENT(IN)      :: jsp 
      INTEGER,INTENT(IN)      :: nbands
      INTEGER,INTENT(IN)      :: it
      INTEGER,INTENT(IN)      ::  n_size 
      INTEGER,INTENT(IN)      :: irank2 ,isize2,comm
      INTEGER,INTENT(IN)      ::  nkpti  ,nk
      INTEGER,INTENT(IN)      ::  maxindxp,maxindxm,maxlcutm,maxbasm
      INTEGER,INTENT(IN)      :: maxbands 
      INTEGER,INTENT(IN)      ::  nbasp
      INTEGER,INTENT(IN)      ::  mnobd
      INTEGER,INTENT(IN)      ::  lmaxcd,maxindxc
      INTEGER,INTENT(IN)      ::  maxfac
      INTEGER,INTENT(IN)      :: gwc


      !     -  arrays -
      INTEGER,INTENT(IN)      ::  nindxm(0:maxlcutm,atoms%ntype)
      INTEGER,INTENT(IN)      ::  lcutm(atoms%ntype) 
      INTEGER,INTENT(IN)      ::  ne_eig(nkpti) 

      INTEGER,INTENT(IN)      ::  kveclo_eig(atoms%nlotot,nkpti)
      INTEGER,INTENT(IN)      ::  nbasm(kpts%nkptf)
      INTEGER,INTENT(IN)      ::  nobd(kpts%nkptf)
      INTEGER,INTENT(IN)      ::  nindxc(0:lmaxcd,atoms%ntype)
      INTEGER,INTENT(IN)      ::  lmaxc(atoms%ntype)
  
      INTEGER,INTENT(IN)      ::  nindxp(0:maxlcutm,atoms%ntype)

      REAL,INTENT(IN)         ::  bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
     &                            bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
      REAL,INTENT(IN)         ::    bas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
     &                            drbas1_MT(hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
      REAL                    ::  basm(atoms%jmtd,maxindxm,0:maxlcutm,atoms%ntype)
      REAL,INTENT(IN)         ::  eig_irr(dimension%neigd,nkpti)
      REAL,INTENT(IN)         ::  fac(0:maxfac),sfac(0:maxfac)
      REAL,INTENT(IN)         ::  gauntarr(2,0:atoms%lmaxd,0:atoms%lmaxd,0:maxlcutm,&
     &                                  -atoms%lmaxd:atoms%lmaxd,-maxlcutm:maxlcutm)
      REAL,INTENT(IN)         ::  core1(atoms%jmtd,maxindxc,0:lmaxcd,atoms%ntype),&
     &                            core2(atoms%jmtd,maxindxc,0:lmaxcd,atoms%ntype)

      REAL,INTENT(IN)      ::  prodm(maxindxm,maxindxp,0:maxlcutm,atoms%ntype)

#ifdef CPP_INVERSION
      REAL,INTENT(INOUT)      ::  a(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#else
      COMPLEX,INTENT(INOUT)   ::  a(dimension%nbasfcn*(dimension%nbasfcn+1)/2)
#endif
#if ( !defined(CPP_INVERSION) )
      COMPLEX,ALLOCATABLE     ::  z(:,:)
#else
      REAL,ALLOCATABLE        ::  z(:,:)
#endif
  
      TYPE(PRODTYPE),INTENT(IN)::  prod(maxindxp,0:maxlcutm,atoms%ntype)

!     - local scalars -
      INTEGER                 ::  i,j,ic,ic1,l,itype
      INTEGER                 ::  iband,iband1,iband2
      INTEGER                 ::  ikpt,ikpt0
      INTEGER                 ::  irec
      INTEGER                 ::  irecl_olap,irecl_z,irecl_vx
      INTEGER                 ::  irecl_gw
      INTEGER                 ::  maxndb
      INTEGER                 ::  nddb,ndb1,ndb2
      INTEGER                 ::  nsymop
      INTEGER                 ::  nkpt_EIBZ
      INTEGER                 ::  ncstd
      INTEGER                 ::  ok

      REAL                    ::  a_ex
!     - local arrays -
      INTEGER                 ::  gpt(3,lapw%nv(jsp))
      INTEGER                 ::  degenerat(ne_eig(nk))
      INTEGER                 ::  nsest(nbands),indx_sest(nbands,nbands)

      INTEGER,ALLOCATABLE     ::  parent(:),symop(:)
      INTEGER,ALLOCATABLE     ::  psym(:)
      INTEGER,ALLOCATABLE     ::  pointer_EIBZ(:)
      INTEGER,ALLOCATABLE     ::  n_q(:)

      REAL                    ::  wl_iks(dimension%neigd,kpts%nkptf)
      REAL                    ::  div_vv(nbands)
      REAL,ALLOCATABLE        ::  gridf(:,:)


#ifdef CPP_INVERSION
      REAL,ALLOCATABLE        ::  olap(:,:),olap_p(:)
#else
      COMPLEX,ALLOCATABLE     ::  olap(:,:),olap_p(:)
#endif

#ifdef CPP_INVERSION
      REAL   ,ALLOCATABLE     ::  trafo(:,:),invtrafo(:,:)
      REAL   ,ALLOCATABLE     ::  ex(:,:),v(:,:),v_x(:),v_xp(:)
      REAL   ,ALLOCATABLE     ::  mat_ex(:)
#else
      COMPLEX,ALLOCATABLE     ::  trafo(:,:),invtrafo(:,:)
      COMPLEX,ALLOCATABLE     ::  ex(:,:),v(:,:),v_x(:),v_xp(:)
      COMPLEX,ALLOCATABLE     ::  mat_ex(:)
#endif
      COMPLEX                 ::  exch(dimension%neigd,dimension%neigd)
      COMPLEX,ALLOCATABLE     ::  carr(:)
      COMPLEX,ALLOCATABLE     ::  rep_c(:,:,:,:,:)

      real :: rdum

      
      CALL timestart("total time hsfock")
      CALL timestart("symm_hf")
    

      !
      ! preparations
      !
      
      ! initialize gridf for radial integration
      CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,gridf)

      !
      ! initialize weighting factor for HF exchange part
      !
      IF     ( xcpot%icorr .eq. icorr_pbe0 ) THEN
        a_ex = amix_pbe0
      ELSE IF( xcpot%icorr .eq. icorr_hf   ) THEN
        a_ex = amix_hf
      ELSE IF (xcpot%icorr .eq. icorr_hse) THEN
        a_ex = aMix_HSE
      ELSE IF (xcpot%icorr .eq. icorr_vhse ) THEN
        a_ex = aMix_VHSE()
      ELSE
        STOP 'mhsfock: xc functional can not be identified'
      END IF


      ! write k1,k2,k3 in gpt
      DO i=1,lapw%nv(jsp)
        gpt(1,i) = lapw%k1(i,jsp)
        gpt(2,i) = lapw%k2(i,jsp)
        gpt(3,i) = lapw%k3(i,jsp)
      END DO

      ! read in lower triangle part of overlap matrix from direct acces file olap
      ALLOCATE( olap_p(dimension%nbasfcn*(dimension%nbasfcn+1)/2), stat=ok)
      IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation olap_p'
      olap_p = 0

#ifdef CPP_INVERSION
      irecl_olap = dimension%nbasfcn*(dimension%nbasfcn+1)*4
#else
      irecl_olap = dimension%nbasfcn*(dimension%nbasfcn+1)*8
#endif
      irec = nkpti*(jsp-1) + nk
      print *, "Olap read:",irec
      OPEN(88,file='olap',form='unformatted',access='direct',&
     &     recl=irecl_olap)
      READ(88,rec=irec) olap_p
      CLOSE(88)

      !unpack olap_p into olap
      ALLOCATE( olap(dimension%nbasfcn,dimension%nbasfcn), stat=ok)
      IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation olap'
      olap = 0
      ic = 0
      DO i=1,dimension%nbasfcn
        DO j=1,i
          ic = ic + 1
          olap(i,j) = olap_p(ic)
#ifdef CPP_INVERSION
          olap(j,i) = olap(i,j)
#else
          olap(j,i) = conjg(olap_p(ic))
#endif
        END DO
      END DO
      DEALLOCATE(olap_p)

      ALLOCATE( mat_ex(dimension%nbasfcn*(dimension%nbasfcn+1)),stat=ok )
      IF ( ok .ne. 0) STOP 'mhsfock: failure allocation mat_ex'
      mat_ex = 0

      IF( hybrid%l_calhf ) THEN
        ncstd = sum( (/ ( (nindxc(l,itype)*(2*l+1)*atoms%neq(itype),&
     &             l=0,lmaxc(itype)), itype = 1,atoms%ntype) /) )
        IF( nk .eq. 1 .and. mpi%irank == 0 )&
     &      WRITE(*,*) 'calculate new HF matrix'
        IF( nk .eq. 1 .and. jsp .eq. 1 .and. input%imix .gt. 10)&
     &      CALL system('rm -f broyd*')

        ! calculate all symmetrie operations, which yield k invariant

        ALLOCATE( parent(kpts%nkptf),symop(kpts%nkptf) ,stat=ok)
        IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation parent/symop'
        parent = 0 ; symop = 0

        CALL symm_hf( kpts,nkpti,nk,sym,&
     &             dimension,ne_eig(nk),eig_irr,nbands,&
     &             atoms,hybrid,bas1,bas2,cell,&
     &             lapw,jsp,&
     &             gpt,&
     &             lmaxcd,&
     &             mpi,irank2,&
     &             nsymop,psym,nkpt_EIBZ,n_q,parent,&
     &             symop,degenerat,pointer_EIBZ,maxndb,nddb,&
     &             nsest,indx_sest,rep_c)

        CALL timestop("symm_hf")
        ! remove weights(wtkpt) in w_iks
        DO ikpt=1,kpts%nkptf
          DO iband=1,dimension%neigd
            ikpt0 = kpts%bkp(ikpt)
            wl_iks(iband,ikpt) = results%w_iks(iband,ikpt0,jsp) /&
     &                           ( kpts%wtkpt(ikpt0) * kpts%nkptf )
          END DO
        END DO

        !
        ! calculate contribution from valence electrons to the
        ! HF exchange
        CALL timestart("valence exchange calculation")

        CALL exchange_valence_hf(&
     &            nk,kpts,nkpti,nkpt_EIBZ,&
     &            sym,atoms,hybrid,&
     &            cell,&
     &            dimension,input,jsp,&
     &            basm,bas1,bas2,bas1_MT,&
     &            drbas1_MT,maxlcutm,lcutm,nindxm,maxindxm,nbasp,&
     &            nbasm,maxbasm,maxindxp,nindxp,&
     &            prod,prodm,mnobd,&
     &            nobd,nbands,ne_eig(nk),lapw,&
     &            eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,&
     &            kveclo_eig,gauntarr,it,xcpot,&
     &            noco,nsest,indx_sest,&
     &            mpi,irank2,isize2,comm,&
     &            div_vv,mat_ex)

        DEALLOCATE ( rep_c )
        CALL timestop("valence exchange calculation")

        CALL timestart("core exchange calculation")
        ! do the rest of the calculation only on master
        IF ( irank2 /= 0 ) RETURN

  
        ! calculate contribution from the core states to the HF exchange
        IF ( xcpot%icorr.eq.icorr_hse .OR. xcpot%icorr.eq.icorr_vhse ) THEN
#ifdef CPP_NEVER           
          CALL exchange_vccvHSE(&
     &                 nk,kpts,nkpti,atoms,&
     &                 hybrid,lmaxc,&
     &                 nindxc,maxindxc,core1,core2,lcutm,&
     &                 bas1,bas2,dimension,jsp,&
     &                 maxfac,fac,sfac,lapw,nbands,&
     &                 gridf,nsymop,nsest,indx_sest,mpi,&
     &                 a_ex,nobd,results,&
     &                 mat_ex%core )
          CALL exchange_ccccHSE(&
     &                 nk,nkpti,obsolete,atoms,lmaxcd,lmaxc,&
     &                 nindxc,maxindxc,ncstd,&
     &                 core1,core2,kpts(:,nk),gridf,&
     &                 sym,maxfac,fac,a_ex,mpi,&
     &                 results%core )
#endif
          STOP "HSE not implemented in hsfock"
        ELSE
          CALL exchange_vccv1(&
     &                 nk,kpts,nkpti,atoms,&
     &                 hybrid,lmaxc,&
     &                 nindxc,maxindxc,core1,core2,lcutm,&
     &                 bas1,bas2,dimension,jsp,&
     &                 maxfac,fac,sfac,lapw,nbands,&
     &                 gridf,nsymop,nsest,indx_sest,mpi,&
     &                 a_ex,nobd,results,&
     &                 mat_ex)
          CALL exchange_cccc(&
     &                 nk,nkpti,atoms,lmaxcd,lmaxc,&
     &                 nindxc,maxindxc,ncstd,&
     &                 core1,core2,gridf,&
     &                 sym,kpts,a_ex,mpi,&
     &                 results )
        END IF

        DEALLOCATE( n_q )
        CALL timestop("core exchange calculation")

        CALL timestart("time for performing T^-1*mat_ex*T^-1*")
        !calculate trafo from wavefunctions to APW basis
        IF( dimension%neigd .lt. nbands ) STOP 'mhsfock:   < nbands; &&
     &trafo from wavefunctions to APW requires at least nbands'

        ALLOCATE(z(dimension%nbasfcn,dimension%neigd),stat=ok)
        IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation z'
        z = 0
#ifdef CPP_INVERSION
        irecl_z   =  dimension%nbasfcn*dimension%neigd*8
#else
        irecl_z   =  dimension%nbasfcn*dimension%neigd*16
#endif

        OPEN(unit=778,file='z',form='unformatted',access='direct',&
     &       recl=irecl_z)
        READ(778,rec=nk) z
        CLOSE(778)

        !unpack mat_ex
        ALLOCATE( ex(nbands,nbands), stat = ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation ex'
        ex = 0

        ic = 0
        DO i = 1,nbands
          DO j = 1,i
            ic = ic + 1
            ex(j,i) = mat_ex(ic)
#ifdef CPP_INVERSION
            ex(i,j) = ex(j,i)
#else
            ex(i,j) = conjg(mat_ex(ic))
#endif
          END DO
        END DO


        ! calculate trafo
        ic = lapw%nv(jsp) + atoms%nlotot

        ALLOCATE( trafo(ic,nbands), stat= ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation trafo'
        trafo  = matmul(olap(:ic,:ic),z(:ic,:nbands))

        ALLOCATE( invtrafo(nbands,ic), stat= ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation invtrafo'

#ifdef CPP_INVERSION
        invtrafo = transpose(trafo)
#else
        invtrafo = conjg(transpose(trafo))
#endif

        ALLOCATE( v(ic,ic), stat = ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation v'

        v = matmul(trafo,matmul(ex,invtrafo))

        DEALLOCATE( mat_ex,ex,trafo,invtrafo,olap )
        CALL timestop("time for performing T^-1*mat_ex*T^-1*")

        !store only lower triangle of v
        ALLOCATE( v_x(dimension%nbasfcn*(dimension%nbasfcn+1)/2), stat = ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation v_x'
        v_x = 0
        ic  = 0
        DO i=1,lapw%nv(jsp)+atoms%nlotot
          DO j = 1,i
            ic      = ic + 1
            v_x(ic) = v(i,j)
          END DO
        END DO

        CALL symmetrizeh(atoms,&
     &                   kpts%bk(:,nk),dimension,jsp,lapw,gpt,&
     &                   sym,kveclo_eig,&
     &                   cell,nsymop,psym,&
     &                   v_x )



        IF( input%imix .ge. 10 ) THEN
#ifdef CPP_INVERSION
          irecl_vx = dimension%nbasfcn*(dimension%nbasfcn+1)*4
#else
          irecl_vx = dimension%nbasfcn*(dimension%nbasfcn+1)*8
#endif

          irec = nkpti*(jsp-1) + nk
          OPEN(778,file='vex',form='unformatted',access='direct',&
     &         recl=irecl_vx)
#ifdef CPP_INVERSION
          WRITE(778,rec=irec) v_x
#else
          WRITE(778,rec=irec) v_x
#endif
          CLOSE(778)
        END IF

      ELSE ! not hybrid%l_calhf

        ALLOCATE(z(dimension%nbasfcn,dimension%neigd),stat=ok)
        IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation z'
        z = 0
#ifdef CPP_INVERSION
        irecl_z   =  dimension%nbasfcn*dimension%neigd*8
#else
        irecl_z   =  dimension%nbasfcn*dimension%neigd*16
#endif

        OPEN(unit=778,file='z',form='unformatted',access='direct',&
     &       recl=irecl_z)
        READ(778,rec=nk) z
        CLOSE(778)

        ALLOCATE( v_x(dimension%nbasfcn*(dimension%nbasfcn+1)/2), stat = ok )
        IF( ok .ne. 0 ) STOP 'mhsfock: error allocation v_x'
        v_x = 0
#ifdef CPP_INVERSION
        irecl_vx = dimension%nbasfcn*(dimension%nbasfcn+1)*4
#else
        irecl_vx = dimension%nbasfcn*(dimension%nbasfcn+1)*8
#endif
        irec = nkpti*(jsp-1) + nk
        OPEN(778,file='vex',form='unformatted',access='direct',&
     &       recl=irecl_vx)
        READ(778,rec=irec) v_x
        CLOSE(778)

      END IF ! hybrid%l_calhf

      ! add non-local x-potential to the hamiltonian a
      a = a - a_ex*v_x


      ! calculate HF energy
      IF( hybrid%l_calhf ) THEN
        WRITE(6,'(A)') new_line('n')//new_line('n')//' ###     '// '        diagonal HF exchange elements (eV)              ###'
        
        WRITE(6,'(A)') new_line('n') // '         k-point      '// 'band       tail      pole     input%total(valence+core)'
        
      END IF
      ic  = lapw%nv(jsp) + atoms%nlotot
      ic1 = ic*(ic+1)/2

#if( !defined CPP_INVERSION )
      v_x(:ic1) = conjg(v_x(:ic1))
#endif

      ! calculate exchange contribution of current k-point nk to total energy (te_hfex)
      ! in the case of a spin-unpolarized calculation the factor 2 is added in eigen_hf.F

      IF( input%gw .eq.2 .and. gwc .eq. 1 ) THEN
        exch = 0
        ALLOCATE( carr(ic),stat=ok )
        IF( ok .ne. 0 ) STOP 'hsfock: error allocation carr'
        DO iband1 = 1,nbands
          carr = matvec(v_x(:ic1),z(:ic,iband1))
          DO iband2 = 1,iband1
            exch(iband2,iband1) = dotprod(z(:ic,iband2),carr(:ic))
            exch(iband1,iband2) = conjg(exch(iband2,iband1))
          END DO
        END DO

        DO iband = 1,nbands
          IF( iband .le. nobd(nk) ) THEN
            results%te_hfex%valence = results%te_hfex%valence -&
     &                        a_ex*results%w_iks(iband,nk,jsp)*exch(iband,iband)
          END IF
          IF(hybrid%l_calhf) THEN
            WRITE(6, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F10.5)')&
     &  kpts%bk(:,nk),iband, (REAL(exch(iband,iband))-div_vv(iband))*(-27.211608),&
     &  div_vv(iband)*(-27.211608),REAL(exch(iband,iband))*(-27.211608)
          END IF
        END DO

        ! write exch(:,:) to file gw_vxnl
        irec     = (jsp-1)*nkpti + nk
        irecl_gw = dimension%neigd*dimension%neigd*16
        OPEN(unit=778,file='gw_vxnl',access='direct',recl=irecl_gw)
        WRITE(778,rec=irec) -exch
        CLOSE(778)
      ELSE
        exch = 0
        DO iband = 1,nbands
          exch(iband,iband) = dotprod(z(:ic,iband),matvec(v_x(:ic1),z(:ic,iband)))
          IF( iband .le. nobd(nk) ) THEN
            results%te_hfex%valence = results%te_hfex%valence -&
     &                        a_ex*results%w_iks(iband,nk,jsp)*exch(iband,iband)
          END IF
          IF(hybrid%l_calhf) THEN
            WRITE(6, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F10.5)')&
     &  kpts%bk(:,nk),iband, (REAL(exch(iband,iband))-div_vv(iband))*(-27.211608),&
     &  div_vv(iband)*(-27.211608),REAL(exch(iband,iband))*(-27.211608)
          END IF
        END DO
      END IF

      DEALLOCATE( z,v_x )

      CALL timestop("total time hsfock")

      END SUBROUTINE hsfock


      END MODULE m_hsfock
