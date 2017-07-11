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
!         kpts%nkptf   :=   number of kpoints                                      c
!         kpts%nkpt   :=   number of irreducible kpoints                          c
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
     &             nk,atoms,hybrid,lapw,dimension,&
     &             kpts,jsp,input,&
     &             hybdat,&
     &             eig_irr,sym,&
     &             cell,noco,&
     &             results,&
     &             it,mnobd,&
     &             xcpot,&
     &             mpi,irank2,isize2,comm,&
     &             hamovlp)

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
      USE m_io_hybrid
      IMPLICIT NONE
      TYPE(t_hybdat),INTENT(IN)       :: hybdat
      TYPE(t_results),INTENT(INOUT)   :: results
      TYPE(t_xcpot),INTENT(IN)        :: xcpot
      TYPE(t_mpi),INTENT(IN)          :: mpi
      TYPE(t_dimension),INTENT(IN)    :: dimension
      TYPE(t_hybrid),INTENT(INOUT)    :: hybrid
      TYPE(t_input),INTENT(IN)        :: input
      TYPE(t_noco),INTENT(IN)         :: noco
      TYPE(t_sym),INTENT(IN)          :: sym
      TYPE(t_cell),INTENT(IN)         :: cell
      TYPE(t_kpts),INTENT(IN)         :: kpts
      TYPE(t_atoms),INTENT(IN)        :: atoms
      TYPE(t_lapw),INTENT(IN)         :: lapw

!     - scalars -
      INTEGER,INTENT(IN)      :: jsp 
      INTEGER,INTENT(IN)      :: it
      INTEGER,INTENT(IN)      :: irank2 ,isize2,comm
      INTEGER,INTENT(IN)      :: nk
      INTEGER,INTENT(IN)      :: mnobd
    

      !     -  arrays -
      REAL,INTENT(IN)         ::  eig_irr(dimension%neigd,kpts%nkpt)
      
      TYPE(t_hamovlp),INTENT(INOUT)::hamovlp

  

!     - local scalars -
      INTEGER                 ::  i,j,ic,ic1,l,itype,n,nn
      INTEGER                 ::  iband,iband1,iband2
      INTEGER                 ::  ikpt,ikpt0
      INTEGER                 ::  irec
      INTEGER                 ::  irecl_olap,irecl_z,irecl_vx
      INTEGER                 ::  maxndb
      INTEGER                 ::  nddb
      INTEGER                 ::  nsymop
      INTEGER                 ::  nkpt_EIBZ
      INTEGER                 ::  ncstd
      INTEGER                 ::  ok

      REAL                    ::  a_ex
!     - local arrays -
      INTEGER                 ::  gpt(3,lapw%nv(jsp))
      INTEGER                 ::  degenerat(hybdat%ne_eig(nk))
      INTEGER                 ::  nsest(hybdat%nbands(nk)),indx_sest(hybdat%nbands(nk),hybdat%nbands(nk))

      INTEGER,ALLOCATABLE     ::  parent(:),symop(:)
      INTEGER,ALLOCATABLE     ::  psym(:)
      INTEGER,ALLOCATABLE     ::  pointer_EIBZ(:)
      INTEGER,ALLOCATABLE     ::  n_q(:)

      REAL                    ::  wl_iks(dimension%neigd,kpts%nkptf)
      REAL                    ::  div_vv(hybdat%nbands(nk))

      TYPE(t_mat)             :: olap,trafo,invtrafo,ex,tmp,v_x,z
      COMPLEX                 ::  exch(dimension%neigd,dimension%neigd)
      COMPLEX,ALLOCATABLE     ::  carr(:)
      COMPLEX,ALLOCATABLE     ::  rep_c(:,:,:,:,:)

      
      CALL timestart("total time hsfock")
      CALL timestart("symm_hf")
    

      !
      ! preparations
      !
      
      ! initialize gridf for radial integration
      !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

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
        gpt(:,i) = (/lapw%k1(i,jsp),lapw%k2(i,jsp),lapw%k3(i,jsp)/)
      END DO

      
      ! read in lower triangle part of overlap matrix from direct acces file olap
      call olap%alloc(hamovlp%l_real,dimension%nbasfcn)
      call read_olap(olap, kpts%nkpt*(jsp-1) + nk)
      if (.not.olap%l_real) olap%data_c=conjg(olap%data_c)

      IF( hybrid%l_calhf ) THEN
        ncstd = sum( (/ ( (hybdat%nindxc(l,itype)*(2*l+1)*atoms%neq(itype),&
     &             l=0,hybdat%lmaxc(itype)), itype = 1,atoms%ntype) /) )
        IF( nk .eq. 1 .and. mpi%irank == 0 )&
     &      WRITE(*,*) 'calculate new HF matrix'
        IF( nk .eq. 1 .and. jsp .eq. 1 .and. input%imix .gt. 10)&
     &      CALL system('rm -f broyd*')
        ! calculate all symmetrie operations, which yield k invariant

        ALLOCATE( parent(kpts%nkptf),symop(kpts%nkptf) ,stat=ok)
        IF( ok .ne. 0 ) STOP 'mhsfock: failure allocation parent/symop'
        parent = 0 ; symop = 0

        CALL symm_hf( kpts,nk,sym,&
     &             dimension,hybdat,eig_irr,&
     &             atoms,hybrid,cell,&
     &             lapw,jsp,&
     &             gpt,&
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
     &            nk,kpts,nkpt_EIBZ, sym,atoms,hybrid,&
     &            cell, dimension,input,jsp, hybdat, mnobd, lapw,&
     &            eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,&
     &            it,xcpot,&
     &            noco,nsest,indx_sest,&
     &            mpi,irank2,isize2,comm,&
     &            div_vv,ex)

        DEALLOCATE ( rep_c )
        CALL timestop("valence exchange calculation")

        CALL timestart("core exchange calculation")
        ! do the rest of the calculation only on master
        IF ( irank2 /= 0 ) RETURN

  
        ! calculate contribution from the core states to the HF exchange
        IF ( xcpot%icorr.eq.icorr_hse .OR. xcpot%icorr.eq.icorr_vhse ) THEN
#ifdef CPP_NEVER           
          CALL exchange_vccvHSE(&
     &                 nk,atoms,&
     &                 hybrid,hybdat,&
     &                 dimension,jsp,&
     &                 lapw,&
     &                 nsymop,nsest,indx_sest,mpi,&
     &                 a_ex,results,&
     &                 mat_ex%core )
          CALL exchange_ccccHSE(&
     &                 nk,obsolete,atoms,hybdat,&
     &                 ncstd,&
     &                 kpts(:,nk),&
     &                 sym,a_ex,mpi,&
     &                 results%core )
#endif
          STOP "HSE not implemented in hsfock"
        ELSE
          CALL exchange_vccv1(&
     &                 nk,atoms,&
     &                 hybrid,hybdat,&
     &                 dimension,jsp,&
     &                 lapw,&
     &                 nsymop,nsest,indx_sest,mpi,&
     &                 a_ex,results,&
     &                 ex)
          CALL exchange_cccc(&
     &                 nk,atoms,hybdat,&
     &                 ncstd,&
     &                 sym,kpts,a_ex,mpi,&
     &                 results )
        END IF

        DEALLOCATE( n_q )
        CALL timestop("core exchange calculation")

        CALL timestart("time for performing T^-1*mat_ex*T^-1*")
        !calculate trafo from wavefunctions to APW basis
        IF( dimension%neigd .lt. hybdat%nbands(nk) ) STOP 'mhsfock: neigd  < nbands(nk) ; &&
     &trafo from wavefunctions to APW requires at least nbands(nk) '

        call z%alloc(olap%l_real,dimension%nbasfcn,dimension%neigd)

        call read_z(z,nk) !what about spin?
      
        ! calculate trafo
        ic = lapw%nv(jsp) + atoms%nlotot
        z%matsize1=ic
        z%matsize2=hybdat%nbands(nk)
        olap%matsize1=ic
        olap%matsize2=ic
        
        call olap%multiply(z,trafo)
        
        call invtrafo%alloc(olap%l_real,hybdat%nbands(nk),ic)
        CALL trafo%transpose(invtrafo)
        
        CALL ex%multiply(invtrafo,tmp)
        CALL trafo%multiply(tmp,v_x)
        
        CALL timestop("time for performing T^-1*mat_ex*T^-1*")

      
        CALL symmetrizeh(atoms,&
     &                   kpts%bkf(:,nk),dimension,jsp,lapw,gpt,&
     &                   sym,hybdat%kveclo_eig,&
     &                   cell,nsymop,psym,&
     &                   v_x )



        IF( input%imix .ge. 10 ) THEN
          CALL write_v_x(v_x,kpts%nkpt*(jsp-1) + nk)
        END IF

      ELSE ! not hybrid%l_calhf
         CALL read_v_x(v_x,kpts%nkpt*(jsp-1) + nk)
     END IF ! hybrid%l_calhf
      ! add non-local x-potential to the hamiltonian a (in packed storage)
     ic=1
     DO n=1,v_x%matsize1
        DO nn=1,n           
           IF (hamovlp%l_real) THEN
              hamovlp%a_r(ic) = hamovlp%a_r(ic) - a_ex*v_x%data_r(n,nn)
           ELSE
              hamovlp%a_c(ic) = hamovlp%a_c(ic) - a_ex*v_x%data_c(n,nn)
           ENDIF
        ENDDO
     END DO
      ! calculate HF energy
      IF( hybrid%l_calhf ) THEN
         WRITE(6,'(A)') new_line('n')//new_line('n')//' ###     '// '        diagonal HF exchange elements (eV)              ###'
        
         WRITE(6,'(A)') new_line('n') // '         k-point      '// 'band       tail      pole     input%total(valence+core)'
        
      END IF
 
     
      ! calculate exchange contribution of current k-point nk to total energy (te_hfex)
      ! in the case of a spin-unpolarized calculation the factor 2 is added in eigen_hf.F
      if (.not.v_x%l_real) v_x%data_c=conjg(v_x%data_c) 
      exch = 0
      call v_x%multiply(z,tmp)
      DO iband = 1,hybdat%nbands(nk)
         if (z%l_real) THEN
            exch(iband,iband) = dot_product(z%data_r(:,iband),tmp%data_r(:,iband))
         else
            exch(iband,iband) = dot_product(z%data_r(:,iband),tmp%data_r(:,iband))
         endif
         IF( iband .le. hybdat%nobd(nk) ) THEN
            results%te_hfex%valence = results%te_hfex%valence -a_ex*results%w_iks(iband,nk,jsp)*exch(iband,iband)
         END IF
         IF(hybrid%l_calhf) THEN
            WRITE(6, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F10.5)')&
                 &  kpts%bkf(:,nk),iband, (REAL(exch(iband,iband))-div_vv(iband))*(-27.211608),&
                 &  div_vv(iband)*(-27.211608),REAL(exch(iband,iband))*(-27.211608)
         END IF
      END DO
     

      CALL timestop("total time hsfock")

      END SUBROUTINE hsfock


      END MODULE m_hsfock
