!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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

SUBROUTINE hsfock(nk,atoms,hybrid,lapw,dimension,kpts,jsp,input,hybdat,eig_irr,sym,cell,noco,&
                  results,it,mnobd,xcpot,mpi,irank2,isize2,comm)

   USE m_types
   USE m_symm_hf
   USE m_util          ,ONLY: intgrf,intgrf_init
   USE m_exchange_valence_hf
   USE m_exchange_core
   USE m_symmetrizeh
   USE m_wrapper
   USE m_hsefunctional ,ONLY: exchange_vccvHSE,exchange_ccccHSE
   USE m_io_hybrid

   IMPLICIT NONE

   TYPE(t_hybdat),        INTENT(IN)    :: hybdat
   TYPE(t_results),       INTENT(INOUT) :: results
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

   ! scalars
   INTEGER,               INTENT(IN)    :: jsp 
   INTEGER,               INTENT(IN)    :: it
   INTEGER,               INTENT(IN)    :: irank2 ,isize2,comm
   INTEGER,               INTENT(IN)    :: nk
   INTEGER,               INTENT(IN)    :: mnobd

   ! arrays
   REAL,                  INTENT(IN)    :: eig_irr(dimension%neigd,kpts%nkpt)

   ! local scalars
   INTEGER                 ::  i,j,ic,ic1,l,itype,n,nn
   INTEGER                 ::  iband,iband1,iband2
   INTEGER                 ::  ikpt,ikpt0
   INTEGER                 ::  irec
   INTEGER                 ::  irecl_olap,irecl_z,irecl_vx
   INTEGER                 ::  maxndb, nbasfcn
   INTEGER                 ::  nddb
   INTEGER                 ::  nsymop
   INTEGER                 ::  nkpt_EIBZ
   INTEGER                 ::  ncstd
   INTEGER                 ::  ok
   REAL                    ::  a_ex

   ! local arrays
   INTEGER                 ::  nsest(hybrid%nbands(nk)),indx_sest(hybrid%nbands(nk),hybrid%nbands(nk))
   INTEGER                 ::  rrot(3,3,sym%nsym)
   INTEGER                 ::  psym(sym%nsym) ! Note: psym is only filled up to index nsymop

   INTEGER,ALLOCATABLE     ::  parent(:)
   INTEGER,ALLOCATABLE     ::  pointer_EIBZ(:)
   INTEGER,ALLOCATABLE     ::  n_q(:)

   REAL                    ::  wl_iks(dimension%neigd,kpts%nkptf)
    
   TYPE(t_mat)             :: olap,trafo,invtrafo,ex,tmp,v_x,z
   COMPLEX                 ::  exch(dimension%neigd,dimension%neigd)
   COMPLEX,ALLOCATABLE     ::  carr(:)
      
   CALL timestart("total time hsfock")
    
   ! preparations

   ! initialize gridf for radial integration
   !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

   ! initialize weighting factor for HF exchange part
   a_ex=xcpot%get_exchange_weight()

   ! read in lower triangle part of overlap matrix from direct acces file olap
   nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
   call olap%alloc(sym%invs,nbasfcn)
   call read_olap(olap, kpts%nkpt*(jsp-1)+nk)
   IF (olap%l_real) THEN
      DO i=1,nbasfcn
        DO j=1,i
          olap%data_r(i,j) = olap%data_r(j,i)
        END DO
      END DO
   ELSE
      DO i=1,nbasfcn
        DO j=1,i
          olap%data_c(i,j) = CONJG(olap%data_c(j,i))
        END DO
      END DO
      olap%data_c=conjg(olap%data_c)
   END IF

   IF(hybrid%l_calhf) THEN
      ncstd = sum( (/ ( (hybdat%nindxc(l,itype)*(2*l+1)*atoms%neq(itype),l=0,hybdat%lmaxc(itype)), itype = 1,atoms%ntype) /) )
      IF( nk .eq. 1 .and. mpi%irank == 0 ) WRITE(*,*) 'calculate new HF matrix'
      IF( nk .eq. 1 .and. jsp .eq. 1 .and. input%imix .gt. 10) CALL system('rm -f broyd*')
      ! calculate all symmetrie operations, which yield k invariant

      ALLOCATE(parent(kpts%nkptf), stat=ok)
      IF(ok.NE.0) STOP 'mhsfock: failure allocation parent'
      parent = 0

      CALL timestart("symm_hf")
      CALL symm_hf_init(sym,kpts,nk,irank2,nsymop,rrot,psym)

      CALL symm_hf(kpts,nk,sym,dimension,hybdat,eig_irr,atoms,hybrid,cell,lapw,jsp,mpi,irank2,&
                   rrot,nsymop,psym,nkpt_EIBZ,n_q,parent,pointer_EIBZ,maxndb,nddb,nsest,indx_sest)
      CALL timestop("symm_hf")

      ! remove weights(wtkpt) in w_iks
      DO ikpt=1,kpts%nkptf
         DO iband=1,dimension%neigd
            ikpt0 = kpts%bkp(ikpt)
            wl_iks(iband,ikpt) = results%w_iks(iband,ikpt0,jsp) / (kpts%wtkpt(ikpt0)*kpts%nkptf)
         END DO
      END DO

      ! calculate contribution from valence electrons to the
      ! HF exchange
      CALL timestart("valence exchange calculation")
      ex%l_real=sym%invs
      CALL exchange_valence_hf(nk,kpts,nkpt_EIBZ, sym,atoms,hybrid,cell,dimension,input,jsp,hybdat,mnobd,lapw,&
                               eig_irr,results,parent,pointer_EIBZ,n_q,wl_iks,it,xcpot,noco,nsest,indx_sest,&
                               mpi,irank2,isize2,comm,ex)
      CALL timestop("valence exchange calculation")

      WRITE(1224,'(a,i7)') 'kpoint: ', nk
      DO i = 1, ex%matsize1
         DO j = 1, i
            IF (ex%l_real) THEN
               WRITE(1224,'(2i7,2f15.8)') i, j, ex%data_r(j,i) !ex%data_r(i,j), ex%data_r(j,i)
            ELSE
               WRITE(1224,'(2i7,4f15.8)') i, j, ex%data_c(j,i) !ex%data_c(i,j), ex%data_c(j,i)
            ENDIF
         END DO
      END DO

      CALL timestart("core exchange calculation")
      ! do the rest of the calculation only on master
      IF (irank2 /= 0) RETURN

      ! calculate contribution from the core states to the HF exchange
      IF (xcpot%is_name("hse").OR.xcpot%is_name("vhse")) THEN
#ifdef CPP_NEVER           
         CALL exchange_vccvHSE(nk,atoms,hybrid,hybdat,dimension,jsp,lapw,nsymop,nsest,indx_sest,mpi,a_ex,results,mat_ex%core)
         CALL exchange_ccccHSE(nk,obsolete,atoms,hybdat,ncstd,kpts(:,nk),sym,a_ex,mpi,results%core)
#endif
         STOP "HSE not implemented in hsfock"
      ELSE
         CALL exchange_vccv1(nk,atoms,hybrid,hybdat,dimension,jsp,lapw,nsymop,nsest,indx_sest,mpi,a_ex,results,ex)
         CALL exchange_cccc(nk,atoms,hybdat,ncstd,sym,kpts,a_ex,mpi,results)
      END IF

      DEALLOCATE(n_q)
      CALL timestop("core exchange calculation")

      CALL timestart("time for performing T^-1*mat_ex*T^-1*")
      !calculate trafo from wavefunctions to APW basis
      IF(dimension%neigd.LT.hybrid%nbands(nk)) STOP 'mhsfock: neigd  < nbands(nk) ; '& 
                                                    'trafo from wavefunctions to APW requires at least nbands(nk) '

      call z%init(olap%l_real,nbasfcn,dimension%neigd)
      call read_z(z,kpts%nkpt*(jsp-1)+nk)
      z%matsize2 = hybrid%nbands(nk) ! reduce "visible matsize" for the following computations

      call olap%multiply(z,trafo)

      CALL invtrafo%alloc(olap%l_real,hybrid%nbands(nk),nbasfcn)
      CALL trafo%TRANSPOSE(invtrafo)
      IF(.NOT.invtrafo%l_real) invtrafo%data_c = CONJG(invtrafo%data_c)

      DO i=1,hybrid%nbands(nk)
         DO j=1,i-1
            IF (ex%l_real) THEN
               ex%data_r(i,j)=ex%data_r(j,i)
            ELSE
               ex%data_c(i,j)=conjg(ex%data_c(j,i))
            END IF
         ENDDO
      ENDDO

      WRITE(1225,'(a,i7)') 'kpoint: ', nk
      DO i = 1, ex%matsize1
         DO j = 1, i
            IF (ex%l_real) THEN
               WRITE(1225,'(2i7,2f15.8)') i, j, ex%data_r(i,j), ex%data_r(j,i)
            ELSE
               WRITE(1225,'(2i7,4f15.8)') i, j, ex%data_c(i,j), ex%data_c(j,i)
            ENDIF
         END DO
      END DO

      CALL ex%multiply(invtrafo,tmp)
      CALL trafo%multiply(tmp,v_x)
        
      CALL timestop("time for performing T^-1*mat_ex*T^-1*")

      WRITE(1231,'(a,i7)') 'kpoint: ', nk
      DO i = 1, v_x%matsize1
         DO j = 1, i
            IF (v_x%l_real) THEN
               WRITE(1231,'(2i7,1f15.8)') i, j, v_x%data_r(i,j)
            ELSE
               WRITE(1231,'(2i7,2f15.8)') i, j, v_x%data_c(i,j)
            ENDIF
         END DO
      END DO

      CALL symmetrizeh(atoms,kpts%bkf(:,nk),dimension,jsp,lapw,sym,hybdat%kveclo_eig,cell,nsymop,psym,v_x)

      WRITE(1232,'(a,i7)') 'kpoint: ', nk
      DO i = 1, v_x%matsize1
         DO j = 1, i
            IF (v_x%l_real) THEN
               WRITE(1232,'(2i7,1f15.8)') i, j, v_x%data_r(j,i) ! Note the different indices in comparison to points above. This is wanted!
            ELSE
               WRITE(1232,'(2i7,2f15.8)') i, j, v_x%data_c(j,i) ! Note the different indices in comparison to points above. This is wanted!
            ENDIF
         END DO
      END DO

      CALL write_v_x(v_x,kpts%nkpt*(jsp-1) + nk)
   END IF ! hybrid%l_calhf

   CALL timestop("total time hsfock")

END SUBROUTINE hsfock

END MODULE m_hsfock
