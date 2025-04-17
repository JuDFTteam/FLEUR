!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnval

   USE m_juDFT
#ifdef CPP_MPI
   use mpi
#endif

CONTAINS

   SUBROUTINE cdnval(eig_id, fmpi, kpts, jspin, noco, nococonv, input, banddos, cell, atoms, enpara, stars, &
                     vacuum, sphhar, sym, vTot, cdnvalJob, den, dos, vacdos, results, &
                     moments, gfinp, hub1inp, hub1data, coreSpecInput, mcd, slab, orbcomp, jDOS, greensfImagPart)

      !************************************************************************************
      !     This is the FLEUR valence density generator
      !******** ABBREVIATIONS *************************************************************
      !     noccbd   : number of occupied bands
      !     pallst   : if set to .true. bands above the Fermi-Energy are taken into account
      !     ener     : band energy averaged over all bands and k-points,
      !                wheighted with the l-like charge of each atom type
      !     sqal     : l-like charge of each atom type. sum over all k-points and bands
      !************************************************************************************

      USE m_types
      USE m_constants
      USE m_eig66_io
      USE m_genMTBasis
      USE m_mcdinit
      USE m_sympsi
      USE m_nmat        ! calculate density matrix for LDA + U
      USE m_nmat21
      USE m_vacden
      USE m_pwden
      USE m_forcea8
      USE m_force_sf ! Klueppelberg (force level 3)
      USE m_checkdopall
      USE m_greensfBZint
      USE m_greensfCalcImagPart
      USE m_local_hamiltonian
      USE m_greensfCalcScalarProducts
      !USE m_cdnmt       ! calculate the density and orbital moments etc.
      !USE m_orbmom      ! coeffd for orbital moments
      USE m_qmtsl       ! These subroutines divide the input%film into banddos%layers
      USE m_qintsl      ! (slabs) and intergate the DOS in these banddos%layers
      USE m_orbcomp     ! calculate orbital composition (like p_x,p_y,p_z)
      USE m_jDOS
      USE m_corespec, only: l_cs    ! calculation of core spectra (EELS)
      USE m_corespec_io, only: corespec_init
      USE m_corespec_eval, only: corespec_gaunt, corespec_rme, corespec_dos, corespec_ddscs
      USE m_xmlOutput
      USE m_types_dos
      USE m_types_mcd
      USE m_types_slab
      USE m_types_jDOS
      USE m_types_vacDOS
      USE m_types_orbcomp
      USE m_denmat_to_full
      USE m_types_denmatrix
      USE m_types_radfun
      use m_l_like
      use m_types_abc
      use m_types_orbmom, only: t_orbmom
      use m_addContribsA21A12
#ifdef CPP_MPI
      USE m_mpi_col_den ! collect density data from parallel nodes
#endif
      USE m_rhonmt
      USE m_nIJmat

      IMPLICIT NONE

      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_mpi), INTENT(IN)    :: fmpi

      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_banddos), INTENT(IN)    :: banddos
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_noco), INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_stars), INTENT(IN)    :: stars
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_gfinp), INTENT(IN)    :: gfinp
      TYPE(t_hub1inp), INTENT(IN)    :: hub1inp
      TYPE(t_potden), INTENT(IN)    :: vTot
      TYPE(t_cdnvalJob), INTENT(IN)    :: cdnvalJob
      TYPE(t_potden), INTENT(INOUT) :: den
      TYPE(t_dos), INTENT(INOUT) :: dos
      TYPE(t_vacdos), INTENT(INOUT) :: vacdos
      TYPE(t_moments), INTENT(INOUT) :: moments
      TYPE(t_hub1data), OPTIONAL, INTENT(INOUT) :: hub1data
      TYPE(t_coreSpecInput), OPTIONAL, INTENT(IN)    :: coreSpecInput
      TYPE(t_mcd), OPTIONAL, INTENT(INOUT) :: mcd
      TYPE(t_slab), OPTIONAL, INTENT(INOUT) :: slab
      TYPE(t_orbcomp), OPTIONAL, INTENT(INOUT) :: orbcomp
      TYPE(t_jDOS), OPTIONAL, INTENT(INOUT) :: jDOS
      TYPE(t_greensfImagPart), OPTIONAL, INTENT(INOUT) :: greensfImagPart

      ! Scalar Arguments
      INTEGER, INTENT(IN)    :: eig_id, jspin

      ! Local Scalars
      INTEGER :: ikpt, ikpt_i, jsp_start, jsp_end, ispin, jsp, max_length_k_list, nk
      INTEGER :: iErr, nbands, noccbd, iType, ispinpr, ispin123
      INTEGER :: skip_t, skip_tt, nbasfcn,abc_itype
      LOGICAL :: l_real, l_corespec, l_empty

      ! Local Arrays
      REAL, ALLOCATABLE  :: we(:), eig(:)
      REAL                  :: bkpt(3)
      INTEGER, ALLOCATABLE  :: ev_list(:)

      TYPE(t_lapw)              :: lapw
      TYPE(t_orbmom)            :: orb
      TYPE(t_denmatrix), allocatable      :: denmatrix(:, :, :)
      TYPE(t_force)             :: force
      TYPE(t_usdus)             :: usdus
      TYPE(t_mat)               :: zMat
      TYPE(t_tlmplm)            :: tlmplm
      TYPE(t_greensfBZintCoeffs):: greensfBZintCoeffs
      TYPE(t_scalarGF), ALLOCATABLE :: scalarGF(:)
      TYPE(t_radfun)             :: radfun
      TYPE(t_abc), allocatable    :: abc(:, :)

      CALL timestart("cdnval")

      call timestart("init")
      l_real = sym%invs .AND. (.NOT. noco%l_soc) .AND. (.NOT. noco%l_noco) .AND. atoms%n_hia == 0

      ! Klueppelberg (force level 3)
      IF (input%l_f .AND. (input%f_level .GE. 3)) THEN
         CALL init_sf(sym, cell, atoms)
      END IF

      IF (noco%l_mperp .OR. banddos%l_jDOS) THEN
         ! when the off-diag. part of the density matrix, i.e. m_x and
         ! m_y, is calculated inside the muffin-tins (l_mperp = T), cdnval
         ! is called only once. therefore, several spin loops have been
         ! added. if l_mperp = F, these loops run only from jspin - jspin.
         jsp_start = 1
         jsp_end = 2
      ELSE
         jsp_start = jspin
         jsp_end = jspin
      END IF

      allocate (denmatrix(jsp_start:jsp_end, jsp_start:jsp_end, atoms%ntype))
      DO ispin = jsp_start, jsp_end
         DO jsp = jsp_start, jsp_end
            DO itype = 1, atoms%ntype
               call denmatrix(ispin, jsp, itype)%init(itype, atoms, input, sphhar)
            end do
         end do
      end do
      allocate (abc(jsp_start:jsp_end, merge(1,atoms%ntype,atoms%n_v==0)))

      !Do we need to consider the unoccupied states
      l_empty = banddos%dos .or. banddos%band
      IF (gfinp%n > 0 .AND. PRESENT(greensfImagPart)) THEN
         l_empty = l_empty .OR. greensfImagPart%l_calc
      END IF

!

      ! Initializations
      CALL usdus%init(atoms, input%jspins)
      CALL force%init1(input, atoms)
      !CALL orb%init(atoms,noco,jsp_start,jsp_end)

      !Greens function always considers the empty states
      IF (gfinp%n > 0 .AND. PRESENT(greensfImagPart)) THEN
         IF (greensfImagPart%l_calc) THEN
            CALL greensfBZintCoeffs%init(gfinp, atoms, noco, SIZE(cdnvalJob%ev_list))
            CALL greensfCalcScalarProducts(gfinp, atoms, input, enpara, noco, sphhar, vTot, fmpi, hub1data=hub1data, &
                                           scalarProducts=scalarGF)
         END IF
      END IF

      IF (banddos%l_mcd .AND. .NOT. PRESENT(mcd)) CALL juDFT_error("mcd is missing", calledby="cdnval")

      ! calculation of core spectra (EELS) initializations -start-
      l_coreSpec = .FALSE.
      IF (PRESENT(coreSpecInput)) THEN
         CALL corespec_init(input, atoms, coreSpecInput)
         IF (l_cs .AND. (fmpi%isize .NE. 1)) CALL juDFT_error('EELS + fmpi not implemented', calledby='cdnval')
         IF (l_cs .AND. jspin .EQ. 1) CALL corespec_gaunt()
         l_coreSpec = l_cs
      END IF
      ! calculation of core spectra (EELS) initializations -end-

      IF (fmpi%irank == 0) THEN
         WRITE (oUnit, FMT=8000) jspin
         CALL openXMLElementPoly('mtCharges', (/'spin'/), (/jspin/))
      END IF
8000  FORMAT(/, /, 10x, 'valence density: spin=', i2)
      BLOCK !TODO This block should be put into subroutine?
         REAL, ALLOCATABLE  :: f(:, :, :, :), g(:, :, :, :), flo(:, :, :, :) ! radial functions
         ALLOCATE (f(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
         ALLOCATE (g(atoms%jmtd, 2, 0:atoms%lmaxd, input%jspins))
         ALLOCATE (flo(atoms%jmtd, 2, atoms%nlod, input%jspins))
         DO iType = 1, atoms%ntype

            DO ispin = 1, input%jspins
               CALL genMTBasis(atoms, enpara, vTot, fmpi, iType, ispin, usdus, f(:, :, 0:, ispin), g(:, :, 0:, ispin), &
                               flo(:, :, :, ispin), hub1data=hub1data)

            END DO
            IF (banddos%l_mcd) CALL mcd_init(atoms, banddos, input, vTot%mt(:, 0, :, :), g, f, mcd, iType, jspin)
            IF (l_coreSpec) CALL corespec_rme(atoms, input, iType, 29, input%jspins, jspin, results%ef, &
                                              atoms%msh, vTot%mt(:, 0, :, :), f, g)

         END DO
         DEALLOCATE (f, g, flo)
      end block

      skip_tt = dot_product(enpara%skiplo(:atoms%ntype, jspin), atoms%neq(:atoms%ntype))
      IF (noco%l_soc .OR. noco%l_noco) skip_tt = 2*skip_tt

      jsp = MERGE(1, jspin, noco%l_noco)
      call timestop("init")

      max_length_k_list = size(cdnvalJob%k_list)
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, max_length_k_list, 1, MPI_INTEGER, MPI_MAX, fmpi%mpi_comm, ierr)
#endif
      DO ikpt_i = 1, size(cdnvalJob%k_list)
         ikpt = cdnvalJob%k_list(ikpt_i)
         bkpt = kpts%bk(:, ikpt)

         CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
         skip_t = skip_tt
         ev_list = cdnvaljob%compact_ev_list(ikpt_i, l_empty)
         noccbd = SIZE(ev_list)
         we = cdnvalJob%weights(ev_list, ikpt)
         eig = results%eig(ev_list, ikpt, jsp)

         IF (cdnvalJob%l_evp) THEN
            IF (minval(ev_list) > skip_tt) skip_t = 0
            IF (maxval(ev_list) <= skip_tt) skip_t = noccbd
            IF ((minval(ev_list) <= skip_tt) .AND. (maxval(ev_list) > skip_tt)) skip_t = mod(skip_tt, noccbd)
         END IF

         nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot, lapw%nv(1) + atoms%nlotot, noco%l_noco)
         CALL zMat%init(l_real, nbasfcn, noccbd)
         CALL read_eig(eig_id, ikpt, jsp, list=ev_list, neig=nbands, zmat=zMat)
#ifdef CPP_MPI
         CALL MPI_BARRIER(fmpi%mpi_comm, iErr) ! Synchronizes the RMA operations
#endif

         IF (noccbd .LE. 0) CYCLE ! Note: This jump has to be after the MPI_BARRIER is called

         ! valence density in the atomic spheres
         DO itype = 1, atoms%ntype
            abc_itype=min(itype,size(abc,2)) !abc might be only needed for a single itype
            call radfun%generate_radial_functions(atoms, input, enpara, hub1data, fmpi, vtot, iType)
            DO ispin = jsp_start, jsp_end
               IF (input%l_f) CALL force%init2(noccbd, input, atoms)
               call abc(ispin, abc_itype)%init(input, atoms, radfun, noccbd, itype)
               call abc(ispin, abc_itype)%calc_abc(input, atoms, sym, cell, lapw, noccbd, usdus, noco, nococonv, ispin, itype, zMat)
               DO ispinpr = jsp_start, ispin
                  ispin123 = merge(ispin, 3, ispin == ispinpr) !sometimes the "3rd" spin is the off-diagonal part
                  !Calculate the density matrix for LDA+U and related methods
                  IF (atoms%n_u + atoms%n_opc .GT. 0) THEN
                     CALL n_mat(atoms, radfun, sym, noccbd, we, abc(ispin, abc_itype), &
                                abc(ispinpr, abc_itype), den%mmpMat(:, :, :, ispin123), itype, ispin, ispinpr)
                  END IF

                  ! Determine weights for DOS and Bandstructures
                  call dos%calc_mt_dos(abc(ispin, abc_itype), abc(ispinpr, abc_itype), banddos, radfun, &
                                       atoms, ev_list, itype, ikpt, ispin, ispinpr)
                  if (ispin == ispinpr) THEN
                     !No off-diagonal contributions yet
                     call mcd%calc_mt_mcd(banddos, atoms, ev_list, abc(ispin, abc_itype), itype, ikpt, ispin)
                     IF (banddos%l_orb) &
                        call orbcomp%calc_orb_comp(atoms, banddos, radfun, abc(ispin, abc_itype), abc(ispin, abc_itype), &
                                                   ev_list, itype, ikpt, ispin, ispin)
                     !IF(l_coreSpec) CALL corespec_dos(atoms,usdus,ispin,atoms%lmaxd*(atoms%lmaxd+2),kpts%nkpt,ikpt,input%neig,&
                     !                                 noccbd,results%ef,banddos%sig_dos,eig,we,eigVecCoeffs) !TODO
                  end if
                  !Decomposition into total angular momentum states
                  IF (banddos%dos .AND. banddos%l_jDOS) THEN
                     IF (PRESENT(jDOS) .AND. ispinpr == jsp_end) THEN
                        call jDOS%calc_jDOS(ikpt, noccbd, ev_list, we, atoms, banddos, input, radfun, abc(1, abc_itype), abc(2, abc_itype))
                     END IF
                  END IF

                  IF (noco%l_soc .and. ispin == ispinpr) CALL orb%calc_orbmom(abc(ispin, abc_itype), atoms, radfun, we, itype, &
                                                                              ispin, moments%clmom(:, itype, ispin))  !TODO MPI collect is missing here

                  !Now calculate the density matrix as needed to construct the charge
                  call denmatrix(ispin, ispinpr, itype)%rhonmt(atoms, sphhar, we, noccbd, itype, &
                                                               sym, abc(ispin, abc_itype), abc(ispinpr, abc_itype))
                  !call denmatrix(ispin,ispin,itype)%mpi_collect() TODO
               end do !loop over ispinpr
               IF (input%l_f) THEN
                  !Calculate force contributions
                  call abc(ispin, abc_itype)%calc_force_abc(input, atoms, sym, cell, lapw, &
                                                        noccbd, usdus, noco, nococonv, ispin, itype, zmat, eig, force)

                  call local_ham(sphhar, atoms, sym, noco, nococonv, enpara, fmpi, vtot, &
                                 vtot, den, input, hub1inp, hub1data, tlmplm, usdus, 0.0)
                  CALL addContribsA21A12(force, input, atoms, sym, cell, enpara, &
                        usdus, tlmplm, vtot, abc(ispin,abc_itype), noccbd, ispin, eig, we, results, jsp_start, jspin, nbasfcn, zMat, lapw, &
                                         sphhar, lapw%gvec(1, :, :), lapw%gvec(2, :, :), lapw%gvec(3, :, :), bkpt, itype)
               END IF

            END DO ! end loop over ispin
         END DO !loop over itypes
         IF (atoms%n_v.GT.0) CALL nIJ_mat(lbound(abc,2),input,atoms,noccbd,usdus,we,abc,cell,kpts,ikpt,den%nIJ_llp_mmp,enpara,vTot) 

         ! valence density in the interstitial and vacuum region has to be called only once (if jspin=1) in the non-collinear case
         IF (.NOT. ((jspin .EQ. 2) .AND. noco%l_noco)) THEN
            ! valence density in the interstitial region
            CALL pwden(stars, kpts, banddos, input, fmpi, noco, nococonv, cell, atoms, sym, ikpt, &
                       jspin, lapw, noccbd, ev_list, we, eig, den, results, force%f_b8, zMat, dos)
            ! charge of each valence state in this k-point of the SBZ in the layer interstitial region of the film
            IF (PRESENT(slab) .AND. banddos%l_slab) &
               CALL q_int_sl(jspin, ikpt, stars, atoms, sym, cell, noccbd, ev_list, lapw, slab, zMat)
            ! valence density in the vacuum region
            IF (input%film) THEN
               CALL vacden(vacuum, stars, input, cell, atoms, noco, nococonv, banddos, &
                           we, ikpt, jspin, REAL(vTot%vac(:, 1, :, :)), noccbd, ev_list, lapw, enpara%evac, den, zMat, vacdos, dos)
            END IF
         END IF

         IF (input%l_sympsi .and. allocated(dos%jsym)) THEN
            CALL sympsi(lapw, jspin, sym, noccbd, cell, eig, noco, dos%jsym(:, ikpt, jspin), zMat)
         END IF
      END DO ! end of k-point loop

#ifdef CPP_MPI
      !print *,"Remaining Barriers:",size(cdnvalJob%k_list)+1,max_length_k_list
      DO nk = size(cdnvalJob%k_list) + 1, max_length_k_list
         CALL MPI_BARRIER(fmpi%MPI_COMM, ierr)
      END DO
      DO ispin = jsp_start, jsp_end
         CALL mpi_col_den(fmpi, sphhar, atoms, stars, vacuum, input, noco, ispin, dos, vacdos, &
                          results, den, mcd, slab, orbcomp, jDOS)
      END DO
#endif

      IF (gfinp%n > 0 .AND. PRESENT(greensfImagPart)) THEN
         IF (greensfImagPart%l_calc) THEN
            call timestart("Green's function: Imag Part collect")
            do ispin = MERGE(1, jsp_start, gfinp%l_mperp), MERGE(3, jsp_end, gfinp%l_mperp)
               CALL greensfImagPart%collect(ispin, fmpi%mpi_comm)
            end do
            call timestop("Green's function: Imag Part collect")
         END IF
      END IF

      IF (fmpi%irank == 0) THEN
         CALL timestart("denmatrix_to_full")
         DO itype = 1, atoms%ntype
            call radfun%generate_radial_functions(atoms, input, enpara, hub1data, fmpi, vtot, iType)
            call print_l_like_charge(lbound(denmatrix, 1), atoms, radfun, denmatrix, itype)
            DO ispin = jsp_start, jsp_end
               do ispinpr=jsp_start, ispin
                  call denmatrix(ispin,ispinpr,itype)%to_full_density(ispin,ispinpr, itype, input, &
                                           sphhar, atoms, noco, sym, radfun,  den%mt, moments=moments)
               enddo
            enddo
         END DO
         CALL timestop("denmatrix_to_full")
         IF (l_coreSpec) CALL corespec_ddscs(jspin, input%jspins)
         DO ispin = jsp_start, jsp_end
            IF (input%cdinf) THEN
               WRITE (oUnit, FMT=8210) ispin
8210           FORMAT(/, 5x, 'check continuity of cdn for spin=', i2)
               CALL checkDOPAll(input, sphhar, stars, atoms, sym, vacuum, cell, den, ispin)
            END IF
            IF (input%l_f) CALL force_a8(input, atoms, sym, sphhar, ispin, vTot%mt(:, :, :, ispin), den%mt, force, fmpi, results)
         END DO
         CALL closeXMLElement('mtCharges')
      END IF

      CALL timestop("cdnval")

   END SUBROUTINE cdnval

!TODO:
! nstm3? functionality in gvacmap
!
! layer charge of each valence state in this k-point of the SBZ from the mt-sphere region of the film
   !IF (PRESENT(slab).and.banddos%l_slab) CALL q_mt_sl(ispin,atoms,sym,noccbd,ev_list,ikpt,noccbd,skip_t,noccbd,eigVecCoeffs,usdus,slab) !TODO
!
!
   !TODO: this has to be added again...
   !IF(gfinp%n>0 .AND. PRESENT(greensfImagPart)) THEN
   !IF(greensfImagPart%l_calc) THEN
   !do ispin = MERGE(1,jsp_start,gfinp%l_mperp),MERGE(3,jsp_end,gfinp%l_mperp)
   !  CALL greensfBZint(ikpt,noccbd,ispin,gfinp,sym,atoms,noco,nococonv,input,kpts,&
   !                    scalarGF,eigVecCoeffs,greensfBZintCoeffs) !TODO
   !  CALL greensfCalcImagPart_single_kpt(ikpt,ikpt_i,ev_list,ispin,gfinp,atoms,input,kpts,noco,fmpi,&
   !                    results,greensfBZintCoeffs,greensfImagPart)
   !enddo
   !ENDIF
   !ENDIF

END MODULE m_cdnval
