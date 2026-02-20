!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt
   USE m_juDFT
   USE m_constants
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho,rho_core, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo)
      USE m_dfpt_check
      USE m_dfpt_sternheimer
      USE m_dfpt_dynmat
      USE m_dfpt_eii2,     only : CalcIIEnerg2, genPertPotDensGvecs, dfpt_e2_madelung
      USE m_dfpt_dynmat_eig
      USE m_juDFT_stop, only : juDFT_error
      USE m_vgen_coulomb
      USE m_dfpt_vgen
      USE m_fleur_init
      USE m_eigen
      USE m_fermie
      USE m_grdchlh
      USE m_dfpt_dynmat_sym
      USE m_types_eigdos
      USE m_make_dos
      USE m_dfpt_gradient
      USE m_dfpt_elph_mat
      USE m_npy
      use m_inv3
      USE m_dfpt_dielecten
      USE m_dfpt_born_effcharge
      USE m_dfpt_vefield
      USE m_dfpt_potdenLocal
      USE m_dfpt_generate_gradient
      USE m_dfpt_phonon
      USE m_dfpt_borncharges
      USE m_dfpt_efield
      USE m_dfpt_interpolation
      USE m_types_BEC


      TYPE(t_mpi),        INTENT(IN)     :: fmpi
      TYPE(t_fleurinput), INTENT(IN)     :: fi
      TYPE(t_sphhar),     INTENT(IN)     :: sphhar
      TYPE(t_stars),      INTENT(IN)     :: stars
      TYPE(t_nococonv),   INTENT(IN)     :: nococonv
      TYPE(t_enpara),     INTENT(INOUT)  :: enpara
      TYPE(t_results),    INTENT(INOUT)  :: results
      TYPE(t_hybdat),     INTENT(INOUT)  :: hybdat
      TYPE(t_mpdata),     INTENT(INOUT)  :: mpdata

      CLASS(t_xcpot),     INTENT(IN)     :: xcpot
      CLASS(t_forcetheo), INTENT(INOUT)  :: forcetheo

      TYPE(t_kpts),       INTENT(IN)  :: qpts !Possibly replace this by fi%kpts [read correctly!]

      TYPE(t_potden),   INTENT(INOUT) :: rho
      TYPE(t_potden),   INTENT(INOUT) :: rho_core
      TYPE(t_potden),   INTENT(IN)    :: vTot, vxc
      INTEGER,          INTENT(IN)    :: eig_id

      TYPE(t_hub1data) :: hub1data
      TYPE(t_potden)                :: grvextdummy, imagrhodummy, vext_dummy, vC_dummy
      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: grgrVext3x3(3,3), grgrvextnum(3,3)
      TYPE(t_potden)                :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, vTot1m, vTot1mIm  ! q-quantities
      TYPE(t_potden), ALLOCATABLE   :: den_elph(:) , denIm_elph(:)
      TYPE(t_results)               :: q_results, results1, qm_results, results1m
      TYPE(t_kpts)                  :: qpts_loc
      TYPE(t_kpts)                  :: kqpts, kqmpts ! basically kpts, but with q added onto each one.
      TYPE(t_dos),TARGET            :: dos
      TYPE(t_BEC)                   :: BEC
      
      TYPE(t_eigdos_list),allocatable :: eigdos(:)

      ! Desymmetrized type variables:
      TYPE(t_fleurinput) :: fiLocal
      TYPE(t_stars)      :: starsq, starsmq


      CLASS(t_xcpot),     ALLOCATABLE :: xcpot_nosym
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo_nosym

      ! Full symmetrized type variables:
      TYPE(t_mpi)        :: fmpi_fullsym
      TYPE(t_fleurinput) :: fi_fullsym
      TYPE(t_sphhar)     :: sphhar_fullsym
      TYPE(t_stars)      :: stars_fullsym
      TYPE(t_nococonv)   :: nococonv_fullsym
      TYPE(t_enpara)     :: enpara_fullsym
      TYPE(t_results)    :: results_fullsym
      TYPE(t_wann)       :: wann_fullsym
      TYPE(t_hybdat)     :: hybdat_fullsym
      TYPE(t_mpdata)     :: mpdata_fullsym

      CLASS(t_xcpot),     ALLOCATABLE :: xcpot_fullsym
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo_fullsym

      ! starsLocal type variables
      TYPE(t_stars) :: starsLocal
      TYPE(t_potden):: imagrhodummyLocal,grvextdummyLocal
      TYPE(t_atoms) :: atomsLocal

      INTEGER,          ALLOCATABLE :: recG(:, :)
      INTEGER                       :: ngdp2km
      complex,           allocatable :: E2ndOrdII(:, :)
      complex,           allocatable :: eigenFreqs(:)
      real,              allocatable :: eigenVals(:), eigenValsFull(:,:,:)
      complex,           allocatable :: eigenVecs(:, :)

      COMPLEX, ALLOCATABLE :: grrhodummy(:, :, :, :, :),grrho_val(:, :, :, :, :)

      COMPLEX, ALLOCATABLE :: dyn_mat(:,:,:), dyn_mat_r(:,:,:), dyn_mat_q_full(:,:,:), dyn_mat_pathq(:,:), sym_dynvec(:,:,:), sym_dyn_mat(:,:,:),dyn_mat_NAC(:,:)
      REAL,    ALLOCATABLE :: e2_vm(:,:,:)

      !For e-field:
      COMPLEX, ALLOCATABLE            :: diel_tensor(:,:) !sdall i put this in an if statement?
      real                  :: qvec_int(3)
      TYPE(t_kpts)         :: qintpts
      COMPLEX,ALLOCATABLE :: born_eff_charge(:,:,:)
      COMPLEX,ALLOCATABLE :: born_eff_charge_contributions(:,:,:,:)

      INTEGER              :: q_start, q_stop
      REAL, ALLOCATABLE                 :: qvecs(:,:)

      INTEGER :: ngdp, iSpin, iQ, iDir, iDtype, nspins, zlim, iVac, lh, iDir2, sym_count
      INTEGER :: iStar, xInd, yInd, zInd, q_eig_id, ikpt, ierr, qm_eig_id, iArray
      INTEGER :: dfpt_eig_id, dfpt_eig_id2, dfpt_eigm_id, dfpt_eigm_id2
      LOGICAL :: l_real, l_minusq, l_dfpt_scf, l_cheated, l_gamma

      LOGICAL :: l_dfpt_band, l_dfpt_dos, l_dfpt_full

      CHARACTER(len=20)  :: dfpt_tag
      CHARACTER(len=4)   :: dynfiletag
      CHARACTER(len=100) :: inp_pref, trash

      INTEGER, ALLOCATABLE :: q_list(:)
      INTEGER, ALLOCATABLE :: sym_list(:) ! For each q: Collect, which symmetries leave q unchanged.

      ! Desym-tests:
      INTEGER :: grid(3), iread
      REAL    :: dr_re(fi%vacuum%nmzd), dr_im(fi%vacuum%nmzd), drr_dummy(fi%vacuum%nmzd), numbers(3*fi%atoms%nat,6*fi%atoms%nat)
      complex                           :: sigma_loc(2), sigma_ext(2), sigma_coul(2), sigma_gext(3,2), constantShift

      l_real = fi%sym%invs.AND.(.NOT.fi%noco%l_soc).AND.(.NOT.fi%noco%l_noco).AND.fi%atoms%n_hia==0

      ! l_minusq is a hard false at the moment. It can be used to ignore +-q symmetries and
      ! run the Sternheimer loop etc etc for both +q and -q.
      ! This was only used for testing but may become relevant again for SOC systems with broken
      ! inversion symmetry
      l_minusq = .FALSE.

      nspins = MERGE(2, 1, fi%noco%l_noco)
      IF (fi%juPhon%l_jpCheck) THEN
          ! This function will be used to check the validity of juPhon's
          ! input. I.e. check, whether all prohibited switches are off and,
          ! once there is more expertise considering this topic, check whether
          ! the cutoffs are chosen appropriately.
          CALL dfpt_check(fi, xcpot)
      END IF

      IF (fi%sym%nop>1) CALL juDFT_error("DFPT breaks the symmetry. Please redo the calculation without symmetry.", calledby="dfpt.F90")

      ! q_results saves the eigen-info for k+q, results1 for the perturbed quantities
      CALL q_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      CALL results1%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      
      IF (l_minusq) THEN
         CALL qm_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
         CALL results1m%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      END IF

      ! The eig_ids, where the stuff of k+q, the perturbed stuff and some extra dynmat stuff will be saved.
      q_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                        .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)
      dfpt_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                             .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
      dfpt_eig_id2 = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                              .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)

      IF (l_minusq) THEN
         qm_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                            .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)
         dfpt_eigm_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                                .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
         dfpt_eigm_id2 = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
                                .NOT.fi%INPUT%eig66(1), .FALSE., fi%noco%l_soc, fi%INPUT%eig66(1), .FALSE., fmpi%n_size)
      END IF

      ! Generate the gradients of the density and the various potentials, that will be used at different points in the programm.
      ! The density gradient is calculated by numerical differentiation, while the potential gradients are constructed (from the
      ! density gradient) by a Weinert construction, just like the potentials are from the density.
      ! This is done to ensure good continuity.
      CALL timestart("Gradient generation")
      call dfpt_generate_gradient(fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
      CALL timestop("Gradient generation")

      if (fi%juPhon%l_scf) then 
         if (fi%juPhon%l_efield) then 
            call timestart("dfpt efield")
            ! Do a scf calculation with an electric field as the external perturbation
            call dfpt_efield(fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,grRho3, &
                             grVtot3,grVext3,results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id,dfpt_eig_id2)
            call timestop("dfpt efield")
         end if 
         if (fi%juPhon%l_borneffcharge) then
            call timestart("dfpt born effective charges") 
            ! Do a scf calculation for phonon for q-Vectors close to Gamma --> calculate the polar response
            call dfpt_borncharges(fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,grRho3,grVtot3,grVC3,grVext3, &
                                    results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id,dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call timestop("dfpt born effective charges")
         end if 
         if (fi%juPhon%l_phonon) then 
            call timestart("dfpt phonons")
            ! Do a scf calculation for a phonon perturbation
            call dfpt_phonon(fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3, &
                             results,q_results,results1,eig_id,q_eig_id,dfpt_eig_id,dfpt_eig_id2,l_minusq,qm_results,results1m,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2)
            call timestop("dfpt phonons")
         end if 
      end if 

#ifdef CPP_MPI
      call MPI_BARRIER(fmpi%mpi_comm,ierr)
#endif 
      if ( fi%juPhon%l_band .or. fi%juPhon%l_dos .or. fi%juPhon%l_intp) then 
         ! Do post processing of converged results 
         ! Interpolate the dynamic matrix with FFT
         call timestart("dfpt interpolation")
         call dfpt_interpolation(fi,fmpi,nococonv,results)
         call timestop("dfpt interpolation")
      end if 


      CALL close_eig(q_eig_id)
      CALL close_eig(dfpt_eig_id)
      CALL close_eig(dfpt_eig_id2)
      IF (l_minusq) THEN
         CALL close_eig(qm_eig_id)
         CALL close_eig(dfpt_eigm_id)
         CALL close_eig(dfpt_eigm_id2)
      END IF


      if (fmpi%irank==0) WRITE (oUnit,*) '------------------------------------------------------'

    END SUBROUTINE dfpt

   SUBROUTINE dfpt_desym(fmpi_nosym,fi_nosym,sphhar_nosym,stars_nosym,nococonv_nosym,enpara_nosym,results_nosym,wann_nosym,hybdat_nosym,mpdata_nosym,xcpot_nosym,forcetheo_nosym,rho_nosym,vTot_nosym,grid,inp_pref,&
                         fi,sphhar,stars,nococonv,enpara,results,rho,vTot)
      USE m_desymmetrizer
      USE m_outcdn
      USE m_plot
      USE m_fleur_init

      TYPE(t_mpi),        INTENT(INOUT) :: fmpi_nosym
      TYPE(t_fleurinput), INTENT(INOUT) :: fi_nosym
      TYPE(t_sphhar),     INTENT(INOUT) :: sphhar_nosym
      TYPE(t_stars),      INTENT(INOUT) :: stars_nosym
      TYPE(t_nococonv),   INTENT(INOUT) :: nococonv_nosym
      TYPE(t_enpara),     INTENT(INOUT) :: enpara_nosym
      TYPE(t_results),    INTENT(INOUT) :: results_nosym
      TYPE(t_wann),       INTENT(INOUT) :: wann_nosym
      TYPE(t_hybdat),     INTENT(INOUT) :: hybdat_nosym
      TYPE(t_mpdata),     INTENT(INOUT) :: mpdata_nosym

      CLASS(t_xcpot),     ALLOCATABLE, INTENT(INOUT) :: xcpot_nosym
      CLASS(t_forcetheo), ALLOCATABLE, INTENT(INOUT) :: forcetheo_nosym

      TYPE(t_potden), INTENT(INOUT):: rho_nosym, vTot_nosym
      TYPE(t_fleurinput), INTENT(IN) :: fi
      TYPE(t_sphhar),     INTENT(IN) :: sphhar
      TYPE(t_stars),      INTENT(IN) :: stars
      TYPE(t_nococonv),   INTENT(IN) :: nococonv
      TYPE(t_enpara),     INTENT(IN) :: enpara
      TYPE(t_results),    INTENT(IN) :: results
      TYPE(t_potden),   INTENT(IN) :: rho, vTot

      INTEGER, INTENT(IN) :: grid(3)

      CHARACTER(len=100), INTENT(IN) :: inp_pref

      INTEGER :: ix, iy, iz, iv_old, iflag_old, iv_new, iflag_new
      INTEGER :: iType_old, iAtom_old, iType_new, iAtom_new!, inversionOp
      REAL    :: old_point(3), new_point(3), pt_old(3), pt_new(3), xdnout_old, xdnout_new!, atom_shift(3)
      LOGICAL :: test_desym
      CALL fleur_init(fmpi_nosym, fi_nosym, sphhar_nosym, stars_nosym, nococonv_nosym, forcetheo_nosym, &
                        enpara_nosym, xcpot_nosym, results_nosym, wann_nosym, hybdat_nosym, mpdata_nosym, &
                        inp_pref)

      CALL rho_nosym%init(stars_nosym,fi_nosym%atoms,sphhar_nosym,fi_nosym%vacuum,fi_nosym%noco,fi%input%jspins,POTDEN_TYPE_DEN)
      CALL vTot_nosym%init(stars_nosym,fi_nosym%atoms,sphhar_nosym,fi_nosym%vacuum,fi_nosym%noco,fi%input%jspins,POTDEN_TYPE_POTTOT)

      ! TODO: Correctly account for such a shift in the desymmetrization.
      ! For now: Just build input, that does not necessitate a shift.
      !        inversionOp = -1
      !        symOpLoop: DO iSym = 1, sym%nop
      !           IF (ALL(sym%mrot(:,:,iSym)==invs_matrix)) THEN
      !              inversionOp = iSym
      !              EXIT symOpLoop
      !           END IF
      !        END DO symOpLoop

      !        atom_shift = 0.0
      !        IF (inversionOp.GT.0) THEN
      !           IF(ANY(ABS(sym%tau(:,inversionOp)).GT.eps7).and..not.(film.and.ABS(sym%tau(3,inversionOp))>eps7)) THEN
      !              atom_shift = 0.5*sym%tau(:,inversionOp)
      !           END IF
      !        END IF

      ALLOCATE(vTot_nosym%pw_w, mold=vTot_nosym%pw)
      vTot_nosym%pw_w = CMPLX(0.0,0.0)

      CALL desymmetrize_pw(fi%sym, stars, stars_nosym, rho%pw, rho_nosym%pw)
      CALL desymmetrize_pw(fi%sym, stars, stars_nosym, vTot%pw, vTot_nosym%pw, vTot%pw_w, vTot_nosym%pw_w)
      CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, rho%mt, rho_nosym%mt)
      CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, vTot%mt, vTot_nosym%mt)

      CALL desymmetrize_types(fi%input, fi_nosym%input, fi%atoms, fi_nosym%atoms, fi%noco, &
                              nococonv, nococonv_nosym, enpara, enpara_nosym, results, results_nosym)

      test_desym = .FALSE.
      IF (test_desym) THEN
         DO iz = 0, grid(3)-1
            DO iy = 0, grid(2)-1
               DO ix = 0, grid(1)-1
                  old_point = fi%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi%cell%amat(:,2)*REAL(iy)/(grid(2)-1) + &
                              fi%cell%amat(:,3)*REAL(iz)/(grid(3)-1)

                  new_point = fi%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi%cell%amat(:,2)*REAL(iy)/(grid(2)-1) + &
                              fi%cell%amat(:,3)*REAL(iz)/(grid(3)-1)! - &
                              !atom_shift

                  ! Set region specific parameters for point
                  ! Get MT sphere for point if point is in MT sphere
                  CALL getMTSphere(fi%input,fi%cell,fi%atoms,old_point,iType_old,iAtom_old,pt_old)
                  CALL getMTSphere(fi_nosym%input,fi_nosym%cell,fi_nosym%atoms,new_point,iType_new,iAtom_new,pt_new)
                  IF (iAtom_old.NE.0) THEN
                     iv_old = 0
                     iflag_old = 1
                  ELSE
                     iv_old = 0
                     iflag_old = 2
                     pt_old(:) = old_point(:)
                  END IF

                  IF (iAtom_new.NE.0) THEN
                     iv_new = 0
                     iflag_new = 1
                  ELSE
                     iv_new = 0
                     iflag_new = 2
                     pt_new(:) = new_point(:)
                  END IF

                  ! Old point:
                  CALL outcdn(pt_old,iType_old,iAtom_old,iv_old,iflag_old,1,.FALSE.,stars,&
                              fi%vacuum,sphhar,fi%atoms,fi%sym,fi%cell ,&
                              rho,xdnout_old)
                  ! New point:
                  CALL outcdn(pt_new,iType_new,iAtom_new,iv_old,iflag_old,1,.FALSE.,stars_nosym,&
                              fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                              rho_nosym,xdnout_new)

                  WRITE(9004,*) xdnout_new-xdnout_old

                  ! Old point:
                  CALL outcdn(pt_old,iType_old,iAtom_old,iv_old,iflag_old,1,.TRUE.,stars,&
                              fi%vacuum,sphhar,fi%atoms,fi%sym,fi%cell ,&
                              vTot,xdnout_old)
                  ! New point:
                  CALL outcdn(pt_new,iType_new,iAtom_new,iv_old,iflag_old,1,.TRUE.,stars_nosym,&
                              fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                              vTot_nosym,xdnout_new)
                  WRITE(9005,*) xdnout_new-xdnout_old
               END DO !x-loop
            END DO !y-loop
         END DO !z-loop

!         CALL save_npy("sym_on_rhopw.npy",rho%pw)
!         CALL save_npy("sym_off_rhopw.npy",rho%pw)
!         CALL save_npy("sym_on_rhomt.npy",rho%mt)
!         CALL save_npy("sym_off_rhomt.npy",rho%mt)
!         CALL save_npy("sym_on_vpw.npy",vTot%pw)
!         CALL save_npy("sym_off_vpw.npy",vTot%pw)
!         CALL save_npy("sym_on_vmt.npy",vTot%mt)
!         CALL save_npy("sym_off_vmt.npy",vTot%mt)
      END IF
   END SUBROUTINE

   SUBROUTINE test_vac_stuff(fi_nosym,stars_nosym,sphhar_nosym,rho_nosym,vTot_nosym,grRho3,grVtot3,grVC3,grVext3,grrhodummy,grid)
      USE m_npy
      USE m_outcdn
      USE m_grdchlh
      USE m_dfpt_gradient

      TYPE(t_fleurinput), INTENT(IN) :: fi_nosym
      TYPE(t_stars),      INTENT(IN) :: stars_nosym
      TYPE(t_sphhar),     INTENT(IN) :: sphhar_nosym
      TYPE(t_potden), INTENT(IN)    :: rho_nosym, vTot_nosym, grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      INTEGER, INTENT(IN) :: grid(3)
      COMPLEX, INTENT(INOUT) :: grrhodummy(:, :, :, :, :)

      INTEGER :: ix, iy, iVac, iStar, iSpin, zlim, xInd, yInd, zInd
      REAL    :: xdnout_grrho_up_pw, xdnout_grrho_up_vac, xdnout_grrho_down_pw, xdnout_grrho_down_vac
      REAL    :: xdnout_grvc_up_pw, xdnout_grvc_up_vac, xdnout_grvc_down_pw, xdnout_grvc_down_vac
      REAL    :: point_plus(3), point_minus(3)
      REAL    :: dr_re(fi_nosym%vacuum%nmzd), dr_im(fi_nosym%vacuum%nmzd), drr_dummy(fi_nosym%vacuum%nmzd)

      COMPLEX, ALLOCATABLE :: grVtotvac(:,:,:,:), grVtotpw(:,:)

      ALLOCATE(grVtotpw(stars_nosym%ng3,3))
      ALLOCATE(grVtotvac(fi_nosym%vacuum%nmz,stars_nosym%ng2,2,3))
      DO iSpin = 1, SIZE(rho_nosym%mt,4)
         CALL mt_gradient_new(fi_nosym%atoms, sphhar_nosym, fi_nosym%sym, vTot_nosym%mt(:, :, :, iSpin), grrhodummy(:, :, :, iSpin, :))
      END DO

      DO zInd = -stars_nosym%mx3, stars_nosym%mx3
         DO yInd = -stars_nosym%mx2, stars_nosym%mx2
            DO xInd = -stars_nosym%mx1, stars_nosym%mx1
               iStar = stars_nosym%ig(xInd, yInd, zInd)
               IF (iStar.EQ.0) CYCLE
               grVtotpw(iStar,1) = vTot_nosym%pw(iStar,1) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grVtotpw(iStar,2) = vTot_nosym%pw(iStar,1) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grVtotpw(iStar,3) = vTot_nosym%pw(iStar,1) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
            END DO
         END DO
      END DO

      IF (fi_nosym%input%film) THEN
         DO yInd = -stars_nosym%mx2, stars_nosym%mx2
            DO xInd = -stars_nosym%mx1, stars_nosym%mx1
               iStar = stars_nosym%ig(xInd, yInd, 0)
               IF (iStar.EQ.0) CYCLE
               iStar = stars_nosym%ig2(iStar)
               grVtotvac(:,iStar,:,1) = vTot_nosym%vac(:,iStar,:,1) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,0]),fi_nosym%cell%bmat)))
               grVtotvac(:,iStar,:,2) = vTot_nosym%vac(:,iStar,:,1) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,0]),fi_nosym%cell%bmat)))
               DO iVac = 1, fi_nosym%vacuum%nvac
                  DO iSpin = 1, SIZE(rho_nosym%vac,4)
                     zlim = MERGE(fi_nosym%vacuum%nmz,fi_nosym%vacuum%nmzxy,iStar==1)
                     CALL grdchlh(fi_nosym%vacuum%delz, REAL(vTot_nosym%vac(:zlim,iStar,iVac,1)),dr_re(:zlim),drr_dummy(:zlim))
                     CALL grdchlh(fi_nosym%vacuum%delz,AIMAG(vTot_nosym%vac(:zlim,iStar,iVac,1)),dr_im(:zlim),drr_dummy(:zlim))
                     grVtotvac(:,iStar,iVac,3) = (3-2*iVac)*(dr_re + ImagUnit * dr_im)
                  END DO
               END DO
            END DO
         END DO
      END IF

      IF (.FALSE.) THEN!!!!! Test grRho/grVC on real space
         DO iy = 0, grid(2)-1
            DO ix = 0, grid(1)-1
               point_plus = fi_nosym%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi_nosym%cell%amat(:,2)*REAL(iy)/(grid(2)-1) + &
                              [0.0,0.0,fi_nosym%cell%z1]

               point_minus = fi_nosym%cell%amat(:,1)*REAL(ix)/(grid(1)-1) + &
                              fi_nosym%cell%amat(:,2)*REAL(iy)/(grid(2)-1) - &
                              [0.0,0.0,fi_nosym%cell%z1]! - &
                              !atom_shift
               
               ! IR rho:
               CALL outcdn(point_plus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           rho_nosym,xdnout_grrho_up_pw)
                           
               ! Vac rho:
               CALL outcdn(point_plus,1,0,1,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           rho_nosym,xdnout_grrho_up_vac)
                           
               ! IR rho:
               CALL outcdn(point_minus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           rho_nosym,xdnout_grrho_down_pw)
                           
               ! Vac rho:
               CALL outcdn(point_minus,1,0,2,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           rho_nosym,xdnout_grrho_down_vac)
                           
               write(5395,*) "Gridx/y:", ix, iy
               write(5395,*) "Upper rho:", xdnout_grrho_up_vac, xdnout_grrho_up_pw
               write(5395,*) "Lower rho:", xdnout_grrho_down_vac,xdnout_grrho_down_pw

               ! IR grrho:
               CALL outcdn(point_plus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grRho3(3),xdnout_grrho_up_pw)
               ! IR grvc:
               CALL outcdn(point_plus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grVC3(3),xdnout_grvc_up_pw)
                           
               ! IR grrho:
               CALL outcdn(point_minus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grRho3(3),xdnout_grrho_down_pw)
               ! IR grvc:
               CALL outcdn(point_minus,1,0,0,2,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grVC3(3),xdnout_grvc_down_pw)

               ! Vac grrho:
               CALL outcdn(point_plus,1,0,1,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grRho3(3),xdnout_grrho_up_vac)
               ! Vac grvc:
               CALL outcdn(point_plus,1,0,1,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grVC3(3),xdnout_grvc_up_vac)
                           
               ! Vac grrho:
               CALL outcdn(point_minus,1,0,2,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grRho3(3),xdnout_grrho_down_vac)
               ! Vac grvc:
               CALL outcdn(point_minus,1,0,2,0,1,.FALSE.,stars_nosym,&
                           fi_nosym%vacuum,sphhar_nosym,fi_nosym%atoms,fi_nosym%sym,fi_nosym%cell ,&
                           grVC3(3),xdnout_grvc_down_vac)

               write(5395,*) "Upper grrho:", xdnout_grrho_up_vac,xdnout_grrho_up_pw
               write(5395,*) "Lower grrho:", xdnout_grrho_down_vac,xdnout_grrho_down_pw
               write(5395,*) "Upper grvc: ", xdnout_grvc_up_vac,xdnout_grvc_up_pw
               write(5395,*) "Lower grvc: ", xdnout_grvc_down_vac,xdnout_grvc_down_pw
            END DO !x-loop
         END DO !y-loop   
      END IF!!!!!
      
      ! IF (fi%input%film)CALL save_npy("rhovac.npy",rho%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("rhogr1vac.npy",grRho3(1)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("rhogr2vac.npy",grRho3(2)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("rhogr3vac.npy",grRho3(3)%vac(:,:,:,1))
      ! CALL save_npy("rhogr3pw.npy",grRho3(3)%pw(:,1))
      ! IF (fi%input%film)CALL save_npy("vcgr1vac.npy",grVC3(1)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vcgr2vac.npy",grVC3(2)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vcgr3vac.npy",grVC3(3)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vextgr1vac.npy",grVext3(1)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vextgr2vac.npy",grVext3(2)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vextgr3vac.npy",grVext3(3)%vac(:,:,:,1))
      ! CALL save_npy("vextgr1pw.npy",grVext3(1)%pw(:,1))
      ! CALL save_npy("vextgr2pw.npy",grVext3(2)%pw(:,1))
      ! CALL save_npy("vextgr3pw.npy",grVext3(3)%pw(:,1))
      ! IF (fi%input%film)CALL save_npy("vtotgr1vac.npy",grVtot3(1)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vtotgr2vac.npy",grVtot3(2)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vtotgr3vac.npy",grVtot3(3)%vac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vtotgr1vacnum.npy",grVtotvac(:,:,:,1))
      ! IF (fi%input%film)CALL save_npy("vtotgr2vacnum.npy",grVtotvac(:,:,:,2))
      ! IF (fi%input%film)CALL save_npy("vtotgr3vacnum.npy",grVtotvac(:,:,:,3))
      ! CALL save_npy("vtotgr1pw.npy",grVtot3(1)%pw(:,1))
      ! CALL save_npy("vtotgr2pw.npy",grVtot3(2)%pw(:,1))
      ! CALL save_npy("vtotgr3pw.npy",grVtot3(3)%pw(:,1))
      ! CALL save_npy("vtotgr1pwnum.npy",grVtotpw(:,1))
      ! CALL save_npy("vtotgr2pwnum.npy",grVtotpw(:,2))
      ! CALL save_npy("vtotgr3pwnum.npy",grVtotpw(:,3))
   END SUBROUTINE

END MODULE m_dfpt
