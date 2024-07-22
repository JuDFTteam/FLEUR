!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
                 & rho, vTot, vxc, eig_id, xcpot, hybdat, mpdata, forcetheo)
      USE m_dfpt_check
      !USE m_dfpt_test
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
      USE m_npy
      USE m_dfpt_dielecten

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

      TYPE(t_kpts),       INTENT(IN)  :: qpts !Possibly replace this by fi_nosym%kpts [read correctly!]

      TYPE(t_potden),   INTENT(INOUT) :: rho
      TYPE(t_potden),   INTENT(IN)    :: vTot, vxc
      INTEGER,          INTENT(IN)    :: eig_id

      TYPE(t_hub1data) :: hub1data
      TYPE(t_potden)                :: grvextdummy, imagrhodummy, rho_nosym, vTot_nosym, vext_dummy, vC_dummy
      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: grgrVC3x3(3,3), grgrvextnum(3,3)
      TYPE(t_potden)                :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, vTot1m, vTot1mIm ! q-quantities
      TYPE(t_results)               :: q_results, results1, qm_results, results1m
      TYPE(t_kpts)                  :: qpts_loc
      TYPE(t_kpts)                  :: kqpts, kqmpts ! basically kpts, but with q added onto each one.
      TYPE(t_dos),TARGET            :: dos
      
      TYPE(t_eigdos_list),allocatable :: eigdos(:)

      ! Desymmetrized type variables:
      TYPE(t_mpi)        :: fmpi_nosym
      TYPE(t_fleurinput) :: fi_nosym
      TYPE(t_sphhar)     :: sphhar_nosym
      TYPE(t_stars)      :: stars_nosym, starsq, starsmq
      TYPE(t_nococonv)   :: nococonv_nosym
      TYPE(t_enpara)     :: enpara_nosym
      TYPE(t_results)    :: results_nosym
      TYPE(t_wann)       :: wann_nosym
      TYPE(t_hybdat)     :: hybdat_nosym
      TYPE(t_mpdata)     :: mpdata_nosym

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

      INTEGER,          ALLOCATABLE :: recG(:, :)
      INTEGER                       :: ngdp2km
      complex,           allocatable :: E2ndOrdII(:, :)
      complex,           allocatable :: eigenFreqs(:)
      real,              allocatable :: eigenVals(:), eigenValsFull(:,:,:)
      complex,           allocatable :: eigenVecs(:, :)

      COMPLEX, ALLOCATABLE :: grrhodummy(:, :, :, :, :)

      COMPLEX, ALLOCATABLE :: dyn_mat(:,:,:), dyn_mat_r(:,:,:), dyn_mat_q_full(:,:,:), dyn_mat_pathq(:,:), sym_dynvec(:,:,:), sym_dyn_mat(:,:,:)
      REAL,    ALLOCATABLE :: e2_vm(:,:,:)

      !For e-field:
      COMPLEX, ALLOCATABLE :: diel_tensor(:,:)
      REAL, ALLOCATABLE    :: diel_tensor_occ1(:,:)!check if this makes sense
      REAL, ALLOCATABLE    :: we1_data(:,:,:,:),eig1_data(:,:,:,:)
      COMPLEX, ALLOCATABLE :: dielten_iden(:,:)   
      INTEGER              :: i_iden

      INTEGER :: ngdp, iSpin, iQ, iDir, iDtype, nspins, zlim, iVac, lh, iDir2, sym_count
      INTEGER :: iStar, xInd, yInd, zInd, q_eig_id, ikpt, ierr, qm_eig_id, iArray
      INTEGER :: dfpt_eig_id, dfpt_eig_id2, dfpt_eigm_id, dfpt_eigm_id2
      LOGICAL :: l_real, l_minusq, l_dfpt_scf, l_cheated

      LOGICAL :: l_dfpt_band, l_dfpt_dos, l_dfpt_full

      CHARACTER(len=20)  :: dfpt_tag
      CHARACTER(len=4)   :: dynfiletag
      CHARACTER(len=100) :: inp_pref, trash

      INTEGER, ALLOCATABLE :: q_list(:)
      INTEGER, ALLOCATABLE :: sym_list(:) ! For each q: Collect, which symmetries leave q unchanged.

      ! Desym-tests:
      INTEGER :: grid(3), iread
      REAL    :: dr_re(fi%vacuum%nmzd), dr_im(fi%vacuum%nmzd), drr_dummy(fi%vacuum%nmzd), numbers(3*fi%atoms%nat,6*fi%atoms%nat)
      complex                           :: sigma_loc(2), sigma_ext(2), sigma_coul(2), sigma_gext(3,2)



      
      ALLOCATE(e2_vm(fi%atoms%nat,3,3))

      l_dfpt_scf   = fi%juPhon%l_scf

      l_dfpt_band  = fi%juPhon%l_band
      l_dfpt_full  = fi%juPhon%l_intp
      l_dfpt_dos   = fi%juPhon%l_dos

      l_cheated = .FALSE.

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
          CALL dfpt_check(fi_nosym, xcpot_nosym)
      END IF
      !print*,"fi%atoms%ntype",fi%atoms%ntype
      !STOP
      IF (fi%sym%nop>1) THEN
         WRITE(*,*) "Desymmetrization needed. Going ahead!"
         ! Grid size for desym quality test:
         grid = 21

         inp_pref = ADJUSTL("desym_")
         fmpi_nosym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
         fmpi_nosym%mpi_comm = fmpi%mpi_comm

         CALL dfpt_desym(fmpi_nosym,fi_nosym,sphhar_nosym,stars_nosym,nococonv_nosym,enpara_nosym,results_nosym,wann_nosym,hybdat_nosym,mpdata_nosym,xcpot_nosym,forcetheo_nosym,rho_nosym,vTot_nosym,grid,inp_pref,&
                         fi,sphhar,stars,nococonv,enpara,results,rho,vTot)
      ELSE
         fmpi_nosym      = fmpi
         fi_nosym        = fi
         sphhar_nosym    = sphhar
         stars_nosym     = stars
         nococonv_nosym  = nococonv
         forcetheo_nosym = forcetheo
         enpara_nosym    = enpara
         xcpot_nosym     = xcpot
         results_nosym   = results
         hybdat_nosym    = hybdat
         mpdata_nosym    = mpdata
         rho_nosym       = rho
         vTot_nosym      = vTot
      END IF

      ! TODO: Maybe rather replace this by switches filtered into dfpt_routines.
      ! IF (fi%juPhon%l_jpTest) THEN
           ! This function will be used to run (parts of) the test suite for
           ! OG juPhon, as provided by CRG.
           !CALL dfpt_test(fi, sphhar, stars, fmpi, rho, grRho, rho0, grRho0, xcpot, ngdp, recG, grVxcIRKern, ylm, dKernMTGPts, gausWts, hybdat)
      !END IF

      ! q_results saves the eigen-info for k+q, results1 for the perturbed quantities
      CALL q_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      CALL results1%init(fi_nosym%input, fi_nosym%atoms, fi_nosym%kpts, fi_nosym%noco)
      
      IF (l_minusq) THEN
         CALL qm_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
         CALL results1m%init(fi_nosym%input, fi_nosym%atoms, fi_nosym%kpts, fi_nosym%noco)
      END IF

      IF (.NOT.fi%juPhon%qmode==0) THEN
         ! Read qpoints from fullsym_inp.xml and fullsym_kpts.xml
         inp_pref = ADJUSTL("fullsym_")
         fmpi_fullsym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
         fmpi_fullsym%mpi_comm = fmpi%mpi_comm
         CALL fleur_init(fmpi_fullsym, fi_fullsym, sphhar_fullsym, stars_fullsym, nococonv_fullsym, forcetheo_fullsym, &
                         enpara_fullsym, xcpot_fullsym, results_fullsym, wann_fullsym, hybdat_fullsym, mpdata_fullsym, &
                         inp_pref)
         qpts_loc = fi_fullsym%kpts

         ALLOCATE(q_list(SIZE(qpts_loc%bk,2)))
         q_list = (/(iArray, iArray=1,SIZE(qpts_loc%bk,2), 1)/)

         ALLOCATE(sym_list(fi_fullsym%sym%nop))
         sym_list = 0
      ELSE
         ! Read qpoints from the juPhon qlist in inp.xml
         qpts_loc = qpts
         qpts_loc%bk(:, :SIZE(fi%juPhon%qvec,2)) = fi%juPhon%qvec

         ALLOCATE(q_list(SIZE(fi%juPhon%qvec,2)))
         q_list = (/(iArray, iArray=1,SIZE(fi%juPhon%qvec,2), 1)/)
      END IF

      ! Generate the gradients of the density and the various potentials, that will be used at different points in the programm.
      ! The density gradient is calculated by numerical differentiation, while the potential gradients are constructed (from the
      ! density gradient) by a Weinert construction, just like the potentials are from the density.
      ! This is done to ensure good continuity.
      CALL timestart("Gradient generation")
      ALLOCATE(grrhodummy(fi_nosym%atoms%jmtd, (fi_nosym%atoms%lmaxd+1)**2, fi_nosym%atoms%nat, SIZE(rho_nosym%mt,4), 3))

      CALL imagrhodummy%copyPotDen(rho_nosym)
      CALL imagrhodummy%resetPotDen()
      CALL grvextdummy%copyPotDen(rho_nosym)

      CALL vext_dummy%copyPotDen(vTot_nosym)
      CALL vext_dummy%resetPotDen()
      CALL vC_dummy%copyPotDen(vTot_nosym)
      CALL vC_dummy%resetPotDen()
      sigma_loc  = cmplx(0.0,0.0)
      sigma_ext  = cmplx(0.0,0.0)
      sigma_coul = cmplx(0.0,0.0)
      sigma_gext = cmplx(0.0,0.0)
      CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, fi%juphon, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .FALSE., imagrhodummy, vext_dummy, sigma_ext)
      CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, fi%juphon, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .FALSE., rho_nosym, vC_dummy, sigma_coul)
      DO iDir = 1, 3
         CALL grRho3(iDir)%copyPotDen(rho_nosym)
         CALL grRho3(iDir)%resetPotDen()
         DO iDir2 = 1, 3
            CALL grgrvextnum(iDir2,iDir)%copyPotDen(vTot_nosym)
            CALL grgrvextnum(iDir2,iDir)%resetPotDen()
            CALL grgrVC3x3(iDir2,iDir)%copyPotDen(vTot_nosym)
            CALL grgrVC3x3(iDir2,iDir)%resetPotDen()
         END DO
         CALL grVext3(iDir)%copyPotDen(vTot_nosym)
         CALL grVext3(iDir)%resetPotDen()
         CALL grVtot3(iDir)%copyPotDen(vTot_nosym)
         CALL grVtot3(iDir)%resetPotDen()
         CALL grVC3(iDir)%copyPotDen(vTot_nosym)
         CALL grVC3(iDir)%resetPotDen()
         ! Generate the external potential gradient.
         write(oUnit, *) "grVext", iDir
         sigma_loc  = cmplx(0.0,0.0)
         !IF (iDir==3) sigma_loc  = sigma_ext
         CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, fi%juphon, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .FALSE., imagrhodummy, grVext3(iDir), sigma_loc, &
                         & dfptdenimag=imagrhodummy, dfptvCoulimag=grvextdummy,dfptden0=imagrhodummy,stars2=stars_nosym,iDtype=0,iDir=iDir)
         IF (iDir==3) sigma_gext(iDir,:) = sigma_loc
      END DO
      !CALL vext_dummy%copyPotDen(vTot_nosym)
      !CALL vext_dummy%resetPotDen()
      ! Density gradient
      DO iSpin = 1, SIZE(rho_nosym%mt,4)
         CALL mt_gradient_new(fi_nosym%atoms, sphhar_nosym, fi_nosym%sym, rho_nosym%mt(:, :, :, iSpin), grrhodummy(:, :, :, iSpin, :))
      END DO
      DO zInd = -stars_nosym%mx3, stars_nosym%mx3
         DO yInd = -stars_nosym%mx2, stars_nosym%mx2
            DO xInd = -stars_nosym%mx1, stars_nosym%mx1
               iStar = stars_nosym%ig(xInd, yInd, zInd)
               IF (iStar.EQ.0) CYCLE
               grRho3(1)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grRho3(2)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grRho3(3)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(1,1)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(1,2)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(1,3)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(2,1)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(2,2)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(2,3)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(3,1)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(3,2)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grgrvextnum(3,3)%pw(iStar,:) = vext_dummy%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat))) &
                                                                   * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
            END DO
         END DO
      END DO

      IF (fi_nosym%input%film) THEN
         DO yInd = -stars_nosym%mx2, stars_nosym%mx2
            DO xInd = -stars_nosym%mx1, stars_nosym%mx1
               iStar = stars_nosym%ig(xInd, yInd, 0)
               IF (iStar.EQ.0) CYCLE
               iStar = stars_nosym%ig2(iStar)
               grRho3(1)%vac(:,iStar,:,:) = rho_nosym%vac(:,iStar,:,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,0]),fi_nosym%cell%bmat)))
               grRho3(2)%vac(:,iStar,:,:) = rho_nosym%vac(:,iStar,:,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,0]),fi_nosym%cell%bmat)))
               DO iVac = 1, fi_nosym%vacuum%nvac
                  DO iSpin = 1, SIZE(rho_nosym%vac,4)
                     zlim = MERGE(fi_nosym%vacuum%nmz,fi_nosym%vacuum%nmzxy,iStar==1)
                     CALL grdchlh(fi_nosym%vacuum%delz, REAL(rho_nosym%vac(:zlim,iStar,iVac,iSpin)),dr_re(:zlim),drr_dummy(:zlim))
                     CALL grdchlh(fi_nosym%vacuum%delz,AIMAG(rho_nosym%vac(:zlim,iStar,iVac,iSpin)),dr_im(:zlim),drr_dummy(:zlim))
                     grRho3(3)%vac(:,iStar,iVac,iSpin) = (3-2*iVac)*(dr_re + ImagUnit * dr_im)
                  END DO
               END DO
            END DO
         END DO
      END IF

      ! Coulomb/Effective potential gradients
      DO iDir = 1, 3
         CALL sh_to_lh(fi_nosym%sym, fi_nosym%atoms, sphhar_nosym, SIZE(rho_nosym%mt,4), 2, grrhodummy(:, :, :, :, iDir), grRho3(iDir)%mt, imagrhodummy%mt)
         CALL imagrhodummy%resetPotDen()
         write(oUnit, *) "grVeff", iDir
         sigma_loc  = cmplx(0.0,0.0)
         IF (iDir==3) sigma_loc  = sigma_coul
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi%juphon, fi_nosym%cell, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVtot3(iDir), .TRUE., grvextdummy, grRho3(iDir), 0, iDir, [0,0], sigma_loc)
         write(oUnit, *) "grVC", iDir
         sigma_loc  = cmplx(0.0,0.0)
         IF (iDir==3) sigma_loc  = sigma_coul
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi%juphon, fi_nosym%cell, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVC3(iDir), .FALSE., grvextdummy, grRho3(iDir), 0, iDir, [0,0], sigma_loc)
      END DO

         DO iDir2 = 1, 3
            DO iDir = 1, 3
               CALL imagrhodummy%resetPotDen()
               sigma_loc = cmplx(0.0,0.0)

               !IF (iDir2==3) sigma_loc = sigma_gext(iDir,:)
               !IF (iDir==3) sigma_loc = sigma_gext(iDir2,:)
               CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, fi%juphon, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .TRUE., imagrhodummy, grgrVC3x3(iDir2,iDir), sigma_loc, &
                         & dfptdenimag=imagrhodummy, dfptvCoulimag=grvextdummy,dfptden0=imagrhodummy,stars2=stars_nosym,iDtype=0,iDir=iDir,iDir2=iDir2, &
                         & sigma_disc2=MERGE(sigma_ext,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir2==3.AND.iDir==3.AND..FALSE.))
               CALL dfpt_e2_madelung(fi_nosym%atoms,fi_nosym%input%jspins,imagrhodummy%mt(:,0,:,:),grgrVC3x3(iDir2,iDir)%mt(:,0,:,1),e2_vm(:,iDir2,iDir))
            END DO
         END DO
      CALL save_npy("radii.npy",fi_nosym%atoms%rmsh(:,1))
      DO iDir2 = 1, 3
         DO iDir = 1, 3
            CALL save_npy("grgrVC_"//int2str(idir2)//int2str(idir)//"_pw.npy",grgrVC3x3(iDir2,iDir)%pw(:,1))
            CALL save_npy("grgrVCnum_"//int2str(idir2)//int2str(idir)//"_pw.npy",grgrvextnum(iDir2,iDir)%pw(:,1))
         END DO
      END DO
      
      CALL grRho3(1)%distribute(fmpi%mpi_comm)
      CALL grRho3(2)%distribute(fmpi%mpi_comm)
      CALL grRho3(3)%distribute(fmpi%mpi_comm)
      CALL grVext3(1)%distribute(fmpi%mpi_comm)
      CALL grVext3(2)%distribute(fmpi%mpi_comm)
      CALL grVext3(3)%distribute(fmpi%mpi_comm)
      CALL timestop("Gradient generation")
      
      CALL test_vac_stuff(fi_nosym,stars_nosym,sphhar_nosym,rho_nosym,vTot_nosym,grRho3,grVtot3,grVC3,grVext3,grrhodummy,grid)

      ! Old CRG-jp Routine to get the vectors G+q for Eii2
      CALL genPertPotDensGvecs( stars_nosym, fi_nosym%cell, fi_nosym%input, ngdp, ngdp2km, [0.0,0.0,0.0], recG )

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

      ALLOCATE(dyn_mat(SIZE(q_list),3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      dyn_mat = cmplx(0.0,0.0)
      ALLOCATE(sym_dyn_mat(SIZE(q_list),3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      sym_dyn_mat = cmplx(0.0,0.0)
      !print*,"dynshape", shape(dyn_mat)
      !print*, "dynmat",dyn_mat
         
      
      !allocations for efield:
      ALLOCATE(we1_data(fi%input%neig,size(fmpi%k_list), MERGE(1,fi%input%jspins,fi%noco%l_noco),3*fi%atoms%ntype))
      ALLOCATE(eig1_data(fi%input%neig,size(fmpi%k_list), MERGE(1,fi%input%jspins,fi%noco%l_noco),3*fi%atoms%ntype))
      we1_data = 0.0
      eig1_data = 0.0
      
      ALLOCATE(dielten_iden(3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      dielten_iden = CMPLX(0,0)

      DO i_iden = 1,3*fi_nosym%atoms%ntype
         print*,i_iden
         dielten_iden(i_iden,i_iden) = CMPLX(1,0)
      END DO
      !print*,"dielten_iden",dielten_iden(:,:)
      !print*,"here"
      !allocate dielectric tensor:
      ALLOCATE(diel_tensor(3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      ALLOCATE(diel_tensor_occ1(3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      diel_tensor = cmplx(0.0,0.0)
      diel_tensor_occ1 = cmplx(0.0,0.0)
      IF (l_dfpt_scf) THEN
         ! Do the self-consistency calculations for each specified q, for all atoms and for
         ! all three cartesian directions.
         ! TODO: The effort here should be greatly reducible by symmetry considerations.
         write(*,*) fi%juPhon%startq/=0, fi%juPhon%stopq, size(q_list)
          
         DO iQ = fi%juPhon%startq, MERGE(fi%juPhon%stopq,SIZE(q_list),fi%juPhon%stopq/=0)
            CALL timestart("q-point")
           ! IF (.NOT.fi%juPhon%qmode==0) THEN
           !    CALL make_sym_list(fi_fullsym%sym, qpts_loc%bk(:,q_list(iQ)),sym_count,sym_list)
           !    ALLOCATE(sym_dynvec(3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype-1,sym_count))
           ! END IF
            kqpts = fi%kpts
            ! Modify this from kpts only in DFPT case.
            DO ikpt = 1, fi%kpts%nkpt
               kqpts%bk(:, ikpt) = kqpts%bk(:, ikpt) + qpts_loc%bk(:,q_list(iQ))
            END DO

            IF (l_minusq) THEN
               kqmpts = fi%kpts
               ! Modify this from kpts only in DFPT case.
               DO ikpt = 1, fi%kpts%nkpt
                  kqmpts%bk(:, ikpt) = kqmpts%bk(:, ikpt) - qpts_loc%bk(:,q_list(iQ))
               END DO
            END IF

            CALL timestart("Eii2")
            CALL CalcIIEnerg2(fi_nosym%atoms, fi_nosym%cell, qpts_loc, stars_nosym, fi_nosym%input, q_list(iQ), ngdp, recG, E2ndOrdII)
            CALL timestop("Eii2")

            write(9991,*) "Eii2 old:", E2ndOrdII
            E2ndOrdII = CMPLX(0.0,0.0)
            DO iDtype = 1, fi_nosym%atoms%ntype
               DO iDir2 = 1, 3
                  DO iDir = 1, 3
                     E2ndOrdII(3*(iDtype-1)+iDir2,3*(iDtype-1)+iDir) = e2_vm(iDtype,iDir2,iDir)
                  END DO
               END DO
            END DO
            CALL timestart("Eigenstuff at k+q")

            ! This was an additional eig_id to test a specific shift from k to k+q in the dynmat setup. We leave it in
            ! for now, as it might be tested again in the future.
            !q_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
            !                  .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)

            ! Get the eigenstuff at k+q
            CALL q_results%reset_results(fi%input)

            CALL eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv, mpdata, &
                     hybdat, 1, q_eig_id, q_results, rho, vTot, vxc, hub1data, &
                     qpts_loc%bk(:,q_list(iQ)))

            ! Fermi level and occupancies
            CALL timestart("determination of fermi energy")
            CALL fermie(q_eig_id, fmpi, kqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, q_results)
            CALL timestop("determination of fermi energy")

#ifdef CPP_MPI
            CALL MPI_BCAST(q_results%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            CALL MPI_BCAST(q_results%w_iks, SIZE(q_results%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

            CALL timestop("Eigenstuff at k+q")

            IF (l_minusq) THEN
               CALL timestart("Eigenstuff at k-q")
               !qm_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
               !                  .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)

               CALL qm_results%reset_results(fi%input)

               CALL eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv, mpdata, &
                        hybdat, 1, qm_eig_id, qm_results, rho, vTot, vxc, hub1data, &
                        -qpts_loc%bk(:,q_list(iQ)))

               ! Fermi level and occupancies
               CALL timestart("determination of fermi energy")
               CALL fermie(qm_eig_id, fmpi, kqmpts, fi%input, fi%noco, enpara%epara_min, fi%cell, qm_results)
               CALL timestop("determination of fermi energy")

#ifdef CPP_MPI
               CALL MPI_BCAST(qm_results%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
               CALL MPI_BCAST(qm_results%w_iks, SIZE(qm_results%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

               CALL timestop("Eigenstuff at k-q")
            END IF

            DO iDtype = 1, fi_nosym%atoms%ntype
               CALL timestart("Typeloop")
               DO iDir = 1, 3
                  CALL timestart("Dirloop")
                 ! IF (.NOT.fi%juPhon%qmode==0.AND.fmpi%irank==0) THEN
                 !    IF (iDtype==1.AND.iDir==2) sym_dyn_mat(iQ, 1, :) = dyn_mat(iQ, 1, :)
                 !    IF (3 *(iDtype-1)+iDir>1) THEN
                 !       CALL cheat_dynmat(fi_fullsym%atoms, fi_fullsym%sym, fi_fullsym%cell%amat, qpts_loc%bk(:,q_list(iQ)), iDtype, iDir, sym_count, sym_list(:sym_count), sym_dynvec, dyn_mat(iQ,:,:), sym_dyn_mat(iQ,:,:), l_cheated)
                 !    END IF
                 !    IF (l_cheated) WRITE(*,*) "Following row was cheated!"
                 !    IF (l_cheated) write(*,*) sym_dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
                 ! END IF
                  dfpt_tag = ''
                  WRITE(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir

                  IF (fmpi%irank==0) THEN
                     WRITE(*,*) 'Starting calculation for:'
                     WRITE(*,*) ' q         = ', qpts_loc%bk(:,q_list(iQ))
                     WRITE(*,*) ' atom      = ', iDtype
                     WRITE(*,*) ' direction = ', iDir
                  END IF

                  !IF (fmpi_nosym%irank==0) THEN
                     CALL starsq%reset_stars()
                     IF (l_minusq) CALL starsmq%reset_stars()
                     CALL denIn1%reset_dfpt()
                     CALL denIn1Im%reset_dfpt()
                     CALL vTot1%reset_dfpt()
                     CALL vTot1Im%reset_dfpt()
                     IF (l_minusq) CALL vTot1m%reset_dfpt()
                     IF (l_minusq) CALL vTot1mIm%reset_dfpt()
                     CALL vC1%reset_dfpt()
                     CALL vC1Im%reset_dfpt()
                     CALL results1%reset_results(fi_nosym%input)
                  !END IF

                  IF (fmpi%irank==0) WRITE(*,*) '-------------------------'
                  ! This is where the magic happens. The Sternheimer equation is solved
                  ! iteratively, providing the scf part of dfpt calculations.
                  IF (l_minusq) THEN
                     CALL timestart("Sternheimer with -q")
                     CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, q_results, enpara_nosym, hybdat_nosym, &
                                          rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                          dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                          denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, MERGE(sigma_ext,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3), &
                                          MERGE(sigma_coul,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3),&
                                          starsmq, qm_results, dfpt_eigm_id, dfpt_eigm_id2, qm_eig_id, results1m, vTot1m, vTot1mIm)
                     CALL timestop("Sternheimer with -q")
                  ELSE
                     CALL timestart("Sternheimer")
                     IF (fi%juPhon%l_efield) THEN
                        CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, q_results, enpara_nosym, hybdat_nosym, &
                                             rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                             dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                             denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, MERGE(sigma_ext,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3), &
                                             MERGE(sigma_coul,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3))
                     ELSE
                        CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, q_results, enpara_nosym, hybdat_nosym, &
                                             rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                             dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                             denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, MERGE(sigma_ext,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3), &
                                             MERGE(sigma_coul,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir==3))
                     END IF
                     CALL timestop("Sternheimer")
                  END IF

                  IF (fmpi%irank==0) WRITE(*,*) '-------------------------'
                  CALL timestart("Dynmat")
                  ! Once the first order quantities are converged, we can construct all
                  ! additional necessary quantities and from that the dynamical matrix.
                  IF (fi%juPhon%l_efield) THEN
                     CALL dfpt_dielecten_row_HF(fi_nosym,stars_nosym,starsq,sphhar_nosym,fmpi_nosym,denIn1,denIn1Im,results_nosym, results1,3 *(iDtype-1)+iDir,diel_tensor(3 *(iDtype-1)+iDir,:))
                     CALL dfpt_dielecten_occ1(fi,fmpi,results1,we1_data,eig1_data,diel_tensor_occ1,3 *(iDtype-1)+iDir)

                     diel_tensor(:,:) = diel_tensor(:,:) + diel_tensor_occ1(:,:)
                     print*, 'diel_tensor(:,:)',diel_tensor(:,:)
                  ELSE
                     IF(.TRUE.) THEN
                        CALL dfpt_dynmat_row(fi_nosym, stars_nosym, starsq, sphhar_nosym, xcpot_nosym, nococonv_nosym, hybdat_nosym, fmpi_nosym, qpts_loc, q_list(iQ), iDtype, iDir, &
                                          eig_id, dfpt_eig_id, dfpt_eig_id2, enpara_nosym, results_nosym, results1, l_real,&
                                          rho_nosym, vTot_nosym, grRho3, grVext3, grVC3, &
                                          denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), E2ndOrdII, sigma_ext, sigma_gext)
                     ELSE
                        CALL dfpt_dynmat_row(fi_nosym, stars_nosym, starsq, sphhar_nosym, xcpot_nosym, nococonv_nosym, hybdat_nosym, fmpi_nosym, qpts_loc, q_list(iQ), iDtype, iDir, &
                                          eig_id, dfpt_eig_id, dfpt_eig_id2, enpara_nosym, results_nosym, results1, l_real,&
                                          rho_nosym, vTot_nosym, grRho3, grVext3, grVC3, &
                                          denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), E2ndOrdII, sigma_ext, sigma_gext, q_eig_id)
                     END IF
                  END IF
                  print*,"diel_tensor(3 *(iDtype-1)+iDir,:)", diel_tensor(3 *(iDtype-1)+iDir,:)
                  !stop
                  CALL timestop("Dynmat")
                  dyn_mat(iQ,3 *(iDtype-1)+iDir,:) = dyn_mat(iQ,3 *(iDtype-1)+iDir,:) + conjg(E2ndOrdII(3 *(iDtype-1)+iDir,:))
                  !IF (.NOT.fi%juPhon%qmode==0) THEN
                  !   CALL make_sym_dynvec(fi_fullsym%atoms, fi_fullsym%sym, fi_fullsym%cell%amat, qpts_loc%bk(:,q_list(iQ)), iDtype, iDir, sym_count, sym_list(:sym_count), dyn_mat(iQ,3 *(iDtype-1)+iDir,:), sym_dynvec)
                  !END IF

                  IF (fmpi%irank==0) write(*,*) "dynmat row for ", dfpt_tag
                  IF (fmpi%irank==0) write(*,*) dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
                  !IF (fmpi%irank==0.AND.l_cheated) write(*,*) "The cheat:"
                  !IF (fmpi%irank==0.AND.l_cheated) write(*,*) sym_dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
                  l_cheated = .FALSE.
                  IF (fmpi%irank==0) WRITE(9339,*) dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
                  CALL timestop("Dirloop")
               END DO
               CALL timestop("Typeloop")

#if defined(CPP_MPI)
               CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif      
            END DO
            diel_tensor(:,:) = dielten_iden(:,:) - (fpi_const/fi%cell%omtil)*diel_tensor(:,:)
            call save_npy("diel_tensor.npy",diel_tensor(:,:))
            STOP
            !IF (fi%juPhon%l_efield) THEN
               !CALL dfpt_dielecten_occ1(3 *(iDtype-1)+iDir)
            !END IF

            IF (fmpi%irank==0) THEN
               WRITE(*,*) '-------------------------'
               CALL timestart("Dynmat diagonalization")
               CALL DiagonalizeDynMat(fi_nosym%atoms, qpts_loc%bk(:,q_list(iQ)), fi%juPhon%calcEigenVec, dyn_mat(iQ,:,:), eigenVals, eigenVecs, q_list(iQ),.TRUE.,"raw",fi_nosym%juphon%l_sumrule)
               CALL timestop("Dynmat diagonalization")

               CALL timestart("Frequency calculation")
               CALL CalculateFrequencies(fi_nosym%atoms, q_list(iQ), eigenVals, eigenFreqs,"raw")
               CALL timestop("Frequency calculation")
               write(9991,*) "Eii2 new:", E2ndOrdII
               DEALLOCATE(eigenVals, eigenVecs, eigenFreqs, E2ndOrdII)
            END IF
            !CALL close_eig(q_eig_id)
            !IF (l_minusq) CALL close_eig(qm_eig_id)
            !IF (.NOT.fi%juPhon%qmode==0) THEN
            !   DEALLOCATE(sym_dynvec)
            !END IF
            CALL timestop("q-point")

         END DO
      END IF

      ! If the Dynmats-Files were already created, we can read them in and do postprocessing.
      ! a) Transform the q-Mesh onto real space.
      ! b) Transform it back onto a dense q-path.
      ! c) Transform it back to a denser grid
      ! d) Perform a DOS calculation for the denser grid.
      IF (fmpi%irank==0) THEN ! Band/Dos stuff
         ! 0) Read
         DO iQ = 1, fi_fullsym%kpts%nkpt ! Loop over dynmat files to read
            IF (iQ<=9) THEN
               OPEN( 3001, file="dynMatq=000"//int2str(iQ), status="old")
            ELSE
               OPEN( 3001, file="dynMatq=00"//int2str(iQ), status="old")
            END IF
            DO iread = 1, 3 + 3*fi%atoms%nat ! Loop over dynmat rows
               IF (iread<4) THEN
                  READ( 3001,*) trash
                  write(*,*) iread, trash
               ELSE
                  READ( 3001,*) numbers(iread-3,:)
                  write(*,*) iread, numbers(iread-3,:)
                  dyn_mat(iQ,iread-3,:) = CMPLX(numbers(iread-3,::2),numbers(iread-3,2::2))
               END IF
            END DO ! iread
            CLOSE(3001)
         END DO ! iQ

         ! a) Real space transformation
         ALLOCATE(dyn_mat_r(fi_fullsym%kpts%nkptf,3*fi%atoms%nat,3*fi%atoms%nat))
         CALL ft_dyn(fi_fullsym%atoms, fi_fullsym%kpts, fi_fullsym%sym, fi_fullsym%cell%amat, dyn_mat, dyn_mat_r, dyn_mat_q_full)
         
         ! b/c) reciprocal space transformation for bands/dense grid
         IF (l_dfpt_band.OR.l_dfpt_full) THEN
            IF (l_dfpt_band) THEN
               dynfiletag = "band"
            ELSE IF (l_dfpt_full) THEN
               dynfiletag = "full"
            ELSE
               dynfiletag = "intp"
            END IF
            IF (l_dfpt_dos) ALLOCATE(eigenValsFull(3*fi%atoms%nat,fi_nosym%kpts%nkpt,fi_nosym%input%jspins))
            DO iQ = 1, fi_nosym%kpts%nkpt
               CALL ift_dyn(fi_fullsym%atoms,fi_fullsym%kpts,fi_fullsym%sym,fi_fullsym%cell%amat,fi_nosym%kpts%bk(:,iQ),dyn_mat_r,dyn_mat_pathq)
               WRITE(*,*) '-------------------------'
               CALL timestart("Dynmat diagonalization")
               CALL DiagonalizeDynMat(fi_nosym%atoms, fi_nosym%kpts%bk(:,iQ), fi%juPhon%calcEigenVec, dyn_mat_pathq, eigenVals, eigenVecs, iQ,.FALSE.,TRIM(dynfiletag),fi_nosym%juphon%l_sumrule)
               CALL timestop("Dynmat diagonalization")

               CALL timestart("Frequency calculation")
               CALL CalculateFrequencies(fi_nosym%atoms, iQ, eigenVals, eigenFreqs,TRIM(dynfiletag))
               CALL timestop("Frequency calculation")

               IF (l_dfpt_dos) eigenValsFull(:,iQ,1) = eigenFreqs(:)

               DEALLOCATE(eigenVals, eigenVecs, eigenFreqs, dyn_mat_pathq)
            END DO ! iQ
         END IF ! bands/interpolation
         IF (l_dfpt_dos) THEN
            fi_nosym%banddos%dos = .TRUE.
            CALL dos%init(fi_nosym%input,fi_nosym%atoms,fi_nosym%kpts,fi_nosym%banddos,eigenValsFull)
            allocate(eigdos(1))
            eigdos(1)%p=>dos
            CALL make_dos(fi_nosym%kpts,fi_nosym%atoms,fi_nosym%vacuum,fi_nosym%input,fi_nosym%banddos,fi_nosym%sliceplot,fi_nosym%noco,fi_nosym%sym,fi_nosym%cell,results,eigdos,fi%juPhon )
         END IF ! dos
      END IF ! Band/Dos stuff

      CALL close_eig(q_eig_id)
      CALL close_eig(dfpt_eig_id)
      CALL close_eig(dfpt_eig_id2)
      IF (l_minusq) THEN
         CALL close_eig(qm_eig_id)
         CALL close_eig(dfpt_eigm_id)
         CALL close_eig(dfpt_eigm_id2)
      END IF

      DEALLOCATE(recG)

      WRITE (oUnit,*) '------------------------------------------------------'

      CALL juDFT_end("Phonon calculation finished.",fmpi%irank)

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

         CALL save_npy("sym_on_rhopw.npy",rho%pw)
         CALL save_npy("sym_off_rhopw.npy",rho_nosym%pw)
         CALL save_npy("sym_on_rhomt.npy",rho%mt)
         CALL save_npy("sym_off_rhomt.npy",rho_nosym%mt)
         CALL save_npy("sym_on_vpw.npy",vTot%pw)
         CALL save_npy("sym_off_vpw.npy",vTot_nosym%pw)
         CALL save_npy("sym_on_vmt.npy",vTot%mt)
         CALL save_npy("sym_off_vmt.npy",vTot_nosym%mt)
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
      
      IF (fi_nosym%input%film)CALL save_npy("rhovac.npy",rho_nosym%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("rhogr1vac.npy",grRho3(1)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("rhogr2vac.npy",grRho3(2)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("rhogr3vac.npy",grRho3(3)%vac(:,:,:,1))
      CALL save_npy("rhogr3pw.npy",grRho3(3)%pw(:,1))
      IF (fi_nosym%input%film)CALL save_npy("vcgr1vac.npy",grVC3(1)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vcgr2vac.npy",grVC3(2)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vcgr3vac.npy",grVC3(3)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vextgr1vac.npy",grVext3(1)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vextgr2vac.npy",grVext3(2)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vextgr3vac.npy",grVext3(3)%vac(:,:,:,1))
      CALL save_npy("vextgr1pw.npy",grVext3(1)%pw(:,1))
      CALL save_npy("vextgr2pw.npy",grVext3(2)%pw(:,1))
      CALL save_npy("vextgr3pw.npy",grVext3(3)%pw(:,1))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr1vac.npy",grVtot3(1)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr2vac.npy",grVtot3(2)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr3vac.npy",grVtot3(3)%vac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr1vacnum.npy",grVtotvac(:,:,:,1))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr2vacnum.npy",grVtotvac(:,:,:,2))
      IF (fi_nosym%input%film)CALL save_npy("vtotgr3vacnum.npy",grVtotvac(:,:,:,3))
      CALL save_npy("vtotgr1pw.npy",grVtot3(1)%pw(:,1))
      CALL save_npy("vtotgr2pw.npy",grVtot3(2)%pw(:,1))
      CALL save_npy("vtotgr3pw.npy",grVtot3(3)%pw(:,1))
      CALL save_npy("vtotgr1pwnum.npy",grVtotpw(:,1))
      CALL save_npy("vtotgr2pwnum.npy",grVtotpw(:,2))
      CALL save_npy("vtotgr3pwnum.npy",grVtotpw(:,3))
   END SUBROUTINE

END MODULE m_dfpt
