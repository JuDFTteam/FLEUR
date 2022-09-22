!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt
   USE m_juDFT
   USE m_constants
   USE m_types
   USE m_dfpt_check
   USE m_dfpt_test
   !USE m_dfpt_init
   USE m_dfpt_sternheimer
   USE m_dfpt_dynmat
   !USE m_jpSternheimer,     only : solveSternheimerSCC
   USE m_jp2ndOrdQuant,     only : CalcIIEnerg2
   !USE m_jpSetupDynMat,     only : SetupDynamicMatrix
   USE m_jpProcessDynMat!,   only : DiagonalizeDynMat, CalculateFrequencies
   USE m_juDFT_stop, only : juDFT_error
   USE m_vgen_coulomb
   USE m_dfpt_vgen
   USE m_convol
   USE m_fleur_init
   USE m_npy
   USE m_desymmetrizer
   USE m_outcdn
   USE m_plot

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho, vTot, vxc, exc, vCoul, eig_id, nvfull, oldmode, xcpot, hybdat, mpdata, forcetheo)

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
      TYPE(t_potden),   INTENT(IN)    :: vTot, vxc, exc, vCoul
      INTEGER,          INTENT(IN)    :: eig_id
      INTEGER,          INTENT(IN)    :: nvfull(:, :)
      LOGICAL,          INTENT(IN)    :: oldmode

      TYPE(t_usdus)                 :: usdus
      TYPE(t_potden)                :: vTotclean, rhoclean, grRho, grvextdummy, imagrhodummy, rho_nosym, vTot_nosym
      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im ! q-quantities
      TYPE(t_jpPotden)              :: rho0, grRho0, vTot0, grVTot0
      TYPE(t_tlmplm)                :: tdHS0
      TYPE(t_results)               :: results1
      TYPE(t_kpts)                  :: qpts_loc

      ! Desymmetrized type variables:
      TYPE(t_mpi)        :: fmpi_nosym
      TYPE(t_fleurinput) :: fi_nosym
      TYPE(t_sphhar)     :: sphhar_nosym
      TYPE(t_stars)      :: stars_nosym, starsq
      TYPE(t_nococonv)   :: nococonv_nosym
      TYPE(t_enpara)     :: enpara_nosym
      TYPE(t_results)    :: results_nosym
      TYPE(t_wann)       :: wann_nosym
      TYPE(t_hybdat)     :: hybdat_nosym
      TYPE(t_mpdata)     :: mpdata_nosym

      CLASS(t_xcpot),     ALLOCATABLE :: xcpot_nosym
      CLASS(t_forcetheo), ALLOCATABLE :: forcetheo_nosym

      !COMPLEX, ALLOCATABLE :: exc_pw_nosym(:,:), vCoul_pw_nosym(:,:)

      !integer                       :: logUnit = 100
      !integer                       :: ngpqdp

      !COMPLEX, ALLOCATABLE          :: loosetd(:, :, :, :)
      !REAL,             ALLOCATABLE :: El(:, :, :, :)
      INTEGER,          ALLOCATABLE :: recG(:, :)!, GbasVec(:, :), ilst(:, :, :)
      INTEGER                       :: ngdp2km
      !INTEGER,          ALLOCATABLE :: gdp2Ind(:, :, :)
      !INTEGER                       :: gdp2iLim(2, 3)
      !INTEGER,          ALLOCATABLE :: nRadFun(:, :), iloTable(:, :, :), ilo2p(:, :)
      !REAL,             ALLOCATABLE :: uuilon(:, :)
      !REAL,             ALLOCATABLE :: duilon(:, :)
      !REAL,             ALLOCATABLE :: ulouilopn(:, :, :)
      !INTEGER,          ALLOCATABLE :: kveclo(:,:)
      !REAL,             ALLOCATABLE :: rbas1(:, :, :, :, :)
      !REAL,             ALLOCATABLE :: rbas2(:, :, :, :, :)
      !REAL,             ALLOCATABLE :: gridf(:, :)
      !COMPLEX,          ALLOCATABLE :: z0(:, :, :, :)
      !complex,           allocatable :: grVxcIRKern(:)
      !real,              allocatable :: dKernMTGPts(:, :, :)
      !real,              allocatable :: gausWts(:)
      !complex,           allocatable :: ylm(:, :)
      !complex,           allocatable :: qpwcG(:, :)
      !complex,           allocatable :: rho1MTCoreDispAt(:, :, :, :)
      !complex,           allocatable :: grVeff0MT_init(:, :, :, :)
      !complex,           allocatable :: grVeff0MT_main(:, :, :, :)
      !complex,           allocatable :: grVext0IR_DM(:, :)
      !complex,           allocatable :: grVext0MT_DM(:, :, :, :)
      !complex,           allocatable :: grVCoul0IR_DM_SF(:, :)
      !complex,           allocatable :: grVCoul0MT_DM_SF(:, :, :, :)
      !complex,           allocatable :: grVeff0IR_DM(:, :)
      !complex,           allocatable :: grVeff0MT_DM(:, :, :, :)
      !complex,           allocatable :: dynMat(:, :)
      complex,           allocatable :: E2ndOrdII(:, :)
      complex,           allocatable :: eigenFreqs(:)
      real,              allocatable :: eigenVals(:)
      complex,           allocatable :: eigenVecs(:, :)
      !integer,           allocatable :: gpqdp(:, :)
      !INTEGER,           ALLOCATABLE :: nocc(:, :)
      !complex,            allocatable :: rho1IR(:, :, :)
      !complex,            allocatable :: rho1MT(:, :, :, :, :)
      !complex,            allocatable :: rho1MTDelta(:, :, :, :, :)
      !complex,            allocatable :: rho1MTz0(:, :, :, :)
      !complex,            allocatable :: vCoul1IRtempNoVol(:, :)
      !complex,            allocatable :: vCoul1MTtempNoVol(:, :, :, :)
      !complex,            allocatable :: vEff1IR(:, :, :)
      !complex,            allocatable :: vEff1MT(:, :, :, :, :)
      !complex,            allocatable :: vEff1MTnoVol(:, :, :, :, :)
      !complex,            allocatable :: vExt1IR_final(:, :, :)
      !complex,            allocatable :: vExt1MT(:, :, :, :, :)
      !complex,            allocatable :: vExt1MTDelta(:, :, :, :, :)
      !complex,            allocatable :: vExt1MTq0(:, :, :, :, :)
      !complex,            allocatable :: vExt1noqIR_final(:, :, :)
      !complex,            allocatable :: vHar1IR_final(:, :, :)
      !complex,            allocatable :: vHar1MTDelta(:, :, :, :, :)
      !complex,            allocatable :: vHar1MTq0(:, :, :, :, :)
      !complex,            allocatable :: vHar1MT_final(:, :, :, :, :)
      !complex,            allocatable :: vXc1MTDelta(:, :, :, :, :)
      !complex,            allocatable :: vXc1MTq0(:, :, :, :, :)
      COMPLEX, ALLOCATABLE :: grrhodummy(:, :, :, :, :)
      !integer,      allocatable :: mapKpq2K(:, :)
      !integer,      allocatable :: kpq2kPrVec(:, :, :)

      COMPLEX, ALLOCATABLE :: dyn_mat(:,:,:)

      INTEGER :: ngdp, iSpin, iType, iR, ilh, iQ, iDir, iDtype
      INTEGER :: iStar, xInd, yInd, zInd
      LOGICAL :: l_real

      CHARACTER(len=20)  :: dfpt_tag
      CHARACTER(len=100) :: inp_pref

      INTEGER, ALLOCATABLE :: q_list(:), dfpt_eig_id_list(:)

      ! Desym-tests:
      INTEGER :: ix, iy, iz, grid(3), iv_old, iflag_old, iv_new, iflag_new
      INTEGER :: iType_old, iAtom_old, iType_new, iAtom_new, inversionOp
      REAL    :: old_point(3), new_point(3), pt_old(3), pt_new(3), xdnout_old, xdnout_new, atom_shift(3)

      !INTEGER, PARAMETER :: invs_matrix(3,3)=RESHAPE([-1,0,0,0,-1,0,0,0,-1],[3,3])
      !REAL,    PARAMETER :: eps7 = 1.0e-7

      l_real = fi%sym%invs.AND.(.NOT.fi%noco%l_soc).AND.(.NOT.fi%noco%l_noco).AND.fi%atoms%n_hia==0

      !IF (fi%juPhon%l_jpCheck) THEN
    !      ! This function will be used to check the validity of juPhon's
    !      ! input. I.e. check, whether all prohibited switches are off and,
    !      ! once there is more expertise considering this topic, check whether
    !      ! the cutoffs are chosen appropriately.
    !      CALL dfpt_check(fi_nosym, xcpot_nosym)
     ! END IF

      IF (fi%sym%nop>1) THEN
         WRITE(*,*) "Desymmetrization needed. Going ahead!"
         ! Grid size for desym quality test:
         grid = 21

         inp_pref = ADJUSTL("desym_")
         fmpi_nosym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
         fmpi_nosym%mpi_comm = fmpi%mpi_comm

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
         !ALLOCATE(exc_pw_nosym(SIZE(vTot_nosym%pw,1),1))
         !ALLOCATE(vCoul_pw_nosym(SIZE(vTot_nosym%pw,1),1))

         CALL desymmetrize_pw(fi%sym, stars, stars_nosym, rho%pw, rho_nosym%pw)
         CALL desymmetrize_pw(fi%sym, stars, stars_nosym, vTot%pw, vTot_nosym%pw, vTot%pw_w, vTot_nosym%pw_w)
         CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, rho%mt, rho_nosym%mt)
         CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, vTot%mt, vTot_nosym%mt)

         !CALL desymmetrize_pw(fi%sym, stars, stars_nosym, exc%pw, exc_pw_nosym)
         !CALL desymmetrize_pw(fi%sym, stars, stars_nosym, vCoul%pw, vCoul_pw_nosym)

         CALL desymmetrize_types(fi%input, fi_nosym%input, fi%atoms, fi_nosym%atoms, fi%noco, &
                                 nococonv, nococonv_nosym, enpara, enpara_nosym, results, results_nosym)

         IF (.FALSE.) THEN
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
            !STOP
         END IF
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

!#ifndef CPP_FFTW
!        call juDFT_error('juPhon is only usable with fftw support.', calledby='dfpt')
!#endif

      CALL results1%init(fi_nosym%input, fi_nosym%atoms, fi_nosym%kpts, fi_nosym%noco)

      !WRITE (oUnit,*) '------------------------------------------------------'
      !WRITE (oUnit,*) 'This output is generated by juPhon, FLEURs DFPT addon.'
      !WRITE (oUnit,*) 'l_dfpt = ', fi%juPhon%l_dfpt
      !WRITE (oUnit,*) 'l_jpCheck = ', fi%juPhon%l_jpCheck
      !WRITE (oUnit,*) 'l_jpTest = ', fi%juPhon%l_jpTest
      !WRITE (oUnit,*) 'l_potout = ', fi%juPhon%l_potout
      !WRITE (oUnit,*) 'l_eigout = ', fi%juPhon%l_eigout
      !WRITE (oUnit,*) 'l_symTsh = ', fi%juPhon%l_symTsh
      !WRITE (oUnit,*) 'l_symTdm = ', fi%juPhon%l_symTdm
      !WRITE (oUnit,*) 'l_bfkq = ', fi%juPhon%l_bfkq
      !WRITE (oUnit,*) 'jplPlus = ', fi%juPhon%jplPlus
      !WRITE (oUnit,*) 'kgqmax = ', fi%juPhon%kgqmax
      !WRITE (oUnit,*) 'gqmax = ', fi%juPhon%gqmax
      !WRITE (oUnit,*) 'eps_pert = ', fi%juPhon%eps_pert
      !WRITE (oUnit,*) 'eDiffCut = ', fi%juPhon%eDiffCut
      !WRITE (oUnit,*) 'qpt_ph = ', fi%juPhon%qpt_ph

      ! TODO: This is a test set of qpoints for a fixed fcc system.
      !       We need to read out actual q-points at some point.
      !       And it needs to be handled properly.
      !ALLOCATE(q_list(5),dfpt_eig_id_list(5))
      !q_list = [1, 10, 19, 28, 37]! 512 k-points: \Gamma to X

      qpts_loc = qpts

      qpts_loc%bk(:,1)  = [0.0,1.0,1.0]*0.025*0.0
      qpts_loc%bk(:,2)  = [0.0,1.0,1.0]*0.025*1.0
      qpts_loc%bk(:,3)  = [0.0,1.0,1.0]*0.025*2.0
      qpts_loc%bk(:,4)  = [0.0,1.0,1.0]*0.025*3.0
      qpts_loc%bk(:,5)  = [0.0,1.0,1.0]*0.025*4.0
      qpts_loc%bk(:,6)  = [0.0,1.0,1.0]*0.025*5.0
      qpts_loc%bk(:,7)  = [0.0,1.0,1.0]*0.025*6.0
      qpts_loc%bk(:,8)  = [0.0,1.0,1.0]*0.025*7.0
      qpts_loc%bk(:,9)  = [0.0,1.0,1.0]*0.025*8.0
      qpts_loc%bk(:,10) = [0.0,1.0,1.0]*0.025*9.0
      qpts_loc%bk(:,11) = [0.0,1.0,1.0]*0.025*10.0
      qpts_loc%bk(:,12) = [0.0,1.0,1.0]*0.025*11.0
      qpts_loc%bk(:,13) = [0.0,1.0,1.0]*0.025*12.0
      qpts_loc%bk(:,14) = [0.0,1.0,1.0]*0.025*13.0
      qpts_loc%bk(:,15) = [0.0,1.0,1.0]*0.025*14.0
      qpts_loc%bk(:,16) = [0.0,1.0,1.0]*0.025*15.0
      qpts_loc%bk(:,17) = [0.0,1.0,1.0]*0.025*16.0
      qpts_loc%bk(:,18) = [0.0,1.0,1.0]*0.025*17.0
      qpts_loc%bk(:,19) = [0.0,1.0,1.0]*0.025*18.0
      qpts_loc%bk(:,20) = [0.0,1.0,1.0]*0.025*19.0
      qpts_loc%bk(:,21) = [0.0,1.0,1.0]*0.025*20.0

      !ALLOCATE(q_list(21),dfpt_eig_id_list(21))
      !ALLOCATE(q_list(1),dfpt_eig_id_list(1))
      ALLOCATE(q_list(5),dfpt_eig_id_list(5))
      !q_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21] ! \Gamma to X in 20 steps.
      !q_list = [21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]
      !q_list = [21,1,11,16,6,19,18,14,13,9,8,4,3,20,17,15,12,10,7,5,2]
      !q_list = [13]
      q_list = [1,6,11,16,21]
      ALLOCATE(grrhodummy(fi_nosym%atoms%jmtd, (fi_nosym%atoms%lmaxd+1)**2, fi_nosym%atoms%nat, SIZE(rho_nosym%mt,4), 3))

      CALL imagrhodummy%copyPotDen(rho_nosym)
      CALL imagrhodummy%resetPotDen()
      CALL grvextdummy%copyPotDen(rho_nosym)
      DO iDir = 1, 3
         CALL grRho3(iDir)%copyPotDen(rho_nosym)
         CALL grRho3(iDir)%resetPotDen()
         CALL grVext3(iDir)%copyPotDen(vTot_nosym)
         CALL grVext3(iDir)%resetPotDen()
         CALL grVtot3(iDir)%copyPotDen(vTot_nosym)
         CALL grVtot3(iDir)%resetPotDen()
         CALL grVC3(iDir)%copyPotDen(vTot_nosym)
         CALL grVC3(iDir)%resetPotDen()
         ! Generate the external potential gradient.
         CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .FALSE., imagrhodummy, grVext3(iDir), &
                         & dfptdenimag=imagrhodummy, dfptvCoulimag=grvextdummy,dfptden0=imagrhodummy,stars2=stars_nosym,iDtype=0,iDir=iDir)
      END DO

      DO iSpin = 1, SIZE(rho_nosym%mt,4)
         CALL mt_gradient_old(fi_nosym%atoms, sphhar_nosym, fi_nosym%sym, sphhar_nosym%clnu, sphhar_nosym%nmem, sphhar_nosym%mlh, rho_nosym%mt(:, :, :, iSpin), grrhodummy(:, :, :, iSpin, :))
      END DO

      DO zInd = -stars_nosym%mx3, stars_nosym%mx3
         DO yInd = -stars_nosym%mx2, stars_nosym%mx2
            DO xInd = -stars_nosym%mx1, stars_nosym%mx1
               iStar = stars_nosym%ig(xInd, yInd, zInd)
               IF (iStar.EQ.0) CYCLE
               grRho3(1)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grRho3(2)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
               grRho3(3)%pw(iStar,:) = rho_nosym%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi_nosym%cell%bmat)))
            END DO
         END DO
      END DO

      DO iDir = 1, 3
         CALL sh_to_lh(fi_nosym%sym, fi_nosym%atoms, sphhar_nosym, SIZE(rho_nosym%mt,4), 2, grrhodummy(:, :, :, :, iDir), grRho3(iDir)%mt, imagrhodummy%mt)
         CALL imagrhodummy%resetPotDen()
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi_nosym%cell, fi_nosym%sliceplot, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVtot3(iDir), .TRUE., grvextdummy, grRho3(iDir), 0, iDir, [0,0])
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi_nosym%cell, fi_nosym%sliceplot, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVC3(iDir), .FALSE., grvextdummy, grRho3(iDir), 0, iDir, [0,0])
      END DO


      ALLOCATE(dyn_mat(SIZE(q_list),3*fi_nosym%atoms%ntype,3*fi_nosym%atoms%ntype))
      DO iQ = 1, SIZE(q_list)
         CALL timestart("Eii2")
         CALL old_get_Gvecs(stars_nosym, fi_nosym%cell, fi_nosym%input, ngdp, ngdp2km, recG, .false.)
         CALL CalcIIEnerg2(fi_nosym%atoms, fi_nosym%cell, qpts_loc, stars, fi_nosym%input, q_list(iQ), ngdp, recG, E2ndOrdII)
         CALL timestop("Eii2")
         DO iDtype = 1, fi_nosym%atoms%ntype
            !DO iDir = 1, 1
            DO iDir = 1, 3
               dfpt_tag = ''
               WRITE(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir
               WRITE(*,*) '-------------------------'
               WRITE(*,*) 'Starting calculation for:'
               WRITE(*,*) ' q         = ', qpts_loc%bk(:,q_list(iQ))
               WRITE(*,*) ' atom      = ', iDtype
               WRITE(*,*) ' direction = ', iDir

               IF (fmpi_nosym%irank==0) THEN
                  CALL starsq%reset_stars()
                  CALL denIn1%reset_dfpt()
                  CALL denIn1Im%reset_dfpt()
                  CALL vTot1%reset_dfpt()
                  CALL vTot1Im%reset_dfpt()
                  CALL vC1%reset_dfpt()
                  CALL vC1Im%reset_dfpt()
               END IF
               ! TODO: Broadcast this.
               WRITE(*,*) '-------------------------'
               ! This is where the magic happens. The Sternheimer equation is solved
               ! iteratively, providing the scf part of dfpt calculations.
               CALL timestart("Sternheimer")
               CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, enpara_nosym, hybdat_nosym, mpdata_nosym, forcetheo_nosym, &
                                     rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                     dfpt_tag, eig_id, l_real, results1, dfpt_eig_id_list(iQ), &
                                     denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im)
               CALL timestop("Sternheimer")

               !CYCLE

               WRITE(*,*) '-------------------------'
               CALL timestart("Dynmat")
               ! Once the first order quantities are converged, we can construct all
               ! additional necessary quantities and from that the dynamical matrix.
               CALL dfpt_dynmat_row(fi_nosym, stars_nosym, starsq, sphhar_nosym, xcpot_nosym, nococonv_nosym, hybdat_nosym, fmpi_nosym, qpts_loc, q_list(iQ), iDtype, iDir, &
                                    eig_id, dfpt_eig_id_list(iQ), enpara_nosym, mpdata_nosym, results_nosym, results1, l_real,&
                                    rho_nosym, vTot_nosym, grRho3, grVext3, grVC3, grVtot3, &
                                    denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:))
               CALL timestop("Dynmat")
               dyn_mat(iQ,3 *(iDtype-1)+iDir,:) = dyn_mat(iQ,3 *(iDtype-1)+iDir,:) + E2ndOrdII(3 *(iDtype-1)+iDir,:)
               write(9989,*) E2ndOrdII(3 *(iDtype-1)+iDir,1)
               write(9989,*) E2ndOrdII(3 *(iDtype-1)+iDir,2)
               write(9989,*) E2ndOrdII(3 *(iDtype-1)+iDir,3)
               write(*,*) "dynmat row for ", dfpt_tag
               write(*,*) dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
            END DO
         END DO
         DEALLOCATE(recG)
         WRITE(*,*) '-------------------------'
         CALL timestart("Dynmat diagonalization")
         CALL DiagonalizeDynMat(fi%atoms, qpts_loc, fi%juPhon%calcEigenVec, dyn_mat(iQ,:,:), eigenVals, eigenVecs, q_list(iQ))
         CALL timestop("Dynmat diagonalization")

         CALL timestart("Frequency calculation")
         CALL CalculateFrequencies(fi%atoms, q_list(iQ), eigenVals, eigenFreqs)
         CALL timestop("Frequency calculation")
         DEALLOCATE(eigenVals, eigenVecs, eigenFreqs)
      END DO

        ! Construct potential without the l=0 prefactor.
        !CALL vTotclean%copyPotDen(vTot)
        !CALL rhoclean%copyPotDen(rho)

        !DO iSpin = 1, fi%input%jspins
         !   DO iType = 1, fi%atoms%ntype
         !       DO ilh = 0, sphhar%nlhd
         !           DO iR = 1, fi%atoms%jri(iType)
         !               IF (ilh.EQ.0) THEN
         !                   vTotclean%mt(iR, 0, iType, iSpin) &
         !               & = vTotclean%mt(iR, 0, iType, iSpin) * sqrt(fpi_const) / fi%atoms%rmsh(iR, iType)
         !               END IF
         !               rhoclean%mt(iR, ilh, iType, iSpin) &
         !           & = rhoclean%mt(iR, ilh, iType, iSpin) / fi%atoms%rmsh(iR, iType) / fi%atoms%rmsh(iR, iType)
         !           END DO
         !       END DO
         !   END DO
        !END DO

        ! This routine will initialize everything for juPhon that isn't already
        ! provided by the FLEUR code/must be provieded in a modified way.
        ! This includes for example the de-symmetrized MT and pw quantities and
        ! their gradients. Notably, q-dependent quantities are initialized and
        ! constructed elsewhere, within the q-loop.
        ! TODO: I ignored the actual significance of clnu_atom etc. They are not type-dependent, but actually
        ! refer to each atom respectively. So this will explode for iatom > 1. This is easily fixed.
       ! CALL timestart("juPhon DFPT initialization")
        !CALL dfpt_init(fi%juPhon, fi%sym, fi%input, fi%atoms, sphhar, stars, fi%cell, fi%noco, nococonv, fi%kpts, &
         !            & fmpi, results, enpara, rho, vTot, eig_id, nvfull, usdus, rho0, grRho0, vTot0, grVTot0, &
         !            & ngdp, El, recG, ngdp2km, gdp2Ind, gdp2iLim, GbasVec, ilst, nRadFun, iloTable, ilo2p, &
         !            & uuilon, duilon, ulouilopn, kveclo, rbas1, rbas2, gridf, z0, grVxcIRKern, dKernMTGPts, &
         !            & gausWts, ylm, qpwcG, rho1MTCoreDispAt, grVeff0MT_init, grVeff0MT_main, grVext0IR_DM, grVext0MT_DM, &
         !            & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM, tdHS0, loosetd, nocc, rhoclean, oldmode, xcpot, grRho)
        !CALL timestop("juPhon DFPT initialization")

        ! < Imagine starting a q-grid-loop here. >
        ! < For now we just select one q-point from the input. >

        !STOP

        !call createkqMapArrays( fi%kpts, qpts, 0, fi%kpts%nkpt3, [0], mapKpq2K, kpq2kPrVec )

        !CALL timestart("juPhon DFPT scf loop")
        !call solveSternheimerSCC( fmpi,  fi%atoms, fi%sym, stars, sphhar, fi%cell, enpara, usdus, fi%input, fi%kpts, qpts, results, usdus,      &
         ! & logUnit, ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
         ! & recG, mapKpq2K, results%neig(:, 1), results%eig, GbasVec, ilst, z0, nvfull, El, nradFun, iloTable, nocc, ilo2p, gdp2Ind,     &
         ! & gdp2iLim, kpq2kPrVec, qpwcG, q_list(2), tdHS0, loosetd, ylm, grRho0%pw(:, 1, 1, :), grRho0%mt(:, :, :, 1, 1, :), grVeff0MT_init, grVeff0MT_main, dKernMTGPts,       &
         ! & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, fi%juPhon%oneSternhCycle, ngpqdp, gpqdp,&
         ! & vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, &
         ! & vXc1MTq0, rho0%pw(:, :, 1, 1), rho0%mt(:, :, :, :, 1, 1), vTot0%pw(:, :, 1, 1), fi%juPhon%noPtsCon, vEff1MTnoVol, vExt1noqIR_final, rho1MTz0, vCoul1IRtempNoVol, vCoul1MTtempNoVol )
        !CALL timestop("juPhon DFPT scf loop")

        !CALL timestart("juPhon DFPT Eii2")
        !CALL CalcIIEnerg2(fi%atoms, fi%cell, qpts, stars, fi%input, q_list(2), ngdp, recG, E2ndOrdII)
        !CALL timestop("juPhon DFPT Eii2")

        !CALL timestart("juPhon DFPT dynmat setup")
        !CALL SetupDynamicMatrix( fmpi, fi%noco, nococonv,  fi%atoms, fi%input, fi%sym, fi%cell, sphhar, stars, fi%kpts, qpts, usdus, results, vTotclean, q_list(2), ngdp, ngpqdp, recG, sphhar%mlh, sphhar%nmem,&
         !   & sphhar%clnu, rho%pw, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, vTot%pw_w, vTotclean%mt(:, 0:, :, 1),&
         !   & rhoclean%mt, E2ndOrdII, El, results%eig, rbas1, rbas2, iloTable, nvfull, nocc, ilst, GbasVec, z0, kveclo, nRadFun, mapKpq2K, kpq2kPrVec,       &
         !   & gpqdp, sphhar%memd, logUnit, vxc%pw, exc%pw(:, 1), vxc%mt, exc%mt(:, 0:, :, 1), vExt1IR_final, vHar1IR_final, vHar1MT_final, grRho0%pw(:, 1, 1, :), grRho0%mt(:, :, :, 1, 1, :), &
         !   & grVext0IR_DM, grVext0MT_DM, grVeff0IR_DM, grVeff0MT_DM, dynMat, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, &
         !   & vXc1MTDelta, vXc1MTq0, vEff1MTnoVol, vExt1noqIR_final, rho1MTz0, &
         !   & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol)
        !CALL timestop("juPhon DFPT dynmat setup")

        !CALL timestart("juPhon DFPT dynmat diagonalization")
        !CALL DiagonalizeDynMat(fi%atoms, qpts, fi%juPhon%calcEigenVec, dynMat, eigenVals, eigenVecs, q_list(2))
        !CALL timestop("juPhon DFPT dynmat diagonalization")

        !CALL timestart("juPhon DFPT frequency calculation")
        !CALL CalculateFrequencies(fi%atoms, q_list(2), eigenVals, eigenFreqs)
        !CALL timestop("juPhon DFPT frequency calculation")

        ! < Imagine ending a q-grid-loop here. >

        WRITE (oUnit,*) '------------------------------------------------------'

        CALL juDFT_end("Phonon calculation finished.")

    END SUBROUTINE dfpt

    subroutine createkqMapArrays( kpts, qpts, nrAddQs, kSetDim, addQnkptis, mapKpq2K, kpq2kPrVec )

      use m_types
      use m_juDFT_stop, only : juDFT_error

      implicit none

      ! Type parameters
      type(t_kpts),              intent(in)  :: kpts
      type(t_kpts),              intent(in)  :: qpts

      ! Array parameter
      integer,                   intent(in)  :: addQnkptis(:)
      integer,                   intent(in)  :: kSetDim(:)
      integer,      allocatable, intent(out) :: mapKpq2K(:, :)
      integer,      allocatable, intent(out) :: kpq2kPrVec(:, :, :)

      ! Scalar parameter
      integer,                   intent(in)  :: nrAddQs

      ! Array variable
      integer,      allocatable              :: mapK2Ind(:, :, :) ! takes components of kpts%bk * lcm and gives its index
      integer,      allocatable              :: mapK2mK(:)
      integer                                :: kpqNomin(3)       ! helps to find k + q mapped back to Brillouin zone
      real                                   :: kpqTemp(3)        ! helps to find k + q mapped back to Brillouin zone
      character(len=1024)                    :: errorMessage      ! stores error message for error output

      ! Scalar variable
      integer                                :: ikpt              ! loop variable
      integer                                :: maxKcomp(3)       ! stores maximal k-point component * lcm
      integer                                :: idir              ! loop variable
      integer                                :: iqpt              ! loop variable
      integer                                :: nkptShift         ! stores shift of shifted k-point set
      integer                                :: ikptSh            ! loop variable
      logical                                :: matchFound        ! is true if index for k + q is found
      real :: lcm

      lcm = real( kgv(kSetDim, 3) )

      ! Determine maximal value of k-vector per direction in internal representation for allocation of mapKpq2K array and allocate it
      maxKcomp = 0
      do idir = 1, 3
        maxKcomp(idir) = maxval( kpts%bk(idir, :kpts%nkpt) * lcm )
      end do
      allocate( mapK2Ind(0:maxKcomp(1), 0:maxKcomp(2), 0:maxKcomp(3)) )
      mapK2Ind = 0
      allocate( mapKpq2K(kpts%nkpt, qpts%nkptf + nrAddQs) )

      ! Fill up array which stores the index of a given k-point
      do ikpt = 1, kpts%nkpt
        mapK2Ind(nint(kpts%bk(1, ikpt) * lcm), nint(kpts%bk(2, ikpt) * lcm), nint(kpts%bk(3, ikpt) * lcm)) = ikpt
      end do

      ! Determine k-point on which k + q can be folded back and determine the respective reciprocal lattice vector.
      ! The absolute value of every coordinate of the reciprocal lattice vector can be 1 maximally.
      allocate( kpq2kPrVec(3, kpts%nkpt, qpts%nkpt) )
      kpq2kPrVec = 0
      do iqpt = 1, qpts%nkpt
        do ikpt = 1, kpts%nkpt
          kpqNomin = 0
          do idir = 1, 3
            kpqNomin(idir) = nint( mod( kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt), 1. ) * lcm )
            !kpqNomin(idir) = nint(lcm*kpts%bk(idir, ikpt) + lcm*qpts%bk(idir, iqpt))
            ! Is in 1st Brillouin zone
            if (abs(real(kpqNomin(idir)) / real(lcm) - (kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt))) < 1e-5) then
            !if (kpqNomin(idir).lt.int(lcm)) then
              kpq2kPrVec(idir, ikpt, iqpt) = 0
            ! Has to be backfolded
            else
              kpq2kPrVec(idir, ikpt, iqpt) = -1
            end if
          end do
          if(.false.) then
            write(1005, '(i5, i5, 3(i5))') iqpt, ikpt, kpq2kPrVec(:, ikpt, iqpt)
          end if
          mapKpq2K(ikpt, iqpt) = mapK2Ind( kpqNomin(1), kpqNomin(2), kpqNomin(3) )
          !mapKpq2K(ikpt, iqpt) = mapK2Ind( mod(kpqNomin(1),int(lcm)), mod(kpqNomin(2),int(lcm)), mod(kpqNomin(3),int(lcm)) )
        end do
      end do

      ! For this found k-vector equals to k + q, determine the index and fill up array which connects the kpts indices of the k and q
      ! qpoint with the index of the k-vector equals to k + q.
      nkptShift = 0
      matchFound = .false.
      do iqpt = qpts%nkpt + 1, qpts%nkpt + nrAddQs
        do ikpt = 1, kpts%nkpt
          kpqTemp(:) = modulo1r( kpts%bk(:, ikpt) + qpts%bk(:, iqpt) )
          do ikptSh = 1, addQnkptis(iqpt - qpts%nkpt)
            if ( norm2( kpts%bk(:, kpts%nkpt + ikptSh + nkptShift) - kpqTemp(:) ) < 1e-7 ) then
              mapKpq2K(ikpt, iqpt) = kpts%nkpt + nkptShift + ikptSh
              matchFound = .true.
              exit
            end if
          end do
          if ( .not.matchFound ) then
            write (errorMessage, '(a,1x,i3)') 'no match for k+q-point', kpqTemp
            call juDFT_error( errorMessage, calledby='createkqMapArrays', hint='Please check k-points, q-points and mapKpq2K array!' )
          else
            matchFound = .false.
          end if
        end do
        nkptShift = nkptShift + addQnkptis(iqpt - qpts%nkpt)
      end do

      if (.false.) then
        ! Finds out which k' results when sending k to -k and then backfolding into 1st Brillouin zone.
        do ikpt = 1, kpts%nkpt
          do idir = 1, 3
            if ( kpts%bk(idir, ikpt)==0 ) then
              kpqTemp(idir) = kpts%bk(idir, ikpt)
            else
              kpqTemp(idir) = -kpts%bk(idir, ikpt) + 1
            end if
          end do ! idir
          mapK2mK(ikpt) = mapK2Ind(int(kpqTemp(1) * lcm), int(kpqTemp(2) * lcm), int(kpqTemp(3) * lcm))
        end do ! ikpt
      end if

    end subroutine createkqMapArrays

    function modulo1r(kpoint)
      implicit none
      real(8), intent(in) :: kpoint(3)
      real(8)             :: modulo1r(3)
      integer             :: i

      modulo1r = modulo (kpoint , 1d0)

      do i = 1,3
        if(abs(1-abs(modulo1r(i))).lt.1d-13) modulo1r(i) = 0d0
      enddo
    end function modulo1r

    function kgv(iarr,n)
      implicit none
      integer              :: kgv
      integer, intent(in)  :: n,iarr(n)
      logical              :: lprim(2:maxval(iarr))
      integer, allocatable :: prim(:),expo(:)
      integer              :: nprim,marr
      integer              :: i,j,ia,k
      ! Determine prime numbers
      marr  = maxval(iarr)
      lprim = .true.
      do i = 2,marr
        j = 2
        do while (i*j.le.marr)
          lprim(i*j) = .false.
          j          = j + 1
        enddo
      enddo
      nprim = count(lprim)
      allocate ( prim(nprim),expo(nprim) )
      j = 0
      do i = 2,marr
        if(lprim(i)) then
          j       = j + 1
          prim(j) = i
        endif
      enddo
      ! Determine least common multiple
      expo = 0
      do i = 1,n
        ia = iarr(i)
        if(ia.eq.0) cycle
        do j = 1,nprim
          k = 0
          do while(ia/prim(j)*prim(j).eq.ia)
            k  = k + 1
            ia = ia / prim(j)
          enddo
          expo(j) = max(expo(j),k)
        enddo
      enddo
      kgv = 1
      do j = 1,nprim
        kgv = kgv * prim(j)**expo(j)
      enddo
      deallocate ( prim,expo )
    end function kgv

    subroutine old_get_Gvecs(starsT, cellT, inputT, ngdp, ngdp2km, gdp, testMode )

      use m_juDFT

      type(t_stars),          intent(in)   :: starsT
      type(t_cell),           intent(in)   :: cellT
      type(t_input),          intent(in)   :: inputT

      integer,                intent(out)  :: ngdp
      integer,                intent(out)  :: ngdp2km
      logical,                intent(in)   :: testMode

      integer,  allocatable,  intent(out)  :: gdp(:, :)

      integer                              :: Gx, Gy, Gz, iG
      integer                              :: ngrest

      integer                              :: gdptemp2kmax(3, (2 * starsT%mx1 + 1) * (2 * starsT%mx2 + 1) * (2 * starsT%mx3 +  1))
      integer                              :: gdptemprest(3, (2 * starsT%mx1 + 1) * (2 * starsT%mx2 + 1) * (2 * starsT%mx3 +  1))
      integer                              :: Gint(3)
      real                                 :: Gext(3)

      ngdp = 0
      gdptemp2kmax = 0
      gdptemprest = 0
      ngdp2km = 0
      ngrest = 0

      do Gx = -starsT%mx1, starsT%mx1
        do Gy = -starsT%mx2, starsT%mx2
          do Gz = -starsT%mx3, starsT%mx3
            Gint = [Gx, Gy, Gz]
            Gext =  matmul(cellT%bmat, Gint)
            if (norm2(Gext) <= inputT%gmax) then
              ngdp = ngdp + 1
              ! Sort G-vectors
              if ( norm2(Gext) <= 2 * inputT%rkmax ) then
                ngdp2km = ngdp2km + 1
                gdptemp2kmax(:, ngdp2km) = Gint(:)
              else
                ngrest = ngrest + 1
                gdptemprest(:, ngrest) = Gint(:)
              end if
            endif
          enddo !Gz
        enddo !Gy
      enddo !Gx
      allocate(gdp(3, ngdp))
      gdp(:, :ngdp2km) = gdptemp2kmax(:, :ngdp2km)
      gdp(:, ngdp2km + 1 : ngdp) = gdptemprest(:, :ngrest)

    end subroutine old_get_Gvecs
END MODULE m_dfpt
