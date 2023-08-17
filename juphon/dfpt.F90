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
   !USE m_dfpt_test
   USE m_dfpt_sternheimer
   USE m_dfpt_dynmat
   USE m_jp2ndOrdQuant,     only : CalcIIEnerg2, genPertPotDensGvecs
   USE m_jpProcessDynMat
   USE m_juDFT_stop, only : juDFT_error
   USE m_vgen_coulomb
   USE m_dfpt_vgen
   USE m_fleur_init
   USE m_npy
   USE m_desymmetrizer
   USE m_outcdn
   USE m_plot
   USE m_eigen
   USE m_fermie

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho, vTot, vxc, exc, vCoul, eig_id, xcpot, hybdat, mpdata, forcetheo)

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

      TYPE(t_usdus)                 :: usdus
      TYPE(t_hub1data) :: hub1data
      TYPE(t_potden)                :: grRho, grvextdummy, imagrhodummy, rho_nosym, vTot_nosym
      TYPE(t_potden)                :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3)
      TYPE(t_potden)                :: denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, vTot1m, vTot1mIm ! q-quantities
      TYPE(t_results)               :: q_results, results1, qm_results, results1m
      TYPE(t_kpts)                  :: qpts_loc
      TYPE(t_kpts)              :: kqpts, kqmpts ! basically kpts, but with q added onto each one.

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

      INTEGER,          ALLOCATABLE :: recG(:, :)
      INTEGER                       :: ngdp2km
      complex,           allocatable :: E2ndOrdII(:, :)
      complex,           allocatable :: eigenFreqs(:)
      real,              allocatable :: eigenVals(:)
      complex,           allocatable :: eigenVecs(:, :)

      COMPLEX, ALLOCATABLE :: grrhodummy(:, :, :, :, :)

      COMPLEX, ALLOCATABLE :: dyn_mat(:,:,:)

      INTEGER :: ngdp, iSpin, iType, iQ, iDir, iDtype, nspins
      INTEGER :: iStar, xInd, yInd, zInd, q_eig_id, ikpt, ierr, qm_eig_id, iArray
      INTEGER :: dfpt_eig_id, dfpt_eig_id2, dfpt_eigm_id, dfpt_eigm_id2
      LOGICAL :: l_real, l_minusq

      CHARACTER(len=20)  :: dfpt_tag
      CHARACTER(len=100) :: inp_pref

      INTEGER, ALLOCATABLE :: q_list(:)!, dfpt_eig_id_list(:), dfpt_eigm_id_list(:), dfpt_eig_id_list2(:), dfpt_eigm_id_list2(:)

      ! Desym-tests:
      INTEGER :: ix, iy, iz, grid(3), iv_old, iflag_old, iv_new, iflag_new
      INTEGER :: iType_old, iAtom_old, iType_new, iAtom_new, inversionOp
      REAL    :: old_point(3), new_point(3), pt_old(3), pt_new(3), xdnout_old, xdnout_new, atom_shift(3)

      l_real = fi%sym%invs.AND.(.NOT.fi%noco%l_soc).AND.(.NOT.fi%noco%l_noco).AND.fi%atoms%n_hia==0

      l_minusq = .FALSE.

      nspins = MERGE(2, 1, fi%noco%l_noco)

      IF (fi%juPhon%l_jpCheck) THEN
          ! This function will be used to check the validity of juPhon's
          ! input. I.e. check, whether all prohibited switches are off and,
          ! once there is more expertise considering this topic, check whether
          ! the cutoffs are chosen appropriately.
          CALL dfpt_check(fi_nosym, xcpot_nosym)
      END IF

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

         CALL desymmetrize_pw(fi%sym, stars, stars_nosym, rho%pw, rho_nosym%pw)
         CALL desymmetrize_pw(fi%sym, stars, stars_nosym, vTot%pw, vTot_nosym%pw, vTot%pw_w, vTot_nosym%pw_w)
         CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, rho%mt, rho_nosym%mt)
         CALL desymmetrize_mt(fi%sym, fi_nosym%sym, fi%cell, fi%atoms, fi_nosym%atoms, sphhar, sphhar_nosym, vTot%mt, vTot_nosym%mt)

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

      CALL q_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
      CALL results1%init(fi_nosym%input, fi_nosym%atoms, fi_nosym%kpts, fi_nosym%noco)
      IF (l_minusq) THEN
         CALL qm_results%init(fi%input, fi%atoms, fi%kpts, fi%noco)
         CALL results1m%init(fi_nosym%input, fi_nosym%atoms, fi_nosym%kpts, fi_nosym%noco)
      END IF

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

      qpts_loc = qpts

      !qpts_loc%bk(:,1)  = [0.0,1.0,1.0]*0.00625*0.0
      !qpts_loc%bk(:,2)  = [0.0,1.0,1.0]*0.00625*1.0
      !qpts_loc%bk(:,3)  = [0.0,1.0,1.0]*0.00625*2.0
      !qpts_loc%bk(:,4)  = [0.0,1.0,1.0]*0.00625*3.0
      !qpts_loc%bk(:,5)  = [0.0,1.0,1.0]*0.00625*4.0
      !qpts_loc%bk(:,6)  = [0.0,1.0,1.0]*0.00625*5.0
      !qpts_loc%bk(:,7)  = [0.0,1.0,1.0]*0.00625*6.0
      !qpts_loc%bk(:,8)  = [0.0,1.0,1.0]*0.00625*7.0
      !qpts_loc%bk(:,9)  = [0.0,1.0,1.0]*0.00625*8.0
      !qpts_loc%bk(:,10) = [0.0,1.0,1.0]*0.00625*9.0
      !qpts_loc%bk(:,11) = [0.0,1.0,1.0]*0.00625*10.0
      !qpts_loc%bk(:,12) = [0.0,1.0,1.0]*0.00625*11.0
      !qpts_loc%bk(:,13) = [0.0,1.0,1.0]*0.00625*12.0
      !qpts_loc%bk(:,14) = [0.0,1.0,1.0]*0.00625*13.0
      !qpts_loc%bk(:,15) = [0.0,1.0,1.0]*0.00625*14.0
      !qpts_loc%bk(:,16) = [0.0,1.0,1.0]*0.00625*15.0
      !qpts_loc%bk(:,17) = [0.0,1.0,1.0]*0.00625*16.0
      !qpts_loc%bk(:,18) = [0.0,1.0,1.0]*0.00625*17.0
      !qpts_loc%bk(:,19) = [0.0,1.0,1.0]*0.00625*18.0
      !qpts_loc%bk(:,20) = [0.0,1.0,1.0]*0.00625*19.0
      !qpts_loc%bk(:,21) = [0.0,1.0,1.0]*0.00625*20.0

      ! Read q-Points from inp.xml!
      qpts_loc%bk(:, :SIZE(fi%juPhon%qvec,2)) = fi%juPhon%qvec

      ALLOCATE(q_list(SIZE(fi%juPhon%qvec,2)))!,dfpt_eig_id_list(SIZE(fi%juPhon%qvec,2)))
      q_list = (/(iArray, iArray=1,SIZE(fi%juPhon%qvec,2), 1)/)

      !ALLOCATE(dfpt_eig_id_list2,mold=dfpt_eig_id_list)
      !IF (l_minusq) ALLOCATE(dfpt_eigm_id_list2,mold=dfpt_eigm_id_list)

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
         write(oUnit, *) "grVext", iDir
         CALL vgen_coulomb(1, fmpi_nosym, fi_nosym%input, fi_nosym%field, fi_nosym%vacuum, fi_nosym%sym, stars_nosym, fi_nosym%cell, &
                         & sphhar_nosym, fi_nosym%atoms, .FALSE., imagrhodummy, grVext3(iDir), &
                         & dfptdenimag=imagrhodummy, dfptvCoulimag=grvextdummy,dfptden0=imagrhodummy,stars2=stars_nosym,iDtype=0,iDir=iDir)
      END DO

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
            END DO
         END DO
      END DO

      CALL grRho3(1)%distribute(fmpi%mpi_comm)
      CALL grRho3(2)%distribute(fmpi%mpi_comm)
      CALL grRho3(3)%distribute(fmpi%mpi_comm)
      CALL grVext3(1)%distribute(fmpi%mpi_comm)
      CALL grVext3(2)%distribute(fmpi%mpi_comm)
      CALL grVext3(3)%distribute(fmpi%mpi_comm)

      DO iDir = 1, 3
         CALL sh_to_lh(fi_nosym%sym, fi_nosym%atoms, sphhar_nosym, SIZE(rho_nosym%mt,4), 2, grrhodummy(:, :, :, :, iDir), grRho3(iDir)%mt, imagrhodummy%mt)
         CALL imagrhodummy%resetPotDen()
         write(oUnit, *) "grVeff", iDir
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi_nosym%cell, fi_nosym%sliceplot, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVtot3(iDir), .TRUE., grvextdummy, grRho3(iDir), 0, iDir, [0,0])
         write(oUnit, *) "grVC", iDir
         CALL dfpt_vgen(hybdat_nosym, fi_nosym%field, fi_nosym%input, xcpot_nosym, fi_nosym%atoms, sphhar_nosym, stars_nosym, fi_nosym%vacuum, fi_nosym%sym, &
                        fi_nosym%cell, fi_nosym%sliceplot, fmpi_nosym, fi_nosym%noco, nococonv_nosym, rho_nosym, vTot_nosym, &
                        stars_nosym, imagrhodummy, grVC3(iDir), .FALSE., grvextdummy, grRho3(iDir), 0, iDir, [0,0])
      END DO

      CALL genPertPotDensGvecs( stars_nosym, fi_nosym%cell, fi_nosym%input, ngdp, ngdp2km, [0.0,0.0,0.0], recG )

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
      DO iQ = 1, SIZE(q_list)
         CALL timestart("q-point")
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

         CALL timestart("Eigenstuff at k+q")
         !q_eig_id = open_eig(fmpi%mpi_comm, lapw_dim_nbasfcn, fi%input%neig, fi%kpts%nkpt, fi%input%jspins, fi%noco%l_noco, &
         !                  .NOT.fi%INPUT%eig66(1), fi%input%l_real, fi%noco%l_soc, fi%input%eig66(1), .FALSE., fmpi%n_size)

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

!#ifdef CPP_MPI
!               CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
!               CALL starsq%mpi_bc(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "starsq OK"
!               CALL denIn1%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "denIn OK"
!               CALL denIn1Im%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "denInIm OK"
!               CALL vTot1%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "vTot OK"
!               CALL vTot1Im%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "vTotIm OK"
!               CALL vC1%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "vC OK"
!               CALL vC1Im%distribute(fmpi%mpi_comm)
!               write(*,*) fmpi_nosym%irank, "vCIm OK"
!               CALL MPI_BARRIER(fmpi%mpi_comm, ierr)
!#endif

               IF (fmpi%irank==0) WRITE(*,*) '-------------------------'
               ! This is where the magic happens. The Sternheimer equation is solved
               ! iteratively, providing the scf part of dfpt calculations.
               IF (l_minusq) THEN
                  CALL timestart("Sternheimer with -q")
                  CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, q_results, enpara_nosym, hybdat_nosym, mpdata_nosym, forcetheo_nosym, &
                                        rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), grVC3(iDir), q_list(iQ), iDtype, iDir, &
                                        dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                        denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, &
                                        starsmq, qm_results, dfpt_eigm_id, dfpt_eigm_id2, qm_eig_id, results1m, vTot1m, vTot1mIm)
                  CALL timestop("Sternheimer with -q")
               ELSE
                  CALL timestart("Sternheimer")
                  CALL dfpt_sternheimer(fi_nosym, xcpot_nosym, sphhar_nosym, stars_nosym, starsq, nococonv_nosym, qpts_loc, fmpi_nosym, results_nosym, q_results, enpara_nosym, hybdat_nosym, mpdata_nosym, forcetheo_nosym, &
                                        rho_nosym, vTot_nosym, grRho3(iDir), grVtot3(iDir), grVext3(iDir), grVC3(iDir), q_list(iQ), iDtype, iDir, &
                                        dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                        denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im)
                  CALL timestop("Sternheimer")
               END IF

               IF (.FALSE.) THEN
                  CALL timestop("Dirloop")
                  CYCLE
               END IF

               IF (fmpi%irank==0) WRITE(*,*) '-------------------------'
               CALL timestart("Dynmat")
               ! Once the first order quantities are converged, we can construct all
               ! additional necessary quantities and from that the dynamical matrix.
               IF (.TRUE.) THEN
                  CALL dfpt_dynmat_row(fi_nosym, stars_nosym, starsq, sphhar_nosym, xcpot_nosym, nococonv_nosym, hybdat_nosym, fmpi_nosym, qpts_loc, q_list(iQ), iDtype, iDir, &
                                       eig_id, dfpt_eig_id, dfpt_eig_id2, enpara_nosym, mpdata_nosym, results_nosym, results1, l_real,&
                                       rho_nosym, vTot_nosym, grRho3, grVext3, grVC3, grVtot3, &
                                       denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), .TRUE., .TRUE., E2ndOrdII)
               ELSE
                  CALL dfpt_dynmat_row(fi_nosym, stars_nosym, starsq, sphhar_nosym, xcpot_nosym, nococonv_nosym, hybdat_nosym, fmpi_nosym, qpts_loc, q_list(iQ), iDtype, iDir, &
                                       eig_id, dfpt_eig_id, dfpt_eig_id2, enpara_nosym, mpdata_nosym, results_nosym, results1, l_real,&
                                       rho_nosym, vTot_nosym, grRho3, grVext3, grVC3, grVtot3, &
                                       denIn1, vTot1, denIn1Im, vTot1Im, vC1, vC1Im, dyn_mat(iQ,3 *(iDtype-1)+iDir,:), .TRUE., .TRUE., E2ndOrdII, q_eig_id)
               END IF
               CALL timestop("Dynmat")
               dyn_mat(iQ,3 *(iDtype-1)+iDir,:) = dyn_mat(iQ,3 *(iDtype-1)+iDir,:) + conjg(E2ndOrdII(3 *(iDtype-1)+iDir,:))
               !IF (fmpi%irank==0) write(9989,*) "Eii2:", E2ndOrdII(3 *(iDtype-1)+iDir,:)
               !IF (fmpi%irank==0) write(9990,*) "Eii2:", E2ndOrdII(3 *(iDtype-1)+iDir,:)
               IF (fmpi%irank==0) write(*,*) "dynmat row for ", dfpt_tag
               IF (fmpi%irank==0) write(*,*) dyn_mat(iQ,3 *(iDtype-1)+iDir,:)
               !STOP
               CALL timestop("Dirloop")
            END DO
            CALL timestop("Typeloop")

#if defined(CPP_MPI)
            CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif
         END DO


         IF (.FALSE.) THEN
            CALL timestop("q-point")
            !CALL close_eig(q_eig_id)
            CYCLE
         END IF

         IF (fmpi%irank==0) THEN
            WRITE(*,*) '-------------------------'
            CALL timestart("Dynmat diagonalization")
            CALL DiagonalizeDynMat(fi_nosym%atoms, qpts_loc, fi%juPhon%calcEigenVec, dyn_mat(iQ,:,:), eigenVals, eigenVecs, q_list(iQ))
            CALL timestop("Dynmat diagonalization")

            CALL timestart("Frequency calculation")
            CALL CalculateFrequencies(fi_nosym%atoms, q_list(iQ), eigenVals, eigenFreqs)
            CALL timestop("Frequency calculation")
            DEALLOCATE(eigenVals, eigenVecs, eigenFreqs, E2ndOrdII)
         END IF

         !CALL close_eig(q_eig_id)
         !IF (l_minusq) CALL close_eig(qm_eig_id)
         CALL timestop("q-point")

      END DO

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

        CALL juDFT_end("Phonon calculation finished.")

    END SUBROUTINE dfpt

    subroutine mt_gradient_new(atoms, sphhar, sym, r2FlhMt, GrFshMt)

      use m_gaunt, only : gaunt1

      type(t_atoms),               intent(in)  :: atoms
      type(t_sphhar),              intent(in)  :: sphhar
      type(t_sym),                 intent(in)  :: sym

      real,                        intent(in)  :: r2FlhMt(:, 0:, :)
      complex,                     intent(out) :: GrFshMt(:, :, :, :)

      real                                     :: pfac
      real                                     :: tGaunt
      integer                                  :: itype
      integer                                  :: imesh
      integer                                  :: mqn_m
      integer                                  :: oqn_l
      integer                                  :: mqn_mpp
      integer                                  :: lm
      integer                                  :: symType
      integer                                  :: ilh
      integer                                  :: imem

      real,           allocatable              :: rDerFlhMt(:)
      complex,        allocatable              :: r2GrFshMtNat(:, :, :, :)

      allocate( r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
      allocate( rDerFlhMt(atoms%jmtd) )
      GrFshMt = cmplx(0., 0.)
      r2GrFshMtNat = cmplx(0., 0.)
      rDerFlhMt = 0.

      pfac = sqrt( fpi_const / 3. )
      do mqn_mpp = -1, 1
        do itype = 1, atoms%ntype
            symType = sym%ntypsy(itype)
            do ilh = 0, sphhar%nlh(symType)
              oqn_l = sphhar%llh(ilh, symType)
              do imem = 1, sphhar%nmem(ilh,symType)
                mqn_m = sphhar%mlh(imem,ilh,symType)

                ! l + 1 block
                ! oqn_l - 1 to l, so oqn_l should be < lmax not <= lmax
                if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                  lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                  call derivative_loc( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                  tGaunt = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                  do imesh = 1, atoms%jri(itype)
                    r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                      &* tGaunt * (rDerFlhMt(imesh) * sphhar%clnu(imem,ilh,symType) &
                      &- ((oqn_l + 2) * r2FlhMt(imesh, ilh, itype) * sphhar%clnu(imem,ilh,symType) / atoms%rmsh(imesh, itype)))
                  end do ! imesh
                end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l )

                ! l - 1 block
                if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                  if ( oqn_l - 1 == -1 ) then
                    write (*, *) 'oqn_l too low'
                  end if
                  lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                  ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                  call derivative_loc( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                  tGaunt = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                  do imesh = 1, atoms%jri(itype)
                    r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                      & * tGaunt * (rDerFlhMt(imesh)  * sphhar%clnu(imem,ilh,symType) &
                      & + ((oqn_l - 1) * r2FlhMt(imesh, ilh, itype) * sphhar%clnu(imem,ilh,symType) / atoms%rmsh(imesh, itype)))
                  end do ! imesh
                end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l )
              end do ! imem
            end do ! ilh
        end do ! itype
      end do ! mqn_mpp

      ! Conversion from natural to cartesian coordinates
      do itype = 1, atoms%ntype
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                grFshMt(imesh, lm, itype, 1:3) = matmul( Tmatrix0(1:3, 1:3), r2GrFshMtNat(imesh, lm, itype, 1:3) ) / atoms%rmsh(imesh, itype)**2
              end do
            end do ! mqn_m
          end do ! oqn_l
      end do ! itype

   end subroutine mt_gradient_new

   subroutine derivative_loc(f, itype, atoms, df)

      integer,       intent(in)  :: itype
      type(t_atoms), intent(in)  :: atoms
      real,          intent(in)  :: f(atoms%jri(itype))
      real,          intent(out) :: df(atoms%jri(itype))
      real                       :: h, r, d21, d32, d43, d31, d42, d41, df1, df2, s
      real                       :: y0, y1, y2
      integer                    :: i, n

      n = atoms%jri(itype)
      h = atoms%dx(itype)
      r = atoms%rmsh(1, itype)

      ! use Lagrange interpolation of 3rd order (and averaging) for points 3 to n
      d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)
      d31 = d21 + d32      ; d42 = d32 + d43
      d41 = d31 + d43
      df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1d0/d21 - 1d0/d31 - 1d0/d41) * f(1)&
     &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)
      df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1d0/d21 - 1d0/d32 - 1d0/d42) * f(2)&
     &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)
      df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +&
     &  ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(3) + d31*d32 / (d41*d42*d43) * f(4)
      do i = 3, n - 2
         d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
         d31 = d42 ; d42 = d42 * exp(h)
         d41 = d41 * exp(h)
         df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1d0/d21 - 1d0/d32 - 1d0/d42) * f(i) + &
     &             d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) * f(i+2)
         df(i) = ( df1 + df2 ) / 2
         df1   = d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42) * f(i) +&
     &    ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(i+1) + d31*d32 / (d41*d42*d43) * f(i+2)
      enddo
      df(n-1) = df1
      df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -&
     &            d41*d42 / (d31*d32*d43) * f(n-1) + ( 1d0/d41 + 1d0/d42 + 1d0/d43 ) * f(n)
      ! for first two points use Lagrange interpolation of second order for log(f(i))
      ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
      s = sign(1d0,f(1))
      if(sign(1d0,f(2)) /= s .or. sign(1d0,f(3))  /= s .or. any(abs(f(:3)) < 1e0)) then
         d21   = r * (exp(h)-1)
         d32   = d21 * exp(h)
         d31   = d21 + d32
         s     = df(3) / (d31*d32) - f(1) / (d21*d31**2) + f(2) / (d21*d32**2) - f(3) / (d31**2*d32) - f(3) / (d31*d32**2)
         df(1) = - (d21+d31) / (d21*d31) * f(1) + d31 / (d21*d32) * f(2) - d21 / (d31*d32) * f(3) + d21*d31 * s

         df(2) = - (d21-d32) / (d21*d32) * f(2) - d32 / (d21*d31) * f(1) + d21 / (d31*d32) * f(3) - d21*d32 * s
      else
         y0    = log(abs(f(1)))
         y1    = log(abs(f(2)))
         y2    = log(abs(f(3)))
         df(1) = ( - 3*y0/2 + 2*y1 - y2/2 ) * f(1) / (h*r)
         df(2) = (y2-y0)/2                  * f(2) / (h*r*exp(h))
      endif
   end subroutine derivative_loc

    SUBROUTINE sh_to_lh(sym, atoms, sphhar, jspins, radfact, rhosh, rholhreal, rholhimag)

        ! WARNING: This routine will not fold back correctly for activated sym-
        !          metry and gradients (rho in l=0 and lattice harmonics do not
        !          allow l=1 --> gradient in l=1 is lost)

        TYPE(t_sym),    INTENT(IN)  :: sym
        TYPE(t_atoms),  INTENT(IN)  :: atoms
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        INTEGER,        INTENT(IN)  :: jspins, radfact
        COMPLEX,        INTENT(IN)  :: rhosh(:, :, :, :)
        REAL,           INTENT(OUT) :: rholhreal(:, 0:, :, :), rholhimag(:, 0:, :, :)

        INTEGER :: iSpin, iType, iEqat, iAtom, ilh, iMem, ilm, iR
        INTEGER :: ptsym, l, m
        REAL    :: factor

        rholhreal = 0.0
        rholhimag = 0.0

        DO iSpin = 1, jspins
            DO iType = 1, atoms%ntype
                iAtom = atoms%firstAtom(iType)
                ptsym = sym%ntypsy(iAtom)
                DO ilh = 0, sphhar%nlh(ptsym)
                    l = sphhar%llh(iLH, ptsym)
                    DO iMem = 1, sphhar%nmem(ilh, ptsym)
                        m = sphhar%mlh(iMem, ilh, ptsym)
                        ilm = l * (l+1) + m + 1
                        DO iR = 1, atoms%jri(iType)
                           IF ((radfact.EQ.0).AND.(l.EQ.0)) THEN
                               factor = atoms%rmsh(iR, iType) / sfp_const
                           ELSE IF (radfact.EQ.2) THEN
                               factor = atoms%rmsh(iR, iType)**2
                           ELSE
                               factor = 1.0
                           END IF
                            rholhreal(iR, ilh, iType, iSpin) = &
                          & rholhreal(iR, ilh, iType, iSpin) + &
                          &  real(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                            rholhimag(iR, ilh, iType, iSpin) = &
                          & rholhimag(iR, ilh, iType, iSpin) + &
                          & aimag(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE sh_to_lh
END MODULE m_dfpt
