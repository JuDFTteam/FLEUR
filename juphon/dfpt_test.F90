!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_test

   USE m_types
   USE m_vgen_coulomb
   USE m_vgen_xcpot
   USE m_jpGrVeff0, ONLY : GenGrVeff0
   USE m_dfpt_init
   USE m_npy

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_test(fi, sphhar, stars, fmpi, rho, grRho, rho0, grRho0, xcpot, ngdp, gdp, grVxcIRKern, ylm, dKernMTGPts, gWghts, hybdat)
      ! Test various of the newly written FLEUR DFPT subroutines for sensibility
      ! and, where applicable, against old juPhon routines.

      TYPE(t_fleurinput), INTENT(IN)  :: fi
      TYPE(t_sphhar),     INTENT(IN)  :: sphhar
      TYPE(t_stars),      INTENT(IN)  :: stars
      TYPE(t_mpi),        INTENT(IN)  :: fmpi
      TYPE(t_potden),     INTENT(IN)  :: rho, grRho
      TYPE(t_jpPotden),   INTENT(IN)  :: rho0, grRho0
      CLASS(t_xcpot),     INTENT(IN)  :: xcpot
      TYPE(t_hybdat),     INTENT(IN)  :: hybdat

      INTEGER,            INTENT(IN)  :: ngdp, gdp(:, :)

      COMPLEX,            INTENT(IN)  :: grVxcIRKern(:), ylm(:, :)
      REAL,               INTENT(IN)  :: dKernMTGPts(:, :, :), gWghts(:)

      LOGICAL :: harSw, extSw, xcSw, testGoldstein, grRhoTermSw

      COMPLEX, ALLOCATABLE :: grVeff0IRDummy(:, :), grVeff0MT_juPhon(:, :, :, :)

      TYPE(t_potden) :: workdenReal, workdenImag, vCoul, dfptvCoulimag, dummyrho, dummyrho2, dummyrho3, dummyrho4
      TYPE(t_jpPotden) :: dummyPot

      CALL workdenReal%copyPotDen(rho)
      CALL workdenImag%copyPotDen(rho)
      CALL vCoul%copyPotDen(rho)
      CALL dfptvCoulimag%copyPotDen(rho)
      dummyPot = rho0
      dummyrho = rho
      dummyrho2 = rho
      dummyrho3 = rho
      dummyrho4 = rho

      workdenReal%mt = 0.0
      workdenReal%pw = CMPLX(0.0,0.0)
      workdenImag%mt = 0.0
      workdenImag%pw = CMPLX(0.0,0.0)
      vCoul%mt = 0.0
      vCoul%pw = CMPLX(0.0,0.0)
      dfptvCoulimag%mt = 0.0
      dfptvCoulimag%pw = CMPLX(0.0,0.0)
      dummyPot%mt = CMPLX(0.0,0.0)
      dummyPot%pw = CMPLX(0.0,0.0)
      dummyrho%mt = 0.0
      dummyrho%pw = CMPLX(0.0,0.0)
      dummyrho2%mt = 0.0
      dummyrho2%pw = CMPLX(0.0,0.0)

      ! Compare the Coulomb potential gradients.
      ! First: gradRho = 0, SF contributions on.
      harSw         = .TRUE.
      extSw         = .TRUE.
      xcSw          = .FALSE.
      testGoldstein = .FALSE.
      grRhoTermSw   = .TRUE.

      CALL GenGrVeff0(fi%atoms, fi%cell, stars, ngdp, harSw, extSw, xcSw, gdp, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                                 & 0.0*grRho0%pw(:, 1, 1 ,:), 0.0*grRho0%mt(:, :, :, 1, 1, :), gWghts, ylm, dKernMTGPts, grVxcIRKern, &
                                 & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MT_juPhon )

      CALL vgen_coulomb(1, fmpi, fi%input, fi%field, fi%vacuum, fi%sym, stars, fi%cell, &
                      & sphhar, fi%atoms, .FALSE., workdenReal, vCoul, &
                      & dfptdenimag=workdenImag, dfptvCoulimag=dfptvCoulimag,dfptden0=rho,stars2=stars,iDtype=0,iDir=1)

      CALL lh_to_sh(fi%sym, fi%atoms, sphhar, 1, 1, vCoul%mt, dummypot%mt)

      !!! Check: dummypot%%mt(:,:,1,1,1,1) should be equal to grVeff0MT_juPhon(:,:,1,1)
      !!! If so: SF contributions and Vext part are correct.
      !CALL save_npy("ir_juph.npy",grVeff0IRdummy(:,1))
      !CALL save_npy("mt_juph.npy",grVeff0MT_juPhon(:,:,1,1))
      !CALL save_npy("ir_maxx.npy",vCoul%pw(:,1))
      !CALL save_npy("sh_maxx.npy",dummypot%mt(:,:,1,1,1,1))
      DEALLOCATE (grVeff0IRdummy)
      DEALLOCATE(grVeff0MT_juPhon)

      ! Now: No SF contributions, gradRho /= 0.
      vCoul%mt = 0.0
      vCoul%pw = CMPLX(0.0,0.0)
      dfptvCoulimag%mt = 0.0
      dfptvCoulimag%pw = CMPLX(0.0,0.0)

      harSw         = .TRUE.
      extSw         = .TRUE.
      xcSw          = .FALSE.
      testGoldstein = .FALSE.
      grRhoTermSw   = .TRUE.

      CALL GenGrVeff0(fi%atoms, fi%cell, stars, ngdp, harSw, extSw, xcSw, gdp, 0.0*rho0%pw(:, :, 1 ,1), 0.0*rho0%mt(:, :, :, :, 1 ,1), &
                                 & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gWghts, ylm, dKernMTGPts, grVxcIRKern, &
                                 & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MT_juPhon )

      workdenReal%pw = grRho%pw
      CALL sh_to_lh(fi%sym, fi%atoms, sphhar, 1, 2, grRho0%mt(:, :, :, :, 1, 1), workdenReal%mt, workdenImag%mt)

      CALL vgen_coulomb(1, fmpi,  fi%input, fi%field, fi%vacuum, fi%sym, stars, fi%cell, &
                      & sphhar, fi%atoms, .FALSE., workdenReal, vCoul, &
                      & dfptdenimag=workdenImag, dfptvCoulimag=dfptvCoulimag,dfptden0=dummyrho,stars2=stars,iDtype=0,iDir=1)

      CALL lh_to_sh(fi%sym, fi%atoms, sphhar, 1, 1, vCoul%mt, dummypot%mt)

      !!! Check: dummypot%%mt(:,:,1,1,1,1) should be equal to grVeff0MT_juPhon(:,:,1,1)
      !!! If so: gradRho contributions and Vext part are correct.
      !CALL save_npy("ir_juph.npy",grVeff0IRdummy(:,1))
      !CALL save_npy("mt_juph.npy",grVeff0MT_juPhon(:,:,1,1))
      !CALL save_npy("ir_maxx.npy",vCoul%pw(:,1))
      !CALL save_npy("sh_maxx.npy",dummypot%mt(:,:,1,1,1,1))
      DEALLOCATE (grVeff0IRdummy)
      DEALLOCATE(grVeff0MT_juPhon)

      ! Everything is fine up to here, so we seem to be able to generate gradVH from gradRho.
      ! Now try xcpot.

      vCoul%mt = 0.0
      vCoul%pw = CMPLX(0.0,0.0)
      dfptvCoulimag%mt = 0.0
      dfptvCoulimag%pw = CMPLX(0.0,0.0)

      harSw         = .FALSE.
      extSw         = .FALSE.
      xcSw          = .TRUE.
      testGoldstein = .FALSE.
      grRhoTermSw   = .TRUE.

      CALL GenGrVeff0(fi%atoms, fi%cell, stars, ngdp, harSw, extSw, xcSw, gdp, rho0%pw(:, :, 1 ,1), rho0%mt(:, :, :, :, 1 ,1), &
                                 & grRho0%pw(:, 1, 1 ,:), grRho0%mt(:, :, :, 1, 1, :), gWghts, ylm, dKernMTGPts, grVxcIRKern, &
                                 & testGoldstein, grRhoTermSw, grVeff0IRdummy, grVeff0MT_juPhon )

      CALL vgen_xcpot(hybdat,fi%input,xcpot,fi%atoms,sphhar,stars,fi%vacuum,fi%sym,&
                     fi%cell,fi%sliceplot,fmpi,fi%noco,rho,rho,dummyrho,vCoul,dummyrho2,dummyrho3,dummyrho4, &
                     & den1Rot=workdenReal, den1Rotimag=workdenImag, dfptvTotimag=dfptvCoulimag,starsq=stars)

      CALL lh_to_sh(fi%sym, fi%atoms, sphhar, 1, 1, vCoul%mt, dummypot%mt)

      !!! Check: dummypot%%mt(:,:,1,1,1,1) should be equal to grVeff0MT_juPhon(:,:,1,1)
      !!! If so: gradRho contributions of Vxc part are correct.
      !CALL save_npy("ir_juph.npy",grVeff0IRdummy(:,1))
      !CALL save_npy("mt_juph.npy",grVeff0MT_juPhon(:,:,1,1))
      !CALL save_npy("ir_maxx.npy",vCoul%pw(:,1))
      !CALL save_npy("sh_maxx.npy",dummypot%mt(:,:,1,1,1,1))
      DEALLOCATE (grVeff0IRdummy)
      DEALLOCATE(grVeff0MT_juPhon)

      ! This also works up to a factor of 1.5 in the xc potential (DUH.)
      ! Note, that the pw part is shifted around due to mxmymz(jp) vs mzmymx(fleur)

      ! Further test ideas:

      ! V1 tests: wait for starsq to be implemented/test with gradRho for now.
      !harSw    = .FALSE.
      !extSw    = .TRUE.
      !xcSw     = .FALSE.
      !vExtFull = .TRUE.
      !vHarNum  = .FALSE.
      !rho1pw   = ???
      !rho1MT   = ???
      !CALL GenVeff1(fi%input, stars, fi%cell, fi%atoms, harSw, extSw, xcSw, vExtFull, &
      !            & ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, grRho0MT, &
      !            & gdp, vEff1IRsh, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, &
      !            & iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )
      !CALL vgen_xcpot(hybdat,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
      !               cell ,sliceplot,fmpi,noco,den,denRot,EnergyDen,dfptvTot,vx,vxc,exc, &
      !               & den1Rotimag=den1imRot, dfptvTotimag=dfptvTotimag,starsq=starsq)

      CALL juDFT_end("Phonon tests finished.")

   END SUBROUTINE dfpt_test
END MODULE m_dfpt_test
