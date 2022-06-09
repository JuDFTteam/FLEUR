!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vgen
   USE m_juDFT

CONTAINS

   SUBROUTINE dfpt_vgen(hybdat,field,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                   cell ,sliceplot,fmpi,noco,nococonv,den,vTot,&
                   &starsq,dfptdenimag,dfptvTot,dfptvTotimag,dfptdenreal,iDtype,iDir)
      !--------------------------------------------------------------------------
      ! FLAPW potential perturbation generator (main routine)
      !
      ! Modification for use with DFPT:
      ! The density variables in the interstitial now live in a G-expansion
      ! shifted by q and those in the Muffin Tin are no longer real. Account for
      ! that with optional arguments: use q-shifted stars and carry the imaginary
      ! part of the MT-density along explicitely. Changes:
      ! - The interstitial part and MT real part is in dfptdenreal
      ! - The MT imaginary part is in dfptdenimag
      ! - vTot will carry the same quantities for V1 as dfptdenreal for rho1
      ! - The MT imaginary part is in dfptvTotimag
      ! - iDtype and iDir tell us where we perturb (atom and direction)
      ! - den is still den; we need it for additional qlm of surface contributions
      !--------------------------------------------------------------------------

      USE m_types
      USE m_constants
      USE m_rotate_int_den_tofrom_local
      USE m_bfield
      USE m_vgen_coulomb
      USE m_vgen_xcpot
      USE m_vgen_finalize
      USE m_rotate_mt_den_tofrom_local
      USE m_get_int_perturbation
      USE m_get_mt_perturbation
      USE m_dfpt_vgen_finalize

      IMPLICIT NONE

      CLASS(t_xcpot),    INTENT(IN)    :: xcpot
      TYPE(t_hybdat),    INTENT(IN)    :: hybdat
      TYPE(t_mpi),       INTENT(IN)    :: fmpi

      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_field),     INTENT(IN)    :: field
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_nococonv),  INTENT(IN)    :: nococonv
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_potden),    INTENT(IN)    :: vTot
      TYPE(t_potden),    INTENT(INOUT) :: den, dfptvTot

      TYPE(t_stars),  OPTIONAL, INTENT(IN)    :: starsq
      TYPE(t_potden), OPTIONAL, INTENT(INOUT) :: dfptdenimag, dfptvTotimag, dfptdenreal

      INTEGER, OPTIONAL, INTENT(IN)           :: iDtype, iDir ! DFPT: Type and direction of displaced atom

      TYPE(t_potden)                   :: workden, denRot, workdenImag, workdenReal, den1Rot, den1imRot
      TYPE(t_potden)                   :: vCoul, dfptvCoulimag, vxc, exc, vx, EnergyDen

      INTEGER :: i, js
      REAL    :: b(3,atoms%ntype), dummy1(atoms%ntype), dummy2(atoms%ntype)

      vCoul = vTot
      vx = vTot
      vxc = vTot
      exc = vTot
      dfptvCoulimag = vTot

      IF (fmpi%irank==0) WRITE (oUnit,FMT=8000)
      IF (fmpi%irank==0) WRITE (oUnit,FMT=8001)
8000  FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
8001  FORMAT (/,/,t10,'          (DFPT edition)              ',/)
      CALL dfptvTot%resetPotDen()
      CALL dfptvTotimag%resetPotDen()
      CALL vCoul%resetPotDen()
      CALL vx%resetPotDen()
      CALL vxc%resetPotDen()
      CALL exc%resetPotDen()
      CALL dfptvCoulimag%resetPotDen()

      ALLOCATE(vx%pw_w,mold=vTot%pw)
      vx%pw_w = 0.0
      ALLOCATE(vxc%pw_w,mold=vTot%pw)
      vxc%pw_w = 0.0
      CALL exc%init(stars, atoms, sphhar, vacuum, noco, 1, 1) !one spin only
      ALLOCATE (exc%pw_w(stars%ng3, 1)); exc%pw_w = 0.0

#ifndef CPP_OLDINTEL
      ALLOCATE(dfptvTot%pw_w,mold=dfptvTot%pw)
      ALLOCATE(dfptvTotimag%pw_w,mold=dfptvTotimag%pw)
#else
      ALLOCATE( dfptvTot%pw_w(size(dfptvTot%pw,1),size(dfptvTot%pw,2)))
      ALLOCATE( dfptvTotimag%pw_w(size(dfptvTotimag%pw,1),size(dfptvTotimag%pw,2)))
#endif

      ALLOCATE(vCoul%pw_w(SIZE(vCoul%pw,1),size(vCoul%pw,2)))
      vCoul%pw_w = CMPLX(0.0,0.0)

        CALL workDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
        CALL workDenReal%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
        CALL workDenImag%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)

        ! a)
        ! Sum up both spins in den into workden:
        CALL den%sum_both_spin(workden)
        CALL dfptdenreal%sum_both_spin(workdenReal)
        CALL dfptdenimag%sum_both_spin(workdenImag)
        ! TODO: Feeding starsq in instead of stars will be meaningless, unless
        !       we also add the q in question at the relevant points.
        ! NOTE: The normal stars are also passed as an optional argument, because
        !       they are needed for surface-qlm.
        CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,starsq,cell,sphhar,atoms,.FALSE.,workdenReal,vCoul,&
                        & dfptdenimag=workdenImag,dfptvCoulimag=dfptvCoulimag,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)

      ! b)
      CALL vCoul%copy_both_spin(dfptvTot)
      CALL dfptvCoulimag%copy_both_spin(dfptvTotimag)

      ! c)
      CALL denRot%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
      denRot=den
      CALL den1Rot%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
      CALL den1imRot%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
      den1Rot=dfptdenreal
      den1imRot=dfptdenimag
      IF (noco%l_noco) THEN
         CALL rotate_int_den_to_local(sym,stars,atoms,sphhar,vacuum,cell,input,noco ,denRot)
         IF (any(noco%l_unrestrictMT)) CALL rotate_mt_den_to_local(atoms,sphhar,sym,noco,denrot)
         !Functions that construct the spin-dependent perturbed densities
         !from the perturbed charge and (vectorial) magnetization density/
         !perturbed density matrix. Also saves the perturbed angles.
          CALL get_int_local_perturbation(sym, stars, atoms, sphhar, input, denRot, den1Rot, den1imRot, starsq)
          IF (any(noco%l_unrestrictMT)) CALL get_mt_local_perturbation(atoms,sphhar,sym,noco,denRot,den1Rot,den1imRot)
      END IF

         ! Skip vxc for rho(1)=0, i.e. starting potential
          IF (ANY(ABS(denRot%pw)>1E-12)) CALL vgen_xcpot(hybdat,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                          cell ,sliceplot,fmpi,noco,den,denRot,EnergyDen,dfptvTot,vx,vxc,exc, &
                          & den1Rot=den1Rot, den1Rotimag=den1imRot, dfptvTotimag=dfptvTotimag,starsq=starsq)

      IF (iDtype/=0) THEN
         ! d)
         ! TODO: This is so different from the base case, that we build a new suboutine.
         CALL dfpt_vgen_finalize(fmpi,atoms,stars,sym,noco,nococonv,input,sphhar,vTot,dfptvTot,dfptvTotimag,denRot,den1Rot,den1imRot,starsq)
         !DEALLOCATE(vcoul%pw_w)
      ELSE
         ! TODO: Write here something for the gradient. It does not need pw-stuff.
      END IF

      CALL dfptvTot%distribute(fmpi%mpi_comm)
      CALL dfptvTotimag%distribute(fmpi%mpi_comm)

  END SUBROUTINE dfpt_vgen

END MODULE m_dfpt_vgen
