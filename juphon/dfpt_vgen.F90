!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vgen
   USE m_juDFT

CONTAINS

   SUBROUTINE dfpt_vgen(hybdat,field,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                   cell,oneD,sliceplot,fmpi,results,noco,nococonv,EnergyDen,den,vTot,vx,vCoul,vxc,exc,&
                   &starsq,dfptdenimag,dfptvTotimag,dfptdenreal,iDtype,iDir)
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
      USE m_magnMomFromDen
      USE m_force_sf ! Klueppelberg (force level 3)
      USE m_fleur_vdW
      IMPLICIT NONE

      TYPE(t_results),   INTENT(INOUT) :: results
      CLASS(t_xcpot),    INTENT(IN)    :: xcpot
      TYPE(t_hybdat),    INTENT(IN)    :: hybdat
      TYPE(t_mpi),       INTENT(IN)    :: fmpi
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_field),     INTENT(IN)    :: field
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_nococonv),  INTENT(INOUT) :: nococonv
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_potden),    INTENT(IN)    :: EnergyDen
      TYPE(t_potden),    INTENT(INOUT) :: den
      TYPE(t_potden),    INTENT(INOUT) :: vTot, vx, vxc, exc

      TYPE(t_stars),  OPTIONAL, INTENT(IN)    :: starsq
      TYPE(t_potden), OPTIONAL, INTENT(INOUT) :: dfptdenimag, dfptvTotimag, dfptdenreal

      INTEGER, OPTIONAL, INTENT(IN)           :: iDtype, iDir ! DFPT: Type and direction of displaced atom

      TYPE(t_potden)                   :: workden, denRot, workdenImag, workdenReal, vCoul, dfptvCoulimag

      INTEGER :: i, js
      REAL    :: b(3,atoms%ntype), dummy1(atoms%ntype), dummy2(atoms%ntype)
      LOGICAL :: l_dfptvgen ! If this is true, we handle things differently!

      l_dfptvgen = PRESENT(starsq)

      IF (fmpi%irank == 0) THEN
         IF (noco%l_sourceFree) THEN
            CALL magnMomFromDen(input,atoms,noco,den,b,dummy1,dummy2)
            DO i=1,atoms%ntype
               WRITE  (oUnit,8025) i,b(1,i),b(2,i),b(3,i),SQRT(b(1,i)**2+b(2,i)**2+b(3,i)**2)
               8025 FORMAT(2x,'Magmom before SF [local frame, atom ',i2,']: ','mx=',f9.5,' my=',f9.5,' mz=',f9.5,' |m|=',f9.5)
            END DO
         END IF
      END IF

      IF (fmpi%irank==0) WRITE (oUnit,FMT=8000)
      IF ((fmpi%irank==0).AND.l_dfptvgen) WRITE (oUnit,FMT=8001)
8000  FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
8001  FORMAT (/,/,t10,'          (DFPT edition)              ',/)
      CALL vTot%resetPotDen()
      CALL vCoul%resetPotDen()
      CALL vx%resetPotDen()
      CALL vxc%resetPotDen()
      CALL exc%resetPotDen()

      ALLOCATE(vx%pw_w,mold=vTot%pw)
      vx%pw_w = 0.0
      ALLOCATE(vxc%pw_w,mold=vTot%pw)
      vxc%pw_w = 0.0
      CALL exc%init(stars, atoms, sphhar, vacuum, noco, 1, 1) !one spin only
      ALLOCATE (exc%pw_w(stars%ng3, 1)); exc%pw_w = 0.0

#ifndef CPP_OLDINTEL
      ALLOCATE(vTot%pw_w,mold=vTot%pw)
#else
      ALLOCATE( vTot%pw_w(size(vTot%pw,1),size(vTot%pw,2)))
#endif

      ALLOCATE(vCoul%pw_w(SIZE(vCoul%pw,1),size(vCoul%pw,2)))
      vCoul%pw_w = CMPLX(0.0,0.0)

      results%force=0.0

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
        CALL vgen_coulomb(1,fmpi,oneD,input,field,vacuum,sym,starsq,cell,sphhar,atoms,.FALSE.,workdenReal,vCoul,&
                        & dfptdenimag=workdenImag,dfptvCoulimag=dfptvCoulimag,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)

      ! b)
      CALL vCoul%copy_both_spin(vTot)
      CALL dfptvCoulimag%copy_both_spin(dfptvTotimag)

      ! c)
      IF (noco%l_noco) THEN
         CALL denRot%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
         denRot=den
         CALL rotate_int_den_to_local(sym,stars,atoms,sphhar,vacuum,cell,input,noco,oneD,denRot)
         IF (any(noco%l_unrestrictMT)) CALL rotate_mt_den_to_local(atoms,sphhar,sym,noco,denrot)
      END IF

          IF (noco%l_noco) THEN
              !Functions that construct the spin-dependent perturbed densities
              !from the perturbed charge and (vectorial) magnetization density/
              !perturbed density matrix. Also saves the perturbed angles.
              !TODO: Calculate in real space but put back onto coefficients, just
              !      like in normal scf.
              !!CALL get_int_local_perturbation()
              !!IF (any(noco%l_unrestrictMT)) CALL get_mt_local_perturbation()
          END IF
          CALL vgen_xcpot(hybdat,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                          cell,oneD,sliceplot,fmpi,noco,den,denRot,EnergyDen,vTot,vx,vxc,exc)

      ! d)
      ! TODO: This is so different from the base case, that we build a new suboutine.
      !!CALL dfpt_vgen_finalize(fmpi,oneD,field,cell,atoms,stars,vacuum,sym,noco,nococonv,input,xcpot,sphhar,vTot,denRot,sliceplot)
      !DEALLOCATE(vcoul%pw_w)

      CALL vTot%distribute(fmpi%mpi_comm)
      CALL vx%distribute(fmpi%mpi_comm)
      CALL vxc%distribute(fmpi%mpi_comm)
      CALL exc%distribute(fmpi%mpi_comm)

      ! Klueppelberg (force level 3)
      IF (input%l_f.AND.(input%f_level.GE.3).AND.(fmpi%irank.EQ.0)) THEN
         DO js = 1,input%jspins
            CALL force_sf_is(atoms,stars,sym,js,cell,den%pw,vTot%pw,exc%pw(:,1),vxc%pw)
         END DO

      END IF

  END SUBROUTINE dfpt_vgen

END MODULE m_dfpt_vgen
