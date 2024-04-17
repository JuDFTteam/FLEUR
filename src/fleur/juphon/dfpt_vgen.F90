!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vgen
   USE m_juDFT

CONTAINS

   SUBROUTINE dfpt_vgen(hybdat,field,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                   juphon,cell,fmpi,noco,nococonv,den,vTot,&
                   &starsq,dfptdenimag,dfptvTot,l_xc,dfptvTotimag,dfptdenreal,iDtype,iDir,killcont,sigma_disc,fi,do_vext,local_dfptvTot,local_dfptvTotimag,local_stars,local_starsq) !qvec
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

      USE m_dfpt_crank_gvecs
      USE m_types_fleurinput
      USE m_convn
      USE m_make_stars
      USE m_constants

      USE m_checkdopall !! Testing

      IMPLICIT NONE

      CLASS(t_xcpot),    INTENT(IN)    :: xcpot
      TYPE(t_hybdat),    INTENT(IN)    :: hybdat
      TYPE(t_mpi),       INTENT(IN)    :: fmpi

      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_field),     INTENT(IN)    :: field
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_nococonv),  INTENT(IN)    :: nococonv
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_juphon),    INTENT(IN)    :: juphon
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_potden),    INTENT(IN)    :: vTot
      TYPE(t_potden),    INTENT(INOUT) :: den, dfptvTot
      TYPE(t_potden),    OPTIONAL, INTENT(OUT) :: local_dfptvTot,local_dfptvTotimag
      TYPE(t_stars), OPTIONAL, INTENT(OUT)     :: local_stars,local_starsq
      !REAL, OPTIONAL, INTENT(IN) :: qvec(3)
      TYPE(t_fleurinput), OPTIONAL,INTENT(IN) :: fi
      TYPE(t_fleurinput) :: local_fi
      TYPE(t_potden) :: local_potden,local_potdenq 
      TYPE(t_potden) :: local_vCoul,local_workdenReal,local_workdenImag,local_workden,local_dfptvCoulimag
      TYPE(t_potden) :: local_dfptdenreal,local_den1Rot,local_den1imRot,local_den,local_denRot,local_dfptdenimag,local_vTot
      TYPE(t_atoms) :: local_atoms
      LOGICAL,OPTIONAL,INTENT(IN) :: do_vext
      LOGICAL :: l_vext
      INTEGER :: ispin
      REAL :: tmp_qvec(3)
      LOGICAL, INTENT(IN) :: l_xc

      TYPE(t_stars),  OPTIONAL, INTENT(IN)    :: starsq
      TYPE(t_potden), OPTIONAL, INTENT(INOUT) :: dfptdenimag, dfptvTotimag, dfptdenreal

      INTEGER, OPTIONAL, INTENT(IN)           :: iDtype, iDir ! DFPT: Type and direction of displaced atom

      INTEGER, OPTIONAL, INTENT(IN)           :: killcont(2)
      complex, OPTIONAL, INTENT(IN)           :: sigma_disc(2)

      TYPE(t_potden)                   :: workden, denRot, workdenImag, workdenReal, den1Rot, den1imRot
      TYPE(t_potden)                   :: vCoul, dfptvCoulimag, vxc, exc, vx, EnergyDen

      complex                           :: sigma_loc(2)

      vCoul = dfptvTot
      vx = vTot
      vxc = vTot
      exc = vTot
      dfptvCoulimag = dfptvTot

      IF (PRESENT(do_vext)) THEN
         l_vext = .TRUE.
      ELSE
         l_vext = .FALSE.
      END IF 
     
      IF (l_vext) THEN
         tmp_qvec = starsq%center
         write(4100,*) "Disc for Vext"
         CALL crank_gvecs(fi,fmpi,fi%sym,fi%cell,fi%input,sphhar, fi%vacuum , fi%noco ,local_stars, local_potden,local_atoms)
         CALL crank_gvecs(fi,fmpi,fi%sym,fi%cell,fi%input,sphhar, fi%vacuum , fi%noco ,local_starsq, local_potdenq,local_atoms, qvec=tmp_qvec ,iDir=iDir,iDtype=iDtype)
      END IF

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
         
        IF (l_vext) THEN
         CALL local_vTot%copyPotDen(local_potden)
         CALL local_vTot%resetPotDen()

         CALL local_dfptvTot%copyPotDen(local_potdenq)
         CALL local_dfptvTot%resetPotDen()

         CALL local_dfptvTotimag%copyPotDen(local_potdenq)
         CALL local_dfptvTotimag%resetPotDen()

         CALL local_den%copyPotDen(local_potden)
         CALL local_den%resetPotDen()
         call cast_onto_larger_grid(local_den,den,stars,input)
         
         CALL local_dfptdenreal%copyPotDen(local_potdenq)
         CALL local_dfptdenreal%resetPotDen()
         call cast_onto_larger_grid(local_dfptdenreal,dfptdenreal,starsq,input)

         CALL local_dfptdenimag%copyPotDen(local_potdenq)
         CALL local_dfptdenimag%resetPotDen()
         call cast_onto_larger_grid(local_dfptdenimag,dfptdenimag,starsq,input)

         CALL local_vCoul%copyPotDen(local_potdenq)
         CALL local_vCoul%resetPotDen()
         call cast_onto_larger_grid(local_vCoul,vCoul,starsq,input)

         CALL local_vCoul%copyPotDen(local_potdenq)
         CALL local_vCoul%resetPotDen()

         CALL local_workdenReal%copyPotDen(local_potdenq)
         CALL local_workdenReal%resetPotDen()
         call cast_onto_larger_grid(local_workdenReal,workdenReal,starsq,input)

         CALL local_workdenImag%copyPotDen(local_potdenq)
         CALL local_workdenImag%resetPotDen()
         call cast_onto_larger_grid(local_workdenImag,workdenImag,starsq,input)

         CALL local_workden%copyPotDen(local_potden)
         CALL local_workden%resetPotDen()
         call cast_onto_larger_grid(local_workden,workden,stars,input)

         CALL local_dfptvCoulimag%copyPotDen(local_potdenq)
         CALL local_dfptvCoulimag%resetPotDen()
         call cast_onto_larger_grid(local_dfptvCoulimag,dfptvCoulimag,starsq,input)

         call local_potdenq%reset_dfpt()
         call local_potden%reset_dfpt()
        END IF 

        ! NOTE: The normal stars are also passed as an optional argument, because
        !       they are needed for surface-qlm.
        sigma_loc = sigma_disc
        IF (l_vext) THEN 
        CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,juphon,local_starsq,cell,sphhar,local_atoms,.TRUE.,local_workdenReal,local_vCoul,sigma_loc,&
                        & dfptdenimag=local_workdenImag,dfptvCoulimag=local_dfptvCoulimag,dfptden0=local_workden,stars2=local_stars,iDtype=iDtype,iDir=iDir)
        ELSE 
         CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,juphon,starsq,cell,sphhar,atoms,.TRUE.,workdenReal,vCoul,sigma_loc,&
         & dfptdenimag=workdenImag,dfptvCoulimag=dfptvCoulimag,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)
        END IF 
      ! b)
      IF (l_vext) THEN
         CALL local_vCoul%copy_both_spin(local_dfptvTot)
         CALL local_dfptvCoulimag%copy_both_spin(local_dfptvTotimag)
      ELSE
         CALL vCoul%copy_both_spin(dfptvTot)
         CALL dfptvCoulimag%copy_both_spin(dfptvTotimag)
      END IF 


      ! c)
      IF (l_vext) THEN
         CALL local_denRot%init(local_stars,atoms,sphhar,vacuum,noco,input%jspins,0)
         local_denRot=local_den
         CALL local_den1Rot%init(local_starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
         CALL local_den1imRot%init(local_starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
         local_den1Rot=local_dfptdenreal
         local_den1imRot=local_dfptdenimag
      ELSE
         CALL denRot%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
         denRot=den
         CALL den1Rot%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
         CALL den1imRot%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0)
         den1Rot=dfptdenreal
         den1imRot=dfptdenimag
      END IF 
      IF (noco%l_noco) THEN
         CALL rotate_int_den_to_local(sym,stars,atoms,sphhar,vacuum,cell,input,noco ,denRot)
         IF (any(noco%l_unrestrictMT)) CALL rotate_mt_den_to_local(atoms,sphhar,sym,noco,denrot)
         !Functions that construct the spin-dependent perturbed densities
         !from the perturbed charge and (vectorial) magnetization density/
         !perturbed density matrix. Also saves the perturbed angles.
         ! TODO: Work on the internal spin logic and add vacuum as well. DFPT_NOCO
          CALL get_int_local_perturbation(sym, stars, atoms, sphhar, input, denRot, den1Rot, den1imRot, starsq)
          IF (any(noco%l_unrestrictMT)) CALL get_mt_local_perturbation(atoms,sphhar,sym,noco,denRot,den1Rot,den1imRot)
      END IF

      ! Skip vxc if we want only vC/vExt
         IF (l_xc) CALL vgen_xcpot(hybdat,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                        cell,fmpi,noco,den,denRot,EnergyDen,dfptvTot,vx,vxc,exc, &
                        & den1Rot=den1Rot, den1Rotimag=den1imRot, dfptvTotimag=dfptvTotimag,starsq=starsq)

      IF (iDtype/=0.AND.ANY(killcont/=0)) THEN
         ! d)
         ! NOTE: This is so different from the base case, that we build a new subroutine.
         IF (.NOT. l_vext) THEN 
            CALL dfpt_vgen_finalize(fmpi,atoms,stars,sym,juphon,noco,nococonv,input,sphhar,vTot,dfptvTot,dfptvTotimag,denRot,den1Rot,den1imRot,starsq,killcont)
         ELSE
            CALL dfpt_vgen_finalize(fmpi,local_atoms,local_stars,sym,juphon,noco,nococonv,input,sphhar,local_vTot,local_dfptvTot,local_dfptvTotimag,&
                                  & local_denRot,local_den1Rot,local_den1imRot,local_starsq,killcont)
         END IF 
         !DEALLOCATE(vcoul%pw_w)
      ELSE
         ! TODO: Write here something for the gradient. It does not need pw(_w)-stuff.
      END IF




      !call local_stars%reset_stars()
      !call local_starsq%reset_stars()
      !call cast_smaller_grid(vCoul,local_vCoul,starsq,input)
      !call cast_smaller_grid(workdenReal,local_workdenReal,starsq,input)
      !call cast_smaller_grid(workdenImag,local_workdenImag,starsq,input)
      !call cast_smaller_grid(workden,local_workden,stars,input)
      !call cast_smaller_grid(dfptvCoulimag,local_dfptvCoulimag,starsq,input)

      !call local_workdenReal%reset_dfpt()
      !call local_vCoul%reset_dfpt()
      !call local_workdenImag%reset_dfpt()
      !call local_dfptvCoulimag%reset_dfpt()
      !call local_workden%reset_dfpt()


      IF (l_vext) THEN
         call local_den%reset_dfpt()
         call local_den1Rot%reset_dfpt()
         call local_denRot%reset_dfpt()
         call local_den1imRot%reset_dfpt()
         call local_workdenReal%reset_dfpt()
         call local_vCoul%reset_dfpt()
         call local_workdenImag%reset_dfpt()
         call local_dfptvCoulimag%reset_dfpt()
         call local_workden%reset_dfpt()


         CALL local_dfptvTot%distribute(fmpi%mpi_comm)
         CALL local_dfptvTotimag%distribute(fmpi%mpi_comm)
      ELSE
         CALL dfptvTot%distribute(fmpi%mpi_comm)
         CALL dfptvTotimag%distribute(fmpi%mpi_comm)
      END IF 
      

      IF (l_vext) THEN 
         DO ispin = 1 , input%jspins
            write(oUnit,*) "I am here calling with q point" , tmp_qvec
            CALL checkDOPALL(input, sphhar, local_starsq ,atoms, sym, vacuum, cell,local_dfptvTot,ispin, local_dfptvTotimag   )
         END DO 
      END IF 

  END SUBROUTINE dfpt_vgen

END MODULE m_dfpt_vgen
