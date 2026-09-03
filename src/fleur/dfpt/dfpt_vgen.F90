!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vgen
   USE m_juDFT

CONTAINS

   SUBROUTINE dfpt_vgen(sternheimerJob,hybdat,field,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                   dfpt, cell,fmpi,noco,nococonv,den,vTot,&
                   &starsq,dfptvTot,l_xc,dfptden1,iDtype,iDir,killcont, sliceplot,l_vextpho)
      !--------------------------------------------------------------------------
      ! FLAPW potential perturbation generator (main routine)
      !
      ! Modification for use with DFPT:
      ! The density variables in the interstitial now live in a G-expansion
      ! shifted by q and those in the Muffin Tin are no longer real. Account for
      ! that with optional arguments: use q-shifted stars and carry the imaginary
      ! part of the MT-density along explicitely. Changes:
      ! - The interstitial part and MT real part is in dfptden1%mt
      ! - The MT imaginary part is in dfptden1%mtIm
      ! - vTot will carry the same quantities for V1 as dfptden1 for rho1
      ! - The MT imaginary part is in dfptvTot%mtIm
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
      USE m_dfpt_int_perturbation
      USE m_dfpt_mt_perturbation
      USE m_dfpt_vgen_finalize
      USE m_dfpt_vefield
      USE m_dfpt_vbfield
      USE m_checkdopall
      USE m_plot
      

      IMPLICIT NONE

      TYPE(t_sternheimerJob),INTENT(IN):: sternheimerJob
      CLASS(t_xcpot),    INTENT(IN)    :: xcpot
      TYPE(t_hybdat),    INTENT(IN)    :: hybdat
      TYPE(t_mpi),       INTENT(IN)    :: fmpi

      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_field),     INTENT(IN)    :: field
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_nococonv),  INTENT(IN)    :: nococonv
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_dfpt),    INTENT(IN)    :: dfpt
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_potden),    INTENT(IN)    :: vTot
      TYPE(t_potden),    INTENT(INOUT) :: den, dfptvTot

      ! for plotting
      TYPE(t_sliceplot), OPTIONAL,INTENT(IN)   :: sliceplot
      LOGICAL                 :: l_plot
      TYPE(t_nococonv)     :: nococonv_int
      TYPE(t_sliceplot)   :: sliceplot_int

      

      LOGICAL, INTENT(IN) :: l_xc

      TYPE(t_stars),  OPTIONAL, INTENT(IN)    :: starsq
      TYPE(t_potden), OPTIONAL, INTENT(INOUT) :: dfptden1

      INTEGER, OPTIONAL, INTENT(IN)           :: iDtype, iDir ! DFPT: Type and direction of displaced atom

      INTEGER, OPTIONAL, INTENT(IN)           :: killcont(2)
      LOGICAL, OPTIONAL, INTENT(IN)           :: l_vextpho

      TYPE(t_potden)                   :: workden, denRot, workden1, den1Rot
      COMPLEX, ALLOCATABLE             :: theta1_mt(:,:), phi1_mt(:,:), theta1_pw(:), phi1_pw(:)
      TYPE(t_potden)                   :: vCoul, vxc, exc, vx, EnergyDen
      TYPE(t_potden)                   :: dfptvefield
      TYPE(t_atoms)                    :: atomsefield 

      COMPLEX :: constantShift
      INTEGER :: ispin 
      
      vCoul = dfptvTot
      vx = vTot
      vxc = vTot
      exc = vTot
      dfptvefield = dfptvTot

      IF (fmpi%irank==0) WRITE (oUnit,FMT=8000)
      IF (fmpi%irank==0) WRITE (oUnit,FMT=8001)
8000  FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
8001  FORMAT (/,/,t10,'          (DFPT edition)              ',/)
      CALL dfptvTot%resetPotDen()
      CALL vCoul%resetPotDen()
      CALL vx%resetPotDen()
      CALL vxc%resetPotDen()
      CALL exc%resetPotDen()
      CALL dfptvefield%resetPotDen()


      ALLOCATE(vx%pw_w,mold=vTot%pw)
      vx%pw_w = 0.0
      ALLOCATE(vxc%pw_w,mold=vTot%pw)
      vxc%pw_w = 0.0
      CALL exc%init(stars, atoms, sphhar, vacuum, noco, 1, 1) !one spin only
      ALLOCATE (exc%pw_w(stars%ng3, 1)); exc%pw_w = 0.0

#ifndef CPP_OLDINTEL
      ALLOCATE(dfptvTot%pw_w,mold=dfptvTot%pw)
      IF (sternheimerJob%l_efield) THEN
         ALLOCATE(dfptvefield%pw_w,mold=dfptvefield%pw)
      end if  
#else
      ALLOCATE( dfptvTot%pw_w(size(dfptvTot%pw,1),size(dfptvTot%pw,2)))

      IF (sternheimerJob%l_efield) THEN
         ALLOCATE( dfptvefield%pw_w(size(dfptvefield%pw,1),size(dfptvefield%pw,2)))
      end if  
         
#endif

      ALLOCATE(vCoul%pw_w(SIZE(vCoul%pw,1),size(vCoul%pw,2)))
      vCoul%pw_w = CMPLX(0.0,0.0)

      CALL workDen%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
      CALL workDen1%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0,l_dfpt=.TRUE.)

      !for plotting
      l_plot = .FALSE.
      IF (PRESENT(sliceplot)) THEN 
         l_plot = .TRUE.
         nococonv_int = nococonv
         sliceplot_int = sliceplot
      END IF 

      ! a)
      ! Sum up both spins in den into workden:
      CALL den%sum_both_spin(workden)
      CALL dfptden1%sum_both_spin(workden1)
      ! NOTE: The normal stars are also passed as an optional argument, because
      !       they are needed for surface-qlm.
      IF (sternheimerJob%l_efield) THEN
         atomsefield = atoms
         atomsefield%zatom(:) = 0.0 ! find out if this is actually needed
         CALL dfpt_vefield(dfpt,starsq,atoms,sym,sphhar,cell,dfptvefield,iDir,1)
         CALL dfptvefield%copy_both_spin(dfptvTot)

         IF (l_xc) THEN
            CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,starsq,cell,sphhar,atomsefield,.TRUE.,workden1,vCoul, sternheimerJob=sternheimerJob,&
                     & dfpt=dfpt,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)
            dfptvTot%pw = dfptvTot%pw + vCoul%pw
            dfptvTot%mt = dfptvTot%mt + vCoul%mt
            dfptvTot%mtIm = dfptvTot%mtIm + vCoul%mtIm
         END IF
      else if (sternheimerJob%l_bfield) then
         atomsefield = atoms
         atomsefield%zatom(:) = 0.0 ! find out if this is actually needed, actually check that
         CALL dfpt_vbfield(input,stars,noco,atoms,dfptvTot)
         !CALL dfptvefield%copy_both_spin(dfptvTot)

         IF (l_xc) THEN
            CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,starsq,cell,sphhar,atomsefield,.TRUE.,workden1,vCoul,sternheimerJob=sternheimerJob,&
                     & dfpt=dfpt,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)
            dfptvTot%pw = dfptvTot%pw + vCoul%pw
            dfptvTot%mt = dfptvTot%mt + vCoul%mt
            dfptvTot%mtIm = dfptvTot%mtIm + vCoul%mtIm
         end if 

      ELSE !standard phonon case
         CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,starsq,cell,sphhar,atoms,.TRUE.,workden1,vCoul, sternheimerJob=sternheimerJob,&
                     & dfpt=dfpt,dfptden0=workden,stars2=stars,iDtype=iDtype,iDir=iDir)
         ! b)

         CALL vCoul%copy_both_spin(dfptvTot)
         !print*,"sum dfptvTot in dfpt_vgen pho", sum(dfptvTot%pw(:,1))
      END IF
      ! c)
      CALL denRot%init(stars,atoms,sphhar,vacuum,noco,input%jspins,0)
      denRot=den
      CALL den1Rot%init(starsq,atoms,sphhar,vacuum,noco,input%jspins,0,l_dfpt=.TRUE.)
      den1Rot=dfptden1
      IF (noco%l_noco) THEN
         CALL rotate_int_den_to_local(sym,stars,atoms,sphhar,vacuum,cell,input,noco ,denRot)
         IF (any(noco%l_unrestrictMT)) CALL rotate_mt_den_to_local(atoms,sphhar,sym,noco,denrot)
         !Functions that construct the spin-dependent perturbed densities
         !from the perturbed charge and (vectorial) magnetization density/
         !perturbed density matrix. Also saves the perturbed angles.
         ! TODO: Work on the internal spin logic and add vacuum as well. DFPT_NOCO
          CALL get_int_local_perturbation(sym, stars, atoms, sphhar, input, denRot, den1Rot, theta1_pw, phi1_pw, starsq)
          IF (any(noco%l_unrestrictMT)) CALL get_mt_local_perturbation(atoms,sphhar,sym,noco,denRot,den1Rot,theta1_mt,phi1_mt)
      END IF

      ! Skip vxc if we want only vC/vExt
      IF (l_xc) CALL vgen_xcpot(hybdat,input,xcpot,atoms,sphhar,stars,vacuum,sym,&
                     cell,fmpi,noco,den,denRot,EnergyDen,dfptvTot,vx,vxc,exc, &
                     & den1Rot=den1Rot,starsq=starsq)
      
      IF (iDtype/=0.AND.ANY(killcont/=0)) THEN
         ! d)
         ! NOTE: This is so different from the base case, that we build a new subroutine.
         CALL dfpt_vgen_finalize(sternheimerJob,fmpi,atoms,stars,sym,noco,nococonv,input,sphhar,vTot,dfptvTot,denRot,den1Rot,theta1_mt,phi1_mt,theta1_pw,phi1_pw,starsq,killcont)
         !DEALLOCATE(vcoul%pw_w)
      ELSE
         ! TODO: Write here something for the gradient. It does not need pw(_w)-stuff.
      END IF
      
      CALL dfptvTot%distribute(fmpi%mpi_comm)

  END SUBROUTINE dfpt_vgen

END MODULE m_dfpt_vgen
