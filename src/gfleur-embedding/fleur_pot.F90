!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_pot
      use m_juDFT
      IMPLICIT NONE
!-------------------------------------------------------------
!Module provides an interface to the potential subroutines of FLEUR
!as used in the GF-code.
!
!In the port to current FLEUR the old visxc/visxcg/vmtxc/vmtxcg calls
!(integer icorr XC ids) are replaced by the standard vis_xc/vmt_xc
!working on t_potden with the CLASS(t_xcpot) from the layer's inp.xml;
!the wrappers keep the array interface of the GF potential generator
!and adapt to t_potden internally. The unwarped interstitial potential
!is returned - the GF code applies its own step functions afterwards.
!------------------------------------------------------------
      CONTAINS
      !<-- S: fleur_xcpot(ld,fmpi,layer,jspins,vr,vpw)

      SUBROUTINE fleur_xcpot(ld,fmpi,layer,jspins,vr,vpw)
!-----------------------------------------------
!Adds the XC potential of the layer's density (loaded from the
!GF density file) to the given Coulomb potential arrays.
!           (last modified: 05-03-23) D. Wortmann
!-----------------------------------------------
      USE m_gf_types
      USE m_gf_iodop
      USE m_vis_xc
      USE m_vmt_xc
      USE m_metagga
      USE m_constants, ONLY: POTDEN_TYPE_POTTOT, POTDEN_TYPE_DEN
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_gflayer),INTENT(IN)   :: ld
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      INTEGER,INTENT(IN)           :: layer
      INTEGER,INTENT(IN)           :: jspins
      REAL   ,INTENT(INOUT)        :: vr(:,0:,:,:)
      COMPLEX,INTENT(INOUT)        :: vpw(:,:)
      !>
      !<-- Locals
      TYPE(t_potden) :: den,vTot,vx,vxc,exc,EnergyDen
      TYPE(t_kinED)  :: kinED
      !>

      IF (ld%fi%noco%l_noco) CALL juDFT_error(                           &
     &     "noco not yet supported in the gfleur port",                  &
     &     calledby="fleur_pot")
      IF (ld%xcpot%exc_is_metagga()) CALL juDFT_error(                   &
     &     "metagga not supported in gfleur",calledby="fleur_pot")

      !<-- load charge density into a t_potden
      CALL den%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,         &
     &             ld%fi%noco,jspins,POTDEN_TYPE_DEN)
      CALL gf_loddop(GF_CDNFILE,layer,jspins,ld%fi%atoms,ld%stars,       &
     &        ld%sphhar,den%mt,den%pw)
      !>

      !<-- Coulomb potential is equal for both spins
      vr(:,:,:,jspins) = vr(:,:,:,1)
      vpw(:,jspins) = vpw(:,1)
      !>

      !<-- adapt the arrays to t_potden objects
      CALL vTot%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,        &
     &             ld%fi%noco,jspins,POTDEN_TYPE_POTTOT)
      CALL vx%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,          &
     &             ld%fi%noco,jspins,POTDEN_TYPE_POTTOT)
      CALL vxc%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,         &
     &             ld%fi%noco,jspins,POTDEN_TYPE_POTTOT)
      CALL exc%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,         &
     &             ld%fi%noco,1,POTDEN_TYPE_POTTOT)
      CALL EnergyDen%init(ld%stars,ld%fi%atoms,ld%sphhar,ld%fi%vacuum,   &
     &             ld%fi%noco,jspins,POTDEN_TYPE_POTTOT)
      IF (.NOT.ALLOCATED(vTot%pw_w)) ALLOCATE(vTot%pw_w,mold=vTot%pw)
      IF (.NOT.ALLOCATED(vx%pw_w)) ALLOCATE(vx%pw_w,mold=vx%pw)
      IF (.NOT.ALLOCATED(vxc%pw_w)) ALLOCATE(vxc%pw_w,mold=vxc%pw)
      IF (.NOT.ALLOCATED(exc%pw_w)) ALLOCATE(exc%pw_w,mold=exc%pw)
      vTot%mt = vr
      vTot%pw = vpw
      vTot%pw_w = 0.0
      !>

      !interstitial region (adds to vTot%pw; the warped pw_w is ignored,
      !the GF step functions are applied later)
      CALL vis_xc(ld%stars,ld%fi%sym,ld%fi%cell,den,ld%xcpot,            &
     &            ld%fi%input,ld%fi%noco,EnergyDen,kinED,vTot,vx,exc,vxc)

      !muffin tin spheres region (adds to vTot%mt)
      CALL vmt_xc(fmpi,ld%sphhar,ld%fi%atoms,den,ld%xcpot,ld%fi%input,   &
     &            ld%fi%sym,EnergyDen,kinED,ld%fi%noco,vTot,vx,exc,vxc)

      vr = vTot%mt
      vpw = vTot%pw

      END SUBROUTINE

      !>
      !<-- S: fleur_vmts(fmpi,ld,vpw,rho,vr)
      SUBROUTINE fleur_vmts(fmpi,ld,vpw,rho,vr)
!-----------------------------------------------
!   interface to vmts of FLEUR: MT Coulomb potential of the (spin-
!   summed) density given the interstitial potential on the boundary
!           (last modified: 05-03-23) D. Wortmann
!-----------------------------------------------
      USE m_vmts
      USE m_gf_types
      USE m_constants, ONLY: POTDEN_TYPE_POTCOUL
      IMPLICIT NONE
      !<-- Arguments
      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_gflayer),INTENT(IN)   :: ld
      COMPLEX, INTENT (IN)         :: vpw(:,:)
      REAL,    INTENT (IN)         :: rho(:,0:,:,:)
      REAL,    INTENT (OUT)        :: vr(:,0:,:,:)
      !>
      vr = 0.0
      CALL vmts(ld%fi%input,fmpi,ld%stars,ld%sphhar,ld%fi%atoms,         &
     &          ld%fi%sym,ld%fi%cell,.FALSE.,vpw(:,1),rho(:,0:,:,1),     &
     &          POTDEN_TYPE_POTCOUL,vr(:,0:,:,1),1)
      END SUBROUTINE
      !>
      END
