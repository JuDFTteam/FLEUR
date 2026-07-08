!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_vgen
    use m_juDFT
    IMPLICIT NONE
    !*****************************************************************
    ! Potential setup for the gfleur code, per layer, on the modern
    ! t_gflayer containers:
    ! gf_vgen      : driver called from gfleur_execute
    ! gf_vgen_make : (re)generate the potential from the charge density
    !                (Coulomb via fleur_vmts, XC via fleur_xcpot)
    ! gf_vgen_read : read gf_pot.hdf and warp the interstitial potential
    !                with the GF step function
    !                         Daniel Wortmann, Juelich 2002
    !*****************************************************************
    PRIVATE
    PUBLIC::gf_vgen
CONTAINS
    SUBROUTINE gf_vgen(gmpi,embinp,gfmix,layers,ld,vchk,jspins,pot_aux)
        !*****************************************************************
        ! Returns the (warped) interstitial potential and the MT potential
        ! of every layer in ld(:)%vTot for the following GF calculation.
        !*****************************************************************
        USE m_constants, ONLY: oUnit
        USE m_gf_plot
        USE m_gf_potdis
        USE m_gf_mix
        USE m_gf_types
        IMPLICIT NONE
        TYPE(t_gfmpi),INTENT(INOUT)  :: gmpi
        TYPE(t_embinp),INTENT(IN)    :: embinp
        TYPE(t_gfmix),INTENT(IN)     :: gfmix
        TYPE(t_layers),INTENT(IN)    :: layers
        TYPE(t_gflayer),INTENT(INOUT):: ld(:)
        LOGICAL,INTENT(IN)           :: vchk
        INTEGER,INTENT(IN)           :: jspins
        COMPLEX,INTENT(OUT)          :: pot_aux(:,:)

        LOGICAL :: l_exist
        INTEGER :: layer,layer_loop
        REAL    :: c_layer

        INQUIRE(FILE='gf_cdn.hdf',EXIST=l_exist)
        l_exist=l_exist.AND.embinp%l_charge
        pot_aux=0.0
        DO layer_loop = 1,gmpi%kl_LayerperPE
            layer = gmpi%kl_layers(layer_loop)
            IF (gfmix%iter>0.AND.l_exist.AND.gmpi%k_kpts(1) == 1) THEN
                CALL gf_vgen_make(embinp,layer,gmpi,ld(layer),vchk,jspins)
                c_layer = MERGE(layers%c(layer),ld(layer)%fi%cell%amat(3,3), &
     &                          layers%c(layer)>0.0)
                CALL gf_potdis(jspins,ld(layer)%fi%atoms,ld(layer)%stars,     &
     &               ld(layer)%sphhar,gmpi%fmpi,ld(layer)%fi%sym,             &
     &               ld(layer)%fi%cell,c_layer,.TRUE.,layer)
                IF (ld(layer)%fi%noco%l_noco) CALL juDFT_warn(               &
     &               "In NoCo-calculation spin rotation of potential is missing")
            ENDIF
        ENDDO
        IF (gfmix%iter>0.AND.l_exist.AND.gmpi%k_kpts(1) == 1.AND.             &
     &      embinp%l_potmix) THEN
            CALL gf_mix(gmpi,embinp,gfmix,layers,ld,jspins)
        ENDIF
        !     Now read the potential from the file!
        IF (gmpi%pe0) WRITE(oUnit,*) "Reading potential from gf_pot.hdf"
        DO layer_loop = 1,gmpi%kl_LayerperPE
            layer = gmpi%kl_layers(layer_loop)
            CALL gf_vgen_READ(layer,embinp,gmpi,ld(layer),jspins,pot_aux)
            CALL gf_plot(layer,ld(layer)%stars,ld(layer)%fi%cell,             &
     &           ld(layer)%fi%atoms,ld(layer)%fi%sym,jspins,ld(layer)%vTot%pw, &
     &           GF_PLOT_TOTPOT,ld(layer)%sphhar,ld(layer)%vTot%mt)
        ENDDO
        IF (gmpi%pe0) WRITE(oUnit,*) "Reading potential from gf_pot.hdf ... done"
    END SUBROUTINE gf_vgen


    SUBROUTINE gf_vgen_read(layer,gfinp,gmpi,ld,jspins,pot_aux)
        !*****************************************************************
        ! Reads the potential from gf_pot.hdf into ld%vTot and applies the
        ! GF step function to the interstitial part.
        !*****************************************************************
        USE m_constants,ONLY:pi_const
        USE m_gf_types
        USE m_gf_iodop
        USE m_gf_stepsanaly
        USE m_gf_checkdop
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: layer
        INTEGER, INTENT (IN)        :: jspins
        TYPE(t_embinp),INTENT(IN)   :: gfinp
        TYPE(t_gfmpi),INTENT(IN)    :: gmpi
        TYPE(t_gflayer),INTENT(INOUT) :: ld
        COMPLEX,INTENT(OUT)         :: pot_aux(:,:)

        REAL,ALLOCATABLE    :: vr(:,:,:,:),vr_c(:,:,:,:)
        COMPLEX,ALLOCATABLE :: vpw(:,:),vpw_w(:,:)
        INTEGER             :: js,n,jspin,ispins

        ispins=jspins
        IF (ld%fi%noco%l_noco) ispins=3
        ALLOCATE(vr(SIZE(ld%vTot%mt,1),0:SIZE(ld%vTot%mt,2)-1,               &
     &              SIZE(ld%vTot%mt,3),jspins))
        ALLOCATE(vpw(ld%stars%ng3,ispins))
        ALLOCATE(vpw_w(ld%stars%ng3,ispins))
        vr=0.0;vpw=0.0

        CALL gf_loddop(GF_POTFILE,layer,jspins,ld%fi%atoms,ld%stars,        &
     &       ld%sphhar,vr,vpw,ld%fi%noco,.FALSE.)

#ifndef CPP_MPI
        ALLOCATE(vr_c(SIZE(vr,1),0:SIZE(vr,2)-1,SIZE(vr,3),SIZE(vr,4)))
        vr_c=vr
        DO jspin=1,jspins
            DO n = 1,ld%fi%atoms%ntype
                vr_c(:ld%fi%atoms%jri(n),0,n,jspin) = SQRT(4*pi_const)       &
     &              *vr(:ld%fi%atoms%jri(n),0,n,jspin)/ld%fi%atoms%rmsh(:,n)
            ENDDO
        ENDDO
        CALL gf_checkdop(ld%fi%atoms,ld%fi%cell,ld%stars,ld%sphhar,          &
     &       ld%fi%sym,.FALSE.,vpw,vr_c)
        DEALLOCATE(vr_c)
#endif

        DO js=1,ispins
            IF (gfinp%npw==0) THEN
                pot_aux(:,js) = CMPLX(0.0,0.0)
            ELSE IF (gfinp%npw<0) THEN
                pot_aux(:,js) = vpw(1,js)
            ELSE
                CALL juDFT_error("aux-potential must be constant",           &
     &               calledby="gf_vgen.F90")
            ENDIF
            CALL gf_initstepsanaly(ld%stars,gfinp%napw(layer))
            CALL gf_gspaceconvolve(layer,ld%stars,REAL(pot_aux(2,js)),        &
     &           vpw(:,js),vpw_w(:,js))
        ENDDO

        ld%vTot%mt = vr
        ld%vTot%pw(:,1:jspins) = vpw_w(:,1:jspins)
    END SUBROUTINE


    SUBROUTINE gf_vgen_make(gfinp,layer,gmpi,ld,vchk,jspins)
        !*****************************************************************
        ! (Re)generate the potential from the charge density. In this
        ! historically frozen build the newly generated potential is only
        ! checked and plotted, not written back - gf_vgen_read reloads the
        ! potential from gf_pot.hdf.
        !*****************************************************************
        USE m_constants, ONLY : oUnit, pi_const
        USE m_gf_types
        USE m_gf_iodop
        USE m_gf_checkdop
        USE m_gf_plot
        USE m_fleur_pot
        USE m_gf_cdntot
        IMPLICIT NONE
        INTEGER, INTENT (IN)         :: layer
        INTEGER, INTENT (IN)         :: jspins
        LOGICAL, INTENT (IN)         :: vchk
        TYPE(t_embinp),INTENT(IN)    :: gfinp
        TYPE(t_gfmpi),INTENT(IN)     :: gmpi
        TYPE(t_gflayer),INTENT(INOUT):: ld

        REAL,ALLOCATABLE    :: vr(:,:,:,:),rho(:,:,:,:)
        COMPLEX,ALLOCATABLE :: vpw(:,:),qpw(:,:)
        REAL                :: q_el,q_nuc
        INTEGER             :: js,n,nlhd,jmtd

        nlhd=SIZE(ld%sphhar%clnu,2)-1
        jmtd=SIZE(ld%fi%atoms%rmsh,1)
        ALLOCATE(vr(jmtd,0:nlhd,ld%fi%atoms%ntype,jspins))
        ALLOCATE(vpw(ld%stars%ng3,jspins))
        ALLOCATE(rho(jmtd,0:nlhd,ld%fi%atoms%ntype,jspins))
        ALLOCATE(qpw(ld%stars%ng3,jspins))

        IF (gmpi%pe0) WRITE (oUnit,FMT = 8000)
8000    FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)

        !<-- load the charge density and check neutrality
        CALL gf_loddop(GF_CDNFILE,layer,jspins,ld%fi%atoms,ld%stars,         &
     &       ld%sphhar,rho,qpw)
        CALL gf_cdntot(layer,gmpi%fmpi,jspins,ld%stars,ld%fi%cell,           &
     &       ld%fi%atoms,rho,qpw,q_el,q_nuc)
        ! spin sum for the Coulomb potential
        IF (jspins==2) THEN
            rho(:,:,:,1) = rho(:,:,:,1) + rho(:,:,:,jspins)
            qpw(:,1) = qpw(:,1) + qpw(:,jspins)
        END IF
        IF (vchk) CALL gf_checkdop(ld%fi%atoms,ld%fi%cell,ld%stars,          &
     &       ld%sphhar,ld%fi%sym,.TRUE.,qpw(:,1:1),rho)

        !<-- interstitial Coulomb potential (filled by gf_makeintcoulpot)
        vpw(:,1) = ld%vTot%pw(:,1)
        IF (SIZE(vpw,2)>1) vpw(:,2:)=0.0
        !<-- Coulomb potential in the MT spheres
        CALL fleur_vmts(gmpi%fmpi,ld,vpw,rho,vr)
        CALL gf_plot(layer,ld%stars,ld%fi%cell,ld%fi%atoms,ld%fi%sym,1,vpw,   &
     &       GF_PLOT_HARTREE,ld%sphhar,vr)
        IF (vchk) CALL gf_checkdop(ld%fi%atoms,ld%fi%cell,ld%stars,          &
     &       ld%sphhar,ld%fi%sym,.FALSE.,vpw(:,1:1),vr)
        !<-- XC potential
        CALL fleur_xcpot(ld,gmpi%fmpi,layer,jspins,vr,vpw)
        IF (vchk) CALL gf_checkdop(ld%fi%atoms,ld%fi%cell,ld%stars,          &
     &       ld%sphhar,ld%fi%sym,.FALSE.,vpw(:,1:1),vr)

        ! store v(l=0) component as r*v(l=0)/sqrt(4pi)
        DO js = 1,jspins
            DO  n = 1,ld%fi%atoms%ntype
                vr(:,0,n,js) = ld%fi%atoms%rmsh(:,n)*vr(:,0,n,js)/SQRT(4*pi_const)
            ENDDO
        ENDDO

        IF (gmpi%pe0) WRITE(*,*) "gf_vgen:no new pot"
    END SUBROUTINE gf_vgen_make
END
