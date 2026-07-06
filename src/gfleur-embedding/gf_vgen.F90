!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_vgen
    use m_juDFT
    IMPLICIT NONE
    !*****************************************************************
    ! This modules contains some subroutines for the gf-potential
    ! setup:
    ! public:
    ! gf_vgen     : General subroutine called to get the potential
    ! private:
    ! gf_vgen_read: read the potential from a file
    ! gf_vgen_make: generate potential from charge density
    !
    !                         Daniel Wortmann, Juelich 2002
    !*****************************************************************
    PRIVATE
    PUBLIC::gf_vgen
CONTAINS
    SUBROUTINE gf_vgen(layers,gfinp,                                  &
    vchk,                                                        &
    stars,sphhar,sym,cell,atoms,mpi,xcpot,mix,                   &
    jspins,noco,enpara,                                                 &
    pot_aux,potential)
        !*****************************************************************
        ! DESC: Potential setup routine based on FLEUR-vgen.
        ! This routine returns the (warped)interstitial potential and the
        ! MT-potential needed for GF-calculations.
        ! The appropriate step functions are applied on the interstitial
        ! potential.
        !                         Daniel Wortmann, Thu Mar  6 13:33:08 2002
        !*****************************************************************
        USE m_gf_plot
        USE m_gf_potdis
        USE m_gf_mix
        USE m_gf_types

        IMPLICIT NONE
        !     C
        !..   Scalar Arguments .. ..
        INTEGER, INTENT (IN)         :: jspins
        LOGICAL, INTENT (IN)         :: vchk
        TYPE(t_gfinp),INTENT(IN)     :: gfinp
        TYPE(t_stars),INTENT(IN)     :: stars(:)
        TYPE(t_sphhar),INTENT(IN)    :: sphhar(:)
        TYPE(t_sym),INTENT(IN)       :: sym
        TYPE(t_cell),INTENT(IN)      :: cell(:)
        TYPE(t_atoms),INTENT(IN)     :: atoms(:)
        TYPE(t_mpi),INTENT(IN)       :: mpi
        TYPE(t_xcpot),INTENT(INOUT)  :: xcpot(:)
        TYPE(t_mix),INTENT(IN)       :: mix
        TYPE(t_layers),INTENT(IN)    :: layers
        TYPE(t_noco),INTENT(IN)      :: noco(:)
        TYPE(t_enpara),INTENT(IN)    :: enpara(:)
        TYPE(t_potential),INTENT(INOUT) :: potential(:)
        !     ..
        !     .. Array Arguments ..
        COMPLEX,INTENT(OUT) :: pot_aux(:,:)
        !
        LOGICAL :: l_exist
        INTEGER :: layer,layer_loop
        INQUIRE(FILE='gf_cdn.hdf',EXIST=l_exist)
        l_exist=l_exist.AND.gfinp%l_charge
        pot_aux=0.0
        DO layer_loop = 1,mpi%kl_LayerperPE
            layer = mpi%kl_layers(layer_loop)
                                                                    !create
            IF (mix%iter>0.AND.l_exist.AND.mpi%k_kpts(1) == 1) THEN
                                                                     !  if (
                CALL gf_vgen_make(gfinp,layer,mpi,noco(layer),                          &
                vchk,                                                  &
                stars(layer),sphhar(layer),sym,cell(layer),atoms(layer)&
                ,xcpot(layer),potential(layer),jspins)
                CALL gf_potdis(jspins,atoms(layer),stars(layer),sphhar(layer),mpi,sym,cell(layer),.TRUE.,layer)
                if (noco(layer)%l_noco) CALL juDFT_warn("In NoCo-calculation spin rotation of potential is missing")
            ENDIF
        ENDDO
        IF (mix%iter>0.AND.l_exist.AND.mpi%k_kpts(1) == 1.AND.             &
        gfinp%l_potmix) THEN
            !call juDFT_error("NOT IMPLEMENTED, potmix=T", calledby="gf_vgen")
            CALL gf_mix(atoms,cell,sphhar,stars,gfinp,noco,mpi,sym,enpara,mix     &
            ,jspins,layers)
        ENDIF
        !     ! Now read the potential from the file!
        WRITE(6,*) "Reading potential from gf_pot.hdf"
        DO layer_loop = 1,mpi%kl_LayerperPE
            layer = mpi%kl_layers(layer_loop)
            CALL gf_vgen_READ(layer,gfinp,                                 &
            stars(layer),sphhar(layer),sym,cell(layer),atoms(layer)   &
            ,mpi,jspins,pot_aux,potential(layer)%vpw,potential(layer  &
            )%vr,noco(1))
            CALL gf_plot(layer,stars(layer),cell(layer),atoms(layer),sym   &
            ,jspins,potential(layer)%vpw,GF_PLOT_TOTPOT,sphhar(layer) &
            ,potential(layer)%vr)
        ENDDO
        WRITE(6,*) "Reading potential from gf_pot.hdf ... done"
                                                                        
                                                                        
 
    END SUBROUTINE gf_vgen
                                                                        
                                                                        
    SUBROUTINE gf_vgen_read(layer,gfinp,                              &
    stars,sphhar,sym,cell,atoms,mpi,                  &
    jspins,                                           &
    pot_aux,vpw,vr,noco)
        !*****************************************************************
        ! DESC: Potential subroutine reads potential from gf_pot.
        ! The appropriate step functions are applied on the interstitial
        ! potential.
        !                         Daniel Wortmann, Thu Mar  6 13:33:08 202
        !*****************************************************************
        !
#include"cpp_double.h"                                                  
        USE m_constants,ONLY:pimach
        USE m_gf_types
        USE m_gf_iodop
        USE m_gf_stepsanaly
        USE m_gf_checkdop
        USE m_gf_fft,ONLY:gf_convol
        IMPLICIT NONE
        !     ..
        !     .. Scalar Arguments ..
        INTEGER,INTENT(IN)        ::layer
        INTEGER, INTENT (IN)      ::jspins
        TYPE(t_gfinp),INTENT(IN)  ::gfinp
        TYPE(t_stars),INTENT(IN)  ::stars
        TYPE(t_sphhar),INTENT(IN) ::sphhar
        TYPE(t_sym),INTENT(IN)    ::sym
        TYPE(t_cell),INTENT(IN)   ::cell
        TYPE(t_atoms),INTENT(IN)  ::atoms
        TYPE(t_mpi),INTENT(IN)    :: mpi
        Type(t_noco),INTENT(IN)   :: noco
        !     ..
        !     .. Array Arguments ..
        REAL,INTENT(OUT)    :: vr(:,0:,:,:)
        COMPLEX,INTENT(OUT) :: pot_aux(:,:)
        COMPLEX,INTENT(OUT) :: vpw(:,:)
                                                                        
        !     locals read  from gf_starstep
        LOGICAL,POINTER  ::step1(:),step2(:),step3(:)
#ifdef CPP_GF_GSTEP                                                     
      REAL,POINTER     ::ufft(:) 
#endif                                                                  
        !     .. Local Scalars ..
        INTEGER          ::js
        !     .. Local Arrays ..
        !for testing of continuity
        REAL,ALLOCATABLE::vr_c(:,:,:,:)
        !     for convolution with step function
        COMPLEX,ALLOCATABLE  :: zero(:)
        COMPLEX, ALLOCATABLE :: vpw_w(:,:)
        !checking of potential
        INTEGER::n,jspin
        INTEGER :: ispins
#ifdef CPP_MPI                                                          
      INCLUDE "mpif.h" 
      CALL MPI_BARRIER(MPI_COMM_WORLD,n) 
#endif                                                                  
                                                                        
        ispins=jspins
        IF (noco%l_noco) ispins=3
        ALLOCATE ( vpw_w(stars%nq3,ispins))
                                                                        
        CALL gf_loddop(GF_POTFILE,layer,jspins,                           &
        atoms,stars,sphhar,                                       &
        vr,vpw,noco,.FALSE.)
                                                                        
#ifndef CPP_MPI                                                         
        ALLOCATE(vr_c(size(vr,1),0:size(vr,2)-1,size(vr,3),size(vr,4)))
        vr_c=vr
        DO jspin=1,jspins
            DO n = 1,atoms%ntype
                vr_c(:atoms%jri(n),0,n,jspin) = SQRT(4*pimach())            &
                *         vr(:atoms%jri(n),0,n,jspin)/atoms%rmsh(:,n)
            ENDDO
        ENDDO
        CALL gf_checkdop(atoms,cell,stars,sphhar,sym,.FALSE.,vpw,vr_c)
        DEALLOCATE(vr_c)
#endif                                                                  
                                                                        
        ALLOCATE(zero(stars%nq3))
        zero=0.0
        DO js=1,ispins
                                   !use zero potential in aux. region   
            IF (gfinp%npw==0) THEN 
                pot_aux(:,js) = CMPLX(0.0,0.0)
              ! pot_aux(2,js) = CMPLX(gfinp%bias,0.0)                   
            ELSE IF (gfinp%npw<0) THEN 
                pot_aux(:,js) = vpw(1,js)
            ELSE 
                CALL juDFT_error("aux-potential must be constant",calledby="gf_vgen.F90")
            ENDIF 
            CALL gf_initstepsanaly(stars,gfinp%napw(layer)) 
            CALL gf_gspaceconvolve(layer,stars            &
            ,REAL(pot_aux(2,js)),vpw(:,js),vpw_w(:,js))
                                                                        
                                                                        
        ENDDO
                                                                        
        DEALLOCATE(zero)
                                                                        
                                                                        
        vpw=vpw_w
                                                                        
                                                                        
                                                                        
        DEALLOCATE (vpw_w)
                                                                        
    END SUBROUTINE
                                                                        
                                                                        
    SUBROUTINE gf_vgen_make(gfinp,layer,mpi,noco,                          &
    vchk,                                             &
    stars,sphhar,sym,cell,atoms,xcpot,potential,      &
    jspins)
        !*****************************************************************
        ! DESC: Potential setup routine based on FLEUR-vgen.
        !
        !*****************************************************************
        !
#include"cpp_double.h"                                                  
        USE m_constants, ONLY : pimach
        USE m_gf_types
        USE m_gf_iodop
        USE m_gf_checkdop
        USE m_gf_plot
        USE m_fleur_pot
        USE m_gf_cdntot
        USE m_gf_intcoul
        USE m_gf_noco
        !      use m_gf_coulombpotential
        IMPLICIT NONE
        !<--Arguments
        INTEGER, INTENT (IN)         :: layer
        INTEGER, INTENT (IN)         :: jspins
        LOGICAL, INTENT (IN)         :: vchk
        TYPE(t_gfinp),INTENT(IN)     :: gfinp
        TYPE(t_stars),INTENT(IN)     :: stars
        TYPE(t_sphhar),INTENT(IN)    :: sphhar
        TYPE(t_sym),INTENT(IN)       :: sym
        TYPE(t_cell),INTENT(IN)      :: cell
        TYPE(t_atoms),INTENT(IN)     :: atoms
        TYPE(t_xcpot),INTENT(INOUT)  :: xcpot
        TYPE(t_mpi),INTENT(IN)       :: mpi
        TYPE(t_potential),INTENT(IN) :: potential
        TYPE(T_noco),INTENT(IN)      :: noco
        LOGICAL                      :: l_nohelpregion
        !>
        !<-- Locals
                                                                        
        REAL,ALLOCATABLE             :: vr(:,:,:,:)
        COMPLEX,ALLOCATABLE          :: vpw(:,:)
        INTEGER             :: js,n
        !dimensions
        INTEGER             :: nlhd,jmtd
        REAL,   ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
        COMPLEX,ALLOCATABLE :: qpw(:,:)
                                                                        
        l_nohelpregion=gfinp%l_nohelpregion
                                                                        
        !>
        nlhd=size(sphhar%clnu,2)-1
        jmtd=size(atoms%rmsh,1)
                                                                        
                                                                        
        ALLOCATE ( vr(jmtd,0:nlhd,atoms%ntype,jspins),vpw(stars%nq3,jspins&
        ),rho(jmtd,0:nlhd,atoms%ntype,jspins),rht(1,2,jspins))
        allocate(qpw(stars%nq3,jspins))

                                                                        
        !<-- Say Hello
        IF (mpi%pe0) WRITE (6,FMT = 8000)
8000    FORMAT (/,/,t10,' p o t e n t i a l   g e n e r a t o r',/)
        !      vpw(:,:) = cmplx(0.,0.)
        !>
                                                                        
        !<-- load the charge density!
                                                                        
        CALL gf_loddop(GF_CDNFILE,layer,jspins,                           &
        atoms,stars,sphhar,                                          &
        rho,qpw)
        !Check for charge neutrality
        CALL gf_cdntot(layer,mpi,jspins,stars,cell,atoms,rho,qpw)
        ! perform spin summation of charge densities
        ! for the calculation of the coulomb potentials
        IF (jspins==2) THEN
            rho(:,:,:,1) = rho(:,:,:,1) + rho(:,:,:,jspins)
            qpw(:,1) = qpw(:,1) + qpw(:,jspins)
        END IF
        !<-- check continuity of the density
        IF (vchk) THEN
            IF (mpi%pe0) WRITE(6,*) "Checking Density..."
            CALL gf_checkdop(atoms,cell,stars,sphhar,sym,.TRUE.,qpw(:,1:1) &
            ,rho)
        END IF
        !>
                                                                        
        !>
                                                                        
        !<--Do Coulomb pot in INT
                                                                        
        IF (mpi%pe0) WRITE (6,FMT=8010) layer
8010    FORMAT (/,5x                                                      &
        ,'Loading interstitial coulomb potential in layer:',i3)
        !      CALL gf_lodcoul(GF_POTFILE,layer,vpw=vpw(:,1))
        vpw(:,1) = potential%vpw(:,1)
        IF (SIZE(vpw,2)>1) vpw(:,2:)=0.0
        !>
                                                                        
        !<-- Coulomb potential in the muffin-tin spheres
        CALL fleur_vmts(mpi%self_SUBCOM,jspins,stars,atoms,sphhar,sym     &
        ,cell,vpw,rho,vr)
        !>
        !<-- plot the Hartree pot
                                                                        
        CALL gf_plot(layer,stars,cell,atoms,sym,1,vpw                     &
        ,GF_PLOT_HARTREE,sphhar,vr)


        !call gf_wrtcoul(gf_potfile,layer,stars,vr=vr)
        !>
                                                                        
        !<-- check continuity of coulomb potential
        IF (vchk) THEN
            IF (mpi%pe0) WRITE(6,*) "Checking Coulomb potential..."
            CALL gf_checkdop(atoms,cell,stars,sphhar,sym,.FALSE.,vpw(:,1:1)&
            ,vr)

        END IF
        !>
        !<-- Calculate XC-potential
        CALL fleur_xcpot(layer,jspins,atoms,stars,                        &
        sphhar,xcpot,sym,cell,vr,vpw)
        !>
        !<-- check continuity of XC-potential
        IF (vchk) THEN
            IF (mpi%pe0) WRITE(6,*) "Checking XC-potential..."
            CALL gf_checkdop(atoms,cell,stars,sphhar,sym,.FALSE.,vpw(:,1:1)&
            ,vr)
        END IF
        !>
                                                                        

                                                                        
        ! store v(l=0) component as r*v(l=0)/sqrt(4pi)
        DO js = 1,jspins
            DO  n = 1,atoms%ntype
                vr(:,0,n,js) = atoms%rmsh(:,n)*vr(:,0,n,js)/sqrt(4*pimach())
            ENDDO
        ENDDO

        !if we are doing a noco-calculation we might have to rotate the interstitial potential
        if (noco%l_noco.and..false.) then
            deallocate(qpw)
            allocate(qpw(stars%nq3,3))
            call gf_loddop(GF_CDNFILE,layer,jspins,atoms,stars,sphhar,rho,qpw)
            call gf_noco_rotate(stars,qpw,vpw)
        endif


        DEALLOCATE (rho,rht)
        DEALLOCATE (qpw)
        print *, "gf_vgen:no new pot"
        return
        !     MOVE potential to old-potential
        CALL gf_renamepot(gf_potfile,mpi%iodop_subcom,layer)
        !     Write this potential to gf_pot_new!
        CALL gf_wrtdop(GF_POTFILE,layer,jspins,                           &
        gfinp,atoms,stars,sphhar,vr,vpw,.FALSE.,mpi%iodop_subcom)
        RETURN
    END SUBROUTINE gf_vgen_make
                                                                        
                                                                        
END
