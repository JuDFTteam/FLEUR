!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_gfleur_compose
    use m_juDFT
    IMPLICIT NONE
CONTAINS
    SUBROUTINE gf_gfleur_compose(layer,noco,gfinp,layers,nk,          &
    jspin,sym,cell,mpi,lapw,                  &
    kpts,pot_aux,charge,                      &
    atoms,stars,sphhar,vr,enpara)
        !****************************************************
        !     num_layers>1 => calculate the properties of the
        !       composed system using the information on the
        !       subsystems.
        !     num_layers=1 => calculate the properties of the
        !       system
        !     Frank Freimuth, November 2007
        !****************************************************
        USE m_gf_energies
        USE m_gf_tmat
        USE m_gf_types
        USE m_gf_gmatfromtmat
        USE m_gf_read_tmat
#include "juDFT_env.h"
        USE m_gf_gfcn
        USE m_gf_landauer
        USE m_gf_propaembcurr
        USE m_gf_propaembcurr2
        USE m_gf_propaembcurr3
        USE m_gf_current4
        USE m_gf_get_spectrum
        USE m_gf_spectrum
        USE m_gf_cdnskval
        USE m_gf_ab_coef
        USE m_gf_nognofinalize
        USE m_gf_dos
        IMPLICIT NONE
        INTEGER,INTENT(IN)::layer
        TYPE(t_noco),INTENT(IN)::noco
        TYPE(t_gfinp),INTENT(IN)::gfinp
        TYPE(t_layers),INTENT(IN)::layers
        TYPE(t_kpts),INTENT(IN) :: kpts
        INTEGER,INTENT(IN)::nk
        INTEGER,INTENT(IN)::jspin
        TYPE(t_sym),INTENT(IN)::sym
        TYPE(t_cell),INTENT(IN)::cell
        TYPE(t_mpi),INTENT(IN)::mpi
        TYPE(t_lapw),INTENT(INOUT)::lapw
        COMPLEX,INTENT(IN)::pot_aux(2,2)
        TYPE(t_potden),INTENT(INOUT)::charge
        TYPE(t_atoms),INTENT(IN)    :: atoms
        TYPE(t_sphhar),INTENT(IN)   :: sphhar
        TYPE(t_stars),INTENT(IN)    :: stars
        REAL, INTENT(INOUT)         :: vr(:,0:,:,:)
        TYPE(t_enpara),INTENT(IN)   :: enpara

        COMPLEX,ALLOCATABLE::g(:,:)
        COMPLEX,ALLOCATABLE::g_sum(:,:)
        COMPLEX,ALLOCATABLE::gij(:,:,:)
                                         !not used
        COMPLEX,ALLOCATABLE::dgij(:,:,:)
        COMPLEX,ALLOCATABLE::tmat(:,:)
        COMPLEX,ALLOCATABLE::embpot_left(:,:)
        COMPLEX,ALLOCATABLE::embpot_right(:,:)
        INTEGER en,ab_spin
        LOGICAL first_g
                                                                        
        first_g=.TRUE.

        ab_spin=1
        if (noco%l_noco) ab_spin=2

        lapw%nmat_sphere=lapw%nmat+lapw%nv_tot_sphere-lapw%nv_tot

        ALLOCATE( gij(2*lapw%nv2_tot,2*lapw%nv2_tot,1) )
        IF(gfinp%l_charge.or.gfinp%l_dos)THEN
            ALLOCATE( g(lapw%nmat_sphere,lapw%nmat_sphere) )
        ENDIF
        IF(gfinp%l_charge.or.gfinp%l_dos)THEN
            ALLOCATE( g_sum(lapw%nmat_sphere,lapw%nmat_sphere) )
            ALLOCATE( embpot_left(lapw%nv2_tot,lapw%nv2_tot) )
            ALLOCATE( embpot_right(lapw%nv2_tot,lapw%nv2_tot) )
            IF(gfinp%l_spectral)THEN
                CPP_juDFT_timestart("gf_get_spectrum_sph")
                CALL gf_spectrum_clean()
                CALL gf_get_spectrum(layer,jspin,gfinp,cell,lapw,.FALSE.,      &
                gfinp%l_fullgreen,gfinp%l_nogno,gfinp%l_nohelpregion,     &
                .false.,nk,.true.)
                !CALL gf_read_spectrum(layer,nk)
                !lapw%nmat=lapw%nmat_sphere
                !lapw%nv=lapw%nv_sphere
                !lapw%nv_tot=lapw%nv_tot_sphere
                CPP_juDFT_timestop("gf_get_spectrum_sph")

            ENDIF 
            g_sum=cmplx(0.0,0.0)
        ENDIF
        IF (gfinp%l_dos.or.gfinp%l_charge) THEN
            CPP_juDFT_timestart("gf_ab_coef")
            CALL gf_ab_coef_calc(noco%l_noco,jspin,kpts%bk(:,nk),sym                        &
            ,enpara%el0,vr(:,0,:,:),atoms,cell         &
            ,lapw)
            CPP_juDFT_timestop("gf_ab_coef")
        ENDIF
        DO en=1,gf_noen()
            IF(gfinp%l_gmat)THEN
                IF(layers%num_layers==1 .OR. gfinp%l_charge.or.gfinp%l_dos)THEN
                    CALL gf_GFCN(.FALSE.,layer,en,nk,jspin,cell,              &
                    lapw,gfinp                                 &
                    ,noco%l_noco,sym%invs,gfinp%l_charge.or.gfinp%l_dos,g,gij)
                    IF(gfinp%l_charge)THEN
                        g_sum = g_sum+gf_weightz(en)*g
                    ENDIF

                    IF (gfinp%l_dos) THEN
                        IF (gfinp%l_nogno) THEN
                        CPP_juDFT_timestart("gf_nognofinalize_single")
                        CALL gf_nognofinalize(g,layer,lapw,en)
                        CPP_juDFT_timestop("gf_nognofinalize_single")
                        endif
                        CPP_juDFT_timestart("gf_dos_mt")
                        CALL gf_dos_mt(layer,atoms,gfinp,                     &
                        g,sym,en,jspin,kpts%weight(nk),nk,noco%l_noco,lapw)
                        IF (gfinp%l_intdos)  CALL gf_dos_INT(gfinp,layer      &
                        ,mpi,lapw,stars,jspin,cell%omtil,G,kpts%weight(nk&
                        ),en,nk,noco%l_noco)
                        CPP_juDFT_timestop("gf_dos_mt")
                    ENDIF
                ENDIF
                  !l_gmat
            ENDIF
            IF (gfinp%curr /= 0)  CPP_juDFT_timestart("current(misc)")
            IF (gfinp%curr ==-1) THEN
                                               !calculate gij from t-matrice
                IF(layers%num_layers>1)THEN
                    ALLOCATE( tmat(2*lapw%nv2_tot,2*lapw%nv2_tot) )
                    CALL gf_read_tmat(                                       &
                    layers,lapw%nv2_tot,en,nk,jspin,        &
                    lapw,                                   &
                    tmat)
                    CALL gf_gmatfromtmat(                                    &
                    lapw%nv2_tot,en,nk,jspin,lapw,                        &
                    .TRUE.,cell,                                          &
                    gij(:,:,1),                                           &
                    tmat)
                    DEALLOCATE( tmat )
                ENDIF
                CALL gf_landauer(lapw%nv2_tot,en,nk,jspin,cell              &
                ,sym,lapw,kpts%bk,gij(:,:,1),mpi)
            ENDIF
            IF (gfinp%curr ==1)THEN
                CALL gf_propaembcurr(layers,lapw%nv2_tot,en,nk,jspin,       &
                lapw,kpts%bk,sym,cell,mpi)
            ENDIF
            IF (gfinp%curr ==2)THEN
                CALL gf_propaembcurr2(layers,lapw%nv2_tot,en,nk,jspin,      &
                lapw,kpts%bk,sym,cell,mpi)
            ENDIF
            IF (gfinp%curr ==3)THEN
                CALL gf_propaembcurr3(layers,                               &
                lapw%nv2_tot,en,nk,jspin,              &
                lapw,kpts%bk,sym,cell,gfinp,mpi)
            ENDIF
            IF (gfinp%curr==4) THEN
                CALL gf_current4 (layers,                                   &
                noco%l_noco,en,nk,jspin,                  &
                lapw,kpts%bk,sym,cell,gfinp,mpi)
            ENDIF
            IF (gfinp%curr /= 0)  CPP_juDFT_timestop("current(misc)")
            IF ( gfinp%l_tmat ) THEN
                CPP_juDFT_timestart("gf_tmat")
                IF (pot_aux(1,jspin) /= pot_aux(2,jspin)) CALL               &
                juDFT_error                                                  &
                ("Cannot treat different aux_potentials")
                                                                        
                CALL gf_tmat(layer,                                          &
                (layers%num_layers>1),layers,en,nk,jspin,sym,cell,      &
                mpi,lapw,kpts%bk,gfinp,pot_aux(1,jspin),gij,dgij)
                CPP_juDFT_timestop("gf_tmat")
            ENDIF

        ENDDO
        IF(allocated(g))DEALLOCATE(g)
        IF(gfinp%l_charge)THEN
            IF (gfinp%l_nogno) THEN
                CPP_juDFT_timestart("gf_nognofinalize")
                CALL gf_nognofinalize(g_sum,layer,lapw)
                CPP_juDFT_timestop("gf_nognofinalize")
            ENDIF
            CPP_juDFT_timestart("gf_cdnskval") 
            CALL gf_cdnskval(noco%l_noco,                                           &
            lapw,jspin,nk,                                        &
            MAXVAL(sym%ntypsy),enpara,vr(:,0,:,:),              &
            G_sum,atoms,stars,sphhar,                             &
            cell,sym,kpts,                                &
            charge%pw(:,:),                                  &
            charge%mt(:,0:,:,:),                             &
            charge%qmtl_new(:,:))

            CPP_juDFT_timestop("gf_cdnskval") 
            DEALLOCATE(g_sum) 
        ENDIF
        call gf_ab_coef_delete()
    END SUBROUTINE gf_gfleur_compose
END
