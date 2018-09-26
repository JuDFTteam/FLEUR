!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen_xcpot

  USE m_juDFT

CONTAINS

  SUBROUTINE vgen_xcpot(hybrid,input,xcpot,dimension,atoms,sphhar,stars,vacuum,sym,&
                        obsolete,cell,oneD,sliceplot,mpi,noco,den,denRot,vTot,vx,results)

    !     ***********************************************************
    !     FLAPW potential generator                           *
    !     ***********************************************************
    !     calculates the density-potential integrals needed for the
    !     total energy
    !     TE_VCOUL  :   charge density-coulomb potential integral
    !     TE_VEFF:   charge density-effective potential integral
    !     TE_EXC :   charge density-ex-corr.energy density integral
    !     ***********************************************************

    USE m_types
    USE m_constants
    USE m_intnv
    USE m_vmt_xc
    USE m_vvacxc
    USE m_vvacxcg
    USE m_vis_xc
    USE m_checkdopall
    USE m_cdn_io
    USE m_convol

    IMPLICIT NONE


    CLASS(t_xcpot),    INTENT(IN)              :: xcpot
    TYPE(t_hybrid),    INTENT(IN)              :: hybrid
    TYPE(t_mpi),       INTENT(IN)              :: mpi
    TYPE(t_dimension), INTENT(IN)              :: dimension
    TYPE(t_oneD),      INTENT(IN)              :: oneD
    TYPE(t_obsolete),  INTENT(IN)              :: obsolete
    TYPE(t_sliceplot), INTENT(IN)              :: sliceplot
    TYPE(t_input),     INTENT(IN)              :: input  
    TYPE(t_vacuum),    INTENT(IN)              :: vacuum
    TYPE(t_noco),      INTENT(IN)              :: noco
    TYPE(t_sym),       INTENT(IN)              :: sym
    TYPE(t_stars),     INTENT(IN)              :: stars
    TYPE(t_cell),      INTENT(IN)              :: cell
    TYPE(t_sphhar),    INTENT(IN)              :: sphhar
    TYPE(t_atoms),     INTENT(IN)              :: atoms 
    TYPE(t_potden),    INTENT(IN)              :: den,denRot
    TYPE(t_potden),    INTENT(INOUT)           :: vTot,vx
    TYPE(t_results),   INTENT(INOUT), OPTIONAL :: results

    ! Local type instances
    TYPE(t_potden) :: workDen,exc,veff
    ! Local Scalars
    INTEGER ifftd,ifftd2,ifftxc3d,ispin
#ifdef CPP_MPI
    include 'mpif.h'
    integer:: ierr
#endif


    CALL exc%init_potden_types(stars,atoms,sphhar,vacuum,1,.false.,1) !one spin only
    ALLOCATE(exc%pw_w(stars%ng3,1));exc%pw_w=0.0
    IF (PRESENT(results)) THEN
       CALL veff%init(stars,atoms,sphhar,vacuum,input%jspins,.FALSE.,1)
#ifndef CPP_OLDINTEL
       ALLOCATE(veff%pw_w,mold=veff%pw)
#else
       ALLOCATE( veff%pw_w(size(veff%pw,1),size(veff%pw,2)) )
#endif
    ENDIF

    ! exchange correlation potential

    ! vacuum region
    IF (mpi%irank == 0) THEN
       IF (input%film) THEN
          CALL timestart("Vxc in vacuum")

          ifftd2 = 9*stars%mx1*stars%mx2
          IF (oneD%odi%d1) ifftd2 = 9*stars%mx3*oneD%odi%M

          IF (.NOT.xcpot%vxc_is_gga()) THEN  ! LDA

             IF (.NOT.oneD%odi%d1) THEN
                CALL vvacxc(ifftd2,stars,vacuum,xcpot,input,noco,Den,vTot,exc)
             ELSE
                CALL judft_error("OneD broken")
                ! CALL vvacxc(stars,oneD%M,vacuum,odi%n2d,dimension,ifftd2,&
                !             xcpot,input,odi%nq2,odi%nst2,den,noco,odi%kimax2%igf,&
                !             odl%pgf,vTot%vacxy,vTot%vacz,excxy,excz)
             END IF
          ELSE      ! GGA
             IF (oneD%odi%d1) THEN
                CALL judft_error("OneD broken")
                ! CALL vvacxcg(ifftd2,stars,vacuum,noco,oneD,&
                !              cell,xcpot,input,obsolete,workDen, ichsmrg,&
                !              vTot%vacxy,vTot%vacz,rhmn, exc%vacxy,exc%vacz)

             ELSE
                CALL vvacxcg(ifftd2,stars,vacuum,noco,oneD,cell,xcpot,input,obsolete,Den,vTot,exc)
             END IF
          END IF
          CALL timestop("Vxc in vacuum")
       END IF

       ! interstitial region
       CALL timestart("Vxc in interstitial")

       
       IF ( (.NOT. obsolete%lwb) .OR. ( .not.xcpot%vxc_is_gga() ) ) THEN
          ! no White-Bird-trick
          CALL vis_xc(stars,sym,cell,den,xcpot,input,noco,vTot,vx,exc)

       ELSE
          ! White-Bird-trick
          WRITE(6,'(a)') "W+B trick cancelled out. visxcwb uses at present common block cpgft3.",&
                         "visxcwb needs to be reprogrammed according to visxcg.f"
          CALL juDFT_error("visxcwb",calledby ="vgen")
       END IF

       CALL timestop("Vxc in interstitial")
    END IF !irank==0


    !
    !     ------------------------------------------
    !     ----> muffin tin spheres region

    IF (mpi%irank == 0) THEN
       CALL timestart ("Vxc in MT")
    END IF

    CALL vmt_xc(DIMENSION,mpi,sphhar,atoms, den,xcpot,input,sym,&
         obsolete, vTot,vx,exc)
    

    !

    ! add MT EXX potential to vr
    IF (mpi%irank == 0) THEN
       CALL timestop ("Vxc in MT")

       ! check continuity of total potential
       IF (input%vchk) CALL checkDOPAll(input,dimension,sphhar,stars,atoms,sym,vacuum,oneD,cell,vTot,1)

       ! TOTAL 
       IF (PRESENT(results)) THEN
          ! CALCULATE THE INTEGRAL OF n1*Veff1 + n2*Veff2
          ! Veff = Vcoulomb + Vxc
          IF (noco%l_noco) THEN 
             workDen = denRot
          ELSE
             workden = den
          END IF

          veff = vTot
          IF(xcpot%is_hybrid().AND.hybrid%l_subvxc) THEN
            DO ispin = 1, input%jspins
                CALL convol(stars,vx%pw_w(:,ispin),vx%pw(:,ispin),stars%ufft)
            END DO
            veff%pw   = vTot%pw   - xcpot%get_exchange_weight() * vx%pw
            veff%pw_w = vTot%pw_w - xcpot%get_exchange_weight() * vx%pw_w
            veff%mt   = vTot%mt   - xcpot%get_exchange_weight() * vx%mt
            exc%pw    = exc%pw    - xcpot%get_exchange_weight() * exc%pw
            exc%pw_w  = exc%pw_w  - xcpot%get_exchange_weight() * exc%pw_w
            exc%mt    = exc%mt    - xcpot%get_exchange_weight() * exc%mt
          END IF

          results%te_veff = 0.0
          DO ispin = 1, input%jspins
             WRITE (6,FMT=8050) ispin
8050         FORMAT (/,10x,'density-effective potential integrals for spin ',i2,/)
             CALL int_nv(ispin,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,veff,workden,results%te_veff)
          END DO

          WRITE (6,FMT=8060) results%te_veff
8060      FORMAT (/,10x,'total density-effective potential integral :', t40,f20.10)

          ! CALCULATE THE INTEGRAL OF n*exc

          ! perform spin summation of charge densities for the calculation of Exc
          CALL  workden%sum_both_spin()

          WRITE (6,FMT=8070)
8070      FORMAT (/,10x,'charge density-energy density integrals',/)

          results%te_exc = 0.0
          CALL int_nv(1,stars,vacuum,atoms,sphhar,cell,sym,input,oneD,exc,workDen,results%te_exc)
          WRITE (6,FMT=8080) results%te_exc

8080      FORMAT (/,10x,'total charge density-energy density integral :', t40,f20.10)
       END IF
    END IF ! mpi%irank == 0

  END SUBROUTINE vgen_xcpot

END MODULE m_vgen_xcpot
