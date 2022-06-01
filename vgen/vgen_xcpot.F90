!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen_xcpot

   USE m_juDFT
#ifdef CPP_MPI
   use mpi
#endif

CONTAINS

   SUBROUTINE vgen_xcpot(hybdat, input, xcpot,  atoms, sphhar, stars, vacuum, sym, &
                          cell,   sliceplot, fmpi, noco, den, denRot, EnergyDen, vTot, vx, vxc, exc, results, &
                          den1Rot, den1Rotimag, dfptvTotimag, starsq)

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
      USE m_vvac_xc
      USE m_vis_xc
      USE m_checkdopall
      USE m_cdn_io
      USE m_convol
      USE m_cdntot
      USE m_intgr
      USE m_metagga
      USE m_dfpt_vmt_xc
      USE m_dfpt_vis_xc

      IMPLICIT NONE

      CLASS(t_xcpot), INTENT(IN)           :: xcpot
      TYPE(t_hybdat), INTENT(IN)              :: hybdat
      TYPE(t_mpi), INTENT(IN)              :: fmpi


      TYPE(t_sliceplot), INTENT(IN)              :: sliceplot
      TYPE(t_input), INTENT(IN)              :: input
      TYPE(t_vacuum), INTENT(IN)              :: vacuum
      TYPE(t_noco), INTENT(IN)              :: noco
      TYPE(t_sym), INTENT(IN)              :: sym
      TYPE(t_stars), INTENT(IN)              :: stars
      TYPE(t_cell), INTENT(IN)              :: cell
      TYPE(t_sphhar), INTENT(IN)              :: sphhar
      TYPE(t_atoms), INTENT(IN)              :: atoms
      TYPE(t_potden), INTENT(IN)              :: den, denRot, EnergyDen
      TYPE(t_potden), INTENT(INOUT)           :: vTot, vx, vxc, exc
      TYPE(t_results), INTENT(INOUT), OPTIONAL :: results
      TYPE(t_potden), INTENT(IN), OPTIONAL     :: den1Rot, den1Rotimag
      TYPE(t_potden), INTENT(INOUT), OPTIONAL  :: dfptvTotimag
      TYPE(t_stars), INTENT(IN), OPTIONAL      :: starsq

      ! Local type instances
      TYPE(t_potden)    :: workDen, veff
      Type(t_kinED)     :: kinED
      REAL, ALLOCATABLE :: tmp_mt(:,:,:), tmp_is(:,:)
      REAL, ALLOCATABLE :: rhoc(:,:,:),rhoc_vx(:)
      REAL, ALLOCATABLE :: tec(:,:), qintc(:,:)
      ! Local Scalars
      INTEGER :: ifftd, ifftd2, ifftxc3d, ispin, i, iType
      REAL    :: dpdot
      LOGICAL :: l_dfptvgen
#ifdef CPP_MPI
      integer:: ierr
#endif

      l_dfptvgen = PRESENT(starsq)

      call set_kinED(fmpi, sphhar, atoms, sym,  xcpot, &
      input, noco, stars,vacuum , cell, Den, EnergyDen, vTot,kinED)

      IF (PRESENT(results)) THEN
         CALL veff%init(stars, atoms, sphhar, vacuum, noco, input%jspins, 1)
#ifndef CPP_OLDINTEL
         ALLOCATE (veff%pw_w, mold=veff%pw)
#else
         ALLOCATE (veff%pw_w(size(veff%pw, 1), size(veff%pw, 2)))
#endif
      ENDIF

      ! exchange correlation potential

      ! vacuum region
      IF (fmpi%irank == 0) THEN
         IF (input%film) THEN
            CALL timestart("Vxc in vacuum")

            ifftd2 = 9*stars%mx1*stars%mx2

            !IF (.NOT. xcpot%needs_grad()) THEN  ! LDA

            CALL vvac_xc(ifftd2, stars, vacuum, noco,   cell, xcpot, input,  Den, vTot, exc)
            CALL timestop("Vxc in vacuum")
         END IF

         ! interstitial region
         CALL timestart("Vxc in interstitial")
         IF (.NOT.l_dfptvgen) THEN
             CALL vis_xc(stars, sym, cell, den, xcpot, input, noco, EnergyDen,kinED, vTot, vx, exc, vxc)
         ELSE
             ! TODO: This is different enough to warrant a separate subroutine, right?
             CALL dfpt_vis_xc(stars, starsq, sym, cell, denRot, den1Rot, xcpot, input, vTot)
         END IF
         CALL timestop("Vxc in interstitial")
      END IF !irank==0

      !
      !     ------------------------------------------
      !     ----> muffin tin spheres region

      IF (fmpi%irank == 0) THEN
         CALL timestart("Vxc in MT")
      END IF

      IF (.NOT.l_dfptvgen) THEN
          CALL vmt_xc(fmpi, sphhar, atoms, den, xcpot, input, sym, &
                      EnergyDen,kinED, noco,vTot, vx, exc, vxc)
      ELSE
          CALL dfpt_vmt_xc(fmpi,sphhar,atoms,denRot,den1Rot,den1Rotimag,xcpot,input,sym,noco,vTot,dfptvTotimag)
      END IF

      ! add MT EXX potential to vr
      IF (fmpi%irank == 0) THEN
         CALL timestop("Vxc in MT")

         ! check continuity of total potential
         IF (input%vchk) CALL checkDOPAll(input,  sphhar, stars, atoms, sym, vacuum,   cell, vTot, 1)

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
            IF (xcpot%is_hybrid() .AND. hybdat%l_subvxc) THEN
               DO ispin = 1, input%jspins
                  CALL convol(stars, vx%pw_w(:, ispin), vx%pw(:, ispin))
               END DO
               veff%pw = vTot%pw - xcpot%get_exchange_weight()*vx%pw
               veff%pw_w = vTot%pw_w - xcpot%get_exchange_weight()*vx%pw_w
               veff%mt = vTot%mt - xcpot%get_exchange_weight()*vx%mt
            END IF

            DO ispin = 1, input%jspins
               DO i = 1, stars%ng3
                  vx%pw(i, ispin) = vx%pw(i, ispin)/stars%nstr(i)
                  vx%pw_w(i, ispin) = vx%pw_w(i, ispin)/stars%nstr(i)
               END DO
            END DO

            results%te_veff = 0.0
            DO ispin = 1, input%jspins
               WRITE (oUnit, FMT=8050) ispin
8050           FORMAT(/, 10x, 'density-effective potential integrals for spin ', i2,/)
               CALL int_nv(ispin, stars, vacuum, atoms, sphhar, cell, sym, input,   veff, workden, results%te_veff)
            END DO

            IF (xcpot%is_hybrid() .AND. hybdat%l_subvxc) THEN

               ALLOCATE(rhoc(atoms%jmtd,atoms%ntype,input%jspins), rhoc_vx(atoms%jmtd))
               ALLOCATE(tec(atoms%ntype,input%jspins),qintc(atoms%ntype,input%jspins))
               CALL readCoreDensity(input,atoms,rhoc,tec,qintc)
               DEALLOCATE(tec,qintc)

               DO ispin = 1, input%jspins
                   DO iType = 1,atoms%ntype
                      rhoc_vx(:atoms%jri(iType)) = rhoc(:atoms%jri(iType),iType,ispin) * &
                                                   vx%mt(:atoms%jri(iType),0,iType,ispin) / sfp_const
                      CALL intgr3(rhoc_vx,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),dpdot)
                      results%te_veff = results%te_veff + xcpot%get_exchange_weight() * dpdot*atoms%neq(iType)
                   END DO
               END DO

               DEALLOCATE(rhoc,rhoc_vx)

            END IF

            WRITE (oUnit, FMT=8060) results%te_veff
8060        FORMAT(/, 10x, 'total density-effective potential integral :', t40, ES20.10)

            ! CALCULATE THE INTEGRAL OF n*exc

            ! perform spin summation of charge densities for the calculation of Exc
            CALL workden%sum_both_spin()

            WRITE (oUnit, FMT=8070)
8070        FORMAT(/, 10x, 'charge density-energy density integrals',/)

            results%te_exc = 0.0
            CALL int_nv(1, stars, vacuum, atoms, sphhar, cell, sym, input,   exc, workDen, results%te_exc)
            WRITE (oUnit, FMT=8080) results%te_exc

8080        FORMAT(/, 10x, 'total charge density-energy density integral :', t40, ES20.10)
         END IF
      END IF ! fmpi%irank == 0

   END SUBROUTINE vgen_xcpot

END MODULE m_vgen_xcpot
