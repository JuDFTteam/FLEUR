!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Assign enpara%vr from the total potential and, for MetaGGA calculations,
!! compute and apply the auxiliary GGA XC correction.
!!
!! For MetaGGA the radial basis functions (u_l, u̇_l) must be generated from
!! a GGA potential rather than the non-multiplicative MetaGGA potential.
!! This routine computes vr = v_tot^{l=0} + (v_xc^{auxGGA} - v_xc^{MGGA})
!! where the auxiliary GGA potential is obtained by calling vmt_xc with a
!! temporary xcpot created from the auxiliary GGA functional IDs.
!!
!! Reference: Doumont et al., Phys. Rev. B 106, 235159 (2022), Section II.B
MODULE m_assign_enpara_potential
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: assign_enpara_potential

CONTAINS

  SUBROUTINE assign_enpara_potential(enpara, fmpi, atoms, input, v, xcpot, vxc, inDen, &
                                     sphhar, sym, stars, vacuum, noco, EnergyDen)
    USE m_types_enpara
    USE m_types_atoms
    USE m_types_input
    USE m_types_potden
    USE m_types_xcpot
    USE m_types_sphhar
    USE m_types_sym
    USE m_types_stars
    USE m_types_vacuum
    USE m_types_noco
    USE m_types_mpi
    USE m_vmt_xc
    USE m_constants
    USE m_judft
    IMPLICIT NONE

    TYPE(t_enpara), INTENT(INOUT)        :: enpara
    TYPE(t_mpi), INTENT(IN)              :: fmpi
    TYPE(t_atoms), INTENT(IN)            :: atoms
    TYPE(t_input), INTENT(IN)            :: input
    TYPE(t_potden), INTENT(IN)           :: v
    CLASS(t_xcpot), INTENT(IN)           :: xcpot
    TYPE(t_potden), INTENT(IN)           :: vxc, inDen, EnergyDen
    TYPE(t_sphhar), INTENT(IN)           :: sphhar
    TYPE(t_sym), INTENT(IN)              :: sym
    TYPE(t_stars), INTENT(IN)            :: stars
    TYPE(t_vacuum), INTENT(IN)           :: vacuum
    TYPE(t_noco), INTENT(IN)             :: noco

    ! Local variables
    INTEGER :: n, jsp
    CLASS(t_xcpot), ALLOCATABLE :: xcpot_aux
    TYPE(t_potden)   :: vTot_aux, vx_aux, vxc_aux, exc_aux, vTau_aux

    ! Assign enpara%vr from the spherical (l=0) component of the total potential
    DO jsp = 1, input%jspins
      DO n = 1, atoms%ntype
        enpara%vr(:, n, jsp) = v%mt(:, 0, n, jsp)
      END DO
    END DO

    ! For MetaGGA: replace MGGA XC potential with auxiliary GGA XC potential
    ! vr = v_Coulomb + v_xc^GGA = (v_Coulomb + v_xc^MGGA) + (v_xc^GGA - v_xc^MGGA)
    IF (xcpot%needs_MetaGGA_ham()) THEN
        
      CALL timestart("Auxiliary GGA for basis")

      ! Create temporary xcpot with auxiliary GGA functional IDs
      CALL xcpot%create_from_aux(xcpot_aux)

      ! Initialize temporary output potentials for vmt_xc
      CALL vTot_aux%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_POTTOT)
      CALL vx_aux%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_POTTOT)
      CALL vxc_aux%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_POTTOT)
      vxc_aux%mt=0.0
      CALL exc_aux%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_POTTOT)
      CALL vTau_aux%init(stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_POTTOT)

      ! Compute auxiliary GGA XC potential via vmt_xc
      CALL vmt_xc(fmpi, sphhar, atoms, inDen, xcpot_aux, input, sym, &
                  EnergyDen, noco, vTot_aux, vx_aux, exc_aux, vxc_aux, vTau_aux)

      ! Apply correction: vr += (v_xc^auxGGA - v_xc^MGGA)
      DO jsp = 1, input%jspins
        DO n = 1, atoms%ntype
          enpara%vr(:, n, jsp) = enpara%vr(:, n, jsp) &
                               + (vxc_aux%mt(:, 0, n, jsp) - vxc%mt(:, 0, n, jsp))*atoms%rmsh(:,n)/sfp_const
          write(77,*) n," vr"
          write(77,*) enpara%vr(:20, n, jsp)                     
          write(77,*) n," v_aux"
          write(77,*) vxc_aux%mt(:20, 0, n, jsp)
          write(77,*) n," v_xc"
          write(77,*) vxc%mt(:20, 0, n, jsp)
        END DO
      END DO

      CALL timestop("Auxiliary GGA for basis")
    END IF

    
     DO n = 1, atoms%ntype
        IF (atoms%l_nonpolbas(n)) THEN
          enpara%vr(:, n, 1) = (v%mt(:, 0, n, 1) + v%mt(:, 0, n, jsp)) / 2
          enpara%vr(:, n, jsp) = v%mt(:, 0, n, 1)
        ENDIF  
     END DO
    
  END SUBROUTINE assign_enpara_potential

END MODULE m_assign_enpara_potential
