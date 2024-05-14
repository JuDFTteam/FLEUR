MODULE m_dftUPotential

   USE m_constants
   USE m_juDFT
   USE m_types
   USE m_uj2f
   USE m_doubleCounting
   USE m_coulombPotential

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dftUPotential(density, ldau, jspins, equivalentAtoms, l_spinoffd, potential, ldaUEnergy, spinavg_dc)

      !This subroutine calculates the DFT+U potential matrix

      COMPLEX,          INTENT(IN)     :: density(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_utype),    INTENT(IN)     :: ldau
      INTEGER,          INTENT(IN)     :: jspins
      INTEGER,          INTENT(IN)     :: equivalentAtoms
      LOGICAL,          INTENT(IN)     :: l_spinoffd
      COMPLEX,          INTENT(INOUT)  :: potential(-lmaxU_const:,-lmaxU_const:,:)
      REAL,             INTENT(INOUT)  :: ldaUEnergy
      LOGICAL, OPTIONAL,INTENT(IN)     :: spinavg_dc


      LOGICAL :: spinavg_dc_local
      INTEGER :: m,ispin,jspin
      REAL    :: u_htr,j_htr,f0,f2,f4,f6,energy_contribution,total_charge, double_counting
      COMPLEX, ALLOCATABLE :: Vdc(:,:,:)

      spinavg_dc_local = .false.
      if(present(spinavg_dc)) spinavg_dc_local = spinavg_dc .and. jspins==2

      call coulombPotential(density, ldau, jspins, l_spinoffd, potential, energy_contribution)

      CALL uj2f(jspins,ldau,f0,f2,f4,f6)
      u_htr = f0/hartree_to_ev_const
      IF (ldau%l.EQ.1) THEN
         j_htr = f2/(5*hartree_to_ev_const)
      ELSE IF (ldau%l.EQ.2) THEN
         j_htr = 1.625*f2/(14*hartree_to_ev_const)
      ELSE IF (ldau%l.EQ.3) THEN
         j_htr = (286.+195.*451./675.+250.*1001./2025.)*f2/(6435*hartree_to_ev_const)
      END IF

      !Add double counting terms
      Vdc = doubleCountingPot(density, ldau, u_htr, j_htr, l_spinoffd,&
                              .FALSE.,spinavg_dc_local,0.0)
      potential = potential - Vdc

      double_counting = doubleCountingEnergy(density, ldau, u_htr, j_htr, l_spinoffd,&
                                             .FALSE.,spinavg_dc_local,0.0)

      if(ldau%l_amf) then
         energy_contribution = energy_contribution - double_counting
      else
         energy_contribution = energy_contribution - double_counting
      endif
      
      total_charge = 0.0
      DO ispin = 1, MIN(2,SIZE(density,3))
         DO m = -ldau%l, ldau%l
            total_charge = total_charge +  REAL(density(m,m,ispin))
         ENDDO
      ENDDO
      energy_contribution = energy_contribution - (u_htr-j_htr)/2.0 * total_charge
      
      ldaUEnergy = ldaUEnergy + energy_contribution * equivalentAtoms

   END SUBROUTINE dftUPotential

END MODULE m_dftUPotential
