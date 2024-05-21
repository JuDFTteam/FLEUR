!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_doubleCounting

   USE m_constants
   USE m_types
   USE m_coulombPotential

   IMPLICIT NONE

   CONTAINS

   FUNCTION doubleCountingPot(density, ldau, U,J, l_spinoffd, l_mix,l_spinAvg, alpha, l_write) RESULT(Vdc)

      !------------------------------------------------------
      ! Calculate the Double Counting Correction in either
      ! the FLL or AMF limit
      ! If l_spinAvg is True the Double counting will be
      ! averaged over the spins
      !------------------------------------------------------

      COMPLEX,             INTENT(IN)  :: density(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_utype),       INTENT(IN)  :: ldau       !LDA+U information
      REAL,                INTENT(IN)  :: U, J
      LOGICAL,             INTENT(IN)  :: l_spinoffd
      LOGICAL,             INTENT(IN)  :: l_mix      !Mix between FLL and AMF
      LOGICAL,             INTENT(IN)  :: l_spinAvg  !Do we want a spin averaged double counting
      REAL,                INTENT(IN)  :: alpha
      LOGICAL, OPTIONAL,   INTENT(IN)  :: l_write

      COMPLEX :: Vdc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(density,3))

      COMPLEX, ALLOCATABLE :: modified_density(:,:,:)
      REAL :: charge, mag(3),mag_m(0:3), tmp, D, Vdcup,Vdcdn
      COMPLEX :: sigma(2,2,3), r21
      INTEGER :: spin_dim, ispin,m, spin1,spin2,mp
      type(t_nococonv) :: nococonv !Used only for the procedures on it

      sigma = cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      spin_dim = SIZE(density,3)
      IF(.NOT.l_spinoffd) spin_dim = MIN(2,spin_dim)

      charge = 0.0
      mag = 0.0

      DO m = -ldau%l, ldau%l
         if (spin_dim==3) then
            r21 = density(m,m,3)
         else
            r21 = 0.0
         endif
         mag_m = nococonv%denmat_to_mag(real(density(m,m,1)),&
                                        real(density(m,m,min(2,spin_dim))),&
                                        r21)

         charge = charge + mag_m(0)
         mag = mag + mag_m(1:)
      ENDDO
      IF(spin_dim == 1) then
         ! In the l_amf case the spin degeneracy will be dealt with in the coulombPotential call.
         IF(.NOT.ldau%l_amf) charge=charge/2.0
         mag = 0.0
      ENDIF

      Vdc = cmplx_0
      IF(ldau%l_amf) THEN

         ALLOCATE(modified_density(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, SIZE(density,3)), source=cmplx_0)
         modified_density = cmplx_0
         D = real(2*(2*ldau%l+1))

         DO ispin = 1, spin_dim

            IF(ispin==3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = ispin
               spin2 = ispin
            ENDIF
         
            DO m = -ldau%l, ldau%l
               IF(spin1==spin2) modified_density(m,m,ispin) = charge/D
               modified_density(m,m,ispin) = modified_density(m,m,ispin) + dot_product(mag,sigma(spin1,spin2,:))/D
            ENDDO

         ENDDO
         
         call coulombPotential(modified_density,ldau, MIN(2,SIZE(density,3)), l_spinoffd,Vdc,tmp)
         
      ELSE
         DO ispin = 1, spin_dim

            IF(ispin==3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = ispin
               spin2 = ispin
            ENDIF

            DO m = -ldau%l, ldau%l
               IF(spin1 == spin2) Vdc(m,m,ispin) = Vdc(m,m,ispin) + U*(charge-0.5) - J*(charge/2.0-0.5)
               Vdc(m,m,ispin) = Vdc(m,m,ispin) - J/2.0*dot_product(mag,sigma(spin1,spin2,:))
            ENDDO
         ENDDO
      ENDIF
      

      IF(PRESENT(l_write)) THEN
         IF(l_write) THEN
            WRITE(oUnit,"(/,A)") 'Double counting chemical potential:'
            IF(ldau%l_amf) THEN
               WRITE(oUnit,9040) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
            ELSE
               WRITE(oUnit,9040) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
            ENDIF
9040        FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
            Vdcup = 0.0
            Vdcdn = 0.0
            do m = -ldau%l, ldau%l
               Vdcup = Vdcup + Vdc(m,m,1)/real(2*ldau%l+1)
               Vdcdn = Vdcdn + Vdc(m,m,min(2,spin_dim))/real(2*ldau%l+1)
            enddo
            WRITE(oUnit,9050) Vdcup,Vdcdn,(Vdcup+Vdcdn)/2.0,Vdcup-Vdcdn
9050        FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
         ENDIF
      ENDIF

      IF(l_spinAvg) THEN
         Vdc(:,:,1) = (Vdc(:,:,1)+Vdc(:,:,min(2,spin_dim)))/2.0
         if(spin_dim>1) Vdc(:,:,2) = Vdc(:,:,1)
         if(spin_dim==3) Vdc(:,:,3) = cmplx_0 !Is this right?
      ENDIF

   END FUNCTION doubleCountingPot


   REAL FUNCTION doubleCountingEnergy(density, ldau, U,J, l_spinoffd, l_mix,l_spinAvg, alpha, l_write)

      !------------------------------------------------------------
      ! Calculate the Double Counting Correction Energy in either
      ! the FLL or AMF limit
      !------------------------------------------------------------

      COMPLEX,             INTENT(IN)  :: density(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_utype),       INTENT(IN)  :: ldau       !LDA+U information
      REAL,                INTENT(IN)  :: U, J
      LOGICAL,             INTENT(IN)  :: l_spinoffd
      LOGICAL,             INTENT(IN)  :: l_mix      !Mix between FLL and AMF
      LOGICAL,             INTENT(IN)  :: l_spinAvg  !Do we want a spin averaged double counting
      REAL,                INTENT(IN)  :: alpha
      LOGICAL, OPTIONAL,   INTENT(IN)  :: l_write

      REAL :: charge, mag(3), mag_m(0:3), D
      INTEGER :: spin_dim, ispin,m, spin1,spin2
      COMPLEX :: tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(density,3)), r21
      COMPLEX, ALLOCATABLE :: modified_density(:,:,:)
      COMPLEX :: sigma(2,2,3)
      type(t_nococonv) :: nococonv !Used only for the procedures on it

      sigma = cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      spin_dim = SIZE(density,3)
      IF(.NOT.l_spinoffd) spin_dim = MIN(2,spin_dim)

      charge = 0.0
      mag = 0.0
      DO m = -ldau%l, ldau%l
         if (spin_dim==3) then
            r21 = density(m,m,3)
         else
            r21 = 0.0
         endif
         mag_m = nococonv%denmat_to_mag(real(density(m,m,1)),&
                                        real(density(m,m,min(2,spin_dim))),&
                                        r21)
         charge = charge + mag_m(0)
         mag = mag + mag_m(1:)
      ENDDO

      IF(spin_dim == 1) then
         charge=charge/2.0
         mag = 0.
      endif
      if (l_spinAvg) mag = 0.

      doubleCountingEnergy = 0.0
      IF(ldau%l_amf) THEN
         ALLOCATE(modified_density(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const, SIZE(density,3)), source=cmplx_0)
         modified_density = cmplx_0
         D = real(2*(2*ldau%l+1)) ! This factor averages over the different m and spin

         DO ispin = 1, spin_dim

            IF(ispin==3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = ispin
               spin2 = ispin
            ENDIF
         
            DO m = -ldau%l, ldau%l
               IF(spin1==spin2) modified_density(m,m,ispin) = charge/D
               modified_density(m,m,ispin) = modified_density(m,m,ispin) + dot_product(mag,sigma(spin1,spin2,:))/D
            ENDDO

         ENDDO
         
         IF(spin_dim == 1) THEN
            modified_density = modified_density * 2.0
         END IF
         
         call coulombPotential(modified_density,ldau, MIN(2,SIZE(density,3)), l_spinoffd,tmp,doubleCountingEnergy)
      ELSE
         doubleCountingEnergy = U/2*charge*(charge-1) -J/2*charge*(charge/2-1)-J*dot_product(mag,mag)/4
      ENDIF

   END FUNCTION doubleCountingEnergy

   FUNCTION doubleCountingMixFactor(mmpmat, l, charge) Result(alpha)

      !---------------------------------------------
      ! Calculate the mixing factor between FLL/AMF
      !---------------------------------------------

      COMPLEX,        INTENT(IN) :: mmpmat(-lmaxU_const:, -lmaxU_const:, :)
      INTEGER,        INTENT(IN) :: l
      REAL,           INTENT(IN) :: charge

      REAL :: alpha

      alpha = 0.0

   END FUNCTION

END MODULE m_doubleCounting
