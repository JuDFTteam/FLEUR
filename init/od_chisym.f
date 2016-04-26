      MODULE m_od_chisym
      use m_juDFT
c*************************************************************
c     establishes the rotational matrices and translational
c     vectors in the case of one-dimensional calculations
c     according to the one-dimensional group of symmetries
c                           Y.Mokrousov 2003
c**********************************************
c     CHIRAL SYMMETRY GROUP
c     The main feature is that this group, been represented in 3D space,
c     is cyclic, with the generator <psi|tau> been simultaneous rotation
c     of the system by angle psi around z-axis and translation by tau 
c     along z-direction. The input parameters, given to the program, are:
c          - M (odd%rot)
c          - N (odd%chi),
c     where N is the number of operations in the symmetry group, so that
c     psi = 2pi/N and M is some integer:tau = M/N in the internal units.
      
      CONTAINS
      SUBROUTINE od_chisym(odd,mrot,tau,zrfs,invs,invs2,amat)

      USE m_constants, ONLY : tpi_const
      USE m_types, ONLY : od_dim
      
      IMPLICIT NONE

      TYPE (od_dim), INTENT (IN) :: odd
      REAL, INTENT (IN) :: amat(3,3)
      LOGICAL, INTENT (IN) :: zrfs,invs,invs2
      REAL, INTENT (OUT) :: mrot(3,3,odd%nop),tau(3,odd%nop)      

      INTEGER n,j,i,half
   
      tau(:,:) = 0.
      mrot(1:3,1:3,1:odd%nop) = 0.
      mrot(3,3,:) = 1.
      mrot(1,1,1) = 1.
      mrot(2,2,1) = 1.

      IF (.NOT.odd%invs .AND. .NOT.odd%zrfs) THEN
c for systems without z-reflection and inversion symmetries
        DO n = 2,odd%nop
         mrot(1,1,n) =  cos(tpi_const*(n-1)/odd%nop)
         mrot(1,2,n) =  sin(tpi_const*(n-1)/odd%nop)
         mrot(2,1,n) = -sin(tpi_const*(n-1)/odd%nop)
         mrot(2,2,n) =  cos(tpi_const*(n-1)/odd%nop)
         IF (odd%chi.NE.1) THEN
            tau(3,n) = tau(3,n-1) +
     +           amat(3,3)*(odd%rot)/(odd%chi)
            IF (abs(tau(3,n)-amat(3,3)).LE.1.e-4) THEN
               tau(3,n) = 0.
            ELSEIF (tau(3,n)-amat(3,3).GT.1.e-4) THEN
               tau(3,n) = tau(3,n) - amat(3,3)
            END IF   
         END IF
        END DO
c for systems with inversion or z-relection symmetries
      ELSE
         half = NINT((odd%nop)/2.)
         IF ((odd%nop)/2. - half > 0.001) THEN
            CALL juDFT_error("nop =/= 2*n & (invs v zrfs) !",calledby
     +           ="od_chisym")
         ENDIF
         DO n = 2,half
            mrot(1,1,n) =  cos(tpi_const*(n-1)/half)
            mrot(1,2,n) =  sin(tpi_const*(n-1)/half)
            mrot(2,1,n) = -sin(tpi_const*(n-1)/half)
            mrot(2,2,n) =  cos(tpi_const*(n-1)/half)
         END DO 
         IF (odd%zrfs) THEN
            mrot(1,1,half+1) =  1.
            mrot(2,2,half+1) =  1.
            mrot(3,3,half+1) = -1. 
         ELSE IF (odd%invs) THEN
            mrot(1,1,half+1) = -1.
            mrot(2,2,half+1) = -1.
            mrot(3,3,half+1) = -1.
         END IF
         DO n = 2,half
           !CALL matmul3r(mrot(:,:,n),mrot(:,:,half+1),mrot(:,:,half+n))
           mrot(:,:,half+n)=matmul(mrot(:,:,n),mrot(:,:,half+1))
         END DO
      END IF
      RETURN
      END SUBROUTINE od_chisym
      END MODULE m_od_chisym
