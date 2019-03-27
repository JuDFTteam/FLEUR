MODULE m_onsite21
   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_onsite21
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module handles the offdiagonal part of the imaginary part of the onsite
   !>  Green's function (spehrically averaged only atm)
   !
   ! REVISION HISTORY:
   ! 14 March 2019 - Initial Version
   !------------------------------------------------------------------------------

   USE m_juDFT
   USE m_constants

   CONTAINS

   SUBROUTINE onsite21(atoms,input,sym,jspins,noccbd,tetweights,wtkpt,eig,usdus,denCoeffsOffDiag,eigVecCoeffs,greensfCoeffs)

      USE m_types
      USE m_differentiate

      IMPLICIT NONE

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs

      INTEGER,                   INTENT(IN)     :: jspins
      INTEGER,                   INTENT(IN)     :: noccbd 

      REAL,                      INTENT(IN)     :: wtkpt
      REAL,                      INTENT(IN)     :: tetweights(:,:)
      REAL,                      INTENT(IN)     :: eig(noccbd)

      INTEGER  i_gf,n,l,natom,nn,i,j,m,lm,mp,lmp
      REAL fac, wk
      REAL, ALLOCATABLE :: dos_weights(:)
      LOGICAL l_zero

      COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      wk = wtkpt/greensfCoeffs%del

      ALLOCATE(dos_weights(greensfCoeffs%ne))

      IF(.NOT.input%onsite_sphavg) CALL juDFT_error("NOCO-offdiagonal + Radial dependence of&
                                                   onsite-GF not implemented",calledby="onsite21")

      DO i_gf = 1, atoms%n_gf
         n = atoms%onsiteGF(i_gf)%atomType
         l = atoms%onsiteGF(i_gf)%l

         natom = SUM(atoms%neq(:n-1))
         DO nn = 1, atoms%neq(n)
            natom = natom + 1
            fac = 1.0 / ( sym%invarind(natom) * atoms%neq(n) )
            
            DO i = 1, noccbd
               l_zero = .true.
               IF(greensfCoeffs%l_tetra) THEN
                  !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                  IF(ANY(tetweights(:,i).NE.0.0)) l_zero = .false.
               ELSE
                  !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
                  j = NINT((eig(i)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
                  IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) ) l_zero = .false.
               END IF

               IF(l_zero) CYCLE

               n_tmp = 0.0
               DO m = -l, l 
                  lm = l*(l+1) + m 
                  DO mp = -l, l 
                     lmp = l*(l+1) + mp 
                     n_tmp(m,mp) = n_tmp(m,mp) - pi_const * ( eigVecCoeffs%acof(i,lm,natom,1) * CONJG(eigVecCoeffs%acof(i,lmp,natom,input%jspins)) * denCoeffsOffdiag%uu21n(l,n)&
                                                            + eigVecCoeffs%bcof(i,lm,natom,1) * CONJG(eigVecCoeffs%bcof(i,lmp,natom,input%jspins)) * denCoeffsOffdiag%dd21n(l,n)&
                                                            + eigVecCoeffs%acof(i,lm,natom,1) * CONJG(eigVecCoeffs%bcof(i,lmp,natom,input%jspins)) * denCoeffsOffdiag%du21n(l,n)&
                                                            + eigVecCoeffs%bcof(i,lm,natom,1) * CONJG(eigVecCoeffs%acof(i,lmp,natom,input%jspins)) * denCoeffsOffdiag%ud21n(l,n))
                  ENDDO
               ENDDO
               DO m = -l,l
                  DO mp = -l,l
                     IF(greensfCoeffs%l_tetra) THEN
                        !We need to differentiate the weights with respect to energy (can maybe be done analytically)
                        CALL diff3(tetweights(:,i),greensfCoeffs%del,dos_weights(:))
                        DO j = 1, greensfCoeffs%ne
                           greensfCoeffs%im_g21(j,i_gf,m,mp) = greensfCoeffs%im_g21(j,i_gf,m,mp) + n_tmp(m,mp) * fac * dos_weights(j)
                        ENDDO
                     ELSE    
                        greensfCoeffs%im_g21(j,i_gf,m,mp) = greensfCoeffs%im_g21(j,i_gf,m,mp) + n_tmp(m,mp) * fac * wk
                     END IF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE onsite21

   SUBROUTINE rot_gf_mat(atoms,noco,greensfCoeffs)

      USE m_types
      USE m_rotdenmat

      IMPLICIT NONE
       
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs

      INTEGER  i_gf,n,l,m,mp,j

      DO i_gf = 1, atoms%n_gf
         n = atoms%onsiteGF(i_gf)%atomType
         l = atoms%onsiteGF(i_gf)%l

         DO m = -l, l 
            DO mp = -l, l 
               DO j = 1, greensfCoeffs%ne 
                  CALL rot_den_mat(noco%alph(n),noco%beta(n),greensfCoeffs%im_g(j,i_gf,m,mp,1),&
                              greensfCoeffs%im_g(j,i_gf,m,mp,2),greensfCoeffs%im_g21(j,i_gf,m,mp))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE rot_gf_mat

END MODULE m_onsite21