MODULE m_onsite21
   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_onsite21
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module handles the spin-offdiagonal part of the imaginary part of the onsite
   !>  Green's function in the noco
   !
   !> In the non-collinear case spin isn't a good quantum number anymore so we get a 2x2 blockatrix
   !
   ! REVISION HISTORY:
   ! 14 March 2019 - Initial Version
   !------------------------------------------------------------------------------

   USE m_types
   USE m_juDFT
   USE m_constants

   CONTAINS

   SUBROUTINE onsite21(atoms,input,nbands,tetweights,ind,wtkpt,eig,denCoeffsOffDiag,eigVecCoeffs,greensfCoeffs)


      IMPLICIT NONE

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs

      INTEGER,                   INTENT(IN)     :: nbands   
      REAL,                      INTENT(IN)     :: wtkpt
      REAL,                      INTENT(IN)     :: tetweights(greensfCoeffs%ne,nbands)
      INTEGER,                   INTENT(IN)     :: ind(nbands,2)
      REAL,                      INTENT(IN)     :: eig(nbands) 

      INTEGER  i_gf,nType,l,natom,ib,j,ie,m,lm,mp,lmp,ispin,spin1,spin2
      REAL     weight
      LOGICAL  l_zero

      IF(.NOT.input%l_gfsphavg) CALL juDFT_error("NOCO-offdiagonal + Radial dependence of onsite-GF not implemented",calledby="onsite21")

      DO i_gf = 1, atoms%n_gf
         nType = atoms%gfelem(i_gf)%atomType
         l = atoms%gfelem(i_gf)%l


         DO natom = SUM(atoms%neq(:nType-1)) + 1, SUM(atoms%neq(:nType))
            
            DO ib = 1, nbands

               l_zero = .true.
               IF(input%tria) THEN
                  !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                  IF(ANY(tetweights(:,ib).NE.0.0)) l_zero = .false.
               ELSE
                  !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
                  j = NINT((eig(ib)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
                  IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) ) l_zero = .false.
               END IF

               IF(l_zero) CYCLE

               DO ie = MERGE(ind(ib,1),j,input%tria), MERGE(ind(ib,2),j,input%tria)

                  weight = 1.0/atoms%neq(nType) * MERGE(tetweights(ie,ib),wtkpt/greensfCoeffs%del,input%tria)
                  DO ispin = 3, 4
                     IF(ispin.EQ.3) THEN
                        spin1 = 2
                        spin2 = 1
                     ELSE
                        spin1 = 1
                        spin2 = 2
                     ENDIF
                     DO m = -l, l 
                        lm = l*(l+1) + m 
                        DO mp = -l, l 
                           lmp = l*(l+1) + mp 
                           greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) - pi_const * weight * &
                                                                     AIMAG( eigVecCoeffs%acof(ib,lm,natom,spin1) * CONJG(eigVecCoeffs%acof(ib,lmp,natom,spin2)) * denCoeffsOffdiag%uu21n(l,nType)&
                                                                     + eigVecCoeffs%bcof(ib,lm,natom,spin1) * CONJG(eigVecCoeffs%bcof(ib,lmp,natom,spin2)) * denCoeffsOffdiag%dd21n(l,nType)&
                                                                     + eigVecCoeffs%acof(ib,lm,natom,spin1) * CONJG(eigVecCoeffs%bcof(ib,lmp,natom,spin2)) * denCoeffsOffdiag%du21n(l,nType)&
                                                                     + eigVecCoeffs%bcof(ib,lm,natom,spin1) * CONJG(eigVecCoeffs%acof(ib,lmp,natom,spin2)) * denCoeffsOffdiag%ud21n(l,nType))                 
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO !ie
            ENDDO !ib
         ENDDO !natom
      ENDDO !i_gf

   END SUBROUTINE onsite21
   SUBROUTINE rot_onsite(atoms,noco,gf)

      USE m_rotdenmat

      IMPLICIT NONE
       
      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_greensf),        INTENT(INOUT)  :: gf

      INTEGER  i_gf,nType,l,m,mp,ie,ipm
      REAL     gf11re,gf22re,gf11imag,gf22imag
      COMPLEX  gf21re,gf21imag

      DO i_gf = 1, atoms%n_gf
         nType = atoms%gfelem(i_gf)%atomType
         l = atoms%gfelem(i_gf)%l

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(i_gf,nType,l) &
         !$OMP SHARED(atoms,noco,gf) &
         !$OMP PRIVATE(ie,m,mp,ipm)&
         !$OMP PRIVATE(gf11re,gf22re,gf11imag,gf22imag,gf21re,gf21imag)
         !$OMP DO
         DO ie = 1, gf%nz
            DO m = -l, l 
               DO mp = -l, l
                  DO ipm = 1, 2 
                     !We need to call rot_den_mat two times for real and imaginary part
                     gf11re = REAL(gf%gmmpMat(ie,i_gf,m,mp,1,ipm))
                     gf22re = REAL(gf%gmmpMat(ie,i_gf,m,mp,2,ipm))
                     CALL rot_den_mat(noco%alph(nType),noco%beta(nType),gf11re,gf22re,gf21re)
                     gf11imag = AIMAG(gf%gmmpMat(ie,i_gf,m,mp,1,ipm))
                     gf22imag = AIMAG(gf%gmmpMat(ie,i_gf,m,mp,2,ipm))
                     CALL rot_den_mat(noco%alph(nType),noco%beta(nType),gf11imag,gf22imag,gf21imag)

                     gf%gmmpMat(ie,i_gf,m,mp,1,ipm) = gf11re + ImagUnit * gf11imag 
                     gf%gmmpMat(ie,i_gf,m,mp,2,ipm) = gf22re + ImagUnit * gf22imag
                     gf%gmmpMat(ie,i_gf,m,mp,3,ipm) = gf21re + ImagUnit * gf21imag
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO
   END SUBROUTINE rot_onsite

END MODULE m_onsite21