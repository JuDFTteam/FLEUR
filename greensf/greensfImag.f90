MODULE m_greensfImag

USE m_juDFT
USE m_types
USE m_constants


CONTAINS

SUBROUTINE greensfImag(atoms,sym,input,ispin,nbands,dosWeights,resWeights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)

   !This Subroutine calculates the contribution to the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !of the current k-Point (it is called in cdnval) inside the MT-sphere
   !and sums over the Brillouin-Zone using the histogram method or linear tetrahedron method
   !It is essentially the l-density of states in a (m,mp) matrix with an additional factor - pi

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
   TYPE(t_usdus),         INTENT(IN)    :: usdus
   TYPE(t_greensfCoeffs), INTENT(INOUT) :: greensfCoeffs
   TYPE(t_input),         INTENT(IN)    :: input

   !-Scalar Arguments
   INTEGER,               INTENT(IN)    :: ispin  !Current spin index
   INTEGER,               INTENT(IN)    :: nbands !Number of bands to be considered
   REAL,                  INTENT(IN)    :: wtkpt  !Weight of the current k-point (not used in tetrahedron method)

   !-Array Arguments
   REAL,                  INTENT(IN)    :: dosWeights(greensfCoeffs%ne,nbands) !Precalculated tetrahedron weights for the current k-point
   REAL,                  INTENT(IN)    :: resWeights(greensfCoeffs%ne,nbands) !Precalculated tetrahedron weights for the current k-point
   INTEGER,               INTENT(IN)    :: ind(nbands,2)                       !Gives the range where the tetrahedron weights are non-zero
   REAL,                  INTENT(IN)    :: eig(nbands)                         !Eigenvalues for the current k-point

   !-Local Scalars
   LOGICAL l_zero,l_tria
   INTEGER i_gf,ib,ie,j,nType,natom,l,m,mp,lm,lmp,ilo,ilop,imat,it,is,isi,i
   REAL    fac
   COMPLEX weight
   COMPLEX im(greensfCoeffs%ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(1,5,input%l_gfsphavg))
   COMPLEX d_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),calc_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   l_tria = (input%tria.OR.input%gfTet).AND..NOT.input%l_hist

   IF(l_tria.AND.(ANY(ind.GT.greensfCoeffs%ne).OR.ANY(ind.LT.1))) THEN
      CALL juDFT_error("Invalid index",calledby="greensfImag")
   ENDIF

   !Loop through the gf elements to be calculated
   DO i_gf = 1, atoms%n_gf

      l     = atoms%gfelem(i_gf)%l
      nType = atoms%gfelem(i_gf)%atomType

      !Loop through equivalent atoms
      DO natom = SUM(atoms%neq(:nType-1)) + 1, SUM(atoms%neq(:nType))
         im = 0.0
         !Loop through bands
         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(natom,l,nType,ispin,wtkpt,i_gf,nbands,l_tria) &
         !$OMP SHARED(atoms,im,input,eigVecCoeffs,usdus,greensfCoeffs,eig,dosWeights,resWeights,ind) &
         !$OMP PRIVATE(ie,m,mp,lm,lmp,ilo,ilop,weight,ib,j,l_zero)
         !$OMP DO
         DO ib = 1, nbands
            !Check wether there is a non-zero weight for the energy window
            l_zero = .true.
            IF(l_tria) THEN
               IF(.NOT.input%l_resolvent) THEN
                  !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                  IF(ANY(dosWeights(ind(ib,1):ind(ib,2),ib).NE.0.0)) l_zero = .false.
               ELSE
                  l_zero = .false.
               ENDIF
            ELSE
               !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
               j = NINT((eig(ib)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
               IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) )         l_zero = .false.
            END IF

            IF(l_zero) CYCLE

            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp
                  !Choose the relevant energy points depending on the bz-integration method
                  DO ie = MERGE(ind(ib,1),j,l_tria), MERGE(ind(ib,2),j,l_tria)
                     !weight for the bz-integration including spin-degeneracy
                     weight = 2.0/input%jspins*(MERGE(resWeights(ie,ib),0.0,input%l_resolvent)&
                                                - ImagUnit * pi_const * MERGE(dosWeights(ie,ib),wtkpt/greensfCoeffs%del,l_tria))
                     !-------------------------
                     !Contribution from states
                     !-------------------------
                     im(ie,m,mp,1) = im(ie,m,mp,1) + weight *&
                                          (conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                          conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                          usdus%ddn(l,nType,ispin))
                     IF(.NOT.input%l_gfsphavg) THEN
                        im(ie,m,mp,2) = im(ie,m,mp,2) + weight *&
                                        conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin)
                        im(ie,m,mp,3) = im(ie,m,mp,3) + weight *&
                                        conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)
                        im(ie,m,mp,4) = im(ie,m,mp,4) + weight *&
                                        conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)
                        im(ie,m,mp,5) = im(ie,m,mp,5) + weight *&
                                        conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin)
                     END IF
                     !------------------------------------------------------------------------------------------------------
                     ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
                     !------------------------------------------------------------------------------------------------------
                     DO ilo = 1, atoms%nlo(nType)
                        IF(atoms%llo(ilo,nType).EQ.l) THEN
                           im(ie,m,mp,1) = im(ie,m,mp,1) + weight * (  usdus%uulon(ilo,nType,ispin) * (&
                                          conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                          conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) )&
                                          + usdus%dulon(ilo,nType,ispin) * (&
                                          conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                          conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)))
                           DO ilop = 1, atoms%nlo(nType)
                              IF (atoms%llo(ilop,nType).EQ.l) THEN
                                 im(ie,m,mp,1) = im(ie,m,mp,1) + weight * usdus%uloulopn(ilo,ilop,nType,ispin) *&
                                                conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,ib,ilo,natom,ispin)
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO! ie
               ENDDO !mp
            ENDDO !m
         ENDDO !ib
         !$OMP END DO
         !$OMP END PARALLEL

         !Rotate the eqivalent atom into the irreducible brillouin zone
         fac = 1.0/(sym%invarind(natom)*atoms%neq(nType))
         IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No symmetry operations available",calledby="greensfImag")
         DO imat = 1, MERGE(1,5,input%l_gfsphavg)
            DO it = 1, sym%invarind(natom)
               is = sym%invarop(natom,it)
               isi = sym%invtab(is)
               d_mat(:,:) = cmplx(0.0,0.0)
               DO m = -l,l
                  DO mp = -l,l
                     d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
                  ENDDO
               ENDDO
               DO ie = 1, greensfCoeffs%ne
                  calc_mat = matmul( transpose( conjg(d_mat) ) , &
                              im(ie,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,imat))
                  calc_mat =  matmul( calc_mat, d_mat )
                  DO m = -l,l
                     DO mp = -l,l
                        IF(imat.EQ.1) THEN
                           greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) - AIMAG(fac * conjg(calc_mat(m,mp)))
                        ELSE IF(imat.EQ.2) THEN
                           greensfCoeffs%uu(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) - AIMAG(fac * conjg(calc_mat(m,mp)))
                        ELSE IF(imat.EQ.3) THEN
                           greensfCoeffs%dd(ie,i_gf,m,mp,ispin) = greensfCoeffs%dd(ie,i_gf,m,mp,ispin) - AIMAG(fac * conjg(calc_mat(m,mp)))
                        ELSE IF(imat.EQ.4) THEN
                           greensfCoeffs%ud(ie,i_gf,m,mp,ispin) = greensfCoeffs%ud(ie,i_gf,m,mp,ispin) - AIMAG(fac * conjg(calc_mat(m,mp)))
                        ELSE IF(imat.EQ.5) THEN
                           greensfCoeffs%du(ie,i_gf,m,mp,ispin) = greensfCoeffs%du(ie,i_gf,m,mp,ispin) - AIMAG(fac * conjg(calc_mat(m,mp)))
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO!it
            ENDDO!ie
         ENDDO!imat
      ENDDO!natom
   ENDDO !i_gf

END SUBROUTINE greensfImag

END MODULE