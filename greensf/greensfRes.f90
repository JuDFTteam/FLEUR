MODULE m_greensfRes



   USE m_juDFT
USE m_types
USE m_constants


CONTAINS

SUBROUTINE greensfRes(atoms,sym,input,ispin,nbands,tetweights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensf)

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
   TYPE(t_greensf), INTENT(INOUT) :: greensf
   TYPE(t_input),         INTENT(IN)    :: input 

   !-Scalar Arguments 
   INTEGER,               INTENT(IN)    :: ispin  !Current spin index
   INTEGER,               INTENT(IN)    :: nbands !Number of bands to be considered
   REAL,                  INTENT(IN)    :: wtkpt  !Weight of the current k-point

   !-Array Arguments
   COMPLEX,               INTENT(IN)    :: tetweights(greensf%nz,nbands) !Precalculated tetrahedron weights for the current k-point
   INTEGER,               INTENT(IN)    :: ind(nbands,2)                       !Gives the range where the tetrahedron weights are non-zero
   REAL,                  INTENT(IN)    :: eig(nbands)                         !Eigenvalues for the current k-point

   !-Local Scalars
   LOGICAL l_zero,l_tria
   INTEGER i_gf,ib,iz,j,nType,natom,l,m,mp,lm,lmp,ilo,ilop,imat,it,is,isi
   REAL    fac
   COMPLEX weight
   COMPLEX g(greensf%nz,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(1,5,input%l_gfsphavg))
   COMPLEX d_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),calc_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   l_tria = input%tria.OR.input%gfTet
   IF(.NOT.l_tria) CALL juDFT_error("Only available with tetra",calledby="greensfRes")
   IF(.NOT.input%l_gfsphavg) CALL juDFT_error("No radial dependence",calledby="greensfRes")
   !Loop through the gf elements to be calculated
   DO i_gf = 1, atoms%n_gf

      l     = atoms%gfelem(i_gf)%l
      nType = atoms%gfelem(i_gf)%atomType

      !Loop through equivalent atoms
      DO natom = SUM(atoms%neq(:nType-1)) + 1, SUM(atoms%neq(:nType))
         g = 0.0
         !Loop through bands
         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(natom,l,nType,ispin,wtkpt,i_gf,nbands,l_tria) &
         !$OMP SHARED(atoms,g,input,eigVecCoeffs,usdus,greensf,eig,tetweights,ind) &
         !$OMP PRIVATE(iz,m,mp,lm,lmp,ilo,ilop,weight,ib,j,l_zero)
         !$OMP DO
         DO ib = 1, nbands
            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp
                  !Choose the relevant energy points depending on the bz-integration method
                  DO iz = 1, greensf%nz
                     !weight for the bz-integration including spin-degeneracy
                     weight = 2.0/input%jspins * tetweights(iz,ib) 
                     !
                     !Contribution from states
                     !
                     g(iz,m,mp,1) = g(iz,m,mp,1) + weight *&
                                    (conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                          conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                          usdus%ddn(l,nType,ispin))
                     !IF(.NOT.input%l_gfsphavg) THEN
                     !   im(ie,m,mp,2) = im(ie,m,mp,2) -  pi_const * weight *&
                     !                   conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                     !   im(ie,m,mp,3) = im(ie,m,mp,3) -  pi_const * weight *&
                     !                   conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                     !   im(ie,m,mp,4) = im(ie,m,mp,4) -  pi_const * weight *&
                     !                   conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                     !   im(ie,m,mp,5) = im(ie,m,mp,5) -  pi_const * weight *&
                     !                   conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                     !END IF
                     !
                     ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
                     !
                     DO ilo = 1, atoms%nlo(nType)
                        IF(atoms%llo(ilo,nType).EQ.l) THEN
                           g(iz,m,mp,1) = g(iz,m,mp,1) + weight * (  usdus%uulon(ilo,nType,ispin) * (&
                                          conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                          conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) )&
                                          + usdus%dulon(ilo,nType,ispin) * (&
                                          conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                          conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)))
                           DO ilop = 1, atoms%nlo(nType)
                              IF (atoms%llo(ilop,nType).EQ.l) THEN
                                 g(iz,m,mp,1) = g(iz,m,mp,1) + weight * usdus%uloulopn(ilo,ilop,nType,ispin) *&
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
            DO iz = 1, greensf%nz 
               DO it = 1, sym%invarind(natom)
                  is = sym%invarop(natom,it)
                  isi = sym%invtab(is)
                  d_mat(:,:) = cmplx(0.0,0.0)
                  DO m = -l,l
                     DO mp = -l,l
                        d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
                     ENDDO
                  ENDDO
                  calc_mat = matmul( transpose( conjg(d_mat) ) , &
                              g(iz,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,imat))
                  calc_mat =  matmul( calc_mat, d_mat )
                  DO m = -l,l
                     DO mp = -l,l
                        IF(imat.EQ.1) THEN
                           greensf%gmmpmat(iz,i_gf,m,mp,ispin,1) = greensf%gmmpmat(iz,i_gf,m,mp,ispin,1) + fac * calc_mat(m,mp)
                           greensf%gmmpmat(iz,i_gf,m,mp,ispin,2) = greensf%gmmpmat(iz,i_gf,m,mp,ispin,2) + fac * conjg(calc_mat(m,mp))
                        !ELSE IF(imat.EQ.2) THEN
                        !   greensfCoeffs%uu(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) + fac * conjg(calc_mat(m,mp))
                        !ELSE IF(imat.EQ.3) THEN
                        !   greensfCoeffs%dd(ie,i_gf,m,mp,ispin) = greensfCoeffs%dd(ie,i_gf,m,mp,ispin) + fac * conjg(calc_mat(m,mp))
                        !ELSE IF(imat.EQ.4) THEN
                        !   greensfCoeffs%ud(ie,i_gf,m,mp,ispin) = greensfCoeffs%ud(ie,i_gf,m,mp,ispin) + fac * conjg(calc_mat(m,mp))
                        !ELSE IF(imat.EQ.5) THEN
                        !   greensfCoeffs%du(ie,i_gf,m,mp,ispin) = greensfCoeffs%du(ie,i_gf,m,mp,ispin) + fac * conjg(calc_mat(m,mp))
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO!it
            ENDDO!ie
         ENDDO!imat
      ENDDO!natom
   ENDDO !i_gf

END SUBROUTINE greensfRes



END MODULE m_greensfRes