MODULE m_greensfImag

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfImag(atoms,gfinp,sym,input,ispin,nbands,dosWeights,resWeights,ind,wtkpt,ef,eig,usdus,eigVecCoeffs,greensfCoeffs)

      !This Subroutine calculates the contribution to the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
      !of the current k-Point (it is called in cdnval) inside the MT-sphere
      !and sums over the Brillouin-Zone using the histogram method or linear tetrahedron method
      !It is essentially the l-density of states in a (m,mp) matrix with an additional factor - pi

      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_gfinp),         INTENT(IN)    :: gfinp
      TYPE(t_sym),           INTENT(IN)    :: sym
      TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
      TYPE(t_usdus),         INTENT(IN)    :: usdus
      TYPE(t_greensfCoeffs), INTENT(INOUT) :: greensfCoeffs
      TYPE(t_input),         INTENT(IN)    :: input
      INTEGER,               INTENT(IN)    :: ispin  !Current spin index
      INTEGER,               INTENT(IN)    :: nbands !Number of bands to be considered
      REAL,                  INTENT(IN)    :: wtkpt  !Weight of the current k-point (not used in tetrahedron method)
      REAL,                  INTENT(IN)    :: ef
      REAL,    ALLOCATABLE,  INTENT(IN)    :: resWeights(:,:)
      REAL,    ALLOCATABLE,  INTENT(IN)    :: dosWeights(:,:) !Precalculated tetrahedron weights for the current k-point
      INTEGER, ALLOCATABLE,  INTENT(IN)    :: ind(:,:)        !Gives the range where the tetrahedron weights are non-zero
      REAL,                  INTENT(IN)    :: eig(:)          !Eigenvalues for the current k-point

      LOGICAL :: l_zero,l_tria
      INTEGER :: i_gf,ib,ie,j,nType,nn,natom
      INTEGER :: l,m,mp,lm,lmp,ilo,ilop
      INTEGER :: ie_start,ie_end
      REAL    :: del,eb
      COMPLEX :: weight
      COMPLEX, ALLOCATABLE :: im(:,:)

      !Temporary until input%tria/input%gauss are sorted out
      !l_tria = (input%tria.OR.input%gfTet).AND..NOT.input%l_hist
      l_tria=.false.

      IF(l_tria.AND.ALLOCATED(ind)) THEN
         IF(ANY(ind.GT.gfinp%ne).OR.ANY(ind.LT.1)) THEN
            CALL juDFT_error("Invalid index",calledby="greensfImag")
         ENDIF
      ENDIF

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del_out=del,eb_out=eb)

      !Loop through the gf elements to be calculated

      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(ispin,wtkpt,nbands,l_tria,del,eb) &
      !$OMP SHARED(atoms,gfinp,input,eigVecCoeffs,usdus,greensfCoeffs,eig,sym) &
      !$OMP SHARED(dosWeights,resWeights,ind) &
      !$OMP PRIVATE(i_gf,natom,l,nType,ie,m,mp,lm,lmp,ilo,ilop,weight,ib,j,l_zero,ie_start,ie_end) &
      !$OMP PRIVATE(im)
      !$OMP DO
      DO i_gf = 1, gfinp%n

         l     = gfinp%elem(i_gf)%l
         nType = gfinp%elem(i_gf)%atomType

         ALLOCATE(im(gfinp%ne,MERGE(1,5,gfinp%l_sphavg)),source=cmplx_0)

         !Loop through equivalent atoms
         DO nn = 1, atoms%neq(nType)
            natom = SUM(atoms%neq(:nType-1)) + nn
            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp
                  im = cmplx_0
                  !Loop through bands
                  DO ib = 1, nbands
                     !Check wether there is a non-zero weight for the energy window
                     l_zero = .true.
                     IF(l_tria) THEN
                        !IF(.NOT.input%l_resolvent) THEN
                           !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                           IF(ANY(dosWeights(ind(ib,1):ind(ib,2),ib).NE.0.0)) l_zero = .false.
                        !ELSE
                        !   l_zero = .false.
                        !ENDIF
                     ELSE
                        !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
                        j = FLOOR((eig(ib)-eb)/del)+1
                        IF((j.LE.gfinp%ne).AND.(j.GE.1)) l_zero = .false.
                     END IF

                     IF(l_zero) CYCLE

                     !Choose the relevant energy points depending on the bz-integration method
                     IF(l_tria) THEN
                        ie_start = ind(ib,1)
                        ie_end = ind(ib,2)
                     ELSE
                        ie_start = j
                        ie_end = j
                     ENDIF
                     DO ie = ie_start, ie_end
                        !weight for the bz-integration including spin-degeneracy
                        IF(l_tria) THEN
                           weight = -2.0/input%jspins * ImagUnit * pi_const * dosWeights(ie,ib)!+resWeights(ie,ib)
                        ELSE
                           weight = -2.0/input%jspins * ImagUnit * pi_const * wtkpt/del
                        ENDIF
                        !-------------------------
                        !Contribution from states
                        !-------------------------
                        im(ie,1) = im(ie,1) + weight *&
                                             (conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                             conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                             usdus%ddn(l,nType,ispin))
                        IF(.NOT.gfinp%l_sphavg) THEN
                           im(ie,2) = im(ie,2) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin)
                           im(ie,3) = im(ie,3) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)
                           im(ie,4) = im(ie,4) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)
                           im(ie,5) = im(ie,5) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin)
                        END IF
                        !------------------------------------------------------------------------------------------------------
                        ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
                        !------------------------------------------------------------------------------------------------------
                        DO ilo = 1, atoms%nlo(nType)
                           IF(atoms%llo(ilo,nType).EQ.l) THEN
                              im(ie,1) = im(ie,1) + weight * (  usdus%uulon(ilo,nType,ispin) * (&
                                             conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                             conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) )&
                                             + usdus%dulon(ilo,nType,ispin) * (&
                                             conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                             conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)))
                              DO ilop = 1, atoms%nlo(nType)
                                 IF (atoms%llo(ilop,nType).EQ.l) THEN
                                    im(ie,1) = im(ie,1) + weight * usdus%uloulopn(ilo,ilop,nType,ispin) *&
                                                   conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,ib,ilo,natom,ispin)
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO!ie
                  ENDDO!ib
                  DO ie = 1, gfinp%ne
                     greensfCoeffs%projdos(ie,m,mp,nn,i_gf,ispin) = greensfCoeffs%projdos(ie,m,mp,nn,i_gf,ispin) + im(ie,1)
                     IF(.NOT.gfinp%l_sphavg) THEN
                        greensfCoeffs%uu(ie,m,mp,nn,i_gf,ispin) = greensfCoeffs%uu(ie,m,mp,nn,i_gf,ispin) + im(ie,2)
                        greensfCoeffs%dd(ie,m,mp,nn,i_gf,ispin) = greensfCoeffs%dd(ie,m,mp,nn,i_gf,ispin) + im(ie,3)
                        greensfCoeffs%ud(ie,m,mp,nn,i_gf,ispin) = greensfCoeffs%ud(ie,m,mp,nn,i_gf,ispin) + im(ie,4)
                        greensfCoeffs%du(ie,m,mp,nn,i_gf,ispin) = greensfCoeffs%du(ie,m,mp,nn,i_gf,ispin) + im(ie,5)
                     ENDIF
                  ENDDO
               ENDDO !mp
            ENDDO !m
         ENDDO!nn
         DEALLOCATE(im)
      ENDDO !i_gf
      !$OMP END DO
      !$OMP END PARALLEL

   END SUBROUTINE greensfImag
END MODULE m_greensfImag