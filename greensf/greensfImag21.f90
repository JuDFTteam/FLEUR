MODULE m_greensfImag21
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
   !> In the non-collinear case spin isn't a good quantum number anymore so we get a 2x2 blockmatrix
   !
   ! REVISION HISTORY:
   ! 14 March 2019 - Initial Version
   !------------------------------------------------------------------------------

   USE m_types
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfImag21(atoms,gfinp,sym,input,nbands,dosWeights,resWeights,ind,wtkpt,ef,eig,denCoeffsOffDiag,eigVecCoeffs,greensfCoeffs)

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs
      INTEGER,                   INTENT(IN)     :: nbands
      REAL,                      INTENT(IN)     :: wtkpt
      REAL,                      INTENT(IN)     :: ef
      REAL,                      INTENT(IN)     :: dosWeights(:,:)
      REAL,                      INTENT(IN)     :: resWeights(:,:)
      INTEGER,                   INTENT(IN)     :: ind(:,:)
      REAL,                      INTENT(IN)     :: eig(:)

      INTEGER  i_gf,nType,l,natom,ib,j
      INTEGER  ie,m,lm,mp,lmp,ilo,ilop,nn
      REAL del,eb
      COMPLEX  weight
      LOGICAL  l_zero,l_tria
      COMPLEX, ALLOCATABLE :: im(:,:)

      IF(.NOT.gfinp%l_sphavg) CALL juDFT_error("NOCO-offdiagonal + Radial dependence of onsite-GF not implemented",calledby="onsite21")

      !Temporary until input%tria/input%gauss are sorted out
      !l_tria = (input%tria.OR.input%gfTet).AND..NOT.input%l_hist
      l_tria=.false.

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb)

      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(wtkpt,nbands,l_tria,del,eb) &
      !$OMP SHARED(atoms,gfinp,input,eigVecCoeffs,greensfCoeffs,denCoeffsOffDiag,eig) &
      !$OMP SHARED(dosWeights,resWeights,ind) &
      !$OMP PRIVATE(i_gf,l,nType,natom,nn,ie,m,mp,lm,lmp,weight,ib,j,l_zero,ilo,ilop) &
      !$OMP PRIVATE(im)
      !$OMP DO
      DO i_gf = 1, gfinp%n
         nType = gfinp%elem(i_gf)%atomType
         l =     gfinp%elem(i_gf)%l

         ALLOCATE(im(gfinp%ne,MERGE(1,5,gfinp%l_sphavg)))

         DO nn = 1, atoms%neq(nType)
            natom = SUM(atoms%neq(:nType-1)) + nn

            DO m = -l, l
               lm = l*(l+1) + m
               DO mp = -l, l
                  lmp = l*(l+1) + mp
                  im = 0.0
                  !Loop through bands
                  DO ib = 1, nbands

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
                        IF( (j.LE.gfinp%ne).AND.(j.GE.1) ) l_zero = .false.
                     END IF

                     IF(l_zero) CYCLE

                     DO ie = MERGE(ind(ib,1),j,l_tria), MERGE(ind(ib,2),j,l_tria)

                        weight = -2.0/input%jspins*ImagUnit * pi_const * MERGE(dosWeights(ie,ib),wtkpt/del,l_tria)
                        !
                        !Contribution from states
                        !
                        im(ie,1) = im(ie,1) + weight *&
                                             (CONJG(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%uu21n(l,nType)&
                                            + CONJG(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%dd21n(l,nType)&
                                            + CONJG(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%ud21n(l,nType)&
                                            + CONJG(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%du21n(l,nType))
                        IF(.NOT.gfinp%l_sphavg) THEN
                           im(ie,2) = im(ie,2) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%uu21n(l,nType)
                           im(ie,3) = im(ie,3) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%dd21n(l,nType)
                           im(ie,4) = im(ie,4) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%ud21n(l,nType)
                           im(ie,5) = im(ie,5) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%du21n(l,nType)
                        ENDIF
                        !
                        !Contribution from local Orbitals
                        !
                        DO ilo = 1, atoms%nlo(nType)
                           IF(atoms%llo(ilo,nType).EQ.l) THEN
                              !TODO: Something is wrong here for noco%l_soc and noco%l_mperp with a local orbital on the same l
                              im(ie,1) = im(ie,1) + weight *&
                              (conjg(eigVecCoeffs%acof(ib,lmp,natom,2))*eigVecCoeffs%ccof(m,ib,ilo,natom,1) *&
                               denCoeffsOffDiag%uulo21n(ilo,nType) +&
                               conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,2))*eigVecCoeffs%acof(ib,lm,natom,1) *&
                               denCoeffsOffDiag%ulou21n(ilo,nType) +&
                               conjg(eigVecCoeffs%bcof(ib,lmp,natom,2))*eigVecCoeffs%ccof(m,ib,ilo,natom,1) *&
                               denCoeffsOffDiag%dulo21n(ilo,nType) +&
                               conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,2))*eigVecCoeffs%bcof(ib,lm,natom,1) *&
                               denCoeffsOffDiag%ulod21n(ilo,nType))
                              DO ilop = 1, atoms%nlo(nType)
                                 IF(atoms%llo(ilop,nType).EQ.l) THEN
                                    im(ie,1) = im(ie,1) + (weight  * denCoeffsOffDiag%uloulop21n(ilo,ilop,nType) *&
                                    conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,2)) *eigVecCoeffs%ccof(m,ib,ilo,natom,1))
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO!local orbitals
                     ENDDO!ie
                  ENDDO!ib
                  DO ie = 1, gfinp%ne
                     greensfCoeffs%projdos(ie,m,mp,nn,i_gf,3) = greensfCoeffs%projdos(ie,m,mp,nn,i_gf,3) + im(ie,1)
                  ENDDO
               ENDDO!mp
            ENDDO!m
         ENDDO !nn
         DEALLOCATE(im)
      ENDDO !i_gf
      !$OMP END DO
      !$OMP END PARALLEL

   END SUBROUTINE greensfImag21
END MODULE m_greensfImag21