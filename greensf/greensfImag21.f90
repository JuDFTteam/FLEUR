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

   CONTAINS

   SUBROUTINE greensfImag21(atoms,sym,angle,input,nbands,dosWeights,resWeights,ind,wtkpt,eig,denCoeffsOffDiag,eigVecCoeffs,greensfCoeffs)


      IMPLICIT NONE

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_sym),               INTENT(IN)     :: sym
      REAL,                      INTENT(IN)     :: angle(sym%nop)
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_greensfCoeffs),     INTENT(INOUT)  :: greensfCoeffs

      INTEGER,                   INTENT(IN)     :: nbands
      REAL,                      INTENT(IN)     :: wtkpt
      REAL,                      INTENT(IN)     :: dosWeights(greensfCoeffs%ne,nbands)
      REAL,                      INTENT(IN)     :: resWeights(greensfCoeffs%ne,nbands)
      INTEGER,                   INTENT(IN)     :: ind(nbands,2)
      REAL,                      INTENT(IN)     :: eig(nbands)

      INTEGER  i_gf,nType,l,natom,ib,j,ie,m,lm,mp,lmp,imat,it,is,isi,ilo,ilop,nn
      REAL     fac
      COMPLEX phase,weight
      LOGICAL  l_zero,l_tria
      COMPLEX, ALLOCATABLE :: im(:,:,:,:)
      COMPLEX d_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),calc_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      IF(.NOT.input%l_gfsphavg) CALL juDFT_error("NOCO-offdiagonal + Radial dependence of onsite-GF not implemented",calledby="onsite21")

      l_tria = (input%tria.OR.input%gfTet).AND..NOT.input%l_hist
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(wtkpt,nbands,l_tria) &
      !$OMP SHARED(atoms,input,eigVecCoeffs,greensfCoeffs,denCoeffsOffDiag,eig) &
      !$OMP SHARED(dosWeights,resWeights,ind) &
      !$OMP PRIVATE(i_gf,imat,l,nType,natom,nn,ie,m,mp,lm,lmp,weight,ib,j,l_zero,ilo,ilop) &
      !$OMP PRIVATE(im)
      !$OMP DO
      DO i_gf = 1, atoms%n_gf
         nType = atoms%gfelem(i_gf)%atomType
         l = atoms%gfelem(i_gf)%l

         ALLOCATE(im(greensfCoeffs%ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(1,5,input%l_gfsphavg)))

         DO nn = 1, atoms%neq(nType)
            natom = SUM(atoms%neq(:nType-1)) + nn

            im = 0.0
            DO m = -l, l
               lm = l*(l+1) + m
               DO mp = -l, l
                  lmp = l*(l+1) + mp
                  !Loop through bands
                  DO ib = 1, nbands

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

                     DO ie = MERGE(ind(ib,1),j,l_tria), MERGE(ind(ib,2),j,l_tria)

                         weight = 2.0/input%jspins*(MERGE(resWeights(ie,ib),0.0,input%l_resolvent)&
                                                - ImagUnit * pi_const * MERGE(dosWeights(ie,ib),wtkpt/greensfCoeffs%del,l_tria))
                        !
                        !Contribution from states
                        !
                        im(ie,m,mp,1) = im(ie,m,mp,1) + weight *&
                                             (CONJG(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%uu21n(l,nType)&
                                            + CONJG(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%dd21n(l,nType)&
                                            + CONJG(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%ud21n(l,nType)&
                                            + CONJG(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%du21n(l,nType))
                        IF(.NOT.input%l_gfsphavg) THEN
                           im(ie,m,mp,2) = im(ie,m,mp,2) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%uu21n(l,nType)
                           im(ie,m,mp,3) = im(ie,m,mp,3) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%dd21n(l,nType)
                           im(ie,m,mp,4) = im(ie,m,mp,4) + weight * conjg(eigVecCoeffs%acof(ib,lmp,natom,2)) * eigVecCoeffs%bcof(ib,lm,natom,1) * denCoeffsOffdiag%ud21n(l,nType)
                           im(ie,m,mp,5) = im(ie,m,mp,5) + weight * conjg(eigVecCoeffs%bcof(ib,lmp,natom,2)) * eigVecCoeffs%acof(ib,lm,natom,1) * denCoeffsOffdiag%du21n(l,nType)
                        ENDIF
                        !
                        !Contribution from local Orbitals
                        !
                        DO ilo = 1, atoms%nlo(nType)
                           IF(atoms%llo(ilo,nType).EQ.l) THEN
                              !TODO: Something is wrong here for noco%l_soc and noco%l_mperp with a local orbital on the same l
                              im(ie,m,mp,1) = im(ie,m,mp,1) + weight *&
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
                                    im(ie,m,mp,1) = im(ie,m,mp,1) + (weight  * denCoeffsOffDiag%uloulop21n(ilo,ilop,nType) *&
                                    conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,2)) *eigVecCoeffs%ccof(m,ib,ilo,natom,1))
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO!local orbitals
                     ENDDO!ie
                  ENDDO!ib
                  DO ie = 1, greensfCoeffs%ne
                     greensfCoeffs%projdos21(ie,m,mp,nn,i_gf) = greensfCoeffs%projdos21(ie,m,mp,nn,i_gf) - AIMAG(im(ie,m,mp,1))
                  ENDDO
               ENDDO!mp
            ENDDO!m

            !Rotate the eqivalent atom into the irreducible brillouin zone
            !fac = 1.0/(sym%invarind(natom)*atoms%neq(nType))
            !IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No symmetry operations",calledby="greensfImag")
            !DO imat = 1, MERGE(1,5,input%l_gfsphavg)
            !   DO ie = 1, greensfCoeffs%ne
            !      DO it = 1, sym%invarind(natom)
            !         is = sym%invarop(natom,it)
            !         isi = sym%invtab(is)
            !         d_mat(:,:) = cmplx(0.0,0.0)
            !         DO m = -l,l
            !            DO mp = -l,l
            !               d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
            !            ENDDO
            !         ENDDO
            !         calc_mat = matmul( transpose( conjg(d_mat) ) , &
            !                     im(ie,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,imat))
            !         calc_mat =  matmul( calc_mat, d_mat )
            !         !phase = exp(ImagUnit * angle(isi))
            !         DO m = -l,l
            !            DO mp = -l,l
            !               IF(imat.EQ.1) THEN

                           !ELSE IF(imat.EQ.2) THEN
                           !   greensfCoeffs%uu(ie,i_gf,m,mp,3) = greensfCoeffs%uu(ie,i_gf,m,mp,3) - AIMAG(fac * conjg(calc_mat(m,mp)))
                           !ELSE IF(imat.EQ.3) THEN
                           !   greensfCoeffs%dd(ie,i_gf,m,mp,3) = greensfCoeffs%dd(ie,i_gf,m,mp,3) - AIMAG(fac * conjg(calc_mat(m,mp)))
                           !ELSE IF(imat.EQ.4) THEN
                           !   greensfCoeffs%ud(ie,i_gf,m,mp,3) = greensfCoeffs%ud(ie,i_gf,m,mp,3) - AIMAG(fac * conjg(calc_mat(m,mp)))
                           !ELSE IF(imat.EQ.5) THEN
                           !   greensfCoeffs%du(ie,i_gf,m,mp,3) = greensfCoeffs%du(ie,i_gf,m,mp,3) - AIMAG(fac * conjg(calc_mat(m,mp)))
                           !ENDIF
                        !ENDDO
                     !ENDDO
                  !ENDDO!it
               !ENDDO!ie
            !ENDDO!imat

         ENDDO !natom
      ENDDO !i_gf
      !$OMP END DO
      !$OMP END PARALLEL

   END SUBROUTINE greensfImag21

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

END MODULE m_greensfImag21