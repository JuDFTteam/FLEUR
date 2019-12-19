MODULE m_nmat21
!     ************************************************************
!     This subroutine calculates the density matrix n^{s}_{m,m'}
!     for a given atom 'n' and l-quantum number 'l'. The l's for
!     all atoms are stored in lda_u(), if lda_u()<0, no +U is used.
!     For details see Eq.(12) of Shick et al. PRB 60, 10765 (1999)
!     Part of the LDA+U package                   G.B., Oct. 2000
!     ************************************************************
   CONTAINS
   SUBROUTINE n_mat21(atoms,sym,angle,ne,we,denCoeffsOffdiag,eigVecCoeffs,n_mmp)

      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffdiag
      REAL,                      INTENT(IN)     :: angle(:)
      INTEGER,                   INTENT(IN)     :: ne
      REAL,                      INTENT(IN)     :: we(:)!(input%neig)
      COMPLEX,                   INTENT(INOUT)  :: n_mmp(-lmaxU_const:,-lmaxU_const:,:)

      INTEGER i,l,m,lp,mp,n,it,is,isi,natom
      INTEGER ilo,ilop,ll1,nn,lmp,lm,i_u,natomTemp
      REAL fac
      COMPLEX c_0,phase

      COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      !
      ! calculate n_mat:
      !

      natom = 0
      i_u = 1
      DO n = 1,atoms%ntype
         DO WHILE (i_u.LE.atoms%n_u)
            IF (atoms%lda_u(i_u)%atomType.GT.n) EXIT
            natomTemp = natom
            n_tmp(:,:) = cmplx_0
            l = atoms%lda_u(i_u)%l
            ll1 = (l+1)*l
            DO nn = 1, atoms%neq(n)
               natomTemp = natomTemp + 1
               !
               !  prepare n_mat in local frame (in noco-calculations this depends
               !                                also on alpha(n) and beta(n) )
               !
               DO m = -l,l
                  lm = ll1+m
                  DO mp = -l,l
                     lmp = ll1+mp
                     c_0 = cmplx_0
                     DO i = 1,ne
                        c_0 = c_0 +  we(i) * ( &
                                    conjg(eigVecCoeffs%acof(i,lmp,natomTemp,2))*eigVecCoeffs%acof(i,lm,natomTemp,1) * denCoeffsOffdiag%uu21n(l,n) &
                                  + conjg(eigVecCoeffs%acof(i,lmp,natomTemp,2))*eigVecCoeffs%bcof(i,lm,natomTemp,1) * denCoeffsOffdiag%ud21n(l,n) &
                                  + conjg(eigVecCoeffs%bcof(i,lmp,natomTemp,2))*eigVecCoeffs%acof(i,lm,natomTemp,1) * denCoeffsOffdiag%du21n(l,n) &
                                  + conjg(eigVecCoeffs%bcof(i,lmp,natomTemp,2))*eigVecCoeffs%bcof(i,lm,natomTemp,1) * denCoeffsOffdiag%dd21n(l,n))
                     ENDDO
                     n_tmp(m,mp) = c_0
                  ENDDO
               ENDDO
               !
               !  add local orbital contribution (if there is one) (untested so far)
               !
               DO ilo = 1, atoms%nlo(n)
                  IF (atoms%llo(ilo,n).EQ.l) THEN
                     DO m = -l,l
                        lm = ll1+m
                        DO mp = -l,l
                           lmp = ll1+mp
                           c_0 = cmplx_0
                           DO i = 1,ne
                              c_0 = c_0 +  we(i) * ( &
                                          conjg(eigVecCoeffs%acof(i,lmp,natomTemp,2))*eigVecCoeffs%ccof(m,i,ilo,natomTemp,1) * denCoeffsOffdiag%uulo21n(l,n) &
                                        + conjg(eigVecCoeffs%ccof(mp,i,ilo,natomTemp,2))*eigVecCoeffs%acof(i,lm,natomTemp,1) * denCoeffsOffdiag%ulou21n(l,n) &
                                        + conjg(eigVecCoeffs%bcof(i,lmp,natomTemp,2))*eigVecCoeffs%ccof(m,i,ilo,natomTemp,1) * denCoeffsOffdiag%dulo21n(l,n) &
                                        + conjg(eigVecCoeffs%ccof(mp,i,ilo,natomTemp,2))*eigVecCoeffs%bcof(i,lm,natomTemp,1) * denCoeffsOffdiag%ulod21n(l,n))
                           ENDDO
                           DO ilop = 1, atoms%nlo(n)
                              IF (atoms%llo(ilop,n).EQ.l) THEN
                                 DO i = 1,ne
                                    c_0 = c_0 +  we(i) * denCoeffsOffdiag%uloulop21n(ilo,ilop,n) *conjg(eigVecCoeffs%ccof(mp,i,ilop,natomTemp,2)) *eigVecCoeffs%ccof(m ,i,ilo ,natomTemp,1)
                                 ENDDO
                              ENDIF
                           ENDDO
                           n_tmp(m,mp) = n_tmp(m,mp) + c_0
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
               !
               !  n_mmp should be rotated by D_mm' ; compare force_a21
               !

               DO it = 1, sym%invarind(natomTemp)

                  fac = 1.0  /  ( sym%invarind(natomTemp) * atoms%neq(n) )
                  is = sym%invarop(natomTemp,it)
                  isi = sym%invtab(is)
                  d_tmp(:,:) = cmplx_0
                  DO m = -l,l
                     DO mp = -l,l
                        d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                     ENDDO
                  ENDDO
                  nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                  n1_tmp =  matmul( nr_tmp, d_tmp )
                  phase = exp(ImagUnit*angle(isi))
                  DO m = -l,l
                     DO mp = -l,l
                        n_mmp(m,mp,i_u) = n_mmp(m,mp,i_u) + conjg(n1_tmp(m,mp)) * fac * phase
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO ! sum  over equivalent
            i_u = i_u +1
         ENDDO
         natom = natom + atoms%neq(n)
      ENDDO     ! loop over atom types

   END SUBROUTINE n_mat21
END MODULE m_nmat21
