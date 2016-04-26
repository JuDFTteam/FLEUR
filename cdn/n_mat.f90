MODULE m_nmat
  !     ************************************************************
  !     This subroutine calculates the density matrix n^{s}_{m,m'}
  !     for a given atom 'n' and l-quantum number 'l'. The l's for
  !     all atoms are stored in lda_u(), if lda_u()<0, no +U is used.
  !     For details see Eq.(12) of Shick et al. PRB 60, 10765 (1999)
  !     Part of the LDA+U package                   G.B., Oct. 2000
  !     ************************************************************
CONTAINS
  SUBROUTINE n_mat(atoms,sym, ne,usdus,jspin,we, acof,bcof,ccof, n_mmp)
    !

    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)   :: usdus
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ne,jspin 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(:)!(dimension%neigd)
    COMPLEX, INTENT (IN) :: acof(:,0:,:)!(nobd,0:atoms%lmaxd*(lmaxd+2) ,natd)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:)!(nobd,0:atoms%lmaxd*(lmaxd+2) ,natd)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%natd)
    COMPLEX, INTENT (INOUT) :: n_mmp(-3:3,-3:3,atoms%n_u)
    !     ..
    !     .. Local Scalars ..
    COMPLEX c_0
    INTEGER i,j,k,l ,mp,n,it,is,isi,natom,n_ldau,lp,m
    INTEGER ilo,ilop,ll1,nn,lmp,lm
    REAL fac
    !     ..
    !     .. Local Arrays ..
    COMPLEX n_tmp(-3:3,-3:3),nr_tmp(-3:3,-3:3),d_tmp(-3:3,-3:3)
    COMPLEX n1_tmp(-3:3,-3:3)
    !     ..
    !
    ! calculate n_mat:
    !
    n_ldau = 0
    natom = 0
    DO n = 1,atoms%ntype
       IF (atoms%lda_u(n)%l.GE.0) THEN
          n_ldau = n_ldau + 1
          n_tmp(:,:) =cmplx(0.0,0.0)
          l = atoms%lda_u(n)%l
          ll1 = (l+1)*l 
          DO nn = 1, atoms%neq(n)
             natom = natom + 1
             !
             !  prepare n_mat in local frame (in noco-calculations this depends 
             !                                also on alpha(n) and beta(n) )
             !
             DO m = -l,l
                lm = ll1+m
                DO mp = -l,l
                   lmp = ll1+mp
                   c_0 = cmplx(0.0,0.0)
                   DO i = 1,ne
                      c_0 = c_0 +  we(i) * ( usdus%ddn(l,n,jspin) *&
                           conjg(bcof(i,lmp,natom))*bcof(i,lm,natom) +&
                           conjg(acof(i,lmp,natom))*acof(i,lm,natom) )
                   ENDDO
                   n_tmp(m,mp) = c_0 
                ENDDO
             ENDDO
             !
             !  add local orbrbital contribution (if there is one) (untested so far)
             !
             DO ilo = 1, atoms%nlo(n)
                IF (atoms%llo(ilo,n).EQ.l) THEN

                   DO m = -l,l
                      lm = ll1+m
                      DO mp = -l,l
                         lmp = ll1+mp
                         c_0 = cmplx(0.0,0.0)
                         DO i = 1,ne
                            c_0 = c_0 +  we(i) * (  usdus%uulon(ilo,n,jspin) * (&
                                 conjg(acof(i,lmp,natom))*ccof(m,i,ilo,natom) +&
                                 conjg(ccof(mp,i,ilo,natom))*acof(i,lm,natom) )&
                                 + usdus%dulon(ilo,n,jspin) * (&
                                 conjg(bcof(i,lmp,natom))*ccof(m,i,ilo,natom) +&
                                 conjg(ccof(mp,i,ilo,natom))*bcof(i,lm,natom)))
                         ENDDO
                         DO ilop = 1, atoms%nlo(n)
                            IF (atoms%llo(ilop,n).EQ.l) THEN
                               DO i = 1,ne
                                  c_0 = c_0 +  we(i) * usdus%uloulopn(ilo,ilop,n,jspin) *&
                                       conjg(ccof(mp,i,ilop,natom)) *ccof(m ,i,ilo ,natom)
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
             DO it = 1, sym%invarind(natom)

                fac = 1.0  /  ( sym%invarind(natom) * atoms%neq(n) )
                is = sym%invarop(natom,it)
                isi = sym%invtab(is)
                d_tmp(:,:) = cmplx(0.0,0.0)
                DO m = -l,l
                   DO mp = -l,l
                      d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                   ENDDO
                ENDDO
                nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                n1_tmp =  matmul( nr_tmp, d_tmp )
                DO m = -l,l
                   DO mp = -l,l
                      n_mmp(m,mp,n_ldau) = n_mmp(m,mp,n_ldau) +conjg(n1_tmp(m,mp)) * fac
                   ENDDO
                ENDDO

             ENDDO

          ENDDO ! sum  over equivalent atoms
       ELSE
          natom = natom + atoms%neq(n)
       ENDIF
    ENDDO     ! loop over atom types

    !     do m=-l,l
    !      write(*,'(14f12.6)') (n_mmp(m,mp),mp=-l,l)
    !     enddo
    !
    RETURN
  END SUBROUTINE n_mat
END MODULE m_nmat
