MODULE m_abcof3
CONTAINS
  SUBROUTINE abcof3(input,atoms,sym,jspin, cell, bkpt,lapw,&
       usdus, kveclo,oneD,a,b,bascof_lo)
    !     ************************************************************
    !     subroutine constructs the a,b coefficients of the linearized
    !     m.t. wavefunctions for each band and atom.       c.l. fu
    !     ************************************************************
#include "cpp_double.h"

    USE m_constants, ONLY : tpi_const
    USE m_setabc1locdn1
    USE m_sphbes
    USE m_dsphbs
    USE m_abclocdn1
    USE m_ylm
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_usdus),INTENT(IN)   :: usdus
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)  :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin 

    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: kveclo(atoms%nlotot)
    REAL,    INTENT (IN) :: bkpt(3)
    COMPLEX, INTENT (OUT):: a(:,0:,:)!(dimension%nvd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (OUT):: b(:,0:,:)!(dimension%nvd,0:dimension%lmd,atoms%natd)
    COMPLEX, INTENT (OUT):: bascof_lo(3,-atoms%llod:atoms%llod,4*atoms%llod+2,atoms%nlod,atoms%natd)
    !     .. Local Scalars ..
    COMPLEX phase,c_0,c_1,c_2,ci
    REAL const,df,r1,s,tmk,wronk
    INTEGER i,j,k,l,ll1,lm ,n,nap,natom,nn,iatom,jatom,lmp,mp
    INTEGER inv_f,ilo,nvmax,lo,n_ldau,inap,iintsp
    INTEGER nk_lo_sv,nk_lo,m
    !     ..
    !     .. Local Arrays ..
    INTEGER kvec(2*(2*atoms%llod+1),atoms%nlod,atoms%natd  )
    INTEGER nbasf0(atoms%nlod,atoms%natd),nkvec(atoms%nlod,atoms%natd)
    REAL dfj(0:atoms%lmaxd),fj(0:atoms%lmaxd),fk(3),fkp(3),fkr(3)
    REAL alo1(atoms%nlod,atoms%ntypd),blo1(atoms%nlod,atoms%ntypd),clo1(atoms%nlod,atoms%ntypd)
    COMPLEX ylm( (atoms%lmaxd+1)**2 )
    LOGICAL enough(atoms%natd),apw(0:atoms%lmaxd,atoms%ntypd)


    !     
    const = 2 * tpi_const/sqrt(cell%omtil)
    !
    a         = cmplx(0.0,0.0)
    b         = cmplx(0.0,0.0)
    bascof_lo = cmplx(0.0,0.0)
    !+APW_LO
    DO n = 1, atoms%ntype
       DO l = 0,atoms%lmax(n)
          apw(l,n) = .false.
          DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(l,n) = .true.
          ENDDO
          IF ((input%l_useapw).AND.(atoms%lapw_l(n).GE.l)) apw(l,n) = .false.

       ENDDO
       DO lo = 1,atoms%nlo(n)
          IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n),n) = .true.
       ENDDO
    ENDDO
    !+APW_LO
    !
    iintsp = 1
 
    CALL setabc1locdn1(jspin, atoms,lapw, sym,usdus,kveclo,enough,nkvec,kvec,&
         nbasf0,alo1,blo1,clo1)


    !---> loop over lapws
    DO  k = 1,nvmax
       !calculate k+G
       fk(1) = bkpt(1) + lapw%k1(k,jspin)
       fk(2) = bkpt(2) + lapw%k2(k,jspin)
       fk(3) = bkpt(3) + lapw%k3(k,jspin)

       !dotirp(f,g,bbmat) calculates the scalar product of f,g in reciprocal space
       s=dot_product(fk,matmul(fk,cell%bbmat))
       s = sqrt(s) ! s=|k+G|

       !--->   loop over atom types
       natom = 0
       DO  n = 1,atoms%ntype
          !calculate R_mt(itype)*|k+G|
          r1 = atoms%rmt(n)*s

          !compute sph. bessel function at r1 up to order lmax(n) stored in fj(0:lmax(n))
          CALL sphbes(atoms%lmax(n),r1, fj)

          !compute derivative of sph. bessel function at r1 up to oder lmax(n) stored in dfj(0:lmax(n))
          CALL dsphbs(atoms%lmax(n),r1,fj, dfj)

          !   ----> construct a and b coefficients
          DO  l = 0,atoms%lmax(n)
             !calculate |k+G|*d/dx j_l(r1)
             df = s*dfj(l)

             wronk = usdus%uds(l,n,jspin)*usdus%dus(l,n,jspin)-usdus%us(l,n,jspin)*usdus%duds(l,n,jspin) !Wronski determinante
             IF (apw(l,n)) THEN
                fj(l) = 1.0*const * fj(l)/usdus%us(l,n,jspin)
                dfj(l) = 0.0d0
             ELSE
                dfj(l) = const* (usdus%dus(l,n,jspin)*fj(l)-df*usdus%us(l,n,jspin))/wronk
                fj(l) = const* (df*usdus%uds(l,n,jspin)-fj(l)*usdus%duds(l,n,jspin))/wronk
             ENDIF
          enddo
          !   ----> loop over equivalent atoms
          DO  nn = 1,atoms%neq(n)
             natom = natom + 1
             !invsat(natom) is 0 if atom natom can't be mapped via inversion symmetrie
             !              is 1 if atom natom can   be mapped via inversion symmetrie and is parent atom
             !              is 2 if atom natom can   be mapped via inversion symmetrie and is second atom

             IF ((atoms%invsat(natom).EQ.0) .OR. (atoms%invsat(natom).EQ.1)) THEN
                tmk = tpi_const* dot_product(fk(:),atoms%taual(:,natom))
                phase = cmplx(cos(tmk),sin(tmk))
                IF (oneD%odi%d1) THEN
                   inap = oneD%ods%ngopr(natom)
                   !                nap = ods%ngopr(natom)
                   !               inap = ods%invtab(nap)
                ELSE
                   nap = atoms%ngopr(natom)
                   inap = sym%invtab(nap)
                END IF
                DO  j = 1,3
                   fkr(j) = 0.
                   DO  i = 1,3
                      IF (oneD%odi%d1) THEN
                         fkr(j) = fkr(j) + fk(i)*oneD%ods%mrot(i,j,inap)
                      ELSE
                         fkr(j) = fkr(j) + fk(i)*sym%mrot(i,j,inap)
                      END IF
                   enddo
                enddo
                !transform fkr from reciprocal internal into reciprocal cartesian coordinates
                fkp=matmul(fkr,cell%bmat)
                !       ----> generate spherical harmonics at fkp up to order lmax(n) stored in ylm(1:(lmax+1)**2)
                CALL ylm4(atoms%lmax(n),fkp,ylm)
                !       ----> loop over l,m
                DO  l = 0,atoms%lmax(n)
                   ll1 = l* (l+1)
                   DO  m = -l,l
                      lm = ll1 + m
                      c_0 = conjg(ylm(lm+1))*phase
                      c_1 = c_0 *  fj(l)
                      c_2 = c_0 * dfj(l)

                      a(k,lm,natom) = c_1
                      b(k,lm,natom) = c_2

                   enddo
                enddo
                IF (.NOT.enough(natom)) THEN
                   CALL abclocdn1(atoms,sym, const,phase,ylm,n,natom,k,s,nvmax,&
                        nbasf0,alo1,blo1,clo1,kvec(1,1,natom), nkvec,enough,bascof_lo )

                ENDIF
             ENDIF    ! invsatom == ( 0 v 1 )
          enddo    ! loop over equivalent atoms
       enddo       ! loop over atom types
    enddo          ! loop over lapws


    iatom = 0
    DO n = 1,atoms%ntype
       DO nn = 1,atoms%neq(n)
          iatom = iatom + 1
          IF (atoms%invsat(iatom).EQ.1) THEN
             jatom = sym%invsatnr(iatom)
             DO ilo = 1,atoms%nlo(n)
                l = atoms%llo(ilo,n)
                DO m = -l,l
                   inv_f = (-1.0)**(m+l)
                   DO i = 1,3
                      bascof_lo(i,m,:,ilo,jatom) = inv_f * conjg(  bascof_lo(i,-m,:,ilo,iatom))
                   ENDDO
                ENDDO
             ENDDO
             DO l = 0,atoms%lmax(n)
                ll1 = l* (l+1)
                DO m =-l,l
                   lm  = ll1 + m
                   lmp = ll1 - m
                   inv_f = (-1.0)**(m+l)
                   DO k = 1,nvmax
                      a(k,lm,jatom) = inv_f *conjg(a(k,lmp,iatom))
                   ENDDO
                   DO k = 1,nvmax
                      b(k,lm,jatom) = inv_f *conjg(b(k,lmp,iatom))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE abcof3
END MODULE m_abcof3
