!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_forcea21soc
   !! Spin-orbit part of the A21 Pulay force for SOC in the first variation (noco).
   !!
   !! The MT SOC matrix of an atom depends on its position only through the phases of
   !! the matching coefficients, dH_SO(G,G')/dR = i(G'-G) H_SO(G,G'). Hence
   !!    F = 2 sum_n w_n Im[ cof_n^* H_SO vcof_n ],
   !! with the MT SOC matrix in the (lm, radial slot, local spin) basis
   !!    H_SO(lm p j1, lm' p' j2) = rso(p,p',l,j1,j2) <l m|L|l m'>.sigma(j1,j2),
   !! i.e. the matrix built in plane-wave form by hsmt_soc_offdiag(_LO).
   USE m_constants, ONLY: tpi_const, ImagUnit
   USE m_juDFT
   USE m_types_atoms
   USE m_types_cell
   USE m_types_sym
   USE m_types_nococonv
   USE m_types_rsoc
   USE m_types_abc
   USE m_types_force
   USE m_types_misc
   USE m_hsmt_spinor
   USE m_forcea21
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: force_a21_soc
CONTAINS
   SUBROUTINE force_a21_soc(atoms,sym,cell,nococonv,rsoc,itype,ne,we,abcUp,abcDn,vecUp,vecDn,f_a21_soc,results)
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_rsoc),     INTENT(IN)    :: rsoc
      INTEGER,          INTENT(IN)    :: itype, ne
      REAL,             INTENT(IN)    :: we(:)
      TYPE(t_abc),      INTENT(IN)    :: abcUp, abcDn   !! coefficients of the local spin components
      TYPE(t_force),    INTENT(IN)    :: vecUp, vecDn   !! their vector coefficients (aveccof, bveccof, cveccof)
      COMPLEX,          INTENT(INOUT) :: f_a21_soc(:,:)
      TYPE(t_results),  INTENT(INOUT) :: results

      INTEGER :: l, mr, mc, lmr, lmc, p1, p2, j1, j2, ie, iat, na, nr, i
      REAL    :: rso
      COMPLEX :: acc(3), t(3), forc(3)
      COMPLEX :: sig(2,2,3)
      COMPLEX, ALLOCATABLE :: lam(:,:,:,:)
      COMPLEX, ALLOCATABLE :: co(:,:,:,:,:), vc(:,:,:,:,:,:)
      REAL,    ALLOCATABLE :: a21(:,:)

      IF (.NOT.atoms%l_geo(itype)) RETURN
      CALL timestart("force_a21_soc")

      nr = MAXVAL(abcUp%n_r(1:atoms%lmax(itype)))
      ALLOCATE(co(ne,0:atoms%lmax(itype)*(atoms%lmax(itype)+2),nr,atoms%neq(itype),2))
      ALLOCATE(vc(3,ne,0:atoms%lmax(itype)*(atoms%lmax(itype)+2),nr,atoms%neq(itype),2),source=CMPLX(0.0,0.0))
      ALLOCATE(lam(-atoms%lmax(itype):atoms%lmax(itype),-atoms%lmax(itype):atoms%lmax(itype),2,2))
      ALLOCATE(a21(3,atoms%nat),source=0.0)

      co(:,:,:,:,1) = abcUp%cof(:ne,:,:nr,:)
      co(:,:,:,:,2) = abcDn%cof(:ne,:,:nr,:)
      DO iat = 1, atoms%neq(itype)
         na = atoms%firstAtom(itype) + iat - 1
         CALL fill_vec(vecUp, na, vc(:,:,:,:,iat,1))
         CALL fill_vec(vecDn, na, vc(:,:,:,:,iat,2))
      END DO

      DO iat = 1, atoms%neq(itype)
         na = atoms%firstAtom(itype) + iat - 1
         CALL spin_ops(na, sig)
         DO l = 1, atoms%lmax(itype)
            CALL l_dot_sigma(l, sig, lam)
            DO j2 = 1, 2
               DO j1 = 1, 2
                  DO p2 = 1, abcUp%n_r(l)
                     DO p1 = 1, abcUp%n_r(l)
                        rso = rsoc%rso(p1,p2,itype,l,j1,j2)
                        IF (ABS(rso) < 1.0e-14) CYCLE
                        DO ie = 1, ne
                           acc = CMPLX(0.0,0.0)
                           DO mr = -l, l
                              lmr = l*(l+1) + mr
                              t = CMPLX(0.0,0.0)
                              DO mc = -l, l
                                 lmc = l*(l+1) + mc
                                 t(:) = t(:) + lam(mr,mc,j1,j2)*vc(:,ie,lmc,p2,iat,j2)
                              END DO
                              acc(:) = acc(:) + CONJG(co(ie,lmr,p1,iat,j1))*t(:)
                           END DO
                           a21(:,na) = a21(:,na) + 2.0*rso*AIMAG(acc(:))*we(ie)/atoms%neq(itype)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL force_a21_symmetrize(atoms,sym,cell,itype,a21,forc)
      DO i = 1, 3
         results%force(i,itype,1) = results%force(i,itype,1) + REAL(forc(i))
         f_a21_soc(i,itype) = f_a21_soc(i,itype) + forc(i)
      END DO

      CALL timestop("force_a21_soc")

   CONTAINS

      SUBROUTINE fill_vec(vec, na, v)
         !! Vector coefficients of atom na in the radial-slot order of abc%cof
         TYPE(t_force), INTENT(IN)    :: vec
         INTEGER,       INTENT(IN)    :: na
         COMPLEX,       INTENT(INOUT) :: v(:,:,0:,:)
         INTEGER :: lo, l, m, lm
         v(:,:,:,1) = vec%aveccof(:,:ne,0:UBOUND(v,3),na)
         IF (nr > 1) v(:,:,:,2) = vec%bveccof(:,:ne,0:UBOUND(v,3),na)
         DO lo = 1, atoms%nlo(itype)
            IF (atoms%l_dulo(lo,itype)) CYCLE ! APW LO: its udot part is in bveccof
            l = atoms%llo(lo,itype)
            IF (l == 0) CYCLE ! no SOC for s states
            DO m = -l, l
               lm = l*(l+1) + m
               v(:,:,lm,atoms%slot_of_lo(lo,itype)) = vec%cveccof(:,m,:ne,lo,na)
            END DO
         END DO
      END SUBROUTINE fill_vec

      SUBROUTINE spin_ops(na, sig)
         !! sigma_a in the local spin frame, expressed in the rotated orbital frame
         !! (r' = R r, see hsmt_ab) in which the coefficients of atom na are given;
         !! L and sigma are axial vectors, so sigma'_a = det(R) sum_b R_ab sigma_b.
         INTEGER, INTENT(IN)  :: na
         COMPLEX, INTENT(OUT) :: sig(2,2,3)
         COMPLEX :: chi(2,2,2,2), isigma(2,2,3)
         REAL    :: rot(3,3), det
         INTEGER :: np, a, b, mr3(3,3)

         CALL hsmt_spinor_soc(itype, nococonv, chi, isigma)
         np = sym%invtab(sym%ngopr(na))
         mr3 = sym%mrot(:,:,np)
         rot = MATMUL(TRANSPOSE(cell%bmat), MATMUL(TRANSPOSE(REAL(mr3)), TRANSPOSE(cell%amat)))/tpi_const
         det = REAL(mr3(1,1)*(mr3(2,2)*mr3(3,3)-mr3(2,3)*mr3(3,2)) &
                   -mr3(1,2)*(mr3(2,1)*mr3(3,3)-mr3(2,3)*mr3(3,1)) &
                   +mr3(1,3)*(mr3(2,1)*mr3(3,2)-mr3(2,2)*mr3(3,1)))
         sig = CMPLX(0.0,0.0)
         DO a = 1, 3
            DO b = 1, 3
               sig(:,:,a) = sig(:,:,a) - ImagUnit*det*rot(a,b)*isigma(:,:,b)
            END DO
         END DO
      END SUBROUTINE spin_ops

      SUBROUTINE l_dot_sigma(l, sig, lam)
         !! lam(mr,mc,j1,j2) = sum_a <l mr|L_a|l mc> sigma_a(j1,j2), Ylm convention of ylm4
         INTEGER, INTENT(IN)  :: l
         COMPLEX, INTENT(IN)  :: sig(2,2,3)
         COMPLEX, INTENT(OUT) :: lam(-atoms%lmax(itype):,-atoms%lmax(itype):,:,:)
         COMPLEX :: lel(3)
         REAL    :: lplus, lminus
         INTEGER :: mr, mc, j1, j2

         DO mc = -l, l
            DO mr = -l, l
               lplus = 0.0 ; lminus = 0.0
               IF (mr == mc+1) lplus  = SQRT(REAL((l-mc)*(l+mc+1)))
               IF (mr == mc-1) lminus = SQRT(REAL((l+mc)*(l-mc+1)))
               lel(1) = CMPLX(0.5*(lplus+lminus), 0.0)
               lel(2) = CMPLX(0.0, -0.5*(lplus-lminus))
               lel(3) = CMPLX(MERGE(REAL(mc), 0.0, mr == mc), 0.0)
               DO j2 = 1, 2
                  DO j1 = 1, 2
                     lam(mr,mc,j1,j2) = SUM(lel(:)*sig(j1,j2,:))
                  END DO
               END DO
            END DO
         END DO
      END SUBROUTINE l_dot_sigma

   END SUBROUTINE force_a21_soc
END MODULE m_forcea21soc
