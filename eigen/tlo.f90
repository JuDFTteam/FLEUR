!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_tlo
   USE m_juDFT
  !***********************************************************************
  !     sets up the extra t-matrix elements due to the local orbitals.
  !     only non=zero elements are calculated
  !
  !     p.kurz jul. 1996
  !***********************************************************************
CONTAINS
   SUBROUTINE tlo(atoms,sym,sphhar,iSpinPr,iSpin,jsp,ntyp,enpara,lh0,input,vr,&
       na,flo,f,g,usdus, tlmplm, one, l_dfpt, l_V1)
      ! Contruct the additional local Hamiltonian matrices
      ! t_{L'L}^{\mu} = \sum_{lh} \int dV u_{l',order'}^{\mu}(r)Y_{l'}^{m'*}(\Omega)
      !                           * V_{lh}(r)Y_{lh}(\Omega)
      !                           * u_{l,order}^{\mu}(r)Y_{l}^{m}(\Omega)
      !                           * i^{l-l'}
      !               + \int dV u_{l',order'}^{\mu}(r) H_{sph}
      !               *         u_{l,order}^{\mu}(r)Y_{l}^{m}(\Omega)
      ! of a real valued potential V(\bm{r}). The superindex L is defined as
      ! L := (l,m,order)
      ! with order = 0 refering to radial functions u and order = 1 denoting
      ! their energy derivatives (udot). This construction is not k-dependent
      ! and therefore executed only once each scf iteration.

      ! Abbreviations:
      ! tuulo:   t-matrix element of an LO and the APW radial fuction
      ! tdulo:   t-matrix element of an LO and the energy derivative of
      !          the APW radial fuction
      ! tulou:   t-matrix element of the APW radial fuction and an LO
      ! tulod:   t-matrix element of the APW radial fuction derivative and an LO
      ! tuloulo: t-matrix element of two LOs

      USE m_intgr, ONLY : intgr3
      USE m_gaunt, ONLY: gaunt1
      USE m_types
      USE m_constants

      IMPLICIT NONE

      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_usdus),  INTENT(IN)    :: usdus
      TYPE(t_tlmplm), INTENT(INOUT) :: tlmplm
      TYPE(t_enpara), INTENT(IN)    :: enpara

      ! Scalar Arguments
      INTEGER, INTENT (IN) :: iSpinPr,iSpin,jsp,ntyp ,lh0,na

      ! Array Arguments
      REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd)
      REAL,    INTENT (IN) :: f(:,:,0:,:),g(:,:,0:,:) !(atoms%jmtd,2,0:atoms%lmaxd,spins)
      REAL,    INTENT (IN) :: flo(:,:,:,:)!(atoms%jmtd,2,atoms%nlod,spins)

      COMPLEX, INTENT(IN) :: one

      LOGICAL, INTENT(IN) :: l_dfpt, l_V1

      ! Local Scalars
      COMPLEX :: cil
      INTEGER :: i,l,lh,lm ,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,lpmin,lpmin0,lpp ,mem,mp,mpp,m,lmx,mlo,mlolo,s
      INTEGER :: loplo_new, mlolo_new

      ! Local Arrays
      REAL :: x(atoms%jmtd),ulovulo(atoms%nlod*(atoms%nlod+1)/2,lh0:sphhar%nlhd)
      REAL :: uvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd),dvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd)

      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         DO lp = 0,atoms%lmax(ntyp)
            lmin = ABS(lp-l)
            lmx = lp + l
            DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
               lpp = sphhar%llh(lh,sym%ntypsy(na))
               IF ((MOD(l+lp+lpp,2).EQ.1).OR.(lpp.LT.lmin).OR.(lpp.GT.lmx)) THEN
                  uvulo(lo,lp,lh) = 0.0
                  dvulo(lo,lp,lh) = 0.0
               ELSE
                  DO i = 1,atoms%jri(ntyp)
                     x(i) = (f(i,1,lp,iSpinPr)*flo(i,1,lo,iSpin)+ f(i,2,lp,iSpinPr)*flo(i,2,lo,iSpin))*vr(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),uvulo(lo,lp,lh))

                  DO i = 1,atoms%jri(ntyp)
                     x(i) = (g(i,1,lp,iSpinPr)*flo(i,1,lo,iSpin)+ g(i,2,lp,iSpinPr)*flo(i,2,lo,iSpin))*vr(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),dvulo(lo,lp,lh))
               END IF
            END DO
         END DO
      END DO
      loplo = 0
      DO lop = 1,atoms%nlo(ntyp)
         lp = atoms%llo(lop,ntyp)
         DO lo = 1,lop
            l = atoms%llo(lo,ntyp)
            loplo = loplo + 1
            IF (loplo>SIZE(ulovulo,1))  CALL juDFT_error("loplo too large!!!" ,calledby ="tlo")
            DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
               lpp = sphhar%llh(lh,sym%ntypsy(na))
               lmin = ABS(lp - l)
               lmx = lp + l
               IF ((MOD(l+lp+lpp,2).EQ.1).OR.(lpp.LT.lmin).OR.(lpp.GT.lmx)) THEN
                  ulovulo(loplo,lh) = 0.0
               ELSE
                  DO i = 1,atoms%jri(ntyp)
                     x(i) = (flo(i,1,lop,iSpinPr)*flo(i,1,lo,iSpin)+flo(i,2,lop,iSpinPr)*flo(i,2,lo,iSpin))*vr(i,lh)
                  END DO
                  CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),ulovulo(loplo,lh))
               END IF
            END DO
         END DO
      END DO

      ! Generate the t-matrices. For optimal performance consider only those
      ! combinations of l,l',l'',m,m',m'' that satisfy the three conditions for
      ! non-zero Gaunt coefficients, i.e.
      !     |l - l''| <= l' <= l + l'' (triangular condition)
      !     m' = m + m''
      !     l + l' + l'' even

      ! Loop over the local orbitals
      mlo=SUM(atoms%nlo(:ntyp-1))
      s=tlmplm%h_loc2_nonsph(ntyp)
      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         DO m = -l,l
            lm=l*(l+1)+m
            ! Loop over the lattice harmonics
            DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
               lpp = sphhar%llh(lh,sym%ntypsy(na))
               lpmin0 = ABS(l-lpp)
               lpmax0 = l + lpp
               ! Check that lpmax is smaller than the max l of the
               ! wavefunction expansion at this atom
               lpmax = MIN(lpmax0,atoms%lnonsph(ntyp))
               ! Make sure that l + l'' + lpmax is even
               lpmax = lpmax - MOD(l+lpp+lpmax,2)
               DO mem = 1,sphhar%nmem(lh,sym%ntypsy(na))
                  mpp = sphhar%mlh(mem,lh,sym%ntypsy(na))
                  mp = m + mpp
                  lpmin = MAX(lpmin0,ABS(mp))
                  !- Make sure that l + l'' + lpmin is even
                  lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                  ! Loop over l'
                  DO lp = lpmin,lpmax,2
                     lmp = lp* (lp+1) + mp
                     cil = ImagUnit**(l-lp) * sphhar%clnu(mem,lh,sym%ntypsy(na)) &
                       & * gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                     tlmplm%tuulo(lmp,m,lo+mlo,iSpinPr,iSpin) = &
                     & tlmplm%tuulo(lmp,m,lo+mlo,iSpinPr,iSpin) + one * cil * uvulo(lo,lp,lh)
                     tlmplm%tdulo(lmp,m,lo+mlo,iSpinPr,iSpin) = &
                     & tlmplm%tdulo(lmp,m,lo+mlo,iSpinPr,iSpin) + one * cil * dvulo(lo,lp,lh)
                     tlmplm%h_LO(lmp,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lmp,m,lo+mlo,iSpinPr,iSpin) + one * cil * uvulo(lo,lp,lh)
                     tlmplm%h_LO(lmp+s,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lmp+s,m,lo+mlo,iSpinPr,iSpin) + one * cil * dvulo(lo,lp,lh)
                       tlmplm%tulou(lmp,m,lo+mlo,iSpinPr,iSpin) = &
                     & tlmplm%tulou(lmp,m,lo+mlo,iSpinPr,iSpin) + one * CONJG(cil*uvulo(lo,lp,lh))
                       tlmplm%tulod(lmp,m,lo+mlo,iSpinPr,iSpin) = &
                     & tlmplm%tulod(lmp,m,lo+mlo,iSpinPr,iSpin) + one * CONJG(cil*dvulo(lo,lp,lh))
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! Generate the t-matrix including two local orbitals for LO <= LO'
      ! Loop over LO'
      mlolo = DOT_PRODUCT(atoms%nlo(:ntyp-1),atoms%nlo(:ntyp-1)+1)/2
      mlolo_new = DOT_PRODUCT(atoms%nlo(:ntyp-1),atoms%nlo(:ntyp-1))
      DO lop = 1,atoms%nlo(ntyp)
         lp = atoms%llo(lop,ntyp)
         DO mp = -lp,lp
            ! Loop over the lattice harmonics
            DO lh = lh0,sphhar%nlh(sym%ntypsy(na))
               lpp = sphhar%llh(lh,sym%ntypsy(na))
               DO mem = 1,sphhar%nmem(lh,sym%ntypsy(na))
                  mpp = sphhar%mlh(mem,lh,sym%ntypsy(na))
                  m = mp - mpp
                  ! Loop over LO
                  DO lo = 1,lop
                     l = atoms%llo(lo,ntyp)
                     loplo = ((lop-1)*lop)/2 + lo
                     loplo_new = (lop-1) * atoms%nlo(ntyp) + lo
                     IF ((ABS(l-lpp).LE.lp).AND.(lp.LE.(l+lpp)).AND.(MOD(l+lp+lpp,2).EQ.0).AND.(ABS(m).LE.l)) THEN
                        cil = ImagUnit**(l-lp) * sphhar%clnu(mem,lh,sym%ntypsy(na)) &
                          & * gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                        tlmplm%tuloulo(mp,m,loplo+mlolo,iSpinPr,iSpin) = &
                      & tlmplm%tuloulo(mp,m,loplo+mlolo,iSpinPr,iSpin) + one * cil * ulovulo(loplo,lh)
                     !   tlmplm%tuloulo_new(mp,m,mlolo_new+loplo_new,iSpinPr,iSpin) = &
                     ! & tlmplm%tuloulo_new(mp,m,mlolo_new+loplo_new,iSpinPr,iSpin) + one * cil * ulovulo(loplo,lh)
                        tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,iSpinPr,iSpin) = &
                      & tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,iSpinPr,iSpin) + one * cil * ulovulo(loplo,lh)
                        IF (lop.NE.lo) THEN
                           !loplo = ((lo-1)*lo)/2 + lop
                           !loplo_new = (lo-1) * atoms%nlo(ntyp) + lop
                        !   tlmplm%tuloulo_new(m,mp,mlolo_new+loplo_new,iSpinPr,iSpin) = &
                        ! & tlmplm%tuloulo_new(m,mp,mlolo_new+loplo_new,iSpinPr,iSpin) + one * CONJG(cil * ulovulo(loplo,lh))
                           tlmplm%tuloulo_newer(m,mp,lo,lop,ntyp,iSpinPr,iSpin) = &
                         & tlmplm%tuloulo_newer(m,mp,lo,lop,ntyp,iSpinPr,iSpin) + one * CONJG(cil * ulovulo(loplo,lh))
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      ! TODO: Do *not* symmetrize the local hamiltonian in dfpt.
      ! Add the diagonal terms from the spherical Hamiltonian. These terms have
      ! to be made Hermitian. If second variation is switched on, the t-matrices
      ! contain only the contributions from the non-spherical Hamiltonian.
      IF (.NOT.input%secvar.AND.iSpinPr==iSpin.AND..NOT.l_V1) THEN
         DO lo = 1,atoms%nlo(ntyp)
            l = atoms%llo(lo,ntyp)
            DO m = -l,l
               lm = l* (l+1) + m
               !IF (.NOT.l_dfpt) THEN
               IF (.TRUE.) THEN
                  tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) &
                         + 0.5 * usdus%uulon(lo,ntyp,iSpinPr) &
                        * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) )
                  tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin) + &
                        + 0.5 * usdus%uulon(lo,ntyp,iSpinPr) &
                        * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) )
                  tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) &
                      + 0.5 * usdus%dulon(lo,ntyp,iSpinPr) &
                      * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) ) &
                      + 0.5 * usdus%uulon(lo,ntyp,iSpinPr)
                  tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin) + &
                      + 0.5 * usdus%dulon(lo,ntyp,iSpinPr) &
                      * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) ) &
                      + 0.5 * usdus%uulon(lo,ntyp,iSpinPr)
                  tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) &
                                                        & + 0.5 * usdus%uulon(lo,ntyp,iSpinPr) &
                                                        & * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) )
                  tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) &
                                                        & + 0.5 * usdus%dulon(lo,ntyp,iSpinPr) &
                                                        & * ( enpara%el0(l,ntyp,iSpinPr)+enpara%ello0(lo,ntyp,iSpinPr) ) &
                                                        & + 0.5 * usdus%uulon(lo,ntyp,iSpinPr)
                  IF (atoms%ulo_der(lo,ntyp).GE.1) THEN
                     tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr)
                     tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr)
                     tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr)
                     tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr)

                     tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr)
                     tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr)
                  END IF
                  !+apw_lo
                  IF (atoms%l_dulo(lo,ntyp)) THEN
                     tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5
                     tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5
                     tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)= 0.0
                     tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = 0.0
                     tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5
                     tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = 0.0
                  END IF
                  !+apw_lo
               ELSE
                  tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) &
                      + usdus%uulon(lo,ntyp,iSpinPr) &
                      * enpara%ello0(lo,ntyp,iSpinPr)
                      tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)&
                      + usdus%uulon(lo,ntyp,iSpinPr) &
                     * enpara%ello0(lo,ntyp,iSpinPr)
                  tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) &
                     + usdus%dulon(lo,ntyp,iSpinPr) &
                     * enpara%ello0(lo,ntyp,iSpinPr) &
                     + 0.0 * usdus%uulon(lo,ntyp,iSpinPr)
                     tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)&
                     + usdus%dulon(lo,ntyp,iSpinPr) &
                     * enpara%ello0(lo,ntyp,iSpinPr) &
                     + 0.0 * usdus%uulon(lo,ntyp,iSpinPr)
                  tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) &
                                                        & + usdus%uulon(lo,ntyp,iSpinPr) &
                                                        & * enpara%el0(l,ntyp,iSpinPr)
                  tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) &
                                                        & + usdus%dulon(lo,ntyp,iSpinPr) &
                                                        & * enpara%el0(l,ntyp,iSpinPr) &
                                                        & + 1.0 * usdus%uulon(lo,ntyp,iSpinPr)
                  ! TODO: Implement boundary term.
                  IF (atoms%ulo_der(lo,ntyp).GE.1) THEN
                     CALL juDFT_error("ulo_der>0 for DFPT" ,calledby ="tlo")
                     tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr) !TODO: 1.0 or 0.0?
                     tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr) !TODO: 1.0 or 0.0?
                     tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr)
                     tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr)
                     tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%uuilon(lo,ntyp,iSpinPr) !TODO: 1.0 or 0.0?
                     tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5 * usdus%duilon(lo,ntyp,iSpinPr) !TODO: 1.0 or 0.0?
                  END IF
                  IF (atoms%l_dulo(lo,ntyp)) THEN
                     CALL juDFT_error("l_dulo for DFPT" ,calledby ="tlo")
                     tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tuulo(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5
                     tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)=tlmplm%h_LO(lm,m,lo+mlo,iSpinPr,iSpin)+0.5
                     tlmplm%h_LO(lm+s,m,lo+mlo,iSpinPr,iSpin)=0.0
                     tlmplm%tdulo(lm,m,lo+mlo,iSpinPr,iSpin) = 0.0
                     tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) = tlmplm%tulou(lm,m,lo+mlo,iSpinPr,iSpin) + 0.5
                     tlmplm%tulod(lm,m,lo+mlo,iSpinPr,iSpin) = 0.0
                  END IF
               END IF
            END DO
         END DO
         DO lop = 1,atoms%nlo(ntyp)
            lp = atoms%llo(lop,ntyp)
            DO lo = atoms%lo1l(lp,ntyp),lop
               loplo = ((lop-1)*lop)/2 + lo
               loplo_new = (lop-1) * atoms%nlo(ntyp) + lo
               DO m = -lp,lp
                  !IF (.NOT.l_dfpt) THEN
                  IF (.TRUE.) THEN
                     tlmplm%tuloulo(m,m,loplo+mlolo,iSpinPr,iSpin) = tlmplm%tuloulo(m,m,loplo+mlolo,iSpinPr,iSpin) &
                                                              & + 0.5 * (enpara%ello0(lop,ntyp,iSpinPr) &
                                                              & +         enpara%ello0(lo,ntyp,iSpinPr)) &
                                                              & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                                                              & + 0.5 * (usdus%ulouilopn(lop,lo,ntyp,iSpinPr) &
                                                              & +        usdus%ulouilopn(lo,lop,ntyp,iSpinPr))
                     !tlmplm%tuloulo_new(m,m,mlolo_new+loplo_new,iSpinPr,iSpin) = tlmplm%tuloulo_new(m,m,mlolo_new+loplo_new,iSpinPr,iSpin) &
                     !                                         & + 0.5 * (enpara%ello0(lop,ntyp,iSpinPr) &
                     !                                         & +         enpara%ello0(lo,ntyp,iSpinPr)) &
                     !                                         & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                     !                                         & + 0.5 * (usdus%ulouilopn(lop,lo,ntyp,iSpinPr) &
                     !                                         & +        usdus%ulouilopn(lo,lop,ntyp,iSpinPr))
                     tlmplm%tuloulo_newer(m,m,lop,lo,ntyp,iSpinPr,iSpin) = tlmplm%tuloulo_newer(m,m,lop,lo,ntyp,iSpinPr,iSpin) &
                                                              & + 0.5 * (enpara%ello0(lop,ntyp,iSpinPr) &
                                                              & +         enpara%ello0(lo,ntyp,iSpinPr)) &
                                                              & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                                                              & + 0.5 * (usdus%ulouilopn(lop,lo,ntyp,iSpinPr) &
                                                              & +        usdus%ulouilopn(lo,lop,ntyp,iSpinPr))
                     IF (.NOT.lop.EQ.lo) THEN
                        !loplo_new = (lo-1) * atoms%nlo(ntyp) + lop
                        !tlmplm%tuloulo_new(m,m,mlolo_new+loplo_new,iSpinPr,iSpin) = tlmplm%tuloulo_new(m,m,mlolo_new+loplo_new,iSpinPr,iSpin) &
                        !                                         & + 0.5 * (enpara%ello0(lo,ntyp,iSpinPr) &
                        !                                         & +         enpara%ello0(lop,ntyp,iSpinPr)) &
                        !                                         & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                        !                                         & + 0.5 * (usdus%ulouilopn(lo,lop,ntyp,iSpinPr) &
                        !                                         & +        usdus%ulouilopn(lop,lo,ntyp,iSpinPr))
                        tlmplm%tuloulo_newer(m,m,lo,lop,ntyp,iSpinPr,iSpin) = tlmplm%tuloulo_newer(m,m,lo,lop,ntyp,iSpinPr,iSpin) &
                                                                 & + 0.5 * (enpara%ello0(lo,ntyp,iSpinPr) &
                                                                 & +         enpara%ello0(lop,ntyp,iSpinPr)) &
                                                                 & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                                                                 & + 0.5 * (usdus%ulouilopn(lo,lop,ntyp,iSpinPr) &
                                                                 & +        usdus%ulouilopn(lop,lo,ntyp,iSpinPr))
                     END IF
                  ELSE
                     !tlmplm%tuloulo(m,m,loplo+mlolo,iSpinPr,iSpin) = tlmplm%tuloulo(m,m,loplo+mlolo,iSpinPr,iSpin) &
                     !                                            & + enpara%ello0(lo,ntyp,iSpinPr) &
                     !                                            & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                     !                                            & + usdus%ulouilopn(lop,lo,ntyp,iSpinPr)
                     tlmplm%tuloulo_newer(m,m,lop,lo,ntyp,iSpinPr,iSpin) = tlmplm%tuloulo_newer(m,m,lop,lo,ntyp,iSpinPr,iSpin) &
                                                                 & + enpara%ello0(lo,ntyp,iSpinPr) &
                                                                 & * usdus%uloulopn(lop,lo,ntyp,iSpinPr) &
                                                                 & + usdus%ulouilopn(lop,lo,ntyp,iSpinPr)
                     IF (.NOT.lop.EQ.lo) THEN
                        tlmplm%tuloulo_newer(m,m,lo,lop,ntyp,iSpinPr,iSpin) = tlmplm%tuloulo_newer(m,m,lo,lop,ntyp,iSpinPr,iSpin) &
                                                                    & + enpara%ello0(lop,ntyp,iSpinPr) &
                                                                    & * usdus%uloulopn(lo,lop,ntyp,iSpinPr) &
                                                                    & + usdus%ulouilopn(lo,lop,ntyp,iSpinPr)
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END IF
     
   END SUBROUTINE tlo
END MODULE m_tlo
