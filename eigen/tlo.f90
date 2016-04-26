MODULE m_tlo
      USE m_juDFT
!***********************************************************************
!     sets up the extra t-matrix elements due to the local orbitals.
!     only non=zero elements are calculated
!
!     p.kurz jul. 1996
!***********************************************************************
      CONTAINS
        SUBROUTINE tlo(atoms,sphhar, jsp,ntyp,enpara,lh0,input,vr,&
             flo,f,g,usdus,uuilon,duilon,ulouilopn, tlmplm )
          !
          !*************** ABBREVIATIONS *****************************************
          ! tuulo      : t-matrix element of the lo and the apw radial fuction
          ! tdulo      : t-matrix element of the lo and the energy derivativ of 
          !              the apw radial fuction
          ! tuloulo    : t-matrix element of two los
          !c***********************************************************************
          !
          USE m_intgr, ONLY : intgr3  
          USE m_gaunt, ONLY: gaunt1
          USE m_types
          IMPLICIT NONE
          TYPE(t_input),INTENT(IN)    :: input
          TYPE(t_sphhar),INTENT(IN)   :: sphhar
          TYPE(t_atoms),INTENT(IN)    :: atoms
          TYPE(t_usdus),INTENT(IN)    :: usdus
          TYPE(t_tlmplm),INTENT(INOUT):: tlmplm
          TYPE(t_enpara),INTENT(IN)   :: enpara
          !     ..
          !     .. Scalar Arguments ..
          INTEGER, INTENT (IN) :: jsp,ntyp ,lh0
          !     ..
          !     .. Array Arguments ..
          REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd)
          REAL,    INTENT (IN) :: f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd)
          REAL,    INTENT (IN) :: flo(atoms%jmtd,2,atoms%nlod)
          REAL,    INTENT (IN) :: uuilon(atoms%nlod,atoms%ntypd),duilon(atoms%nlod,atoms%ntypd)
          REAL,    INTENT (IN) :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntypd)
          !     ..
          !     .. Local Scalars ..
          COMPLEX ci,cil
          INTEGER i,l,lh,lm ,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,lpmin,lpmin0,lpp ,mem,mp,mpp,m,lmx
          !     ..
          !     .. Local Arrays ..
          REAL x(atoms%jmtd),ulovulo(atoms%nlod*(atoms%nlod+1)/2,lh0:sphhar%nlhd)
          REAL uvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd),dvulo(atoms%nlod,0:atoms%lmaxd,lh0:sphhar%nlhd)
          !     ..

          ci = CMPLX(0.0,1.0)
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             DO lp = 0,atoms%lmax(ntyp)
                lmin = ABS(lp-l)
                !               lmin = lp - l
                lmx = lp + l
                DO lh = lh0,sphhar%nlh(atoms%ntypsy(atoms%nat))
                   lpp = sphhar%llh(lh,atoms%ntypsy(atoms%nat))
                   IF ((MOD(l+lp+lpp,2).EQ.1) .OR. (lpp.LT.lmin) .OR.&
                        (lpp.GT.lmx)) THEN
                      uvulo(lo,lp,lh) = 0.0
                      dvulo(lo,lp,lh) = 0.0
                   ELSE
                      DO i = 1,atoms%jri(ntyp)
                         x(i) = (f(i,1,lp)*flo(i,1,lo)+ f(i,2,lp)*flo(i,2,lo))*vr(i,lh)
                      END DO
                      CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),uvulo(lo,lp,lh))
                      DO i = 1,atoms%jri(ntyp)
                         x(i) = (g(i,1,lp)*flo(i,1,lo)+ g(i,2,lp)*flo(i,2,lo))*vr(i,lh)
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
                IF (loplo>size(ulovulo,1))  CALL juDFT_error("loplo too large!!!" ,calledby ="tlo")
                DO lh = lh0,sphhar%nlh(atoms%ntypsy(atoms%nat))
                   lpp = sphhar%llh(lh,atoms%ntypsy(atoms%nat))
                   lmin = ABS(lp - l)
                   lmx = lp + l
                   IF ((MOD(l+lp+lpp,2).EQ.1).OR.(lpp.LT.lmin).OR.(lpp.GT.lmx)) THEN
                      ulovulo(loplo,lh) = 0.0
                   ELSE
                      DO i = 1,atoms%jri(ntyp)
                         x(i) = (flo(i,1,lop)*flo(i,1,lo)+flo(i,2,lop)*flo(i,2,lo))*vr(i,lh)
                      END DO
                      CALL intgr3(x,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),ulovulo(loplo,lh))
                   END IF
                END DO
             END DO
          END DO
          !---> generate the different t matrices
          !---> but first initialize them ( done in eigen )
          !     
          !---> generate the t-matrices. for optimal performance consider only
          !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
          !---> conditions for non-zero gaunt-coeff. i.e.
          !---> |l - l''| <= l' <= l + l'' (triangular condition)
          !---> m' = m + m'' and l + l' + l'' even
          !---> loop over the local orbitals
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             DO m = -l,l
                !--->       loop over the lattice harmonics
                DO lh = lh0,sphhar%nlh(atoms%ntypsy(atoms%nat))
                   lpp = sphhar%llh(lh,atoms%ntypsy(atoms%nat))
                   lpmin0 = ABS(l-lpp)
                   lpmax0 = l + lpp
                   !--->          check that lpmax is smaller than the max l of the
                   !--->          wavefunction expansion at this atom
                   lpmax = MIN(lpmax0,atoms%lmax(ntyp))
                   !--->          make sure that l + l'' + lpmax is even
                   lpmax = lpmax - MOD(l+lpp+lpmax,2)
                   DO mem = 1,sphhar%nmem(lh,atoms%ntypsy(atoms%nat))
                      mpp = sphhar%mlh(mem,lh,atoms%ntypsy(atoms%nat))
                      mp = m + mpp
                      lpmin = MAX(lpmin0,ABS(mp))
                      !--->             make sure that l + l'' + lpmin is even
                      lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                      !--->             loop over l'
                      DO lp = lpmin,lpmax,2
                         lmp = lp* (lp+1) + mp
                         cil = ((ci** (l-lp))*sphhar%clnu(mem,lh,atoms%ntypsy(atoms%nat)))* gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                         tlmplm%tuulo(lmp,m,lo,jsp) = tlmplm%tuulo(lmp,m,lo,jsp) + cil*uvulo(lo,lp,lh)
                         tlmplm%tdulo(lmp,m,lo,jsp) = tlmplm%tdulo(lmp,m,lo,jsp) + cil*dvulo(lo,lp,lh)
                      END DO
                   END DO
                END DO
             END DO
          END DO
          !---> generate the t-matrix including two local orbitals for lo' >= lo
          !---> loop over lo'
          loplo = 0
          DO lop = 1,atoms%nlo(ntyp)
             lp = atoms%llo(lop,ntyp)
             DO mp = -lp,lp
                !--->       loop over the lattice harmonics
                DO lh = lh0,sphhar%nlh(atoms%ntypsy(atoms%nat))
                   lpp = sphhar%llh(lh,atoms%ntypsy(atoms%nat))
                   DO mem = 1,sphhar%nmem(lh,atoms%ntypsy(atoms%nat))
                      mpp = sphhar%mlh(mem,lh,atoms%ntypsy(atoms%nat))
                      m = mp - mpp
                      !--->             loop over lo
                      DO lo = 1,lop
                         l = atoms%llo(lo,ntyp)
                         loplo = ((lop-1)*lop)/2 + lo
                         IF ((ABS(l-lpp).LE.lp) .AND. (lp.LE. (l+lpp)) .AND.&
                              (MOD(l+lp+lpp,2).EQ.0) .AND. (ABS(m).LE.l)) THEN
                            cil = ((ci** (l-lp))*sphhar%clnu(mem,lh,atoms%ntypsy(atoms%nat)))* gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                            tlmplm%tuloulo(mp,m,loplo,jsp) = tlmplm%tuloulo(mp,m,loplo,jsp) + cil*ulovulo(loplo,lh)
                         END IF
                      END DO
                   END DO
                END DO
             END DO
          END DO
          !---> add the diagonal terms from the muffin-tin hamiltonian. these
          !---> terms have to be made hermitian. if second variation is switched
          !---> on, the t-matrices contain only the contributions from the
          !---> non-spherical hamiltonian.
          IF (.NOT.input%secvar) THEN
             DO lo = 1,atoms%nlo(ntyp)
                l = atoms%llo(lo,ntyp)
                DO m = -l,l
                   lm = l* (l+1) + m
                   tlmplm%tuulo(lm,m,lo,jsp) = tlmplm%tuulo(lm,m,lo,jsp) + 0.5 * usdus%uulon(lo,ntyp,jsp) *&
                        ( enpara%el0(l,ntyp,jsp)+enpara%ello0(lo,ntyp,jsp) )
                   tlmplm%tdulo(lm,m,lo,jsp) = tlmplm%tdulo(lm,m,lo,jsp) + 0.5 * usdus%dulon(lo,ntyp,jsp) *&
                        ( enpara%el0(l,ntyp,jsp)+enpara%ello0(lo,ntyp,jsp) ) + 0.5 * usdus%uulon(lo,ntyp,jsp)
                   IF (atoms%ulo_der(lo,ntyp).GE.1) THEN
                      tlmplm%tuulo(lm,m,lo,jsp) = tlmplm%tuulo(lm,m,lo,jsp) + 0.5 * uuilon(lo,ntyp)
                      tlmplm%tdulo(lm,m,lo,jsp) = tlmplm%tdulo(lm,m,lo,jsp) + 0.5 * duilon(lo,ntyp)
                   ENDIF
                   !+apw_lo
                   IF (atoms%l_dulo(lo,ntyp)) THEN         
                      tlmplm%tuulo(lm,m,lo,jsp) = tlmplm%tuulo(lm,m,lo,jsp) + 0.5
                      tlmplm%tdulo(lm,m,lo,jsp) = 0.0
                   ENDIF
                   !+apw_lo
                 END DO
             END DO
             DO lop = 1,atoms%nlo(ntyp)
                lp = atoms%llo(lop,ntyp)
                DO lo = atoms%lo1l(lp,ntyp),lop
                   loplo = ((lop-1)*lop)/2 + lo
                   DO m = -lp,lp
                      tlmplm%tuloulo(m,m,loplo,jsp) = tlmplm%tuloulo(m,m,loplo,jsp) + 0.5* (enpara%ello0(lop,ntyp,jsp)+&
                           enpara%ello0(lo,ntyp,jsp))* usdus%uloulopn(lop,lo,ntyp,jsp) + 0.5* (ulouilopn(lop,lo,ntyp) +&
                           ulouilopn(lo,lop,ntyp))
                   END DO
                END DO
             END DO
          END IF

        END SUBROUTINE tlo
      END MODULE m_tlo
