      MODULE m_force_kin_sf_mt

      CONTAINS

      SUBROUTINE sf_ham(
     >                  )

      USE m_constants, ONLY : pimach
      USE m_gaunt
      USE m_cotra

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,ntype,ne

!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),invarind(natd)
      INTEGER, INTENT (IN) :: invarop(natd,???),invtab(???)
      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL   , INTENT (IN) :: el(0:lmaxd,ntypd),eig(nobd),rmt(ntypd)
      REAL   , INTENT (IN) :: us(0:lmaxd,ntypd),uds(0:lmaxd,ntypd)
      REAL   , INTENT (IN) :: amat(3,3),bmat(3,3)
      COMPLEX, INTENT (IN) :: acof(nobd,0:(lmaxd+1)*lmaxd,natd)
      COMPLEX, INTENT (IN) :: bcof(nobd,0:(lmaxd+1)*lmaxd,natd)
      COMPLEX, INTENT (INOUT) :: force(3,ntypd),fsfmt(3,ntypd)

!     .. Local Scalars ..
      INTEGER :: natom,itype,natrun,t,l,iband,lp,m,lm,mp,lmp,it,is,isinv
      REAL    :: tpi,sfp
      COMPLEX :: img,lsum,afac,bfac,asum,bsum,fac

!     .. Local Arrays ..
      COMPLEX :: coeff(3,-1:1),tsum(3),vec(3),gvint(3),vecsum(3)
      COMPLEX :: starsum(3),force_sf_mt(3)


      tpi = 2. * pimach()
      sfp = sqrt(2.*tpi)
      img = cmplx(0.0,1.0)
      force_sf_mt = 0.0

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi/3.)
      coeff(1, 1) =    -sqrt(tpi/3.)
      coeff(2,-1) = img*sqrt(tpi/3.)
      coeff(2, 1) = img*sqrt(tpi/3.)
      coeff(3, 0) =  sqrt(2.*tpi/3.)

      natom = 1
      DO itype = 1,ntype
      DO natrun = natom,natom+neq(itype)-1
        force_sf_mt = 0.0

!         nd = ntypsy(natom)

        tsum = 0.0
        DO t = -1,1

          lsum = 0.0
          DO l = 0,lmax(itype)

            DO iband = 1,ne
              afac = (el(l,itype,jsp,nwd)-eig(iband))*us(l,itype)
              bfac = (el(l,itype,jsp,nwd)-eig(iband))*uds(l,itype)
     +                                               +us(l,itype)

              asum = 0.0
              bsum = 0.0
              DO lp = iabs(l-1),min(l+1,lmax(itype)),2
              DO m = -l,l
                lm = l*(l+1) + m + 1

                mp = m+t
                lmp = lp*(lp+1) + mp
                IF (lp.ge.abs(mp)) THEN
                  fac =  us(lp,itype)*conjg(acof(iband,lmp,natrun))
     *                + uds(lp,itype)*conjg(bcof(iband,lmp,natrun))
                  fac = fac * gaunt1(lp,l,1,mp,m,t,lmaxd)
                  asum = asum + fac *       acof(iband,lm,natrun)
                  bsum = bsum + fac *       bcof(iband,lm,natrun)
                END IF

                mp = m-t
                lmp = lp*(lp+1) + mp
                IF (lp.ge.abs(mp)) THEN
                  fac =  us(lp,itype)*      acof(iband,lmp,natrun)
     *                + uds(lp,itype)*      bcof(iband,lmp,natrun)
                  fac = fac * gaunt1(lp,l,1,mp,m,t,lmaxd)
                  asum = asum + fac * conjg(acof(iband,lm,natrun))
                  bsum = bsum + fac * conjg(bcof(iband,lm,natrun))
                END IF

              END DO ! m
              END DO ! lp

              lsum = lsum + afac * asum + bfac * bsum

            END DO ! iband

          END DO ! l

          tsum(:) = tsum(:) + coeff(:,t)*lsum

        END DO ! t

        tsum(:) = tsum(:)*0.5*rmt(itype)**2

        vec(i) = tsum(:)/neq(itype)
        CALL cotra1(vec,gvint,bmat)
        vecsum(i) = 0.0

        DO it = 1,invarind(natom)
          is =invarop(natom,it)
          isinv = invtab(is)
!           IF (odi%d1) isinv = ods%ngopr(natom)

!           IF (.NOT.odi%d1) THEN
            vec = matmul(    mrot(:,:,isinv),gvint)
!           ELSE
!             vec = matmul(ods%mrot(:,;,isinv),gvint)
!           END IF
          vecsum(:) = vecsum(:) + vec(:)

! !           DO i = 1,3
! !             vec(i) = zero
! !             DO j = 1,3
! !               IF (.NOT.odi%d1) THEN
! !                  vec(i) = vec(i) + mrot(i,j,isinv)*gvint(j)
! !               ELSE
! !                  vec(i) = vec(i) + ods%mrot(i,j,isinv)*gvint(j)
! !               END IF
! !             END DO
! !           END DO
! !           DO i = 1,3
! !             vecsum(i) = vecsum(i) + vec(i)
! !           END DO

        END DO

        CALL cotra0(vecsum,starsum,amat)
        force_sf_mt(:) = force_sf_mt(:) + starsum(:)/invarind(natrun)

      END DO ! natrun

        force(:,itype) = force(:,itype) + real(force_sf_mt(:))
        fsfmt(:,itype) = real(force_sf_mt(:))

        natom = natom + neq(itype)
      END DO ! itype

      END SUBROUTINE sf_ham



      SUBROUTINE sf_potden(
     >                     )

      USE

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: natd,ntypd,ntype,ntypsd,nlhd,memd,jmtd
      INTEGER, INTENT (IN) :: jspd

!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ntypsy(natd),neq(ntypd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: llh(0:nlhd,ntypsd),nmem(nlhd,0:ntypsd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd)
      REAL   , INTENT (IN) :: vr(jmtd,ntypd,jspd),rmt(ntypd)
      REAL   , INTENT (IN) :: rho(jmtd,0:nlhd,ntypsd,jspd)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      REAL   , INTENT (INOUT) :: force(3,ntypd),fsfmt(3,ntypd)

!     .. Local Scalars ..
      INTEGER :: natom,itype,nd,lh,t,jmem
      REAL    :: tpi,sfp,vr0
      COMPLEX :: img

!     .. Local Arrays ..
      COMPLEX :: force_sf_mt(3),coeff(3,-1:1),tsum(3)

      tpi = 2. * pimach()
      sfp = sqrt(2.*tpi)
      img = cmplx(0.0,1.0)
      force_sf_mt = 0.0

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi/3.)
      coeff(1, 1) =    -sqrt(tpi/3.)
      coeff(2,-1) = img*sqrt(tpi/3.)
      coeff(2, 1) = img*sqrt(tpi/3.)
      coeff(3, 0) =  sqrt(2.*tpi/3.)

!     load core density to subtract it from total density
!     not necessary since core density is l=0 and we need
!     the l=1 component of the density => no core contained
!     in there

      natom = 1
      DO itype = 1,ntype

        vr0 = sfp * vr(jri(itype),itype,jsp) / rmt(itype)

        nd = ntypsy(natom)

        tsum = 0.0
        DO lh = 0,nlh(nd)
          IF (llh(lh,nd).ne.1) CYCLE
        DO t = -1,1
        DO jmem = 1,nmem(lh,nd)
          IF (mlh(jmem,lh,nd).ne.t) CYCLE

          tsum(:) = tsum(:) + conjg(coeff(:,t)) * vr0
     *                      * rho(jri(itype),lh,itype,jsp)
     *                      * clnu(jmem,lh,nd)

        END DO ! jmem
        END DO ! t
        END DO ! lh


        force_sf_mt(:) = force_sf_mt(:) + tsum(:)

        force(:,itype) = force(:,itype) + real(force_sf_mt(:))
        fsfmt(:,itype) = real(force_sf_mt(:))

        natom = natom + neq(itype)
      END DO ! itype

      END SUBROUTINE sf_potden

      END MODULE m_force_kin_sf_mt