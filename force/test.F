      MODULE m_test

      CONTAINS
      SUBROUTINE test(
     >                nvd,jspd,ntypd,nop,natd,neigd,nobd,nbasfcn,jsp,
     >                ntype,ne,nk,omtil,k1,k2,k3,nv,neq,mrot,invarop,
     >                invarind,invtab,bkpt,amat,bmat,taual,eig,rmt,we,
     >                l_geo,z,odi,ods,
     X                force,f_a12)

      USE m_od_types
      USE m_constants, ONLY : pimach
      USE m_sphbes

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nvd,jspd,ntypd,nop,natd,neigd
      INTEGER, INTENT (IN) :: nobd,nbasfcn,jsp,ntype,ne,nk
      REAL   , INTENT (IN) :: omtil

!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: nv(jspd),neq(ntypd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invarop(natd,nop),invarind(natd)
      INTEGER, INTENT (IN) :: invtab(nop)
      REAL   , INTENT (IN) :: bkpt(3),amat(3,3),bmat(3,3),taual(3,natd)
      REAL   , INTENT (IN) :: eig(neigd),rmt(ntypd),we(nobd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL   , INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
!-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
!+odim
      COMPLEX, INTENT (INOUT) :: f_a12(3,ntypd)
      REAL   , INTENT (INOUT) :: force(3,ntypd,jspd)

!     .. Local Scalars ..
      INTEGER :: kn,km,iband,iatom,itype,ieq,i,it,is,isinv
      REAL    :: tpi,norm,kGnm,r
      COMPLEX :: img,bandcontribution,expf

!     .. Local Arrays ..
      REAL :: force_a12(3,ntypd),G(3,nvd),kG(3),kGreal(3,nvd),dG(3)
      REAL :: dGreal(3),fj(0:1),gv(3),vec(3),gvint(3),vecsum(3)
      REAL :: starsum(3),ab(3,3)

      IF (.not.any(l_geo)) RETURN ! if no forces on any atom shall be calculated, skip the routine

      tpi = 2 * pimach()
      img = cmplx(0.0,1.0)

      IF (nk.eq.1) THEN
        WRITE (3001,'(2f15.10)') tpi,1./tpi
        ab = matmul(amat,bmat)
        WRITE (3001,'(3f15.10)') ab(1,1),ab(1,2),ab(1,3)
        WRITE (3001,'(3f15.10)') ab(2,1),ab(2,2),ab(2,3)
        WRITE (3001,'(3f15.10)') ab(3,1),ab(3,2),ab(3,3)
        ab = matmul(bmat,amat)
        WRITE (3001,'(3f15.10)') ab(1,1),ab(1,2),ab(1,3)
        WRITE (3001,'(3f15.10)') ab(2,1),ab(2,2),ab(2,3)
        WRITE (3001,'(3f15.10)') ab(3,1),ab(3,2),ab(3,3)
      END IF

      DO kn = 1,nv(jsp)
        G(1,kn) = k1(kn,jsp)
        G(2,kn) = k2(kn,jsp)
        G(3,kn) = k3(kn,jsp)
        kG = bkpt + G(:,kn)
        kGreal(:,kn) = matmul(kG,bmat)
      END DO !

      force_a12 = 0
      DO kn = 1,nv(jsp)
!         G(1,kn) = k1(kn,jsp)
!         G(2,kn) = k2(kn,jsp)
!         G(3,kn) = k3(kn,jsp)
!         kG = bkpt + G(:,kn)
!         kGreal(:,kn) = matmul(kG,bmat)

        IF (nk.eq.1) THEN
          WRITE (3000,'(i4,3f15.10)') kn,G(1,kn),G(2,kn),G(3,kn)
        END IF

        DO km = 1,nv(jsp)!1,kn-1
          IF (kn.eq.km) CYCLE
          dG = G(:,kn) - G(:,km)
          dGreal = matmul(dG,bmat)
          norm = sqrt(dot_product(dGreal,dGreal))
!           kGnm = dot_product(kGreal(:,kn),kGreal(:,km))
          kGnm = 1.0*dot_product(kGreal(:,km),kGreal(:,km))
!      +         + 0.5*dot_product(kGreal(:,kn),kGreal(:,kn))

          bandcontribution = 0
          DO iband = 1,nobd
            bandcontribution = bandcontribution
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     +      + we(iband)*z(kn,iband)*z(km,iband)*(0.5*kGnm-eig(iband))
#else
     +      + we(iband)*conjg(z(kn,iband))*z(km,iband)
     *        *(0.5*kGnm-eig(iband))
#endif
          END DO ! iband
          IF (nk.eq.1) THEN
            WRITE (3000,'(2i4,2f15.10)') kn,km,bandcontribution
          END IF

          iatom = 0
          DO itype = 1,ntype
            IF (.not.l_geo(itype)) THEN ! only calculate forces for atoms that require it
              iatom = iatom + neq(itype)
              CYCLE
            END IF
            r = norm*rmt(itype)
            CALL sphbes(1,r,fj)
          DO ieq = 1,neq(itype)
            iatom = iatom + 1

            expf = exp(-tpi*img*dot_product(dG,taual(:,iatom)))
     *             * 2*tpi*img/omtil
            DO i = 1,3
              gv(i) = -1.*real(dGreal(i)*expf*bandcontribution*fj(1)
     *                                              *rmt(itype)**2/norm)
              vec(i) = real(gv(i))/neq(itype)
            END DO ! i
            IF (nk.eq.1) THEN
              WRITE (3000,'(3i4,3f15.10)') kn,km,iatom,gv
            END IF

!           starsumgedoens
            gvint = matmul(bmat,vec)/tpi
            vecsum = 0

            DO it = 1,invarind(iatom) ! loop over invariant operations
              is = invarop(iatom,it)
              isinv = invtab(is)
              IF (odi%d1) isinv = ods%ngopr(iatom)

              IF (.NOT.odi%d1) THEN
                vecsum = vecsum+matmul(mrot(:,:,isinv),gvint(:))
              ELSE
                vecsum = vecsum+matmul(ods%mrot(:,:,isinv),gvint(:))
              END IF
            END DO ! it

            starsum = matmul(amat,vecsum)
            force_a12(:,itype) = force_a12(:,itype)
     +                         + starsum(:)/invarind(iatom)

            IF (nk.eq.1) THEN
              WRITE (3000,'(4i4,3f15.10)')
     &          kn,km,itype,iatom,force_a12(:,itype)
            END IF

          END DO ! ieq
          END DO ! itype

        END DO ! km

      END DO ! kn
      WRITE (850,*) 'old test:'
      DO itype = 1,ntype ! loop over atom types
!       sum up forces
        DO i = 1,3 ! directions
          force(i,itype,jsp) = force(i,itype,jsp) + force_a12(i,itype)
          f_a12(i,itype)     = f_a12(i,itype)     + force_a12(i,itype)
        END DO ! i
        WRITE (850,*) (force_a12(i,itype),i=1,3)
      END DO ! itype atom types

      END SUBROUTINE test

      SUBROUTINE test2(
     >                 ntypd,ntype,jspd,jsp,nvd,nobd,neigd,natd,omtil,
     >                 lmax,nv,k1,k2,k3,bkpt,bmat,eig,rmt,taual,we,z,
     >                 neq,nbasfcn,nop,invarind,invarop,invtab,mrot,
     >                 amat )

!     does the same as above, just using numerical integration
!     right now, it can only be used with kinetic energy as Laplacian
!     nabla(left)nabla(right) needs to be implemented

      USE m_constants, ONLY : pimach
      USE m_lebMesh

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,ntype,jspd,jsp,nvd,nobd,neigd,natd
      INTEGER, INTENT (IN) :: nbasfcn
      REAL   , INTENT (IN) :: omtil

!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: lmax(ntypd),nv(jspd),neq(ntypd)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      REAL   , INTENT (IN) :: bkpt(3),bmat(3,3),eig(neigd),rmt(ntypd)
      REAL   , INTENT (IN) :: taual(3,natd),we(nobd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL   , INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif

      INTEGER, INTENT (IN) :: nop
      INTEGER, INTENT (IN) :: invarind(natd),invarop(natd,nop)
      INTEGER, INTENT (IN) :: invtab(nop),mrot(3,3,nop)
      REAL   , INTENT (IN) :: amat(3,3)

!     .. Local Scalars ..
      INTEGER :: meshptmax,itype,kn,iband,iatom,i,jatom,ieq
      REAL    :: normsq,prefactor,tpi
      COMPLEX :: ppw,fpw,expf,img

!     .. Local Arrays ..
      REAL   , ALLOCATABLE :: meshpt(:,:,:),weight(:,:)
      INTEGER :: nummeshpts(ntypd)
      REAL    :: G(3,nvd),kG(3,nvd),kGreal(3,nvd)
      REAL    :: kineticfactor(nvd,nobd),gv(3),r(3)
      COMPLEX :: force_a12(3,ntypd)

      INTEGER :: it,is,isinv
      COMPLEX :: gvint(3),vecsum(3),starsum(3)

      meshptmax = 1202
      ALLOCATE ( meshpt(3,meshptmax,ntypd),weight(meshptmax,ntypd) )

      tpi = 2. * pimach()
      img = cmplx(0.0,1.0)

!     obtain mesh points and weights
      DO itype = 1,ntype
! !         CALL getNumLebMeshPts(lmax(itype)+1,nummeshpts(itype))
! !         CALL genLebMesh(lmax(itype)+1,nummeshpts(itype),
! !      <                  meshpt(:,:,itype),weight(:,itype))
        CALL getNumLebMeshPts(58,nummeshpts(itype))
        CALL genLebMesh(58,nummeshpts(itype),
     <                  meshpt(:,:,itype),weight(:,itype))
      END DO ! itype

!     loop over plane waves
      DO kn = 1,nv(jsp)
        G(1,kn) = k1(kn,jsp)
        G(2,kn) = k2(kn,jsp)
        G(3,kn) = k3(kn,jsp)
        kG(:,kn) = bkpt(:) + G(:,kn)
        kGreal(:,kn) = matmul(kG(:,kn),bmat)
!       this is only useable using the Laplacian for kinetic energy
        normsq = dot_product(kGreal(:,kn),kGreal(:,kn))
        DO iband = 1,nobd
          kineticfactor(kn,iband) = 0.5*normsq - eig(iband)
        END DO ! iband
      END DO ! kn

      force_a12 = 0.0
      WRITE (851,*) 'new test:'
      iatom = 1
      DO itype = 1,ntype
      DO ieq = 1,neq(itype)
        prefactor = 2.*tpi*rmt(itype)**2/omtil
        gv = 0.0
        DO i = 1,nummeshpts(itype)
          r = rmt(itype)*matmul(bmat,meshpt(:,i,itype))/tpi
          r = r + taual(:,iatom)
          DO iband = 1,nobd
            ppw = 0.0
            fpw = 0.0
            DO kn = 1,nv(jsp)
!               expf = 0.0
!               jatom = iatom
!               DO ieq = 1,neq(itype)
!                 jatom = jatom + 1
!                 expf = expf + exp(tpi*img
!      *               * dot_product(kG(:,kn),r+taual(:,jatom)))
                expf = exp(tpi*img*dot_product(kG(:,kn),r))
!               END DO
              ppw = ppw + z(kn,iband) * expf
              fpw = fpw + z(kn,iband) * expf * kineticfactor(kn,iband)
            END DO ! kn
            gv(:) = gv(:) + weight(i,itype) * meshpt(:,i,itype)
     *                    * prefactor * we(iband) * conjg(ppw) * fpw
          END DO ! iband
        END DO ! i

!       starsumgedoens
        gvint = matmul(bmat,gv)/tpi/neq(itype)
        vecsum = 0

        DO it = 1,invarind(iatom) ! loop over invariant operations
          is = invarop(iatom,it)
          isinv = invtab(is)
!           IF (odi%d1) isinv = ods%ngopr(iatom)

!           IF (.NOT.odi%d1) THEN
            vecsum = vecsum+matmul(mrot(:,:,isinv),gvint(:))
!           ELSE
!             vecsum = vecsum+matmul(ods%mrot(:,:,isinv),gvint(:))
!           END IF
        END DO ! it

        starsum = matmul(amat,vecsum)
        force_a12(:,itype) = force_a12(:,itype)
     +                     + starsum(:)/invarind(iatom)

!         force_a12(:,itype) = gv(:)
!         WRITE (851,*) (real(force_a12(i,itype)),i=1,3)
        iatom = iatom + 1
!         iatom = iatom + neq(itype)
      END DO ! ieq
        WRITE (851,*) (real(force_a12(i,itype)),i=1,3)
      END DO ! itype

      DEALLOCATE ( meshpt,weight )

      END SUBROUTINE test2

      SUBROUTINE test3(
     >                 jsp,jspd,nobd,neigd,ntypd,ntype,natd,nbasfcn,
     >                 nop,nvd,lmaxd,omtil,
     >                 nv,neq,k1,k2,k3,invarind,invarop,invtab,mrot,
     >                 ngopr,amat,bmat,eig,rmt,taual,we,bkpt,z,
     X                 f_a12,force )
!     same as above, but this time, the plane waves are expanded in Y

      USE m_constants, ONLY : pimach
      USE m_ylm
      USE m_sphbes
      USE m_gaunt

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jsp,jspd,nobd,neigd,ntypd,ntype,natd
      INTEGER, INTENT (IN) :: nbasfcn,nop,nvd,lmaxd
      REAL   , INTENT (IN) :: omtil

!     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nv(jspd),neq(ntypd)
      INTEGER, INTENT (IN) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      INTEGER, INTENT (IN) :: invarind(natd),invarop(natd,nop)
      INTEGER, INTENT (IN) :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      REAL   , INTENT (IN) :: amat(3,3),bmat(3,3),eig(neigd),rmt(ntypd)
      REAL   , INTENT (IN) :: taual(3,natd),we(nobd),bkpt(3)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL   , INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
      COMPLEX, INTENT (INOUT) :: f_a12(3,ntypd)
      REAL   , INTENT (INOUT) :: force(3,ntypd,jspd)

!     .. Local Scalars ..
      INTEGER :: lmax,kn,iband,iatom,itype,ieq,l,m,lm,t,lp,mp,lmp,it,is
      INTEGER :: isinv,i
      REAL    :: tpi,normsq,r,r2vol
      COMPLEX :: img,noband,fpil

!     .. Local Arrays ..
!       REAL    :: G(3,nvd),kG(3,nvd),kGreal(3,nvd)
!       REAL    :: kineticfactor(nobd,nvd)
      COMPLEX :: force_a12(3,ntypd),gv(3),Ygaunt(3)!,expf(nvd,natd)
      COMPLEX :: coeff(3,-1:1),gvint(3),vecsum(3),starsum(3)
      REAL   , ALLOCATABLE :: bsl(:,:,:),G(:,:),kG(:,:),kGreal(:,:)
      REAL   , ALLOCATABLE :: kineticfactor(:,:)
      COMPLEX, ALLOCATABLE :: ylm(:,:),ppw(:,:),fpw(:,:),expf(:,:)

      lmax = 2*lmaxd
      tpi = 2.0 * pimach()
      img = cmplx(0.0,1.0)

      ALLOCATE ( bsl(nvd,0:lmax,ntypd),ylm((lmax+1)**2,nvd) )
      ALLOCATE ( ppw(nobd,(lmax+1)**2),fpw(nobd,(lmax+1)**2) )
      ALLOCATE ( G(3,nvd),kG(3,nvd),kGreal(3,nvd),expf(nvd,natd) )
      ALLOCATE ( kineticfactor(nobd,nvd) )

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi/3.)
      coeff(1, 1) =    -sqrt(tpi/3.)
      coeff(2,-1) = img*sqrt(tpi/3.)
      coeff(2, 1) = img*sqrt(tpi/3.)
      coeff(3, 0) =  sqrt(2.*tpi/3.)

!       WRITE (5555,*) 'kpt',bkpt(1),bkpt(2),bkpt(3)
!       WRITE (5556,*) 'kpt',bkpt(1),bkpt(2),bkpt(3)

!     loop over plane waves
      DO kn = 1,nv(jsp)
        G(1,kn) = k1(kn,jsp)
        G(2,kn) = k2(kn,jsp)
        G(3,kn) = k3(kn,jsp)
!         WRITE (5555,'(4i10)') kn,k1(kn,jsp),k2(kn,jsp),k3(kn,jsp)
!         kG(:,kn) = matmul(G(:,kn),bmat)
!         WRITE (5556,'(i10,3f10.6)') kn,kG(1,kn),kG(2,kn),kG(3,kn)
        kG(:,kn) = bkpt(:) + G(:,kn)
        kGreal(:,kn) = matmul(kG(:,kn),bmat)
!       this is only useable using the Laplacian for kinetic energy
        normsq = dot_product(kGreal(:,kn),kGreal(:,kn))
        DO iband = 1,nobd
          kineticfactor(iband,kn) = 0.5*normsq - eig(iband)
        END DO ! iband
        CALL ylm4(lmax,kGreal(:,kn),ylm(:,kn))
!         CALL ylm4(lmax,matmul(G(:,kn),bmat),ylm(:,kn))
        iatom = 1
        DO itype = 1,ntype
          r = sqrt(normsq)*rmt(itype)
!           r = sqrt(dot_product(matmul(G(:,kn),bmat),
!      ,                                 matmul(G(:,kn),bmat)))*rmt(itype)
          CALL sphbes(lmax,r,bsl(kn,:,itype))
          DO ieq = 1,neq(itype)
            expf(kn,iatom) =
     =                 exp(tpi*img*dot_product(kG(:,kn),taual(:,iatom)))
!      =                 exp(tpi*img*dot_product(G(:,kn),taual(:,iatom)))
            iatom = iatom + 1
          END DO ! ieq
        END DO ! itype
      END DO ! kn

      force_a12 = 0.0
      WRITE (852,*) 'new test3'

      iatom = 1
      DO itype = 1,ntype
        r2vol = rmt(itype)**2/omtil
      DO ieq = 1,neq(itype)

        gv  = 0.0
        ppw = 0.0
        fpw = 0.0

        DO l = 0,lmax-1 ! (arrays run to lmax, l+1 can only be supported til lmax-1)
          fpil = 2.*tpi*img**l
!         Calculate ppw and fpw
          DO kn = 1,nv(jsp)
            DO m = -l,l
              lm = l*(l+1) + m + 1
              noband = expf(kn,iatom)*conjg(ylm(lm,kn))*bsl(kn,l,itype)
     *               * fpil
              DO iband = 1,nobd
                ppw(iband,lm) = ppw(iband,lm) + z(kn,iband) * noband
                IF (l.gt.0) CYCLE
                fpw(iband,lm) = fpw(iband,lm) + z(kn,iband) * noband
     *                                        * kineticfactor(iband,kn)
              END DO ! iband
            END DO ! m
!           If ppw is used with l, one needs fpw with l-1 (already calculated) and l+1 (calculated now)
            DO m = -l-1,l+1
              lm = (l+1)*(l+2) + m + 1
              noband =expf(kn,iatom)*conjg(ylm(lm,kn))*bsl(kn,l+1,itype)
     *               * fpil * img
              DO iband = 1,nobd
                fpw(iband,lm) = fpw(iband,lm) + z(kn,iband) * noband
     *                                        * kineticfactor(iband,kn)
              END DO ! iband
            END DO ! m
          END DO ! kn
!         Calculate forces
          DO m = -l,l
            lm = l*(l+1) + m + 1
            DO lp = abs(l-1),l+1,2
              DO t = -1,1
                mp = m-t
                IF (lp.lt.abs(mp)) CYCLE
                lmp = lp*(lp+1) + mp + 1
                Ygaunt(:) = gaunt1(l,1,lp,m,t,mp,lmax)*coeff(:,t)
                DO iband = 1,nobd
                  gv(:) = gv(:) + we(iband) * r2vol * Ygaunt(:)
     *                  * conjg(ppw(iband,lm)) * fpw(iband,lmp)
                END DO ! iband
              END DO ! t
            END DO ! lp
          END DO ! m
        END DO ! lmax

        WRITE (852,*) 'atom',iatom,'(',itype,neq(itype),')'
        WRITE (852,*) (gv(i),i=1,3)

!       starsumgedoens
!       k summation is only on IBZ. Expansion to FBZ can be achieved by
!       applying the rotation of the k-points to the atomic positions of
!       equivalent atoms.
!       Until now, the forces for one k-point are calculated for one
!       specific atom, not considering such a rotation. To do so, the
!       resulting forces have to be rotated back to the representative
!       atom. This is done in the next 6 lines.
        vecsum = matmul(bmat,gv)/tpi/neq(itype)

!         is = ngopr(iatom)
!         it = invtab(is)
! 
! !         gvint = matmul(mrot(:,:,it),vecsum)
        gvint = matmul(mrot(:,:,1),vecsum)

        vecsum = 0

        DO it = 1,invarind(iatom) ! loop over invariant operations
          is = invarop(iatom,it)
          isinv = invtab(is)

!           IF (odi%d1) isinv = ods%ngopr(iatom)

!           IF (.NOT.odi%d1) THEN
            vecsum = vecsum + matmul(mrot(:,:,isinv),gvint(:))
!           ELSE
!             vecsum = vecsum+matmul(ods%mrot(:,:,isinv),gvint(:))
!           END IF
        END DO ! it

        is = ngopr(iatom)
        it = invtab(is)

        vecsum = matmul(mrot(:,:,is),vecsum)
!         starsum = matmul(amat,
!      &    matmul(mrot(:,:,invtab(ngopr(iatom))),vecsum))
        starsum = matmul(amat,vecsum)
        force_a12(:,itype) = force_a12(:,itype)
     +                     + real(starsum(:))/invarind(iatom)

        iatom = iatom + 1
      END DO ! ieq
        force(:,itype,jsp) = force(:,itype,jsp) + force_a12(:,itype)
        f_a12(:,itype)     = f_a12(:,itype)     + force_a12(:,itype)
        WRITE (852,*) (force_a12(i,itype),i=1,3)
      END DO ! itype

      DEALLOCATE ( bsl,ylm,ppw,fpw )
      DEALLOCATE ( G,kG,kGreal,expf,kineticfactor )
!       STOP 'debug'
      END SUBROUTINE test3

      END MODULE m_test
