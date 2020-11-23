MODULE m_force_a12_lv2 ! Klueppelberg (force level 2)
CONTAINS
   SUBROUTINE force_a12_lv2(jsp,jspd,nobd,neigd,ntypd,ntype,natd,nbasfcn,nop,nvd,lmaxd,omtil,nv,neq,k1,k2,k3,invarind,invarop,invtab,mrot,ngopr,amat,bmat,eig,rmt,taual,we,bkpt,zMat,f_a12,force )
      USE m_constants
      USE m_ylm
      USE m_sphbes
      USE m_gaunt
      USE m_types_mat

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
      TYPE(t_mat), INTENT(IN)             :: zMat
      COMPLEX, INTENT (INOUT) :: f_a12(3,ntypd)
      REAL   , INTENT (INOUT) :: force(3,ntypd,jspd)

!     .. Local Scalars ..
      INTEGER :: lmax,kn,iband,iatom,itype,ieq,l,m,lm,t,lp,mp,lmp,it,is
      INTEGER :: isinv,i
      REAL    :: normsq,r,r2vol
      COMPLEX :: img,noband,fpil

!     .. Local Arrays ..
!       REAL    :: G(3,nvd),kG(3,nvd),kGreal(3,nvd)
!       REAL    :: kineticfactor(nobd,nvd)
      COMPLEX :: force_a12(3,ntypd),gv(3),Ygaunt(3)!,expf(nvd,natd)
      COMPLEX :: coeff(3,-1:1),gvint(3),vecsum(3),starsum(3)
      REAL   , ALLOCATABLE :: bsl(:,:,:),G(:,:),kG(:,:),kGreal(:,:)
      REAL   , ALLOCATABLE :: kineticfactor(:,:)
      COMPLEX, ALLOCATABLE :: ylm(:,:),ppw(:,:),fpw(:,:),expf(:,:)
      REAL :: zr(nbasfcn,nobd)
      COMPLEX :: zc(nbasfcn,nobd)

      IF (zMat%l_real) THEN
        zr=zMat%data_r 
      ELSE
        zc=zMat%data_c 
      END IF

      lmax = 2*lmaxd
      img = cmplx(0.0,1.0)

      ALLOCATE ( bsl(0:lmax,nvd,ntypd),ylm((lmax+1)**2,nvd) )
      ALLOCATE ( ppw(nobd,(lmax+1)**2),fpw(nobd,(lmax+1)**2) )
      ALLOCATE ( G(3,nvd),kG(3,nvd),kGreal(3,nvd),expf(nvd,natd) )
      ALLOCATE ( kineticfactor(nobd,nvd) )

      coeff(:, :) =   cmplx(0.0,0.0)
      coeff(1,-1) =     sqrt(tpi_const/3.)
      coeff(1, 1) =    -sqrt(tpi_const/3.)
      coeff(2,-1) = img*sqrt(tpi_const/3.)
      coeff(2, 1) = img*sqrt(tpi_const/3.)
      coeff(3, 0) =  sqrt(2.*tpi_const/3.)

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
          CALL sphbes(lmax,r,bsl(:,kn,itype))
          DO ieq = 1,neq(itype)
            expf(kn,iatom) = exp(tpi_const*img*dot_product(kG(:,kn),taual(:,iatom)))
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
          fpil = 2.*tpi_const*img**l
!         Calculate ppw and fpw
          DO kn = 1,nv(jsp)
            DO m = -l,l
              lm = l*(l+1) + m + 1
              noband = expf(kn,iatom)*conjg(ylm(lm,kn))*bsl(l,kn,itype) * fpil
              DO iband = 1,nobd
                IF (zMat%l_real) THEN
                   ppw(iband,lm) = ppw(iband,lm) + zr(kn,iband) * noband
                   IF (l.gt.0) CYCLE
                   fpw(iband,lm) = fpw(iband,lm) + zr(kn,iband) * noband * kineticfactor(iband,kn)
                ELSE
                   ppw(iband,lm) = ppw(iband,lm) + zc(kn,iband) * noband
                   IF (l.gt.0) CYCLE
                   fpw(iband,lm) = fpw(iband,lm) + zc(kn,iband) * noband * kineticfactor(iband,kn)                   
                END IF
              END DO ! iband
            END DO ! m
!           If ppw is used with l, one needs fpw with l-1 (already calculated) and l+1 (calculated now)
            DO m = -l-1,l+1
              lm = (l+1)*(l+2) + m + 1
              noband =expf(kn,iatom)*conjg(ylm(lm,kn))*bsl(l+1,kn,itype) * fpil * img
              DO iband = 1,nobd
                IF (zMat%l_real) THEN
                   fpw(iband,lm) = fpw(iband,lm) + zr(kn,iband) * noband * kineticfactor(iband,kn)
                ELSE
                   fpw(iband,lm) = fpw(iband,lm) + zc(kn,iband) * noband * kineticfactor(iband,kn)
                END IF
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
                Ygaunt(:) = gaunt2(l,1,lp,m,t,mp,lmax)*coeff(:,t)
                DO iband = 1,nobd
                  gv(:) = gv(:) + we(iband) * r2vol * Ygaunt(:) * conjg(ppw(iband,lm)) * fpw(iband,lmp)
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
        vecsum = matmul(bmat,gv)/tpi_const/neq(itype)

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
        force_a12(:,itype) = force_a12(:,itype) + real(starsum(:))/invarind(iatom)

        iatom = iatom + 1
      END DO ! ieq
        force(:,itype,jsp) = force(:,itype,jsp) + force_a12(:,itype)
        f_a12(:,itype)     = f_a12(:,itype)     + force_a12(:,itype)
        WRITE (852,*) (force_a12(i,itype),i=1,3)
      END DO ! itype

      DEALLOCATE ( bsl,ylm,ppw,fpw )
      DEALLOCATE ( G,kG,kGreal,expf,kineticfactor )
!       STOP 'debug'
   END SUBROUTINE force_a12_lv2
END MODULE m_force_a12_lv2
