MODULE m_phasy1
   !-----------------------------------------------------------------------------
   ! Calculate 4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) }
   !     e. wimmer   oct.1984
   !-----------------------------------------------------------------------------
   USE m_constants
   USE m_ylm
   USE m_spgrot
   USE m_types

   IMPLICIT NONE

CONTAINS
   SUBROUTINE phasy1(atoms,stars,sym, cell,k, pylm)

!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN)::atoms
      TYPE(t_stars),INTENT(IN)::stars
      TYPE(t_sym),INTENT(IN)  ::sym
      TYPE(t_cell),INTENT(IN) ::cell
      INTEGER, INTENT (IN) :: k

!     .. Array Arguments ..
      COMPLEX, INTENT (OUT):: pylm(:,:)

!     .. Local Scalars ..
      COMPLEX sf,csf
      REAL x
      INTEGER iOp,l,m,iType,iAtom,lm,ll1

!     .. Local Arrays ..
      COMPLEX ciall(0:atoms%lmaxd)
      COMPLEX phas(sym%nop)
      REAL rg(3)
      INTEGER kr(3,sym%nop)
      COMPLEX, ALLOCATABLE :: ylm(:,:)

      ciall(0) = fpi_const/sym%nop
      DO l = 1,atoms%lmaxd
         ciall(l) = ciall(0)*ImagUnit**l
      ENDDO

      pylm = CMPLX(0.0,0.0)

      CALL spgrot(sym%nop, sym%symor, sym%mrot, sym%tau, sym%invtab, &
                  stars%kv3(:,k), kr, phas)

      ALLOCATE ( ylm( (atoms%lmaxd+1)**2, sym%nop ) )
      DO iOp = 1,sym%nop !center/=0 only works for sym = 1
          rg=matmul(real(kr(:,iOp))+stars%center,cell%bmat)
          CALL ylm4(atoms%lmaxd, rg, ylm(:,iOp))!keep
      ENDDO
      ylm = conjg( ylm )

      DO iType = 1,atoms%ntype
         iAtom = atoms%firstAtom(iType)
         DO iOp = 1,sym%nop
            x = tpi_const* dot_product(real(kr(:,iOp))+stars%center,atoms%taual(:,iAtom))
            sf = cmplx(cos(x),sin(x))*phas(iOp)
            DO l = 0,atoms%lmax(iType)
               ll1 = l*(l+1) + 1
               csf = ciall(l)*sf
               DO m = -l,l
                  lm = ll1 + m
                  pylm(lm,iType) = pylm(lm,iType) + csf*ylm(lm,iOp)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE ( ylm )

   END SUBROUTINE phasy1

   SUBROUTINE phasy2(atoms, stars, sym, cell, k, iType, iAtom, pylm2)
      ! phasy2 has i*RG in the sum of phasy1 and produces a vector
      ! routine built to be called with a specific atom (type)

!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN)::atoms
      TYPE(t_stars),INTENT(IN)::stars
      TYPE(t_sym),INTENT(IN)  ::sym
      TYPE(t_cell),INTENT(IN) ::cell
      INTEGER, INTENT (IN) :: k, iType, iAtom

!     .. Array Arguments ..
      COMPLEX, INTENT (OUT):: pylm2(:,:,:)

!     .. Local Scalars ..
      COMPLEX sf,csf
      REAL x
      INTEGER iOp,l,m,lm,ll1,dir

!     .. Local Arrays ..
      COMPLEX ciall(0:atoms%lmaxd)
      COMPLEX phas(sym%nop)
      REAL phasr(sym%nop)
      REAL rg(3)
      INTEGER kr(3,sym%nop)
      COMPLEX, ALLOCATABLE :: ylm(:)

      ciall(0) = fpi_const/sym%nop
      DO l = 1, atoms%lmax(iType)
         ciall(l) = ciall(0)*ImagUnit**l
      ENDDO

      pylm2= CMPLX(0.0,0.0)

      CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),kr,phas)
      phasr=REAL(phas)

      ALLOCATE (ylm( (atoms%lmaxd+1)**2))
      DO iOp = 1,sym%nop
         ylm = cmplx(0.0,0.0)
         rg(:)=matmul(kr(:,iOp),cell%bmat)
         CALL ylm4(atoms%lmaxd, rg(:), ylm(:))!keep
         ylm = conjg(ylm)
         x = tpi_const* dot_product(real(kr(:,iOp)),atoms%taual(:,iAtom))
         DO dir = 1,3
            sf = cmplx(cos(x),sin(x))*phasr(iOp)*ImagUnit*rg(dir)
            DO l = 0,atoms%lmax(iType)
               ll1 = l*(l+1) + 1
               csf = ciall(l)*sf
               DO m = -l,l
                  lm = ll1 + m
                  pylm2(lm,dir,iOp) = pylm2(lm,dir,iOp) + csf*ylm(lm) !shouldn't iOp be iType in the first 2 terms?
               ENDDO
            ENDDO
         END DO ! direction
      END DO
      DEALLOCATE ( ylm )

   END SUBROUTINE phasy2
END MODULE m_phasy1
