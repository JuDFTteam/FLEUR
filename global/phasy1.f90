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
      INTEGER j,l,m,n,na,lm,ll1

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

      CALL spgrot(sym%nop, sym%symor, sym%mrot, sym%tau, sym%invtab, &
                  stars%kv3(:,k), kr, phas)

      ALLOCATE ( ylm( (atoms%lmaxd+1)**2, sym%nop ) )
      DO j = 1,sym%nop !center/=0 only works for sym = 1
          rg=matmul(real(kr(:,j))+stars%center,cell%bmat)
          CALL ylm4(atoms%lmaxd, rg, ylm(:,j))!keep
      ENDDO
      ylm = conjg( ylm )

      na = 1
      DO n = 1,atoms%ntype
         DO lm = 1, (atoms%lmax(n)+1)**2
               pylm(lm,n) = cmplx(0.,0.)
         ENDDO
         DO j = 1,sym%nop
            x = tpi_const* dot_product(real(kr(:,j))+stars%center,atoms%taual(:,na))
            sf = cmplx(cos(x),sin(x))*phas(j)
            DO l = 0,atoms%lmax(n)
               ll1 = l*(l+1) + 1
               csf = ciall(l)*sf
               DO m = -l,l
                  lm = ll1 + m
                  pylm(lm,n) = pylm(lm,n) + csf*ylm(lm,j)
               ENDDO
            ENDDO
         ENDDO
         na = na + atoms%neq(n)
      ENDDO
      DEALLOCATE ( ylm )

   END SUBROUTINE phasy1

   SUBROUTINE phasy2(atoms, stars, sym, cell, k, n, na, pylm2)
      ! phasy2 has i*RG in the sum of phasy1 and produces a vector
      ! routine built to be called with a specific atom (type)

!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN)::atoms
      TYPE(t_stars),INTENT(IN)::stars
      TYPE(t_sym),INTENT(IN)  ::sym
      TYPE(t_cell),INTENT(IN) ::cell
      INTEGER, INTENT (IN) :: k, n, na

!     .. Array Arguments ..
      COMPLEX, INTENT (OUT):: pylm2(:,:,:)

!     .. Local Scalars ..
      COMPLEX sf,csf
      REAL x
      INTEGER j,l,m,lm,ll1,dir

!     .. Local Arrays ..
      COMPLEX ciall(0:atoms%lmaxd)
      COMPLEX phas(sym%nop)
      REAL phasr(sym%nop)
      REAL rg(3,sym%nop)
      INTEGER kr(3,sym%nop)
      COMPLEX, ALLOCATABLE :: ylm(:,:)

      ciall(0) = fpi_const/sym%nop
      DO l = 1, atoms%lmax(n)
         ciall(l) = ciall(0)*ImagUnit**l
      ENDDO

      CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),kr,phas)
      phasr=REAL(phas)

      ALLOCATE ( ylm( (atoms%lmaxd+1)**2, sym%nop ) )
      ylm = cmplx(0.,0.)
      DO j = 1,sym%nop
          rg(:,j)=matmul(kr(:,j),cell%bmat)
          CALL ylm4(atoms%lmaxd, rg(:,j), ylm(:,j))!keep
      ENDDO
      ylm = conjg( ylm )

      pylm2 = cmplx(0.,0.)
      DO j = 1, sym%nop
        x = tpi_const* dot_product(real(kr(:,j)),atoms%taual(:,na))
        DO dir = 1,3
          sf = cmplx(cos(x),sin(x))*phasr(j)*ImagUnit*rg(dir,j)
          DO l = 0,atoms%lmax(n)
            ll1 = l*(l+1) + 1
            csf = ciall(l)*sf
            DO m = -l,l
              lm = ll1 + m
              pylm2(lm,dir,j) = pylm2(lm,dir,j) + csf*ylm(lm,j)!shouldn't j be n in the first 2 terms?
            ENDDO
          ENDDO
        END DO ! direction
      ENDDO
      DEALLOCATE ( ylm )

   END SUBROUTINE phasy2
END MODULE m_phasy1
