!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_phasy1
   !-----------------------------------------------------------------------------
   ! Calculate 4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) }
   !     e. wimmer   oct.1984
   !-----------------------------------------------------------------------------
   USE m_constants
   USE m_ylm
   USE m_spgrot
   USE m_types
   implicit none

   private
   public :: phasy1, phasy2, phasy1nSym

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
      REAL rg(3)
      INTEGER kr(3,sym%nop)
      COMPLEX, ALLOCATABLE :: ylm(:)

      ciall(0) = fpi_const/sym%nop
      DO l = 1, atoms%lmax(iType)
         ciall(l) = ciall(0)*ImagUnit**l
      ENDDO

      pylm2= CMPLX(0.0,0.0)

      CALL spgrot(sym%nop,sym%symor,sym%mrot,sym%tau,sym%invtab,stars%kv3(:,k),kr,phas)

      ALLOCATE (ylm( (atoms%lmaxd+1)**2))
      DO iOp = 1,sym%nop
         ylm = cmplx(0.0,0.0)
         rg(:)=matmul(kr(:,iOp),cell%bmat)
         CALL ylm4(atoms%lmaxd, rg(:), ylm(:))!keep
         ylm = conjg(ylm)
         x = tpi_const* dot_product(real(kr(:,iOp)),atoms%taual(:,iAtom))
         DO dir = 1,3
            sf = cmplx(cos(x),sin(x))*phas(iOp)*ImagUnit*rg(dir)
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

   subroutine phasy1nSym(atoms, cell, Gvec, qptn, pylm)
      !Routine by C. Gerhorst to calculate phasefactors for dfpt
      use m_ylm
      use m_types_atoms
      use m_types_cell
    
      implicit none
    
      ! Scalar Type Arguments
      type(t_atoms),  intent(in)  ::  atoms
      type(t_cell),   intent(in)  ::  cell
    
      ! Array Arguments
      integer,        intent(in)  ::  Gvec(:)
      real,           intent(in)  ::  qptn(:)
      complex,        intent(out) ::  pylm(:, :)
    
      !-------------------------------------------------------------------------------------------------------------------------------
      ! Local Scalar Variables
      ! iatom : runs over all atoms
      ! itype : runs over all types
      ! ineq  : runs over all equivalent atoms of one atom type
      ! lm    : encodes oqn_l and mqn_m
      ! sf    : stores exponential function
      ! csf   : stores exponential function times 4 π i^l
      ! x     : stores argument of exponential function
      ! mqn_m : magnetic quantum number m
      ! oqn_l : orbital quantum number l
      ! ll1   : auxiliary variable to calculate lm
      !-------------------------------------------------------------------------------------------------------------------------------
      integer                     ::  iatom
      integer                     ::  itype
      integer                     ::  ineq
      integer                     ::  lm
      complex                     ::  sf
      complex                     ::  csf
      real                        ::  x
      integer                     ::  mqn_m
      integer                     ::  oqn_l
      integer                     ::  ll1
    
      !-------------------------------------------------------------------------------------------------------------------------------
      ! Local Array Variables
      ! fpiul: stores 4 π i^l
      ! Gqext: stores G + q in external coordinates
      ! ylm  : stores Y_lm
      !-------------------------------------------------------------------------------------------------------------------------------
      complex, allocatable        ::  fpiul(:)
      complex, allocatable        ::  ylm(:)
      real                        ::  Gqext(3)
    
      allocate(fpiul(0:atoms%lmaxd))
      fpiul(:) = cmplx(0., 0.)
      ! calculates 4 π i^l resolved for every l, not divided by nop because no loop over symmetry operations
      do oqn_l = 0, atoms%lmaxd
         fpiul(oqn_l) = fpi_const * ImagUnit**oqn_l
      enddo
    
    
      ! calculates Y*_lm(\vec{G} + \vec{q}) for every l and m. The argument Gqext must be in external coordinates.
      allocate(ylm((atoms%lmaxd + 1)**2))
      ylm(:) = cmplx(0., 0.)
      Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), real(Gvec(1:3) + qptn(1:3)))
      call ylm4(atoms%lmaxd, Gqext(1:3), ylm)
      ylm = conjg(ylm)
    
    
      ! calculates first exp(i (G + q) tau)  and multiplies recent factors before storing the final result to pylm
      iatom = 1
      pylm(:, :) = cmplx(0.,0.)
      do itype = 1, atoms%ntype
         do ineq = 1, atoms%neq(itype)
            x = tpi_const * dot_product(real(Gvec(1:3) + qptn(1:3)), atoms%taual(1:3, iatom))
            sf = exp(ImagUnit *  x)
            do oqn_l = 0, atoms%lmax(itype)
               ll1 = oqn_l * (oqn_l + 1) + 1
               csf = fpiul(oqn_l) * sf
               do mqn_m = -oqn_l, oqn_l
                  lm = ll1 + mqn_m
                  pylm(lm, iatom) = csf * ylm(lm)
               enddo ! mqn_m
            enddo ! oqn_l
            iatom = iatom + 1
         enddo ! ineq
      enddo ! itype
    
    end subroutine phasy1nSym
END MODULE m_phasy1
