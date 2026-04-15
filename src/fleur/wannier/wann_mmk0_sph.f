!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_mmk0_sph
c***********************************************************************
c   computes the Mmn(K) matrix elements which are the overlaps
c   between the Bloch wavefunctions, in the spheres
c   a modification of the eparas.F routine, so go there
c   and to wannier.F for more information on variables
c                                Y.Mokrousov 15.6.06
c***********************************************************************
      use m_juDFT
      CONTAINS
      SUBROUTINE wann_mmk0_sph(
     >                  atoms,noccbd,lmd,acof,bcof,ccof,
     >                  ddn,uulon,dulon,uloulopn,
     =                  mmn)
      use m_types
      implicit none
      TYPE(t_atoms), INTENT(IN) :: atoms
c     .. scalar arguments ..
      integer, intent (in) :: lmd
      integer, intent (in) :: noccbd
c     .. array arguments ..
      real,    intent (in)  :: ddn(0:,:)
      real,    intent (in)  :: uloulopn(:,:,:)
      real,    intent (in)  :: uulon(:,:),dulon(:,:)
      complex, intent (in)  ::
     >   ccof(-atoms%llod:,:,:,:)
      complex, intent (in)  :: acof(:,0:,:)
      complex, intent (in)  :: bcof(:,0:,:)
      complex, intent (inout) :: mmn(:,:)
c     .. local scalars ..
      integer i,j,l,lo,lop,m,natom,nn,ntyp
      integer nt1,nt2,lm,n,ll1
      complex suma,sumb
C     ..
C     .. local arrays ..
      complex, allocatable :: qlo(:,:,:,:,:)
      complex, allocatable :: qaclo(:,:,:,:),qbclo(:,:,:,:)
C     ..
C     .. intrinsic functions ..
      intrinsic conjg

      call timestart("wann_mmk0_sph")

      allocate (qlo(noccbd,noccbd,atoms%nlod,atoms%nlod,atoms%ntype), 
     +          qaclo(noccbd,noccbd,atoms%nlod,atoms%ntype),
     +          qbclo(noccbd,noccbd,atoms%nlod,atoms%ntype) )
c---> performs summations of the overlaps of the wavefunctions
      do 140 i = 1,noccbd            
       do 145 j = 1,noccbd
         nt1 = 1
         do 130 n = 1,atoms%ntype
            nt2 = nt1 + atoms%neq(n) - 1
            do 120 l = 0,atoms%lmax(n)
               suma = cmplx(0.,0.)
               sumb = cmplx(0.,0.)
               ll1 = l* (l+1)
               do 110 m = -l,l
                  lm = ll1 + m
                  do natom = nt1,nt2
                    suma = suma + acof(i,lm,natom)*
     +                     conjg(acof(j,lm,natom))
                    sumb = sumb + bcof(i,lm,natom)*
     +                     conjg(bcof(j,lm,natom))
                  enddo
 110          continue
               mmn(i,j) = mmn(i,j) + (suma+sumb*ddn(l,n))
  120       continue
            nt1 = nt1 + atoms%neq(n)
  130    continue
  145  continue   ! cycle by j-band
  140 continue  !  cycle by i-band
c---> initialize qlo arrays
      qlo(:,:,:,:,:) = 0.0
      qaclo(:,:,:,:) = 0.0
      qbclo(:,:,:,:) = 0.0
c---> prepare the coefficients
      natom = 0
      do ntyp = 1,atoms%ntype
         do nn = 1,atoms%neq(ntyp)
            natom = natom + 1
            do lo = 1,atoms%nlo(ntyp)
               l = atoms%llo(lo,ntyp)
               ll1 = l* (l+1)
               do m = -l,l
                  lm = ll1 + m
                  do i = 1,noccbd
                   do j = 1,noccbd
                     qbclo(i,j,lo,ntyp) = qbclo(i,j,lo,ntyp) + 
     +                      bcof(i,lm,natom)*conjg(ccof(m,j,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(bcof(j,lm,natom)) 
                     qaclo(i,j,lo,ntyp) = qaclo(i,j,lo,ntyp) + 
     +                      acof(i,lm,natom)*conjg(ccof(m,j,lo,natom)) +
     +                      ccof(m,i,lo,natom)*conjg(acof(j,lm,natom)) 
                   enddo
                  enddo
               enddo
               do lop = 1,atoms%nlo(ntyp)
                 if (atoms%llo(lop,ntyp).eq.l) then
                   do m = -l,l
                     do i = 1,noccbd
                      do j = 1,noccbd
                       qlo(i,j,lop,lo,ntyp) = qlo(i,j,lop,lo,ntyp) + 
     +                        conjg(ccof(m,j,lop,natom))
     *                                  *ccof(m,i,lo,natom)
                      enddo
                     enddo
                   enddo
                 endif
               enddo
            enddo
         enddo
      enddo

c---> perform summation of the coefficients with the integrals
c---> of the radial basis functions
      do ntyp = 1,atoms%ntype
         do lo = 1,atoms%nlo(ntyp)
            l = atoms%llo(lo,ntyp)
            do i = 1,noccbd
             do j = 1,noccbd
               mmn(i,j)= mmn(i,j)  + 
     +                      ( qaclo(i,j,lo,ntyp)*uulon(lo,ntyp)     +
     +                        qbclo(i,j,lo,ntyp)*dulon(lo,ntyp)     )
             enddo
            enddo 
            do lop = 1,atoms%nlo(ntyp)
               if (atoms%llo(lop,ntyp).eq.l) then
               do i = 1,noccbd
                do j = 1,noccbd
                 mmn(i,j) = mmn(i,j)  + 
     +                      qlo(i,j,lop,lo,ntyp)*uloulopn(lop,lo,ntyp)
                enddo
               enddo
               endif
            enddo
         enddo 
      enddo 
      deallocate ( qlo,qaclo,qbclo )

      call timestop("wann_mmk0_sph")
      END SUBROUTINE wann_mmk0_sph
      END MODULE m_wann_mmk0_sph
