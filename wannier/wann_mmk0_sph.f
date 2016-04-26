      MODULE m_wann_mmk0_sph
c***********************************************************************
c   computes the Mmn(K) matrix elements which are the overlaps
c   between the Bloch wavefunctions, in the spheres
c   a modification of the eparas.F routine, so go there
c   and to wannier.F for more information on variables
c                                Y.Mokrousov 15.6.06
c***********************************************************************
      CONTAINS
      SUBROUTINE wann_mmk0_sph(
     >                  llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >                  ntype,neq,nlo,llo,acof,bcof,ccof,
     >                  ddn,uulon,dulon,uloulopn,
     =                  mmn)
      implicit none
c     .. scalar arguments ..
      integer, intent (in) :: llod,nlod,natd,ntypd,lmaxd,lmd
      integer, intent (in) :: ntype,noccbd
c     .. array arguments ..
      integer, intent (in)  :: neq(ntypd)
      integer, intent (in)  :: nlo(ntypd),llo(nlod,ntypd)
      real,    intent (in)  :: ddn(0:lmaxd,ntypd)
      real,    intent (in)  :: uloulopn(nlod,nlod,ntypd)
      real,    intent (in)  :: uulon(nlod,ntypd),dulon(nlod,ntypd)
      complex, intent (in)  :: ccof(-llod:llod,noccbd,nlod,natd)
      complex, intent (in)  :: acof(:,0:,:) !acof(noccbd,0:lmd,natd)
      complex, intent (in)  :: bcof(:,0:,:) !bcof(noccbd,0:lmd,natd)
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
      allocate (qlo(noccbd,noccbd,nlod,nlod,ntypd), 
     +          qaclo(noccbd,noccbd,nlod,ntypd),
     +          qbclo(noccbd,noccbd,nlod,ntypd) )
c---> performs summations of the overlaps of the wavefunctions
      do 140 i = 1,noccbd            
       do 145 j = 1,noccbd
         nt1 = 1
         do 130 n = 1,ntype
            nt2 = nt1 + neq(n) - 1
            do 120 l = 0,lmaxd
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
            nt1 = nt1 + neq(n)
  130    continue
  145  continue   ! cycle by j-band
  140 continue  !  cycle by i-band
c---> initialize qlo arrays
      qlo(:,:,:,:,:) = 0.0
      qaclo(:,:,:,:) = 0.0
      qbclo(:,:,:,:) = 0.0
c---> prepare the coefficients
      natom = 0
      do ntyp = 1,ntype
         do nn = 1,neq(ntyp)
            natom = natom + 1
            do lo = 1,nlo(ntyp)
               l = llo(lo,ntyp)
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
               do lop = 1,nlo(ntyp)
                 if (llo(lop,ntyp).eq.l) then
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
      do ntyp = 1,ntype
         do lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            do i = 1,noccbd
             do j = 1,noccbd
               mmn(i,j)= mmn(i,j)  + 
     +                      ( qaclo(i,j,lo,ntyp)*uulon(lo,ntyp)     +
     +                        qbclo(i,j,lo,ntyp)*dulon(lo,ntyp)     )
             enddo
            enddo 
            do lop = 1,nlo(ntyp)
               if (llo(lop,ntyp).eq.l) then
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

      END SUBROUTINE wann_mmk0_sph
      END MODULE m_wann_mmk0_sph
