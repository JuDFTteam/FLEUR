      MODULE m_wann_anglmom
c***********************************************************************
c     Compute matrix elements of angular momentum operator 
c     in the muffin-tin spheres.
c
c     Frank Freimuth
c***********************************************************************
      CONTAINS
      SUBROUTINE wann_anglmom(
     >                  llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >                  ntype,neq,nlo,llo,acof,bcof,ccof,
     >                  ddn,uulon,dulon,uloulopn,
     =                  mmn)
      implicit none
c     .. scalar arguments ..
      integer, intent (in) :: llod,nlod,natd,ntypd,lmaxd,lmd
      integer, intent (in) :: ntype,noccbd
c     .. array arguments ..
      integer, intent (in)  :: neq(:)!neq(ntypd)
      integer, intent (in)  :: nlo(:)!nlo(ntypd)
      integer, intent (in)  :: llo(:,:)!llo(nlod,ntypd)
      real,    intent (in)  :: ddn(0:,:)!ddn(0:lmaxd,ntypd)
      real,    intent (in)  :: uloulopn(:,:,:)!uloulopn(nlod,nlod,ntypd)
      real,    intent (in)  :: uulon(:,:)!uulon(nlod,ntypd)
      real,    intent (in)  :: dulon(:,:)!dulon(nlod,ntypd)
      complex, intent (in)  :: ccof(-llod:,:,:,:) !ccof(-llod:llod,noccbd,nlod,natd)
      complex, intent (in)  :: acof(:,0:,:)!acof(noccbd,0:lmd,natd)
      complex, intent (in)  :: bcof(:,0:,:)!bcof(noccbd,0:lmd,natd)
      complex, intent (inout) :: mmn(:,:,:)!mmn(3,noccbd,noccbd)
c     .. local scalars ..
      integer :: i,j,l,lo,lop,m,natom,nn,ntyp
      integer :: nt1,nt2,lm,n,ll1
      complex :: suma_z,sumb_z
      complex :: suma_p,sumb_p
      complex :: suma_m,sumb_m
      complex :: suma_x,sumb_x
      complex :: suma_y,sumb_y
      real    :: lplus,lminus
C     ..
C     .. local arrays ..
      complex, allocatable :: qlo_z(:,:,:,:,:)
      complex, allocatable :: qlo_p(:,:,:,:,:)
      complex, allocatable :: qlo_m(:,:,:,:,:)

      complex, allocatable :: qaclo_z(:,:,:,:),qbclo_z(:,:,:,:)
      complex, allocatable :: qaclo_p(:,:,:,:),qbclo_p(:,:,:,:)
      complex, allocatable :: qaclo_m(:,:,:,:),qbclo_m(:,:,:,:)
C     ..
C     .. intrinsic functions ..
      intrinsic conjg
      allocate (qlo_z(noccbd,noccbd,nlod,nlod,ntypd), 
     +          qaclo_z(noccbd,noccbd,nlod,ntypd),
     +          qbclo_z(noccbd,noccbd,nlod,ntypd) )

      allocate (qlo_p(noccbd,noccbd,nlod,nlod,ntypd), 
     +          qaclo_p(noccbd,noccbd,nlod,ntypd),
     +          qbclo_p(noccbd,noccbd,nlod,ntypd) )

      allocate (qlo_m(noccbd,noccbd,nlod,nlod,ntypd), 
     +          qaclo_m(noccbd,noccbd,nlod,ntypd),
     +          qbclo_m(noccbd,noccbd,nlod,ntypd) )

c-----> lapw-lapw-Terms
      do i = 1,noccbd            
       do j = 1,noccbd
         nt1 = 1
         do n = 1,ntype
            nt2 = nt1 + neq(n) - 1
            do l = 0,lmaxd
               suma_z = cmplx(0.,0.); sumb_z = cmplx(0.,0.)
               suma_m = cmplx(0.,0.); sumb_m = cmplx(0.,0.)
               suma_p = cmplx(0.,0.); sumb_p = cmplx(0.,0.)
               ll1 = l* (l+1)
               do m = -l,l
                  lm = ll1 + m
                  lplus=sqrt(real( (l-m)*(l+m+1) ) )
                  lminus=sqrt(real( (l+m)*(l-m+1) ) )
                  do natom = nt1,nt2
                    suma_z = suma_z + acof(i,lm,natom)*
     +                          conjg(acof(j,lm,natom))*real(m)
                    sumb_z = sumb_z + bcof(i,lm,natom)*
     +                          conjg(bcof(j,lm,natom))*real(m)
                    if(m+1.le.l)then
                     suma_p = suma_p + acof(i,lm,natom)*
     +                          conjg(acof(j,lm+1,natom))*lplus
                     sumb_p = sumb_p + bcof(i,lm,natom)*
     +                          conjg(bcof(j,lm+1,natom))*lplus
                    endif
                    if(m-1.ge.-l)then
                     suma_m = suma_m + acof(i,lm,natom)*
     +                          conjg(acof(j,lm-1,natom))*lminus
                     sumb_m = sumb_m + bcof(i,lm,natom)*
     +                          conjg(bcof(j,lm-1,natom))*lminus
                    endif
                  enddo
               enddo
               mmn(3,j,i) = mmn(3,j,i) + (suma_z+sumb_z*ddn(l,n))

               suma_x=0.5*(suma_p+suma_m)
               sumb_x=0.5*(sumb_p+sumb_m)
               mmn(1,j,i) = mmn(1,j,i) + (suma_x+sumb_x*ddn(l,n))

               suma_y=cmplx(0.0,-0.5)*(suma_p-suma_m)
               sumb_y=cmplx(0.0,-0.5)*(sumb_p-sumb_m)
               mmn(2,j,i) = mmn(2,j,i) + (suma_y+sumb_y*ddn(l,n))
            enddo ! l
            nt1 = nt1 + neq(n)
         enddo ! n
       enddo ! j
      enddo ! i


c---> Terms involving local orbitals.
      qlo_z = 0.0; qlo_p = 0.0; qlo_m = 0.0
      qaclo_z = 0.0; qaclo_p = 0.0; qaclo_m = 0.0
      qbclo_z = 0.0; qbclo_p = 0.0; qbclo_m = 0.0

      natom = 0
      do ntyp = 1,ntype
       do nn = 1,neq(ntyp)
         natom = natom + 1
         do lo = 1,nlo(ntyp)
           l = llo(lo,ntyp)
           ll1 = l* (l+1)
           do m = -l,l
            lm = ll1 + m
            lplus=sqrt(real( (l-m)*(l+m+1) ) )
            lminus=sqrt(real( (l+m)*(l-m+1) ) )
            do i = 1,noccbd
             do j = 1,noccbd
                qbclo_z(j,i,lo,ntyp) = qbclo_z(j,i,lo,ntyp) + (
     +         bcof(i,lm,natom) * conjg(ccof(m,j,lo,natom)) +
     +         ccof(m,i,lo,natom)*conjg(bcof(j,lm,natom)) )*real(m)

                qaclo_z(j,i,lo,ntyp) = qaclo_z(j,i,lo,ntyp) + (
     +         acof(i,lm,natom) * conjg(ccof(m,j,lo,natom)) +
     +         ccof(m,i,lo,natom)*conjg(acof(j,lm,natom)) )*real(m)
                if(m+1.le.l)then
                 qbclo_p(j,i,lo,ntyp) = qbclo_p(j,i,lo,ntyp) + (
     +           bcof(i,lm,natom) * conjg(ccof(m+1,j,lo,natom)) +
     +           ccof(m,i,lo,natom)*conjg(bcof(j,lm+1,natom)) )*lplus

                 qaclo_p(j,i,lo,ntyp) = qaclo_p(j,i,lo,ntyp) + (
     +           acof(i,lm,natom) * conjg(ccof(m+1,j,lo,natom)) +
     +           ccof(m,i,lo,natom)*conjg(acof(j,lm+1,natom)) )*lplus
                endif
                if(m-1.ge.-l)then
                 qbclo_m(j,i,lo,ntyp) = qbclo_m(j,i,lo,ntyp) + (
     +           bcof(i,lm,natom) * conjg(ccof(m-1,j,lo,natom)) +
     +           ccof(m,i,lo,natom)*conjg(bcof(j,lm-1,natom)) )*lminus

                 qaclo_m(j,i,lo,ntyp) = qaclo_m(j,i,lo,ntyp) + (
     +           acof(i,lm,natom) * conjg(ccof(m-1,j,lo,natom)) +
     +           ccof(m,i,lo,natom)*conjg(acof(j,lm-1,natom)) )*lminus
                endif

             enddo !j
            enddo !i
           enddo !m
           do lop = 1,nlo(ntyp)
              if (llo(lop,ntyp).eq.l) then
                do m = -l,l
                  lplus=sqrt(real( (l-m)*(l+m+1) ) )
                  lminus=sqrt(real( (l+m)*(l-m+1) ) )
                  do i = 1,noccbd
                   do j = 1,noccbd
                    qlo_z(j,i,lop,lo,ntyp) = qlo_z(j,i,lop,lo,ntyp) + 
     +                    conjg(ccof(m,j,lop,natom))
     *                               *ccof(m,i,lo,natom)*real(m)
                    if(m+1.le.l)then
                       qlo_p(j,i,lop,lo,ntyp) = 
     +                   qlo_p(j,i,lop,lo,ntyp) + 
     +                    conjg(ccof(m+1,j,lop,natom))
     *                         *ccof(m,i,lo,natom)*lplus

                    endif
                    if(m-1.ge.-l)then
                       qlo_m(j,i,lop,lo,ntyp) = 
     +                   qlo_m(j,i,lop,lo,ntyp) + 
     +                    conjg(ccof(m-1,j,lop,natom))
     *                         *ccof(m,i,lo,natom)*lminus
                    endif
                   enddo ! j
                  enddo ! i
                enddo ! m
              endif
           enddo ! lop
         enddo ! lo
       enddo ! nn
      enddo ! ntyp
c---> perform summation of the coefficients with the integrals
c---> of the radial basis functions
      do ntyp = 1,ntype
         do lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            do j = 1,noccbd
             do i = 1,noccbd
               mmn(3,i,j)= mmn(3,i,j)  + 
     +               qaclo_z(i,j,lo,ntyp)*uulon(lo,ntyp) +
     +               qbclo_z(i,j,lo,ntyp)*dulon(lo,ntyp)  

               suma_p=qaclo_p(i,j,lo,ntyp)*uulon(lo,ntyp) +
     +                qbclo_p(i,j,lo,ntyp)*dulon(lo,ntyp)

               suma_m=qaclo_m(i,j,lo,ntyp)*uulon(lo,ntyp) +
     +                qbclo_m(i,j,lo,ntyp)*dulon(lo,ntyp)

               suma_x=            0.5*(suma_p+suma_m)
               suma_y=cmplx(0.0,-0.5)*(suma_p-suma_m)

               mmn(1,i,j)= mmn(1,i,j)  + suma_x
               mmn(2,i,j)= mmn(2,i,j)  + suma_y 

             enddo !i
            enddo !j 
            do lop = 1,nlo(ntyp)
              if (llo(lop,ntyp).eq.l) then
               do j = 1,noccbd
                do i = 1,noccbd
                 mmn(3,i,j) = mmn(3,i,j)  + 
     +                  qlo_z(i,j,lop,lo,ntyp)*uloulopn(lop,lo,ntyp)
                 suma_p=qlo_p(i,j,lop,lo,ntyp)*uloulopn(lop,lo,ntyp)
                 suma_m=qlo_m(i,j,lop,lo,ntyp)*uloulopn(lop,lo,ntyp)
                 mmn(1,i,j) = mmn(1,i,j) + 0.5*(suma_p+suma_m)
                 mmn(2,i,j) = mmn(2,i,j) + 
     +                 cmplx(0.0,-0.5)*(suma_p-suma_m)
                enddo ! i
               enddo ! j
              endif
            enddo !lop
         enddo !lo 
      enddo !ntyp 
      deallocate ( qlo_z,qaclo_z,qbclo_z )
      deallocate ( qlo_m,qaclo_m,qbclo_m )
      deallocate ( qlo_p,qaclo_p,qbclo_p )

      END SUBROUTINE wann_anglmom
      END MODULE m_wann_anglmom
