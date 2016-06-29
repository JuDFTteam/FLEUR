!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_mmk0_updown_sph
      CONTAINS
      SUBROUTINE wann_mmk0_updown_sph(
     >               l_noco,alph,beta,
     >               llod,noccbd,nlod,natd,ntypd,lmaxd,lmd,
     >               ntype,neq,nlo,llo,
     >               radial1_ff,radial1_gg,
     >               radial1_fg,radial1_gf,
     >               radial1_flo,radial1_glo,
     >               radial1_lof,radial1_log, 
     >               radial1_lolo,
     >               acof,bcof,ccof,
     >               ddn,uulon,dulon,uloulopn,
     =               mmn)
c************************************************************
c     Overlaps of the spin-down parts of the Bloch functions
c     with the spin-up parts.
c                           Frank Freimuth
c************************************************************
      implicit none
      logical, intent (in)  :: l_noco
      integer, intent (in)  :: llod,nlod,natd,ntypd,lmaxd,lmd
      integer, intent (in)  :: ntype,noccbd
      REAL,    INTENT (IN)  :: alph(ntypd),beta(ntypd)
      integer, intent (in)  :: neq(ntypd)
      integer, intent (in)  :: nlo(ntypd),llo(nlod,ntypd)
      real,    intent (in)  :: radial1_ff(:,:,0:,0:,:)
      real,    intent (in)  :: radial1_gg(:,:,0:,0:,:)      
      real,    intent (in)  :: radial1_fg(:,:,0:,0:,:)
      real,    intent (in)  :: radial1_gf(:,:,0:,0:,:)
      real,intent(in)       :: radial1_flo(:,:,0:,:,:)
      real,intent(in)       :: radial1_glo(:,:,0:,:,:)
      real,intent(in)       :: radial1_lof(:,:,:,0:,:)
      real,intent(in)       :: radial1_log(:,:,:,0:,:)
      real,intent(in)       :: radial1_lolo(:,:,:,:,:)
      real,    intent (in)  :: ddn(0:lmaxd,ntypd,2)
      real,    intent (in)  :: uloulopn(nlod,nlod,ntypd,2)
      real,    intent (in)  :: uulon(nlod,ntypd,2),dulon(nlod,ntypd,2)
      complex, intent (in)  :: ccof(-llod:llod,noccbd,nlod,natd,2)
      complex, intent (in)  :: acof(noccbd,0:lmd,natd,2)
      complex, intent (in)  :: bcof(noccbd,0:lmd,natd,2)
      complex, intent (inout) :: mmn(noccbd,noccbd)

      integer           :: i,j,l,lo,lop,m,natom,nn,ntyp
      integer           :: nt1,nt2,lm,n,ll1,i1spin,i2spin
      complex           :: suma,sumb,sumc,sumd
      complex           :: suma12(2,2),sumb12(2,2)
      complex           :: sumc12(2,2),sumd12(2,2)
      real, allocatable :: qlo(:,:,:,:,:)
      real, allocatable :: qaclo(:,:,:,:),qbclo(:,:,:,:)
      COMPLEX           :: ccchi(2,2),ci

      ci = cmplx(0.0,1.0)
      allocate (qlo(noccbd,noccbd,nlod,nlod,ntypd), 
     +          qaclo(noccbd,noccbd,nlod,ntypd),
     +          qbclo(noccbd,noccbd,nlod,ntypd) )
c---> performs summations of the overlaps of the wavefunctions
      do i = 1,noccbd            
       do j = 1,noccbd
         nt1 = 1
         do n = 1,ntype
            if(l_noco)then
               ccchi(1,1) = conjg( exp( ci*alph(n)/2)*cos(beta(n)/2))
               ccchi(1,2) = conjg(-exp( ci*alph(n)/2)*sin(beta(n)/2))
               ccchi(2,1) = conjg( exp(-ci*alph(n)/2)*sin(beta(n)/2))
               ccchi(2,2) = conjg( exp(-ci*alph(n)/2)*cos(beta(n)/2))
            endif
            nt2 = nt1 + neq(n) - 1
            do l = 0,lmaxd
             if(.not.l_noco)then  
               suma = cmplx(0.,0.)
               sumb = cmplx(0.,0.)
               sumc = cmplx(0.,0.)
               sumd = cmplx(0.,0.)
               ll1 = l* (l+1)
               do m = -l,l
                  lm = ll1 + m
                  do natom = nt1,nt2
                    suma = suma + acof(i,lm,natom,1)*
     +                      conjg(acof(j,lm,natom,2))
                    sumb = sumb + bcof(i,lm,natom,1)*
     +                      conjg(bcof(j,lm,natom,2))
                    sumc = sumc + acof(i,lm,natom,1)*
     +                      conjg(bcof(j,lm,natom,2))
                    sumd = sumd + bcof(i,lm,natom,1)*
     +                      conjg(acof(j,lm,natom,2))
                  enddo !natom
               enddo !m      
               mmn(i,j) = mmn(i,j) + ( suma*radial1_ff(1,2,l,l,n)+
     +                                 sumb*radial1_gg(1,2,l,l,n)+
     +                                 sumc*radial1_fg(1,2,l,l,n)+
     +                                 sumd*radial1_gf(1,2,l,l,n)  )      
             else
               suma12 = cmplx(0.,0.)
               sumb12 = cmplx(0.,0.)
               sumc12 = cmplx(0.,0.)
               sumd12 = cmplx(0.,0.)
               ll1 = l* (l+1)
               do i1spin=1,2
                do i2spin=1,2
                 do m = -l,l
                  lm = ll1 + m
                  do natom = nt1,nt2
                    suma12(i1spin,i2spin) = suma12(i1spin,i2spin) 
     +                      + acof(i,lm,natom,i1spin)*
     +                      conjg(acof(j,lm,natom,i2spin))
                    sumb12(i1spin,i2spin) = sumb12(i1spin,i2spin) 
     +                      + bcof(i,lm,natom,i1spin)*
     +                      conjg(bcof(j,lm,natom,i2spin))
                    sumc12(i1spin,i2spin) = sumc12(i1spin,i2spin) 
     +                      + acof(i,lm,natom,i1spin)*
     +                      conjg(bcof(j,lm,natom,i2spin))
                    sumd12(i1spin,i2spin) = sumd12(i1spin,i2spin) 
     +                      + bcof(i,lm,natom,i1spin)*
     +                      conjg(acof(j,lm,natom,i2spin))
                  enddo !natom
                 enddo !m
                 mmn(i,j) =    mmn(i,j)
     &                           +
     &             suma12(i1spin,i2spin)*radial1_ff(i1spin,i2spin,l,l,n)
     &                    *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &             sumb12(i1spin,i2spin)*radial1_gg(i1spin,i2spin,l,l,n)
     &                    *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &             sumc12(i1spin,i2spin)*radial1_fg(i1spin,i2spin,l,l,n)
     &                    *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &             sumd12(i1spin,i2spin)*radial1_gf(i1spin,i2spin,l,l,n)
     &                    *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))
                enddo !i2spin
               enddo !i1spin 
             endif   

            enddo !l
            nt1 = nt1 + neq(n)
         enddo !n   
       enddo   ! cycle by j-band
      enddo  !  cycle by i-band

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
                     qbclo(i,j,lo,ntyp) = qbclo(i,j,lo,ntyp) + real(
     +                  bcof(i,lm,natom,1)*conjg(ccof(m,j,lo,natom,2)) +
     +                  ccof(m,i,lo,natom,1)*conjg(bcof(j,lm,natom,2)) )
                     qaclo(i,j,lo,ntyp) = qaclo(i,j,lo,ntyp) + real(
     +                  acof(i,lm,natom,1)*conjg(ccof(m,j,lo,natom,2)) +
     +                  ccof(m,i,lo,natom,1)*conjg(acof(j,lm,natom,2)) )
                   enddo
                  enddo
               enddo
               do lop = 1,nlo(ntyp)
                 if (llo(lop,ntyp).eq.l) then
                   do m = -l,l
                     do i = 1,noccbd
                      do j = 1,noccbd
                       qlo(i,j,lop,lo,ntyp) = qlo(i,j,lop,lo,ntyp) + 
     +                        real(conjg(ccof(m,j,lop,natom,2))
     *                                  *ccof(m,i,lo,natom,1))
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
            stop 'not yet finished'
            l = llo(lo,ntyp)
            do i = 1,noccbd
             do j = 1,noccbd
               mmn(i,j)= mmn(i,j)  + 
     +                      ( qaclo(i,j,lo,ntyp)*uulon(lo,ntyp,2)     +
     +                        qbclo(i,j,lo,ntyp)*dulon(lo,ntyp,2)     )
             enddo
            enddo 
            do lop = 1,nlo(ntyp)
               if (llo(lop,ntyp).eq.l) then
               do i = 1,noccbd
                do j = 1,noccbd
                 mmn(i,j) = mmn(i,j)  + 
     +                      qlo(i,j,lop,lo,ntyp)*uloulopn(lop,lo,ntyp,2)
                enddo
               enddo
               endif
            enddo
         enddo 
      enddo 
      deallocate ( qlo,qaclo,qbclo )

      END SUBROUTINE wann_mmk0_updown_sph
      END MODULE m_wann_mmk0_updown_sph
