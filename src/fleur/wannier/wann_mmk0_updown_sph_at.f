      MODULE m_wann_mmk0_updown_sph_at
        use m_juDFT
      CONTAINS
      SUBROUTINE wann_mmk0_updown_sph_at(
     >               l_noco,alph,beta,
     >               llod,noccbd,nlod,natd,ntypd,lmaxd,lmax,lmd,
     >               ntype,neq,nlo,llo,
     >               radial1_ff,radial1_gg,
     >               radial1_fg,radial1_gf,
     >               acof,bcof,ccof,
     >               ddn,uulon,dulon,uloulopn,
     >               atomlist_num,atomlist,
     =               mmn)
c**************************************************************
c     Overlaps of the spin-down parts of the Bloch functions
c     with the spin-up parts in the MT-spheres. Atom-resolved.
c                           Frank Freimuth
c**************************************************************
      implicit none
      logical, intent (in)  :: l_noco
      integer, intent (in)  :: llod,nlod,natd,ntypd,lmaxd,lmd
      integer, intent (in)  :: lmax(:) !(ntypd)
      integer, intent (in)  :: ntype,noccbd
      REAL,    INTENT (IN)  :: alph(ntypd),beta(ntypd)
      integer, intent (in)  :: neq(ntypd)
      integer, intent (in)  :: nlo(ntypd),llo(nlod,ntypd)
      real,    intent (in)  :: radial1_ff(:,:,0:,0:,:)
      real,    intent (in)  :: radial1_gg(:,:,0:,0:,:)      
      real,    intent (in)  :: radial1_fg(:,:,0:,0:,:)
      real,    intent (in)  :: radial1_gf(:,:,0:,0:,:)
      real,    intent (in)  :: ddn(0:lmaxd,ntypd,2)
      real,    intent (in)  :: uloulopn(nlod,nlod,ntypd,2)
      real,    intent (in)  :: uulon(nlod,ntypd,2),dulon(nlod,ntypd,2)
      complex, intent (in)  :: ccof(-llod:llod,noccbd,nlod,natd,2)
      complex, intent (in)  :: acof(noccbd,0:lmd,natd,2)
      complex, intent (in)  :: bcof(noccbd,0:lmd,natd,2)
      integer, intent(in)   :: atomlist_num
      integer, intent(in)   :: atomlist(:)

      complex, intent (inout) :: mmn(:,:,:) !mmn(noccbd,noccbd,natd)

      integer           :: i,j,l,lo,lop,m,natom,nn,ntyp
      integer           :: nt1,nt2,lm,n,ll1,i1spin,i2spin
      complex           :: suma(natd),sumb(natd)
      complex           :: sumc(natd),sumd(natd)
      complex           :: suma12(2,2),sumb12(2,2)
      complex           :: sumc12(2,2),sumd12(2,2)
      real, allocatable :: qlo(:,:,:,:,:)
      real, allocatable :: qaclo(:,:,:,:),qbclo(:,:,:,:)
      COMPLEX           :: ccchi(2,2),ci
      integer           :: nat2
      logical           :: l_inthelist

      call timestart("wann_mmk0_updown_sph_at")
      ci = cmplx(0.0,1.0)
      allocate (qlo(noccbd,noccbd,nlod,nlod,natd), 
     +          qaclo(noccbd,noccbd,nlod,natd),
     +          qbclo(noccbd,noccbd,nlod,natd) )
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
            do l = 0,lmax(n)
             if(.not.l_noco)then  
               suma = cmplx(0.,0.)
               sumb = cmplx(0.,0.)
               sumc = cmplx(0.,0.)
               sumd = cmplx(0.,0.)
               ll1 = l* (l+1)
               do m = -l,l
                  lm = ll1 + m
                  do natom = nt1,nt2
                    suma(natom) = suma(natom) + acof(i,lm,natom,1)*
     +                      conjg(acof(j,lm,natom,2))
                    sumb(natom) = sumb(natom) + bcof(i,lm,natom,1)*
     +                      conjg(bcof(j,lm,natom,2))
                    sumc(natom) = sumc(natom) + acof(i,lm,natom,1)*
     +                      conjg(bcof(j,lm,natom,2))
                    sumd(natom) = sumd(natom) + bcof(i,lm,natom,1)*
     +                      conjg(acof(j,lm,natom,2))
                  enddo !natom
               enddo !m      
               do natom=nt1,nt2
                 l_inthelist=.false. 
                 do nat2=1,atomlist_num
                   if(atomlist(nat2).eq.natom)then
                      l_inthelist=.true.
                      exit
                   endif    
                 enddo !nat2
                 if(l_inthelist)then
                   mmn(j,i,nat2) = mmn(j,i,nat2) + 
     +                     ( suma(natom)*radial1_ff(1,2,l,l,n)+
     +                       sumb(natom)*radial1_gg(1,2,l,l,n)+
     +                       sumc(natom)*radial1_fg(1,2,l,l,n)+
     +                       sumd(natom)*radial1_gf(1,2,l,l,n)  )      
                 endif
               enddo !natom
             else
               stop 'not yet finished' 
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
                 do natom=nt1,nt2
                  mmn(i,j,natom) =    mmn(i,j,natom)
     &                           +
     &         suma12(i1spin,i2spin)*radial1_ff(i1spin,i2spin,l,l,n)
     &               *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &         sumb12(i1spin,i2spin)*radial1_gg(i1spin,i2spin,l,l,n)
     &               *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &         sumc12(i1spin,i2spin)*radial1_fg(i1spin,i2spin,l,l,n)
     &               *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))

     &                           +
     &         sumd12(i1spin,i2spin)*radial1_gf(i1spin,i2spin,l,l,n)
     &               *ccchi(1,i2spin)*conjg(ccchi(2,i1spin))
                 enddo !natom 
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
              if(l_noco)then
                 stop 'not yet finished'
              else
               l = llo(lo,ntyp)
               ll1 = l* (l+1)
               do m = -l,l
                  lm = ll1 + m
                  do i = 1,noccbd
                   do j = 1,noccbd
                 qbclo(j,i,lo,natom) = qbclo(j,i,lo,natom) + real(
     +               bcof(i,lm,natom,1)*conjg(ccof(m,j,lo,natom,2)) +
     +               ccof(m,i,lo,natom,1)*conjg(bcof(j,lm,natom,2)) )
                 qaclo(j,i,lo,natom) = qaclo(j,i,lo,natom) + real(
     +               acof(i,lm,natom,1)*conjg(ccof(m,j,lo,natom,2)) +
     +               ccof(m,i,lo,natom,1)*conjg(acof(j,lm,natom,2)) )
                   enddo
                  enddo
               enddo
               do lop = 1,nlo(ntyp)
                 if (llo(lop,ntyp).eq.l) then
                   do m = -l,l
                     do i = 1,noccbd
                      do j = 1,noccbd
                   qlo(j,i,lop,lo,natom) = qlo(j,i,lop,lo,natom)+ 
     +                     real(conjg(ccof(m,j,lop,natom,2))
     *                               *ccof(m,i,lo,natom,1))
                      enddo
                     enddo
                   enddo
                 endif
               enddo !lop
              endif !l_noco 
            enddo !lo
         enddo !nn
      enddo !ntyp
c---> perform summation of the coefficients with the integrals
c---> of the radial basis functions
      natom=0
      do ntyp = 1,ntype
        do nn=1,neq(ntyp)
          natom=natom+1

          l_inthelist=.false.
          do nat2=1,atomlist_num
              if(atomlist(nat2).eq.natom)then
                 l_inthelist=.true.
                 exit
              endif
          enddo !nat2   
          if(.not.l_inthelist) cycle

          do lo = 1,nlo(ntyp)
            l = llo(lo,ntyp)
            do i = 1,noccbd
             do j = 1,noccbd
               mmn(j,i,nat2)= mmn(j,i,nat2)  + 
     +                      ( qaclo(j,i,lo,natom)*uulon(lo,ntyp,2) +
     +                        qbclo(j,i,lo,natom)*dulon(lo,ntyp,2)  )
             enddo
            enddo 
            do lop = 1,nlo(ntyp)
               if (llo(lop,ntyp).eq.l) then
               do i = 1,noccbd
                do j = 1,noccbd
                 mmn(j,i,nat2) = mmn(j,i,nat2)  + 
     +                  qlo(j,i,lop,lo,natom)*uloulopn(lop,lo,ntyp,2)
                enddo
               enddo
               endif
            enddo
          enddo !lo 
        enddo !nn  
      enddo !ntyp
      deallocate ( qlo,qaclo,qbclo )

      call timestop("wann_mmk0_updown_sph_at")
      END SUBROUTINE wann_mmk0_updown_sph_at
      END MODULE m_wann_mmk0_updown_sph_at