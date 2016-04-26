      module m_wann_projgen
      contains
      subroutine wann_projgen(
     >               ntype,neq,natd,zatom,l_nocosoc,wann)
c*****************************************
c     Generate the proj-file.
c     Frank Freimuth
c*****************************************
      use m_types, only: t_wann
      implicit none
      integer, intent(in) :: ntype
      integer, intent(in) :: neq(ntype)
      integer, intent(in) :: natd
      real,intent(in)     :: zatom(ntype)
      logical,intent(in)  :: l_nocosoc
      type(t_wann), intent(in) :: wann

      integer             :: j,nat,n,at,k,num_proj
      integer             :: num_wann,kk
      character(len=2)    :: symbol(100*ntype)
      character(len=2)    :: symb(natd)
      integer             :: lwf(100*ntype)
      integer             :: mrwf(100*ntype)
      integer             :: rwf(100*ntype)
      real                :: alpha(100*ntype),beta(100*ntype)
      real                :: gamma(100*ntype),zona(100*ntype)
      real                :: regio(100*ntype)      
      character*2         :: namat(0:103)
      integer             :: projections(-5:5)
      integer             :: projspin

      DATA namat/'va',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     +     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     +     ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     +     'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     +     'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',
     +     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     +     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',
     +     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw'/

      projections(:)=0
      projections(-5)=6
      projections(-4)=5
      projections(-3)=4
      projections(-2)=3
      projections(-1)=2
      projections(0)=1
      projections(1)=3
      projections(2)=5
      projections(3)=7


      nat=0
      do n=1,ntype
         do at=1,neq(n)
            nat=nat+1
            symb(nat)=namat(nint(zatom(n)))
         enddo
      enddo

      open(200,file='projgen_inp',status='old')
      j=0
      do 
         j=j+1
         read(200,*,end=299,err=299)symbol(j),lwf(j),mrwf(j),rwf(j)
         num_proj=j
      enddo   
 299  continue
      print*,"found ",num_proj,"projections"
      write(6,*)"found ",num_proj,"projections"
      close(200)

      open(300,file='proj')
      alpha=0.0
      beta=0.0
      gamma=0.0
      zona=0.0
      regio=1.0
      num_wann=0
      do j=1,nat
         do k=1,num_proj
            print*,"symbols:",symbol(k)," ",symb(j)
            write(6,*)"symbols:",symbol(k)," ",symb(j)
            if(trim(adjustl(symbol(k))).eq.
     &             trim(adjustl(symb(j))))then
               if(mrwf(k).ne.0)then
                  num_wann=num_wann+1
               else
                  num_wann=num_wann+projections(lwf(k))
               endif
            endif
         enddo
      enddo
      if(l_nocosoc)num_wann=num_wann*2
      print*,num_wann," wannier functions have been defined"
      if(l_nocosoc)then
         write(300,*) num_wann, 
     &                MAX(wann%band_max(1)-wann%band_min(1)+1,num_wann),
     &                "  t  "
      else
         write(300,*) num_wann, 
     &                MAX(wann%band_max(1)-wann%band_min(1)+1,num_wann)
      endif
      if(l_nocosoc)then
       do projspin=1,-1,-2  
        do j=1,nat
         do k=1,num_proj
            if(
     &          trim(adjustl(symbol(k)))
     &                    .eq.
     &          trim(adjustl(symb(j)))
     &                    )then
               if(mrwf(k).ne.0)then
                 write (300,'(i3,1x,i2,1x,i2,1x,i2,1x,i2)') 
     &                j,lwf(k),mrwf(k),rwf(k),projspin
                 write (300,'(2x,4f10.6,1x,f4.2)') 
     &                alpha(k),beta(k),gamma(k),
     &                zona(k),regio(k)
               else
                  do kk=1,projections(lwf(k))
                     write (300,'(i3,1x,i2,1x,i2,1x,i2,1x,i2)') 
     &                j,lwf(k),kk,rwf(k),projspin
                     write (300,'(2x,4f10.6,1x,f4.2)') 
     &                alpha(k),beta(k),gamma(k),
     &                zona(k),regio(k)
                  enddo
               endif
            endif
         enddo
        enddo
       enddo 
      else
        do j=1,nat
         do k=1,num_proj
            if(
     &          trim(adjustl(symbol(k)))
     &                    .eq.
     &          trim(adjustl(symb(j)))
     &                    )then
               if(mrwf(k).ne.0)then
                 write (300,'(i3,1x,i2,1x,i2,1x,i2)') 
     &                j,lwf(k),mrwf(k),rwf(k)
                 write (300,'(2x,4f10.6,1x,f4.2)') 
     &                alpha(k),beta(k),gamma(k),
     &                zona(k),regio(k)
               else
                  do kk=1,projections(lwf(k))
                     write (300,'(i3,1x,i2,1x,i2,1x,i2)') 
     &                j,lwf(k),kk,rwf(k)
                     write (300,'(3x,4f10.6,1x,f4.2)') 
     &                alpha(k),beta(k),gamma(k),
     &                zona(k),regio(k)
                  enddo
               endif
            endif
         enddo
        enddo
      endif
      close(300)

      end subroutine wann_projgen
      end module m_wann_projgen
