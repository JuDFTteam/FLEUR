      module m_wann_optional
      use m_juDFT
      contains
      subroutine wann_optional(
     >               jspins,ntype,neq,natd,zatom,
     >               nop,mrot,bmat,amat,tau,
     >               taual,film,
     >               l_onedimens,l_soc,l_noco,
     >               omtil,pos)
c**************************************************
c     Make preparations for the calculation of 
c     Wannier functions.
c     Frank Freimuth
c**************************************************
      use m_types
      use m_wann_read_inp
      use m_wann_projgen
      use m_wann_kpointgen
      use m_wann_w90kpointgen
      use m_wann_kptsreduc
      use m_wann_kptsreduc2
      use m_wann_wan90prep
      use m_wann_dipole3
      use m_wann_dipole

      implicit none
      integer, intent(in) :: jspins
      integer, intent(in) :: ntype
      integer, intent(in) :: neq(ntype)
      integer, intent(in) :: natd
      real,intent(in)     :: zatom(ntype)
      integer,intent(in)  :: nop
      integer,intent(in)  :: mrot(3,3,nop)
      real,intent(in)     :: bmat(3,3),amat(3,3)
      real,intent(in)     :: tau(3,nop)
      real,intent(in)     :: taual(3,natd)
      logical,intent(in)  :: film
      logical,intent(in)  :: l_onedimens
      logical,intent(in)  :: l_soc
      logical,intent(in)  :: l_noco
      real,intent(in)     :: omtil
      real,intent(in)     :: pos(3,natd)

      type(t_wann)  :: wann
      integer       :: num_wann(2)
      logical       :: l_nocosoc

      l_nocosoc=l_noco.or.l_soc

c-----read the input file to determine what to do
      call wann_read_inp(
     >         .true.,wann)

c-----generate projection-definition-file
      if(wann%l_projgen) then
         call wann_projgen(
     >       ntype,neq,natd,zatom,l_nocosoc, wann)
         wann%l_stopopt=.true.
      endif

c-----generate k-point-files
      if(wann%l_kpointgen) then
             call wann_kpointgen()
             wann%l_stopopt=.true.
      endif
      if(wann%l_w90kpointgen) then
             call wann_w90kpointgen()
             wann%l_stopopt=.true.
      endif

c-----find Wannier-irreducible part of BZ
      if(wann%l_kptsreduc)then
         call wann_kptsreduc(
     >            nop,mrot,bmat,tau,film,
     >            l_onedimens,(l_soc.or.l_noco))
         wann%l_stopopt=.true.
      endif

c-----find Wannier-irreducible part of BZ
      if(wann%l_kptsreduc2)then
         call wann_kptsreduc2(
     >            wann%mhp,
     >            nop,mrot,bmat,tau,film,
     >            l_onedimens,(l_soc.or.l_noco))
         wann%l_stopopt=.true.
      endif

c-----generate WF1.win and bkpts
      if(wann%l_prepwan90)then
         call wann_wan90prep(
     >            jspins,amat,bmat,
     >            natd,taual,zatom,ntype,
     >            ntype,neq,wann%l_bzsym,film,l_onedimens)
      endif

c-----calculate polarization, if not wannierize
c-----if wannierize, then calculate polarization later (after wannierize)
      if(wann%l_dipole3.and..not.wann%l_wannierize)then
         num_wann(1)=wann%band_max(1)-wann%band_min(1)+1
         num_wann(2)=wann%band_max(2)-wann%band_min(2)+1
         call wann_dipole3(
     >            jspins,omtil,natd,pos,
     >            amat,bmat,taual,num_wann,
     >            ntype,neq,zatom,(l_soc.or.l_noco))
         wann%l_stopopt=.true.
      endif

c-----calculate polarization, if not wannierize
c-----if wannierize, then calculate polarization later (after wannierize)
      if(wann%l_dipole.and..not.wann%l_wannierize)then
         call wann_dipole(
     >            jspins,omtil,natd,pos,
     >            amat,ntype,neq,zatom)
         wann%l_stopopt=.true.
      endif

      IF(wann%l_stopopt)  CALL juDFT_end("wann_optional done")

      end subroutine wann_optional
      end module m_wann_optional
