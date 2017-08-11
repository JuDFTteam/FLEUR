!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_optional
      use m_juDFT
      contains
      subroutine wann_optional(
     >               input,atoms,sym,cell,oneD,noco,wann,
     >               l_ms,l_sgwf,l_socgwf,
     >               aux_latt_const,param_file,l_dim)
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

      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_wann),      INTENT(INOUT) :: wann

      real,intent(in)     :: aux_latt_const
      character(len=20),intent(in) :: param_file

      integer       :: num_wann(2)
      logical       :: l_nocosoc

      logical, intent(in) :: l_ms,l_sgwf,l_socgwf
      logical,intent(in) :: l_dim(3)

      l_nocosoc=noco%l_noco.or.noco%l_soc

c-----read the input file to determine what to do
      call wann_read_inp(
     >         .true.,wann)

c-----generate projection-definition-file
      if(wann%l_projgen) then
         call wann_projgen(
     >      atoms%ntype,atoms%neq,atoms%nat,atoms%zatom,l_nocosoc,wann)
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
     >            sym%nop,sym%mrot,cell%bmat,sym%tau,input%film,
     >            oneD%odi%d1,l_nocosoc)
         wann%l_stopopt=.true.
      endif

c-----find Wannier-irreducible part of BZ
      if(wann%l_kptsreduc2)then
         call wann_kptsreduc2(
     >            wann%mhp,
     >            sym%nop,sym%mrot,cell%bmat,sym%tau,input%film,
     >            oneD%odi%d1,l_nocosoc)
         wann%l_stopopt=.true.
      endif

c-----generate WF1.win and bkpts
      if(wann%l_prepwan90)then
         call wann_wan90prep(
     >            input%jspins,cell%amat,cell%bmat,
     >            atoms%nat,atoms%taual,atoms%zatom,atoms%ntype,
     >            atoms%ntype,atoms%neq,wann%l_bzsym,input%film,
     >            oneD%odi%d1,l_ms,l_sgwf,l_socgwf,aux_latt_const,
     >            param_file,l_dim)
      endif

c-----calculate polarization, if not wannierize
c-----if wannierize, then calculate polarization later (after wannierize)
      if(wann%l_dipole3.and..not.wann%l_wannierize)then
         num_wann(1)=wann%band_max(1)-wann%band_min(1)+1
         num_wann(2)=wann%band_max(2)-wann%band_min(2)+1
         call wann_dipole3(
     >            input%jspins,cell%omtil,atoms%nat,atoms%pos,
     >            cell%amat,cell%bmat,atoms%taual,num_wann,
     >            atoms%ntype,atoms%neq,atoms%zatom,l_nocosoc)
         wann%l_stopopt=.true.
      endif

c-----calculate polarization, if not wannierize
c-----if wannierize, then calculate polarization later (after wannierize)
      if(wann%l_dipole.and..not.wann%l_wannierize)then
         call wann_dipole(
     >            input%jspins,cell%omtil,atoms%nat,atoms%pos,
     >            cell%amat,atoms%ntype,atoms%neq,atoms%zatom)
         wann%l_stopopt=.true.
      endif

      IF(wann%l_stopopt)  CALL juDFT_end("wann_optional done",1) ! The 1 is temporarily. Should be mpi%irank.

      end subroutine wann_optional
      end module m_wann_optional
