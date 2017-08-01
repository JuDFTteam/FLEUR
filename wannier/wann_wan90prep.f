!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_wan90prep
      use m_juDFT
c**********************************************************
c         prepares WF//spin12(jspin)//.win for 
c         input to wannier90 and creates bkpts file
c                FF, September 2006
c**********************************************************
      CONTAINS
      SUBROUTINE wann_wan90prep(
     >      jspins,amat,bmat,natd,taual_inp,zatom,ntype,
     >      ntypd,neq,l_bzsym,film,l_onedimens,l_ms,l_sgwf,l_socgwf,
     >      aux_latt_const,param_file,l_dim)

      USE m_wann_postproc_setup
      USE m_wann_postproc_setup4
      USE m_wann_postproc_setup5
      use m_wann_get_mp
      use m_wann_get_kpts
      use m_wann_get_qpts
      use m_wann_gwf_tools, only : get_index_kq,get_dimension,get_shift
      use m_wann_gwf_auxbrav

      IMPLICIT NONE

      integer,intent(in)  :: jspins
      integer,intent(in)  :: ntype
      integer,intent(in)  :: ntypd
      logical,intent(in)  :: l_bzsym
      logical,intent(in)  :: film
      real,intent(in)     :: amat(3,3),aux_latt_const
      real,intent(in)     :: bmat(3,3)
      integer,intent(in)  :: natd
      real,intent(in)     :: taual_inp(3,natd)
      real,intent(in)     :: zatom(ntype)
      integer,intent(in)  :: neq(ntypd)
      logical,intent(in)  :: l_onedimens

      character(len=20),intent(in) :: param_file

      real             :: scale,scale_q
      real             :: taual(3,natd)
      real             :: amat_ang(3,3)
      real             :: bmat_ang(3,3)
      integer          :: at,j,n,bb
      integer          :: nkpts,iter,len,num_wann,num_bands,nn,i
      integer          :: num(3),dim,jspin,num_iter,nntot
      real             :: weight,maxi,mini,increm,compare
      real,allocatable :: kpoints(:,:)
      logical          :: l_file
      character(len=1) :: spin12(2)
      data spin12/'1','2'/
      integer          :: nkptd,ikpt_help,ikpt,iqpt
      character*2      :: namat(0:103)
      real,parameter   :: bohr=0.5291772108
      character(len=2) :: spin012(0:2)
      data spin012/'  ', '.1', '.2'/
      character(len=6) :: filename


      real,allocatable :: qpoints(:,:)
      real             :: amat_ang_q(3,3),bmat_ang_q(3,3)
      real             :: amat_q(3,3),bmat_q(3,3)
      real             :: taual_q(3,natd)
      integer          :: nqpts=1,nqptd=1
      integer          :: numq(3)
      logical, intent(in) :: l_ms,l_sgwf,l_socgwf
      logical, intent(in) :: l_dim(3)
      character(len=20) :: win_filename

      real,allocatable :: kqpoints(:,:)
      integer          :: nkqpts,ikqpt
      integer,allocatable :: numkq(:)
      real,allocatable :: amat_kq(:,:),bmat_kq(:,:)
      real,allocatable :: taual_kq(:,:)
      real             :: tpi,dummy_omtil
      integer :: arr_len,shift(3)

      DATA namat/'va',' h','he','li','be',' b',' c',' n',' o',' f','ne',
     +     'na','mg','al','si',' p',' s','cl','ar',' k','ca','sc','ti',
     +     ' v','cr','mn','fe','co','ni','cu','zn','ga','ge','as','se',
     +     'br','kr','rb','sr',' y','zr','nb','mo','tc','ru','rh','pd',
     +     'ag','cd','in','sn','sb','te',' j','xe','cs','ba','la','ce',
     +     'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     +     'lu','hf','ta',' w','re','os','ir','pt','au','hg','tl','pb',
     +     'bi','po','at','rn','fr','ra','ac','th','pa',' u','np','pu',
     +     'am','cm','bk','cf','es','fm','md','no','lw'/
      num_iter=100

      call get_dimension(l_dim,arr_len)
      call get_shift(l_dim,shift)

c**********************************************************
c     Read in kpoints from kpts/w90kpts file.
c**********************************************************
      call wann_get_kpts(
     >         l_bzsym,film,l_onedimens,.false.,
     <         nkpts,kpoints)
      allocate(kpoints(3,nkpts))
      call wann_get_kpts(
     >         l_bzsym,film,l_onedimens,.true.,
     <         nkpts,kpoints)


c**********************************************************
c     Read in qpoints from qpts file.
c**********************************************************
      IF(l_sgwf.or.l_socgwf)THEN
      call wann_get_qpts(
!     >         l_bzsym,film,l_onedimens,.false.,
     >         .false.,.false.,.false.,.false.,
     <         nqpts,qpoints,param_file)
      allocate(qpoints(3,nqpts))
      call wann_get_qpts(
!     >         l_bzsym,film,l_onedimens,.true.,
     >         .false.,.false.,.false.,.true.,
     <         nqpts,qpoints,param_file)
      ELSE
         allocate(qpoints(3,1))
         qpoints=0.
      ENDIF
c*********************************************************
c     Find out the structure of k-point set.
c*********************************************************
      call wann_get_mp(
     >         nkpts,kpoints,
     <         num)


      taual=taual_inp
      if(film.and.l_onedimens)then
         do n=1,natd
           taual(1,n)=taual_inp(1,n)+0.5
           taual(2,n)=taual_inp(2,n)+0.5
         enddo
      elseif(film)then
         do n=1,natd
           taual(3,n)=taual_inp(3,n)+0.5
         enddo
      endif

c*********************************************************
c     Find out the structure of q-point set.              
c*********************************************************
      IF(l_sgwf.or.l_socgwf)THEN
      call wann_get_mp(
     >         nqpts,qpoints,
     <         numq)
      ELSE
         numq=1
      ENDIF

c*******************************************************
c     Write information to WF//spin12(jspin)//.win
c*******************************************************
      if(l_ms) nqptd = nqpts
         
      do jspin=1,jspins
c        proj file provides num_wann and num_bands
         l_file=.false.
         do j=jspin,0,-1
           inquire(file=trim('proj'//spin012(j)),exist=l_file)
           if(l_file)then
             filename='proj'//spin012(j)
             exit
           endif
         enddo
         if(l_file)then
            open(712,file=trim(filename),form='formatted',status='old')
            rewind(712)
         else
            CALL juDFT_error("no proj/proj.1/proj.2",calledby
     +           ="wann_wan90prep")
         endif
         read(712,*)num_wann,num_bands
         close(712)

         do iqpt=1,nqptd
            win_filename = 'WF'//spin12(jspin)      
            if(l_ms) then
              write(win_filename,'("WF",a1,"_",i4.4)')spin12(jspin),iqpt
            endif

            call wann_write_win(win_filename,film,num_wann,num_bands,
     >               num_iter,jspin,ntype,ntypd,natd,nkpts,neq,num,amat,
     >               kpoints,zatom,taual,namat,3)
         enddo
            
      enddo!jspin

      amat_ang=amat*bohr
      bmat_ang=bmat/bohr

c******************************************************
c     call wannier90 routines to get bkpts
c******************************************************
      call wann_postproc_setup(
     >         natd,nkpts,kpoints,amat_ang,bmat_ang,
     >         num,num_bands,ntype,neq,
     >         zatom,taual,namat,win_filename,'bkpts')

c******************************************************
c     call wannier90 routines to get bqpts
c******************************************************
      if(.not.(l_sgwf.or.l_socgwf)) goto 8765

      call wann_gwf_auxbrav(aux_latt_const,l_sgwf,l_socgwf,
     >                      amat_q,bmat_q,l_dim)
      amat_ang_q = amat_q*bohr
      bmat_ang_q = bmat_q/bohr

      if(l_sgwf.or.l_socgwf) then
         qpoints = qpoints/2.0
      endif

      taual_q=0.5
      if(l_dim(1))taual_q(1,:)=0.0
      if(l_dim(2))taual_q(2,:)=0.0
      if(l_dim(3))taual_q(3,:)=0.0

      call wann_write_win('WF1_q',film,num_wann,num_bands,
     >      num_iter,1,ntype,ntypd,natd,nqpts,neq,numq,amat_q,
     >      qpoints,zatom,taual_q,namat,3)

      call wann_postproc_setup(
     >         natd,nqpts,qpoints,amat_ang_q,bmat_ang_q,
     >         numq,num_bands,ntype,neq,
     >         zatom,taual_q,namat,'WF1_q','bqpts') 

      !bkqpts part below
      nkqpts = nqpts*nkpts
      allocate(kqpoints(arr_len,nkqpts))
      kqpoints = 0.
      do iqpt=1,nqpts
         do ikpt=1,nkpts 
            ikqpt = get_index_kq(ikpt,iqpt,nkpts)
            kqpoints(1,ikqpt) = kpoints(1,ikpt)
            kqpoints(2,ikqpt) = kpoints(2,ikpt)
            kqpoints(3,ikqpt) = kpoints(3,ikpt)
            if(l_dim(1)) kqpoints(4+shift(1),ikqpt)=qpoints(1,iqpt)
            if(l_dim(2)) kqpoints(4+shift(2),ikqpt)=qpoints(2,iqpt)
            if(l_dim(3)) kqpoints(4+shift(3),ikqpt)=qpoints(3,iqpt)
         enddo
      enddo

      allocate(numkq(arr_len))
      allocate(taual_kq(arr_len,natd))
      allocate(amat_kq(arr_len,arr_len))
      allocate(bmat_kq(arr_len,arr_len))

      numkq(1)=num(1); numkq(2)=num(2); numkq(3)=num(3)
      if(l_dim(1))numkq(4+shift(1))=numq(1)
      if(l_dim(2))numkq(4+shift(2))=numq(2)
      if(l_dim(3))numkq(4+shift(3))=numq(3)

      taual_kq(1,:) = taual(1,:); taual_kq(2,:) = taual(2,:)
      taual_kq(3,:) = taual(3,:); taual_kq(4:arr_len,:) = 0.0!0.0

      ! set up amat_kq
      amat_kq = 0.0
      bmat_kq = 0.0
      amat_kq(1:3,1:3) = amat(1:3,1:3)
      bmat_kq(1:3,1:3) = bmat(1:3,1:3)
      if(l_dim(1)) amat_kq(4+shift(1),4+shift(1)) = amat_q(1,1)
      if(l_dim(1)) bmat_kq(4+shift(1),4+shift(1)) = bmat_q(1,1)
      if(l_dim(2)) amat_kq(4+shift(2),4+shift(2)) = amat_q(2,2)
      if(l_dim(2)) bmat_kq(4+shift(2),4+shift(2)) = bmat_q(2,2)
      if(l_dim(3)) amat_kq(4+shift(3),4+shift(3)) = amat_q(3,3)
      if(l_dim(3)) bmat_kq(4+shift(3),4+shift(3)) = bmat_q(3,3)

      call wann_write_win('WF1_gwf',film,num_wann,num_bands,
     >      0,1,ntype,ntypd,natd,nkqpts,neq,numkq,amat_kq,
     >      kqpoints,zatom,taual_kq,namat,arr_len)


      if(arr_len.eq.4) then
         call wann_postproc_setup4(
     >         natd,nkqpts,kqpoints,amat_kq*bohr,bmat_kq/bohr,
     >         numkq,num_bands,ntype,neq,
     >         zatom,taual_kq,namat,'WF1_gwf','bkqpts') 
      elseif(arr_len.eq.5) then
         call wann_postproc_setup5(
     >         natd,nkqpts,kqpoints,amat_kq*bohr,bmat_kq/bohr,
     >         numkq,num_bands,ntype,neq,
     >         zatom,taual_kq,namat,'WF1_gwf','bkqpts') 
!         call juDFT_error("arr_len.eq.5 not yet implemented",
!     >                  calledby='wann_wan90prep')
      elseif(arr_len.eq.6) then
         call juDFT_error("arr_len.eq.6 not yet implemented",
     >                  calledby='wann_wan90prep')
      else
         call juDFT_error("Dimension arr_len not recognized",
     >                  calledby='wann_wan90prep')
      endif

      deallocate(kqpoints)
      deallocate(amat_kq,bmat_kq,numkq,taual_kq)

8765  continue

      deallocate(kpoints)
      deallocate(qpoints)

      END SUBROUTINE wann_wan90prep


      subroutine wann_write_win(win_filename,film,num_wann,num_bands,
     >               num_iter,jspin,ntype,ntypd,natd,nkpts,neq,num,amat,
     >               kpoints,zatom,taual,namat,rdim)
      implicit none
      character(len=*),intent(in) ::win_filename
      logical,intent(in) :: film
      integer,intent(in) :: num_wann,num_iter,num_bands
      integer,intent(in) :: rdim
      integer,intent(in) :: jspin,ntype,natd,nkpts,ntypd
      integer,intent(in) :: neq(ntypd),num(rdim)
      real,intent(in) :: amat(rdim,rdim),kpoints(rdim,nkpts)
      real,intent(in) :: zatom(ntype),taual(rdim,natd)
      character*2,intent(in) :: namat(0:103)

      integer :: dim,nn,iter,at,i,search_shells
      logical :: l_exist

         open(911,file=trim(win_filename)//'.win')
         write(911,*)"length_unit=Bohr"
         write(911,*)"num_wann=",num_wann
         write(911,*)"num_iter=",num_iter
         write(911,*)"num_bands=",num_bands
         write(911,*)"          "

	 if(rdim.gt.3) then
            search_shells=200
            inquire(file='searchshells_inp',exist=l_exist)
            if(l_exist) then
             open(777,file='searchshells_inp')
             read(777,*)search_shells
             close(777)
             write(*,*)'search_shells=',search_shells
            endif
            inquire(file=trim(win_filename)//'.kshell',exist=l_exist)
            if(l_exist) then
             write(*,*)'found .kshell file; set devel_flag'
             write(911,*)'devel_flag=kmesh_degen'
            endif
            write(911,*)"!iprint=",5
            write(911,*)"search_shells=",search_shells
         else
            write(911,*)"search_shells=",200
         endif

         write(911,*)"!optional parameters for wannierization"
         write(911,*)"!num_cg_steps="
         write(911,*)"!trial_step="
         write(911,*)"!fixed_step="
         write(911,*)"!restart=wannierise"
         write(911,*)"         "

         if(num_bands.ne.num_wann)then
            write(911,*)"! optional parameters for disentangling"
            write(911,*)"!dis_win_min="
            write(911,*)"!dis_win_max="
            write(911,*)"!dis_froz_min="
            write(911,*)"!dis_froz_max="
            write(911,*)"dis_num_iter=10000"
            write(911,*)"!dis_mix_ratio="
            write(911,*)"!dis_conv_tol="
            write(911,*)"!dis_conv_window="
            write(911,*)"            "
         endif

         write(911,*)"! optional parameters for plotting"
         if(jspin.eq.1)then
            write(911,*)"spin=up"
         else   
            write(911,*)"spin=down"
         endif   
         write(911,*)"!restart=plot"
         write(911,*)"!wannier_plot=true"
         write(911,*)"!wannier_plot_supercell=3"
         write(911,*)"!bands_plot=true"


         write(911,*)"!fermi_surface_plot=true"
         write(911,*)"            "

         write(911,*)"!options for Hamiltonian in Wannier basis"
         write(911,*)"!HR_PLOT=true"
         write(911,*)"!DIST_CUTOFF=3.0"
         write(911,*)"           "

         write(911,*)"!some more options"
         write(911,*)"!WRITE_R2MN=true"
         write(911,*)"!NUM_PRINT_CYCLES=10"

         write(911,*)"begin unit_cell_cart"
         write(911,*)"bohr"
         if(rdim.eq.3) then
            do dim=1,3
               write(911,*)amat(:,dim)
            enddo   
         elseif(rdim.eq.4)then
            do dim=1,4
               write(911,'(4f14.8)')amat(:,dim)
            enddo             
         elseif(rdim.eq.5)then
            do dim=1,5
               write(911,'(5f14.8)')amat(:,dim)
            enddo             
         elseif(rdim.eq.6)then
            do dim=1,6
               write(911,'(6f13.7)')amat(:,dim)
            enddo             
         endif
         write(911,*)"end unit_cell_cart"
         write(911,*)

         write(911,*)"begin atoms_frac"
         if(film)then
           write(911,*)!for reasons of plotting: shift the
           write(911,*)!coordinates in the 1d and 2d case
         endif
         nn=0
         do iter=1,ntype
            at=nint(zatom(iter))
            do i=1,neq(iter)
               nn=nn+1
             if(rdim.eq.3) write(911,*)namat(at),taual(:,nn)   
             if(rdim.eq.4) write(911,'(1x,a2,2x,4f12.6)')
     >                     namat(at),taual(:,nn)
             if(rdim.eq.5) write(911,'(1x,a2,2x,5f12.6)')
     >                     namat(at),taual(:,nn)
             if(rdim.eq.6) write(911,'(1x,a2,2x,6f12.6)')
     >                     namat(at),taual(:,nn)
            enddo   
         enddo   
         write(911,*)"end atoms_frac"
         write(911,*)

         if(rdim.eq.3) write(911,*)"mp_grid",(num(dim),dim=1,3)
         if(rdim.eq.4) write(911,'(1x,a7,2x,4(i4,1x))')
     >                 "mp_grid",(num(dim),dim=1,4)
         if(rdim.eq.5) write(911,'(1x,a7,2x,5(i4,1x))')
     >                 "mp_grid",(num(dim),dim=1,5)
         if(rdim.eq.6) write(911,'(1x,a7,2x,6(i4,1x))')
     >                 "mp_grid",(num(dim),dim=1,6)
         write(911,*)

         write(911,*)"begin kpoints"
         if(rdim.eq.3) then
            do iter=1,nkpts
               write(911,*)kpoints(:,iter)
            enddo
         elseif(rdim.eq.4)then
            do iter=1,nkpts
               write(911,'(4f19.15)')kpoints(:,iter)
            enddo
         elseif(rdim.eq.5)then
            do iter=1,nkpts
               write(911,'(5f15.11)')kpoints(:,iter)
            enddo
         elseif(rdim.eq.6)then
            do iter=1,nkpts
               write(911,'(6f13.9)')kpoints(:,iter)
            enddo
         endif
         write(911,*)"end kpoints"
         write(911,*)

         write(911,*)"!begin kpoint_path"
         if(rdim.eq.3) then
            write(911,*)"!X -0.5   0.0   0.0   G  0.0   0.0   0.0"
            write(911,*)"!G  0.0   0.0   0.0   X  0.5   0.0   0.0"
         elseif(rdim.eq.4) then
            write(911,*)"!X  0.0   0.0  -0.5   0.0    ",
     >                   "G  0.0   0.0   0.0   0.0"
            write(911,*)"!G  0.0   0.0   0.0   0.0    ",
     >                   "X  0.0   0.0   0.5   0.0"
         elseif(rdim.eq.5) then
            write(911,*)"!X  0.0  0.0 -0.5  0.0  0.0   ",
     >                   "G  0.0  0.0  0.0  0.0  0.0"
            write(911,*)"!G  0.0  0.0  0.0  0.0  0.0   ",
     >                   "X  0.0  0.0  0.5  0.0  0.0"
         elseif(rdim.eq.6) then
            write(911,*)"!X  0.0  0.0 -0.5  0.0  0.0  0.0   ",
     >                   "G  0.0  0.0  0.0  0.0  0.0  0.0"
            write(911,*)"!G  0.0  0.0  0.0  0.0  0.0  0.0   ",
     >                   "X  0.0  0.0  0.5  0.0  0.0  0.0"
         endif
         write(911,*)"!end kpoint_path"
         write(911,*)

         write(911,*)"wvfn_formatted=.true."
        
         close(911)


      end subroutine wann_write_win


      END MODULE m_wann_wan90prep
