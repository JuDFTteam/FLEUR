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
     >               jspins,amat,bmat,natd,taual_inp,zatom,ntype,
     >               ntypd,neq,l_bzsym,film,l_onedimens)

      USE m_wann_postproc_setup
      use m_wann_get_mp
      use m_wann_get_kpts
      IMPLICIT NONE

      integer,intent(in)  :: jspins
      integer,intent(in)  :: ntype
      integer,intent(in)  :: ntypd
      logical,intent(in)  :: l_bzsym
      logical,intent(in)  :: film
      real,intent(in)     :: amat(3,3)
      real,intent(in)     :: bmat(3,3)
      integer,intent(in)  :: natd
      real,intent(in)     :: taual_inp(3,natd)
      real,intent(in)     :: zatom(ntype)
      integer,intent(in)  :: neq(ntypd)
      logical,intent(in)  :: l_onedimens

      real             :: scale
      real             :: taual(3,natd)
      real             :: amat_ang(3,3)
      real             :: bmat_ang(3,3)
      integer          :: at,j,n
      integer          :: nkpts,iter,len,num_wann,num_bands,nn,i
      integer          :: num(3),dim,jspin,num_iter,nntot
      real             :: weight,maxi,mini,increm,compare
      real,allocatable :: kpoints(:,:)
      logical          :: l_file
      character(len=1) :: spin12(2)
      data spin12/'1','2'/
      integer          :: nkptd,ikpt_help,ikpt
      character*2      :: namat(0:103)
      real,parameter   :: bohr=0.5291772108
      character(len=2) :: spin012(0:2)
      data spin012/'  ', '.1', '.2'/
      character(len=6) :: filename

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

c*******************************************************
c     Write information to WF//spin12(jspin)//.win
c*******************************************************
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
         print*,"num_wann=",num_wann
         print*,"num_bands=",num_bands

         open(911,file='WF'//spin12(jspin)//'.win')
         write(911,*)"length_unit=Bohr"
         write(911,*)"num_wann=",num_wann
         write(911,*)"num_iter=",num_iter
         write(911,*)"num_bands=",num_bands
         write(911,*)"          "

	 write(911,*)"search_shells=",200

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
         do dim=1,3
            write(911,*)amat(:,dim)
         enddo   
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
               write(911,*)namat(at),taual(:,nn)   
            enddo   
         enddo   
         write(911,*)"end atoms_frac"
         write(911,*)

         write(911,*)"mp_grid",(num(dim),dim=1,3)
         write(911,*)

         write(911,*)"begin kpoints"
         do iter=1,nkpts
            write(911,*)kpoints(:,iter)
         enddo
         write(911,*)"end kpoints"
         write(911,*)

         write(911,*)"wvfn_formatted=.true."
         close(911)
      enddo

      amat_ang=amat*bohr
      bmat_ang=bmat/bohr

c******************************************************
c     call wannier90 routines to get bkpts
c******************************************************
      call wann_postproc_setup(
     >         natd,nkpts,kpoints,amat_ang,bmat_ang,
     >         num,num_bands,ntype,neq,
     >         zatom,taual,namat)

      deallocate(kpoints)

      END SUBROUTINE wann_wan90prep
      END MODULE m_wann_wan90prep
