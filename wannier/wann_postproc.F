!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_postproc
      contains 
      subroutine wann_postproc(
     >          DIMENSION,stars,vacuum,atoms,sphhar,input,sym,mpi,
     >          lapw,oneD,noco,cell,eig_id,
     >          mpi_comm,wann,l_p0,film,jspins,n2d,
     >          natd,pos,amat,bmat,ntype,neq,zatom,
     >          omtil,l_soc,l_noco,neigd,fullnkpts,
     >          odi,l_proj,
     >          l_dulo,l_ss,lmaxd,ntypd,nop,nvd,jspd,nbasfcn,
     >          llod,nlod,nwdd,nlo,llo,lapw_l,invtab,mrot,ngopr,
     >          lmax,invsat,invsatnr,nkpt,taual,rmt,bbmat,alph,
     >          beta,qss,sk2,phi2,ods,
     >          irank,isize,n3d,nmzxyd,nmzd,jmtd,nlhd,nq3,nvac,
     >          invs,invs2,nlh,jri,ntypsy,ntypsd,nkptd,
     >          dx,rmsh,
     >          e1s,e2s,ulo_der,ustep,k1d,k2d,k3d,nop2,
     >          ig,rgphs,slice,kk,nnne,irecl,nv2d,nmzxy,nmz,nq2,
     >          delz,area,z1,volint,ig2,tau,symor,ef,l_sgwf,fullnqpts)
c     <          fermi_weights)

c***********************************************
c     Collection of those subroutines which may
c     be called after the Wannier functions have
c     been computed.
c     Frank Freimuth
c***********************************************
      use m_types
      use m_wann_dipole
      use m_wann_dipole2
      use m_wann_wannierize
      use m_wann_hopping
      use m_wann_plot_um_dat
      use m_wannier_to_lapw
      use m_wann_plot_from_lapw
      use m_wann_nabla_rs
      use m_wann_pauli_rs
      use m_wann_nabla_pauli_rs
      use m_wann_socmat_rs
      use m_wann_perpmag_rs
      use m_types
      use m_wann_wigner_seitz
      use m_wann_get_mp
      use m_wann_get_kpts
      use m_wann_fft4
      use m_wann_fft5
      use m_wann_rmat
      use m_wann_nocoplot
#ifdef CPP_TOPO
      use m_wann_torque_rs
      use m_wann_offdiposop_rs
      use m_wann_fft6
      use m_wann_fft3
      use m_wann_nedrho
#endif

#ifdef DCPP_WANN_EXT
      use m_wannier_to_lapw_kpts
      use m_wannier_lapw_gfleur
#endif
      implicit none

      TYPE(t_dimension),INTENT(IN) :: DIMENSION
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_lapw),INTENT(IN)      :: lapw
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_cell),INTENT(IN)      :: cell

      type(t_wann), intent(in) :: wann
      logical,      intent(in) :: l_p0
      logical,      intent(in) :: film
      integer,      intent(in) :: jspins,n2d,mpi_comm, eig_id, irecl
      integer,      intent(in) :: natd
      real,         intent(in) :: pos(3,natd)
      real,         intent(in) :: amat(3,3),bmat(3,3)
      integer,      intent(in) :: ntype
      integer,      intent(in) :: neq(ntype)
      real,         intent(in) :: zatom(ntype)
      real,         intent(in) :: omtil
      logical,      intent(in) :: l_soc,l_noco
      integer,      intent(in) :: neigd
      integer,      intent(in) :: fullnkpts,fullnqpts
      type(od_inp), intent(in) :: odi
      logical,      intent(in) :: l_proj
      integer, intent (in)     :: lmaxd,ntypd,nop,nvd,jspd
      integer, intent (in)     :: nbasfcn,llod,nlod,nwdd
      integer, intent (in)     :: nlo(ntypd),llo(nlod,ntypd)
      integer, intent (in)     :: lapw_l(ntypd)
      logical, intent (in)     :: l_dulo(nlod,ntypd),l_ss
      integer, intent (in)     :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      integer, intent (in)     :: lmax(ntypd)
      integer, intent (in)     :: invsat(natd),invsatnr(natd),nkpt(nwdd)
      real,    intent (in)     :: taual(3,natd),rmt(ntypd)
      real,    intent (in)     :: bbmat(3,3)
      real, intent (inout)     :: alph(ntypd),beta(ntypd)
      real, intent (inout)     :: qss(3,fullnqpts)
      real,    intent (in)     :: sk2(n2d),phi2(n2d)
      type (od_sym),intent(in) :: ods
      integer, intent (in)     :: irank,isize,n3d
      integer, intent (in) :: nmzxyd,nmzd,jmtd,nlhd,nq3,nvac
      logical, intent (in) :: invs,invs2
      integer, intent (in) :: ntypsd,nkptd
      integer, intent (in) :: nlh(ntypsd),jri(ntypd),ntypsy(natd)
      real,    intent (in) :: dx(ntypd)
      real,    intent (in) :: rmsh(jmtd,ntypd)
      real,    intent (in) :: e1s,e2s
      integer, intent (in) :: ulo_der(nlod,ntypd)
      complex, intent (in) :: ustep(n3d)
      integer, intent (in) :: k1d,k2d,k3d,nop2
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      real,    intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      logical, intent (in) :: slice
      integer, intent (in) :: kk,nnne
      integer, intent (in) :: nv2d,nmzxy,nmz,nq2
      real,    intent (in) :: delz,area,z1,volint
      integer, intent (in) :: ig2(n3d)
      real,    intent (in) :: tau(3,nop)
      logical, intent (in) :: symor
      real,    intent (in) :: ef
      logical, intent (in) :: l_sgwf
c      real,intent(inout)   :: fermi_weights(:,:,:) !(neigd,nkptd,jspd)

      character(len=12) :: fending
      integer :: i,kpts,ikpt,nkqpts,iqpt
      real(4) :: delta3,time_lapw_expand,delta2,time_lapw_plot
      logical :: l_need_fft,l_file
      integer :: hopmin_z,hopmax_z
      integer :: hopmin_y,hopmax_y
      integer :: hopmin_x,hopmax_x
      integer             :: ii
      integer             :: rvecnum,rvecind,num(3),int_dummy
      integer,allocatable :: rvec(:,:),ndegen(:)
      real                :: scale
      real                :: kpoints(3,fullnkpts)
      integer             :: r1,r2,r3,spinspin,num_angl
      logical             :: l_nocosoc

      l_nocosoc=.false.
      if(l_soc)l_nocosoc=.true.
      if(l_noco)l_nocosoc=.true.

c***************************************************
c     Read in the kpoints from w90kpts or kpts.
c***************************************************      
      call wann_get_kpts(
     >         wann%l_bzsym,film,odi%d1,.false.,
     <         kpts,kpoints)
      if (kpts /= fullnkpts)
     +  CALL juDFT_error ('mismatch in number of kpoints',
     +                  calledby="wann_postproc")

      call wann_get_kpts(
     >         wann%l_bzsym,film,odi%d1,.true.,
     <         kpts,kpoints)
c*********************************************************
c     Find out the structure of k-point set.
c*********************************************************
      call wann_get_mp(
     >         fullnkpts,kpoints,
     <         num)

      if(wann%l_wannierize.and.l_p0)then
#ifndef CPP_WANN
         WRITE(*,*) 'At this point a wannierization has to be performed'
         WRITE(*,*) 'but the Wannier90 library is not linked!'
         CALL juDFT_error ('Wannierization without Wannier90 library',
     +                     calledby="wann_postproc")
#else
         call wann_wannierize(
     >            film,wann%l_bzsym,jspins,
     >            natd,pos,amat,bmat,ntype,neq,zatom)
#endif
      endif

      if(wann%l_dipole2.and.l_p0)then
         call wann_dipole2(
     >            jspins,pos,omtil,natd,
     >            (l_soc.or.l_noco))
      endif   

      if(wann%l_dipole.and.l_p0)then
         call wann_dipole(
     >            jspins,omtil,natd,pos,amat,
     >            ntype,neq,zatom)
      endif

#ifdef CPP_TOPO
      if(wann%l_nedrho.and.l_p0)then
         call wann_nedrho(
     >       fullnkpts,area,ef,
     <       fermi_weights)
      endif
#endif

      if(wann%l_byenergy) i=3
      if(wann%l_bynumber) i=2
      if(wann%l_byindex)  i=1

      if(wann%l_plotw90)then
         CALL juDFT_error("not in this version",
     +                    calledby ="wann_postproc")

c         call wann_plotw90(i,wann%band_min,wann%band_max,numbands,nwfs,
c     >   l_dulo,l_noco,l_ss,lmaxd,ntypd,
c     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
c     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
c     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
c     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
c     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
c     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
c     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
c     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
c     >      volint,symor,pos,ef,wann%l_bzsym,irecl)

      endif

      l_need_fft=.false.
      if(wann%l_hopping.or.wann%l_nablars.or.
     &   wann%l_nablapaulirs.or.wann%l_pauli.or.
     &   wann%l_perpmagrs.or.wann%l_socmatrs.or.
     &   wann%l_torquers.or.wann%l_offdiposoprs.or.
     &   wann%l_socspicomrs.or.wann%l_spindisprs.or.
     &   wann%l_anglmomrs .or.wann%l_perturbrs.or.
     &   wann%l_orbcomprs.or.wann%l_rmat)l_need_fft=.true.

      if(l_need_fft.and.l_p0)then

       if(.false.)then !specify r-mesh by its boundaries
         hopmin_z=-5;hopmax_z=5
         hopmin_x=0;hopmax_x=0
         hopmin_y=0;hopmax_y=0
         rvecnum=(hopmax_z-hopmin_z+1)
         if(.not.odi%d1.and.film)then
           hopmin_x=-5;hopmax_x=5
           hopmin_y=-5;hopmax_y=5
           hopmin_z=0;     hopmax_z=0
         else
           hopmin_x=-5;hopmax_x=5
           hopmin_y=-5;hopmax_y=5
         endif
         rvecnum=        (hopmax_z-hopmin_z+1)
         rvecnum=rvecnum*(hopmax_y-hopmin_y+1)
         rvecnum=rvecnum*(hopmax_x-hopmin_x+1)

         allocate(rvec(3,rvecnum))
         rvecind=0
         do r3=hopmin_z,hopmax_z
          do r2=hopmin_y,hopmax_y
           do r1=hopmin_x,hopmax_x
            rvecind=rvecind+1
            if (rvecind > rvecnum)
     +        CALL juDFT_error ('mismatch in number of kpoints',
     +                         calledby="wann_postproc")

            rvec(1,rvecind)=r1
            rvec(2,rvecind)=r2
            rvec(3,rvecind)=r3
           enddo !r1
          enddo !r2
         enddo !r3
       else !determine optimal r-mesh
         call wann_wigner_seitz(
     >       .true.,num,amat,0,
     <       rvecnum,rvec,ndegen)
         allocate(rvec(3,rvecnum))
         allocate(ndegen(rvecnum))
         call wann_wigner_seitz(
     >       .false.,num,amat,rvecnum,
     <       int_dummy,rvec,ndegen)

         open(333,file='wig_vectors',recl=1000)
         do ii=1,rvecnum
            write(333,*)ii,rvec(1,ii),rvec(2,ii),rvec(3,ii),
     &             ndegen(ii)
         enddo
         close(333)

       endif

      endif

      if(wann%l_hopping.and.l_p0) then
         call wann_hopping(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_nocosoc,wann%band_min,wann%band_max,
     >       neigd,wann%l_socmmn0,wann%l_ndegen,ndegen,
     >       wann%wan90version,wann%l_unformatted)
      endif

      if(wann%l_nablars.and.l_p0) then
         call wann_nabla_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,
     >       neigd,.false.,
     >       wann%wan90version)
      endif  

      if(wann%l_orbcomprs.and.l_p0)then
         num_angl=9
         if(wann%l_oc_f)num_angl=16
         call wann_fft5(
     >          .false.,
     >          num_angl,wann%oc_num_orbs,
     >          'WF1.orbcomp','orbcomp.1',
     >          rvecnum,rvec,kpoints,
     >          1,fullnkpts,film,
     >          l_soc,wann%band_min,wann%band_max,neigd,
     >          wann%wan90version)
         if( jspins.eq.2 )then
            spinspin=2
            if(l_soc.or.l_noco)spinspin=1
            call wann_fft5(
     >             .false.,
     >             num_angl,wann%oc_num_orbs,
     >             'WF2.orbcomp','orbcomp.2',
     >             rvecnum,rvec,kpoints,
     >             spinspin,fullnkpts,film,
     >             l_soc,wann%band_min,wann%band_max,neigd,
     >             wann%wan90version)
         endif   
      endif
 
      if(wann%l_rmat)then
        call wann_rmat(
     >      bmat,amat,
     >      rvecnum,rvec,kpoints,
     >      jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >      l_nocosoc,wann%band_min,wann%band_max,
     >      neigd,.false.,wann%wan90version)
      endif


       if(wann%l_anglmomrs.and.l_p0)then
        call wann_fft4(
     >       'WF1.anglmom',
     >       'anglmomrs.1',.false.,
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.,wann%wan90version)
       endif
 
#ifdef CPP_TOPO
      if(wann%l_offdiposoprs.and.l_p0)then
         call wann_offdiposop_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.)
      endif

      if(wann%l_spindisprs.and.l_p0)then
          call wann_fft5(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.)
      endif   

      if(wann%l_perturbrs.and.l_p0)then
          call wann_fft3(
     >       'WF1.perturb'  ,'perturbrs.1'  ,.false.,
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.)
      endif

      if(wann%l_socspicomrs.and.l_p0)then
        if(l_soc)then
          call wann_fft4(
     >       'WF1.socspicom','socspicomrs.1',.true.,
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.)
        else
          call wann_fft6(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,neigd,
     >       .false.)
        endif  
      endif



      if(wann%l_torquers.and.l_p0) then
         call wann_torque_rs(
     >       ntype,neq,rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,
     >       neigd,.false.)
      endif  
#endif
      if (wann%l_nablapaulirs.and.l_p0)then
         call wann_nabla_pauli_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,
     >       neigd,.false.,wann%wan90version)
      endif   

      if (wann%l_pauli.and.l_p0)then
         call wann_pauli_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_nocosoc,wann%band_min,wann%band_max,
     >       neigd,.false.,wann%l_ndegen,ndegen,
     >       wann%wan90version,wann%l_unformatted)
      endif

      if (wann%l_perpmagrs.and.l_p0)then
         call wann_perpmag_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,
     >       neigd,.false.,wann%l_ndegen,ndegen,wann%wan90version,
     >       wann%l_unformatted)
      endif

      if (wann%l_socmatrs.and.l_p0)then
         call wann_socmat_rs(
     >       rvecnum,rvec,kpoints,
     >       jspins,fullnkpts,wann%l_bzsym,film,odi%d1,
     >       l_soc,wann%band_min,wann%band_max,
     >       neigd,.false.,wann%wan90version)
      endif

      if(wann%l_plot_umdat)then
         call wann_plot_um_dat(
     >      DIMENSION,stars,vacuum,atoms,sphhar,input,sym,mpi,
     >      lapw,oneD,noco,cell,eig_id,
     >      mpi_comm,i,wann%band_min,wann%band_max,l_soc,
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss(:,1),sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,wann%l_bzsym,wann%l_proj_plot,
     >      wann%wan90version)
      endif


      call cpu_time(delta2)
      if(wann%l_lapw.and.l_p0)then
         call wannier_to_lapw(
     >      mpi_comm,eig_id,
     >      input,lapw,oneD,noco,sym,cell,atoms,stars,vacuum,sphhar,
     >      l_soc,wann%unigrid,i,wann%band_min,wann%band_max,
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss(:,1),sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,wann%l_bzsym,wann%l_proj_plot,
     >      wann%wan90version)
      endif

      if(wann%l_lapw_kpts)then
#ifdef DCPP_WANN_EXT
         call wannier_to_lapw_kpts(
     >      unigrid,i,wann%band_min,wann%band_max,
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss(:,1),sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,wann%l_bzsym,wann%l_proj_plot,irecl)
#else
         CALL juDFT_error("not yet tested in this release",calledby
     +        ="wann_postproc")
#endif
      endif
      if(wann%l_lapw_gfleur)then
#ifdef DCPP_WANN_EXT
         call wannier_lapw_gfleur(
     >      gfthick,gfcut,i,wann%band_min,wann%band_max,
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss(:,1),sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,wann%l_bzsym,wann%l_proj_plot,irecl)
#else
         CALL juDFT_error("not yet tested in this release",calledby
     +        ="wann_postproc")
#endif
      endif

      call cpu_time(delta3)
      time_lapw_expand=delta3-delta2

      call cpu_time(delta2)
      if(wann%l_plot_lapw.and.l_p0)then
         CALL juDFT_error("not yet tested in this release",calledby
     +        ="wann_postproc")
c         call wann_plot_from_lapw(
c     >     nv2d,jspins,odi,ods,n3d,nmzxyd,n2d,
c     >     ntypsd,
c     >     ntype,lmaxd,jmtd,natd,nmzd,neq,nq3,nvac,
c     >     nmz,nmzxy,nq2,nop,nop2,volint,film,slice,symor,
c     >     invs,invs2,z1,delz,ngopr,ntypsy,jri,pos,zatom,
c     >     lmax,mrot,tau,rmsh,invtab,amat,bmat,bbmat,nnne,kk,
c     >     nlod,llod,lmd,omtil,nlo,llo)
      endif
      call cpu_time(delta3)
      time_lapw_plot=delta3-delta2

      write(6,*)"time_lapw_expand=",time_lapw_expand
      write(6,*)"time_lapw_plot=",time_lapw_plot

      if(wann%l_finishnocoplot.and.l_p0) then
c         write(*,*)'doing the UNK mixing'
c         write(*,*)'alph',alph
c         write(*,*)'beta',beta
c         nkqpts=fullnkpts
c         if(l_sgwf) nkqpts=fullnkpts*fullnqpts 
         call wann_nocoplot(atoms,slice,nnne,!nslibd
     >                      amat,bmat,fullnkpts,odi,film,
     >                      natd,ntypd,jmtd,ntype,neq,pos,
     >                      jri,rmsh,alph,beta,fullnqpts,qss,
     >                      z1,zatom)
      endif

      end subroutine wann_postproc
      end module m_wann_postproc
