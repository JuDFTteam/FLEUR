!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_plotw90
      use m_juDFT

c******************************************************************
c       plot wannierfunctions directly within fleur
c       based on wann_plot
c       FF, September 2006
c       last update: 27.10.2006 by FF
c******************************************************************
      CONTAINS
      SUBROUTINE wann_plotw90(bandswitch,band_min,band_max,num_bands,
     >  num_wann,l_dulo,l_noco,l_ss,lmaxd,
     >      ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,l_bzsym,irecl)

      use m_abcof
      use m_wann_2dvacabcof
      use m_dotir
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0, cdn_read
      use m_od_types, only : od_inp, od_sym
      use m_loddop
      use m_constants, only : pimach
      use m_wann_real
      use m_xsf_io
      use m_cotra,only:cotra1
      use m_wann_plot_vac
      USE m_wann_plot_wannier
      IMPLICIT NONE
 
      integer,intent(in)   :: bandswitch
      integer,intent(in)   :: band_min
      integer,intent(in)   :: band_max
      integer, intent (in) :: num_wann
      logical, intent (in) :: invs,invs2,film,slice,symor
      integer, intent (in) :: lmaxd,ntypd,neigd,nkptd,kk,nnne
      integer, intent (in) :: natd,nop,nvd,jspd,nbasfcn,nq2,nop2
      integer, intent (in) :: llod,nlod,ntype,nwdd,n3d,n2d
      integer, intent (in) :: nmzxyd,nmzd,jmtd,nlhd,nq3,nvac 
      integer, intent (in) :: ntypsd,jspins,k1d,k2d,k3d
      real,    intent (in) :: omtil,e1s,e2s,delz,area,z1,volint
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      real,    intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: nlh(ntypsd),jri(ntypd),ntypsy(natd)
      integer, intent (in) :: nlo(ntypd),llo(nlod,ntypd),lapw_l(ntypd)
      integer, intent (in) :: invtab(nop),mrot(3,3,nop),ngopr(natd)
      integer, intent (in) :: neq(ntypd),lmax(ntypd)
      integer, intent (in) :: invsat(natd),invsatnr(natd),nkpt(nwdd)
      integer, intent (in) :: irank,isize,nv2d,nmzxy,nmz
      integer, intent (in) :: ulo_der(nlod,ntypd),ig2(n3d)
      real,    intent (in) :: taual(3,natd),rmt(ntypd),dx(ntypd) 
      real,    intent (in) :: amat(3,3),bmat(3,3),bbmat(3,3)
      real,    intent (in) :: rmsh(jmtd,ntypd),tau(3,nop),zatom(ntype)
      real, intent (inout) :: alph(ntypd),beta(ntypd),qss(3)  
      real,    intent (in) :: pos(3,natd),ef
      complex, intent (in) :: ustep(n3d)
      logical, intent (in) :: l_dulo(nlod,ntypd),l_noco,l_ss,l_bzsym
c-odim
      real,    intent (in) :: sk2(n2d),phi2(n2d)
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
c+odim
      integer, intent (in) :: irecl
cccccccccccccccccc   local variables   cccccccccccccccccccc
      integer lmd,nlotot,n,nmat,nw,ispin,iter,ikpt,ilo
      integer noccbd,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer jsp_start,jsp_end,nrec,nrec1,nbands
      integer nodeu,noded,n_size,na,n_rank,nbnd,nkbnd
      integer i1,i2,i3,in,ikpt_k,lda,num_bands
      integer n_bands(0:neigd),nslibd
      character*8 dop,iop,name(10)
      real bkpt(3),sfp,tpi,wronk,wk,wk_b,phase
      real eig(neigd),cp_time(9)
      logical l_p0,l_bkpts,l_proj,l_amn,l_mmn,l_eig,wann,wann_plott
!!! energy window boundaries
      integer, allocatable :: kveclo(:),nv(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      real,    allocatable :: we(:)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, allocatable :: z(:,:),zz(:,:)
#else
      real,    allocatable :: z(:,:),zz(:,:)
#endif
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)
#endif
      real,    allocatable :: eigg(:)
      real kpoints(nkptd)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:)
      complex, allocatable :: bcof(:,:,:)
      complex, allocatable :: ccof(:,:,:,:)
!!! the parameters for the number of wfs
      integer :: nwfs
!!! the potential in the spheres and the vacuum
      real, allocatable :: vr(:,:,:),vz(:,:,:),vrf(:,:,:,:)
!!! auxiliary potentials
      complex, allocatable :: vpw(:,:),vzxy(:,:,:,:)
!!! bkpts data
      integer nntot,ikpt_help
      integer, allocatable :: gb(:,:,:),bpt(:,:)
!!! radial wavefunctions in the muffin-tins and more ...
      real,    allocatable :: flo(:,:,:,:)
      real,    allocatable :: ff(:,:,:,:),gg(:,:,:,:)
      real,    allocatable :: us(:,:),uds(:,:),ulos(:,:)
      real,    allocatable :: dus(:,:),duds(:,:),dulos(:,:)
      real,    allocatable :: ddn(:,:),uulon(:,:),dulon(:,:)
      real,    allocatable :: uloulopn(:,:,:)
      real     uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real     ulouilopn(nlod,nlod,ntypd)
!!! energy parameters
      real    ello(nlod,ntypd,jspd),evac(2,jspd)
      real    epar(0:lmaxd,ntypd,jspd),evdu(2,jspd)
      character(len=3) :: spin12(2)
      data spin12/'WF1' , 'WF2'/
      character(len=2)::spinspin12(2)
      data spinspin12/'.1','.2'/
      complex,allocatable::wannierfunc(:,:,:)
      integer num_kpts,num_nnmax
      integer posi
      integer n2
      real tmp_omi
      integer kpt,oper
      real poinint(3)
      real phas,tmax
      real bkrot(3)
      integer j1,j2,j3
      real v(3),scale,ev
      complex av,bv
      logical have_disentangled,l_chk,l_umdat
      integer :: chk_unit,nkp,ntmp,fullnkpts
      integer,allocatable::irreduc(:),mapkoper(:)
      character(len=20)::checkpoint
      real :: tmp_latt(3,3)
      real,allocatable:: tmp_kpt_latt(:,:)
      real omega_invariant
      complex,allocatable::u_matrix_opt(:,:,:)
      logical l_file
      logical,allocatable:: inc_band(:)
      integer num_inc,counter,kptibz
      logical l_bandsort


      REAL          :: s
      COMPLEX       :: xdnout,factor
      INTEGER       :: ii3,ix,iy,iz,nplo,nbn
      INTEGER       :: nbmin,nbmax
      INTEGER       :: nplot,nq,nt,jm,ii1,ii2
      LOGICAL       :: twodim
      logical l_bynumber,l_byindex,l_byenergy

c     ..basis wavefunctions in the vacuum
      complex, allocatable :: ac(:,:,:),bc(:,:,:)
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      real,    allocatable :: u(:,:,:),ue(:,:,:)
      real vz0(2)

      integer ik,nv2,ivac,jvac,symvac,symvacvac
      real evacp,sign,arg,const
      complex c_1
      integer kvac1(nv2d),kvac2(nv2d),map2(nvd)
      real fas,zks
      complex d_wgn(-3:3,-3:3,3,nop)
      integer pos1,pos2,ato,loc,invop


      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      INTEGER :: grid(3),k
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(jmtd)
      REAL    :: point(3),post(3,natd)
      CHARACTER(len=30):: filename
      CHARACTER(len=20):: name1,name2,name3
      CHARACTER(len=10):: vandername
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename 
c      external dotirp !module now
c      real dotirp

      l_byindex=.false.
      l_byenergy=.false.
      l_bynumber=.false.
      if (bandswitch==1)l_byindex=.true.
      if(bandswitch==2)l_bynumber=.true.
      if(bandswitch==3)l_byenergy=.true.
      


c     read in plot_inp

      INQUIRE(file ="plot_inp",exist= twodim)
      IF (.NOT.twodim) THEN !no input file exists, create a template and
                            !exit
         OPEN(20,file ="plot_inp")
         WRITE(20,'(i2,a5,l1)') 1,",xsf=",.false.
c         WRITE(20,*) "&PLOT twodim=t,cartesian=t"
c         WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
c         WRITE(20,*) "  filename='plot1' /"
         WRITE(20,*) "&PLOT twodim=f,cartesian=f"
         WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
         WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
         WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
         WRITE(20,*) "  grid(1)=20  grid(2)=20  grid(3)=20  "
c         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.5 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No plot_inp file found. Created a template"
         CALL juDFT_error("Missing input for plot; modify plot_inp"
     +        ,calledby ="wann_plotw90")
      ENDIF


      


      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)') nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (nplot.ge.2) 
     &     CALL juDFT_error
     +     ("plots one by one, please, this is not charge density"
     +     ,calledby ="wann_plotw90")
      twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
      vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
      zero = (/0.,0.,0./);filename="default"
      READ(18,plot)
      IF (twodim.AND.ANY(grid(1:2)<1)) 
     +     CALL juDFT_error("Illegal grid size in plot",calledby
     +     ="wann_plotw90")
      IF (.NOT.twodim.AND.ANY(grid<1)) 
     +     CALL juDFT_error("Illegal grid size in plot",calledby
     +     ="wann_plotw90")
      IF (twodim) grid(3) = 1
      !calculate cartesian coordinates if needed
      IF (.NOT.cartesian) THEN
         vec1=matmul(amat,vec1)
         vec2=matmul(amat,vec2)
         vec3=matmul(amat,vec3)
         zero=matmul(amat,zero)
      ENDIF
      Close(18)
    


      sfp = 2* sqrt( pimach() )
      tpi = 2* pimach()
      lmd = lmaxd*(lmaxd+2)
      nkpts = maxval(nkpt(:))

      nrec = 0
      nlotot = 0
      do n = 1, ntype
        do l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        enddo
      enddo

c      if((nlotot.ne.0).and.l_bzsym) then
c         call d_wigner(nop,mrot,bmat,3,d_wgn)
c      endif   



cccccccccccccccc   initialize the potential   cccccccccccc

      allocate ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      allocate ( vz(nmzd,2,jspd),vrf(jmtd,0:nlhd,ntypd,jspd) )
      allocate ( vr(jmtd,ntypd,jspd) )

      open (8,file='pottot',form='unformatted',status='old')
      rewind (8)

      call loddop(
     >      jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd,
     >      jspins,nq3,odi%nq2,nvac,ntype,invs,invs2,film,
     >      nlh,jri,ntypsd,ntypsy,8,natd,neq,
     <      iop,dop,iter,vrf,vpw,vz,vzxy,name)

      close (8)

      do jspin = 1,jspins
        do n = 1, ntype
          do j = 1,jri(n)
            vr(j,n,jspin) = vrf(j,0,n,jspin)
          enddo
        enddo
      enddo

      deallocate ( vpw,vzxy,vrf )
 
cccccccccccccccc   end of the potential part  ccccccccccc

      allocate ( kveclo(nlotot),nv(jspd) )
      allocate ( k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )
      allocate ( ff(ntypd,jmtd,2,0:lmaxd) )
      allocate ( gg(ntypd,jmtd,2,0:lmaxd) )
      allocate ( us(0:lmaxd,ntypd) )
      allocate ( uds(0:lmaxd,ntypd) )
      allocate ( dus(0:lmaxd,ntypd) )
      allocate ( duds(0:lmaxd,ntypd) )
      allocate ( ddn(0:lmaxd,ntypd) )
      allocate ( ulos(nlod,ntypd) )
      allocate ( dulos(nlod,ntypd) )
      allocate ( uulon(nlod,ntypd) )
      allocate ( dulon(nlod,ntypd) )
      allocate ( uloulopn(nlod,nlod,ntypd) )

c cycle by spins starts! NO NON-COLLINEAR WHATSOEVER!

      do 110 jspin=1,jspd   ! cycle by spins
         print*,"spin=",jspin
      jsp_start = jspin ; jsp_end = jspin

c**************************************************************
c   for bzsym=.true.: determine mapping between kpts and w90kpts
c**************************************************************
      if (l_bzsym) then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby
     +        ="wann_plotw90")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts
c         if(l_matrixmmn)then
c            allocate(kpoi(3,fullnkpts))
c            do k=1,fullnkpts
c               read(412,*)kpoi(:,k)
c            enddo   
c         endif   
         close(412)
         print*,"fullnkpts=",fullnkpts
         IF(fullnkpts<=nkpts) CALL juDFT_error("fullnkpts.le.nkpts"
     +        ,calledby ="wann_plotw90")
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         l_file=.false.
         inquire(file='kptsmap',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("kptsmap not found, needed if bzsym",calledby
     +        ="wann_plotw90")
         open(713,file='kptsmap')
         do i=1,fullnkpts
            read(713,*)kpt,irreduc(i),mapkoper(i)
            IF(kpt/=i) CALL juDFT_error("kpt.ne.i",calledby
     +           ="wann_plotw90")
            print*,i,irreduc(i),mapkoper(i)
         enddo   
         close(713)
         IF(MAXVAL(irreduc(:))/=nkpts) CALL juDFT_error
     +        ("max(irreduc(:))/=nkpts",calledby ="wann_plotw90")
      else   
         fullnkpts=nkpts
c         if(l_matrixmmn)then
c            allocate(kpoi(3,fullnkpts))
c            open(412,file='kpts',form='formatted')
c            read(412,*)
c            do k=1,fullnkpts
c               read(412,*)kpoi(:,k)
c            enddo   
c         endif   
      endif   

c      if(l_matrixmmn) then
c         print*,kpoi(:,:)
c      endif
c**************************************************************








      if(l_byindex) num_bands=band_max-band_min+1

      if(l_byenergy)then   !determine number of bands
         l_file=.false.
         inquire(file=spin12(jspin)//'.eig',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error("WF1(2).eig not found",calledby
     +        ="wann_plotw90")
         open(396,file=spin12(jspin)//'.eig',form='formatted',
     &    status='old')
         k=0
         do i=1,(fullnkpts*num_wann)
            read(396,*)j
            if(j.gt.k)then
               k=j
            elseif(k.gt.j)then
               num_bands=i-1
               exit
            else
               CALL juDFT_error("k == j in wann_plotw90.F",calledby
     +              ="wann_plotw90")
            endif   
         enddo   
         close(396)
      endif   

      print*,"number of bands=",num_bands
   
      allocate(wannierfunc(
     &  (grid(1))*(grid(2))*(grid(3)),num_bands,fullnkpts))
      allocate(inc_band(num_bands))

cccccccccccc   read in the eigenvalues and vectors   cccccc

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.

      call cdn_read0(
     >               lmaxd,ntypd,nlod,neigd,jspd,
     >               irank,isize,jspin,jsp_start,jsp_end,
     >               l_noco,nrec,66,
     <               ello,evac,epar,bkpt,wk,n_bands,nrec1,n_size)


      allocate ( flo(ntypd,jmtd,2,nlod) )

      na = 1
      do 40 n = 1,ntype
       do 30 l = 0,lmax(n)
c...compute the l-dependent, k-independent radial MT- basis functions
         call radfun(
     >              l,epar(l,n,jspin),vr(1,n,jspin),jri(n),rmsh(1,n),
     >              dx(n),jmtd,
     <              ff(n,:,:,l),gg(n,:,:,l),us(l,n),
     <              dus(l,n),uds(l,n),duds(l,n),
     <              ddn(l,n),nodeu,noded,wronk)
   30  continue
c...and the local orbital radial functions
c       do ilo = 1, nlo(n)
         call radflo(
     >             ntypd,nlod,jspd,jmtd,lmaxd,n,jspin,
     >             ello(1,1,jspin),vr(1,n,jspin),
     >             jri(n),rmsh(1,n),dx(n),ff(n,1:,1:,0:),
     >             gg(n,1:,1:,0:),llo,nlo,l_dulo(1,n),irank,ulo_der,
     <             ulos(1,1),dulos(1,1),uulon(1,1),dulon(1,1),
     <             uloulopn(1,1,1),uuilon,duilon,ulouilopn,flo(n,:,:,:))
c       enddo
c       na = na + neq(n)
   40 continue
      i_rec = 0 ; n_rank = 0

      do ikpt = 1,fullnkpts  ! loop by k-points starts
         print*,"k-point=",ikpt
        kptibz=ikpt 
        if(l_bzsym) kptibz=irreduc(ikpt)
        if(l_bzsym) oper=mapkoper(ikpt)





      i_rec = i_rec + 1

      if (mod(i_rec-1,isize).eq.irank) then

      allocate ( z(nbasfcn,neigd),we(neigd),
     &                   zz(nbasfcn,neigd),eigg(neigd) )

      call cdn_read(
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,kptibz,66,neigd,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eigg,zz,cp_time)

      nslibd = 0

c...we work only within the energy window

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      z(:,:) = cmplx(0.,0.)
#else
      z(:,:) = 0.
#endif
      eig(:) = 0.

      print*,"bands used"
      do i = 1,nbands
         if ((eigg(i).ge.e1s.and.nslibd.lt.num_bands.and.l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.l_byenergy)
     &.or.(i.ge.band_min.and.i.le.band_max.and.l_byindex))then
            print*,i
        nslibd = nslibd + 1
        eig(nslibd) = eigg(i)
        we(nslibd) = we(i)
        do  j = 1,nv(jspin) + nlotot
         z(j,nslibd) = zz(j,i)
        enddo
       endif 
      enddo 
c***********************************************************
c              rotate the wavefunction
c***********************************************************
      if (l_bzsym.and.oper.ne.1) then  !rotate bkpt
         bkrot(:)=0.0
         do k=1,3
           bkrot(:)=bkrot(:)+mrot(k,:,oper)*bkpt(k)
         enddo
         bkpt(:)=bkrot(:)

         jloop:do j=1,nv(jspin)
               j1=mrot(1,1,oper)*k1(j,jspin)+
     +             mrot(2,1,oper)*k2(j,jspin)+mrot(3,1,oper)*k3(j,jspin)
               j2=mrot(1,2,oper)*k1(j,jspin)+
     +             mrot(2,2,oper)*k2(j,jspin)+mrot(3,2,oper)*k3(j,jspin)  
               j3=mrot(1,3,oper)*k1(j,jspin)+
     +             mrot(2,3,oper)*k2(j,jspin)+mrot(3,3,oper)*k3(j,jspin)
               k1(j,jspin)=j1
               k2(j,jspin)=j2
               k3(j,jspin)=j3
         enddo jloop  

      endif   
      print*,"bkpt=",bkpt(:)
c******************************************************************


      noccbd = nslibd

      allocate ( acof(noccbd,0:lmd,natd),
     &           bcof(noccbd,0:lmd,natd),
     &           ccof(-llod:llod,noccbd,nlod,natd))

      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:) = cmplx(0.,0.)

c...generation of the A,B,C coefficients in the spheres 
c...for the lapws and local orbitals, summed by the basis functions 

      call abcof(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,jspd,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,noccbd,z,
     >           us(0,1),dus(0,1),uds(0,1),
     >           duds(0,1),ddn(0,1),invsat,invsatnr,
     >           ulos(1,1),uulon(1,1),dulon(1,1),
     >           dulos(1,1),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,jspin,alph,beta,qss,kveclo,odi,ods,
     <           acof(1,0,1),bcof(1,0,1),
     <           ccof(-llod,1,1,1))

      call abcrot(
     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >        noccbd,lmax,nlo,llo,nop,ngopr,mrot,invsat,invsatnr,
     >                 bmat,odi,ods,
     X                 acof,bcof,ccof)


c***********************************************************************
c make preparations for plotting in vacuum
c***********************************************************************
      if (film)then
         allocate ( ac(nv2d,nslibd,2),bc(nv2d,nslibd,2),
     +         u(nmzd,nv2d,nvac),ue(nmzd,nv2d,nvac))
         call wann_2dvacabcof(
     >         nv2d,nslibd,nvac,nmzd,nmz,omtil,vz(:,:,jspin),
     >         nv(jspin),bkpt,z1,
     >         nvd,k1(:,jspin),k2(:,jspin),k3(:,jspin),evac(:,jspin),
     >         bbmat,delz,bmat,nbasfcn,
     >         neigd,z,
     <         ac,bc,u,ue)
      endif !preparations for vacuum
c**************************************************************************
c**************************************************************************
      nbmin=1
      nbmax=nslibd


      band:DO nbn = nbmin,nbmax

          

          DO iz = 0,grid(3)-1
          DO iy = 0,grid(2)-1
           xloop:DO ix = 0,grid(1)-1
            posi=ix+1+iy*(grid(1))+iz*(grid(1))*(grid(2))
            point = zero+vec1*(ix+0.0)/(grid(1))+vec2*(iy+0.0)
     $                 /(grid(2))
            IF (.NOT.twodim) point = point+vec3*(iz+0.0)/(grid(3))

!Check if the point is in MT-sphere
            if(film)then
               fas=-bkpt(3)*bmat(3,3)*z1
               factor=cmplx(cos(fas),sin(fas))
             else
                factor=cmplx(1.0,0.0)
             endif   
             ii1 = 3
             ii2 = 3
             ii3 = 3
             IF (film .AND. .NOT.odi%d1) ii3 = 0
             IF (odi%d1) THEN
                ii1 = 0 ; ii2 = 0
             END IF
             DO  i1 = -ii1,ii1
              DO  i2 = -ii2,ii2
               DO  i3 = -ii3,ii3
                pt = point+MATMUL(amat,(/i1,i2,i3/))
                na = 0
                DO nt = 1,ntype
                 DO nq = 1,neq(nt)
                  na   = na + 1
                  s  = SQRT(dot_PRODUCT(pos(:,na)-pt,pos(:,na)-pt))
                  IF (s<rmsh(jri(nt),nt)) THEN
                    CALL wann_real(
     >                   pt,nt,na,0,1,bkpt,
     >                   n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >                   natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >                   nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >                   lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >                   omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >                   ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >                   ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >                   lmd,nbasfcn,
     <                   xdnout)
                    wannierfunc(posi,nbn,ikpt)=xdnout*factor


                   CYCLE xloop
                  ENDIF
                 ENDDO
                ENDDO !nt
               ENDDO
              ENDDO
             ENDDO !i1
!Check for point in vacuum
             IF (film.AND..NOT.odi%d1.AND.ABS(point(3))>=z1) THEN
                ivac=1
                if (point(3).lt. 0.0)ivac=2
                jvac=ivac
                if(nvac==1)jvac=1
                call wann_plot_vac(point,z1,nmzd,nv2d,n3d,nvac,
     >            nmz,delz,bmat,bbmat,evac(:,jspin),bkpt,vz,jspin,
     >            k1(:,jspin),k2(:,jspin),k3(:,jspin),nvd, 
     >            nbasfcn,neigd,z(:,nbn),nv(jspin),omtil,nslibd,
     >            ac(:,nbn,ivac),
     &             bc(:,nbn,ivac),
     &            u(:,:,jvac),ue(:,:,jvac),xdnout)
                xdnout=xdnout*factor
                if(real(xdnout).gt.9.0 .or.real(xdnout).lt.-9.0
     &        .or.imag(xdnout).gt.9.0 .or. imag(xdnout).lt.-9.0)then
                xdnout=cmplx(0.0,0.0)
                print*,"vac-problem at z=",point(3)
                endif
c               CALL wann_real(
c     >              point,0,0,1,0,bkpt,
c     >              n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
c     >              natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
c     >              nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
c     >              lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
c     >              omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
c     >              ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
c     >              ccof(:,nbn,:,:),z(:,nbn),
c     >              nv,k1,k2,k3,lmd,nbasfcn,
c     <              xdnout)
                wannierfunc(posi,nbn,ikpt)=xdnout*factor
              CYCLE xloop
             END IF
            
             IF (odi%d1) THEN
              IF (SQRT((pt(1))**2 + (pt(2))**2)>=z1) THEN
                CALL wann_real(
     >               pt,0,0,1,0,bkpt,
     >               n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >               natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >               nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >               lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >               omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >               ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >               ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >                lmd,nbasfcn,
     <               xdnout)
                wannierfunc(posi,nbn,ikpt)=xdnout

                CYCLE xloop
              END IF
             END IF
             CALL wann_real(
     >             point,0,0,0,2,bkpt,
     >             n3d,nmzxyd,n2d,ntypsd,lmaxd,jmtd,
     >             natd,ntypd,nmzd,nop,nop2,mrot,tau,invtab,
     >             nq3,nvac,invs,z1,delz,nmz,nmzxy,nq2,
     >             lmax,rmsh,jri,pos,ngopr,ntypsy,nvd,
     >             omtil,amat,bmat,odi,ods,nlod,llod,nlo,llo,
     >             ff,gg,flo,acof(nbn,:,:),bcof(nbn,:,:),
     >             ccof(:,nbn,:,:),z(:,nbn),
     >               nv(jspin),k1(:,jspin),k2(:,jspin),k3(:,jspin),
     >             lmd,nbasfcn,
     <             xdnout)
             wannierfunc(posi,nbn,ikpt)=xdnout*factor


            ENDDO xloop
           ENDDO
          ENDDO !z-loop

c..end of the loop by the bands



      ENDDO band   

      deallocate ( acof,bcof,ccof,z,we,zz,eigg )

      write (*,*) 'nslibd=',nslibd

      if(film)then
         deallocate(ac,bc,u,ue)
      endif   

      endif!processors

      enddo !loop over k-points

      deallocate(flo)
      if(l_bzsym)deallocate(mapkoper,irreduc)
      call wann_plot_wannier(wannierfunc(1,1,1),grid(1),
     &    jspin,fullnkpts,num_bands)
      
      deallocate(wannierfunc)
      nrec=nrec+nkpts
      deallocate(inc_band)

110   continue ! end of cycle by spins


      deallocate ( vr,vz,kveclo,nv,k1,k2,k3 )
      deallocate ( ff,gg,us,dus,duds,uds,ddn)
      deallocate ( ulos,dulos,uulon,dulon,uloulopn )


      END SUBROUTINE wann_plotw90
      END MODULE m_wann_plotw90
