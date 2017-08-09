!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

c******************************************************************
c       Calculate lapw-representation of wannierfunctions.
c       Lapw means here: lapw-like, but not identical to
c       fleur-lapw. Wannierfunctions can not be expanded
c       in fleur-lapw-set of a single k-point in general.
c
c       Frank Freimuth, November 2006
c******************************************************************
      module m_wannier_to_lapw
      use m_juDFT
      contains
      subroutine wannier_to_lapw(
     >      mpi_comm,eig_id,
     >      input,lapw,oneD,noco,sym,cell,atoms,stars,vacuum,sphhar,
     >      l_soc,unigrid,sortrule,band_min,band_max,
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,l_bzsym,l_proj_plot,wan90version)
      use m_wann_rw_eig
      use m_wann_read_umatrix
      use m_abcof
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0
      use m_types
      use m_loddop
      use m_constants
      use m_wann_real
      use m_xsf_io
      use m_wann_plot_vac
      USE m_abcrot


      implicit none
#include "cpp_double.h"

      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_lapw),INTENT(IN)   :: lapw
      TYPE(t_oneD),INTENT(IN)   :: oneD
      TYPE(t_noco),INTENT(IN)   :: noco
      TYPE(t_sym),INTENT(IN)    :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_stars),INTENT(IN)  :: stars
      TYPE(t_vacuum),INTENT(IN) :: vacuum
      TYPE(t_sphhar),INTENT(IN) :: sphhar

#ifdef CPP_MPI
      include 'mpif.h'
      integer mpiierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)
#endif
      logical,intent(in):: l_soc
      integer,intent(in)::unigrid(4),mpi_comm,eig_id
      integer,intent(in)::band_min(2),band_max(2)
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
      logical,intent(in)::l_proj_plot
      integer,intent(in)::sortrule
      integer, intent(in):: wan90version
c-odim
      real,    intent (in) :: sk2(n2d),phi2(n2d)
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
c+odim
      logical l_spreadcal
      complex, allocatable::spreads(:,:)
      real,allocatable:: centers(:,:)
cccccccccccccccccc   local variables   cccccccccccccccccccc
      integer lmd,nlotot,n,nmat,nw,ispin,iter,ikpt,ilo
      integer :: wannierspin,jspin2
      integer noccbd,nn,nkpts,i,jspin,j,l,i_rec,m,nwf,nwfp
      integer jsp_start,jsp_end,nrec,nrec1,nbands
      integer nodeu,noded,n_size,na,n_rank,nbnd,nkbnd
      integer i1,i2,i3,in,ikpt_k,lda
      integer n_bands(0:neigd),nslibd
      character*8 dop,iop,name(10)
      real bkpt(3),wronk,wk,wk_b,phase
      real eig(neigd),cp_time(9)
      logical l_p0,l_bkpts,l_proj,l_amn,l_mmn,l_eig,wann,wann_plott
!!! energy window boundaries
      integer, allocatable :: kveclo(:),nv(:)
      integer, allocatable :: k1(:,:),k2(:,:),k3(:,:)
      real,    allocatable :: we(:)

      real,    allocatable :: eigg(:)
      real kpoints(nkptd)
!!! a and b coeff. constructed for each k-point
      complex, allocatable :: acof(:,:,:)
      complex, allocatable :: bcof(:,:,:)
      complex, allocatable :: ccof(:,:,:,:)
      complex, allocatable :: wann_acof(:,:,:,:,:,:)
      complex, allocatable :: wann_bcof(:,:,:,:,:,:)
      complex, allocatable :: wann_ccof(:,:,:,:,:,:,:)
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

      real     uuilon(nlod,ntypd),duilon(nlod,ntypd)
      real     ulouilopn(nlod,nlod,ntypd)
!!! energy parameters
      real    ello(nlod,ntypd,max(2,jspd)),evac(2,max(2,jspd))
      real    epar(0:lmaxd,ntypd,max(2,jspd)),evdu(2,max(jspd,2))
      character(len=3) :: spin12(2)
      data spin12/'WF1' , 'WF2'/
      complex,allocatable::wannierfunc(:,:)
      complex,allocatable::wannierfunc_temp(:,:)
      integer posi
      complex,allocatable::u_matrix(:,:,:)
      integer kpt,oper
      real poinint(3)
      real phas,tmax
      real bkrot(3)
      integer j1,j2,j3
      logical um_format
      logical have_disentangled
      integer,allocatable :: ndimwin(:)
      logical,allocatable :: lwindow(:,:)
      integer :: chk_unit,nkp,ntmp,ierr,fullnkpts
      integer,allocatable::irreduc(:),mapkoper(:)
      logical l_file
      logical,allocatable::inc_band(:)
      integer num_inc,counter,kptibz
      logical l_byindex, l_byenergy, l_bynumber
      integer num_wann,num_bands,kpun
      complex d_wgn(-3:3,-3:3,3,nop)
      integer pos1,pos2,ato,loc,invop
      real vz0(2)
      integer ik,nv2,ivac,jvac,symvac,symvacvac
      real evacp,sign,arg
      complex c_1
      integer kvac1(nv2d),kvac2(nv2d),map2(nvd)
      real fas,zks
      integer mesh
      integer n2
      real v(3),scale,ev
      complex av,bv
      real volume
      REAL          :: s,const
      COMPLEX       :: xdnout,factor
      INTEGER       :: ii3,ix,iy,iz,nplo,nbn
      INTEGER       :: nbmin,nbmax
      INTEGER       :: nplot,nq,nt,jm,ii1,ii2
      LOGICAL       :: twodim
      real,allocatable::knorm(:,:)
      real,allocatable::wfnorm(:)
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      INTEGER :: grid(3),k
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(jmtd)
      REAL    :: point(3),post(3,natd)
      CHARACTER(len=30):: filename
      CHARACTER(len=20):: name1,name2,name3
      CHARACTER(len=10):: vandername
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename 
      integer cell1,cell2,cell3,pw1,pw2,pw3
      complex ci
      complex,allocatable::wannint(:,:,:,:)
      complex,allocatable::wannz(:,:),wannz2(:,:)
      real denom

      TYPE(t_zmat)  :: zzMat, zMat
      TYPE(t_usdus) :: usdus

      ci=cmplx(0.,1.)
c specify number of unit-cells that are calculated
      cell1=3
      cell2=3
      cell3=3

      um_format=.true.
      l_byindex=.false.
      l_byenergy=.false.
      l_bynumber=.false.
      if(sortrule==1)l_byindex=.true.
      if(sortrule==2)l_bynumber=.true.
      if(sortrule==3)l_byenergy=.true.


      lmd = lmaxd*(lmaxd+2)
      nkpts = maxval(nkpt(:))

      nrec = 0
      nlotot = 0
      do n = 1, ntype
        do l = 1,nlo(n)
          nlotot = nlotot + neq(n) * ( 2*llo(l,n) + 1 )
        enddo
      enddo


cccccccccccccccc   initialize the potential   cccccccccccc

      allocate ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
      allocate ( vz(nmzd,2,jspd),vrf(jmtd,0:nlhd,ntypd,jspd) )
      allocate ( vr(jmtd,ntypd,jspd) )

      open (8,file='pottot',form='unformatted',status='old')
      rewind (8)

      call loddop(stars,vacuum,atoms,sphhar,input,sym,
     >            8,
     <            iter,vrf,vpw,vz,vzxy)

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
      wannierspin=jspd
      if(l_soc) wannierspin=2

      allocate ( kveclo(nlotot),nv(jspd) )
      allocate ( k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd) )
      allocate ( ff(ntypd,jmtd,2,0:lmaxd) )
      allocate ( gg(ntypd,jmtd,2,0:lmaxd) )
      allocate ( usdus%us(0:lmaxd,ntypd,jspins) )
      allocate ( usdus%uds(0:lmaxd,ntypd,jspins) )
      allocate ( usdus%dus(0:lmaxd,ntypd,jspins) )
      allocate ( usdus%duds(0:lmaxd,ntypd,jspins) )
      allocate ( usdus%ddn(0:lmaxd,ntypd,jspins) )
      allocate ( usdus%ulos(nlod,ntypd,jspins) )
      allocate ( usdus%dulos(nlod,ntypd,jspins) )
      allocate ( usdus%uulon(nlod,ntypd,jspins) )
      allocate ( usdus%dulon(nlod,ntypd,jspins) )
      allocate ( usdus%uloulopn(nlod,nlod,ntypd,jspins) )
c****************************************************
c cycle by spins starts! 
c****************************************************
      do 110 jspin=1,wannierspin   ! cycle by spins

       jspin2=jspin
       if(l_soc .and. jspins.eq.1)jspin2=1
       jsp_start = jspin ; jsp_end = jspin

      jsp_start = jspin ; jsp_end = jspin

c*******************************************************
c      get num_bands and num_wann from WF1.amn (WF2.amn)
c*******************************************************
      l_file=.false.
      inquire(file=spin12(jspin)//'.amn',exist=l_file)
      open(355,file=spin12(jspin)//'.amn')
      read(355,*)
      read(355,*)num_bands,kpun,num_wann
      close(355)
      if(l_byindex.and.
     &     .not.((1+band_max(jspin)-band_min(jspin)).eq.num_bands))
     &     CALL juDFT_error("1+band_max-band_min/=num_bands",calledby
     +     ="wannier_to_lapw")

c**************************************************************
!   for bzsym = .true.: determine mapping between kpts and w90kpts
c**************************************************************
      if (l_bzsym) then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("w90kpts not found, needed if bzsym",calledby
     +        ="wannier_to_lapw")
         open(412,file='w90kpts',form='formatted')
         read(412,*)fullnkpts
         close(412)
         print*,"fullnkpts=",fullnkpts
         IF(fullnkpts<=nkpts) CALL juDFT_error("fullnkpts.le.nkpts"
     +        ,calledby ="wannier_to_lapw")
         allocate(irreduc(fullnkpts),mapkoper(fullnkpts))
         l_file=.false.
         inquire(file='kptsmap',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error
     +        ("kptsmap not found, needed if bzsym",calledby
     +        ="wannier_to_lapw")
         open(713,file='kptsmap')
         do i=1,fullnkpts
            read(713,*)kpt,irreduc(i),mapkoper(i)
            IF(kpt/=i) CALL juDFT_error("kpt.ne.i",calledby
     +           ="wannier_to_lapw")
            print*,i,irreduc(i),mapkoper(i)
         enddo   
         close(713)
         IF(MAXVAL(irreduc(:))/=nkpts) CALL juDFT_error
     +        ("max(irreduc(:))/=nkpts",calledby ="wannier_to_lapw")
      else   
         fullnkpts=nkpts

      endif   

      IF(kpun/=fullnkpts) CALL juDFT_error
     +     ("mismatch in kpun and fullnkpts",calledby ="wannier_to_lapw"
     +     )

c**************************************************************
c          read in u_matrix
c**************************************************************

      if(.not.l_proj_plot)then
         allocate(lwindow(num_bands,fullnkpts))
         allocate(ndimwin(fullnkpts))
         allocate(u_matrix(num_bands,num_wann,fullnkpts))
         call wann_read_umatrix(fullnkpts,num_wann,num_bands,
     >              um_format,jspin,wan90version,
     <              have_disentangled,lwindow,ndimwin,
     <              u_matrix)
         if(.not.have_disentangled)
     &       deallocate(lwindow,ndimwin)
         if(have_disentangled)allocate(inc_band(num_bands))

      else
c**************************************************************
c             read WF1.umn (WF2.umn) (if projmethod)
c**************************************************************
         l_file=.false.
         inquire(file=spin12(jspin)//'.umn',exist=l_file)
         IF(.NOT.l_file)  CALL juDFT_error("no umn file found" ,calledby
     +        ="wannier_to_lapw")
         open(419,file=spin12(jspin)//'.umn')
         read(419,*)     !num_wann,num_bands
         allocate(u_matrix(num_bands,num_wann,fullnkpts))
         do ikpt=1,fullnkpts
            do j=1,num_wann
               do i=1,num_bands
                  read(419,*)u_matrix(i,j,ikpt)
               enddo
            enddo   
         enddo   
         close(419)
      endif   
      
***********************************************************
***********************************************************

      print*,"num_wann=",num_wann
      print*,"num_bands=",num_bands

cccccccccccc   read in the eigenvalues and vectors   cccccc

      l_p0 = .false.
      if (irank.eq.0) l_p0 = .true.
#ifdef CPP_NEVER
      call cdn_read0(
     >               lmaxd,ntypd,nlod,neigd,wannierspin,
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
#else
      call judft_error("NOT implemented")
#endif


      allocate(knorm(fullnkpts,num_bands))
      allocate(wann_acof(num_wann,0:lmd,natd,-cell1:cell1,-cell2:cell2,
     & -cell3:cell3))
      allocate(wann_bcof(num_wann,0:lmd,natd,-cell1:cell1,-cell2:cell2,
     & -cell3:cell3))
      allocate(wann_ccof(num_wann,-llod:llod,nlod,natd,-cell1:cell1,
     & -cell2:cell2,-cell3:cell3))
      allocate(wannint(-unigrid(4):unigrid(4),-unigrid(4):unigrid(4),
     , -unigrid(4):unigrid(4),num_wann))
      wann_acof(:,:,:,:,:,:)=0.0
      wann_bcof(:,:,:,:,:,:)=0.0
      wann_ccof(:,:,:,:,:,:,:)=0.0
      wannint(:,:,:,:)=0.0

      print*,"num_bands=",num_bands
      print*,"num_wann=",num_wann
      knorm(:,:)=0.0

c******************************************************************
c          beginning of k-point loop,each may be a separate task
c******************************************************************
      i_rec = 0 ; n_rank = 0
      do ikpt = 1,fullnkpts  ! loop by k-points starts

        i_rec = i_rec + 1
        if (mod(i_rec-1,isize).eq.irank) then
        print*,"k-point=",ikpt
        kptibz=ikpt 
        if(l_bzsym) kptibz=irreduc(ikpt)
        if(l_bzsym) oper=mapkoper(ikpt)

       if(have_disentangled) then
          inc_band(:)=lwindow(:,ikpt)
          num_inc=ndimwin(ikpt)
       end if

      allocate (we(neigd),eigg(neigd))



      call wann_read_eig(
     >              eig_id,
     >              lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >              irank,isize,kptibz,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,
     <              nmat,nv,ello,evdu,epar,kveclo,
     <              k1,k2,k3,bkpt,wk,nbands,eigg,zzMat,
     >              .false.,1)


      zMat%l_real = zzMat%l_real
      zMat%nbasfcn = zzMat%nbasfcn
      zMat%nbands = zzMat%nbands
      IF (zzMat%l_real) THEN
         ALLOCATE (zMat%z_r(zMat%nbasfcn,zMat%nbands))
         zMat%z_r = 0.0
      ELSE
         ALLOCATE (zMat%z_c(zMat%nbasfcn,zMat%nbands))
         zMat%z_c = CMPLX(0.0,0.0)
      END IF

      nslibd = 0

c...we work only within the energy window

      eig(:) = 0.

      print*,"bands used"
      do i = 1,nbands
       if ((eigg(i).ge.e1s .and. nslibd.lt.num_bands.and.l_bynumber)
     &.or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.l_byenergy)
     &.or.(i.ge.band_min(jspin)
     & .and.i.le.band_max(jspin).and.l_byindex))then
          print*,i
        nslibd = nslibd + 1
        eig(nslibd) = eigg(i)
        we(nslibd) = we(i)
        IF(zzMat%l_real) THEN
          do j = 1,nv(jspin) + nlotot
             zMat%z_r(j,nslibd) = zzMat%z_r(j,i)
          end do
        ELSE
          do j = 1,nv(jspin) + nlotot
             zMat%z_c(j,nslibd) = zzMat%z_c(j,i)
          end do
        END IF
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
c          calculate a-, b-, and c-coefficients
c******************************************************************

      noccbd = nslibd

      allocate ( acof(noccbd,0:lmd,natd),
     &           bcof(noccbd,0:lmd,natd),
     &           ccof(-llod:llod,noccbd,nlod,natd))

      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      ccof(:,:,:,:) = cmplx(0.,0.)

c...generation of the A,B,C coefficients in the spheres 
c...for the lapws and local orbitals, summed by the basis functions 

      CALL abcof(input,atoms,noccbd,sym,cell,bkpt,lapw,noccbd,usdus,
     +           noco,jspin,kveclo,oneD,acof,bcof,ccof,zMat)

      call abcrot(
     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >        noccbd,lmax,nlo,llo,nop,ngopr,mrot,invsat,invsatnr,
     >                 bmat,odi,ods,
     X                 acof,bcof,ccof)

c***************************************************************
c        calculate wannierfunctions' a-,b-, and c-coefficients
c***************************************************************
      do i1=-cell1,cell1
       do i2=-cell2,cell2
        do i3=-cell3,cell3
         factor=exp(ci*tpi_const*(i1*bkpt(1)+i2*bkpt(2)+i3*bkpt(3)))  
         do i=1,num_wann
          do n=1,noccbd
            wann_acof(i,0:lmd,1:natd,i1,i2,i3)=
     =          wann_acof(i,0:lmd,1:natd,i1,i2,i3)+
     +          u_matrix(n,i,ikpt)*acof(n,0:lmd,1:natd)*factor 
            
            wann_bcof(i,0:lmd,1:natd,i1,i2,i3)=
     =          wann_bcof(i,0:lmd,1:natd,i1,i2,i3)+
     +          u_matrix(n,i,ikpt)*bcof(n,0:lmd,1:natd)*factor     

            wann_ccof(i,-llod:llod,1:nlod,1:natd,i1,i2,i3)=
     =          wann_ccof(i,-llod:llod,1:nlod,1:natd,i1,i2,i3)+
     +    u_matrix(n,i,ikpt)*ccof(-llod:llod,n,1:nlod,1:natd)*factor
          enddo
         enddo
        enddo
       enddo   
      enddo
c***************************************************************
c       calculate wannierfunctions' planewave-expansion
c***************************************************************   
      allocate(wannz(nv(jspin),num_wann))
      allocate(wannz2(nv(jspin),num_wann))
      wannz(:,:)=cmplx(0.0,0.0)
      do n=1,noccbd
       do m=1,num_wann
        IF(zMat%l_real) THEN
          do j=1, nv(jspin)
            wannz(j,m)=wannz(j,m)+zMat%z_r(j,n)*
     +                            u_matrix(n,m,ikpt)/sqrt(omtil)
          enddo  
        ELSE
          do j=1, nv(jspin)
            wannz(j,m)=wannz(j,m)+zMat%z_c(j,n)*
     +                            u_matrix(n,m,ikpt)/sqrt(omtil)
          enddo  
        END IF

       enddo   
      enddo
      print*,"unigrid=",unigrid(:)
      do j=1,nv(jspin)  
       do pw1=-unigrid(4),unigrid(4)
        do pw2=-unigrid(4),unigrid(4)
         do pw3=-unigrid(4),unigrid(4)
          denom=-pw1+unigrid(1)*(k1(j,jspin)+bkpt(1))  
          wannz2(j,:)=wannz(j,:)
          if(abs(denom).gt.1e-5)then
             denom=denom*tpi_const/2
             factor=cmplx(cos(denom),sin(denom))
             wannz2(j,:)=wannz2(j,:)*(factor-conjg(factor))/
     /              (ci*denom*2)
          endif   
          denom=-pw2+unigrid(2)*(k2(j,jspin)+bkpt(2))  
          if(abs(denom).gt.1e-5)then
             denom=denom*tpi_const/2
             factor=cmplx(cos(denom),sin(denom))
             wannz2(j,:)=wannz2(j,:)*(factor-conjg(factor))/
     /              (ci*denom*2)
          endif   
          denom=-pw3+unigrid(3)*(k3(j,jspin)+bkpt(3))  
          if(abs(denom).gt.1e-5)then
             denom=denom*tpi_const/2
             factor=cmplx(cos(denom),sin(denom))
             wannz2(j,:)=wannz2(j,:)*(factor-conjg(factor))/
     /              (ci*denom*2)
          endif   
             wannint(pw1,pw2,pw3,:)=wannint(pw1,pw2,pw3,:)+
     +                wannz2(j,:)
         enddo    
        enddo
       enddo
      enddo 
      deallocate(wannz,wannz2)

      deallocate ( acof,bcof,ccof,we,eigg )

      write (*,*) 'nslibd=',nslibd

      endif!processors

      enddo !loop over k-points
c****************************************************************
c        end of k-loop
c****************************************************************

c****************************************************************
c     write radial wavefunctions to file
c****************************************************************
      open(344,file=spin12(jspin)//'.lapw',form='unformatted')
      write(344)num_wann,cell1,cell2,cell3
      write(344)amat(1:3,1:3)
      write(344)ntype,jmtd,lmaxd,nlod,natd,lmd,llod
      write(344)pos(1:3,1:natd)
      write(344)neq(1:ntype)
      write(344)zatom(1:ntype)
      write(344)dx(1:ntype)
      write(344)rmt(1:ntype)
      write(344)jri(1:ntype)
      write(344)rmsh(1:jmtd,1:ntype)
      write(344)ff(1:ntype,1:jmtd,1:2,0:lmaxd)
      write(344)gg(1:ntype,1:jmtd,1:2,0:lmaxd)
      write(344)flo(1:ntype,1:jmtd,1:2,1:nlod)
c****************************************************************
c      write a-,b-, and c-coefficients to file
c**************************************************************** 

      write(344)wann_acof(1:num_wann,0:lmd,1:natd,-cell1:cell1,
     &                    -cell2:cell2,-cell3:cell3)
      write(344)wann_bcof(1:num_wann,0:lmd,1:natd,-cell1:cell1,
     &                    -cell2:cell2,-cell3:cell3)
      write(344)wann_ccof(1:num_wann,-llod:llod,1:nlod,1:natd,
     &                -cell1:cell1,-cell2:cell2,-cell3:cell3)
c****************************************************************
c        write planewave expansion to file
c****************************************************************
      write(344)unigrid(1:4)
      write(344)wannint(:,:,:,:)
      close(344)






      deallocate(flo)


      if(have_disentangled)deallocate(lwindow,ndimwin,inc_band)


      deallocate(u_matrix)

#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,mpiierr)
#endif      

      nrec=nrec+nkpts
110   continue ! end of cycle by spins


      deallocate ( vr,vz,kveclo,nv,k1,k2,k3 )
      deallocate ( ff,gg)
 
#ifdef CPP_MPI
      call MPI_BARRIER(mpi_comm,mpiierr)
#endif      

      end subroutine wannier_to_lapw
      end module m_wannier_to_lapw
