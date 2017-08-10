      MODULE m_wann_orbmag
      USE m_fleurenv
      CONTAINS
      SUBROUTINE wann_orbmag(
     >      l_dulo,l_noco,l_ss,lmaxd,ntypd,
     >      neigd,natd,nop,nvd,jspd,nbasfcn,llod,nlod,ntype,
     >      nwdd,omtil,nlo,llo,lapw_l,invtab,mrot,ngopr,neq,lmax,
     >      invsat,invsatnr,nkpt,taual,rmt,amat,bmat,bbmat,alph,
     >      beta,qss,sk2,phi2,odi,ods,irank,isize,n3d,nmzxyd,nmzd,
     >      jmtd,nlhd,nq3,nvac,invs,invs2,film,nlh,jri,ntypsd,
     >      ntypsy,jspins,nkptd,dx,n2d,rmsh,e1s,e2s,ulo_der,
     >      ustep,ig,k1d,k2d,k3d,rgphs,slice,kk,nnne,
     >      z1,nv2d,nmzxy,nmz,delz,ig2,area,tau,zatom,nq2,nop2,
     >      volint,symor,pos,ef,irecl,l_soc,
     >      memd,lnonsph,clnu,lmplmd,mlh,nmem,llh,lo1l,
     >      theta,phi,soc_opt,
     >      ikpt,nntot,fullnkpts,bpt,wannierspin,jspin,jspin_b,
     >      nlotot,nrec,nrec_b,
     >      l_gwf,iqpt,numbands,wann,lmd,ngopr1,
     >      us,uds,dus,duds,ddn,ulos,dulos,uulon,dulon,
     >      qpt_i,alph_i,beta_i,k1,k2,k3,nv,nmat,noccbd,nslibd,
     >      kveclo,bkpt,chi,
     >      flo,ff,gg,
     >      gb,kdiff,nbnd,uHu)

      use m_types
      use m_wann_rw_eig
      use m_abcof_small
      use m_radfun
      use m_radflo
      use m_cdnread, only : cdn_read0, cdn_read
      use m_od_types, only : od_inp, od_sym
      use m_wann_abinv
      use m_wann_kptsrotate
      use m_matmul,only : matmul3,matmul3r
      use m_wann_uHu_int2
      use m_wann_mmkb_sph
      use m_wann_ujugaunt

      IMPLICIT NONE
#include "cpp_double.h"
#ifdef CPP_MPI
      include 'mpif.h'
      integer ierr(3)
      integer cpu_index
      integer stt(MPI_STATUS_SIZE)

#endif
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
      real,    intent (in) :: alph(ntypd),beta(ntypd),qss(3)
      real,    intent (in) :: pos(3,natd)
      real, intent(in)     :: ef

      complex, intent (in) :: ustep(n3d)
      logical, intent (in) :: l_dulo(nlod,ntypd),l_noco,l_ss
      integer, intent (in) :: irecl
      logical, intent (in) :: l_soc
      integer, intent (in) :: memd
      integer, intent (in) :: lnonsph(ntypd)
      complex, intent (in) :: clnu(memd,0:nlhd,ntypsd)
      integer, intent (in) :: lmplmd
      integer, intent (in) :: mlh(memd,0:nlhd,ntypsd)
      integer, intent (in) :: nmem(0:nlhd,ntypsd)
      integer, intent (in) :: llh(0:nlhd,ntypsd)
      integer, intent (in) :: lo1l(0:llod,ntypd)
      real,    intent (in) :: sk2(n2d),phi2(n2d)
      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods
      real, intent(in) :: theta,phi
      logical,intent(in) :: soc_opt(ntype+2)

      integer, intent(in) :: ikpt,nntot,fullnkpts,lmd
      integer, intent(in) :: bpt(nntot,fullnkpts)
      integer, intent(in) :: ngopr1(natd)
      integer, intent(in) :: wannierspin,jspin,jspin_b,nbnd,nslibd
      integer, intent(in) :: nv(wannierspin)
      integer, intent(in) :: nlotot,nrec,nrec_b,iqpt,numbands
      integer, intent(in) :: nmat,noccbd,kveclo(nlotot)
      real, intent(in) :: bkpt(3)

      integer,intent(in) :: k1(nvd,wannierspin)
      integer,intent(in) :: k2(nvd,wannierspin)
      integer,intent(in) :: k3(nvd,wannierspin)
      logical, intent(in) :: l_gwf
      type(t_wann), intent(in) :: wann
      real, intent(in) :: us(0:lmaxd,ntypd,2)
      real, intent(in) :: uds(0:lmaxd,ntypd,2)
      real, intent(in) :: dus(0:lmaxd,ntypd,2)
      real, intent(in) :: duds(0:lmaxd,ntypd,2)
      real, intent(in) :: ddn(0:lmaxd,ntypd,2)
      real, intent(in) :: ulos(nlod,ntypd,2)
      real, intent(in) :: uulon(nlod,ntypd,2)
      real, intent(in) :: dulon(nlod,ntypd,2)
      real, intent(in) :: dulos(nlod,ntypd,2)
      real, intent(in) :: qpt_i(3),alph_i(ntypd),beta_i(ntypd)

      real, intent(in) :: ff(ntypd,jmtd,2,0:lmaxd,2)
      real, intent(in) :: gg(ntypd,jmtd,2,0:lmaxd,2)
      real, intent(in) :: flo(ntypd,jmtd,2,nlod,2)

      complex, intent(in) :: chi(ntypd)
      complex  :: ujug(0:lmd,0:lmd,ntype,nntot,4),
     >           ujdg(0:lmd,0:lmd,ntype,nntot,4), 
     >           djug(0:lmd,0:lmd,ntype,nntot,4),
     >           djdg(0:lmd,0:lmd,ntype,nntot,4), 
     >           ujulog(0:lmd,nlod,-llod:llod,ntype,nntot,4), 
     >           djulog(0:lmd,nlod,-llod:llod,ntype,nntot,4), 
     >           ulojug(0:lmd,nlod,-llod:llod,ntype,nntot,4), 
     >           ulojdg(0:lmd,nlod,-llod:llod,ntype,nntot,4),
     >           ulojulog(nlod,-llod:llod,nlod,-llod:llod,ntype,nntot,4) 
      integer, intent (in)  :: gb(3,nntot,fullnkpts)
      real, intent (in)     :: kdiff(3,nntot)


      complex, intent(inout) :: uHu(nbnd,nbnd,nntot,nntot)

      integer :: i,j,g1,g2,g3,g4,n
      integer :: funbas
      integer :: info
      integer, allocatable :: ipiv(:)
      complex, allocatable :: work(:)
      complex, allocatable :: tmp_mat1(:,:),tmp_mat2(:,:)

      integer :: ikpt_b,n_start,n_end,kptibz_b
      integer :: jspin2,jspin3,jspin4,jspin5
      integer :: addnoco_kb,addnoco_kb2,addnoco,addnoco2
      integer :: nmat_b,nbands_b,nslibd_b,noccbd_b
      integer :: nv_b(wannierspin)
      integer :: kveclo_b(nlotot)
      integer :: k1_b(nvd,wannierspin)
      integer :: k2_b(nvd,wannierspin)
      integer :: k3_b(nvd,wannierspin)
      real :: eig_b(neigd),we_b(neigd)

      integer :: ikpt_b2,kptibz_b2
      integer :: nmat_b2,nbands_b2,nslibd_b2,noccbd_b2
      integer :: nv_b2(wannierspin)
      integer :: kveclo_b2(nlotot)
      integer :: k1_b2(nvd,wannierspin)
      integer :: k2_b2(nvd,wannierspin)
      integer :: k3_b2(nvd,wannierspin)
      real :: eig_b2(neigd),we_b2(neigd)

      complex, allocatable :: acof(:,:,:,:)
      complex, allocatable :: bcof(:,:,:,:)
      complex, allocatable :: ccof(:,:,:,:,:)

      complex, allocatable :: acof_b(:,:,:)
      complex, allocatable :: bcof_b(:,:,:)
      complex, allocatable :: ccof_b(:,:,:,:)

      complex, allocatable :: acof_b2(:,:,:)
      complex, allocatable :: bcof_b2(:,:,:)
      complex, allocatable :: ccof_b2(:,:,:,:)

      complex, allocatable :: acof_tmp(:,:,:)
      complex, allocatable :: bcof_tmp(:,:,:)
      complex, allocatable :: ccof_tmp(:,:,:,:)

      complex :: I1(nvd,nvd),II1(neigd,nvd),III1(neigd,nvd)
      complex :: I2(nvd,nvd),II2(nvd,neigd)

      real :: cp_time(9)
      real :: bkpt_b(3),wk_b
      real :: bkpt_b2(3),wk_b2
      real :: ello(nlod,ntypd,max(2,jspd)),evac(2,max(2,jspd))
      real :: epar(0:lmaxd,ntypd,max(2,jspd)),evdu(2,max(jspd,2))
      real :: eigg(neigd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex :: zz(nbasfcn,neigd)
      complex :: z_b(nbasfcn,neigd),z_tmp(nbasfcn,neigd)
      complex :: z_b2(nbasfcn,neigd)
#else
      real :: zz(nbasfcn,neigd)
      real :: z_b(nbasfcn,neigd),ztemp(nbasfcn,neigd)
      real :: z_b2(nbasfcn,neigd)
#endif

#ifdef CPP_INVERSION
      real,allocatable :: hpack(:), ham(:,:)  
      real,allocatable :: Spack(:), Sinv(:,:)
#else
      complex,allocatable :: hpack(:), ham(:,:)
      complex,allocatable :: Spack(:), Sinv(:,:)
      complex,allocatable :: S(:,:)
#endif
      integer :: matsize,nv1,nv2,nbasfcn2,lwork
      integer :: ind_uu(2),ind_dd(2),ind_du(2)
      real :: diff,mdiff
      character(len=20)::fn

      info=0
      addnoco =0
      addnoco2=0

      addnoco_kb =0
      addnoco_kb2=0

      nbasfcn2 = nbasfcn
      nbasfcn2 = nv(1)+nv(2)+nlotot
!**************************************!
!      read in Hamiltonian at k        !
!**************************************!
      jspin2 = jspin
      if(l_noco) jspin2 = 1
      write(*,'(a,f8.5,1x,f8.5,1x,f8.5)')'k   : ',bkpt
      write(fn,'("ham_",i4.4,".",i1)')ikpt,jspin2
      open(879,file=trim(fn),form='unformatted')
      read(879)matsize,nv1,nv2
      if(nv1.ne.nv(1)) stop 'nv1.ne.nv(1)'
      if(nv2.ne.nv(2)) stop 'nv2.ne.nv(2)'
      allocate(hpack(matsize))
      allocate(Spack(matsize))
      read(879)hpack(1:matsize)
      read(879)Spack(1:matsize)
      close(879)
      write(*,*)'matsize',matsize

      ! indices for checks of up-up, down-down, and down-up blocks
      !ind_uu(1)=1
      !ind_uu(2)=nv1*(nv1+1)/2
      !ind_dd(1)=ind_uu(2)+nv(1)+1
      !ind_dd(2)=(nv1+nv2)*(nv1+nv2+1)/2
      !do i=ind_uu(1),ind_uu(2)
      !   write(*,*)'up-up',i,Spack(i)
      !enddo
      !do i=ind_dd(1),ind_dd(2)
      !   write(*,*)'down-down',i,Spack(i)
      !enddo

      ! are overlaps real?
      do i=1,matsize
            if(aimag(Spack(i)).ne.0.0) write(*,*)'imag',i
      enddo

      ! set up S for later check 1=SS^{-1}
      allocate(S(nbasfcn2,nbasfcn2))
      do j=1,nbasfcn2
         do i=1,j
            n = i+(j-1)*j/2
            S(i,j)= Spack(n)
         enddo
      enddo
      do j=1,nbasfcn2
         do i=j+1,nbasfcn2
            S(i,j) = conjg(S(j,i))
         enddo
      enddo

      ! write out diagonal of S
      !do i=1,nbasfcn2
      !   write(*,*)'S(i,i)',i,S(i,i)
      !enddo

      ! S^{-1} for indefinite Hermitian S (packed)
      !allocate(ipiv(nbasfcn2), work(nbasfcn2))
      !call zhptrf('U',nbasfcn2,Spack,ipiv,info)
      !if(info.ne.0) then stop 'info1'
      !call zhptri('U',nbasfcn2,Spack,ipiv,work,info)
      !if(info.ne.0) then stop 'info2'
      !deallocate(ipiv,work)

      ! set up H and S^{-1} resp. S matrices
      allocate(ham(nbasfcn2,nbasfcn2))
      allocate(Sinv(nbasfcn2,nbasfcn2))
      do j=1,nbasfcn2   ! FOR 'U'
         do i=1,j
            n = i+(j-1)*j/2
            ham(i,j) = hpack(n)
            Sinv(i,j)= Spack(n)
         enddo
      enddo
      !do i=1,nbasfcn2   ! FOR 'L'
      !   do j=1,i
      !      n = i+(j-1)*(2*nbasfcn2-j)/2
      !      ham(i,j) = hpack(n)
      !      Sinv(i,j)= Spack(n)
      !   enddo
      !enddo
      deallocate(hpack)
      deallocate(Spack)

      ! complete Hermitian S^{-1}
      do j=1,nbasfcn2       ! FOR 'U'
         do i=j+1,nbasfcn2
            ham(i,j) = conjg(ham(j,i))
            Sinv(i,j) = conjg(Sinv(j,i))
         enddo
      enddo
      !do i=1,nbasfcn2       ! FOR 'L'
      !   do j=i+1,nbasfcn2
      !      ham(i,j) = conjg(ham(j,i))
      !      Sinv(i,j) = conjg(Sinv(j,i))
      !   enddo
      !enddo

      ! is S^{-1} real?
      do i=1,nbasfcn2
         do j=1,i
            if(aimag(Sinv(i,j)).ne.0.0) write(*,*)'imag',i,j
         enddo
      enddo

      ! S^{-1} from positive definite Hermitian S
      call zpotrf('U',nbasfcn2,Sinv,nbasfcn2,info)
      if(info.ne.0) stop 'info1'
      call zpotri('U',nbasfcn2,Sinv,nbasfcn2,info)
      if(info.ne.0) stop 'info2'

      ! S^{-1} from indefinite Hermitian S
      !lwork = nbasfcn2
      !allocate(ipiv(nbasfcn2),work(lwork))
      !call zhetrf('U',nbasfcn2,Sinv,nbasfcn2,
     >!            ipiv,work,lwork,info)
      !if(info.ne.0) stop 'info1'
      !call zhetri('U',nbasfcn2,Sinv,nbasfcn2,
     >!            ipiv,work,info)
      !call zgetrf(nbasfcn2,nbasfcn2,Sinv,nbasfcn2,ipiv,info)
      !if(info.ne.0) stop 'info1'
      !call zgetri(nbasfcn2,Sinv,nbasfcn2,ipiv,work,lwork,info)
      !deallocate(ipiv,work)
      !if(info.ne.0) stop 'info2'

      do j=1,nbasfcn2
         do i=j+1,nbasfcn2
            Sinv(i,j) = conjg(Sinv(j,i))
         enddo
      enddo

      ! write out some elements of S^{-1} and H
      do i=1,2
         do j=1,2
            write(*,'(a,i2,i2,1x,2f24.18)')'Sinv',i,j,Sinv(i,j)
            write(*,'(a,i2,i2,1x,2f24.18)')'H   ',i,j,ham(i,j)
         enddo
      enddo
      write(*,*)

      ! check whether spin-offdiagonal blocks are zero for overlaps
      do j=1+nv(1)+nlotot,nbasfcn2
         do i=1,nv(1)
            if(S(i,j).ne.cmplx(0.,0.)) write(*,*)'S up-down',i,j
            if(S(i,j).ne.cmplx(0.,0.)) write(*,*)'S up-down',i,j
            if(Sinv(i,j).ne.cmplx(0.,0.)) write(*,*)'Sinv up-down',i,j
            if(Sinv(j,i).ne.cmplx(0.,0.)) write(*,*)'Sinv down-up',i,j
         enddo
      enddo

      ! check whether resulting matrices are Hermitian
      do i=1,nbasfcn2
         do j=1,i!nbasfcn2
            if(abs(ham(j,i)-conjg(ham(i,j))).gt.1e-14) then
               write(*,*)'H not Hermitian',j,i,ikpt,jspin2
            endif
            if(abs(Sinv(j,i)-conjg(Sinv(i,j))).gt.1e-14) then
               write(*,*)'inv. S not Hermitian',j,i,ikpt,jspin2
            endif
         enddo
      enddo

      allocate(tmp_mat1(nbasfcn2,nbasfcn2))
      ! compute S^{-1} S = 1
      call zgemm('N','N',nbasfcn2,nbasfcn2,nbasfcn2,
     >           cmplx(1.0,0.0),Sinv,nbasfcn2,
     >           S,nbasfcn2,cmplx(0.0,0.0),tmp_mat1,nbasfcn2)
      do i=1,nbasfcn2
         do j=1,nbasfcn2
            if(i.eq.j) then
               if(abs(tmp_mat1(i,j)-1.0).gt.1e-10)
     >            write(*,*)'diag',i,tmp_mat1(i,j)
            else
               if(abs(tmp_mat1(i,j)).gt.1e-10)
     >            write(*,*)'offd',i,j,tmp_mat1(i,j)
            endif
         enddo
      enddo

      ! replace H with S^{-1} H S^{-1}
      call zgemm('N','N',nbasfcn2,nbasfcn2,nbasfcn2,
     >           cmplx(1.0,0.0),Sinv,nbasfcn2,
     >           ham,nbasfcn2,cmplx(0.0,0.0),tmp_mat1,nbasfcn2)
      call zgemm('N','N',nbasfcn2,nbasfcn2,nbasfcn2,
     >           cmplx(1.0,0.0),tmp_mat1,nbasfcn2,
     >           Sinv,nbasfcn2,cmplx(0.0,0.0),ham,nbasfcn2)
      deallocate(tmp_mat1)
      deallocate(Sinv)
      deallocate(S)

      ! write out elements of S^{-1} H S^{-1}
      do i=1,2
         do j=1,2
            !write(*,'(i2,i2,2(f24.18))')j,i,real(ham(j,i)),
     >      !                               aimag(ham(j,i))
            write(*,'(i2,i2,1x,2(f24.12))')j,i,real(ham(j,i)),
     >                                         aimag(ham(j,i))
         enddo
      enddo

      !mdiff = 0.0
      !do i=1,nbasfcn2
      !   do j=1,i
      !      diff = abs(ham(j,i)-conjg(ham(i,j)))
      !      !if(diff.gt.1e-8) write(*,*)j,i,diff
      !      if(diff.gt.mdiff) mdiff=diff
      !      !if(abs(ham(j,i)-conjg(ham(i,j))).gt.1e-14) then
      !      !   write(*,*)'S^-1 H S^-1 not Hermitian',j,i,ikpt,jspin2
      !      !endif
      !   enddo
      !enddo
      !write(*,*)'max. diff',mdiff


      ! Htr -> eV
      ham = 27.2 * ham  !not yet tested


      allocate(tmp_mat2(nbnd,nbnd))
!**************************************!
!      a,b,c coefficients at k         !
!**************************************!
      allocate ( acof(nvd,0:lmd,natd,2),
     &           bcof(nvd,0:lmd,natd,2),
     &           ccof(-llod:llod,nbasfcn,nlod,natd,2))
      acof = cmplx(0.,0.)
      bcof = cmplx(0.,0.)
      ccof = cmplx(0.,0.)
      allocate ( acof_tmp(noccbd,0:lmd,natd),
     &           bcof_tmp(noccbd,0:lmd,natd),
     &           ccof_tmp(-llod:llod,noccbd,nlod,natd))
      do jspin3 = 1,2
         addnoco=0
         if(l_noco.and.jspin3.eq.2) addnoco = nv(1)+nlotot
      do i=1,nv(jspin3)
         z_tmp = 0.0
         z_tmp(i+addnoco,1)=1.0
         call abcof_small(
     >           lmaxd,ntypd,neigd,noccbd,natd,nop,nvd,wannierspin,
     >           lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >           ntype,mrot,ngopr1,taual,neq,lmax,rmt,omtil,
     >           bmat,bbmat,bkpt,k1,k2,k3,nv,nmat,noccbd,z_tmp,
     >           us(0,1,jspin3),dus(0,1,jspin3),uds(0,1,jspin3),
     >           duds(0,1,jspin3),ddn(0,1,jspin3),invsat,invsatnr,
     >           ulos(1,1,jspin3),uulon(1,1,jspin3),dulon(1,1,jspin3),
     >           dulos(1,1,jspin3),llo,nlo,l_dulo,lapw_l,
     >           l_noco,l_ss,jspin3,alph_i,beta_i,qpt_i,
     >           kveclo,odi,ods,
     <           acof_tmp(1,0,1),bcof_tmp(1,0,1),
     <           ccof_tmp(-llod,1,1,1),i)

         call wann_abinv(
     >        ntypd,natd,noccbd,lmaxd,lmd,llod,nlod,ntype,neq,
     >        noccbd,lmax,nlo,llo,invsat,invsatnr,bkpt,taual,
     X        acof_tmp,bcof_tmp,ccof_tmp)
         acof(i,:,:,jspin3) = acof_tmp(1,:,:)
         bcof(i,:,:,jspin3) = bcof_tmp(1,:,:)
         ccof(:,i,:,:,jspin3)=ccof_tmp(:,1,:,:)
      enddo
      enddo!jspin3
      deallocate(acof_tmp,bcof_tmp,ccof_tmp)

      ! compute necessary ujugaunt(spin3,spin4)
      ! TODO: ujugaunt are k-point independent->do computation before k-loop
      do jspin4=1,2
         do jspin3=1,2
            jspin5=jspin3+2*(jspin4-1)
         call wann_ujugaunt(
     >            llod,nntot,kdiff,lmax,ntype,ntypd,bbmat,
     >            bmat,nlod,nlo,llo,flo(:,:,:,:,jspin3),
     >            flo(:,:,:,:,jspin4),
     >            ff(:,:,:,:,jspin3),
     >            ff(:,:,:,:,jspin4),
     >            gg(:,:,:,:,jspin3),
     >            gg(:,:,:,:,jspin4),jri,rmsh,dx,jmtd,
     >            lmaxd,lmd,
     <            ujug(:,:,:,:,jspin5),
     >            ujdg(:,:,:,:,jspin5),
     >            djug(:,:,:,:,jspin5),
     >            djdg(:,:,:,:,jspin5),
     <            ujulog(:,:,:,:,:,jspin5),
     >            djulog(:,:,:,:,:,jspin5),
     >            ulojug(:,:,:,:,:,jspin5),
     >            ulojdg(:,:,:,:,:,jspin5),
     >            ulojulog(:,:,:,:,:,:,jspin5),.false.,1)
         enddo
      enddo


!**************************************!
!      loop over neighbors b1          !
!**************************************!
      do ikpt_b = 1,2!nntot
         kptibz_b=bpt(ikpt_b,ikpt)
         n_start=1
         n_end=neigd

   !**************************************!
   !      read eigenvectors at k+b1       !
   !**************************************!   
#if(!defined(CPP_HDF) && defined(CPP_MPI))
         call wann_read_eig(
     >              lmaxd,ntypd,nlod,neigd,nvd,wannierspin,
     >              irank,isize,kptibz_b,jspin,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec,irecl,
     <              nmat_b,nv_b,ello,evdu,epar,kveclo_b,
     <              k1_b,k2_b,k3_b,bkpt_b,wk_b,nbands_b,
     <              eigg,zz,cp_time,66,l_gwf,iqpt)
#else
         call cdn_read(
     >                lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >                irank,isize,kptibz_b,jspin,nbasfcn,nlotot,
     >                l_ss,l_noco,nrec,kptibz_b,66,
     >                neigd,n_start,n_end,
     <                nmat_b,nv_b,ello,evdu,epar,kveclo_b,
     <                k1_b,k2_b,k3_b,bkpt_b,wk_b,nbands_b,
     <                eigg,zz,cp_time)
#endif
         !write(*,'(a,f8.5,1x,f8.5,1x,f8.5)')'k+b1: ',bkpt_b

         nslibd_b = 0
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
         z_b(:,:) = cmplx(0.,0.)
#else
         z_b(:,:) = 0.
#endif
         eig_b(:) = 0.

         do i = 1,nbands_b
          if((eigg(i).ge.e1s.and.nslibd_b.lt.numbands
     &       .and.wann%l_bynumber)
     &       .or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &       .or.(i.ge.wann%band_min(jspin)
     &       .and.(i.le.wann%band_max(jspin))
     &       .and.wann%l_byindex))then
            nslibd_b = nslibd_b + 1
            eig_b(nslibd_b) = eigg(i)
            we_b(nslibd_b) = we_b(i)
            if(l_noco)then
               funbas =        nv_b(1) + nlotot
               funbas = funbas+nv_b(2) + nlotot
            else
               funbas = nv_b(jspin) + nlotot
            endif
            do  j = 1,funbas
              z_b(j,nslibd_b) = zz(j,i)
            enddo
          endif
         enddo
         noccbd_b = nslibd_b

         addnoco_kb=0
         if(l_noco.and.jspin.eq.2) addnoco_kb =nv_b(1)+nlotot
   !**************************************!
   !      a,b,c coefficients at k+b1      !
   !**************************************!
         allocate ( acof_b(nv_b(jspin),0:lmd,natd),
     &             bcof_b(nv_b(jspin),0:lmd,natd),
     &             ccof_b(-llod:llod,nv_b(jspin),nlod,natd))

         allocate ( acof_tmp(noccbd_b,0:lmd,natd),
     &           bcof_tmp(noccbd_b,0:lmd,natd),
     &           ccof_tmp(-llod:llod,noccbd_b,nlod,natd))     
 
         do i=1,nv_b(jspin)
            z_tmp = 0.0
            z_tmp(i+addnoco_kb,1)=1.0
            call abcof_small(
     >             lmaxd,ntypd,neigd,noccbd_b,natd,nop,nvd,wannierspin,
     >             lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >             ntype,mrot,ngopr1,taual,neq,lmax,rmt,omtil,
     >             bmat,bbmat,bkpt_b,
     >             k1_b,k2_b,k3_b,nv_b,nmat_b,noccbd_b,z_tmp,
     >             us(0,1,jspin),dus(0,1,jspin),uds(0,1,jspin),
     >             duds(0,1,jspin),ddn(0,1,jspin),invsat,invsatnr,
     >             ulos(1,1,jspin),uulon(1,1,jspin),dulon(1,1,jspin),
     >             dulos(1,1,jspin),llo,nlo,l_dulo,lapw_l,
     >             l_noco,l_ss,jspin,alph_i,beta_i,qpt_i,
     >             kveclo_b,odi,ods,
     <             acof_tmp(1,0,1),bcof_tmp(1,0,1),
     <             ccof_tmp(-llod,1,1,1),i)

            call wann_abinv(
     >         ntypd,natd,noccbd_b,lmaxd,lmd,llod,nlod,ntype,neq,
     >         noccbd_b,lmax,nlo,llo,invsat,invsatnr,bkpt_b,taual,
     X         acof_tmp,bcof_tmp,ccof_tmp)

            acof_b(i,:,:)  = acof_tmp(1,:,:)
            bcof_b(i,:,:)  = bcof_tmp(1,:,:)
            ccof_b(:,i,:,:)= ccof_tmp(:,1,:,:)

         enddo
         deallocate(acof_tmp,bcof_tmp,ccof_tmp)


      ! inner spin loop 1
      do jspin3=1,2
         addnoco = 0
         if(l_noco .and. jspin3.eq.2) addnoco =nv(1)+nlotot

         I1 = cmplx(0.,0.)
   !**************************************!
   !       INT contribution at k+b1       !
   !**************************************! 
         call wann_uHu_int2( nvd,k1d,k2d,k3d,
     >                   n3d,k1_b(:,jspin),k2_b(:,jspin),k3_b(:,jspin),
     >                   nv_b(jspin),
     >                   k1(:,jspin3),k2(:,jspin3),k3(:,jspin3),
     >                   nv(jspin3),
     >                   rgphs,ustep,ig,-gb(:,ikpt_b,ikpt),I1)

            !write(*,'(a,3f8.5,3f8.5,3i1)')'k,k+b ',bkpt,bkpt_b,
     >      !                               jspin,jspin3,jspin4
            !write(*,*)'I1 int' 
            !do i=1,2
            !   do j=1,2
            !      write(*,*)I1(j,i)
            !   enddo
            !enddo

         !do i=1,nv_b(jspin)
         !   if(all(acof_b(i,:,:).eq.0.0)) write(*,*)'acof_b',i
         !   if(all(bcof_b(i,:,:).eq.0.0)) write(*,*)'bcof_b',i
         !enddo
         !do i=1,nv(jspin3)
         !   if(all(acof(i,:,:,jspin3).eq.0.0)) write(*,*)'acof',i
         !   if(all(bcof(i,:,:,jspin3).eq.0.0)) write(*,*)'bcof',i
         !enddo
       
 
   !**************************************!
   !       MT contribution at k+b1        !
   !**************************************!
         jspin5=jspin+2*(jspin3-1)
         call wann_mmkb_sph(nvd,llod,nv_b(jspin),nv(jspin3),nlod,
     >                      natd,ntypd,lmd,jmtd,taual,nop,lmax,
     >                      ntype,neq,nlo,llo,
     >                      acof_b,bcof_b,
     >                      ccof_b,bkpt,
     >                      acof(1:nv(jspin3),:,:,jspin3),
     >                      bcof(1:nv(jspin3),:,:,jspin3),
     >                      ccof(:,1:nv(jspin3),:,:,jspin3),
     >                      -gb(:,ikpt_b,ikpt),bkpt_b,
     <            ujug(:,:,:,:,jspin5),
     >            ujdg(:,:,:,:,jspin5),
     >            djug(:,:,:,:,jspin5),
     >            djdg(:,:,:,:,jspin5),
     <            ujulog(:,:,:,:,:,jspin5),
     >            djulog(:,:,:,:,:,jspin5),
     >            ulojug(:,:,:,:,:,jspin5),
     >            ulojdg(:,:,:,:,:,jspin5),
     >            ulojulog(:,:,:,:,:,:,jspin5),
     >                      kdiff,nntot,chi,I1)
       
            !write(*,*)'I1 tot' 
            !do i=1,2
            !   do j=1,2
            !      write(*,*)I1(j,i)
            !   enddo
            !enddo

         !I1 = conjg(I1)

         ! II1 = z_b^\dag . I1
         call zgemm('C','N',nslibd_b,nv(jspin3),nv_b(jspin),
     >              cmplx(1.0,0.0),
     >              z_b(1+addnoco_kb,1),nbasfcn,
     >              I1,nvd,cmplx(0.0,0.0),II1,neigd)


   !**************************************!
   !      loop over neighbors b2          !
   !**************************************!
         do ikpt_b2 = 1,2!nntot
            kptibz_b2=bpt(ikpt_b2,ikpt)
            n_start=1
            n_end=neigd

      !**************************************!
      !      read eigenvectors at k+b2       !
      !**************************************!   
#if(!defined(CPP_HDF) && defined(CPP_MPI))
            call wann_read_eig(
     >              lmaxd,ntypd,nlod,neigd,nvd,wannierspin,
     >              irank,isize,kptibz_b2,jspin_b,nbasfcn,nlotot,
     >              l_ss,l_noco,nrec_b,irecl,
     <              nmat_b2,nv_b2,ello,evdu,epar,kveclo_b2,
     <              k1_b2,k2_b2,k3_b2,bkpt_b2,wk_b2,nbands_b2,
     <              eigg,zz,cp_time,66,l_gwf,iqpt)
#else
            call cdn_read(
     >                lmaxd,ntypd,nlod,neigd,nvd,jspd,
     >                irank,isize,kptibz_b2,jspin_b,nbasfcn,nlotot,
     >                l_ss,l_noco,nrec_b,kptibz_b2,66,
     >                neigd,n_start,n_end,
     <                nmat_b2,nv_b2,ello,evdu,epar,kveclo_b2,
     <                k1_b2,k2_b2,k3_b2,bkpt_b2,wk_b2,nbands_b2,
     <                eigg,zz,cp_time)
#endif
!            write(*,'(a,f8.5,1x,f8.5,1x,f8.5)')'k+b2: ',bkpt_b2

            nslibd_b2 = 0
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
            z_b2(:,:) = cmplx(0.,0.)
#else
            z_b2(:,:) = 0.
#endif
            eig_b2(:) = 0.

            do i = 1,nbands_b2
               if((eigg(i).ge.e1s.and.nslibd_b2.lt.numbands
     &       .and.wann%l_bynumber)
     &       .or.(eigg(i).ge.e1s.and.eigg(i).le.e2s.and.wann%l_byenergy)
     &       .or.(i.ge.wann%band_min(jspin_b)
     &       .and.(i.le.wann%band_max(jspin_b))
     &       .and.wann%l_byindex))then
             nslibd_b2 = nslibd_b2 + 1
             eig_b2(nslibd_b2) = eigg(i)
             we_b2(nslibd_b2) = we_b2(i)
             if(l_noco)then
               funbas =        nv_b2(1) + nlotot
               funbas = funbas+nv_b2(2) + nlotot
             else
               funbas = nv_b2(jspin_b) + nlotot
             endif
             do  j = 1,funbas
              z_b2(j,nslibd_b2) = zz(j,i)
             enddo
               endif
            enddo
            noccbd_b2 = nslibd_b2

            addnoco_kb2=0 
            if(l_noco.and.jspin_b.eq.2) addnoco_kb2=nv_b2(1)+nlotot
      !**************************************!
      !      a,b,c coefficients at k+b2      !
      !**************************************!
            allocate ( acof_b2(nv_b2(jspin_b),0:lmd,natd),
     &             bcof_b2(nv_b2(jspin_b),0:lmd,natd),
     &             ccof_b2(-llod:llod,nv_b2(jspin_b),nlod,natd))

            allocate ( acof_tmp(noccbd_b2,0:lmd,natd),
     &           bcof_tmp(noccbd_b2,0:lmd,natd),
     &           ccof_tmp(-llod:llod,noccbd_b2,nlod,natd))     
 
            do i=1,nv_b2(jspin_b)
               z_tmp = 0.0
               z_tmp(i+addnoco_kb2,1)=1.0
               call abcof_small(
     >             lmaxd,ntypd,neigd,noccbd_b2,natd,nop,nvd,wannierspin,
     >             lmd,nbasfcn,llod,nlod,nlotot,invtab,
     >             ntype,mrot,ngopr1,taual,neq,lmax,rmt,omtil,
     >             bmat,bbmat,bkpt_b2,
     >             k1_b2,k2_b2,k3_b2,nv_b2,nmat_b2,noccbd_b2,z_tmp,
     >             us(0,1,jspin_b),dus(0,1,jspin_b),uds(0,1,jspin_b),
     >             duds(0,1,jspin_b),ddn(0,1,jspin_b),invsat,invsatnr,
     >             ulos(1,1,jspin_b),uulon(1,1,jspin_b),
     >             dulon(1,1,jspin_b),
     >             dulos(1,1,jspin_b),llo,nlo,l_dulo,lapw_l,
     >             l_noco,l_ss,jspin_b,alph_i,beta_i,qpt_i,
     >             kveclo_b2,odi,ods,
     <             acof_tmp(1,0,1),bcof_tmp(1,0,1),
     <             ccof_tmp(-llod,1,1,1),i)

               call wann_abinv(
     >         ntypd,natd,noccbd_b2,lmaxd,lmd,llod,nlod,ntype,neq,
     >         noccbd_b2,lmax,nlo,llo,invsat,invsatnr,bkpt_b2,taual,
     X         acof_tmp,bcof_tmp,ccof_tmp)

               acof_b2(i,:,:) = acof_tmp(1,:,:)
               bcof_b2(i,:,:) = bcof_tmp(1,:,:)
               ccof_b2(:,i,:,:)=ccof_tmp(:,1,:,:)
            enddo
            deallocate(acof_tmp,bcof_tmp,ccof_tmp)


         ! inner spin loop 2
         do jspin4=1,2
            addnoco2 = 0
            if(l_noco .and. jspin4.eq.2) addnoco2=nv(1)+nlotot

            I2 = cmplx(0.,0.)
      !**************************************!
      !       INT contribution at k+b2       !
      !**************************************! 
            call wann_uHu_int2( nvd,k1d,k2d,k3d,
     >                 n3d,k1(:,jspin4),k2(:,jspin4),k3(:,jspin4),
     >                 nv(jspin4),
     >                 k1_b2(:,jspin_b),k2_b2(:,jspin_b),
     >                 k3_b2(:,jspin_b),nv_b2(jspin_b),
     >                 rgphs,ustep,ig,gb(:,ikpt_b2,ikpt),I2)

            !write(*,'(a,3f8.5,3f8.5,3i1)')'k,k+b2 ',bkpt,bkpt_b2,
     >      !                               jspin,jspin3,jspin4
            !write(*,*)'I2 int' 
            !do i=1,2
            !   do j=1,2
            !      write(*,*)I2(j,i)
            !   enddo
            !enddo

         !do i=1,nv_b2(jspin)
         !   if(all(acof_b2(i,:,:).eq.0.0)) write(*,*)'acof_b',i
         !   if(all(bcof_b2(i,:,:).eq.0.0)) write(*,*)'bcof_b',i
         !enddo
         !do i=1,nv(jspin4)
         !   if(all(acof(i,:,:,jspin4).eq.0.0)) write(*,*)'acof',i
         !   if(all(bcof(i,:,:,jspin4).eq.0.0)) write(*,*)'bcof',i
         !enddo
        
      !**************************************!
      !       MT contribution at k+b2        !
      !**************************************!
            jspin5=jspin4+2*(jspin_b-1)
            call wann_mmkb_sph(nvd,llod,nv(jspin4),nv_b2(jspin_b),nlod,
     >                      natd,ntypd,lmd,jmtd,taual,nop,lmax,
     >                      ntype,neq,nlo,llo,
     >                      acof(1:nv(jspin4),:,:,jspin4),
     >                      bcof(1:nv(jspin4),:,:,jspin4),
     >                      ccof(:,1:nv(jspin4),:,:,jspin4),bkpt_b2,
     >                      acof_b2,bcof_b2,
     >                      ccof_b2,gb(:,ikpt_b2,ikpt),bkpt,
     <            ujug(:,:,:,:,jspin5),
     >            ujdg(:,:,:,:,jspin5),
     >            djug(:,:,:,:,jspin5),
     >            djdg(:,:,:,:,jspin5),
     <            ujulog(:,:,:,:,:,jspin5),
     >            djulog(:,:,:,:,:,jspin5),
     >            ulojug(:,:,:,:,:,jspin5),
     >            ulojdg(:,:,:,:,:,jspin5),
     >            ulojulog(:,:,:,:,:,:,jspin5),
     >                      kdiff,nntot,chi,I2)                 
       
            !write(*,*)'I2 tot'
            !do i=1,2
            !   do j=1,2
            !      write(*,*)I2(j,i)
            !   enddo
            !enddo

            !I2 = conjg(I2)

            ! II2 = I2 . z_b2
            call zgemm('N','N',nv(jspin4),nslibd_b2,nv_b2(jspin_b),
     >              cmplx(1.0,0.0),I2,nvd,
     >              z_b2(1+addnoco_kb2,1),nbasfcn,cmplx(0.0,0.0),
     >              II2,nvd)

            ! III1 = II1 . H
            call zgemm('N','N',nslibd_b,nv(jspin4),nv(jspin3),
     >              cmplx(1.0,0.0),II1,neigd,
     >              ham(1+addnoco,1+addnoco2),nbasfcn2,
     >              cmplx(0.0,0.0),III1,neigd)

            !write(*,*)'III1'
            !do i=1,2
            !   do j=1,2
            !      write(*,*)III1(j,i)
            !   enddo
            !enddo

            tmp_mat2 = cmplx(0.,0.)
            ! uHu = III1 . II2 = z_b^\dag . I1 . H . I2 . z_b2
            call zgemm('N','N',nslibd_b,nslibd_b2,nv(jspin4),
     >              cmplx(1.0,0.0),III1,neigd,
     >              II2,nvd,cmplx(0.0,0.0),
     >              tmp_mat2,nbnd)
            uHu(:,:,ikpt_b2,ikpt_b) = uHu(:,:,ikpt_b2,ikpt_b)
     >                              + tmp_mat2
            !write(*,'(a,i3,i3,i3)')'(k,k+b1,k+b2): ',ikpt,ikpt_b,ikpt_b2
            !do i=1,2
            !   do j=1,2
            !      write(*,'(i2,i2,2(f20.14))')j,i,tmp_mat2(j,i)
            !   enddo
            !enddo
            ! uHu(n,m,i,j) = < u_{jn} | H | u_{im} >

         enddo!jspin4

            deallocate(acof_b2,bcof_b2,ccof_b2)
         enddo ! loop b2
      enddo!jspin3

         deallocate(acof_b,bcof_b,ccof_b)
      enddo ! loop b1

      !write(*,'(a,i2,i2)')'spins:',jspin,jspin_b
      !do ikpt_b=1,2
      !   do ikpt_b2=1,2
      !do i=1,2
      !   do j=1,2
      !      write(*,'(i2,i2,i2,i2,2f20.14)')j,i,ikpt_b2,ikpt_b,
     >!                                  uHu(j,i,ikpt_b2,ikpt_b)
      !   enddo
      !enddo
      !   enddo
      !enddo

      deallocate(acof,bcof,ccof)
      deallocate(tmp_mat2)
      deallocate(ham)

      END SUBROUTINE wann_orbmag

 
      END MODULE m_wann_orbmag


