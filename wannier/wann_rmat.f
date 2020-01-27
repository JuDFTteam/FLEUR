!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_rmat
      USE m_juDFT
      contains
      subroutine wann_rmat(
     >          bmat,amat, 
     >          rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >          l_nocosoc,band_min,band_max,neigd,
     >          l_socmmn0,wan90version)
c***********************************************************
c     Calculate the matrix elements of the position operator
c     between wannier functions.
c
c     FF, 2009
c***********************************************************
      use m_constants, only:pimach, ImagUnit
      use m_wann_read_umatrix

      implicit none
      real, intent(in) :: bmat(:,:)
      real, intent(in) :: amat(:,:)

      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)
      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical, intent(in) :: l_bzsym
      logical, intent(in) :: film
      logical, intent(in) :: l_onedimens

      logical, intent(in) :: l_nocosoc
      integer, intent(in) :: band_min(2),band_max(2),neigd

      logical, intent(in) :: l_socmmn0
      integer, intent(in) :: wan90version

      real,allocatable    :: bpunkt(:,:)
      real                :: kdiff(3)
      real,allocatable    :: wb(:)
      integer,allocatable :: gb(:,:,:),bpt(:,:)
      integer             :: ikpt,jspins,ikpt_help
      integer             :: kpts,nntot,nn
      logical             :: l_file
      integer             :: num_wann,num_kpts,num_nnmax,jspin
      integer             :: kspin,kkspin
      integer             :: num_wann2
      integer             :: i,j,k,m,info,r1,r2,r3,dummy1
      integer             :: dummy2,dummy3,dummy4
      integer             :: hopmin,hopmax,counter,m1,m2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:),m_matrix(:,:,:,:)
      complex,allocatable :: posop(:,:,:,:),posop2(:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      real                :: rvalue
      logical             :: um_format
      logical             :: repro_eig
      logical             :: l_chk,l_proj
      logical             :: have_disentangled
      integer,allocatable :: ndimwin(:)
      logical,allocatable :: lwindow(:,:)
      integer             :: chk_unit,nkp,ntmp,ierr
      character(len=33)   :: header
      character(len=20)   :: checkpoint
      real                :: tmp_latt(3,3), tmp_kpt_latt(3,nkpts)
      real                :: omega_invariant
      complex,allocatable :: u_matrix_opt(:,:,:)
      integer             :: num_bands,nn2
      logical             :: l_umdat
      real,allocatable    :: eigval2(:,:)
      real,allocatable    :: eigval_opt(:,:)
      real                :: scale,a,b
      character(len=2)    :: spinspin12(0:2)
      character(len=3)    :: spin12(2)
      character(len=6)    :: filename
      integer             :: jp,mp,kk
      integer             :: rvecind,rvecind_0
      real                :: realvec(3)
      complex,allocatable :: mmnk(:,:,:,:)
      logical             :: l_worksout,l_she

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      tpi=2*pimach()

      jspins=jspins_in
      if(l_nocosoc)jspins=1

      write(6,*)"nkpts=",nkpts

      do jspin=1,jspins  !spin loop
c*****************************************************
c        get num_bands and num_wann from the proj file
c*****************************************************
         do j=jspin,0,-1
          inquire(file=trim('proj'//spinspin12(j)),exist=l_file)
          if(l_file)then
            filename='proj'//spinspin12(j)
            exit
          endif
         enddo
         if(l_file)then
          open (203,file=trim(filename),status='old')
          rewind (203)
         else
            CALL judft_error("no proj/proj.1/proj.2",calledby
     +           ="wann_rmat")
         endif
         read (203,*) num_wann,num_bands
         close (203)
         write(6,*)'According to proj there are ',num_bands,' bands'
         write(6,*)"and ",num_wann," wannier functions."

c****************************************************************
c        get nntot and bk and wb from bkpts file
c****************************************************************
         inquire (file='bkpts',exist=l_file)
         if (.not.l_file)  CALL judft_error("need bkpts"
     +      ,calledby ="wann_rmat")
         open (202,file='bkpts',form='formatted',status='old')
         read (202,'(i4)') nntot
         allocate ( gb(3,nntot,nkpts) )
         allocate ( bpt(nntot,nkpts) )
         do ikpt=1,nkpts
          do nn=1,nntot
           read (202,'(2i6,3x,3i4)')
     &      ikpt_help,bpt(nn,ikpt),(gb(i,nn,ikpt),i=1,3)
          enddo !nn
         enddo !ikpt
         allocate( wb(nntot) )
         allocate( bpunkt(3,nntot) )
         do nn=1,nntot
            read(202,*)bpunkt(:,nn),wb(nn)
         enddo
         close (202)

         do ikpt=1,nkpts
          do nn=1,nntot
           kdiff(1)=kpoints(1,bpt(nn,ikpt))+
     +            gb(1,nn,ikpt) - kpoints(1,ikpt)
           kdiff(2)=kpoints(2,bpt(nn,ikpt))+
     +            gb(2,nn,ikpt) - kpoints(2,ikpt)
           kdiff(3)=kpoints(3,bpt(nn,ikpt))+
     +            gb(3,nn,ikpt) - kpoints(3,ikpt)
           kdiff=matmul(transpose(bmat),kdiff)
           l_worksout=.false.
           do nn2=1,nntot
             if (abs(kdiff(1)-bpunkt(1,nn2)).lt.1.e-5)then
               if (abs(kdiff(2)-bpunkt(2,nn2)).lt.1.e-5)then
                 if (abs(kdiff(3)-bpunkt(3,nn2)).lt.1.e-5)then
                     l_worksout=.true. 
                     exit
                 endif          
               endif
             endif 
           enddo !nn2
           if(.not.l_worksout)then
              write(*,*)"ikpt,nn=",ikpt,nn
              write(*,*)"kdiff(1)=",kdiff(1)
              write(*,*)"kdiff(2)=",kdiff(2)
              write(*,*)"kdiff(3)=",kdiff(3)
              stop 'worksout'
           endif   
          enddo !nn
         enddo !ikpt 

c****************************************************************
c        read in chk
c****************************************************************
         num_kpts=nkpts
         allocate( u_matrix_opt(num_bands,num_wann,nkpts) )
         allocate( u_matrix(num_wann,num_wann,nkpts) )
         allocate( lwindow(num_bands,nkpts) )
         allocate( ndimwin(nkpts) )
         allocate( m_matrix(num_wann,num_wann,nntot,num_kpts) )
         call wann_read_umatrix2(
     >            nkpts,num_wann,num_bands,
     >            um_format,jspin,wan90version,
     <            have_disentangled,
     <            lwindow,ndimwin,
     <            u_matrix_opt,u_matrix,m_matrix)
c****************************************************************
c        read in eig-file
c****************************************************************
         write(6,*)"read in eig-file"
         allocate(energy(num_bands,num_kpts))
         inquire(file=spin12(jspin)//'.eig',exist=l_umdat)
         IF(.NOT.l_umdat)  CALL judft_error
     +        ("Thou shall not hide your eig file",calledby
     +        ="wann_hopping")
         open(300,file=spin12(jspin)//'.eig',form='formatted')
         do i=1,num_kpts
           do j=1,num_bands
              read(300,*)band,kpt,energy(j,i)
           enddo
         enddo
         close(300)

         minenerg=minval(energy(:,:))
         maxenerg=maxval(energy(:,:))
         write(6,*)"minenerg=",minenerg
         write(6,*)"maxenerg=",maxenerg


         allocate(eigval_opt(num_bands,nkpts))
         allocate(eigval2(num_wann,nkpts))
         eigval_opt=0.0
         eigval2=0.0

         if(have_disentangled) then

           do nkp=1,num_kpts
            counter=0
            do j=1,num_bands
              if(lwindow(j,nkp)) then
                counter=counter+1
                eigval_opt(counter,nkp)=energy(j,nkp)
              end if
            end do
           end do
       
           do nkp=1,num_kpts
            do j=1,num_wann
             do m=1,ndimwin(nkp)
                eigval2(j,nkp)=eigval2(j,nkp)+eigval_opt(m,nkp)* 
     &    real(conjg(u_matrix_opt(m,j,nkp))*u_matrix_opt(m,j,nkp))
             enddo
            enddo
           enddo

         else
           eigval2 = energy
         end if                    !have_disentangled

         deallocate(eigval_opt)
         deallocate(energy)

c****************************************************************
c        Set up posop.
c****************************************************************
         write(6,*)"Set up posop."
         allocate( posop2(3,num_wann,num_wann,num_kpts) )
         posop2=cmplx(0.0,0.0)
         
         do ikpt=1,nkpts
          do nn=1,nntot
           kdiff(1)=kpoints(1,bpt(nn,ikpt))+
     +            gb(1,nn,ikpt) - kpoints(1,ikpt)
           kdiff(2)=kpoints(2,bpt(nn,ikpt))+
     +            gb(2,nn,ikpt) - kpoints(2,ikpt)
           kdiff(3)=kpoints(3,bpt(nn,ikpt))+
     +            gb(3,nn,ikpt) - kpoints(3,ikpt)
           kdiff=matmul(transpose(bmat),kdiff)          
           do i=1,num_wann
            do j=1,num_wann
             if(j.eq.i)then
              do kk=1,3
               posop2(kk,j,i,ikpt)=posop2(kk,j,i,ikpt)+ImagUnit*
     &             wb(nn)*kdiff(kk)*(m_matrix(j,i,nn,ikpt)-1.0)  
              enddo !kk
             else
              do kk=1,3
               posop2(kk,j,i,ikpt)=posop2(kk,j,i,ikpt)+ImagUnit*
     &             wb(nn)*kdiff(kk)*m_matrix(j,i,nn,ikpt)  
              enddo !kk
             endif
            enddo !j
           enddo !i
          enddo !nn   
         enddo !ikpt

         do ikpt=1,nkpts
          do i=1,num_wann
           do j=1,i
            do kk=1,3
              posop2(kk,j,i,ikpt)=(posop2(kk,j,i,ikpt)+
     &                        conjg(posop2(kk,i,j,ikpt)))/2.0
              posop2(kk,i,j,ikpt)=posop2(kk,j,i,ikpt)
            enddo !kk
           enddo !j
          enddo !i
         enddo !ikpt

         allocate( posop(3,num_wann,num_wann,rvecnum) )
         posop=cmplx(0.0,0.0)
         do rvecind=1,rvecnum
            do k=1,nkpts  
              rdotk=tpi*( kpoints(1,k)*rvec(1,rvecind)+
     &                    kpoints(2,k)*rvec(2,rvecind)+
     &                    kpoints(3,k)*rvec(3,rvecind) )
              fac=cmplx(cos(rdotk),-sin(rdotk))
              do i=1,num_wann
               do j=1,num_wann
                do kk=1,3
                 posop(kk,j,i,rvecind)=
     &                   posop(kk,j,i,rvecind)+
     &                   fac*posop2(kk,j,i,k)
                enddo !kk
               enddo !j
              enddo !i
            enddo !k
         enddo !rvecind
         posop=posop/cmplx(real(nkpts),0.0)
         deallocate( posop2 )

         do rvecind=1,rvecnum
          if(rvec(1,rvecind).eq.0)then
           if(rvec(2,rvecind).eq.0)then
            if(rvec(3,rvecind).eq.0)then
               rvecind_0=rvecind
               goto 123
            endif
           endif
          endif
         enddo !rvecind
         stop 'Ou est ce point-la ?'
 123     continue

         do i=1,num_wann
          do kk=1,3 
             posop(kk,i,i,rvecind_0)=0.0
          enddo !kk
         enddo !i

         do ikpt=1,nkpts
          do nn=1,nntot  
           kdiff(1)=kpoints(1,bpt(nn,ikpt))+
     +            gb(1,nn,ikpt) - kpoints(1,ikpt)
           kdiff(2)=kpoints(2,bpt(nn,ikpt))+
     +            gb(2,nn,ikpt) - kpoints(2,ikpt)
           kdiff(3)=kpoints(3,bpt(nn,ikpt))+
     +            gb(3,nn,ikpt) - kpoints(3,ikpt)
           kdiff=matmul(transpose(bmat),kdiff)/real(nkpts)             
           do i=1,num_wann
            do kk=1,3 
             posop(kk,i,i,rvecind_0)=posop(kk,i,i,rvecind_0)-
     &        wb(nn)*kdiff(kk)*aimag(log(m_matrix(i,i,nn,ikpt)))
            enddo !kk
           enddo !i
          enddo !nn
         enddo !ikpt

c********************************************
c        Print posop.
c********************************************
         open(321,file='posop'//spinspin12(jspin),form='formatted')
         do rvecind=1,rvecnum
            r3=rvec(3,rvecind)
            r2=rvec(2,rvecind)
            r1=rvec(1,rvecind)             
            do j=1,num_wann
             do i=1,num_wann
              do kk=1,3  
               write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &                 1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &                 r1,r2,r3,i,j,kk,posop(kk,i,j,rvecind) 
              enddo !kk  
             enddo !j
            enddo !i
         enddo !rvecind  
         close(321)

         deallocate( lwindow,u_matrix_opt,ndimwin )
         deallocate( u_matrix,m_matrix )
         deallocate( posop,eigval2 )
         deallocate( gb,bpt,wb,bpunkt )
      enddo !jspin

      end subroutine wann_rmat
      end module m_wann_rmat
