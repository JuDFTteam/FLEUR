      module m_wann_hopping
      use m_juDFT
      contains 
      subroutine wann_hopping(
     >          rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >          l_soc,band_min,band_max,neigd,
     >          l_socmmn0,l_ndegen,ndegen,wan90version,
     >          l_unformatted)
c****************************************************
c     Calculate the hoppings using the information
c     in file WF1.chk produced by wannier90.
c
c     Frank Freimuth
c****************************************************
      use m_constants, only:pimach
      use m_wann_read_umatrix
c      use m_wann_get_mp
c      use m_wann_get_kpts
c      use m_wann_wigner_seitz

      implicit none
      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)
      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical, intent(in) :: l_bzsym,l_soc
      logical, intent(in) :: film,l_onedimens
      integer, intent(in) :: band_min(2),band_max(2),neigd
      logical, intent(in) :: l_socmmn0
      logical, intent(in) :: l_ndegen
      integer, intent(in) :: ndegen(:)
      integer, intent(in) :: wan90version
      logical, intent(in) :: l_unformatted	

      integer             :: ikpt,jspins
      integer             :: kpts
      logical             :: l_file
c      real                :: kpoints(3,nkpts)
      integer             :: num_wann,num_kpts,num_nnmax,jspin
      integer             :: kspin,kkspin
      integer             :: wann_shift,num_wann2
      integer             :: i,j,k,m,info,r1,r2,r3,dummy1
      integer             :: dummy2,dummy3
      integer             :: counter,m1,m2
      integer             :: num_bands2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:),hwann(:,:,:)
      complex,allocatable :: hreal(:,:,:)
      complex,allocatable :: hrealsoc(:,:,:,:,:)
      complex,allocatable :: hwannsoc(:,:,:,:,:)
      complex,allocatable :: mmn0(:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      complex             :: value(4)
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
      integer             :: num_bands,num(3)
      logical             :: l_umdat
      real,allocatable    :: eigval2(:,:)
      real,allocatable    :: eigval_opt(:,:)
      real                :: scale,a,b
      character(len=2)    :: spinspin12(0:2)
      character(len=3)    :: spin12(2)
      character(len=6)    :: filename
      logical             :: l_socham
!      integer             :: hopmin_x,hopmin_y,hopmin_z
!      integer             :: hopmax_x,hopmax_y,hopmax_z
      integer             :: rvecind!,rvecnum
!      integer,allocatable :: rvec(:,:)
      integer             :: ii,int_dummy

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      um_format=.false.   !if you would like to get a formatted um_dat
      repro_eig=.false.   !if you would like to check unitarity of u_matrix

      tpi=2*pimach()

      jspins=jspins_in
      if(l_soc)jspins=1
      l_socham=.false.

      write(6,*)"nkpts=",nkpts

      do jspin=1,jspins  !spin loop
c*****************************************************
c     get num_bands and num_wann from the proj file
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
            CALL juDFT_error("no proj/proj.1/proj.2",calledby
     +           ="wann_hopping")
         endif
         read (203,*) num_wann,num_bands
         close (203)
         write(6,*)'According to proj there are ',num_bands,' bands'
         write(6,*)"and ",num_wann," wannier functions."
c****************************************************************
c     read in chk
c****************************************************************
         num_kpts=nkpts
         allocate( u_matrix_opt(num_bands,num_wann,nkpts) )
         allocate( u_matrix(num_wann,num_wann,nkpts) )
         allocate( lwindow(num_bands,nkpts) )
         allocate( ndimwin(nkpts) )
         call wann_read_umatrix2(
     >               nkpts,num_wann,num_bands,
     >               um_format,jspin,wan90version,
     <               have_disentangled,
     <               lwindow,ndimwin,
     <               u_matrix_opt,u_matrix)
c****************************************************************
c     read in eig-file
c****************************************************************
         num_bands2=num_bands
         if(l_soc.and.l_socmmn0)then
          num_bands2=neigd
         endif
         write(6,*)"read in eig-file"
         allocate(energy(num_bands2,num_kpts))
         inquire(file=spin12(jspin)//'.eig',exist=l_umdat)
         IF(.NOT.l_umdat)  CALL juDFT_error
     +        ("Thou shall not hide your eig file",calledby
     +        ="wann_hopping")
         open(300,file=spin12(jspin)//'.eig',form='formatted')
         do i=1,num_kpts
           do j=1,num_bands2
              read(300,*)band,kpt,energy(j,i)
           enddo
         enddo
         close(300)

         minenerg=minval(energy(:,:))
         maxenerg=maxval(energy(:,:))
         write(6,*)"minenerg=",minenerg
         write(6,*)"maxenerg=",maxenerg
c*********************************************
c       Preparations for spin-orbit coupling
c*********************************************
         if(l_socham)then
           allocate(mmn0(nkpts,num_bands2,num_bands2,2))
           open(304,file='spinup.mmn0',form='formatted')
           read(304,*)
           read(304,*)
           do nkp=1,num_kpts
            do i=1,num_bands2
             do j=1,num_bands2
                read(304,*)dummy1,dummy2,dummy3,a,b
c                mmn0(nkp,j,i,1)=cmplx(a,b)
                mmn0(nkp,i,j,1)=cmplx(a,-b)
             enddo !j
            enddo !i
           enddo !nkp
           close(304)
           open(304,file='spindown.mmn0',form='formatted')
           read(304,*)
           read(304,*)
           do nkp=1,num_kpts
            do i=1,num_bands2
             do j=1,num_bands2
                read(304,*)dummy1,dummy2,dummy3,a,b
c                mmn0(nkp,j,i,2)=cmplx(a,b)
                mmn0(nkp,i,j,2)=cmplx(a,-b)
             enddo !j
            enddo !i
           enddo !nkp
           close(304)
         endif!l_soc
c****************************************************************
c               calculate matrix elements of hamiltonian
c****************************************************************

         write(6,*)"calculate matrix elements of hamiltonian
     &   between wannier orbitals"


         allocate(eigval_opt(num_bands,nkpts))
         allocate(eigval2(num_bands2,nkpts))
         eigval_opt=0.0
         eigval2=0.0


         if(have_disentangled) then

         ! slim down eigval to contain states within the outer window

          do nkp=1,num_kpts
            counter=0
            do j=1,num_bands
              if(lwindow(j,nkp)) then
                counter=counter+1
                eigval_opt(counter,nkp)=energy(j,nkp)
              end if
            end do
          end do
       
          ! rotate eigval into the optimal subspace
          ! in general eigval would be a matrix at each kpoints
          ! but we choose u_matrix_opt such that the Hamiltonian is
          ! diagonal at each kpoint. (I guess we should check it here)
       
          do nkp=1,num_kpts
           do j=1,num_wann
             do m=1,ndimwin(nkp)
                eigval2(j,nkp)=eigval2(j,nkp)+eigval_opt(m,nkp)* 
     &    real(conjg(u_matrix_opt(m,j,nkp))*u_matrix_opt(m,j,nkp))
             enddo
           enddo
          enddo

         else
          eigval2(1:num_bands2,:)=energy(1:num_bands2,:)
         end if                    !have_disentangled

         deallocate(eigval_opt)
         deallocate(energy)

         allocate(hwann(num_wann,num_wann,num_kpts))
         hwann=cmplx(0.0,0.0)
         wann_shift=0
         if(l_socmmn0)then
            wann_shift=band_min(jspin)-1
         endif

         do k=1,num_kpts
          do i=1,num_wann
           do j=1,num_wann
            do m=1,num_wann
                    hwann(i,j,k)=hwann(i,j,k)
     +    +eigval2(m+wann_shift,k)*u_matrix(m,j,k)*
     *     conjg(u_matrix(m,i,k))
            enddo
           enddo
          enddo
         enddo

         if(l_socham)then

         num_wann2=num_wann
         wann_shift=0
         if(l_socmmn0)then
          num_wann2=neigd
          wann_shift=band_min(jspin)-1
         endif
         write(6,*)"num_wann2=",num_wann2
         write(6,*)"num_wann=",num_wann
         write(6,*)"wann_shift=",wann_shift
         allocate(hwannsoc(num_wann,num_wann,num_kpts,2,2))
         hwannsoc=cmplx(0.0,0.0)
         do kspin=1,2
          do kkspin=1,2
           do i=1,num_wann
            do j=1,num_wann
             do k=1,num_kpts
              do m=1,num_wann2
              do m1=1,num_wann
              do m2=1,num_wann
               hwannsoc(i,j,k,kkspin,kspin)=
     =          hwannsoc(i,j,k,kkspin,kspin)
     +           +eigval2(m,k)*u_matrix(m2,j,k)*conjg(u_matrix(m1,i,k))*
     *           mmn0(k,m1+wann_shift,m,kkspin)*
     *           mmn0(k,m,m2+wann_shift,kspin)
              enddo
              enddo
              enddo !m
             enddo !k
            enddo !j
           enddo !i
          enddo !kkspin
         enddo !kspin
         deallocate(mmn0)


         do i=1,num_wann
          do j=1,num_wann
           do k=1,num_kpts
             eulav=hwannsoc(i,j,k,1,1)+hwannsoc(i,j,k,1,2)+
     &                   hwannsoc(i,j,k,2,1)+hwannsoc(i,j,k,2,2)
             eulav=eulav-hwann(i,j,k)
             if(abs(eulav).gt.1.e-6)then
                 write(6,*)"soc-hop: something wrong:",eulav
             endif
           enddo
          enddo
         enddo



         endif !l_soc




c***************************************************************
c        repro_eig
c***************************************************************

         if(repro_eig)then
   
         write(6,*)"As a check, try to reproduce the old eigenvalues"

c       Note, that this check will give positive result,
c       if u_matrix is unitary

        allocate (vec(num_wann,num_wann),ei(num_wann))
        allocate(eigw(num_kpts,num_wann))

        lee = log( dble(num_wann) )/log(2.d0) + 1
        lwork = 1 + 5*num_wann + 2*num_wann*lee + 3*(num_wann**2)
        lrwork = 1 + 4*num_wann + 2*num_wann*lee + 3*(num_wann**2)
        liwork = 2 + 5*num_wann +1

        allocate (work(lwork),rwork(lrwork),iwork(liwork))



        do i=1,num_kpts
          do j = 1,num_wann
            do k = 1,num_wann
               vec(j,k) = hwann(j,k,i)
            enddo
          enddo


          jobz = 'V' ; uplo = 'L' ; n = num_wann ; lda = num_wann
          call zheevd(jobz,uplo,n,vec,lda,ei,work,lwork,
     &             rwork,lrwork,iwork,liwork,info)

          if (info.lt.0) write (*,*)
     &               'ith argument had an illegal value ',info
          IF (info>0)  CALL juDFT_error("not converged diagonalization"
     +         ,calledby ="wann_hopping")

          do j = 1,num_wann
              eigw(i,j) = ei(j)
              if(abs(eigw(i,j)-eigval2(j,i)).gt.0.0001)then
                 write(6,*)"found different eigenvalues:"
                 write(6,*)"kpt=",i
                 write(6,*)"band=",j
                 write(6,*)"eig=",eigval2(j,i)
                 write(6,*)"neweig=",eigw(i,j)
              endif
          enddo

        enddo
        deallocate(work,rwork,iwork)
        deallocate(vec,ei)
        open(500,file='reeig'//spin12(jspin),form='formatted')
        do i=1,num_kpts
          do j=1,num_wann
             write(500,'(i5,3x,i5,3x,f20.16)')i,j,eigw(i,j)
          enddo
        enddo
        close(500)
        deallocate(eigw)
      endif !repro_eig
c************************************************************
c     Calculate hoppings.
c***********************************************************      
      write(6,*)"calculate hoppings"
      allocate( hreal(num_wann,num_wann,rvecnum) )
      hreal=cmplx(0.0,0.0)
      do rvecind=1,rvecnum
       do k=1,nkpts
          rdotk=tpi*(  kpoints(1,k)*rvec(1,rvecind)+
     +                 kpoints(2,k)*rvec(2,rvecind)+
     +                 kpoints(3,k)*rvec(3,rvecind)  )
          fac=cmplx(cos(rdotk),-sin(rdotk))
          do m2=1,num_wann
           do m1=1,num_wann
            hreal(m1,m2,rvecind)=hreal(m1,m2,rvecind)+
     &            fac*hwann(m1,m2,k)
           enddo !m1
          enddo !m2
       enddo !k
      enddo !rvecind
      hreal=hreal/cmplx(real(nkpts),0.0)

      if(l_ndegen)then
        do rvecind=1,rvecnum
          do m2=1,num_wann
            do m1=1,num_wann
              hreal(m1,m2,rvecind)=
     & hreal(m1,m2,rvecind)/cmplx(real(ndegen(rvecind)),0.0)
            enddo !m1
          enddo !m2
        enddo !rvecind  
      endif !l_ndegen

c************************************************************
c              calculate hoppings for soc
c***********************************************************   
      if(l_socham)then
       allocate(hrealsoc(num_wann,num_wann,rvecnum,
     &         2,2))
       hrealsoc=cmplx(0.0,0.0)
       do rvecind=1,rvecnum
         do k=1,nkpts
           rdotk=tpi*(  kpoints(1,k)*rvec(1,rvecind)+
     +                  kpoints(2,k)*rvec(2,rvecind)+
     +                  kpoints(3,k)*rvec(3,rvecind)  )
           fac=cmplx(cos(rdotk),-sin(rdotk))/nkpts
           do kspin=1,2
             do kkspin=1,2
               hrealsoc(:,:,rvecind,kspin,kkspin)=
     =           hrealsoc(:,:,rvecind,kspin,kkspin)+
     +           fac*hwannsoc(:,:,k,kspin,kkspin)
             enddo
           enddo
         enddo !k
       enddo !rvecnum
      endif !l_soc

c****************************************************************
c           make the hoppings real
c****************************************************************
c$$$      if(.false.)then
c$$$      do i=1,num_wann
c$$$       do j=i+1,num_wann
c$$$        eulav=cmplx(0.0,0.0)
c$$$        do r1=hopmin,hopmax
c$$$         do r2=hopmin,hopmax
c$$$          do r3=hopmin,hopmax
c$$$            eulav1=hreal(i,j,r1,r2,r3)
c$$$            if(abs(eulav1).gt.abs(eulav))then
c$$$               eulav=eulav1
c$$$            endif
c$$$          enddo !r1
c$$$         enddo !r2
c$$$        enddo !r3
c$$$        if(abs(eulav).gt.1e-6)then
c$$$         eulav=eulav/abs(eulav)
c$$$         do r1=hopmin,hopmax
c$$$          do r2=hopmin,hopmax
c$$$           do r3=hopmin,hopmax
c$$$             hreal(i,j,r1,r2,r3)=hreal(i,j,r1,r2,r3)/eulav
c$$$             hreal(j,i,-r1,-r2,-r3)=conjg(hreal(i,j,r1,r2,r3))
c$$$!             hreal(j,i,r1,r2,r3)=hreal(j,i,r1,r2,r3)/conjg(eulav)
c$$$             if(l_soc)then
c$$$               hrealsoc(i,j,r1,r2,r3,:,:)=
c$$$     &                  hrealsoc(i,j,r1,r2,r3,:,:)/eulav
c$$$               hrealsoc(j,i,-r1,-r2,-r3,:,:)=
c$$$     &                  conjg( hrealsoc(i,j,r1,r2,r3,:,:) )
c$$$!               hrealsoc(j,i,r1,r2,r3,:,:)=
c$$$!     &                  hrealsoc(j,i,r1,r2,r3,:,:)/conjg(eulav)
c$$$             endif
c$$$           enddo !r3
c$$$          enddo !r2   
c$$$         enddo !r1      
c$$$        endif !eulav
c$$$       enddo !j
c$$$      enddo !i
c$$$      endif

      if(l_socham)then
       do kspin=1,2
       do kkspin=1,2
       open(321,file='hopping'//spinspin12(kspin)//spinspin12(kkspin),
     &                        form='formatted')
       do rvecind=1,rvecnum
          r1=rvec(1,rvecind)
          r2=rvec(2,rvecind)
          r3=rvec(3,rvecind)
          do i=1,num_wann
           do j=1,num_wann
            write(321,'(i3,i3,i3,i3,i3,f20.8,f20.8)')
     &        r1,r2,r3,i,j,hrealsoc(i,j,rvecind,kspin,kkspin) 
           enddo !j
          enddo !i
       enddo !rvecind   
       close(321)       
       enddo  !kkspin
       enddo  !kspin
       deallocate(hrealsoc,hwannsoc)
      endif !l_soc  

      if(l_unformatted)then
       if(l_soc)then  
        if(l_ndegen)then		
         open(321,file='hunfndegen',form='unformatted'
#ifdef CPP_INTEL
     &                  ,convert='BIG_ENDIAN'
#endif
     &                   )
        else
         open(321,file='hunf',form='unformatted'
#ifdef CPP_INTEL
     &                  ,convert='BIG_ENDIAN'
#endif
     &                   )
        endif	
       else
        if(l_ndegen)then		
         open(321,file='hunfndegen'//spinspin12(jspin),
     &                  form='unformatted'
#ifdef CPP_INTEL
     &                  ,convert='BIG_ENDIAN'
#endif
     &                   )
        else
         open(321,file='hunf'//spinspin12(jspin),
     &                  form='unformatted'
#ifdef CPP_INTEL
     &                  ,convert='BIG_ENDIAN'
#endif
     &                   )
        endif	
       endif
       write(321)num_wann
       write(321)rvecnum
       write(321)rvec
       write(321)hreal	
      else	
       if(l_ndegen)then
        open(321,file='hopping_ndegen'//spinspin12(jspin),
     &            form='formatted')
       else
        open(321,file='hopping'//spinspin12(jspin),
     &            form='formatted')
       endif   
       do rvecind=1,rvecnum
         r3=rvec(3,rvecind)
         r2=rvec(2,rvecind)
         r1=rvec(1,rvecind)
         do j=1,num_wann
           do i=1,num_wann
            write(321,'(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,hreal(i,j,rvecind) 
           enddo
         enddo   
       enddo !rvecnum
      endif !l_unformatted 	 
      close(321)

      deallocate(lwindow,u_matrix_opt,ndimwin)
      deallocate(eigval2,u_matrix,hwann,hreal)
!      deallocate(rvec)
      enddo !jspin



      end subroutine wann_hopping
      end module m_wann_hopping
