      module m_wann_torque_rs
        use m_juDFT
      contains 
      subroutine wann_torque_rs( 
     >          ntype,neq,rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,
     >          l_soc,band_min,band_max,neigd,
     >          l_ndegen,ndegen,wan90version,
     >          l_unformatted)
c*************************************************
c     Calculate the matrix elements of the 
c     muffin-tin spin current operator
c     in real space from the
c     files WF1.chk (and WF1_um.dat) produced
c     by wannier90.
c
c     Frank Freimuth, January 2011
c*************************************************
      use m_constants, only:pimach
      use m_wann_read_umatrix

      implicit none
      integer, intent(in) :: ntype
      integer, intent(in) :: neq(:)
      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)

      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical, intent(in) :: l_bzsym,l_soc
      logical, intent(in) :: film
      

      integer, intent(in) :: band_min(2),band_max(2),neigd

      logical, intent(in) :: l_ndegen
      integer, intent(in) :: ndegen(:)
      integer, intent(in) :: wan90version
      logical, intent(in) :: l_unformatted

      integer             :: nbnd,fullnkpts
      integer             :: ikpt,jspins
      integer             :: kpts
      logical             :: l_file
      integer             :: num_wann,num_kpts,num_nnmax,jspin
      integer             :: kspin,kkspin
      integer             :: num_wann2
      integer             :: i,j,k,m,info,r1,r2,r3,dummy1
      integer             :: dummy2,dummy3
      integer             :: counter,m1,m2
      integer             :: num_bands2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:,:),hwann(:,:,:,:)
      complex,allocatable :: hwannmix(:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:)
      complex,allocatable :: hrealsoc(:,:,:,:,:,:,:)
      complex,allocatable :: hwannsoc(:,:,:,:,:)
      complex,allocatable :: paulimat_opt(:,:,:,:)
      complex,allocatable ::    nablamat(:,:,:,:)
      complex,allocatable ::   nablamat2(:,:,:,:)
      complex,allocatable :: nablamatmix(:,:,:,:)
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
      integer,allocatable :: ndimwin(:,:)
      logical,allocatable :: lwindow(:,:,:)
      integer             :: chk_unit,nkp,ntmp,ierr
      character(len=33)   :: header
      character(len=20)   :: checkpoint
      real                :: tmp_latt(3,3), tmp_kpt_latt(3,nkpts)
      real                :: omega_invariant
      complex,allocatable :: u_matrix_opt(:,:,:,:)
      integer             :: num_bands,counter1,counter2
      logical             :: l_umdat
      real,allocatable    :: eigval2(:,:)
      real,allocatable    :: eigval_opt(:,:)
      real                :: scale,a,b
      character(len=2)    :: spinspin12(0:2)
      character(len=3)    :: spin12(2)
      character(len=6)    :: filename
      integer             :: jp,mp,kk
      integer             :: rvecind
      integer             :: spin1,spin2
      character(len=10)   :: torquename
      character(len=9)    :: torquename2
      character(len=15)   :: torquename3
      integer             :: at,nn,itype
      complex,allocatable :: perpmag(:,:,:,:)
      logical             :: l_yesss

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      call timestart("wann_torque_rs")

      tpi=2*pimach()

      jspins=jspins_in
      if(l_soc)jspins=1

      write(6,*)"nkpts=",nkpts

c*****************************************************
c     get num_bands and num_wann from the proj file
c*****************************************************
      do j=jspins,0,-1
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
          stop 'no proj/proj.1/proj.2' 
      endif
      read (203,*) num_wann,num_bands
      close (203)
      write(6,*)'According to proj there are ',num_bands,' bands'
      write(6,*)"and ",num_wann," wannier functions."

      allocate( u_matrix_opt(num_bands,num_wann,nkpts,jspins) )
      allocate( u_matrix(num_wann,num_wann,nkpts,jspins) )
      allocate( lwindow(num_bands,nkpts,jspins) )
      allocate( ndimwin(nkpts,jspins) )

      do jspin=1,jspins  !spin loop
c****************************************************************
c        read in chk
c****************************************************************
         num_kpts=nkpts

         call wann_read_umatrix2(
     >            nkpts,num_wann,num_bands,
     >            um_format,jspin,wan90version,
     <            have_disentangled,
     <            lwindow(:,:,jspin),
     <            ndimwin(:,jspin),
     <            u_matrix_opt(:,:,:,jspin),
     <            u_matrix(:,:,:,jspin) )

      enddo !jspin

      jspin=1

      num_bands2=jspins*num_bands
      num_wann2=jspins*num_wann

c************************************************
c     Read the file torque_...
c************************************************
      allocate( nablamat(3,num_bands2,num_bands2,nkpts) )
      at=0  
      do itype=1,ntype
       do nn=1,neq(itype)
        at=at+1
        write(torquename,fmt='("torque_",i3.3)')at
        if(l_unformatted)then  

       inquire(file=trim(torquename)//'_unf',
     &                  exist=l_yesss)
       if(.not.l_yesss)cycle

         open(304,file=trim(torquename)//'_unf',
     &              form='unformatted')
         read(304)nbnd,nbnd,fullnkpts
         if(nbnd.ne.num_bands2)stop'problem1'
         if(fullnkpts.ne.nkpts)stop'problem2'
         read(304)nablamat(1:3,1:nbnd,1:nbnd,1:fullnkpts)
        else   
         open(304,file=torquename,form='formatted')
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1,num_bands2
           do j=1,num_bands2
            do k=1,3  
             read(304,*)dummy1,dummy2,dummy3,a,b
             nablamat(k,j,i,nkp)=cmplx(a,b)
            enddo !k 
           enddo !j
          enddo !i
         enddo !nkp
        endif !l_unformatted 
        close(304)

c********************************************************************
c        Calculate matrix elements of MT spin current in the basis of
c        rotated Bloch functions.
c********************************************************************
         allocate( nablamat2(3,num_wann2,num_wann2,nkpts) )
         write(6,*)"calculate matrix elements of MT spin current 
     &   operator between wannier orbitals"

         if(have_disentangled) then       

          allocate( paulimat_opt(3,num_bands2,num_bands2,nkpts) )
          allocate( nablamatmix(3,num_bands2,num_wann2,nkpts) )

          do nkp=1,num_kpts
           counter1=0
           do m=1,num_bands2
            spin1=(m-1)/num_bands 
            if(lwindow(m-spin1*num_bands,nkp,spin1+1))then
             counter1=counter1+1  
             counter2=0
             do mp=1,num_bands2
              spin2=(mp-1)/num_bands 
              if(lwindow(mp-spin2*num_bands,nkp,spin2+1))then
               counter2=counter2+1
               do k=1,3
                paulimat_opt(k,counter2,counter1,nkp)=
     &          nablamat(k,mp,m,nkp)  
               enddo
              endif
             enddo !mp
            endif
           enddo !m 
          enddo

          nablamatmix=0.0  
          do spin1=0,jspins-1
           do spin2=0,jspins-1
            do nkp=1,num_kpts
             do j=1,num_wann
              do m=1,ndimwin(nkp,spin2+1)
               do mp=1,ndimwin(nkp,spin1+1)  
                do k=1,3  
      nablamatmix(k,mp+spin1*ndimwin(nkp,1),j+spin2*num_wann,nkp)=
     &nablamatmix(k,mp+spin1*ndimwin(nkp,1),j+spin2*num_wann,nkp)+ 
     &          paulimat_opt(k,mp+spin1*ndimwin(nkp,1),
     &                         m +spin2*ndimwin(nkp,1),nkp)*
     &                 u_matrix_opt(m,j,nkp,spin2+1)
                enddo !k 
               enddo !mp  
              enddo !m
             enddo !j
            enddo !nkp
           enddo !spin2
          enddo !spin1

          deallocate(paulimat_opt)

          nablamat2=0.0
          do spin1=0,jspins-1
           do spin2=0,jspins-1
            do nkp=1,num_kpts
             do j=1,num_wann
              do jp=1,num_wann  
               do mp=1,ndimwin(nkp,spin1+1)  
                do k=1,3  
                nablamat2(k,jp+spin1*num_wann,j+spin2*num_wann,nkp)=
     &          nablamat2(k,jp+spin1*num_wann,j+spin2*num_wann,nkp)+ 
     &                 conjg(u_matrix_opt(mp,jp,nkp,spin1+1))*
     & nablamatmix(k,mp+spin1*ndimwin(nkp,1),j+spin2*num_wann,nkp)
                enddo !k 
               enddo !mp  
              enddo !jp 
             enddo !j
            enddo !nkp
           enddo !spin2
          enddo !spin1
          deallocate( nablamatmix )
         else
          nablamat2(:,:,:,:)=nablamat(:,:,:,:)
         end if                    !have_disentangled

         allocate(hwann(3,num_wann2,num_wann2,num_kpts))
         hwann=cmplx(0.0,0.0)

         allocate(hwannmix(3,num_wann2,num_wann2,num_kpts))
         hwannmix=cmplx(0.0,0.0)

         do spin2=0,jspins-1
          do spin1=0,jspins-1
           do k=1,num_kpts
            do m=1,num_wann
             do j=1,num_wann
              do mp=1,num_wann
               do kk=1,3  
       hwannmix(kk,mp+spin1*num_wann,m+spin2*num_wann,k)=
     & hwannmix(kk,mp+spin1*num_wann,m+spin2*num_wann,k)+
     *              conjg(u_matrix(j,mp,k,spin1+1))*
     *                nablamat2(kk,j+spin1*num_wann,
     *                             m+spin2*num_wann,k)
               enddo !kk
              enddo !mp   
             enddo !j
            enddo !m
           enddo !k
          enddo !spin2
         enddo !spin1 


         do spin2=0,jspins-1
          do spin1=0,jspins-1
           do k=1,num_kpts
            do m=1,num_wann
             do i=1,num_wann
              do mp=1,num_wann
               do kk=1,3  
                hwann(kk,mp+spin1*num_wann,m+spin2*num_wann,k)=
     *          hwann(kk,mp+spin1*num_wann,m+spin2*num_wann,k)+
     *               hwannmix(kk,mp+spin1*num_wann,
     *                            i+spin2*num_wann,k)*
     *                    u_matrix(i,m,k,spin2+1)
               enddo !kk
              enddo !mp   
             enddo !i
            enddo !m
           enddo !k
          enddo !spin1
         enddo !spin2
         deallocate(hwannmix)
         deallocate(nablamat2)
c************************************************************
c        Calculate matrix elements in real space.
c***********************************************************      
         write(6,*)"calculate nabla-mat in rs"

         allocate(hreal(3,num_wann2,num_wann2,rvecnum))
         hreal=cmplx(0.0,0.0)

         if(l_ndegen)then
          do rvecind=1,rvecnum
            do k=1,nkpts  
              rdotk=tpi*(kpoints(1,k)*rvec(1,rvecind)+
     &                   kpoints(2,k)*rvec(2,rvecind)+
     &                   kpoints(3,k)*rvec(3,rvecind) )
              fac=cmplx(cos(rdotk),-sin(rdotk))/
     /         cmplx(real(ndegen(rvecind)),0.0)
              do i=1,num_wann2
               do j=1,num_wann2
                do kk=1,3
                 hreal(kk,j,i,rvecind)=
     &                   hreal(kk,j,i,rvecind)+
     &                   fac*hwann(kk,j,i,k)
                enddo !kk
               enddo !j
              enddo !i
            enddo !k
          enddo !rvecind
         else
          do rvecind=1,rvecnum
            do k=1,nkpts  
              rdotk=tpi*(kpoints(1,k)*rvec(1,rvecind)+
     &                   kpoints(2,k)*rvec(2,rvecind)+
     &                   kpoints(3,k)*rvec(3,rvecind) )
              fac=cmplx(cos(rdotk),-sin(rdotk))
              do i=1,num_wann2
               do j=1,num_wann2
                do kk=1,3
                 hreal(kk,j,i,rvecind)=
     &                   hreal(kk,j,i,rvecind)+
     &                   fac*hwann(kk,j,i,k)
                enddo !kk
               enddo !j
              enddo !i
            enddo !k
          enddo !rvecind
         endif
         hreal=hreal/cmplx(real(nkpts),0.0)
         deallocate(hwann)
         if(l_unformatted)then
           if(l_ndegen)then
             write(torquename3,fmt='("rstorndegen_",i3.3)')at
             open(321,file=trim(torquename3)//'_unf',             
     &            form='unformatted'
#ifdef CPP_INTEL
     &         ,convert='BIG_ENDIAN'
#endif
     &                   )             
           else ! l_ndegen
             write(torquename2,fmt='("rstor_",i3.3)')at
             open(321,file=trim(torquename2)//'_unf',
     &            form='unformatted'
#ifdef CPP_INTEL
     &         ,convert='BIG_ENDIAN'
#endif             
     &                   )
           endif ! l_ndegen 
           allocate(perpmag(num_wann2,num_wann2,3,rvecnum))           
           do rvecind=1,rvecnum
            do j=1,num_wann2
             do i=1,num_wann2
              do kk=1,3
                perpmag(i,j,kk,rvecind)=hreal(kk,i,j,rvecind)
              enddo !kk
             enddo !i
            enddo !j
           enddo !rvecind
           write(321)perpmag
           deallocate(perpmag)
         else   
           if(l_ndegen)then
             write(torquename3,fmt='("rstorndegen_",i3.3)')at
             open(321,file=torquename3)
           else
             write(torquename2,fmt='("rstor_",i3.3)')at
             open(321,file=torquename2)
           endif  
           do rvecind=1,rvecnum
             r3=rvec(3,rvecind)
             r2=rvec(2,rvecind)
             r1=rvec(1,rvecind)             
             do j=1,num_wann2
              do i=1,num_wann2
               do kk=1,3  
                write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &                 1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &                 r1,r2,r3,i,j,kk,hreal(kk,i,j,rvecind) 
               enddo !kk  
              enddo !j
             enddo !i
           enddo !rvecind  
         endif !l_unformatted  
         close(321)
         deallocate(hreal)
       enddo !nn  
      enddo !itype

      call timestop("wann_torque_rs")
      end subroutine wann_torque_rs
      end module m_wann_torque_rs
