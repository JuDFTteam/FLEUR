      module m_wann_pauli_rs
      USE m_juDFT
      contains 
      subroutine wann_pauli_rs(
     >          rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >          l_soc,band_min,band_max,neigd,
     >          l_socmmn0,l_ndegen,ndegen,wan90version,
     >          l_unformatted)
c*************************************************
c     Calculate the matrix elements of the 
c     Pauli matrix in real space from the
c     files WF1.chk (and WF1_um.dat) produced
c     by wannier90.
c     FF, January 2009
c*************************************************
      use m_constants, only:pimach
      use m_wann_read_umatrix

      implicit none
      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)
      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical,intent (in) :: l_bzsym,l_soc
      logical,intent(in)  :: film
      integer,intent(in)  :: band_min(2),band_max(2),neigd
      logical, intent(in) :: l_socmmn0
      logical, intent(in) :: l_ndegen
      integer, intent(in) :: ndegen(:)	
      integer, intent(in) :: wan90version
      logical, intent(in) :: l_unformatted

      integer             :: ikpt,jspins
      integer             :: kpts
      logical             :: l_file
      integer             :: num_wann,num_kpts,num_nnmax,jspin
      integer             :: kspin,kkspin
      integer             :: wann_shift,num_wann2
      integer             :: i,j,k,m,info,r1,r2,r3,dummy1
      integer             :: dummy2,dummy3
      integer             :: hopmin,hopmax,counter,m1,m2
      integer             :: num_bands2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:,:),hwann(:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:)
      complex,allocatable :: hrealsoc(:,:,:,:,:,:,:)
      complex,allocatable :: hwannsoc(:,:,:,:,:)
      complex,allocatable :: paulimat(:,:,:,:)
      complex,allocatable :: paulimat_opt(:,:,:,:)
      complex,allocatable :: paulimat2(:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      complex             :: value(4)
      logical             :: um_format
      logical             :: repro_eig
      logical             :: l_chk,l_proj,l_onedimens
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
      complex,parameter   :: ci=(0.0,1.0)
      integer             :: dir,rvecind
      integer             :: spin1,spin2
      complex,allocatable :: amn(:,:,:)
      complex,allocatable :: paulimatmix(:,:,:,:)
      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/
      integer             :: nbnd,nwfs,fullnkpts
      complex,allocatable :: hwannmix(:,:,:,:)
      complex,allocatable :: perpmag(:,:,:,:)

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
             CALL fleur_err("no proj/proj.1/proj.2",calledby
     +           ="wann_pauli_rs")
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
     <            u_matrix(:,:,:,jspin))

 
      enddo   

      jspin=1

      num_bands2=jspins*num_bands
      num_wann2=jspins*num_wann

      allocate( paulimat(3,num_bands2,num_bands2,nkpts) )
      paulimat=0.0

c****************************************************
c        Read the file "WF1.socmmn0".
c**************************************************** 
      if(l_unformatted)then
       if(l_soc)then
         open(304,file='WF1.socmmn0_unf',form='unformatted')
       else
         open(304,file='WF1.mmn0_unf',form='unformatted')
       endif
       read(304)nbnd,fullnkpts,nwfs
       allocate(amn(nbnd,nwfs,fullnkpts))
       read(304)amn(1:nbnd,1:nwfs,1:fullnkpts)
       do nkp=1,num_kpts
        do i=1,num_bands
         do j=1,num_bands
          paulimat(3,j,i,nkp)=amn(i,j,nkp)
         enddo !j
        enddo !i
       enddo !nkp
       close(304)
       deallocate(amn)
      else !l_unformatted
         if(l_soc)then
            open(304,file='WF1.socmmn0',form='formatted')
         else
            open(304,file='WF1.mmn0',form='formatted')
         endif
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1,num_bands
           do j=1,num_bands  
             read(304,*)dummy1,dummy2,dummy3,a,b
             paulimat(3,i,j,nkp)=cmplx(a,b)
           enddo !j
          enddo !i
         enddo !nkp
         close(304)
      endif   
c****************************************************
c        Read the file "WF2.socmmn0".
c****************************************************
      if(l_unformatted)then
       if(l_soc)then
         open(304,file='WF2.socmmn0_unf',form='unformatted')
       else
         open(304,file='WF2.mmn0_unf',form='unformatted')
       endif
       read(304)nbnd,fullnkpts,nwfs
       allocate(amn(nbnd,nwfs,fullnkpts))
       read(304)amn(1:nbnd,1:nwfs,1:fullnkpts)
       do nkp=1,num_kpts
        do i=1+(jspins-1)*num_bands,num_bands2
         do j=1+(jspins-1)*num_bands,num_bands2
           paulimat(3,i,j,nkp)=
     &      paulimat(3,i,j,nkp)-
     &      amn(j-(jspins-1)*num_bands,
     &          i-(jspins-1)*num_bands,nkp)
         enddo !j
        enddo !i
       enddo !nkp
       close(304)
       deallocate(amn)
      else !l_unformatted
         if(l_soc)then
            open(304,file='WF2.socmmn0',form='formatted')
         else
            open(304,file='WF2.mmn0',form='formatted')            
         endif
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1+(jspins-1)*num_bands,num_bands2
           do j=1+(jspins-1)*num_bands,num_bands2
             read(304,*)dummy1,dummy2,dummy3,a,b
             paulimat(3,i,j,nkp)=
     &          paulimat(3,i,j,nkp)-cmplx(a,b)
           enddo !j
          enddo !i
         enddo !nkp
         close(304)
      endif !l_unformatted   

c****************************************************
c        Read the file "updown.mmn0".
c****************************************************
      if(l_unformatted)then
         open(304,file='updown.mmn0_unf',form='unformatted')
         read(304)nbnd,fullnkpts,nwfs
         allocate(amn(nbnd,nwfs,fullnkpts))
         read(304)amn(1:nbnd,1:nwfs,1:fullnkpts)
         close(304)
         deallocate(amn)
         do nkp=1,fullnkpts
          do i=1,num_bands
           do j=1,num_bands
             paulimat(1,j,i+num_bands*(jspins-1),nkp)=
     &             conjg(amn(j,i,nkp))
           enddo !j
          enddo !i
          paulimat(2,:,:,nkp)= -ci*paulimat(1,:,:,nkp)
          paulimat(1,:,:,nkp)=paulimat(1,:,:,nkp)+
     &       transpose(conjg( paulimat(1,:,:,nkp) ))
          paulimat(2,:,:,nkp)=paulimat(2,:,:,nkp)+
     &       transpose(conjg( paulimat(2,:,:,nkp) ))
         enddo !nkp

      else !l_unformatted   
         open(304,file='updown.mmn0',form='formatted')
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1,num_bands
           do j=1,num_bands
             read(304,*)dummy1,dummy2,dummy3,a,b
             paulimat(1,j,i+num_bands*(jspins-1),nkp)=cmplx(a,-b)
           enddo !j
          enddo !i
          paulimat(2,:,:,nkp)= -ci*paulimat(1,:,:,nkp)
          paulimat(1,:,:,nkp)=paulimat(1,:,:,nkp)+
     &        transpose(conjg( paulimat(1,:,:,nkp) ))
          paulimat(2,:,:,nkp)=paulimat(2,:,:,nkp)+
     &        transpose(conjg( paulimat(2,:,:,nkp) ))
         enddo !nkp
         close(304)
      endif !l_unformatted   

c****************************************************************
c        Calculate matrix elements of Pauli in the basis of
c        rotated Bloch functions.
c****************************************************************
         allocate( paulimat2(3,num_wann2,num_wann2,nkpts) )
         write(6,*)"calculate matrix elements of spin operator
     &   between wannier orbitals"

         if(have_disentangled) then       

          allocate( paulimat_opt(3,num_bands2,num_bands2,nkpts) )
          allocate( paulimatmix(3,num_wann2,num_bands2,nkpts))

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
     &          paulimat(k,mp,m,nkp)  
               enddo
              endif
             enddo !mp
            endif
           enddo !m 
          enddo

          paulimatmix=0.0  
          do spin1=0,jspins-1
            do spin2=0,jspins-1
             do nkp=1,num_kpts
              do jp=1,num_wann  
               do m=1,ndimwin(nkp,spin2+1)
                do mp=1,ndimwin(nkp,spin1+1)
                 do k=1,3  

          paulimatmix(k,jp+spin1*num_wann,m+spin2*ndimwin(nkp,1),nkp)=
     =    paulimatmix(k,jp+spin1*num_wann,m+spin2*ndimwin(nkp,1),nkp)+
     &                 conjg(u_matrix_opt(mp,jp,nkp,spin1+1))*
     &                        paulimat_opt(k,mp+spin1*ndimwin(nkp,1),
     &                                   m+spin2*ndimwin(nkp,1),nkp)
                 enddo !k 
                enddo !mp   
               enddo !m
              enddo !jp 
             enddo !nkp
            enddo !spin2
          enddo !spin1

          paulimat2=0.0  
          do spin1=0,jspins-1
           do spin2=0,jspins-1
            do nkp=1,num_kpts
             do j=1,num_wann
              do jp=1,num_wann  
               do m=1,ndimwin(nkp,spin2+1)

                 do k=1,3  


                paulimat2(k,jp+spin1*num_wann,j+spin2*num_wann,nkp)=
     =          paulimat2(k,jp+spin1*num_wann,j+spin2*num_wann,nkp)+ 
     &                        paulimatmix(k,jp+spin1*num_wann,

     &                                   m+spin2*ndimwin(nkp,1),nkp)*
     &                       u_matrix_opt(m,j,nkp,spin2+1)
                 enddo !k 
               enddo !m
              enddo !jp 
             enddo !j
            enddo !nkp
           enddo !spin2
          enddo !spin1   


        deallocate(paulimatmix)
          deallocate(paulimat_opt)

         else
          paulimat2=paulimat
         end if !have_disentangled

         allocate(hwann(3,num_wann2,num_wann2,num_kpts))
         hwann=cmplx(0.0,0.0)
         wann_shift=0

      allocate(hwannmix(3,num_wann2,num_wann2,num_kpts))
      hwannmix=cmplx(0.0,0.0)

      do spin1=0,jspins-1
       do spin2=0,jspins-1
        do k=1,num_kpts
         do i=1,num_wann
          do mp=1,num_wann
           do j=1,num_wann

             do kk=1,3  
              hwannmix(kk,mp+spin1*num_wann,i+spin2*num_wann,k)=
     =         hwannmix(kk,mp+spin1*num_wann,i+spin2*num_wann,k)+
     *           conjg(u_matrix(j,mp,k,spin1+1))*
     *               paulimat2(kk,j+spin1*num_wann,
     *                            i+spin2*num_wann,k)

             enddo !kk
           enddo !j
          enddo !mp

         enddo !i     
        enddo !k
       enddo !spin2
      enddo !spin1 

      do spin1=0,jspins-1
       do spin2=0,jspins-1
        do k=1,num_kpts
         do m=1,num_wann
          do i=1,num_wann
           do mp=1,num_wann
            do kk=1,3  
              hwann(kk,mp+spin1*num_wann,m+spin2*num_wann,k)=
     =         hwann(kk,mp+spin1*num_wann,m+spin2*num_wann,k)+
     *               hwannmix(kk,mp+spin1*num_wann,
     *                            i+spin2*num_wann,k)*
     *                 u_matrix(i,m,k,spin2+1)
            enddo !kk
           enddo !mp
          enddo !i     
         enddo !m
        enddo !k
       enddo !spin2
      enddo !spin1

      deallocate(hwannmix)



c************************************************************
c        Calculate matrix elements in real space.
c***********************************************************      
         write(6,*)"calculate pauli-mat in rs"
         allocate(hreal(3,num_wann2,num_wann2,rvecnum))
         hreal = cmplx(0.0,0.0)
         if(l_ndegen)then
          do rvecind=1,rvecnum
           do k=1,nkpts  
            rdotk=tpi*( kpoints(1,k)*rvec(1,rvecind)+
     &                  kpoints(2,k)*rvec(2,rvecind)+
     &                  kpoints(3,k)*rvec(3,rvecind) )
            fac=cmplx(cos(rdotk),-sin(rdotk))/real(ndegen(rvecind))
            do m2=1,num_wann2
             do m1=1,num_wann2
              do dir=1,3  
               hreal(dir,m1,m2,rvecind)=
     &                   hreal(dir,m1,m2,rvecind)+
     &                   fac*hwann(dir,m1,m2,k)
              enddo !dir 
             enddo !m1  
            enddo !m2  
           enddo !k
          enddo !rvecind
         else
          do rvecind=1,rvecnum
           do k=1,nkpts  
            rdotk=tpi*( kpoints(1,k)*rvec(1,rvecind)+
     &                  kpoints(2,k)*rvec(2,rvecind)+
     &                  kpoints(3,k)*rvec(3,rvecind) )
            fac=cmplx(cos(rdotk),-sin(rdotk))
            do m2=1,num_wann2
             do m1=1,num_wann2
              do dir=1,3  
               hreal(dir,m1,m2,rvecind)=
     &                   hreal(dir,m1,m2,rvecind)+
     &                   fac*hwann(dir,m1,m2,k)
              enddo !dir 
             enddo !m1  
            enddo !m2  
           enddo !k
          enddo !rvecind
         endif !l_ndegen 
         hreal=hreal/cmplx(real(nkpts),0.0)


       if(l_unformatted)then
	 allocate(perpmag(num_wann2,num_wann2,3,rvecnum))

         if(l_ndegen)then
           open(321,file='rspauli_ndegen_unf'//spinspin12(1),
     &       form='unformatted'
#ifdef CPP_INTEL
     &       ,convert='BIG_ENDIAN'
#endif
     &             )
         else
           open(321,file='rspauli_unf'//spinspin12(1),
     &       form='unformatted'
#ifdef CPP_INTEL
     &       ,convert='BIG_ENDIAN'
#endif
     &             )
         endif

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
         close(321)

       else !l_unformatted 

         if(l_ndegen)then
           open(321,file='rspauli_ndegen'//spinspin12(1),
     &            form='formatted')
         else
           open(321,file='rspauli'//spinspin12(1),
     &            form='formatted')
         endif
         do rvecind=1,rvecnum
           r3=rvec(3,rvecind)
           r2=rvec(2,rvecind)
           r1=rvec(1,rvecind)
           do j=1,num_wann2
            do i=1,num_wann2
             do kk=1,3   
              write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &           1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,kk,hreal(kk,i,j,rvecind) 
             enddo !kk 
            enddo !i
           enddo !j   
         enddo !rvecnum        

         close(321)

      endif    

      deallocate(lwindow,u_matrix_opt,ndimwin)
      deallocate(u_matrix,hwann,hreal)

      end subroutine wann_pauli_rs
      end module m_wann_pauli_rs
