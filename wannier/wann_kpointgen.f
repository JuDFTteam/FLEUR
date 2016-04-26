      module m_wann_kpointgen
      use m_juDFT
      private
      public::wann_kpointgen
      contains
      subroutine wann_kpointgen()
c***********************************************************
c     Generate the k-point file 'kpts' appropriate for the 
c     calculation of Wannier Functions with the Fleur-code.
c     In order to use the symmetry of the Brillouin zone,
c     use the program w90kpointgen instead.
c     Frank Freimuth, June 2007
c***********************************************************
      implicit none
      integer num,dim
      integer c1,c2,c3
      real i1,i2,i3
      logical l_shift
      integer Nx,Ny,Nz
      integer nnx,nny,nnz
      real ilen1,ilen2,ilen3
      integer limit1, limit2, limit3
      real scale,shift1,shift2,shift3
      
      print*,"specify dimension"
      read(*,*)dim
      print*,"symmetric to origin?"
      read(*,*)l_shift
      Nx=1;Ny=1;Nz=1
      if(dim==3)then
         print*,"Creates three dimensional k-point set."
      elseif(dim==2)then
         print*,"Create two-dimensional k-point set."
      elseif(dim==1)then
         print*,"Create one-dimensional k-point set."
      else
         CALL juDFT_error("unknown dimension",calledby ="wann_kpointgen")
      endif
      print*,"Specify the number of k-point steps"
      print*,"for each direction."
      if(.not.dim==2)then
         print*,"Nz="
         read(*,*)Nz
      endif
      if(.not.dim==1)then
         print*,"Ny="
         read(*,*)Ny
         print*,"Nx="
         read(*,*)Nx
      endif   
      num=Nx*Ny*Nz
      print*,"Number of k-points: ",num
      call findkgt(Nx,nnx)
      call findkgt(ny,nny)
      call findkgt(nz,nnz)
      scale=nnx
      if(nny.ne.nnx)scale=scale*nny
      if(nnz.ne.nny .and. nnz.ne.nnx)
     &      scale=scale*nnz
      print*,"scale=",scale
      ilen1=1.0/nx; ilen2=1.0/ny; ilen3=1.0/nz
      limit1=nx-1; limit2=ny-1; limit3=nz-1   
      if(l_shift)then
         shift1=limit1*ilen1/2.0
         shift2=limit2*ilen2/2.0
         shift3=limit3*ilen3/2.0
      endif
      open(100,file='kpts',form='formatted',status='unknown')
      if(dim.ne.2)then !not the film-format of kpts
       write(100,'(i5,f20.10)')num,scale
       do c1=0,limit1
         do c2=0,limit2
            do c3=0,limit3
               i1=(ilen1*c1-shift1)*scale
               i2=(ilen2*c2-shift2)*scale
               i3=(ilen3*c3-shift3)*scale
               write(100,'(4f10.5)')i1,i2,i3,1.0
            enddo
         enddo
       enddo
      else !film-format
         write(100,'(i5,f20.10,3x,l1)')num,scale,.false.
         do c1=0,limit1
           do c2=0,limit2
               i1=(ilen1*c1-shift1)*scale
               i2=(ilen2*c2-shift2)*scale
               i3=(ilen3*c3-shift3)*scale
               write(100,'(3f10.5)')i1,i2,1.0
           enddo
         enddo
      endif
      close(100)
      end subroutine wann_kpointgen

      subroutine findkgt(nu,sc)
      implicit none
      integer,intent(out)::sc
      integer,intent(in)::nu
      integer k,nnu
      nnu=nu
      IF(nnu==0)  CALL juDFT_error("nnu.eq.0",calledby ="wann_kpointgen")
      do k=1,3
         if(nnu.eq.1)exit
         if(mod(nnu,5).ne.0)exit
         nnu=nnu/5
      enddo
      do k=1,3
         if(nnu.eq.1)exit
         if(mod(nnu,2).ne.0)exit
         nnu=nnu/2
      enddo
      sc=nnu
      end subroutine findkgt

      end module m_wann_kpointgen
