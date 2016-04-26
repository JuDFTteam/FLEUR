c************************************************
c     generate equidistant k-point mesh as needed
c     for wannier functions
c     Frank Freimuth, September 2006
c************************************************
      program kpunktgen
      implicit none
      integer num,numkpt,dim
      integer c1,c2,c3,limit
      real ilen,shift
      real i1,i2,i3
      logical l_shift


      print*,"specify dimension"
      read(*,*)dim
      print*,"symmetric to origin?"
      read(*,*)l_shift
      if(dim==3)then
      print*,"create three dimensional k-point set"   
      print*,"Divisions per direction"
      read(*,*)num
      print*,"Divisions per direction: ",num
      numkpt=num**3
      print*,"Number of k-points: ",numkpt
      ilen=1.0/num
      print*,"length of intervals: ",ilen
c      open(100,file='kpts',form='formatted',status='unknown')
      open(200,file='w90kpts',form='formatted',status='unknown')
c      write(100,'(2x,i3,8x,f7.5)')numkpt,1.0
      write(200,*)numkpt
      limit=num-1
      if (l_shift) then
         shift=limit*ilen/2.0        
      endif

      do c1=0,limit
         do c2=0,limit
            do c3=0,limit
               i1=ilen*c1-shift
               i2=ilen*c2-shift
               i3=ilen*c3-shift
c               write(100,'(3x,f7.5,3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,i3,1.00000
c               write(200,'(3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,i3
               write(200,*)i1,i2,i3
            enddo
         enddo
      enddo
      close(100)
      close(200)

      elseif(dim==2)then
      print*,"create two dimensional k-point set"
         print*,"Divisions per direction"
      read(*,*)num
      print*,"Divisions per direction: ",num
      numkpt=num**2
      print*,"Number of k-points: ",numkpt
      ilen=1.0/num
      print*,"Length of intervals: ",ilen
c      open(100,file='kpts',form='formatted',status='unknown')
      open(200,file='w90kpts',form='formatted',status='unknown')
c      write(100,'(2x,i3,8x,f7.5,8x,1a)')numkpt,1.0,"F"
      write(200,*)numkpt
      limit=num-1
      if(l_shift)then
         shift=limit*ilen/2.0
      endif   
      do c1=0,limit
         do c2=0,limit
               i1=ilen*c1-shift
               i2=ilen*c2-shift
c               write(100,'(3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,1.00000
c               write(200,'(3x,f7.5,3x,f7.5,3x,f6.5)')i1,i2,0.00000
               write(200,*)i1,i2,0.00000
         enddo
      enddo
      close(100)
      close(200)


      elseif(dim==1)then
      print*,"create one dimensional k-point set"   
      print*,"Divisions per direction"
      read(*,*)num
      print*,"Divisions per direction: ",num
      numkpt=num
      print*,"Number of k-points: ",numkpt
      ilen=1.0/num
      print*,"length of intervals: ",ilen
c      open(100,file='kpts',form='formatted',status='unknown')
      open(200,file='w90kpts',form='formatted',status='unknown')
c      write(100,'(2x,i3,8x,f7.5)')numkpt,1.0
      write(200,*)numkpt
      limit=num-1
      if (l_shift) then
         shift=limit*ilen/2.0        
      endif

      do c1=0,limit

               i1=ilen*c1-shift
 
c               write(100,'(3x,f7.5,3x,f7.5,3x,f7.5,3x,f7.5)')i1,i2,i3,1.00000
c               write(200,'(3x,f7.5,3x,f7.5,3x,f7.5)')i1,0.0000,0.00000
               write(200,*)i1,0.0000,0.00000

      enddo
      close(100)
      close(200)



      endif
      end program kpunktgen
