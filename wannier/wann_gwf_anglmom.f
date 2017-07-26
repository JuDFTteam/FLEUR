c*************************c
c  routine to set up the  c
c  composite matrix anglmom c
c*************************c
      module m_wann_gwf_anglmom
      USE m_fleurenv
      implicit none
      contains

      subroutine wann_gwf_anglmom(nkpts,nqpts,l_unformatted)
      use m_wann_gwf_tools

      implicit none
      
      ! input parameters
      logical,intent(in) :: l_unformatted
      integer,intent(in) :: nkpts,nqpts

      integer :: nwfs,nbnd,nkqpts
      integer :: d,i,j,k,q,kq,t1,t2,t3,t4
      real :: tempr,tempi 

      complex,allocatable :: anglmom(:,:,:,:)
      character(len=12) :: fending

      nkqpts = nkpts*nqpts

      ! get number of bands and wfs from proj
      open(405,file='proj',status='old')
      read(405,*)nwfs,nbnd
      close(405)
      write(*,*)'nbnd=',nbnd
      write(*,*)'nwfs=',nwfs

      allocate(anglmom(3,nbnd,nbnd,nkqpts))

      ! read in separate files
      do q=1,nqpts
       WRITE(fending,'("_",i4.4)')q
       open(405,file='WF1'//trim(fending)//'.anglmom')
       read(405,*) !header
       read(405,*) !nbnd and nkpts
       do k=1,nkpts
        kq=get_index_kq(k,q,nkpts)
        do i=1,nbnd
         do j=1,nbnd
          do d=1,3
           read(405,*)t1,t2,t3,t4,tempr,tempi
           anglmom(d,j,i,kq) = cmplx(tempr,tempi)
          enddo
         enddo
        enddo
       enddo
       close(405,status='delete')    
      enddo 

      ! write file either unformatted or formatted
      if(l_unformatted) then
       open(405,file='WF1_gwf.anglmom',form='unformatted')
       write(405)anglmom
       close(405)
      else
       open(405,file='WF1_gwf.anglmom')
       write(405,*)'Matrix elements of angular momentum'
       write(405,'(3i5)')nbnd,nbnd,nkqpts
       do kq=1,nkqpts
        do i=1,nbnd
         do j=1,nbnd
          do d=1,3
           write(405,'(4i5,3x,2f18.12)')d,j,i,kq,
     >                                  anglmom(d,j,i,kq)
          enddo
         enddo
        enddo
       enddo
       close(405)
      endif

      deallocate(anglmom)

      end subroutine wann_gwf_anglmom
      end module m_wann_gwf_anglmom
