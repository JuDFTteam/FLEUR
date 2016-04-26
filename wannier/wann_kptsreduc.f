      module m_wann_kptsreduc
      use m_juDFT
      contains
      subroutine wann_kptsreduc(
     >               nop,mrot,bmat,tau,film,
     >               l_onedimens,l_nocosoc)
c*****************************************************************
c     Apply the symmetries to reduce the number of k-points.
c     Frank Freimuth
c*****************************************************************
      implicit none
      integer, intent(in) :: nop
      logical,intent(in)  :: film
      real, intent(in)    :: bmat(3,3)
      real,intent(in)     :: tau(3,nop)
      integer, intent(in) :: mrot(3,3,nop)
      logical,intent(in)  :: l_onedimens
      logical,intent(in)  :: l_nocosoc

      real,allocatable    :: weight(:)
      integer             :: reduznumk
      integer,allocatable :: mapk(:),mapkoper(:)
      integer             :: ikpt,nmat,k,j
      real                :: bbkpt(3)
      real                :: bkpt(3),wk,brot(3)
      real,allocatable    :: kpoints(:,:)
      integer             :: oper,iter
      logical             :: l_file
      real                :: scale
      integer             :: nkpts
      integer             :: sumweights
      integer,allocatable :: irreduc(:)
      integer,allocatable :: shiftkpt(:,:)
      real                :: bbmat(3,3)
      real,allocatable    :: kptslen(:)
      logical             :: l_onlysymor

      inquire(file='onlysymor',exist=l_onlysymor)

      bbmat=matmul(bmat,transpose(bmat))
      write(6,*) "Apply the symmetries to w90kpts"
c**********************************************************
c     read in kpoints from w90kpts file
c**********************************************************
      l_file=.false.
      inquire(file='w90kpts',exist=l_file)
      IF(.NOT.l_file) CALL juDFT_error("where is w90kpts?",calledby
     +     ="wann_kptsreduc")
      open(987,file='w90kpts',status='old',form='formatted')
      read(987,*)nkpts,scale
      print*,"nkpts=",nkpts
      allocate(mapk(nkpts),mapkoper(nkpts),kpoints(3,nkpts))
      allocate(weight(nkpts),irreduc(nkpts))
      allocate(shiftkpt(3,nkpts))
      do iter=1,nkpts
         read(987,*)kpoints(:,iter)
      enddo
      close(987)
      write(6,*) "kpoints read from w90kpts:"
      do iter=1,nkpts
         write(6,'(3f10.5)')kpoints(:,iter)/scale
      enddo

c*********************************************************
c     determine lengths of kpoints
c*********************************************************
c      allocate( kptslen(nkpts) )
c      do iter=1,nkpts
c         kptslen(iter)=dot_product( axmintx(kpoints(:,iter)),
c     &                matmul(bbmat,axmintx(kpoints(:,iter))) )
c      enddo

c***********************************************
c            determine irreducible part
c***********************************************
      weight(:)=1.0
      mapk(:)=0
      reduznumk=0
      mapkoper(:)=0
      irreduc(:)=0
      shiftkpt(:,:)=0
      do ikpt=1,nkpts
         if(mapk(ikpt).ne.0) cycle
         reduznumk=reduznumk+1
         irreduc(ikpt)=reduznumk
         bkpt(:)=kpoints(:,ikpt)         
         do oper=1,nop
           if (all(abs(tau(:,oper)).lt.1e-10).or..not.l_onlysymor)then
             brot(:)=0.0
             do k=1,3
               brot(:)=brot(:)+mrot(k,:,oper)*bkpt(k)
             enddo
             do j=ikpt+1,nkpts
              if(mapk(j).ne.0)cycle
c              if(abs(kptslen(j)-kptslen(ikpt)).gt.1.e-8)cycle
              if( all( axmintx( (brot(:)-kpoints(:,j))/scale ) 
     &                            .lt.1e-6) ) then
                    weight(ikpt)=weight(ikpt)+1.0
                    mapk(j)=ikpt
                    mapkoper(j)=oper
                    shiftkpt(:,j)=
     &                 nint( (brot(:)-kpoints(:,j))/scale )
              endif
             enddo
           endif
         enddo
         if(.not.l_nocosoc)then
          do oper=1,nop
           if (all(abs(tau(:,oper)).lt.1e-10).or..not.l_onlysymor)then
             brot(:)=0.0
             do k=1,3
               brot(:)=brot(:)-mrot(k,:,oper)*bkpt(k)
             enddo
             do j=ikpt+1,nkpts
              if(mapk(j).ne.0)cycle
c              if(abs(kptslen(j)-kptslen(ikpt)).gt.1.e-8)cycle
              if( all( axmintx( (brot(:)-kpoints(:,j))/scale) 
     &                            .lt.1e-6) ) then
                    weight(ikpt)=weight(ikpt)+1.0
                    mapk(j)=ikpt
                    mapkoper(j)=-oper
                    shiftkpt(:,j)=nint( (brot(:)-kpoints(:,j))/scale )
              endif
             enddo
           endif
          enddo
         endif 
      enddo   

c****************************************************
c         write results to files
c         w90kpts: whole Brillouin zone (for w90)
c         kpts: irreducible part (for fleur)
c         kptsmap: mapping from w90kpts to kpts
c****************************************************
      open(117,file='kptsmap',form='formatted',recl=1000)
      do ikpt=1,nkpts
         if(mapk(ikpt)==0)then
            write(117,*)ikpt,irreduc(ikpt),1,0,0,0
         else
            write(117,*)ikpt,irreduc(mapk(ikpt)),mapkoper(ikpt),
     &                  shiftkpt(:,ikpt)
         endif  
      enddo   
      close(117)
      open(119,file='kpts',form='formatted')
      if (film.and..not.l_onedimens)then
        write(119,'(i5,f20.10,3x,l1)')reduznumk,scale,.false.         
      else
        write(119,'(i5,f20.10)')reduznumk,scale
      endif   
      sumweights=0
      do ikpt=1,nkpts
         bkpt(:)=kpoints(:,ikpt)
         if (mapk(ikpt)==0)then
            sumweights=sumweights+weight(ikpt)
            write(6,*)"ikpt=",ikpt
            write(6,*)"irreducible"
            write(6,fmt='(a10,3f9.6)')"internal: ",bkpt(:)/scale

            if (film.and..not.l_onedimens)then
               write(119,'(3f10.5)')bkpt(1:2),weight(ikpt)
            else
               write(119,'(4f10.5)')bkpt(:),weight(ikpt)
            endif   
      
         elseif(mapkoper(ikpt).gt.0)then
            write(6,*)"ikpt=",ikpt
            write(6,*)"reducible"
            write(6,*)"map=",mapk(ikpt)
            brot(:)=0.0
            do k=1,3
               brot(:)=brot(:)+
     + mrot(k,:,mapkoper(ikpt))*kpoints(k,mapk(ikpt))
            enddo   
            write(6,'(a19,3f9.6)')"rotated internal: ",brot(:)/scale
         elseif(mapkoper(ikpt).lt.0)then
            write(6,*)"ikpt=",ikpt
            write(6,*)"reducible"
            write(6,*)"map=",mapk(ikpt)
            brot(:)=0.0
            do k=1,3
               brot(:)=brot(:)-
     + mrot(k,:,-mapkoper(ikpt))*kpoints(k,mapk(ikpt))
            enddo   
            write(6,'(a19,3f9.6)')"rotated internal: ",brot(:)/scale
         endif   
      enddo   
      close(119)

      write(6,*)"reduznumk=",reduznumk     
      write(6,*)"nkpts=",nkpts
      write(6,*)"sumweights=",sumweights
      
      IF(sumweights/=nkpts) CALL juDFT_error
     +     ("sum of weights differs from nkpts",calledby
     +     ="wann_kptsreduc")
      deallocate(kpoints,mapk,mapkoper,weight,irreduc)

      contains
      real elemental function axmintx(x)
      implicit none
      real,intent(in) :: x
      axmintx = abs(x-nint(x))
      end function     

      end subroutine wann_kptsreduc
      end module m_wann_kptsreduc
