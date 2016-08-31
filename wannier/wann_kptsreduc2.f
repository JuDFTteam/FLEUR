!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_kptsreduc2
      use m_juDFT
      contains
      subroutine wann_kptsreduc2(
     >               mhp,
     >               nop,mrot,bmat,tau,film,
     >               l_onedimens,l_nocosoc)
c*****************************************************************
c     Apply the symmetries to reduce the number of k-points.
c     Frank Freimuth
c*****************************************************************
      implicit none
      integer, intent(in) :: mhp(3)
      integer, intent(in) :: nop
      logical,intent(in)  :: film
      real, intent(in)    :: bmat(3,3)
      real,intent(in)     :: tau(3,nop)
      integer, intent(in) :: mrot(3,3,nop)
      logical,intent(in)  :: l_onedimens
      logical,intent(in)  :: l_nocosoc

      real,allocatable    :: weight(:)
      integer             :: reduznumk,ikpt2
      integer,allocatable :: mapk(:),mapkoper(:)
      integer             :: ikpt,nmat,k,j
      integer             :: k1,k2,k3
      real                :: bbkpt(3)
      real                :: wk
      integer             :: brot(3)
      integer             :: bkpt(3),bkpt2(3)
      real                :: kpoint(3)
      integer             :: oper,iter
      logical             :: l_file
      real                :: scale
      integer             :: nkpts
      integer             :: sumweights
      integer,allocatable :: irreduc(:)
      integer,allocatable :: shiftkpt(:,:)
      real                :: bbmat(3,3)
      real,allocatable    :: kptslen(:)
      integer             :: i,kk1,kk2,kk3


      bbmat=matmul(bmat,transpose(bmat))
      write(6,*) "Apply the symmetries to w90kpts"
c**************************************************************
c     The array 'mhp' specifies the Monkhorst-Pack mesh.
c**************************************************************
      nkpts = mhp(1) * mhp(2) * mhp(3)

      allocate(mapk(nkpts),mapkoper(nkpts))
      allocate(weight(nkpts),irreduc(nkpts))

c*********************************************************
c     determine lengths of kpoints
c*********************************************************
      allocate( kptslen(nkpts) )
      iter=0
      do k3=0,mhp(3)-1
       do k2=0,mhp(2)-1
        do k1=0,mhp(1)-1  
         iter=iter+1
         kpoint(1)=real(k1)/real(mhp(1))
         kpoint(2)=real(k2)/real(mhp(2))
         kpoint(3)=real(k3)/real(mhp(3))
         kptslen(iter)=dot_product( axmintx(kpoint(:)),
     &                matmul(bbmat,axmintx(kpoint(:))) )
        enddo !k1 
       enddo !k2
      enddo !k3

c***********************************************
c     determine irreducible part
c***********************************************
      weight(:)=1.0
      mapk(:)=0
      reduznumk=0
      mapkoper(:)=0
      irreduc(:)=0
      shiftkpt(:,:)=0
      ikpt=0
      do k3=0,mhp(3)-1
       do k2=0,mhp(2)-1
        do k1=0,mhp(1)-1
         ikpt=ikpt+1
         if(mapk(ikpt).ne.0) cycle
         reduznumk=reduznumk+1
         irreduc(ikpt)=reduznumk
c         bkpt(:)=kpoints(:,ikpt)         
         bkpt(1)=k1*mhp(2)*mhp(3)
         bkpt(2)=k2*mhp(1)*mhp(3)
         bkpt(3)=k3*mhp(1)*mhp(2)
         do oper=1,nop
          do i=1,3
           brot(i)=0.0
           do k=1,3
            brot(i)=brot(i)+mrot(k,i,oper)*bkpt(k)
           enddo
          enddo
          ikpt2=0
          do kk3=0,mhp(3)-1
           do kk2=0,mhp(2)-1
            do kk1=0,mhp(1)-1
             ikpt2=ikpt2+1  
             if(ikpt2.le.ikpt)cycle  
             bkpt2(1)=kk1*mhp(2)*mhp(3)
             bkpt2(2)=kk2*mhp(1)*mhp(3)
             bkpt2(3)=kk3*mhp(1)*mhp(2)
             if(mapk(ikpt2).ne.0)cycle
             if(abs(kptslen(ikpt2)-kptslen(ikpt)).gt.1.e-8)cycle
             kpoint(1)=real(bkpt(1)-bkpt2(1))/real(mhp(1)*mhp(2)*mhp(3))
             kpoint(2)=real(bkpt(2)-bkpt2(2))/real(mhp(1)*mhp(2)*mhp(3))
             kpoint(3)=real(bkpt(3)-bkpt2(3))/real(mhp(1)*mhp(2)*mhp(3))
             if( all(   axmintx(kpoint(:)) .lt.1e-6     ) ) then
                    weight(ikpt)=weight(ikpt)+1.0
                    mapk(ikpt2)=ikpt
                    mapkoper(ikpt2)=oper
             endif
            enddo !kk1
           enddo !kk2
          enddo !kk3
         enddo !oper
         if(.not.l_nocosoc)then
          do oper=1,nop
           do i=1,3
            brot(i)=0.0
            do k=1,3
             brot(i)=brot(i)+mrot(k,i,oper)*bkpt(k)
            enddo
           enddo
           ikpt2=0
           do kk3=0,mhp(3)-1
            do kk2=0,mhp(2)-1
             do kk1=0,mhp(1)-1
              ikpt2=ikpt2+1  
              if(ikpt2.le.ikpt)cycle  
              bkpt2(1)=kk1*mhp(2)*mhp(3)
              bkpt2(2)=kk2*mhp(1)*mhp(3)
              bkpt2(3)=kk3*mhp(1)*mhp(2)
              if(mapk(ikpt2).ne.0)cycle
              if(abs(kptslen(ikpt2)-kptslen(ikpt)).gt.1.e-8)cycle
              kpoint(1)=
     &           real(bkpt(1)-bkpt2(1))/real(mhp(1)*mhp(2)*mhp(3))
              kpoint(2)=
     &           real(bkpt(2)-bkpt2(2))/real(mhp(1)*mhp(2)*mhp(3))
              kpoint(3)=
     &           real(bkpt(3)-bkpt2(3))/real(mhp(1)*mhp(2)*mhp(3))
              if( all(  axmintx(kpoint(:)).lt.1e-6  ) ) then
                    weight(ikpt)=weight(ikpt)+1.0
                    mapk(ikpt2)=ikpt
                    mapkoper(ikpt2)=oper
              endif
             enddo !kk1
            enddo !kk2
           enddo !kk3
          enddo !oper
         endif !l_nocosoc 
        enddo !k1
       enddo !k2 
      enddo !k3  

c****************************************************
c     write results to files
c     w90kpts: whole Brillouin zone (for w90)
c     kpts: irreducible part (for fleur)
c     kptsmap: mapping from w90kpts to kpts
c****************************************************
c      open(117,file='kptsmap',form='formatted',recl=1000)
c      do ikpt=1,nkpts
c         if(mapk(ikpt)==0)then
c            write(117,*)ikpt,irreduc(ikpt),1,0,0,0
c         else
c            write(117,*)ikpt,irreduc(mapk(ikpt)),mapkoper(ikpt),
c     &                  shiftkpt(:,ikpt)
c         endif  
c      enddo   
c      close(117)

      scale=1.000

      open(119,file='kpts',form='formatted')
      if (film.and..not.l_onedimens)then
        write(119,'(i5,f20.10,3x,l1)')reduznumk,scale,.false.         
      else
        write(119,'(i5,f20.10)')reduznumk,scale
      endif   
      sumweights=0

      ikpt=0
      do k3=0,mhp(3)-1
       do k2=0,mhp(2)-1
        do k1=0,mhp(1)-1  
         ikpt=ikpt+1
         kpoint(1)=real(k1)/real(mhp(1))
         kpoint(2)=real(k2)/real(mhp(2))
         kpoint(3)=real(k3)/real(mhp(3))


         if (mapk(ikpt)==0)then
            sumweights=sumweights+weight(ikpt)
            write(6,*)"ikpt=",ikpt
            write(6,*)"irreducible"
            write(6,fmt='(a10,3f9.6)')"internal: ",kpoint(:)/scale

            if (film.and..not.l_onedimens)then
               write(119,'(3f10.5)')kpoint(1:2),weight(ikpt)
            else
               write(119,'(4f10.5)')kpoint(:),weight(ikpt)
            endif   
      
         elseif(mapkoper(ikpt).gt.0)then
c            write(6,*)"ikpt=",ikpt
c            write(6,*)"reducible"
c            write(6,*)"map=",mapk(ikpt)
c            brot(:)=0.0
c            do k=1,3
c               brot(:)=brot(:)+
c     + mrot(k,:,mapkoper(ikpt))*kpoints(k,mapk(ikpt))
c            enddo   
c            write(6,'(a19,3f9.6)')"rotated internal: ",brot(:)/scale
         elseif(mapkoper(ikpt).lt.0)then
c            write(6,*)"ikpt=",ikpt
c            write(6,*)"reducible"
c            write(6,*)"map=",mapk(ikpt)
c            brot(:)=0.0
c            do k=1,3
c               brot(:)=brot(:)-
c     + mrot(k,:,-mapkoper(ikpt))*kpoints(k,mapk(ikpt))
c            enddo   
c            write(6,'(a19,3f9.6)')"rotated internal: ",brot(:)/scale
         endif   
        enddo !k1 
       enddo !k2 
      enddo !k3  
      close(119)

      write(6,*)"reduznumk=",reduznumk     
      write(6,*)"nkpts=",nkpts
      write(6,*)"sumweights=",sumweights
      
      IF(sumweights/=nkpts) CALL juDFT_error
     +     ("sum of weights differs from nkpts",calledby
     +     ="wann_kptsreduc2")
      deallocate(mapk,mapkoper,weight,irreduc)

      contains
      real elemental function axmintx(x)
      implicit none
      real,intent(in) :: x
      axmintx = abs(x-nint(x))
      end function     

      end subroutine wann_kptsreduc2
      end module m_wann_kptsreduc2
