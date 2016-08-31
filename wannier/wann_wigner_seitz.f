!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_wigner_seitz
      use m_juDFT
      contains
      subroutine wann_wigner_seitz(
     >          l_get_rvecnum,num,amat,
     >          rvecnum_in,
     <          rvecnum,rvec,ndegen)

      implicit none
      logical, intent(in) :: l_get_rvecnum
      integer, intent(in) :: num(:)
      real, intent(in)    :: amat(:,:)
      integer,intent(in)  :: rvecnum_in
      integer, intent(out):: rvecnum
      integer, intent(out):: rvec(:,:)
      integer, intent(out):: ndegen(:)	

      integer             :: idist (3)
      real                :: dist(125),summa,dist_min
      integer             :: k1,k2,k3,i1,i2,i3,count,i,j
      real                :: eps7,eps8
      real                :: metric(3,3)

      eps7=1.e-7
      eps8=1.e-8
      rvecnum=0

      metric=matmul(transpose(amat),amat)

      rvecnum = 0  
      do k1=-num(1),num(1)  
       do k2=-num(2),num(2)  
        do k3=-num(3),num(3)  
         count=0  
         do i1=-2,2  
          do i2=-2,2  
           do i3=-2,2  
            count=count+1  
            ! Get |r-R|^2
            idist(1)=k1-i1*num(1)  
            idist(2)=k2-i2*num(2)  
            idist(3)=k3-i3*num(3)  
            dist(count)=0.0 
            do i=1,3  
             do j=1,3  
              dist(count)=dist(count)+ 
     +         real(idist(i))*metric(i,j)*real(idist(j))
             enddo !i
            enddo !j
           enddo !i3
          enddo !i2
         enddo !i1
         dist_min=minval(dist)
         if (abs(dist(63) - dist_min ) .lt. eps7 ) then
                rvecnum = rvecnum + 1  
                if(.not. l_get_rvecnum) then
c                   if(.not.allocated(ndegen))
c     &                 allocate(ndegen(rvecnum_in))
                   ndegen(rvecnum)=0
                   do i=1,125
                      if (abs (dist (i) - dist_min) .lt. eps7 ) 
     &                   ndegen(rvecnum)=ndegen(rvecnum)+1
                   end do
                   rvec(1,rvecnum) = k1  
                   rvec(2,rvecnum) = k2   
                   rvec(3,rvecnum) = k3   
                endif
         endif
        enddo !k3
       enddo !k2
      enddo !k1
      if(l_get_rvecnum) return

c------ Consistency Check.
      summa = 0.0
      do i = 1, rvecnum  
       summa = summa + 1.0/real(ndegen(i))  
      enddo
      if (abs (summa - real(num(1)*num(2)*num(3)) ) > eps8) then
         CALL juDFT_error("problem finding Wigner-Seitz points",calledby
     +        ="wann_wigner_seitz")
      endif

      end subroutine wann_wigner_seitz

      end module m_wann_wigner_seitz
