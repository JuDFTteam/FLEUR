      MODULE m_wann_plot_symm
      CONTAINS
      SUBROUTINE wann_plot_symm(jspin,mrot,ikpt,kplot,l_conjugate)

      implicit none
      integer,intent(in)::jspin,ikpt,kplot
      logical,intent(in)::l_conjugate
      integer,intent(in)::mrot(3,3)
      character(len=10)::vandernameold
      character(len=10)::vandernamenew
      integer::grid(3),nslibd,ikpt_copy,nbnd
      integer::ix,iy,iz,ixx,iyy,izz
      complex::xdnout
      real,allocatable::xdnreal(:,:,:,:),xdnimag(:,:,:,:)
      real::point(3)
      real::rotpoint(3)
      integer::opoint1,opoint2,opoint3

      WRITE (vandernameold,201) ikpt,jspin
      write (vandernamenew,201) kplot,jspin
 201  FORMAT ('UNK',i5.5,'.',i1)      
      open(55,file=vandernameold)
      read(55,*)grid(1),grid(2),grid(3),ikpt_copy,nslibd
      allocate(xdnreal(0:grid(1)-1,0:grid(2)-1,0:grid(3)-1,nslibd))
      allocate(xdnimag(0:grid(1)-1,0:grid(2)-1,0:grid(3)-1,nslibd))
      do nbnd=1,nslibd
       do iz=0,grid(3)-1
        do iy=0,grid(2)-1
         do ix=0,grid(1)-1
            read(55,*)xdnreal(ix,iy,iz,nbnd),xdnimag(ix,iy,iz,nbnd)
         enddo
        enddo
       enddo
      enddo
      close(55)
      open(666,file=vandernamenew,form='formatted')
      write(666,'(5i4)')grid(1),grid(2),grid(3),kplot,nslibd
      do nbnd=1,nslibd
       do iz=0,grid(3)-1
        do iy=0,grid(2)-1
         do ix=0,grid(1)-1
          point(1)=real(ix)/grid(1)
          point(2)=real(iy)/grid(2)
          point(3)=real(iz)/grid(3)
          rotpoint(:)=matmul(point(1:3),1.0*mrot(1:3,1:3))
          where(rotpoint.lt.0)
             rotpoint=rotpoint+1.0
          endwhere
          opoint1=int(grid(1)*rotpoint(1)+0.01)
          opoint2=int(grid(2)*rotpoint(2)+0.01)
          opoint3=int(grid(3)*rotpoint(3)+0.01)
          if(l_conjugate)then
            write(666,*)xdnreal(opoint1,opoint2,opoint3,nbnd),
     &               -xdnimag(opoint1,opoint2,opoint3,nbnd)
          else
            write(666,*)xdnreal(opoint1,opoint2,opoint3,nbnd),
     &               xdnimag(opoint1,opoint2,opoint3,nbnd)
          endif
         enddo
        enddo
       enddo
      enddo
      close(666)

      END SUBROUTINE wann_plot_symm
      END MODULE m_wann_plot_symm
