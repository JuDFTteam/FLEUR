      module m_wann_gwf_auxovlp
      implicit none
      contains

      subroutine wann_gwf_auxovlp(q1,q2,a,mmnaux)
      use m_constants, only : pimach

      implicit none
      real,    intent(in)  :: q1(3),q2(3),a
      complex, intent(out) :: mmnaux

      integer :: i
      real :: pi, a_aux, b_aux, tau(3)
      complex :: ov(3)

      pi = pimach()
      a_aux = a
      b_aux = 0.98*a
      tau = 2.0*pi*(q1-q2)/a_aux
      ov = cmplx(1.,0.)

      do i=1,3
         if(abs(tau(i)).lt.1e-10) cycle
         ov(i) = 8.0*pi*pi*sin(tau(i)*b_aux/2.0)
         ov(i) = ov(i) / ( 4.0*pi*pi - (tau(i)*b_aux)**2 )
         ov(i) = ov(i) / ( tau(i)*b_aux )
      enddo
      mmnaux = ov(1)*ov(2)*ov(3)

      write(*,*)'ov(1) =',ov(1)
      write(*,*)'ov(2) =',ov(2)
      write(*,*)'ov(3) =',ov(3)
      write(*,*)'mmnaux=',mmnaux

      end subroutine wann_gwf_auxovlp

      end module m_wann_gwf_auxovlp
