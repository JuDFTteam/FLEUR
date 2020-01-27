c*****************************c
c  routine to define Bravais  c
c  in the auxiliary space     c
c*****************************c
      module m_wann_gwf_auxbrav
      implicit none
      contains

      subroutine wann_gwf_auxbrav(a,l_sgwf,l_socgwf,amat,bmat,l_dim)
      use m_inv3
      use m_constants, only : pimach

      implicit none
      real,    intent(in)  :: a
      logical, intent(in)  :: l_sgwf,l_socgwf,l_dim(3)
      real,    intent(out) :: amat(3,3),bmat(3,3)
      real :: tpi, omtil,kmesh_tol

      write(*,*)'auxiliary lattice constant a=',a

      tpi = 2.0*pimach()      
      kmesh_tol = 0.001 
      amat = 0.0
      bmat = 0.0

      amat(1,1)=1.1*a; amat(2,2)=1.1*a; amat(3,3)=1.1*a
      if(l_dim(1))amat(1,1)=a
      if(l_dim(2))amat(2,2)=a
      if(l_dim(3))amat(3,3)=a
      if(l_sgwf.or.l_socgwf) then
         amat = amat/2.0
c         amat(1,1)=amat(1,1)-kmesh_tol
c         amat(2,2)=amat(2,2)-kmesh_tol
c         amat(3,3)=amat(3,3)-kmesh_tol
      endif
      call inv3(amat,bmat,omtil)
      bmat = tpi*bmat

      !write(*,*)'aux. amat',amat
      !write(*,*)'aux. bmat',bmat
      end subroutine wann_gwf_auxbrav

      end module m_wann_gwf_auxbrav
