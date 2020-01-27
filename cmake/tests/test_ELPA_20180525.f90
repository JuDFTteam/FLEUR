      program test
      use elpa

      implicit none
      class(elpa_t),pointer :: elpa1
      integer :: success,na
      REAL :: H(2,2),s(2,2),eig(2),ev(2,2)
      
      call elpa1%set("na",na,success) 
      CALL elpa1%generalized_eigenvectors(h,s,eig, ev, .FALSE.)
      end
