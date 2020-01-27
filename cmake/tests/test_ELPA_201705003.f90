      program test
      use elpa

      implicit none
      class(elpa_t),pointer :: elpa1
      integer :: success,na

      call elpa1%set("na",na,success) 

      end
