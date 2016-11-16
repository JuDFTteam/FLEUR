      program test
      real a(2,2),work(2),w(2)
      integer info

      call ssyev('N','U',2,a,2,w,work,2,info)
      end
