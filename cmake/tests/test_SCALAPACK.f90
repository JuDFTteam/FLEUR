      program test

      !call a blacs routine
      CALL BLACS_PINFO(1,1)
      

      CALL pzhegvx(1,'V','I','U',1,1.0,1,1,(/1,1/),1.0)
      end
