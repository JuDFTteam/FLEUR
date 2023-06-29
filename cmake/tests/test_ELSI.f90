      program test
      use elsi
      type(elsi_handle):: eh
      INTEGER,PARAMETER::BLACS_DENSE=0
      call elsi_init(eh,1,1,BLACS_DENSE,100,50.,50)
      call elsi_finalize(eh)
      end program
