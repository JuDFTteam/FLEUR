      program test
      use elsi
      use elsi_precision, only: r8
      type(elsi_handle):: eh
      INTEGER,PARAMETER::BLACS_DENSE=0
      REAL(r8)::elec
      elec=50
      call elsi_init(eh,1,1,BLACS_DENSE,100,elec,50)
      call elsi_finalize(eh)
      end program
