      program test  !Basically test if the elpa module is available and if it can diagonalize a single precision matrix
      use elpa
      type(elpa_t),pointer :: elpa
      integer:: err
      integer, parameter:: sp = selected_real_kind(6)
      real(kind=sp):: h(2,2),e(2),z(2,2)

      call elpa%eigenvalues(h,e,z,err)
      
      end program
