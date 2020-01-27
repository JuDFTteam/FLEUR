      program test
      use hdf5
      integer:: error
      CALL h5open_f(error)
      CALL h5close_f(error)

      end
