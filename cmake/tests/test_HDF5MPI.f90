      program test
      use hdf5
      integer:: error
      integer(hid_t)   :: access_prp
      CALL h5open_f(error)
      CALL h5pset_fapl_mpio_f(access_prp, 1,1,error)
      CALL h5close_f(error)

      end
