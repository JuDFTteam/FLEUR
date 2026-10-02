!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hdf_tools
!-----------------------------------------------
!     major rewrite of hdf_tools
!     this module only collects the features
!
!-----------------------------------------------
   USE m_hdf_tools1
   USE m_hdf_tools2
   USE m_hdf_tools3
   USE m_hdf_tools4
   USE m_hdf_tools5
   USE m_hdf_tools6
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: io_write_attchar0, io_read_attchar0, io_write_attlog0, io_read_attlog0, io_write_attreal0, &
      io_write_attreal1, io_write_attreal2, io_write_attreal3, io_read_attreal0, io_read_attreal1, io_read_attreal2, &
      io_read_attreal3, io_write_attint0, io_write_attint1, io_write_attint2, io_write_attint3, io_read_attint0, &
      io_read_attint1, io_read_attint2, io_read_attint3, io_write_att, io_read_att, io_read_real0, io_read_real1, &
      io_read_real2, io_read_real3, io_read_real4, io_read_real5, io_read_real6, io_write_real0, io_write_real1, &
      io_write_real2, io_write_real3, io_write_real4, io_write_real5, io_write_real6, io_read_integer0, &
      io_read_integer1, io_read_integer2, io_read_integer3, io_read_integer4, io_read_integer5, io_read_integer6, &
      io_write_integer0, io_write_integer1, io_write_integer2, io_write_integer3, io_write_integer4, io_write_integer5, &
      io_write_integer6, io_read_complex0, io_read_complex1, io_read_complex2, io_read_complex3, io_read_complex4, &
      io_read_complex5, io_write_complex0, io_write_complex1, io_write_complex2, io_write_complex3, io_write_complex4, &
      io_write_complex5, io_read, io_write, io_dataexists, io_attexists, io_groupexists, io_hdfopen, io_hdfclose, &
      io_gopen, io_gcreate, io_gdelete, io_gclose, io_dopen, io_dclose, io_layername, hdf_init, hdf_close, &
      io_createvar, gettransprop, cleartransprop, io_check, checklib, hdf_err, io_datadim, rkind, ckind, &
      io_write_real3s, io_write_real1s, io_write_real2s, io_write_var, io_read_var
END MODULE m_hdf_tools
