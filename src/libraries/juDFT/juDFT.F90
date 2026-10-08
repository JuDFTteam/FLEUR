!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_juDFT
  USE m_juDFT_logging
  USE m_juDFT_stop
  USE m_juDFT_string
  USE m_juDFT_time
  USE m_juDFT_init
  USE m_judft_args
  USE m_judft_info
  USE m_judft_usage
  use m_judft_para
  use m_juDFT_round
  use m_npy
  implicit none
  private
  public :: log_start, log_stop, t_log_message, logmode_status, logmode_info, logmode_warning, logmode_error, &
     logmode_bug, judft_init_errormessages, judft_free_errormessages, judft_file_readable, judft_bug, judft_error, &
     judft_warn, judft_end, strip, str2int, int2str_int4, int2str_int8, float2str, gen_filename, replace_text, &
     get_byte_str, int2str, whitespaces, addtime, timestart, timestop, writetimes, writetimesxml, &
     check_time_for_next_iteration, resetiterationdependenttimers, cputime, judft_time_lastlocation, judft_init, &
     judft_was_argument, judft_string_for_argument, judft_info, judft_write_infos, send_usage_data, add_usage_data, &
     judft_check_para, check_omp_para, round_to_deci, run_sys, addrpl_cmplx_sng_vec, addrpl_cmplx_sng_mtx, &
     addrpl_cmplx_dbl_vec, addrpl_cmplx_dbl_mtx, addrpl_dbl_vec, addrpl_dbl_mtx, addrpl_sng_vec, addrpl_sng_mtx, &
     addrpl_int8_vec, addrpl_int8_mtx, addrpl_int16_vec, addrpl_int16_mtx, addrpl_int32_vec, addrpl_int32_mtx, &
     addrpl_int64_vec, addrpl_int64_mtx, write_cmplx_sgn_mtx, write_cmplx_sgn_vec, write_cmplx_dbl_6dt, &
     write_cmplx_dbl_5dt, write_cmplx_dbl_4dt, write_cmplx_dbl_3dt, write_cmplx_dbl_mtx, write_cmplx_dbl_vec, &
     write_sng_3dt, write_sng_mtx, write_sng_vec, write_dbl_3dt, write_dbl_4dt, write_dbl_5dt, write_dbl_mtx, &
     write_dbl_vec, write_int64_mtx, write_int64_vec, write_int32_mtx, write_int32_3d, write_int32_vec, &
     write_int16_mtx, write_int16_vec, write_int8_mtx, write_int8_3d, write_int8_vec, dict_str, shape_str, save_npy, &
     add_npz, p_un, magic_num, major, minor, zip_flag, magic_str
!  use m_npy
END MODULE m_juDFT
