module m_ex_to_vx
   USE m_judft
   USE m_types
   USE m_symmetrizeh

contains
   subroutine ex_to_vx(fi, nk, jsp, nsymop, psym, hybdat, lapw, z, ex, v_x)
      use m_juDFT
      use m_eig66_io
      implicit none

      type(t_fleurinput), intent(in)    :: fi

      TYPE(t_hybdat), INTENT(inout) :: hybdat
      type(t_lapw), intent(in)      :: lapw
      integer, intent(in)           :: nk, jsp, nsymop, psym(fi%sym%nsym)
      TYPE(t_mat), intent(in)       :: z
      type(t_mat), intent(inout)    :: ex
      type(t_mat), intent(inout)    :: v_x

      integer     :: nbasfcn
      type(t_mat) :: trafo, tmp, olap

      CALL timestart("T^-1*mat_ex*T^-1*")
      nbasfcn = lapw%hyb_num_bas_fun(fi)

      !calculate trafo from wavefunctions to APW basis
      IF (fi%input%neig < hybdat%nbands(nk)) call judft_error(' mhsfock: neigd  < nbands(nk) ;trafo from wavefunctions to APW requires at least nbands(nk)')
      
      call ex%u2l()

      call ex%save_npy("ex_in_x2v.npy")

      call olap%init(z%l_real, z%matsize1, z%matsize1)
      CALL read_eig(hybdat%eig_id,nk,jsp, smat=olap)
      call olap%save_npy("olap.npy")
      call olap%u2l()
      call olap%conjugate()

      !call judft_error("stop after olap")

      call olap%multiply(z, res=trafo)
      call trafo%save_npy("trafo.npy")

      CALL ex%multiply(trafo, res=tmp, transB="C")
      call tmp%save_npy("tmp.npy")

      CALL trafo%multiply(tmp, res=v_x)
      call v_x%save_npy("v_x_pre_sym.npy")

      CALL timestop("T^-1*mat_ex*T^-1*")

      call timestart("symmetrizeh")
      CALL symmetrizeh(fi%atoms, fi%kpts%bkf(:, nk), jsp, lapw, fi%sym, hybdat%kveclo_eig, fi%cell, nsymop, psym, v_x)
      call v_x%save_npy("v_x_post_sym.npy")
      call timestop("symmetrizeh")
   end subroutine ex_to_vx
end module
