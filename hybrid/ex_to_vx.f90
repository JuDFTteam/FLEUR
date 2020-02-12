module m_ex_to_vx

contains
   subroutine ex_to_vx()
      use m_juDFT
      implicit none
      CALL timestart("time for performing T^-1*mat_ex*T^-1*")
      !calculate trafo from wavefunctions to APW basis
      IF (input%neig < hybdat%nbands(nk)) call judft_error(' mhsfock: neigd  < nbands(nk) ;trafo from wavefunctions to APW requires at least nbands(nk)')

      call olap%multiply(z, trafo)

      CALL invtrafo%alloc(olap%l_real, hybdat%nbands(nk), nbasfcn)
      CALL trafo%TRANSPOSE(invtrafo)

      DO i = 1, hybdat%nbands(nk)
         DO j = 1, i - 1
            IF (ex%l_real) THEN
               ex%data_r(i, j) = ex%data_r(j, i)
            ELSE
               ex%data_c(i, j) = conjg(ex%data_c(j, i))
            END IF
         ENDDO
      ENDDO

      CALL ex%multiply(invtrafo, tmp)
      CALL trafo%multiply(tmp, v_x)

      CALL timestop("time for performing T^-1*mat_ex*T^-1*")

      call timestart("symmetrizeh")
      CALL symmetrizeh(atoms, kpts%bkf(:, nk), jsp, lapw, sym, hybdat%kveclo_eig, cell, nsymop, psym, v_x)
      call timestop("symmetrizeh")

      CALL write_v_x(v_x, kpts%nkpt*(jsp - 1) + nk)
   end subroutine ex_to_vx
end module
