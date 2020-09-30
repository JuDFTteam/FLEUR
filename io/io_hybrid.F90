!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_io_hybinp
   use m_io_matrix
   use m_judft
   use m_types
   implicit none
   !private
   integer, save :: id_olap, id_z, id_v_x
   !public:: open_hybinp_io,read_cmt,write_cmt
contains
   subroutine read_z(atoms, cell, hybdat, kpts, sym, noco,nococonv, input, ik,&
                     jsp, z_out, parent_z, c_phase, p_list)
      USE m_eig66_io
      use m_types
      use m_trafo
      implicit none
      type(t_atoms), intent(in)    :: atoms
      type(t_cell), intent(in)     :: cell
      type(t_hybdat), intent(in)   :: hybdat
      type(t_kpts), intent(in)     :: kpts
      type(t_sym), intent(in)      :: sym
      type(t_noco), intent(in)     :: noco
      TYPE(t_nococonv),INTENT(IN)  :: nococonv
      type(t_input), intent(in)    :: input
      integer, intent(in)          :: ik, jsp
      TYPE(t_mat), INTENT(INOUT)   :: z_out

      type(t_mat), intent(inout), target, optional :: parent_z
      complex, intent(inout), optional             :: c_phase(:)
      integer, intent(in), optional                :: p_list(:)

      INTEGER              :: ikp, iop, i
      type(t_mat), pointer :: ptr_mat
      type(t_mat), target  :: tmp_mat
      type(t_lapw)         :: lapw_ik, lapw_ikp
      integer, allocatable :: list(:)

      call timestart("read_z")

      if(present(p_list)) then 
         list = p_list 
      else
         list = [(i, i=1,hybdat%nbands(ik))]
      endif

      if(ik <= kpts%nkpt) then
         call read_eig(hybdat%eig_id,ik,jsp,zmat=z_out, list=list)
         if(size(list) /= z_out%matsize2) then
            write (*,*)  size(list), z_out%matsize1, z_out%matsize2
            call judft_error("this doesn't match")
         endif
         if(present(parent_z)) then
            call parent_z%copy(z_out,1,1)
         endif
      else
         if(present(parent_z)) then
            ptr_mat => parent_z
         else
            call tmp_mat%init(z_out)
            ptr_mat => tmp_mat
         endif

         ikp = kpts%bkp(ik) ! parrent k-point
         iop = kpts%bksym(ik) ! connecting symm
         call read_eig(hybdat%eig_id,ikp, jsp,zmat=ptr_mat, list=list)
         if(size(list) /= ptr_mat%matsize2) then
            write (*,*) "list:", size(list)
            write (*,*) "ptr_mat", ptr_mat%matsize1, ptr_mat%matsize2
            write (*,*) "z_out", z_out%matsize1, z_out%matsize2
            call judft_error("this doesn't match ptr mat")
         endif 
         
         CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, sym%zrfs)
         CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell, sym%zrfs)
         call waveftrafo_gen_zmat(ptr_mat, ikp, iop, kpts, sym, jsp, input, &
                                  size(list), lapw_ikp, lapw_ik, z_out, &
                                  c_phase)
      endif
      call timestop("read_z")
   END subroutine read_z
end module m_io_hybinp
