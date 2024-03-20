!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_io_hybrid
   use m_io_matrix
   use m_judft
   use m_types
   use m_unify_zmat
   implicit none
   !private
   integer, save :: id_olap, id_z, id_v_x
   !public:: open_hybinp_io,read_cmt,write_cmt
contains
   subroutine read_z(atoms, cell, hybdat, kpts, sym, noco,nococonv, input, ik,&
                     jsp, z_out, parent_z, c_phase, list)
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
      integer, intent(in), optional                :: list(:)

      INTEGER              :: ikp, iop, i
      type(t_mat), pointer :: ptr_mat
      type(t_mat), target  :: tmp_mat
      type(t_lapw)         :: lapw_ik, lapw_ikp
      integer, allocatable :: p_list(:)
      logical, parameter   :: unify_z = .False.

      REAL :: eig(input%neig)

      call timestart("read_z")

      if(present(list)) then 
         p_list = list 
      else
         p_list = [(i, i=1,hybdat%nbands(ik,jsp))]
      endif

      if(ik <= kpts%nkpt) then
         call read_eig(hybdat%eig_id,ik,jsp,zmat=z_out, list=p_list, eig=eig)
         if(unify_z) then 
            call check_p_list(p_list, eig)
            call unify_zmat(eig, z_out)
         endif 

         if(size(p_list) /= z_out%matsize2) then
            write (*,*)  size(p_list), z_out%matsize1, z_out%matsize2
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

         call read_eig(hybdat%eig_id,ikp, jsp,zmat=ptr_mat, list=p_list, eig=eig)
         if(unify_z) then 
            call check_p_list(p_list, eig)
            call unify_zmat(eig, ptr_mat)
         endif

         if(size(p_list) /= ptr_mat%matsize2) then
            write (*,*) "list:", size(p_list)
            write (*,*) "ptr_mat", ptr_mat%matsize1, ptr_mat%matsize2
            write (*,*) "z_out", z_out%matsize1, z_out%matsize2
            call judft_error("this doesn't match ptr mat")
         endif 
         
         CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell)
         CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell)
         call waveftrafo_gen_zmat(ptr_mat, ikp, iop, kpts, sym, jsp, &
                                  size(p_list), lapw_ikp, lapw_ik, z_out, &
                                  c_phase)
      endif
      call timestop("read_z")
   END subroutine read_z

   subroutine check_p_list(p_list, eig)
      implicit none 
      integer, intent(in) :: p_list(:)
      real, intent(in)    :: eig(:) 

      integer, allocatable :: groups(:)
      logical :: succ
      integer :: n_g, end_group

      groups = make_groups(eig)
      succ = .false.

      do n_g = 1, size(groups)
         end_group = sum(groups(1:n_g))
         if(end_group == maxval(p_list)) then 
            succ = .True.
         endif 
      enddo
      if(.not. succ) call judft_error("You can't cut in the middle of deg eigenvals")
   end subroutine check_p_list
end module m_io_hybrid
