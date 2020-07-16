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
   subroutine read_olap(olap, rec, nbasfcn)
      implicit none
      TYPE(t_mat), INTENT(INOUT):: olap
      INTEGER, INTENT(IN)           :: rec, nbasfcn
      integer :: i, j, id

      id = open_matrix(olap%l_real, olap%matsize1, 1, 1, "olap." // int2str(rec))
      CALL read_matrix(olap, 1, id)
      call close_matrix(id)

      IF(olap%l_real) THEN
         DO i = 1, nbasfcn
            DO j = 1, i
               olap%data_r(i, j) = olap%data_r(j, i)
            END DO
         END DO
      ELSE
         DO i = 1, nbasfcn
            DO j = 1, i
               olap%data_c(i, j) = CONJG(olap%data_c(j, i))
            END DO
         END DO
         olap%data_c = conjg(olap%data_c)
      END IF
   END subroutine read_olap

   subroutine write_olap(mat, rec)
      implicit none
      TYPE(t_mat), INTENT(IN)   :: mat
      INTEGER, INTENT(IN)           :: rec
      integer :: id

      id = open_matrix(mat%l_real, mat%matsize1, 1, 1, "olap." // int2str(rec))
      CALL write_matrix(mat, 1, id)
      call close_matrix(id)
   END subroutine write_olap

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
      complex              :: cmt(input%neig,hybdat%maxlmindx,atoms%nat)
      complex              :: cmthlp(input%neig,hybdat%maxlmindx,atoms%nat)
      type(t_lapw)         :: lapw_ik, lapw_ikp

      cmt=0;cmthlp=0

      call timestart("read_z")

      if(ik <= kpts%nkpt) then
         call read_eig(hybdat%eig_id,ik,jsp,zmat=z_out, list=list)
         z_out%matsize2 = hybdat%nbands(ik)
         if(present(parent_z)) then
            call parent_z%copy(z_out,1,1)
            parent_z%matsize2=z_out%matsize2
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
         ptr_mat%matsize2 = hybdat%nbands(ik)

         CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, sym%zrfs)
         CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell, sym%zrfs)
         call waveftrafo_gen_zmat(ptr_mat, ikp, iop, kpts, sym, jsp, input, &
                                  hybdat%nbands(ikp), lapw_ikp, lapw_ik, z_out, &
                                  c_phase)
      endif
      call timestop("read_z")
   END subroutine read_z
end module m_io_hybinp
