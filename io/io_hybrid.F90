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
   integer, save :: id_olap, id_z, id_v_x, id_coulomb, id_coulomb_spm
   !public:: open_hybinp_io,read_cmt,write_cmt
contains
   SUBROUTINE open_hybinp_io2(mpdata, hybinp, hybdat, input, atoms, l_real)
      IMPLICIT NONE
      type(t_mpdata), intent(in) :: mpdata
      TYPE(t_hybinp), INTENT(IN)  :: hybinp
      TYPE(t_hybdat), INTENT(IN)  :: hybdat
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_atoms), INTENT(IN)   :: atoms
      LOGICAL, INTENT(IN)         :: l_real
      INTEGER:: irecl_coulomb

      ! if the sparse matrix technique is used, several entries of the
      ! matrix vanish so that the size of each entry is smaller
      irecl_coulomb = (atoms%ntype*(maxval(hybinp%lcutm1) + 1)*(maxval(mpdata%num_radbasfn) - 1)**2 &
                       + atoms%nat*(maxval(hybinp%lcutm1) + 2)*(2*maxval(hybinp%lcutm1) + 1)*(maxval(mpdata%num_radbasfn) - 1) &
                       + (maxval(mpdata%num_radbasfn) - 1)*atoms%nat**2 &
                       + ((maxval(hybinp%lcutm1) + 1)**2*atoms%nat + maxval(mpdata%n_g)) &
                       *((maxval(hybinp%lcutm1) + 1)**2*atoms%nat + maxval(mpdata%n_g) + 1)/2)*8
      if(.not. l_real) irecl_coulomb = irecl_coulomb*2
      OPEN(unit=778, file='coulomb1', form='unformatted', access='direct', recl=irecl_coulomb)
      id_coulomb_spm = 778
   END SUBROUTINE open_hybinp_io2

   subroutine close_hybinp_io2
      implicit none

      close(778)
   end subroutine close_hybinp_io2

   subroutine write_coulomb(nk, l_real, coulomb)
      implicit none
      complex, intent(in) :: coulomb(:)
      integer, intent(in) :: nk
      logical, intent(in) :: l_real

      if(l_real) THEN
         write(id_coulomb, rec=nk) real(coulomb)
      else
         write(id_coulomb, rec=nk) coulomb
      end if
   end subroutine write_coulomb

   subroutine write_coulomb_spm_r(nk, coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir)
      implicit none
      real, intent(in)    :: coulomb_mt1(:, :, :, :)
      real, intent(in) :: coulomb_mt2(:, :, :, :), coulomb_mt3(:, :, :)
      real, intent(in) :: coulomb_mtir(:)
      integer, intent(in) :: nk
      write(id_coulomb_spm, rec=nk) coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir
   end subroutine write_coulomb_spm_r

   subroutine write_coulomb_spm_c(nk, coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir)
      implicit none
      real, intent(in)    :: coulomb_mt1(:, :, :, :)
      complex, intent(in) :: coulomb_mt2(:, :, :, :), coulomb_mt3(:, :, :)
      complex, intent(in) :: coulomb_mtir(:)
      integer, intent(in) :: nk

      write(id_coulomb_spm, rec=nk) coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir
   end subroutine write_coulomb_spm_c

   subroutine read_coulomb_spm_r(nk, coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir)
      implicit none
      real, intent(out)    :: coulomb_mt1(:, :, :, :)
      real, intent(out) :: coulomb_mt2(:, :, :, :), coulomb_mt3(:, :, :)
      real, intent(out) :: coulomb_mtir(:)
      integer, intent(in) :: nk

      !print *, "read coulomb",nk,size(coulomb_mt1),size(coulomb_mt2),size(coulomb_mt3),size(coulomb_mtir)
      read(id_coulomb_spm, rec=nk) coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir
   end subroutine read_coulomb_spm_r

   subroutine read_coulomb_spm_c(nk, coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir)
      implicit none
      real, intent(out)    :: coulomb_mt1(:, :, :, :)
      complex, intent(out) :: coulomb_mt2(:, :, :, :), coulomb_mt3(:, :, :)
      complex, intent(out) :: coulomb_mtir(:)
      integer, intent(in) :: nk
      read(id_coulomb_spm, rec=nk) coulomb_mt1, coulomb_mt2, coulomb_mt3, coulomb_mtir
   end subroutine read_coulomb_spm_c

   subroutine read_coulomb_r(nk, coulomb)
      implicit none
      real, intent(out) :: coulomb(:)
      integer, intent(in) :: nk

      read(id_coulomb, rec=nk) coulomb
   end subroutine read_coulomb_r

   subroutine read_coulomb_c(nk, coulomb)
      implicit none
      complex, intent(out) :: coulomb(:)
      integer, intent(in) :: nk

      read(id_coulomb, rec=nk) coulomb
   end subroutine read_coulomb_c

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
                     jsp, z_out, parent_z, c_phase)
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

      INTEGER              :: ikp, iop
      type(t_mat), pointer :: ptr_mat
      type(t_mat), target  :: tmp_mat
      complex              :: cmt(input%neig,hybdat%maxlmindx,atoms%nat)
      complex              :: cmthlp(input%neig,hybdat%maxlmindx,atoms%nat)
      type(t_lapw)         :: lapw_ik, lapw_ikp

      cmt=0;cmthlp=0

      call timestart("read_z")
      if(ik <= kpts%nkpt) then
         call read_eig(hybdat%eig_id,ik,jsp,zmat=z_out)
         ! z_out%matsize2 = hybdat%nbands(ik)
         ! call z_out%save_npy("z_ik=" // int2str(ik) // ".npy")
         if(present(parent_z)) call parent_z%copy(z_out,1,1)
      else
         if(present(parent_z)) then
            ptr_mat => parent_z
         else
            call tmp_mat%init(z_out)
            ptr_mat => tmp_mat
         endif

         ikp = kpts%bkp(ik) ! parrent k-point
         iop = kpts%bksym(ik) ! connecting symm

         call read_eig(hybdat%eig_id,ikp, jsp,zmat=ptr_mat)
         ! ptr_mat%matsize2 = hybdat%nbands(ik)
         ! call ptr_mat%save_npy("z_ik=" // int2str(ik) // "_ikp=" // int2str(ikp) // ".npy")

         CALL lapw_ik%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, sym%zrfs)
         CALL lapw_ikp%init(input, noco, nococonv, kpts, atoms, sym, ikp, cell, sym%zrfs)
         call waveftrafo_gen_zmat(ptr_mat, ikp, iop, kpts, sym, jsp, input, &
                                  hybdat%nbands(ikp), lapw_ikp, lapw_ik, z_out, &
                                  c_phase)
      endif
      call timestop("read_z")
   END subroutine read_z

   subroutine read_v_x(mat, rec)
      implicit none
      TYPE(t_mat), INTENT(INOUT)  :: mat
      INTEGER, INTENT(IN)         :: rec
      integer                     :: id

      id = open_matrix(mat%l_real, mat%matsize1, 1, 1, "v_x." // int2str(rec))
      CALL read_matrix(mat, 1, id)
      call close_matrix(id)
   END subroutine read_v_x

   subroutine write_v_x(mat, rec)
      use m_juDFT
      implicit none
      TYPE(t_mat), INTENT(IN)   :: mat
      INTEGER, INTENT(IN)       :: rec
      integer :: id

      id = open_matrix(mat%l_real, mat%matsize1, 1, 1, "v_x." // int2str(rec))
      CALL write_matrix(mat, 1, id)
      call close_matrix(id)
   END subroutine write_v_x

end module m_io_hybinp
