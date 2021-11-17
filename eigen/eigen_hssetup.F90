!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_hssetup
CONTAINS
   !> The setup of the Hamiltonian and Overlap matrices are performed here
  !!
  !! The following steps are executed:
  !! 1. The matrices are a allocated (in the fi%noco-case these are 2x2-arrays of matrices)
  !! 2. The Interstitial contribution is calculated (in hs_int())
  !! 3. The MT-part is calculated (in hsmt() )
  !! 4. The vacuum part is added (in hsvac())
  !! 5. The matrices are copied to the final matrix, in the fi%noco-case the full matrix is constructed from the 4-parts.

   SUBROUTINE eigen_hssetup(isp, fmpi, fi, mpdata, results, vx, xcpot, enpara, nococonv, stars, sphhar, hybdat, &
                            ud, td, v, lapw, nk, smat_final, hmat_final)
      USE m_types
      USE m_types_mpimat
      USE m_types_gpumat
      USE m_hs_int
      USE m_hsvac
      USE m_od_hsvac
      USE m_hsmt
      USE m_eigen_redist_matrix
      USE m_add_vnonlocal
      USE m_subvxc
      USE m_eig66_io, ONLY: open_eig, write_eig, read_eig
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: isp
      TYPE(t_mpi), INTENT(IN)       :: fmpi
      type(t_fleurinput), intent(in)    :: fi
      type(t_mpdata), intent(inout):: mpdata
      type(t_results), intent(inout):: results
      class(t_xcpot), intent(in)   :: xcpot
      TYPE(t_stars), INTENT(IN)     :: stars
      TYPE(t_enpara), INTENT(IN)    :: enpara
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      type(t_hybdat), intent(inout):: hybdat
      TYPE(t_usdus), INTENT(INout)  :: ud
      TYPE(t_tlmplm), INTENT(IN)    :: td
      TYPE(t_lapw), INTENT(IN)      :: lapw
      TYPE(t_potden), INTENT(IN)    :: v, vx
      integer, intent(in)          :: nk
      CLASS(t_mat), ALLOCATABLE, INTENT(INOUT)   :: smat_final, hmat_final

      CLASS(t_mat), ALLOCATABLE :: smat(:, :), hmat(:, :)
      INTEGER :: i, j, nspins
      complex, allocatable :: vpw_wTemp(:,:)

      !Matrices for Hamiltonian and Overlapp
      !In fi%noco case we need 4-matrices for each spin channel
      nspins = MERGE(2, 1, fi%noco%l_noco)
      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::smat(nspins, nspins), hmat(nspins, nspins))
      ELSE
         ALLOCATE (t_mpimat::smat(nspins, nspins), hmat(nspins, nspins))
      END IF
      DO i = 1, nspins
         DO j = 1, nspins
            CALL smat(i, j)%init(fi%input%l_real, lapw%nv(i) + fi%atoms%nlotot, lapw%nv(j) + fi%atoms%nlotot, fmpi%sub_comm, .false.)
            CALL hmat(i, j)%init(smat(i, j))
         END DO
      END DO

      CALL timestart("Interstitial part")
      !Generate interstitial part of Hamiltonian
      ALLOCATE(vpw_wTemp(SIZE(v%pw_w,1),SIZE(v%pw_w,2)))
      vpw_wTemp = merge(v%pw_w - xcpot%get_exchange_weight() * vx%pw_w, v%pw_w, hybdat%l_subvxc)
      CALL hs_int(fi%input, fi%noco, stars, lapw, fmpi, fi%cell, isp, vpw_wTemp, smat, hmat)
      DEALLOCATE(vpw_wTemp)

      CALL timestop("Interstitial part")
      CALL timestart("MT part")
      !MT-part of Hamiltonian. In case of fi%noco, we need an loop over the local spin of the fi%atoms
      DO i = 1, nspins; DO j = 1, nspins
            !$acc enter data copyin(hmat(i,j),smat(i,j))
            !$acc enter data copyin(hmat(i,j)%data_r,smat(i,j)%data_r,hmat(i,j)%data_c,smat(i,j)%data_c)
         END DO; END DO
      CALL hsmt(fi%atoms, fi%sym, enpara, isp, fi%input, fmpi, fi%noco, nococonv, fi%cell, lapw, ud, td, smat, hmat)
      DO i = 1, nspins; DO j = 1, nspins; if (hmat(1, 1)%l_real) THEN
            !$acc exit data copyout(hmat(i,j)%data_r,smat(i,j)%data_r) delete(hmat(i,j)%data_c,smat(i,j)%data_c)
            !$acc exist data delete(hmat(i,j),smat(i,j))
         ELSE
            !$acc exit data copyout(hmat(i,j)%data_c,smat(i,j)%data_c) delete(hmat(i,j)%data_r,smat(i,j)%data_r)
            !$acc exist data delete(hmat(i,j),smat(i,j))
         END IF; END DO; END DO
      CALL timestop("MT part")

      !Vacuum contributions
      IF (fi%input%film) THEN
         CALL timestart("Vacuum part")
         CALL hsvac(fi%vacuum, stars, fmpi, isp, fi%input, v, enpara%evac, fi%cell, &
                    lapw, fi%sym, fi%noco, nococonv, hmat, smat)
         CALL timestop("Vacuum part")
      END IF

      !Deal with hybrid code
      IF (fi%hybinp%l_hybrid .OR. fi%input%l_rdmft) THEN
         if (any(shape(smat) /= 1)) then
            call judft_error("Hybrid doesn't do noco.")
         end if

         CALL write_eig(hybdat%eig_id, nk, isp, smat=smat(1, 1), n_start=fmpi%n_size, n_end=fmpi%n_rank)
      END IF

      IF (fi%hybinp%l_hybrid) THEN
         IF (hybdat%l_addhf) CALL add_Vnonlocal(nk, lapw, fi, hybdat, isp, results, xcpot, fmpi, nococonv, hmat(1, 1))

         !IF (hybdat%l_subvxc) THEN
         !   CALL subvxc(lapw, fi%kpts%bk(:, nk), fi%input, isp, v%mt(:, 0, :, :), fi%atoms, ud, &
         !               mpdata, hybdat, enpara%el0, enpara%ello0, fi%sym, &
         !               fi%cell, sphhar, stars, xcpot, fmpi, fi%oneD, hmat(1, 1), vx)
         !END IF
      END IF ! fi%hybinp%l_hybrid

      !Now copy the data into final matrix
      ! Collect the four fi%noco parts into a single matrix
      ! In collinear case only a copy is done
      ! In the parallel case also a redistribution happens
      ALLOCATE (smat_final, mold=smat(1, 1))
      ALLOCATE (hmat_final, mold=smat(1, 1))
      CALL timestart("Matrix redistribution")
      CALL eigen_redist_matrix(fmpi, lapw, fi%atoms, smat, smat_final)
      CALL eigen_redist_matrix(fmpi, lapw, fi%atoms, hmat, hmat_final, smat_final)
      CALL timestop("Matrix redistribution")

   END SUBROUTINE eigen_hssetup
END MODULE m_eigen_hssetup
