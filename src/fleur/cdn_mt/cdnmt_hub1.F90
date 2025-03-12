! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_cdnmt_hub1

contains

   subroutine cdnmt_hub1(atoms, jsp_start, jsp_end, vr, denCoeffs, usdus, hub1inp, hub1data)
      !!
      !! Calculate the atomic density needed for Hubbard-1 calculations

      use m_types
      use m_constants
      use m_radfun

      implicit none

      type(t_usdus), intent(INOUT) :: usdus !in fact only the lo part is intent(in)
      type(t_atoms), intent(IN)    :: atoms
      type(t_hub1inp), intent(IN)    :: hub1inp
      type(t_denCoeffs), intent(IN)    :: denCoeffs
      real, intent(IN)                          :: vr(:, :, :)
      type(t_hub1data), intent(INOUT) :: hub1data
      integer, intent(in):: jsp_start, jsp_end

      integer :: itype, l, llp, ispin, noded, nodeu
      logical :: l_hia, l_performSpinavg

      real              :: vrTmp(atoms%jmtd)

      real, allocatable :: Rf(:, :, :, :, :)

      do itype = 1, atoms%ntype
         do ispin = jsp_start, jsp_end
            do l = 0, lmaxU_const
               !Check if the orbital is treated with Hubbard 1
               l_hia = any(atoms%lda_u(atoms%n_u + 1:atoms%n_u + atoms%n_hia)%atomType == itype &
                           .and. atoms%lda_u(atoms%n_u + 1:atoms%n_u + atoms%n_hia)%l == l)

               !In the case of a spin-polarized calculation with Hubbard 1 we want to treat
               !the correlated orbitals with a non-spin-polarized basis
               if (l_hia .and. atoms%l_nonpolbas(itype) .or. hub1data%l_performSpinavg) then
                  vrTmp = (vr(:, itype, 1) + vr(:, itype, 2))/2
               else
                  vrTmp = vr(:, itype, ispin)
               end if

               call radfun(l, itype, ispin, enpara%el0(l, itype, ispin), vrTmp, atoms, &
                           Rf(1, 1, 0, l, ispin), Rf(1, 1, 1, l, ispin), usdus, nodeu, noded, wronk)

               llp = l*(atoms%lmax(itype) + 1) + l
               hub1data%cdn_atomic(j, l, itype, ispin) = hub1data%cdn_atomic(j, l, itype, ispin) &
                                                         + real(denCoeffs%nmt_coeff(llp, 0, itype, 0, 0, ispin, ispin)) &
                                                         *(Rf(j, 1, 0, l, ispin)*Rf(j, 1, 0, l, ispin) &
                                                           + Rf(j, 2, 0, l, ispin)*Rf(j, 2, 0, l, ispin)) &
                                                         *1.0/(atoms%neq(itype)*sfp_const)

            end do
         end do
      end do

   end subroutine cdnmt_hub1
end module m_cdnmt_hub1
