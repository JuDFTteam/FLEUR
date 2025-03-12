!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_radfun
   use m_judft
   implicit none
   private
   type:: t_radfun
      integer, private:: itype = 0
      integer, allocatable:: n_r(:) !! number of radial functions per l-chanel
      real, allocatable:: r(:, :, :, :,:) !!radial function (index,jmtd,1:2,l,spin)
   contains
      procedure, pass :: init
      procedure, pass :: generate_radial_functions
   end type
   public:: t_radfun
contains
   subroutine init(this, atoms, input, itype)
      use m_types
      implicit none
      class(t_radfun), intent(inout):: this
      type(t_atoms), intent(IN)   :: atoms
      type(t_input), intent(in)     :: input
      integer, intent(in)           :: itype

      integer:: l, lo
      this%itype = itype
      do l = 1, atoms%lmax(itype)
         this%n_r(l) = 2
         if (input%l_useapw) call judft_bug("APW not implemented")
         do lo = 1, atoms%nlo(itype)
            if (l /= atoms%llo(lo, itype)) cycle !no LO for this l
            this%n_r(l) = this%n_r(l) + 1
         end do
      end do
   end subroutine

   subroutine generate_radial_functions(this, atoms, input, enpara, hub1data, fmpi, vtot, iType)
      use m_genMTBasis
      use m_types
      implicit none
      class(t_radfun), intent(inout)         ::this
      type(t_atoms), intent(IN)      :: atoms
      type(t_input), intent(IN)    :: input
      type(t_enpara), intent(IN)   :: enpara
      type(t_hub1data), intent(IN) :: hub1data
      type(t_mpi), intent(IN)     :: fmpi
      type(t_potden), intent(IN)   :: vtot
      integer, intent(in)                    :: itype

      !temp variables not really used but required by genMTBasis
      type(t_usdus) :: usdus
      !radial functions to copy into type
      real            :: f(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: g(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: flo(atoms%jmtd, 2, atoms%nlod)

      integer:: ispin, l, n, lo

      !check if data is already available
      if (this%itype == itype .and. allocated(this%r)) return
      !init type
      if (allocated(this%r)) deallocate (this%r)
      call this%init(atoms, input, itype)

      allocate (this%r(maxval(this%n_r), atoms%jmtd, 2, atoms%lmaxd, input%jspins))

      do ispin = 1, input%jspins
         call genMTBasis(atoms, enpara, vTot, fmpi, iType, ispin, usdus, f, g, flo, hub1data, l_writeArg=.false.)
         do l = 1, atoms%lmax(itype)
            this%R(1, 1:atoms%jri(itype), 1:2, l, ispin) = f(1:atoms%jri(itype), 1:2, l)
            this%R(2, 1:atoms%jri(itype), 1:2, l, ispin) = g(1:atoms%jri(itype), 1:2, l)
            n = 2
            if (input%l_useapw) call judft_bug("APW not implemented")
            do lo = 1, atoms%nlo(itype)
               if (l /= atoms%llo(lo, itype)) cycle !no LO for this l
               n = n + 1
               this%R(n, 1:atoms%jri(itype), 1:2, l, ispin) = flo(1:atoms%jri(itype), 1:2, lo)
            end do
         end do
      end do
   end subroutine

   

end module m_types_radfun
