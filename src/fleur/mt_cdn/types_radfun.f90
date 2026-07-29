!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
      real, allocatable:: integral(:,:,:,:,:) !! int of r over MT sphere (index,index,l,spin,spin)
   contains
      procedure, pass :: init
      procedure, pass :: generate_radial_functions
   end type
   public:: t_radfun
contains
   pure subroutine init(this, atoms, input, itype)
      use m_types_atoms
      use m_types_input
      implicit none
      class(t_radfun), intent(inout):: this
      type(t_atoms), intent(IN)   :: atoms
      type(t_input), intent(in)     :: input
      integer, intent(in)           :: itype

      integer:: l, lo
      this%itype = itype
      if (allocated(this%n_r)) deallocate(this%n_r)
      allocate(this%n_r(0:atoms%lmaxd))
      this%n_r(0:)=atoms%num_radial_functions_per_l(itype)
   end subroutine

   subroutine generate_radial_functions(this, atoms, input, enpara, fmpi, vtot, iType, hub1data, usdus_out)
      use m_genMTBasis
      use m_types_atoms
      use m_types_input
      use m_types_enpara
      use m_types_hub1data
      use m_types_mpi
      use m_types_potden
      use m_types_usdus
      use m_intgr
      implicit none
      class(t_radfun), intent(inout)         ::this
      type(t_atoms), intent(IN)      :: atoms
      type(t_input), intent(IN)    :: input
      type(t_enpara), intent(IN)   :: enpara
      type(t_hub1data), intent(IN),optional :: hub1data
      type(t_mpi), intent(IN)     :: fmpi
      type(t_potden), intent(IN)             :: vtot
      integer, intent(in)                    :: itype
      !accumulates data for all itype this routine is called for; on a cache
      !hit (data for itype already generated) it is left untouched
      type(t_usdus), intent(INOUT),optional  :: usdus_out

      !temp variables not really used but required by genMTBasis
      type(t_usdus) :: usdus_tmp
      !radial functions to copy into type
      real            :: f(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: g(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: flo(atoms%jmtd, 2, atoms%nlod)

      integer:: ispin, jspin, i,j, l, n, lo
      real,allocatable:: rf(:)
      real :: ovlp
      call timestart("generate radial functions")
      if (present(usdus_out)) then
         if (.not.allocated(usdus_out%us)) call usdus_out%init(atoms, input%jspins)
      else
         call usdus_tmp%init(atoms,input%jspins)
      end if


      !check if data is already available
      if (this%itype /= itype .or. .not.allocated(this%r)) THEN
         !init type
         call this%init(atoms, input, itype)

         if (allocated(this%r)) deallocate (this%r)
         allocate (this%r( atoms%jmtd, 2,maxval(this%n_r),0:atoms%lmaxd, input%jspins),source=0.0)
         if (allocated(this%integral)) deallocate (this%integral)
         allocate (this%integral(maxval(this%n_r), maxval(this%n_r),0:atoms%lmaxd, input%jspins,input%jspins),source=0.0)

         do ispin = 1, input%jspins
            if (present(usdus_out)) then
               call genMTBasis(atoms, enpara, vTot, fmpi, iType, ispin, usdus_out, f, g, flo, hub1data, l_writeArg=.false.)
            else
               call genMTBasis(atoms, enpara, vTot, fmpi, iType, ispin, usdus_tmp, f, g, flo, hub1data, l_writeArg=.false.)
            end if
            do l = 0, atoms%lmax(itype)
               this%R( 1:atoms%jri(itype), 1:2, 1,l, ispin) = f(1:atoms%jri(itype), 1:2, l)
               this%R( 1:atoms%jri(itype), 1:2, 2,l, ispin) = g(1:atoms%jri(itype), 1:2, l)
               n = 2
               if (input%l_useapw) call judft_bug("APW not implemented")
               do lo = 1, atoms%nlo(itype)
                  if (l /= atoms%llo(lo, itype)) cycle !no LO for this l
                  n = n + 1
                  this%R( 1:atoms%jri(itype), 1:2, n,l, ispin) = flo(1:atoms%jri(itype), 1:2, lo)
               end do
            end do
         end do


         !Calculate the overlaps
         DO ispin=1,input%jspins
            DO jspin=1,input%jspins
               DO l=0,atoms%lmax(itype)
                  DO i=1,this%n_r(l)
                     DO j=1,i
                        rf=this%r(1:atoms%jri(itype),1,i,l,ispin)*this%r(1:atoms%jri(itype),1,j,l,jspin)&
                        +this%r(1:atoms%jri(itype),2,i,l,ispin)*this%r(1:atoms%jri(itype),2,j,l,jspin)
                        CALL intgr0(rf,atoms%rmsh(1,itype),atoms%dx(itype),atoms%jri(itype),ovlp)
                        this%integral(i,j,l,ispin,jspin)=ovlp
                        this%integral(j,i,l,ispin,jspin)=ovlp
                     enddo
                  enddo
               ENDDO
            ENDDO
         ENDDO         
      end if
      call timestop("generate radial functions")

   end subroutine

   

end module m_types_radfun
