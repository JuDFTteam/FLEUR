!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_types_radfun
   !! Radial basis functions of one atom type in the unified slot basis:
   !! per l, slot 1=u, 2=udot, 3.. = LOs of that l in atoms%llo order (atoms%slot_of_lo).
   !! Each slot solves H r_i = e_i r_i + P r_i inside the sphere (P u = 0, P udot = u,
   !! P ulo = filo for energy-derivative LOs), except relLOs which are Dirac solutions.
   use m_judft
   implicit none
   private
   type:: t_radfun
      integer, private:: itype = 0
      real, private   :: rmt = 0.0
      integer, allocatable:: n_r(:) !! number of radial functions per l-chanel
      real, allocatable:: r(:, :, :, :,:) !!radial function (index,jmtd,1:2,l,spin)
      real, allocatable:: integral(:,:,:,:,:) !! int of r over MT sphere (index,index,l,spin,spin)
      real, allocatable:: bnd(:,:,:,:)   !! value and radial derivative at R_MT (1:2,index,l,spin)
      real, allocatable:: e(:,:,:)        !! energy parameter (index,l,spin)
      real, allocatable:: ipred(:,:,:,:)  !! <r_i|P r_j> (index,index,l,spin)
      logical, allocatable:: l_sra(:,:)   !! scalar-relativistic solution, false for relLOs (index,l)
   contains
      procedure, pass :: init
      procedure, pass :: generate_radial_functions
      procedure, pass :: reset
      procedure, pass :: hsph
      procedure, pass :: to_usdus
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

      this%itype = itype
      this%rmt = atoms%rmt(itype)
      if (allocated(this%n_r)) deallocate(this%n_r)
      allocate(this%n_r(0:atoms%lmaxd))
      this%n_r(0:)=atoms%num_radial_functions_per_l(itype)
   end subroutine

   subroutine generate_radial_functions(this, atoms, input, enpara, fmpi, vtot, iType, hub1data,usdus_out)
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
      !receives the data of itype; entries of other types are left untouched
      type(t_usdus), intent(INOUT),optional  :: usdus_out

      type(t_usdus) :: usdus
      real            :: f(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: g(atoms%jmtd, 2, 0:atoms%lmaxd)
      real            :: flo(atoms%jmtd, 2, atoms%nlod)

      integer:: ispin, jspin, i,j, l, lo, jlo, nd
      real,allocatable:: rf(:)
      real :: ovlp
      call timestart("generate radial functions")

      !data is cached per itype only; call reset when the potential or energy parameters change
      if (this%itype /= itype .or. .not.allocated(this%r)) THEN
         if (input%l_useapw) call judft_bug("APW not implemented")
         if (any(atoms%l_dulo(:atoms%nlo(itype),itype))) call judft_bug("l_dulo not implemented in t_radfun")
         call this%init(atoms, input, itype)
         call usdus%init(atoms,input%jspins)
         nd = maxval(this%n_r)

         if (allocated(this%r)) deallocate (this%r, this%integral, this%bnd, this%e, this%ipred, this%l_sra)
         allocate (this%r( atoms%jmtd, 2,nd,0:atoms%lmaxd, input%jspins),source=0.0)
         allocate (this%integral(nd, nd,0:atoms%lmaxd, input%jspins,input%jspins),source=0.0)
         allocate (this%bnd(2,nd,0:atoms%lmaxd,input%jspins),source=0.0)
         allocate (this%e(nd,0:atoms%lmaxd,input%jspins),source=0.0)
         allocate (this%ipred(nd,nd,0:atoms%lmaxd,input%jspins),source=0.0)
         allocate (this%l_sra(nd,0:atoms%lmaxd),source=.true.)

         do ispin = 1, input%jspins
            call genMTBasis(atoms, enpara, vTot, fmpi, iType, ispin, usdus, f, g, flo, hub1data, l_writeArg=.false.)
            do l = 0, atoms%lmax(itype)
               this%R( 1:atoms%jri(itype), 1:2, 1,l, ispin) = f(1:atoms%jri(itype), 1:2, l)
               this%R( 1:atoms%jri(itype), 1:2, 2,l, ispin) = g(1:atoms%jri(itype), 1:2, l)
            end do
            do lo = 1, atoms%nlo(itype)
               this%R( 1:atoms%jri(itype), 1:2, atoms%slot_of_lo(lo,itype),atoms%llo(lo,itype), ispin) = flo(1:atoms%jri(itype), 1:2, lo)
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
                        this%integral(j,i,l,jspin,ispin)=ovlp
                     enddo
                  enddo
               ENDDO
            ENDDO
         ENDDO

         !Boundary values, energies and predecessor overlaps per slot
         do ispin = 1, input%jspins
            do l = 0, atoms%lmax(itype)
               this%bnd(:,1,l,ispin) = [usdus%us(l,itype,ispin), usdus%dus(l,itype,ispin)]
               this%bnd(:,2,l,ispin) = [usdus%uds(l,itype,ispin), usdus%duds(l,itype,ispin)]
               this%e(1:2,l,ispin) = enpara%el0(l,itype,ispin)
               this%ipred(:this%n_r(l),2,l,ispin) = this%integral(:this%n_r(l),1,l,ispin,ispin)
            end do
            do lo = 1, atoms%nlo(itype)
               l = atoms%llo(lo,itype)
               i = atoms%slot_of_lo(lo,itype)
               this%bnd(:,i,l,ispin) = [usdus%ulos(lo,itype,ispin), usdus%dulos(lo,itype,ispin)]
               this%e(i,l,ispin) = enpara%ello0(lo,itype,ispin)
               this%l_sra(i,l) = .not.atoms%l_relLO(lo,itype)
               if (atoms%ulo_der(lo,itype) < 1) cycle
               this%ipred(1,i,l,ispin) = usdus%uuilon(lo,itype,ispin)
               this%ipred(2,i,l,ispin) = usdus%duilon(lo,itype,ispin)
               do jlo = 1, atoms%nlo(itype)
                  if (atoms%llo(jlo,itype) /= l) cycle
                  this%ipred(atoms%slot_of_lo(jlo,itype),i,l,ispin) = usdus%ulouilopn(jlo,lo,itype,ispin)
               end do
            end do
         end do
      end if
      if (present(usdus_out)) then
         if (.not.allocated(usdus_out%us)) call usdus_out%init(atoms, input%jspins)
         call this%to_usdus(atoms, usdus_out)
      end if
      call timestop("generate radial functions")

   end subroutine

   function hsph(this, l, ispin) result(h)
      !! symmetrized spherical Hamiltonian <r_i|H_sph|r_j> in slot space.
      !! <H r_i|r_j> = e_i O_ij + <P r_i|r_j>; the hermitian part carries half the
      !! kinetic Wronskian W_ij on each side. For a pair with one non-SRA function
      !! (relLO) H is applied to the SRA partner only.
      class(t_radfun), intent(in) :: this
      integer, intent(in)         :: l, ispin
      real, allocatable           :: h(:,:)

      real, allocatable :: a(:,:), w(:,:), wi(:,:)
      integer :: n, i, j

      n = this%n_r(l)
      associate(o => this%integral(:n,:n,l,ispin,ispin), b => this%bnd(:,:n,l,ispin), sra => this%l_sra(:n,l))
         a = spread(this%e(:n,l,ispin),2,n)*o + transpose(this%ipred(:n,:n,l,ispin))
         w = -0.5*this%rmt**2*(spread(b(1,:),2,n)*spread(b(2,:),1,n) - spread(b(2,:),2,n)*spread(b(1,:),1,n))
         wi = reshape([((merge(0.5, merge(1.0,0.0,sra(i)), sra(i).eqv.sra(j)), i=1,n), j=1,n)], [n,n])
      end associate
      h = wi*(a + 0.5*w) + transpose(wi)*(transpose(a) - 0.5*w)
   end function

   subroutine to_usdus(this, atoms, usdus)
      !! fill the legacy t_usdus entries of this atom type
      use m_types_atoms
      use m_types_usdus
      class(t_radfun), intent(in)  :: this
      type(t_atoms), intent(in)    :: atoms
      type(t_usdus), intent(inout) :: usdus

      integer :: n, ispin, l, lo, jlo, i

      n = this%itype
      do ispin = 1, size(this%e,3)
         do l = 0, atoms%lmax(n)
            usdus%us(l,n,ispin)   = this%bnd(1,1,l,ispin)
            usdus%dus(l,n,ispin)  = this%bnd(2,1,l,ispin)
            usdus%uds(l,n,ispin)  = this%bnd(1,2,l,ispin)
            usdus%duds(l,n,ispin) = this%bnd(2,2,l,ispin)
            usdus%ddn(l,n,ispin)  = this%integral(2,2,l,ispin,ispin)
         end do
         do lo = 1, atoms%nlo(n)
            l = atoms%llo(lo,n)
            i = atoms%slot_of_lo(lo,n)
            usdus%ulos(lo,n,ispin)   = this%bnd(1,i,l,ispin)
            usdus%dulos(lo,n,ispin)  = this%bnd(2,i,l,ispin)
            usdus%uulon(lo,n,ispin)  = this%integral(1,i,l,ispin,ispin)
            usdus%dulon(lo,n,ispin)  = this%integral(2,i,l,ispin,ispin)
            usdus%uuilon(lo,n,ispin) = this%ipred(1,i,l,ispin)
            usdus%duilon(lo,n,ispin) = this%ipred(2,i,l,ispin)
            do jlo = 1, atoms%nlo(n)
               if (atoms%llo(jlo,n) == l) then
                  usdus%uloulopn(lo,jlo,n,ispin)  = this%integral(i,atoms%slot_of_lo(jlo,n),l,ispin,ispin)
                  usdus%ulouilopn(jlo,lo,n,ispin) = this%ipred(atoms%slot_of_lo(jlo,n),i,l,ispin)
               else
                  usdus%uloulopn(lo,jlo,n,ispin)  = 0.0
                  usdus%ulouilopn(jlo,lo,n,ispin) = 0.0
               end if
            end do
         end do
      end do
   end subroutine

   subroutine reset(this)
      !! invalidate cached radial functions
      class(t_radfun), intent(inout):: this
      this%itype = 0
      if (allocated(this%r)) deallocate(this%r, this%integral, this%bnd, this%e, this%ipred, this%l_sra)
   end subroutine

end module m_types_radfun
