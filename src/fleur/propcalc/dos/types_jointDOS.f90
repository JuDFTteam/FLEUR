!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_jointdos
   use m_judft
   use m_types_eigdos
   implicit none
   PRIVATE
   public t_jointdos
   TYPE, extends(t_eigdos):: t_jointDOS
      REAL, ALLOCATABLE :: qis(:, :, :)  !i,kpt,spin
      REAL, ALLOCATABLE :: qal(:, :, :, :, :) !l,ndos,i,kpt,spin
      REAL, ALLOCATABLE :: qTot(:, :, :)  !i,kpt,spin
      CHARACTER(LEN=20), ALLOCATABLE :: weight_names(:)

   CONTAINS
      PROCEDURE, PASS :: init => jointDOS_init
      PROCEDURE      :: get_weight_eig
      PROCEDURE      :: get_num_weights
      PROCEDURE      :: get_weight_name
      procedure      :: postprocessing
   END TYPE t_jointDOS
CONTAINS
   !NOTE: charge_mag is a module procedure on purpose. It is only used by
   !      postprocessing and would naturally be an internal procedure there,
   !      but nvfortran (23.9) cannot compile this module if postprocessing has
   !      a CONTAINS section: it emits the internal subprogram and the static
   !      block holding the host association under the same mangled name and
   !      then dies with "cannot be a common block and a subprogram", or, when
   !      -g is added, with "Internal compiler error. flowgraph: node is zero".
   !      Do not turn this back into an internal procedure.
   function charge_mag(vec1,vec2)
      real :: charge_mag(2)
      real, intent(in):: vec1(:)
      real, intent(in):: vec2(:)
      real :: rho1,mz1,rho2,mz2
      real, parameter :: eps_rho = 1.0e-10

      !distribution into charge and magnetisation parts (mz only, other components are directly given in density matrix vector), factor 1/2 included later
      rho1=vec1(1)+vec1(2)
      mz1=vec1(1)-vec1(2)
      rho2=vec2(1)+vec2(2)
      mz2=vec2(1)-vec2(2)
      if (abs(rho1) < eps_rho.or. abs(rho2) < eps_rho) THEN
         charge_mag(1)=0.0
         charge_mag(2)=0.0
         return
      endif
      charge_mag(1)= 0.25*rho1*rho2  !charge part
      charge_mag(2)= 0.25*mz1*mz2+vec1(3)*vec2(3)+vec1(4)*vec2(4) !mag part
      !Convert to sum and difference
      rho1=charge_mag(1)
      charge_mag(1) = 0.5*(rho1+charge_mag(2))
      charge_mag(2) = 0.5*(rho1-charge_mag(2))
   end function charge_mag

   subroutine postprocessing(this, noco,nococonv, banddos, alldos, ef)
      use m_types_atoms
      use m_types_noco
      use m_types_nococonv
      use m_types_banddos
      use m_types_dos
      use m_sort 
      class(t_jointDOS), intent(inout):: this
      class(t_eigdos_list), intent(in),optional    :: alldos(:)
      TYPE(t_noco), INTENT(IN)        :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_banddos), INTENT(IN)    :: banddos
      real, intent(in),optional       :: ef
      type(t_dos),pointer :: dos
      integer :: ikpt, ispin, jspin, ii, ntype, l, i, j, iispin,n,n_dos
      integer :: idx(size(this%eig,1))
      !find a DOS of type t_dos from the given alldos for the postprocessing
      n_dos=0
      NULLIFY(dos)
      DO n=1,size(alldos)
         if (.not. associated(alldos(n)%p)) cycle
         select type(d=>alldos(n)%p)
            type is (t_dos)
               dos=>d
               n_dos=n
               exit
         end select
      end do
      if (n_dos==0) then
         call judft_error("No eigdos of type t_dos found for jointDOS postprocessing")
      end if

      !calculate the eigdos from the DOS
      DO ikpt=1,size(dos%qis,2)
         if (size(dos%qis,3)<3) then !collinear case
            DO ispin=1,size(dos%qis,3)
               DO jspin=1,ispin
                  ii=0
                  DO i=1,size(dos%qis,1)
                     if (dos%eig(i,ikpt,ispin) < ef) then  !valid initial state
                        DO j=1,size(dos%qis,1)
                           if (dos%eig(j,ikpt,jspin) > ef) then !valid final state
                              ii=ii+1
                              if (ispin==jspin) then
                                 iispin=ispin
                              else
                                 iispin=3
                              end if   
                              this%eig(ii,ikpt,iispin)=this%eig(ii,ikpt,iispin)+dos%eig(j,ikpt,jspin)-dos%eig(i,ikpt,ispin)
                              !add contribution to jointDOS
                              this%qis(ii,ikpt,iispin)=this%qis(ii,ikpt,iispin)+dos%qis(i,ikpt,ispin)*dos%qis(j,ikpt,jspin)
                              this%qTot(ii,ikpt,iispin)=this%qTot(ii,ikpt,iispin)+dos%qTot(i,ikpt,ispin)*dos%qTot(j,ikpt,jspin)
                              DO ntype=1,size(banddos%dos_typelist)
                                 DO l=0,3
                                    this%qal(l,ntype,ii,ikpt,iispin)=this%qal(l,ntype,ii,ikpt,iispin)+&
                                                      dos%qal(l,ntype,i,ikpt,ispin)*dos%qal(l,ntype,j,ikpt,jspin)
                                 ENDDO
                              ENDDO
                           endif
                        enddo
                     endif  
                  enddo
                  !Sort the results according to excitation energy
                  CALL sort(idx(:ii),this%eig(:ii,ikpt,iispin))
                  this%eig(:ii,ikpt,iispin)=this%eig(idx(:ii),ikpt,iispin)
                  this%qis(:ii,ikpt,iispin)=this%qis(idx(:ii),ikpt,iispin)
                  this%qTot(:ii,ikpt,iispin)=this%qTot(idx(:ii),ikpt,iispin)
                  DO ntype=1,size(banddos%dos_typelist)
                     DO l=0,3
                        this%qal(l,ntype,:ii,ikpt,iispin)=this%qal(l,ntype,idx(:ii),ikpt,iispin)
                     ENDDO
                  ENDDO
               enddo   
            enddo   
         else  !noncollinear case
            ii=0
            DO i=1,size(dos%qis,1)  
               if (dos%eig(i,ikpt,1) < ef) then !valid initial state
                  DO j=i,size(dos%qis,1)
                     if (dos%eig(j,ikpt,1) >ef) then !valid final state
                        ii=ii+1
                        this%eig(ii,ikpt,:)=dos%eig(j,ikpt,1)-dos%eig(i,ikpt,1)
                        print *, "jDOS:",ii,"=", i, "->", j,this%eig(ii,ikpt,1)
                        !add contribution to jointDOS (charge and magnetisation components)
                        this%qis(ii,ikpt,1:2)=this%qis(ii,ikpt,1:2)+charge_mag(dos%qis(i,ikpt,:),dos%qis(j,ikpt,:))
                        this%qTot(ii,ikpt,1:2)=this%qTot(ii,ikpt,1:2)+charge_mag(dos%qTot(i,ikpt,:),dos%qTot(j,ikpt,:))
                        DO ntype=1,size(banddos%dos_typelist)
                           DO l=0,3
                              this%qal(l,ntype,ii,ikpt,1:2)=this%qal(l,ntype,ii,ikpt,1:2)+&
                                                charge_mag(dos%qal(l,ntype,i,ikpt,:),dos%qal(l,ntype,j,ikpt,:))
                           ENDDO
                        ENDDO
                     endif
                  enddo
               end if
            ENDDO
            !Sort the results according to excitation energy
            CALL sort(idx(:ii),this%eig(:ii,ikpt,1))
            DO iispin=1,size(this%eig,3)
               this%eig(:ii,ikpt,iispin)=this%eig(idx(:ii),ikpt,iispin)
               this%qis(:ii,ikpt,iispin)=this%qis(idx(:ii),ikpt,iispin)
               this%qTot(:ii,ikpt,iispin)=this%qTot(idx(:ii),ikpt,iispin)
               DO ntype=1,size(banddos%dos_typelist)
                  DO l=0,3
                     this%qal(l,ntype,:ii,ikpt,iispin)=this%qal(l,ntype,idx(:ii),ikpt,iispin)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
            
      ENDDO

   end subroutine postprocessing
   
   
   integer function get_num_weights(this)
      class(t_jointDOS), intent(in):: this
      get_num_weights = 0
      if (allocated(this%weight_names)) get_num_weights = size(this%weight_names)
   end function

   character(len=20) function get_weight_name(this, id) !TODO
      class(t_jointDOS), intent(in):: this
      INTEGER, intent(in)         :: id
      if (.not. allocated(this%weight_names)) call judft_error("No weight names in t_jointDOS")
      if (id > size(this%weight_names)) call judft_error("Not enough weight names in t_jointDOS")
      get_weight_name = this%weight_names(id)
   end function

   integer function get_num_spins(this)
      class(t_jointDOS), intent(in):: this
      get_num_spins = size(this%qis, 3)
   end function

   function get_weight_eig(this, id)
      class(t_jointDOS), intent(in):: this
      INTEGER, intent(in)     :: id
      real, allocatable:: get_weight_eig(:, :, :)

      INTEGER :: ind, l, ntype, i
      allocate (get_weight_eig, mold=this%qis)

      if (id == 1) THEN
         get_weight_eig = this%qTot
         if (all(this%qis == 0.0)) then
            get_weight_eig = 1.0
         END IF
      END IF
      if (id == 2) THEN
         get_weight_eig = this%qis
         if (all(get_weight_eig == 0.0)) then
            !No INT dos calculated so far...
            get_weight_eig = 1.0
            DO ntype = 1, size(this%qal, 2)
               DO l = 0, 3
                  get_weight_eig = get_weight_eig - this%qal(l, ntype, :, :, :)
               END DO
            END DO
         end if
      end if
      ind = 2
      DO ntype = 1, size(this%qal, 2)
         DO l = 0, 3
            ind = ind + 1
            if (ind == id) get_weight_eig = this%qal(l, ntype, :, :, :)
         END DO
      END DO
   end function

   SUBROUTINE jointDOS_init(thisDOS, input, atoms, kpts, banddos, l_noco, eig)
      USE m_types_input
      USE m_types_atoms
      USE m_types_banddos
      USE m_types_kpts
      USE m_types_noco
      IMPLICIT NONE
      CLASS(t_jointDOS), INTENT(INOUT) :: thisDOS
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_atoms), INTENT(IN)    :: atoms
      LOGICAL, INTENT(IN)          :: l_noco
      TYPE(t_kpts), INTENT(IN)    :: kpts
      TYPE(t_banddos), INTENT(IN)    :: banddos
      real, intent(in)                       :: eig(:, :, :)

      INTEGER :: ntype, l, i, ind,ispin
      character :: spdfg(0:4) = ["s", "p", "d", "f", "g"]
      thisDOS%name_of_dos = "jointDOS"
      
      ispin= merge(2,input%jspins,l_noco)
      ALLOCATE (thisDOS%qis((input%neig*input%neig)/4+1, kpts%nkpt, ispin))
      ALLOCATE (thisDOS%qal(0:3, size(banddos%dos_typelist), (input%neig*input%neig)/4+1, kpts%nkpt,ispin))
      ALLOCATE (thisDOS%qTot((input%neig*input%neig)/4+1, kpts%nkpt, ispin))
      ALLOCATE (thisDOS%eig((input%neig*input%neig)/4+1, kpts%nkpt, ispin))
      thisDOS%qis = 0.0
      thisDOS%qal = 0.0
      thisDOS%qTot = 0.0

      allocate (thisDOS%weight_names(2 + 4*size(banddos%dos_typelist)))
      thisDOS%weight_names(1) = "Total"
      thisDOS%weight_names(2) = "INT"
      ind = 2
      DO ntype = 1, size(banddos%dos_typelist)
         DO l = 0, 3
            ind = ind + 1
            write (thisDOS%weight_names(ind), "(a,i0,a)") "MT:", banddos%dos_typelist(ntype), spdfg(l)
         END DO
      END DO

   END SUBROUTINE jointDOS_init
end MODULE m_types_jointdos
