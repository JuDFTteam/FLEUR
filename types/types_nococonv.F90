!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_nococonv
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  TYPE:: t_nococonv
    REAL   :: theta=0.0
    REAL   :: phi=0.0
    REAL   :: qss(3)=[0.,0.,0.]
    REAL, ALLOCATABLE    :: alph(:)
    REAL, ALLOCATABLE    :: beta(:)
    REAL, ALLOCATABLE    :: alphRlx(:)
    REAL, ALLOCATABLE    :: betaRlx(:)
    REAL, ALLOCATABLE    :: betaPrev(:)
    REAL, ALLOCATABLE    :: alphPrev(:)
    REAL, ALLOCATABLE    :: b_con(:,:)
  CONTAINS
    procedure:: init=>t_nococonv_init
    procedure:: init_ss=>t_nococonv_initss
    procedure:: chi
    procedure:: mpi_bc => nococonv_mpi_bc
  end TYPE
  public :: t_nococonv
CONTAINS
  function chi(nococonv,n)
    use m_constants
    IMPLICIT NONE
    CLASS(t_nococonv),INTENT(IN)  :: nococonv
    INTEGER,INTENT(IN)           :: n
    COMPLEX                      :: chi(2,2)

    chi(1,1) =  exp(ImagUnit*nococonv%alph(n)/2)*cos(nococonv%beta(n)/2)
    chi(2,1) = -EXP(-ImagUnit*nococonv%alph(n)/2)*SIN(nococonv%beta(n)/2)
    chi(1,2) = EXP(ImagUnit*nococonv%alph(n)/2)*SIN(nococonv%beta(n)/2)
    chi(2,2) =  EXP(-ImagUnit*nococonv%alph(n)/2)*COS(nococonv%beta(n)/2)
  end function

  subroutine t_nococonv_init(this,noco)
    use m_types_noco
    class(t_nococonv),INTENT(OUT):: This
    type(t_noco),INTENT(IN)      :: noco

    this%theta=noco%theta_inp
    this%phi=noco%phi_inp
    this%alph=noco%alph_inp
    this%beta=noco%beta_inp
    if (noco%l_ss) THEN
      this%qss=noco%qss_inp
    else
      this%qss=0.0
    endif
    if (allocated(this%b_con)) deallocate(this%b_con)
    allocate(this%b_con(2,size(this%alph)))
    this%b_con=0.0
    allocate(this%alphprev(size(this%alph)),this%betaprev(size(this%beta)))
  end subroutine

  subroutine t_nococonv_initss(nococonv,noco,atoms,qss)
      use m_types_noco
      use m_types_atoms
      use m_constants
      CLASS(t_nococonv),INTENT(inout):: nococonv
      TYPE(t_noco),INTENT(IN) :: noco
      TYPE(t_atoms),INTENT(IN):: atoms
      REAL,INTENT(IN),OPTIONAL :: qss(3)


      integer :: na,itype
      if (noco%l_ss) THEN
          nococonv%qss=noco%qss_inp
          if (present(qss)) nococonv%qss=qss
      endif
      ! Check noco stuff and calculate missing noco parameters
      IF (noco%l_noco) THEN
         IF (noco%l_ss) THEN
            !--->    the angle beta is relative to the spiral in a spin-spiral
            !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
            !--->    that means that the moments are "in line" with the spin-spiral
            !--->    (beta = qss * taual). note: this means that only atoms within
            !--->    a plane perpendicular to qss can be equivalent!
            na = 1
            DO iType = 1,atoms%ntype
               nococonv%alph(iType) = noco%alph_inp(iType) + tpi_const*dot_product(nococonv%qss,atoms%taual(:,na))
               na = na + atoms%neq(iType)
            END DO
         END IF
      ELSE

         IF (noco%l_ss) THEN
            CALL judft_error("l_noco=F and l_ss=T is meaningless.")
         END IF
      END IF
    end subroutine

    SUBROUTINE nococonv_mpi_bc(this, mpi_comm, irank)
      USE m_mpi_bc_tool
      CLASS(t_nococonv), INTENT(INOUT)::this
      INTEGER, INTENT(IN):: mpi_comm
      INTEGER, INTENT(IN), OPTIONAL::irank
      INTEGER ::rank,myrank,n,ierr
      IF (PRESENT(irank)) THEN
         rank = irank
      ELSE
         rank = 0
      END IF
      CALL mpi_bc(this%theta,rank,mpi_comm)
      CALL mpi_bc(this%phi,rank,mpi_comm)
      CALL mpi_bc(rank,mpi_comm,this%qss)
      CALL mpi_bc(this%alph,rank,mpi_comm)
      CALL mpi_bc(this%beta,rank,mpi_comm)
      CALL mpi_bc(this%alphRlx,rank,mpi_comm)
      CALL mpi_bc(this%betaRlx,rank,mpi_comm)
      CALL mpi_bc(this%alphPrev,rank,mpi_comm)
      CALL mpi_bc(this%betaPrev,rank,mpi_comm)
      CALL mpi_bc(this%b_con,rank,mpi_comm)

   END SUBROUTINE nococonv_mpi_bc

end module
