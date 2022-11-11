!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_nococonv
   USE m_judft
   Use m_constants

   IMPLICIT NONE
   PRIVATE
   TYPE:: t_nococonv
      REAL   :: theta = 0.0
      REAL   :: phi = 0.0
      REAL   :: qss(3) = [0., 0., 0.]
      REAL, ALLOCATABLE    :: alph(:)
      REAL, ALLOCATABLE    :: beta(:)
      REAL, ALLOCATABLE    :: alphRlx(:)
      REAL, ALLOCATABLE    :: betaRlx(:)
      REAL, ALLOCATABLE    :: betaPrev(:)
      REAL, ALLOCATABLE    :: alphPrev(:)
      REAL, ALLOCATABLE    :: b_con(:, :)
   CONTAINS
      procedure:: init => t_nococonv_init
      procedure:: init_ss => t_nococonv_initss
      !Routines to obtain chi transformation matrix
      procedure:: chi_pass
      procedure:: chi_explicit
      generic :: chi => chi_pass, chi_explicit
      generic :: umat => chi_pass, chi_explicit
      !Routines to rotate density matrix
      procedure:: rotdenmat_mat, rotdenmat_denmat
      procedure:: rotdenmat_explicit_mat, rotdenmat_explicit_denmat
      generic  :: rotdenmat => rotdenmat_mat, rotdenmat_denmat, rotdenmat_explicit_mat, rotdenmat_explicit_denmat
      !Functions to get magnetiszation vector from density matrix
      procedure :: denmat_to_mag_mat, denmat_to_mag_denmat
      generic   :: denmat_to_mag => denmat_to_mag_mat, denmat_to_mag_denmat
      !function to construct density matrix from magnetisaztion vector
      procedure:: mag_to_denmat
      !Rotate magnetisation vector
      procedure :: rot_magvec
      procedure :: avg_moments
      procedure :: mpi_bc => mpi_bc_nococonv
   end TYPE
   public :: t_nococonv
CONTAINS

SUBROUTINE mpi_bc_nococonv(this,mpi_comm,irank)
   USE m_mpi_bc_tool
   CLASS(t_nococonv),INTENT(INOUT)::this
   INTEGER,INTENT(IN):: mpi_comm
   INTEGER,INTENT(IN),OPTIONAL::irank
   INTEGER ::rank
   IF (PRESENT(irank)) THEN
      rank=irank
   ELSE
      rank=0
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

 END SUBROUTINE mpi_bc_nococonv

   function chi_pass(nococonv, n)
      CLASS(t_nococonv), INTENT(IN)  :: nococonv
      INTEGER, INTENT(IN)           :: n
      COMPLEX                      :: chi_pass(2, 2)
      chi_pass = nococonv%chi_explicit(nococonv%alph(n), nococonv%beta(n))
   end function

   pure function chi_explicit(nococonv, alpha, beta) result(chi)
      class(t_nococonv), intent(in) :: nococonv
      REAL, INTENT(IN) :: alpha, beta
      COMPLEX         :: chi(2, 2)
      chi(1, 1) =  EXP( ImagUnit*alpha/2)*COS(beta/2)
      chi(2, 1) = -EXP( ImagUnit*alpha/2)*SIN(beta/2)
      chi(1, 2) =  EXP(-ImagUnit*alpha/2)*SIN(beta/2)
      chi(2, 2) =  EXP(-ImagUnit*alpha/2)*COS(beta/2)
      chi=transpose(conjg(chi))
   end function

   function denmat_to_mag_mat(nococonv, mat) result(mag)
      class(t_nococonv), intent(in) :: nococonv
      complex, intent(in):: mat(2, 2)
      real :: mag(0:3)
      mag = nococonv%denmat_to_mag_denmat(real(mat(1, 1)), real(mat(2, 2)), mat(2, 1))
   end function

   function mag_to_denmat(nococonv, mag) result(mat)
      class(t_nococonv), intent(in) :: nococonv
      complex:: mat(2, 2)
      real, intent(in) :: mag(0:3)
      mat(1, 1) = 0.5*(mag(3) + mag(0))
      mat(2, 2) = 0.5*(mag(0) - mag(3))
      mat(2, 1) = cmplx(mag(1), mag(2))*0.5
      mat(1, 2) = cmplx(mag(1), -mag(2))*0.5
   end function

   function denmat_to_mag_denmat(nococonv, r11, r22, r21) result(mag)
      class(t_nococonv), intent(in) :: nococonv
      real, INTENT(IN)   :: r11, r22
      complex, intent(in):: r21
      real :: mag(0:3)
      mag(0) = r11 + r22
      mag(1) = 2*Real(r21)
      mag(2) = 2*Aimag(r21)
      mag(3) = r11 - r22
   end function

   subroutine rot_magvec(nococonv, n, mag, toGlobal)
      CLASS(t_nococonv), INTENT(IN) :: nococonv
      INTEGER, INTENT(IN)           :: n
      REAL, INTENT(INOUT)      :: mag(0:3)
      LOGICAL, INTENT(IN), OPTIONAL  :: toGlobal

      complex :: mat(2, 2)

      mat = nococonv%mag_to_denmat(mag)
      call nococonv%rotdenmat(n, mat, toGlobal)
      mag = nococonv%denmat_to_mag(mat)
   end subroutine

   subroutine rotdenmat_mat(nococonv, n, mat, toGlobal)
      CLASS(t_nococonv), INTENT(IN) :: nococonv
      INTEGER, INTENT(IN)           :: n
      COMPLEX, INTENT(INOUT) :: mat(2, 2)
      LOGICAL, INTENT(IN), OPTIONAL:: toGlobal

      real :: r11, r22
      r11 = real(mat(1, 1)); r22 = real(mat(2, 2))
      call nococonv%rotdenmat_explicit_denmat(nococonv%alph(n), nococonv%beta(n), r11, r22, mat(2, 1), toGlobal)
      mat(1, 1) = r11
      mat(2, 2) = r22
      mat(1, 2) = conjg(mat(2, 1))
   end subroutine

   subroutine rotdenmat_denmat(nococonv, n, rho11, rho22, rho21, toGlobal)
      CLASS(t_nococonv), INTENT(IN) :: nococonv
      INTEGER, INTENT(IN)           :: n
      REAL, INTENT(INOUT) :: rho11
      REAL, INTENT(INOUT) :: rho22
      COMPLEX, INTENT(INOUT) :: rho21
      LOGICAL, INTENT(IN), OPTIONAL:: toGlobal
      call nococonv%rotdenmat_explicit_denmat(nococonv%alph(n), nococonv%beta(n), rho11, rho22, rho21, toGlobal)
   end subroutine

   subroutine rotdenmat_explicit_mat(nococonv, alph, beta, mat, toGlobal)
      CLASS(t_nococonv), INTENT(IN) :: nococonv
      REAL, INTENT(IN) :: alph, beta
      COMPLEX, INTENT(INOUT) :: mat(2, 2)
      LOGICAL, INTENT(IN), OPTIONAL:: toGlobal
      real :: r11, r22
      r11 = real(mat(1, 1)); r22 = real(mat(2, 2))
      call nococonv%rotdenmat_explicit_denmat(alph, beta, r11, r22, mat(2, 1), toGlobal)
      mat(1, 1) = r11
      mat(2, 2) = r22
      mat(1, 2) = conjg(mat(2, 1))
   end subroutine

   SUBROUTINE rotdenmat_explicit_denmat(nococonv, alph, beta, rho11, rho22, rho21, toGlobal)
      use m_constants
      IMPLICIT NONE

      CLASS(t_nococonv), INTENT(IN) :: nococonv
      REAL, INTENT(IN) :: alph, beta
      REAL, INTENT(INOUT) :: rho11
      REAL, INTENT(INOUT) :: rho22
      COMPLEX, INTENT(INOUT) :: rho21
      LOGICAL, INTENT(IN), OPTIONAL:: toGlobal
      REAL r11n, r22n
      COMPLEX r21n
      if (present(toGlobal)) THEN
         if (toGlobal) THEN
            r11n = 0.5*(1.0 + cos(beta))*rho11 - sin(beta)*real(rho21) + 0.5*(1.0 - cos(beta))*rho22
            r22n = 0.5*(1.0 - cos(beta))*rho11 + sin(beta)*real(rho21) + 0.5*(1.0 + cos(beta))*rho22
            r21n = CMPLX(cos(alph), sin(alph))*(0.5*sin(beta)*(rho11 - rho22) + cos(beta)*real(rho21) + cmplx(0.0, aimag(rho21)))
            rho11 = r11n
            rho22 = r22n
            rho21 = r21n

            RETURN
         end if
      end if
      r11n = sin(beta)*(cos(alph)*real(rho21) + sin(alph)*AIMAG(rho21)) + (rho11 - rho22)*0.5*(1 + cos(beta)) + rho22
      r22n = -sin(beta)*(cos(alph)*real(rho21) + sin(alph)*AIMAG(rho21)) + (rho22 - rho11)*0.5*(1 + cos(beta)) + rho11
      r21n = (cos(alph)*real(rho21) + sin(alph)*AIMAG(rho21))*(1 + cos(beta)) - 0.5*sin(beta)*(rho11 - rho22) - cmplx(cos(alph), sin(alph))*conjg(rho21)
      rho11 = r11n
      rho22 = r22n
      rho21 = r21n

   end subroutine

   subroutine t_nococonv_init(this, noco)
      use m_types_noco
      class(t_nococonv), INTENT(OUT):: This
      type(t_noco), INTENT(IN)      :: noco

      this%theta = noco%theta_inp
      this%phi = noco%phi_inp
      this%alph = noco%alph_inp
      this%beta = noco%beta_inp
      if (noco%l_ss) THEN
         this%qss = noco%qss_inp
      else
         this%qss = 0.0
      end if
      if (allocated(this%b_con)) deallocate (this%b_con)
      allocate (this%b_con(2, size(this%alph)))
      this%b_con = 0.0
      allocate (this%alphprev(size(this%alph)), this%betaprev(size(this%beta)))

   end subroutine

   subroutine t_nococonv_initss(nococonv, noco, atoms, qss)
      use m_types_noco
      use m_types_atoms
      use m_constants
      CLASS(t_nococonv), INTENT(inout):: nococonv
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_atoms), INTENT(IN):: atoms
      REAL, INTENT(IN), OPTIONAL :: qss(3)

      integer :: na, itype
      if (noco%l_ss) THEN
         nococonv%qss = noco%qss_inp
         if (present(qss)) nococonv%qss = qss
      end if
      ! Check noco stuff and calculate missing noco parameters
      IF (noco%l_noco) THEN
         IF (noco%l_ss) THEN
            !--->    the angle beta is relative to the spiral in a spin-spiral
            !--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
            !--->    that means that the moments are "in line" with the spin-spiral
            !--->    (beta = qss * taual). note: this means that only atoms within
            !--->    a plane perpendicular to qss can be equivalent!
            na = 1
            DO iType = 1, atoms%ntype
               nococonv%alph(iType) = noco%alph_inp(iType) + tpi_const*dot_product(nococonv%qss, atoms%taual(:, na))
               na = na + atoms%neq(iType)
            END DO
         END IF
      ELSE

         IF (noco%l_ss) THEN
            CALL judft_error("l_noco=F and l_ss=T is meaningless.")
         END IF
      END IF
   end subroutine

   subroutine avg_moments(nococonv, den, atoms, magm, theta, phi)
      use m_types_atoms
      use m_types_potden
      use m_polangle
      use m_intgr
      class(t_nococonv), intent(in) :: nococonv
      class(t_potden), INTENT(IN):: den
      type(t_atoms), INTENT(IN)  :: atoms
      real, INTENT(OUT)          :: magm(3, atoms%ntype)
      real, INTENT(OUT), OPTIONAL          :: theta(atoms%ntype)
      real, INTENT(OUT), OPTIONAL          :: phi(atoms%ntype)

      integer:: i, j
      real:: integral(4)
      magm = 0.0
      DO i = 1, atoms%ntype
         integral = 0.0
         DO j = 1, size(den%mt, 4)
            call intgr3(den%mt(:, 0, i, j), atoms%rmsh(:, i), atoms%dx(i), atoms%jri(i), integral(j))
         END DO
         magm(3, i) = (integral(1) - integral(2))*sfp_const
         if (size(den%mt, 4) > 2) THEN
            magm(1, i) = -2*integral(3)*sfp_const
            magm(2, i) = 2*integral(4)*sfp_const
         end if
         if (present(theta)) THEN
            CALL pol_angle(magm(1, i), magm(2, i), magm(3, i), theta(i), phi(i), .true.)
         end if
      end do
   END subroutine

end module
