#define POTENTIAL
      MODULE param
      USE m_constants
      IMPLICIT NONE
!

      REAL,PARAMETER   :: au2A =0.529177249
      REAL,PARAMETER   :: Ha2eV=hartree_to_ev_const
      REAL,PARAMETER   :: pi=3.14159265358979323846
!
      INTEGER, PARAMETER   :: jnlout=6
      INTEGER, PARAMETER   :: ftable=20,fdebug=21,fxcrysd=22,fwarn=23
!
! parameters used to calculate the vdW-DF functional
!
      REAL,PARAMETER   ::  Zab_v1=-0.8491        ! vdW-DF1
      REAL,PARAMETER   ::  Zab_v2=-1.887         ! vdW-DF2
!
! parameters used to calculate the LDA correlation energy
!
      REAL,PARAMETER   ::  A       = 0.031091    ! some fitted parameters for the
      REAL,PARAMETER   ::  alpha_1 = 0.2137      ! evaluation of LDA correlation
      REAL,PARAMETER   ::  beta_1  = 7.5957      ! energy according to:
      REAL,PARAMETER   ::  beta_2  = 3.5876      !
      REAL,PARAMETER   ::  beta_3  = 1.6382      ! J. P. Perdew and Yue Wang,
      REAL,PARAMETER   ::  beta_4  = 0.49294     ! Phys. Rev. B 45, 13244 (1992)
      REAL,PARAMETER   ::  beta_H  = 0.066725    ! for H function from PRL 77, 3865 (1996)
      REAL,PARAMETER   ::  gamma_H = 0.031091    ! for H function
!
! parameters used to calculate the PBE exchange energy
!
      REAL,PARAMETER   ::  kappa_revPBE = 1.245   ! revPBE
      REAL,PARAMETER   ::  kappa_PBE    = 0.804   ! PBE from PRL 77, 3865 (1996)
      REAL,PARAMETER   ::  mu = 0.21951           ! PBE from PRL 77, 3865 (1996)
!
      END MODULE param
!
! module containing the main non-local variables and subroutines
!
      MODULE nonlocal_data

      IMPLICIT NONE
!
      REAL             ::  Zab             ! This is the Zab really used. Zab_v1
                                               ! is the default. Can be switched to
                                               ! Zab_v2 by invoking the program with
                                               ! command line option 'vdw2'
!
      INTEGER              :: nx, ny, nz       ! number of grid points in x, y, z direction
      INTEGER              :: n_grid           ! total number of grid points
      INTEGER              :: n_k              ! number of k-points for which the kernel was
                                               ! tabulated
!
      REAL             :: r_max            ! maximum r for which the kernel has been generated
!
      REAL,ALLOCATABLE :: q_alpha(:)       ! q mesh for the interpolation
      REAL,ALLOCATABLE :: phi(:,:,:)       !
      REAL,ALLOCATABLE :: d2phi_dk2(:,:,:) !
!
      INTEGER                 :: n_gvectors    ! number of G vectors
      REAL                :: G_cut         ! radius of cutoff sphere
      REAL,   ALLOCATABLE :: G(:,:)        ! the G vectors
      INTEGER,ALLOCATABLE :: G_ind(:)      ! the index of the G vectors
!
      REAL             :: a1(3),a2(3),a3(3)  ! lattice vectors
      REAL             :: b1(3),b2(3),b3(3)  ! reciprocal lattice vectors
!
      INTEGER              :: n_alpha            ! number of q points
      REAL             :: q_cut              ! maximum q value
      INTEGER              :: m_c                ! maximum m for the saturation of q_0
!
      REAL             :: omega              ! volume of the unit cell
      REAL             :: tpibya             ! (2*pi/omega)
!
      REAL             :: lambda             ! parameter for the logarithmic q mesh
      REAL             :: dk                 ! spacing of the uniform radial k grid
!
      INTEGER              :: time1(8),time2(8)  ! arrays for measuring the execution time
!
      END MODULE nonlocal_data
!
      MODULE driver_fft
!
 CONTAINS
!
!==========================================================================================
!
!     this subroutine is an interface for the fft. It will fourier transform in place the
!     complex array fftin. idir = -1 means forward and idir = +1 means backward fourier
!     transform. At the moment it is designed to use fftw3.

      SUBROUTINE inplfft( fftin, idir )
      USE m_juDFT
      USE param,        ONLY: jnlout
      USE nonlocal_data,ONLY: nx,ny,nz,n_grid

      implicit none


      complex, intent(inout) :: fftin(n_grid)

      integer                    :: idir

      integer*8                  :: plan
#ifdef CPP_FFTW
      include 'fftw3.f'  ! some definitions for fftw3

      if ( idir == -1 ) then

         CALL dfftw_plan_dft_3d(plan,nz,ny,nx,fftin,fftin,FFTW_FORWARD,FFTW_ESTIMATE)
         CALL dfftw_execute_dft(plan,fftin,fftin)
         CALL dfftw_destroy_plan(plan)

         fftin = fftin / real(n_grid) ! rescale as fftw3 puts a factor n_grid on forward FFT

      elseif( idir == 1 ) then

         CALL dfftw_plan_dft_3d(plan,nz,ny,nx,fftin,fftin,FFTW_BACKWARD,FFTW_ESTIMATE)
         CALL dfftw_execute_dft(plan,fftin,fftin)
         CALL dfftw_destroy_plan(plan)

      else

        write(jnlout,*) 'ERROR during FFT: neither FORWARD &
                         nor BACKWARD FFT was chosen.'
        STOP 'Error in FFT'

      endif
#else
	CALL judft_error("VdW functionals only available if compiled with CPP_FFTW")
#endif
      END SUBROUTINE inplfft
!
      END MODULE driver_fft
!
      MODULE functionals
!
 CONTAINS
!
!================================================================================
!
      SUBROUTINE calc_PBE_correlation(n,grad_n,Ec_PBE,Ec_LDA,e_cLDA,e_cSL)
!
      USE param,        ONLY: Ha2eV,pi,        &
                              jnlout,             &
                              beta_H,gamma_H,     &                        ! Parameters for LDA_c
                              A,alpha_1,beta_1,beta_2,beta_2,beta_3,beta_4 ! Parameters for LDA_c
      USE nonlocal_data,ONLY: n_grid,omega
!
      IMPLICIT NONE
!
      REAL,INTENT(in)  :: n(n_grid)          ! charge density
      REAL,INTENT(in)  :: grad_n(n_grid,3)   ! gradient of charge density
      REAL,INTENT(out) :: Ec_PBE             ! PBE-correlation energy
      REAL,INTENT(out) :: Ec_LDA             ! LDA-correlation energy
      REAL,INTENT(out) :: e_cLDA(n_grid)     ! LDA correlation energy density
      REAL,INTENT(out) :: e_cSL(n_grid)      ! SL correlation energy density in Ha
!
      REAL             :: r_s                ! intermediate values needed for the
      REAL             :: k_F                ! formulas given below.
      REAL             :: zeta               !
      REAL             :: phi_zeta           !
      REAL             :: grad_n_squ         !
      REAL             :: t                  !             -- "" --
      REAL             :: k_s                !
      REAL             :: AA                 !
      REAL             :: H                  !
!
      REAL             :: sqrt_r_s
      REAL             :: eLDA_c
      REAL             :: LDA_dummy_1
      REAL             :: LDA_dummy_2
!
      integer              :: i1
!
      e_cLDA(:)=0.0
      e_cSL(:) =0.0
!
      Ec_PBE=0.0
      Ec_LDA=0.0
!
!     Formula for PBE correlation ( Perdew, Burke and Ernzerhof PRL 77, 18 (1996) eqn. [3]):
!     Ec_PBE = int d^3 r n(r) (ec_LSDA (r_s, zeta) + H (r_s, zeta, t))
!
!     r_s = 3 / (4* pi*n)**(1/3) Seitz radius.
!     zeta is relative spin polarization will be set to 0 as vdW-DF only works without spin
!     t = | grad_n | / (2*phi(zeta)*k_s*n) is a dimensionless gradient
!     phi(zeta) = [ (1 + zeta)**(2/3) + (1 - zeta)**(2/3) ] / 2
!     k_s = sqrt( 4*k_F / (pi*a_0) )
!     a_0 = ( hbar / (m*e^2) ) = 1 in a. u.
!
!     H = (e**2/a_0) g * phi**3
!          * ln(1 + beta_H/gamma_H * t**2 * [1 + AA*t**2 / (1 + AA*t**2 + AA**2 * t**4 )] )
!     AA = beta_H/gamma_H * [ exp( -ec_LSDA / ( gamma_H * phi**3 * e**2 / a_0)) - 1 ]**(-1)
!
!
      zeta = 0.0     ! non-spin polarized case
!
      phi_zeta = ( (1.0 + zeta)**(2.0/3.0) + &
                   (1.0 - zeta)**(2.0/3.0))/2.0
!
      do i1=1,n_grid
!
         if(n(i1) < 1.0d-12) cycle

         r_s = (3.0/(4.0*pi*n(i1)))**(1.0/3.0)

!     we also need LDA correlation per particle so I repeat here the formulas from calc_q0

         sqrt_r_s = sqrt(r_s)

         LDA_dummy_1 = (8.0*pi/3.0)*A*(1.0+ alpha_1*r_s)
         LDA_dummy_2 =  2.0*A*(beta_1*sqrt_r_s     + beta_2*r_s + &
                                  beta_3*sqrt_r_s*r_s + beta_4*r_s*r_s)

         eLDA_c = (-3.0/(4.0*pi))*LDA_dummy_1*log(1.0+1.0/LDA_dummy_2)

!     now start calculation of PBE correlation. first check wether density is very small (approx
!     zero) or negative which is also not correct
!
         k_F = (3.0*pi*pi*n(i1))**(1.0/3.0)

         k_s = sqrt( 4.0 * k_F / pi)

         grad_n_squ = grad_n(i1,1)**2 + grad_n(i1,2)**2 + grad_n(i1,3)**2

         t = sqrt(grad_n_squ)/(2.0*phi_zeta*k_s*n(i1))

         AA = (beta_H/gamma_H)* &
              (1.0/(exp(-1.0*eLDA_c/(gamma_H*phi_zeta**3))- 1.0 ))

         H = gamma_H*phi_zeta**3 *         &
             log(1.0+                  &
                  (beta_H/gamma_H)*(t**2)* &
                  (1.0+AA*t**2)/(1.0+AA*t**2+AA**2 * t**4))
!
!     LDA and SL correlation energy density
!
         e_cLDA(i1) = eLDA_c*n(i1)
         e_cSL(i1)  =      H*n(i1)
!
#if 0
         Ec_PBE = Ec_PBE + n(i1)*eLDA_c + n(i1)*H
         Ec_LDA = Ec_LDA + n(i1)*eLDA_c
#else
         Ec_PBE = Ec_PBE + e_cLDA(i1) + e_cSL(i1)
         Ec_LDA = Ec_LDA + e_cLDA(i1)
#endif
      enddo
!
      Ec_PBE = Ec_PBE*omega/real(n_grid)
      Ec_LDA = Ec_LDA*omega/real(n_grid)
!
      END SUBROUTINE calc_PBE_correlation
!
!================================================================================
!
      SUBROUTINE calc_GGA_exchange(n,grad_n,Ex_PBE,Ex_revPBE,Ex_PW86,Ex_LDA)
!
      USE param,        ONLY: pi, &
                              mu,kappa_PBE,kappa_revPBE  ! Param. for PBE_ex
      USE nonlocal_data,ONLY: n_grid,omega
!
      IMPLICIT NONE
!
      REAL,INTENT(IN)  :: n(n_grid)          ! charge density
      REAL,INTENT(IN)  :: grad_n(n_grid,3)   ! gradient of charge density
      REAL,INTENT(OUT) :: Ex_PBE             !    PBE-exchange energy in Ha
      REAL,INTENT(OUT) :: Ex_revPBE          ! revPBE-exchange energy in Ha
      REAL,INTENT(OUT) :: Ex_PW86            ! refit PW86-exchange energy
                                                 ! in Ha
      REAL,INTENT(OUT) :: Ex_LDA             ! LDA part of the PBE exchange
                                                 ! energy Ex_PBE
!
      INTEGER              :: i1
      REAL             :: k_F                ! formulas given below.
      REAL             :: ss                 ! reduced density gradient
      REAL             :: sqrt_grad_n
!
! ss is actually the reduced density gradient
!                      ss = |\grad n|/2(3\pi^2)^1/3n^4/3
!
      Ex_PBE   =0.0
      Ex_revPBE=0.0
      Ex_PW86  =0.0
      Ex_LDA   =0.0
!
! here the charge density has already the correct indices
!
      do i1=1,n_grid
!
      if(n(i1) < 1.0d-12) cycle
!
        k_F = (3.0*pi*pi*n(i1))**(1.0/3.0)
!
        sqrt_grad_n = sqrt(grad_n(i1,1)**2 + grad_n(i1,2)**2 + grad_n(i1,3)**2)
!
        ss=sqrt_grad_n/(2.0*k_F*n(i1))
!
! LDA part of the PBE exchange energy density
!
        Ex_LDA = Ex_LDA +                                            &
                 (-3.0/4.0)*(3.0/pi)**(1.0/3.0)*      &
                 (n(i1)**(4.0/3.0))
!
! PBE exchange energy
!
        Ex_PBE = Ex_PBE +                                            &
                 (-3.0/4.0)*(3.0/pi)**(1.0/3.0)*      &
                 (n(i1)**(4.0/3.0))*                           &
                 (1.0+kappa_PBE-                                  &
                         kappa_PBE/(1.0+mu*ss**2.0/kappa_PBE))
!
! revPBE exchange energy
!
        Ex_revPBE = Ex_revPBE +                                      &
                 (-3.0/4.0)*(3.0/pi)**(1.0/3.0)*      &
                 (n(i1)**(4.0/3.0))*                           &
                 (1.0+kappa_revPBE-                               &
                         kappa_revPBE/(1.0+mu*ss**2.0/kappa_revPBE))
!
! use refit PW86
!
        Ex_PW86 = Ex_PW86 +                                          &
                 (-3.0/4.0)*(3.0/pi)**(1.0/3.0)*      &
                 (n(i1)**(4.0/3.0))*                           &
                 (1.0+15.0*0.1234*ss**2.0 +              &
                  17.33*ss**4.0+0.163*ss**6.0)**(1.0/15.0)
      enddo
!
      Ex_PBE    = Ex_PBE   *omega/real(n_grid)
      Ex_revPBE = Ex_revPBE*omega/real(n_grid)
      Ex_PW86   = Ex_PW86  *omega/real(n_grid)
      Ex_LDA    = Ex_LDA   *omega/real(n_grid)
!
      END SUBROUTINE calc_GGA_exchange
!
!==========================================================================================
!
      SUBROUTINE calc_ehartree(n)
      USE param,        ONLY: Ha2ev,pi,       &
                              jnlout
      USE nonlocal_data,ONLY: n_grid,n_gvectors, &
                              omega,             &
                              G,G_ind
      USE driver_fft
!
      IMPLICIT NONE
!
      REAL,INTENT(in)     :: n(n_grid)
!
      INTEGER                 :: i1,ind
      REAL                :: E_Hartree
      REAL                :: fac,size
      COMPLEX,ALLOCATABLE :: n_cmplx(:)
!
! E_Hartree= 2*pi*omega \sum_{G/=0} n(G)^2/G^2
!
      allocate(n_cmplx(n_grid))
!
      n_cmplx = (0.0,0.0)
      n_cmplx(1:n_grid) = cmplx( n(1:n_grid), 0.0 )
!
      call inplfft( n_cmplx, -1 )
!
      E_Hartree=0.0
!
      do i1 = 1, n_gvectors
!
      ind = G_ind(i1)
      size=sum(G(i1,1:3)**2)
      if (size > 1.d-12) then
         fac=1.0/size
      else
      print'(A,F10.3)','Total nr. of electrons: ',omega*real(n_cmplx(ind))
         fac=0.0
      endif
!
      E_Hartree=E_Hartree+ &
                (real(n_cmplx(ind))**2+aimag(n_cmplx(ind))**2)*fac
!
      enddo
!
      deallocate(n_cmplx)
!
      E_Hartree = E_Hartree*omega*2.0*pi
!
      write(jnlout,'(/,A)') '-------------- Hartree energy  --------------'
      write(jnlout,'(A)')      'E_Hartree (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') E_Hartree,E_Hartree*2.0,E_Hartree*Ha2eV
!
      END SUBROUTINE calc_ehartree
!
      END MODULE functionals
!
      MODULE plot_functions
!
 CONTAINS
!
!==========================================================================================
!
      SUBROUTINE write_xcrysden_LDA(file_name,ene_dens)
!
      USE param,        ONLY: au2A,fxcrysd
      USE nonlocal_data,ONLY: nx,ny,nz,n_grid, &
                              a1,a2,a3
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*),INTENT(IN) :: file_name
      REAL,        INTENT(IN) :: ene_dens(n_grid)
      INTEGER                     :: i1,ind,ix,iy,iz
      REAL,ALLOCATABLE        :: ene_temp(:,:,:)
!
      allocate(ene_temp(nx,ny,nz))
!
      open(fxcrysd,FILE=trim(file_name))
!
      write(fxcrysd,'(A)') 'CRYSTAL'
      write(fxcrysd,'(A)') 'PRIMVEC'
      write(fxcrysd,'(3F16.10)') a1(:)*au2A
      write(fxcrysd,'(3F16.10)') a2(:)*au2A
      write(fxcrysd,'(3F16.10)') a3(:)*au2A
!
      write(fxcrysd,'(A)')       'PRIMCOORD'
      write(fxcrysd,'(A,2x,I1)') 'XXX',1
      write(fxcrysd,'(A)')       'Atomic_positions'
!
      write(fxcrysd,'(A)') 'BEGIN_BLOCK_DATAGRID_3D'
      write(fxcrysd,'(A)') '3D_PWSCF'
      write(fxcrysd,'(A)') 'DATAGRID_3D_UNKNOWN'
      write(fxcrysd,'(3(3x,I5))') nx,ny,nz
      write(fxcrysd,'(3F16.10)') 0.0,0.0,0.0 ! origin of the system
      write(fxcrysd,'(3F16.10)') a1(:)*au2A
      write(fxcrysd,'(3F16.10)') a2(:)*au2A
      write(fxcrysd,'(3F16.10)') a3(:)*au2A
!
      do ix=1,nx
        do iy=1,ny
          do iz=1,nz
            ind = ix + nx*(iy - 1) + nx*ny*(iz - 1)
            ene_temp(ix,iy,iz)=ene_dens(ind)
            ! in the case of the semi-local corr. ene. dens.
            !   ene_temp(ix,iy,iz) can be positive
            !if(ene_temp(ix,iy,iz) > 0.0) then
            !    print'(A)','WARNING: ene_temp(ix,iy,iz) > 0.0!'
            !endif
          enddo
        enddo
      enddo
!
      do i1=1,nz
        write(fxcrysd,'(5(1x,E15.8))') ene_temp(1:nx,1:ny,i1)
      enddo
!
! write the end of the XSF file
!
      write(fxcrysd,'(A)') 'END_DATAGRID_3D'
      write(fxcrysd,'(A)') 'END_BLOCK_DATAGRID_3D'
!
      close(fxcrysd)
!
      deallocate(ene_temp)
!
      END SUBROUTINE write_xcrysden_LDA
!
!==========================================================================================
!
      SUBROUTINE write_xcrysden_NL(file_name,ene_dens)
!
      USE param,        ONLY: au2A,fxcrysd
      USE nonlocal_data,ONLY: nx,ny,nz,       &
                              n_grid,n_alpha, &
                              a1,a2,a3
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*),INTENT(IN) :: file_name
      REAL,        INTENT(IN) :: ene_dens(n_grid,n_alpha)
      INTEGER                     :: i1,ind,ix,iy,iz
      REAL,ALLOCATABLE        :: ene_temp(:,:,:)
!
      allocate(ene_temp(nx,ny,nz))
!
      open(fxcrysd,FILE=trim(file_name))
!
      write(fxcrysd,'(A)') 'CRYSTAL'
      write(fxcrysd,'(A)') 'PRIMVEC'
      write(fxcrysd,'(3F16.10)') a1(:)*au2A
      write(fxcrysd,'(3F16.10)') a2(:)*au2A
      write(fxcrysd,'(3F16.10)') a3(:)*au2A
!
      write(fxcrysd,'(A)')       'PRIMCOORD'
      write(fxcrysd,'(A,2x,I1)') 'XXX',1
      write(fxcrysd,'(A)')       'Atomic_positions'
!
      write(fxcrysd,'(A)') 'BEGIN_BLOCK_DATAGRID_3D'
      write(fxcrysd,'(A)') '3D_PWSCF'
      write(fxcrysd,'(A)') 'DATAGRID_3D_UNKNOWN'
      write(fxcrysd,'(3(3x,I5))') nx,ny,nz
      write(fxcrysd,'(3F16.10)') 0.0,0.0,0.0 ! origin of the system
      write(fxcrysd,'(3F16.10)') a1(:)*au2A
      write(fxcrysd,'(3F16.10)') a2(:)*au2A
      write(fxcrysd,'(3F16.10)') a3(:)*au2A
!
      do ix=1,nx
        do iy=1,ny
          do iz=1,nz
            ind = ix + nx*(iy - 1) + nx*ny*(iz - 1)
            ene_temp(ix,iy,iz)=sum(ene_dens(ind,1:n_alpha))
          enddo
        enddo
      enddo
!
      do i1=1,nz
        write(fxcrysd,'(5(1x,E15.8))') ene_temp(1:nx,1:ny,i1)
      enddo
!
! write the end of the XSF file
!
      write(fxcrysd,'(A)') 'END_DATAGRID_3D'
      write(fxcrysd,'(A)') 'END_BLOCK_DATAGRID_3D'
!
      close(fxcrysd)
!
      deallocate(ene_temp)
!
      END SUBROUTINE write_xcrysden_NL
!
      END MODULE plot_functions
!
      MODULE nonlocal_funct
      PRIVATE
      PUBLIC :: soler
!
 CONTAINS
!
!========================================================================================
!
      SUBROUTINE soler(n, Ecnl, v_nl)
      USE param,        ONLY: Ha2eV,au2A,jnlout
      USE nonlocal_data,ONLY: n_grid,n_alpha, &
                              omega,          &
                              q_alpha,phi,    &
                              G,G_ind,        &
                              time1,time2
      USE functionals
      USE driver_fft
      USE plot_functions
!
      IMPLICIT NONE
!
      REAL,          INTENT(INOUT) :: n(n_grid)       ! charge density
      REAL,          INTENT(INOUT) :: v_nl(n_grid)
      REAL                         :: Ecnl            ! the non-local correlation

      CHARACTER(LEN=132)               :: option
!
      REAL,ALLOCATABLE    :: grad_n(:,:)     ! gradient of charge density
      REAL,ALLOCATABLE    :: e_cLDA(:)       ! LDA correlation energy density
      REAL,ALLOCATABLE    :: e_cSL(:)        ! semi local correlation energy density
      REAL,ALLOCATABLE    :: e_cNL(:,:)      ! non local correlation energy density
      REAL,ALLOCATABLE    :: q_0(:)          ! array holding the q0 on the grid

#ifdef POTENTIAL

      REAL,ALLOCATABLE    :: dq0_dn(:)
      REAL,ALLOCATABLE    :: dq0_dgrad_n(:)

#endif

      COMPLEX,ALLOCATABLE :: theta_alpha(:,:)! array holding theta_i
      COMPLEX,ALLOCATABLE :: u_a(:,:)        ! quantity needed for the non local
                                                 ! correlation energy density
      REAL                :: Ecnl_rsp        ! Ecnl calculated in real space
      REAL                :: Ec_LDA_rsp      ! Ec_LDA calculated in real space
      REAL                :: Ec_LDA          ! LDA correlation
      REAL                :: Ec_PBE          ! PBE correlation energy. will be replaced by
                                                 ! LDA correlation + non local correlation
      REAL                :: Ex_PBE          ! PBE exchange energy
      REAL                :: Ex_revPBE       ! revPBE echange energy
      REAL                :: Ex_PW86         ! refit PW86 exchange energy
      REAL                :: Ex_LDA          ! LDA part of the PBE exchange
                                                 ! energy Ex_PBE
      INTEGER                 :: i, j, k
      INTEGER                 :: ind ,mem_size
      REAL                :: sum_test
!
!     the actual program starts here. allocation for n and grad_n have to be reconsidered when
!     putting this as a subroutine into a different program.
!
!     read the pregenerated kernel from file 'vdW_kernel_table'. This will set variables n_alpha
!     n_k and dk. The arrays phi(n_alpha,n_alpha,n_k) and d2phi_dk2(n_alpha,n_alpha,n_k) needed for the
!     later interpolation of the kernel will be allocated and read.
!
      call read_kernel()
!
      allocate(theta_alpha(n_grid,n_alpha))
      allocate(u_a(n_grid,n_alpha))
      allocate(grad_n(n_grid,3),q_0(n_grid))
      allocate(e_cLDA(n_grid),e_cSL(n_grid),e_cNL(n_grid,n_alpha))
!
      mem_size=2*n_grid*n_alpha + 2*n_grid*n_alpha + n_grid*3 + n_grid + &
               n_grid + n_grid + n_grid*n_alpha
      write(jnlout,'(A,F12.3,A,/)')  &
                   'Memory required: ',real(mem_size)*8.0/(1024.0*1024.0),' MB'
!
!     setup the G vectors which we need for the calculation of the gradient, the non local correlation
!     energy and the non local contribution to the potential.
!
      call DATE_AND_TIME(values=time1)
!
      call setup_g_vectors()
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time to setup G vectors:'
      call timing( time2 - time1 )
!
      write(jnlout,*)
!
!     calculate the gradient of the density numerically by FFT.
!
      call DATE_AND_TIME(values=time1)
!
      call calc_gradient( n, grad_n )
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time to calculate gradient of n numerically:'
      call timing( time2 - time1 )
!
!     calculate q_0 for every grid point according to equations (11),(12) from Dion and (7)
!     from Soler. The latter cares for the saturation.
!
      call DATE_AND_TIME(values=time1)
!
      call calc_q0(n, grad_n, q_0)
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time to calculate q0:'
      call timing( time2 - time1 )
!
!     calculate theta_i defined as theta_i = n(r_i) * p_alpha(q_i) where p_alpha are the
!     polynomials obtained by interpolating dirac delta. These polynomials will also be
!     derived in this subroutine
!
      call DATE_AND_TIME(values=time1)
!
      call calc_theta_i(n, q_alpha, q_0, theta_alpha)
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time to setup Thetas:'
      call timing( time2 - time1 )
!
!     to calculate the non local energy density we need theta_alpha_i in real space. As they are
!     overwritten during the fourier transform we save them here.
!
      e_cNL = theta_alpha
!
!     Fourier transform the theta_alpha_i in order to get theta_alpha_k.
!
      call DATE_AND_TIME(values=time1)
!
      do i=1,n_alpha

         call inplfft( theta_alpha(:,i), -1 )

      enddo
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time spent to Fourier transform Thetas:'
      call timing( time2 - time1 )
!
!     now we can evaluate the integral in eqn. (8) of Soler.
!
      call calc_ecnl(Ecnl, theta_alpha, u_a, e_cNL)
!
#ifdef POTENTIAL

      allocate(dq0_dn(n_grid))
      allocate(dq0_dgrad_n(n_grid))

      call calc_dq0(n, grad_n, dq0_dn, dq0_dgrad_n)

      call calc_potential(v_nl, u_a, q_0, dq0_dn, dq0_dgrad_n, n, grad_n)

      deallocate( dq0_dn, dq0_dgrad_n)

#endif

!     In vdW-DF GGA-correlation is been replaced by LDA-correlation + non local correlation. Thus
!     we have to calculate Ec_LDA and Ec_PBE now to output th energy term which has to be added
!     to the total energy. We already have the LDA correlation energy density from calc_q0 so we
!     just have to integrate that array and afterwards call a subroutine calculating the PBE.
!     omega/n_grid is the according volume element.
!
      call calc_PBE_correlation(n,grad_n,Ec_PBE,Ec_LDA,e_cLDA,e_cSL)
!
      write(jnlout,'(/,A)')     'PBE_c (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ec_PBE,Ec_PBE*2.0,Ec_PBE*Ha2eV
!
      write(jnlout,'(/,A)')     'LDA_c (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ec_LDA,Ec_LDA*2.0,Ec_LDA*Ha2eV
!
!     for testing write the sum over LDA energy density.
!
      Ec_LDA_rsp = sum(e_cLDA)*omega/real(n_grid)
      write(jnlout,'(/,A)')     &
                   '(check) LDA_c evaluated in real space (Ha/Ry):'
      write(jnlout,'(2F24.15)') Ec_LDA_rsp,Ec_LDA_rsp*2.0
!
      write(jnlout,'(/,A)')     'E_c,nl (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ecnl,Ecnl*2.0,Ecnl*Ha2eV
!
!     for testing write the sum over non local energy density.
!     it should be equal to Ec_nl
!
      Ecnl_rsp=sum(e_cNL)*omega/real(n_grid)
      write(jnlout,'(/,A)')     &
                   '(check) E_c,nl evaluated in real space (Ha/Ry):'
      write(jnlout,'(2F24.15)') 0.5*Ecnl_rsp,Ecnl_rsp
!
      call write_xcrysden_LDA('Ec_LDA.xsf',e_cLDA)
      call write_xcrysden_LDA('Ec_SL.xsf' ,e_cSL)
      call write_xcrysden_NL ('Ec_NL.xsf' ,e_cNL)
!
      call calc_GGA_exchange(n,grad_n,Ex_PBE,Ex_revPBE,Ex_PW86,Ex_LDA)
!
      write(jnlout,'(/,A)')     'Ex_LDA (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ex_LDA,Ex_LDA*2.0,Ex_LDA*Ha2eV
!
      write(jnlout,'(/,A)') '-------------- PBE exchange --------------'
      write(jnlout,'(A)')      'Ex_PBE (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ex_PBE,Ex_PBE*2.0,Ex_PBE*Ha2eV
!
      write(jnlout,'(A)') 'Ex_PBE+Ec_LDA+Ecnl (Ha/Ry/eV):'
      sum_test=Ex_PBE+Ec_LDA+Ecnl
      write(jnlout,'(3F24.15)') sum_test,sum_test*2.0,sum_test*Ha2eV
!
      write(jnlout,'(/,A)') '-------------- revPBE exchange --------------'
      write(jnlout,'(A)')      'Ex_revPBE (Ha/Ry/eV):'
      write(jnlout,'(3F24.15)') Ex_revPBE,Ex_revPBE*2.0,Ex_revPBE*Ha2eV
!
      write(jnlout,'(A)') 'Ex_revPBE+Ec_LDA+Ecnl (Ha/Ry/eV):'
      sum_test=Ex_revPBE+Ec_LDA+Ecnl
      write(jnlout,'(3F24.15)') sum_test,sum_test*2.0,sum_test*Ha2eV
!
      if (trim(option) == 'vdw2' .or. trim(option) == 'vdW2') then
!
! in the case of vdW-DF2 use refit PW86 exchange
!
        write(jnlout,'(/,A)') '-------------- for vdW-DF2 --------------'
        write(jnlout,'(A)') 'refit PW86 exchange:'
        write(jnlout,'(A)')      'Ex_PW86 (Ha/Ry/eV):'
        write(jnlout,'(3F24.15)') Ex_PW86,Ex_PW86*2.0,Ex_PW86*Ha2eV
!
        write(jnlout,'(A)') 'Ex_PW86+Ec_LDA+Ecnl (Ha/Ry/eV):'
        sum_test=Ex_PW86+Ec_LDA+Ecnl
        write(jnlout,'(3F24.15)') sum_test,sum_test*2.0,sum_test*Ha2eV
!
      endif
!
      call calc_ehartree(n)
!
      deallocate(grad_n,G,G_ind)
      deallocate(e_cNL,e_cSL,e_cLDA)
      deallocate(q_0,theta_alpha,u_a)
!
      END SUBROUTINE soler
!
!==========================================================================================
!
!     We take the kernel from the Thonhauser implementation (Espresso). This subroutine to
!     read the kernel thus is also taken from there.

      SUBROUTINE read_kernel()
      USE param,        ONLY: pi,ftable,jnlout
      USE nonlocal_data,ONLY: n_alpha,n_k,          &
                              r_max,q_cut,dk,       &
                              q_alpha,phi,d2phi_dk2
      implicit none

      character(len=30)    :: double_format = '(1p4e23.14)'
      integer              :: q1, q2
      if (allocated(q_alpha)) return
      write(jnlout,*) 'reading kernel table'

      open(ftable, file='vdW_kernel_table')

      read(ftable, '(2i5)' ) n_alpha, n_k
      read(ftable, double_format) r_max

      allocate(q_alpha(n_alpha), phi(0:n_k,n_alpha,n_alpha), d2phi_dk2(0:n_k,n_alpha,n_alpha))

      read(ftable, double_format) q_alpha

      q_cut = q_alpha(n_alpha)

      write(jnlout,*)'n_a:', n_alpha
      write(jnlout,*)'q_c:', q_cut
      write(jnlout,*)

      do q1 = 1, n_alpha
         do q2 = 1, q1

            read(ftable, double_format ) phi(0:n_k, q1, q2)
            phi(0:n_k,q2, q1) = phi(0:n_k, q1, q2)

         end do
      end do

      do q1 = 1, n_alpha
         do q2 = 1, q1

            read(ftable, double_format ) d2phi_dk2(0:n_k,q1, q2)
            d2phi_dk2(0:n_k,q2, q1) = d2phi_dk2(0:n_k,q1, q2)

         end do
      end do

      dk = 2.0*pi/r_max
!
      close(ftable)
!
      END SUBROUTINE read_kernel
!
!==========================================================================================
!
      SUBROUTINE setup_g_vectors()
      USE param,        ONLY: jnlout
      USE nonlocal_data,ONLY: nx,ny,nz,n_gvectors,n_grid, &
                              G_cut,b1,b2,b3,             &
                              G,G_ind
      implicit none
!
      integer             :: igx,igy,igz
      integer             :: tmpx,tmpy,tmpz
      integer             :: ind
!
      real            :: G_squ
      real            :: gx, gy, gz

!     We need the G vectors in a number of subroutines of this code. So the idea is to set them up
!     once and for all to minimize possible sources of errors. The array G(n_gvectors,3) will hold
!     the components of the g_vectors and G_ind their index on the fft mesh. We have to do the loop
!     over all G vectors twice. First time to count number of G_vectors inside cut off sphere and
!     second time to set them up.

      n_gvectors = 0

      do igz = -(nz-1)/2,(nz-1)/2
         do igy = -(ny-1)/2,(ny-1)/2
            do igx = -(nx-1)/2,(nx-1)/2

               gx = igx * b1(1) + igy * b2(1) + igz * b3(1)
               gy = igx * b1(2) + igy * b2(2) + igz * b3(2)
               gz = igx * b1(3) + igy * b2(3) + igz * b3(3)

               G_squ =  gx**2 + gy**2 + gz**2

               if ( G_squ .le. G_cut ) n_gvectors = n_gvectors + 1

            enddo
         enddo
      enddo
      if (allocated(g)) return
      allocate(G(n_gvectors, 3), G_ind(n_gvectors))

      write(jnlout,*) "#g-vectors:",n_gvectors," outof ",nz*ny*nx

      G     = 0.0
      G_ind = 0

      n_gvectors = 0

      do igz = -(nz-1)/2,(nz-1)/2
         do igy = -(ny-1)/2,(ny-1)/2
            do igx = -(nx-1)/2,(nx-1)/2

               gx = igx * b1(1) + igy * b2(1) + igz * b3(1)
               gy = igx * b1(2) + igy * b2(2) + igz * b3(2)
               gz = igx * b1(3) + igy * b2(3) + igz * b3(3)

               G_squ =  gx**2 + gy**2 + gz**2

               tmpx = 0
               tmpy = 0
               tmpz = 0

               IF (igx .LT. 0) tmpx = nx
               IF (igy .LT. 0) tmpy = ny
               IF (igz .LT. 0) tmpz = nz

               ind = 1 + ( igx + tmpx ) + &
                      nx*( igy + tmpy ) + &
                   nx*ny*( igz + tmpz )

               if ( G_squ .le. G_cut ) then

                   n_gvectors = n_gvectors + 1

                   G(n_gvectors, 1) = gx
                   G(n_gvectors, 2) = gy
                   G(n_gvectors, 3) = gz

                   G_ind(n_gvectors) = ind

               endif

            enddo
         enddo
      enddo

      end subroutine setup_g_vectors
!
!==========================================================================================
!
      SUBROUTINE calc_gradient(n,grad_n)
      USE param,        ONLY: jnlout
      USE nonlocal_data,ONLY: n_grid,n_gvectors,nx,ny,nz, &
                              a1,a2,a3,G,G_ind
      USE driver_fft
!
      implicit none
!
      real, intent(in)     :: n(n_grid)
      real, intent(inout)  :: grad_n(n_grid,3)

      complex, allocatable :: n_cmplx(:)
      complex, allocatable :: grad_n_cmplx(:,:)

      integer                  :: i

      integer                  :: ind

      integer*8                :: plan

      grad_n = 0.0

!     A derivative d/dx f(x) in real space corresponds to i*G*f(G) in reciprocal space. So we have
!     to fourier transform the density, find the according G vectors, calculate the product and back
!     fourier transform the gradient

      allocate(n_cmplx(n_grid), grad_n_cmplx(n_grid,3))

      n_cmplx      = (0.0, 0.0)
      grad_n_cmplx = (0.0, 0.0)

      do i = 1, n_grid

         n_cmplx(i) = cmplx( n(i), 0.0)

      enddo

      call inplfft( n_cmplx, -1 )

      do i = 1, n_gvectors

         ind = G_ind(i)

         grad_n_cmplx(ind,:) = (0.0,1.0) * G(i,:) * n_cmplx(ind)

      enddo

      do i=1,3

        call inplfft( grad_n_cmplx(:,i), 1)

      enddo
!
      grad_n(:,:) = real(grad_n_cmplx(:,:))

      deallocate(n_cmplx, grad_n_cmplx)
!
      END SUBROUTINE calc_gradient
!
!==========================================================================================
!
!     This subroutine calculates q_0 following eqns. 11 and 12 of Dion and saturates it
!     following eqn. 7 of Soler
!
      SUBROUTINE calc_q0(n, grad_n, q_0)
      USE param,        ONLY: pi, &
                              A,alpha_1,beta_1,beta_2,beta_2,beta_3,beta_4 ! Parameters for LDA_c
      USE nonlocal_data,ONLY: n_grid,q_cut,m_c,Zab
!
      implicit none
!
      real, intent(in)    ::  n(n_grid)        ! charge density
      real, intent(in)    ::  grad_n(n_grid,3) ! gradient of charge density
      real, intent(out)   ::  q_0(n_grid)      ! array holding g_0 on the grid

      real                ::  q                ! q before saturation

      real                ::  k_F                ! fermi wave vector
      real                ::  r_s                !
      real                ::  sqrt_r_s           !
      real                ::  eLDA_c             ! LDA correlation
      real                ::  eLDA_x             ! LDA exchange
      real                ::  grad_n_squ         ! squared gradient of n

      real                ::  LDA_dummy_1
      real                ::  LDA_dummy_2

      real                ::  sum_m

      integer                 ::  i, m         ! counters
!
!     loop over all grid points
!
      q_0(:) = q_cut
!
      do i=1,n_grid

!     check if charge density is negative. If so treat it like zero charge density. Zero
!     density corresponds to high q_0. Thats why we set it to q_cut the highest possible
!     q_0 and continue with the next grid point.

         if ( n(i) < 1.d-12 ) then
            q_0(i) = q_cut
            cycle
         end if
!
!     calculate q according to eqns (11) and (12) from Dion.
!
         k_F = (3.0*(pi**2.0)*n(i))**(1.0/3.0)
         r_s = (3.0/(4.0*pi*n(i)))**(1.0/3.0)
         sqrt_r_s = r_s**(1.0/2.0)

         grad_n_squ = grad_n(i,1)**2.0 + grad_n(i,2)**2.0 + grad_n(i,3)**2.0

         LDA_dummy_1 = (8.0*pi/3.0)*A*(1.0+ alpha_1*r_s)
         LDA_dummy_2 =  2.0*A*(beta_1*sqrt_r_s     + beta_2*r_s + &
                                  beta_3*sqrt_r_s*r_s + beta_4*r_s*r_s)

         eLDA_x = (-3.0/(4.0*pi))*k_F
         eLDA_c = (-3.0/(4.0*pi))*LDA_dummy_1*log(1.0+1.0/LDA_dummy_2)
!
!     LDA correlation energy density is epsilon_c[n(r)]*n(r)
!
!         e_cLDA(i) = eLDA_c*n(i)
!
         q = (1.0 + ( eLDA_c / eLDA_x) - (Zab/9.0)*grad_n_squ &
                     /(4.0*(n(i)**2.0)*(k_F**2.0)))*k_F

!     now saturate q according to eq. (7) from Soler.

         sum_m = 0.0

         do m=1,m_c

            sum_m = sum_m + (q/q_cut)**m / real(m)

         enddo

         q_0(i) = q_cut*(1.0 - exp(-sum_m))

      enddo

      END SUBROUTINE calc_q0
!
!==========================================================================================
!
!
      SUBROUTINE calc_theta_i(n, q_alpha, q_0, theta_alpha)

      USE nonlocal_data,ONLY: n_grid, n_alpha
!
      implicit none
!
      real, intent(in)      ::  n(n_grid)        ! charge density
      real, intent(in)      ::  q_alpha(n_alpha) ! q-mesh for interpolation
      real, intent(in)      ::  q_0(n_grid)      ! q values on the grid

      complex, intent(inout) ::  theta_alpha(n_grid,n_alpha) ! thetas on real space grid

      integer                    ::  i               ! counter

!     this call will interpolate dirac delta which gives p_alpha(q_i) according to a method
!     from numerical recipes. Whithin this subroutine also the initial setup of the second
!     derivatives needed for spline interpolation will be done once.

      call splint( q_alpha, q_0, theta_alpha )

      do i=1,n_grid

!     theta_alpha will hold at this stage the p_alpha(q_i)

            theta_alpha(i,:) = n(i)*theta_alpha(i,:)

      enddo

      END SUBROUTINE calc_theta_i
!
!==========================================================================================
!
      SUBROUTINE calc_ecnl(Ecnl, theta_alpha, u_a, e_cNL)
      USE param        ,ONLY: pi,jnlout
      USE nonlocal_data,ONLY: nx,ny,nz,n_grid,n_alpha,n_gvectors, &
                              omega,time1,time2,                  &
                              G, G_ind
      USE driver_fft
!
      implicit none
!
      real,   intent(out)   :: Ecnl
      complex,intent(in)    :: theta_alpha(n_grid,n_alpha)
      complex,intent(out)   :: u_a(n_grid,n_alpha)   ! function u_a(r) needed for the construction
                                                         ! of the potential and also for non local
      real,   intent(inout) :: e_cNL(n_grid,n_alpha) ! non local correlation energy density                                                      ! correlation energy density.
!
      integer                 :: i, ind
      integer                 :: alpha,beta
      real                :: k
      real                :: phi_k(n_alpha,n_alpha)
      complex             :: integral
!
      u_a = (0.0, 0.0)
      Ecnl = 0.0
!
!     integration of theta_a*theta_b*phi_ab
!
      write (jnlout,*) 'calculating E_c,nl:'
!
      call DATE_AND_TIME(values=time1)
!
!     The difference to older versions of this code is the way how the G-vectors are set up. The
!     index is not any more calculated on the fly but stored in G_ind.
!
      do i = 1, n_gvectors

          ind = G_ind(i)

          k = sqrt(dble( G(i,1)**2 + G(i,2)**2 + G(i,3)**2))

          call interpolate_kernel(k, phi_k)

          integral = (0.0, 0.0)

          do alpha = 1,n_alpha
             do beta = 1,n_alpha

                integral = integral + conjg(theta_alpha(ind,alpha)) * &
                    theta_alpha(ind,beta)*cmplx(phi_k(alpha,beta),0.0) !&
!
!    the array u_a(k,a) = sum_b theta_b(k) * phi_ab(k) equals the fourier transform
!    of the function u_a(r) = sum_b int d^3 r' theta_b(r') phi_ab(r - r') which we need for
!    the calculation of the potential and the non local correlation energy density.
!
                u_a(ind,alpha) = u_a(ind,alpha) + &
                    theta_alpha(ind,beta)*cmplx(phi_k(alpha,beta),0.0)

             enddo
          enddo

          Ecnl = Ecnl + real(integral)

      enddo
!
!     backward fourier transform of u_a(k,a) in order to get u_a(r) the non local correlation
!     energy density.
!
      do i=1,n_alpha

         call inplfft( u_a(:,i), 1 )

      enddo
!
!     at this moment e_cNL(i,alpha) holds the theta_alpha_i which we have to multiply by the convolution
!     u_a(r) in order to get the non local energy density as a function of alpha
!
      e_cNL = e_cNL*real(u_a)
!
      call DATE_AND_TIME(values=time2)
!
      write(jnlout,'(A)') 'time for integration:'
      call timing( time2 - time1 )
!
!     taking care of the units. (2*pi)^3/omega is the volume element of the reciprocal unit
!     cell. omega^2 we get when replacing the integrals by sums over a real space grid. The
!     (2*pi)^3 cancels with an according factor from the radial fourier transform on the kernel
!
      Ecnl = Ecnl*0.5*omega
!
      END SUBROUTINE calc_ecnl
!
!==========================================================================================
!
      SUBROUTINE make_q_alpha(q_alpha)


      USE nonlocal_data,ONLY: n_alpha,q_cut,lambda ! n_alpha: number of q points, q_cut: maximum
                                                   ! q value, lambda: parameter for logarithmic
                                                   ! mesh.
      implicit none
!
      real, intent(out) ::  q_alpha(n_alpha) ! array holding the qs

      real              ::  q1               ! auxiliary vairable holding the first
                                                 ! q which can be calculated explicitely

      integer               ::  i                ! counter

!     first we have to setup the qs starting from the maximum value q_cut on a logarithmic mesh
!     which satisfies ( q_a+1 - q_a ) = lambda*( q_a - q_a-1 ). Lambda > 1 is a parameter.
!     In GPAW lambda = 1.2 is used.

      q1 = q_cut * ( lambda - 1.0) / (lambda**( n_alpha - 1.0 ) - 1.0 )

      do i = 1,n_alpha

         q_alpha(i) = q1 * (lambda**( i - 1.0 ) - 1.0 ) / ( lambda - 1.0 )

      enddo

      END SUBROUTINE make_q_alpha
!
!==========================================================================================
!
!     From Numerical Recipes

      SUBROUTINE splint( x_i, x, p_iofx  )


      USE nonlocal_data,ONLY: n_grid, n_alpha
!
      implicit none
!
      real, intent(in)  ::  x_i(n_alpha) ! the x_i values where the function is known
      real, intent(in)  ::  x(n_grid)    ! the x values for which the function
                                             ! shall be interpolated

      complex, intent(inout) ::  p_iofx(n_grid,n_alpha) ! the function values for each x

      real, allocatable :: d2y_dx2(:,:) ! second derivatives needed for the splines

      real, allocatable :: y(:)

      real              :: a, b, c, d, h ! some intermediate variables for the interpolation

      integer               :: i, j               ! counters
      integer               :: u_bound, l_bound   ! variables for finding the section of x_i
                                                  ! in which x is located.
      integer               :: alpha              ! index of the found section

      allocate(y(n_alpha))

      allocate(d2y_dx2(n_alpha,n_alpha))

      call setup_spline(x_i,d2y_dx2)

      do i=1,n_grid

!     first find the correct section of x_i in which x is located by bisectioning. According to
!     numerical recipes this is efficient if the x values are random. In our case there might
!     be some correlation as the density and its gradient are smooth.

         l_bound = 1
         u_bound = n_alpha

         do while ((u_bound - l_bound) > 1)

            alpha = (u_bound + l_bound) / 2

            if ( x(i) > x_i(alpha)) then
               l_bound = alpha
            else
               u_bound = alpha
            endif
         enddo

         h = x_i(u_bound) - x_i(l_bound)

         a = (x_i(u_bound) - x(i))/h
         b = (x(i) - x_i(l_bound))/h
         c = ((a**3.0 - a)*h**2.0)/6.0
         d = ((b**3.0 - b)*h**2.0)/6.0

         do j=1,n_alpha

            y    = 0.0
            y(j) = 1.0

            p_iofx(i,j) = a*y(l_bound) + b*y(u_bound) + &
               c*d2y_dx2(j,l_bound) + d*d2y_dx2(j,u_bound)

         enddo
      enddo

      deallocate(y, d2y_dx2)

      END SUBROUTINE splint
!
!==========================================================================================
!
!     From Numerical Recipes

      SUBROUTINE setup_spline(x_i,d2y_dx2)


      USE nonlocal_data,ONLY: n_alpha
!
      implicit none
!
      real, intent(in)     ::  x_i(n_alpha)
      real, intent(inout)  ::  d2y_dx2(n_alpha,n_alpha)

      real, allocatable    ::  y(:)       ! this array holds the function values at x_i which
                                              ! are going to be interpolated.
      real, allocatable    ::  dy_dx(:)   ! temporary array holding the first derivative

      real                 ::  sig, p     ! temporary values for the interpolation

      integer                  ::  i, j       ! counter

      allocate(y(n_alpha), dy_dx(n_alpha))

      do i=1,n_alpha

!     as we are interpolating dirac delta set the according function values:

         y=0.0
         y(i)=1.0

!     now according to numerical recipes we set the boundary conditions which will give
!     so called "natural" splines.

         d2y_dx2(i,1) = 0.0
         dy_dx(1) = 0.0

         do j=2,n_alpha-1

            sig = (x_i(j) - x_i(j-1)) / (x_i(j+1) - x_i(j-1))
            p = sig*d2y_dx2(i,j-1) + 2.0
            d2y_dx2(i,j) = (sig - 1.0)/p
            dy_dx(j) = (y(j+1) - y(j))/(x_i(j+1) - x_i(j)) - (y(j) - y(j-1))/(x_i(j) - x_i(j-1))
            dy_dx(j) = (6.0*dy_dx(j)/(x_i(j+1) - x_i(j-1)) - sig*dy_dx(j-1))/p

         enddo

         d2y_dx2(i,n_alpha) = 0.0

         do j=n_alpha - 1, 1, -1

            d2y_dx2(i,j) = d2y_dx2(i,j)*d2y_dx2(i,j+1) + dy_dx(j)

         enddo
      enddo

      deallocate(y, dy_dx)

      END SUBROUTINE setup_spline
!
!===========================================================================================
!
!     similar to the Thonhauser implementation

      SUBROUTINE interpolate_kernel(k, phi_k)


      USE nonlocal_data,ONLY: n_alpha, dk, n_k, phi, d2phi_dk2
!
      implicit none
!
      real, intent(in)    :: k
      real, intent(inout) :: phi_k(n_alpha,n_alpha)

      real                :: a, b, c, d

      integer                 :: q1, q2, k_i

!     the algorithm for interpolation will be more or less the same as in splint().
      phi_k = 0.0

!     find the interval in which our k lies

      k_i = int(k/dk)

!     simple case when k equals one of the values we have tabulated the kernel for

      if (mod(k,dk) == 0.0) then

         do q1=1,n_alpha
            do q2=1,q1

               phi_k(q1, q2) = phi(k_i, q1, q2)
               phi_k(q2, q1) = phi(k_i, q2, q1)

            enddo
         enddo

         return

      endif

      a = (dk*(k_i+1.0) - k)/dk
      b = (k - dk*k_i)/dk
      c = (a**3.0-a)*dk**2.0/6.0
      d = (b**3.0-b)*dk**2.0/6.0

      do q1 = 1, n_alpha
         do q2 = 1, q1

            phi_k(q1, q2) = a*phi(k_i, q1, q2) + b*phi(k_i+1, q1, q2) &
            +(c*d2phi_dk2(k_i, q1, q2) + d*d2phi_dk2(k_i+1,q1, q2))

            phi_k(q2, q1) = phi_k(q1, q2)

         end do
      end do

      END SUBROUTINE interpolate_kernel
!
!=================================================================================================
!
#ifdef POTENTIAL

      subroutine calc_dq0(n, grad_n, dq0_dn, dq0_dgrad_n)

      USE param,        ONLY: pi, &
                              A,alpha_1,beta_1,beta_2,beta_2,beta_3,beta_4 ! Parameters for LDA_c
      USE nonlocal_data,ONLY: n_grid,q_cut,m_c,Zab


      implicit none

      real, intent(in)    ::  n(n_grid)        ! charge density
      real, intent(in)    ::  grad_n(n_grid,3) ! gradient of charge density

      real, intent(out)   ::  dq0_dn(n_grid)   ! derivative of q0 w.r.t. the density
      real, intent(out)   ::  dq0_dgrad_n(n_grid) ! derivative of q0 w.r.t. the gradient

      real                ::  q                ! q before saturation

      real                ::  dq0_dq           ! derivative needed for dq0_dn and dq0_dgrad_n

      real                ::  k_F                ! fermi wave vector
      real                ::  r_s                !
      real                ::  sqrt_r_s           !
      real                ::  eLDA_c             ! LDA correlation
      real                ::  eLDA_x             ! LDA exchange
      real                ::  grad_n_squ         ! squared gradient of n

      real                ::  LDA_dummy_1
      real                ::  LDA_dummy_2

      real                ::  sum_m

      integer                 ::  i, m         ! counters

!     loop over all grid points

      dq0_dn(:) = 0.0
      dq0_dgrad_n(:) = 0.0

      do i=1,n_grid

         if ( n(i) < 1.d-12 ) then
            cycle     ! NOTE: The derivatives will be zero at these points
         end if
!
!     calculate q according to eqns (11) and (12) from Dion.
!
         k_F = (3.0*(pi**2.0)*n(i))**(1.0/3.0)
         r_s = (3.0/(4.0*pi*n(i)))**(1.0/3.0)
         sqrt_r_s = r_s**(1.0/2.0)

         grad_n_squ = grad_n(i,1)**2.0 + grad_n(i,2)**2.0 + grad_n(i,3)**2.0

         LDA_dummy_1 = (8.0*pi/3.0)*A*(1.0+ alpha_1*r_s)
         LDA_dummy_2 =  2.0*A*(beta_1*sqrt_r_s     + beta_2*r_s + &
                                  beta_3*sqrt_r_s*r_s + beta_4*r_s*r_s)

         eLDA_x = (-3.0/(4.0*pi))*k_F
         eLDA_c = (-3.0/(4.0*pi))*LDA_dummy_1*log(1.0+1.0/LDA_dummy_2)

         q = (1.0 + ( eLDA_c / eLDA_x) - (Zab/9.0)*grad_n_squ &
                     /(4.0*(n(i)**2.0)*(k_F**2.0)))*k_F

!     now saturate q according to eq. (7) from Soler.

         sum_m = 0.0
         dq0_dq = 0.0

         do m=1,m_c

            sum_m = sum_m + (q/q_cut)**m / real(m)

            dq0_dq = dq0_dq + ((q/q_cut)**(m-1))

         enddo

         dq0_dq = dq0_dq * exp(-sum_m)

!     here we calculate the derivatives needed for the potential later. i.e. we calculate
!     n*dq0/dn and n*dq0/dgrad_n so that we do not have to divide by n which might be very
!     small at some points and thus produce numerical errors.

         dq0_dn(i) = dq0_dq * ( k_F/3.0 &
                     + k_F*7.0/3.0*(Zab/9.0)*grad_n_squ &
                     /(4.0*(n(i)**2.0)*(k_F**2.0)) &
                     - (8.0*pi/9.0)*A*alpha_1*r_s*log(1.0+1.0/LDA_dummy_2) &
                     + LDA_dummy_1/(LDA_dummy_2*(1.0 + LDA_dummy_2)) &
                     * (2.0*A*(beta_1/6.0*sqrt_r_s + beta_2/3.0*r_s &
                     + beta_3/2.0*r_s*sqrt_r_s + 2.0*beta_4/3.0*r_s**2)))

         dq0_dgrad_n(i) = dq0_dq * 2.0 * (-1.0*Zab/9.0) &
                          * sqrt(grad_n_squ) / (4.0*n(i)*k_F)

      enddo

      end subroutine calc_dq0
!
!=======================================================================================
!
!     similar to the Thonhauser implementation

      subroutine calc_potential(v_nl, u_a, q_0, dq0_dn, dq0_dgrad_n, n, grad_n)


      use nonlocal_data, only : n_grid, n_alpha, q_alpha, q_cut, nx, ny, nz, G, G_ind, n_gvectors

      USE driver_fft

      implicit none

      real, intent(out)    :: v_nl(n_grid) ! non local part of the potential

      real, intent(in)     :: q_0(n_grid)

      real, intent(in)     :: dq0_dn(n_grid)      ! n*dq0/dn
      real, intent(in)     :: dq0_dgrad_n(n_grid) ! n*dq0/dgrad_n

      real, intent(in)     :: n(n_grid)
      real, intent(in)     :: grad_n(n_grid,3)

      complex, intent(in)  :: u_a(n_grid, n_alpha)

      complex, allocatable :: h_j(:,:)

      real, allocatable    :: d2y_dx2(:,:) ! second derivatives needed for the splines

      real, allocatable    :: y(:)

      real                 :: a, b, c, d, d1, d2, h ! some intermediate variables
                                                        ! for the interpolation of p_a

      real                 :: p_a
      real                 :: dpa_dq0

      real                 :: grad_n_abs ! | grad_n |

      integer                  :: i, alpha

      integer                  :: l_bound, u_bound

      integer                  :: ind

      allocate(h_j(n_grid,3))
      allocate(y(n_alpha))

      h_j = (0.0, 0.0)

      allocate( d2y_dx2(n_alpha, n_alpha) )

      call setup_spline( q_alpha, d2y_dx2 )

!     The potential will be calculated using FFT following White and Bird. PRB 50, 4954
!     (1994) eqn. (10). v^xc(r) = d f_xc / d n(r) - (1/N) * sum_G,r' i * G * ( grad_n(r')
!     / |grad_n(r')| ) * ( d f_xc / d |grad_n(r')| ) * e^(i*G*(r - r')) with E_c,NL =
!     int d^3 r f_xc .

      do i=1,n_grid

!     first we have to interpolate the polynomials p_a and their derivatives with respect
!     to q0 again as we need them for the calculation of the potential. This will be the
!     same procedure as in splint.

         l_bound = 1
         u_bound = n_alpha

         do while ((u_bound - l_bound) > 1)

            ind = (u_bound + l_bound) / 2

            if ( q_0(i) > q_alpha(ind)) then
               l_bound = ind
            else
               u_bound = ind
            endif
         enddo

         h = q_alpha(u_bound) - q_alpha(l_bound)

         a = (q_alpha(u_bound) - q_0(i))/h
         b = (q_0(i) - q_alpha(l_bound))/h
         c = ((a**3.0 - a)*h**2.0)/6.0
         d = ((b**3.0 - b)*h**2.0)/6.0
         d1 = (3.0*a**2.0 - 1.0)*h/6.0
         d2 = (3.0*b**2.0 - 1.0)*h/6.0

         do alpha = 1, n_alpha

            y = 0.0
            y(alpha) = 1.0

            p_a = a*y(l_bound) + b*y(u_bound) + &
                  c*d2y_dx2(alpha,l_bound) + d*d2y_dx2(alpha,u_bound)

            dpa_dq0 = (y(u_bound) - y(l_bound))/h &
                      - d1*d2y_dx2(alpha,l_bound) + d2*d2y_dx2(alpha,u_bound)

!     first term sum_a u_ai * ( p_ai + n_i * dpai/dqi * dqi/dni ). the factor n(i) is already
!     contained in dq0_dn

            v_nl(i) = v_nl(i) + u_a(i,alpha) * (p_a + dpa_dq0 * dq0_dn(i)  )

!     the following sum we need for h_j which will be fourier transformed later. The IF
!     condition excludes the contributions of high q values.

            grad_n_abs = sqrt(grad_n(i,1)**2.0 + grad_n(i,2)**2.0 + grad_n(i,3)**2.0)

            if ( q_0(i) .ne. q_cut ) then

               h_j(i,:) = h_j(i,:) + u_a(i,alpha)*cmplx((grad_n(i,:) / grad_n_abs) &
                                   * dpa_dq0 * dq0_dgrad_n(i), 0.0)

            endif
         enddo
      enddo

!     now fourier transform h_j and carry out the sum over G.

      do i = 1,3

          call inplfft( h_j(:,i), -1 )

      enddo

      do i = 1, n_gvectors

         ind = G_ind(i)

         h_j(ind,:) = (0.0,1.0) * G(i,:) * h_j(ind,:)

      enddo

!     back fourier transform h_j and add it to the potential

      do i = 1,3

         call inplfft( h_j(:,i), 1 )

      enddo

      v_nl(:) = v_nl(:) - real( h_j(:,1) + h_j(:,2) + h_j(:,3) )

      deallocate(h_j, y,d2y_dx2)

      end subroutine
!
!==========================================================================================
!
#endif

      SUBROUTINE timing ( time )
      USE param,ONLY: jnlout
!
      implicit none
!
      integer, intent(in)     :: time(8)
      integer                 :: tmp(8), time_out(8)

      integer                 :: i

!     In the older versions of this program the time written to output file could be negative. This
!     subroutine will fix this problem.

      tmp = 0

      tmp(8) = 1000
      tmp(7) = 60
      tmp(6) = 60
      tmp(5) = 24

      do i=1,8

         if (time(i) .lt. 0) then
            time_out(i) = time(i) + tmp(i)
            time_out(i-1) = time(i-1) - 1
         else
            time_out(i) = time(i)
         endif

      enddo
!
      write(jnlout, '(I4,A,I4,A,I4,A,I4,A)') time_out(5)," h  ", &
                                             time_out(6)," min", &
                                             time_out(7)," s  ", &
                                             time_out(8)," ms "
!
      END SUBROUTINE timing
!
!     As the hole thing is intended to be interfaced with different codes I put here just an
!     example program for testing, which is capable of reading an input file, reading a
!     charge density, which otherwise has to be given by the calling program, and then call
!     the subroutine soler(n, ... ).
!
      END MODULE nonlocal_funct
