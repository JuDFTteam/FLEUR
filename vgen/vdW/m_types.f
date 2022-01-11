      MODULE m_types
       IMPLICIT NONE
!
! Types for orbital moment calculation:
!
      TYPE t_orb                           ! 'normal' contributions
        SEQUENCE
        REAL :: uu,dd                      ! z   component
        COMPLEX :: uup,uum,ddp,ddm         ! +/- component
      END TYPE t_orb

      TYPE t_orbl                          ! local orbitals & (u,d)
        SEQUENCE
        REAL :: uulo,dulo
        COMPLEX :: uulop,uulom,dulop,dulom
      END TYPE t_orbl

      TYPE t_orblo                         ! lo,lo' contributions
        SEQUENCE
        REAL :: z
        COMPLEX :: p,m
      END TYPE t_orblo
!
! Types for spin-off-diagonal charge density:
!
      TYPE t_mt21                          ! 'normal' contributions
        SEQUENCE
        REAL ::  uun2,udn2,dun2,ddn2           ! normes of radial overlaps
        COMPLEX :: uu2,ud2,du2,dd2             ! values
      END TYPE t_mt21

      TYPE t_lo21                          ! local orbitals & (u,d)
        SEQUENCE
        REAL ::  uulon2,dulon2,uloun2,ulodn2   ! normes of radial overlaps
        COMPLEX :: uulo2,dulo2,ulou2,ulod2     ! values
      END TYPE t_lo21
!
! Type for LDA+U:
!
      TYPE t_utype
        SEQUENCE
        REAL u,j
        INTEGER l
        LOGICAL :: l_amf
      END TYPE t_utype
!
! Type for vacuum stars, etc.
!
      TYPE t_2dstar
        COMPLEX, ALLOCATABLE :: r2gphs(:,:)
        INTEGER, ALLOCATABLE :: i2g(:,:)
        INTEGER, ALLOCATABLE :: igvac(:,:)
        INTEGER, ALLOCATABLE :: ig2(:)
      END TYPE t_2dstar
!
! Type for the electric field
!
      TYPE t_efield
        REAL    :: zsigma  = 10.0  ! Distance to the charged plates
        REAL    :: sigma   =  0.0  ! charge at the plates
        REAL    :: sig_b(2)=  0.0  ! Extra charge for the top/bottom plate
        REAL,    ALLOCATABLE :: sigEF(:,:,:) ! (nx, ny, nvac)
        COMPLEX, ALLOCATABLE :: rhoEF(:,:)   ! (g_||, nvac)
        LOGICAL :: l_segmented = .false.
        LOGICAL :: plot_charge = .false. ! Plot charge as inputted
        LOGICAL :: plot_rho    = .false. ! Plot Fourier-transformed charge
        LOGICAL :: autocomp    = .true.  ! Auto-compensate film charge
      END TYPE t_efield
!
! type for wannier-functions
!
      TYPE t_wann
        integer :: wan90version
        integer :: oc_num_orbs
        integer,allocatable :: oc_orbs(:)
        logical :: l_unformatted
        logical :: l_oc_f
        logical :: l_ndegen
        logical :: l_orbitalmom
        logical :: l_orbcomp
        logical :: l_orbcomprs
        logical :: l_denmat
        logical :: l_perturbrs
        logical :: l_perturb
        logical :: l_nedrho
        logical :: l_anglmomrs
        logical :: l_anglmom
        logical :: l_spindisp
        logical :: l_spindisprs
        logical :: l_socspicom
        logical :: l_socspicomrs
        logical :: l_offdiposoprs
        logical :: l_offdiposop
        logical :: l_torque
        logical :: l_torquers
        logical :: l_atomlist
        integer :: atomlist_num
        integer,allocatable :: atomlist(:)
        logical :: l_berry
        logical :: l_perpmagrs
        logical :: l_perpmag
        logical :: l_perpmagat
        logical :: l_perpmagatrs
        logical :: l_socmatrs
        logical :: l_socmat
        logical :: l_soctomom
        logical :: l_kptsreduc2
        logical :: l_nablapaulirs
        logical :: l_nablars
        logical :: l_surfcurr
        logical :: l_updown
        logical :: l_ahe
        logical :: l_she
        logical :: l_rmat
        logical :: l_nabla
        logical :: l_socodi
        logical :: l_pauli
        logical :: l_pauliat
        logical :: l_potmat
        logical :: l_projgen
        logical :: l_plot_symm
        logical :: l_socmmn0
        logical :: l_bzsym
        logical :: l_hopping
        logical :: l_kptsreduc
        logical :: l_prepwan90
        logical :: l_plot_umdat
        logical :: l_wann_plot
        logical :: l_bynumber
        logical :: l_stopopt
        logical :: l_matrixmmn
        logical :: l_matrixamn
        logical :: l_projmethod
        logical :: l_wannierize
        logical :: l_plotw90
        logical :: l_byindex
        logical :: l_byenergy
        logical :: l_proj_plot
        logical :: l_bestproj
        logical :: l_ikptstart
        logical :: l_lapw
        logical :: l_plot_lapw
        logical :: l_fermi
        logical :: l_dipole
        logical :: l_dipole2
        logical :: l_dipole3
        logical :: l_mmn0
        logical :: l_mmn0at
        logical :: l_manyfiles
        logical :: l_collectmanyfiles
        logical :: l_ldauwan
        logical :: l_lapw_kpts
        logical :: l_lapw_gfleur
        logical :: l_kpointgen
        logical :: l_w90kpointgen
        integer :: ikptstart
        integer :: band_min(1:2)
        integer :: band_max(1:2)
        integer :: gfthick
        integer :: gfcut
        integer :: unigrid(6)
        integer :: mhp(3)
      END TYPE t_wann
!
! Types for vdW-forces
!
      TYPE :: sc_data
        REAL                :: alat
        REAL,DIMENSION(3,3) :: bravais
        REAL,DIMENSION(3,3) :: bravais_inv
      END TYPE sc_data
!
      TYPE atom_data
        INTEGER, DIMENSION(:),  POINTER :: nr_atom_type
        INTEGER, DIMENSION(:),  POINTER :: atomic_number
        REAL,    DIMENSION(:,:),POINTER :: coord_bravais
        REAL,    DIMENSION(:,:),POINTER :: coord_cart
      END TYPE atom_data

      END MODULE m_types
