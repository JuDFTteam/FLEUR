!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! module m_jpSternheimer: Self-consistent solution of the Sternheimer cycle
!
!> @author
!> Christian-Roman Gerhorst
!>
!> @brief
!> Self-consistent solution of the Sternheimer cycle. This task is operated by the subroutine m_jpsternheimer::solvesternheimerscc.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpSternheimer.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpSternheimer

    USE m_constants

  implicit none

  contains

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> This routine initializes all quantities which are required for the Sternheimer SCC and can be shifted from the scope of the
  !> k-point as well as atom loops of the Sternheimer SCC.
  !>
  !> @details
  !> A shifted G-set fulfilling |G + q | <= Gmax, the gradient of the full charge density, the functional derivatives of the
  !> xc-kernel, at the moment only x-alpha, the gradients of the external and effective potential, with and without the term
  !> containing the gradient of the density, which is canceling within the Sternheimer SCC, the t-matrices of the non-spherical
  !> potential for the Pulay contributions to the Sternheimer SCC and the gradient of the Gauss-curve pseudo core-density in the
  !> displaced atom as well as the interstitial are determined. See also Section 7.4.5 (dissertation CRG).
  !>
  !> @param[in]   atoms            : Atoms type, see types.f90
  !> @param[in]   sym              : Symmetries type, see types.f90
  !> @param[in]   cell             : Unit cell type, see types.f90
  !> @param[in]   stars            : Stars type, see types.f90
  !> @param[in]   dimens           : Dimension type, see types.f90
  !> @param[in]   lathar           : Lattice harmonics type, see types.f90
  !> @param[in]   enpara           : Energy parameter type, see types.f90
  !> @param[in]   uds              : Type containing quantities consisting of the radial solutions, see types.f90
  !> @param[in]   input            : Input type, see types.f90
  !> @param[in]   qpts             : Q-point set represented by a k-points type, see types.f90
  !> @param[out]  tdHS0            : Tlmplm matrix type for Sternheimer Pulay matrix elements, see types.f90
  !> @param[in]   logUnit          : Unit number for juPhon.log
  !> @param[in]   ngdp             : Number of G-vectors for potentials and densities
  !> @param[in]   iqpt             : Index of q-point in kpts file which is evaluated ! DEPRECATED shifted G-set
  !> @param[in]   gdp              : G-vectors of potentials and densities
  !> @param[in]   rho0IR           : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]   rho0MT           : Lattice harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[in]   mlh_atom         : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[in]   nmem_atom        : Number of lattice harmonic members for every atom.
  !> @param[in]   clnu_atom        : Phase mediating between stars and plane waves.
  !> @param[in]   rbas1            : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]   rbas2            : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]   uuilon           : Overlap integral between the radial functions of the integral (multiplied by ulo_der) of a
  !>                                 local orbital and the flapw radial function with the same l
  !> @param[in]   duilon           : Overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                                 derivative of the flapw radial function with the same l
  !> @param[in]   ulouilopn        : Overlap integral between the radial functions of the integral of a local orbital and another
  !>                                 local orbital with the same l.
  !> @param[in]   ilo2p            : Mapping array giving the p value for given number of LO and itype
  !> @param[in]   rho0IRpw         : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]   rho0MTsh         : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from
  !>                                 Fleur
  !> @param[in]   vEff0MTsh        : Spherical harmonic coefficients of unperturbed and converged muffin-tin effective potential
  !>                                 parsed from Fleur
  !> @param[out]  qpwcG            : Plane-wave expansion coefficients of the FFT of the gradient of the pseudo core-density ( Gauss
  !>                                 curve )
  !> @param[out]  rho1MTCoreDispAt : Spherical-harmonic expansion coefficients of the gradient of the pseudo core-density ( Gauss
  !>                                 curve ) at the displaced atom.
  !> @param[out]  grVxcIRKern      : Plane-wave interstitial coefficients of the functional derivative of the xc-kernel with respect
  !>                                 to the density
  !> @param[out]  dKernMTGPts      : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                                 respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[out]  grVeff0MT_init   : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed external
  !>                                 potential (part of the effective) for initial Sternheimer iteration created in
  !>                                 m_jpGrVeff0::gengrveff0
  !> @param[out]  grVeff0MT_main   : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective
  !>                                 potential for a regular and the final Sternheimer iteration created in m_jpGrVeff0::gengrveff0
  !> @param[out]  grVext0IR_DM     : Plane-wave interstitial coefficients of the gradient of the unperturbed external potential for
  !>                                 the dynamical matrix created in m_jpGrVeff0::gengrveff0
  !> @param[out]  grVext0MT_DM     : Spherical harmonic muffin-tin coefficients of the complete gradient of the unperturbed external
  !>                                 potential for the dynamical matrix created in m_jpGrVeff0::gengrveff0
  !> @param[out]  grVeff0IR_DM     : Plane-wave interstitial coefficients of the gradient of the unperturbed effective potential for
  !>                                 the dynamical matrix created in m_jpGrVeff0::gengrveff0
  !> @param[out]  grVeff0MT_DM     : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective
  !>                                 potential for the dynamical matrix created in m_jpGrVeff0::gengrveff0
  !> @param[out]  grRho0IR         : Plane-wave expansion coefficients of the gradient of the charge density
  !> @param[out]  grRho0MT         : Muffin-tin expansion coefficients of the gradient of the charge density created in
  !>                                 mod_juphonutils::calcgrr2finlh
  !> @param[out]  gausWts          : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[out]  ylm              : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to
  !>                                 lmax
  !>
  !> @todo Insert spin, replace the 1
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
!!! Removed initSCC !!!
  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Main subroutine of this module. It operates the self-consistent solution of the Sternheimer cycle/equations.
  !>
  !> @param[in]   atoms            : Atoms type, see types.f90
  !> @param[in]   sym              : Symmetries type, see types.f90
  !> @param[in]   stars            : Stars type, see types.f90
  !> @param[in]   lathar           : Lattice harmonics type, see types.f90
  !> @param[in]   dimens           : Dimension type, see types.f90
  !> @param[in]   cell             : Unit cell type, see types.f90
  !> @param[in]   enpara           : Energy parameter type, see types.f90
  !> @param[in]   uds              : Type containing quantities consisting of the radial solutions, see types.f90
  !> @param[in]   input            : Input type, see types.f90
  !> @param[in]   kpts             : K-points type, see types.f90
  !> @param[in]   qpts             : Q-point set represented by a k-points type, see types.f90
  !> @param[in]   results          : Results type, see types.f90
  !> @param[in]   usdus            : Type containing quantities consisting of the radial solutions, see types.f90
  !> @param[in]   tdHS0            : Tlmplm matrix type for Sternheimer Pulay matrix elements, see types.f90
  !> @param[in]   logUnit          : Unit number for juPhon.log
  !> @param[in]   ngdp             : Number of G-vectors for potentials and densities
  !> @param[in]   iqpt             : Index of q-point in kpts file which is evaluated
  !> @param[in]   oneSternhCycle   : Switch to start the Sternheimer SCC from the final iteration (An emtpy file start_cdn1 should
  !>                                 exist, so that the converged density is read in)
  !> @param[in]   noPts2chCont     : Number of points for which continuity test should be performed.
  !> @param[in]   rbas1            : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]   rbas2            : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]   El               : Contains LAPW and LO energy parameters
  !> @param[in]   kveclo           : Basis G-vectors of local orbitals
  !> @param[in]   uuilon           : overlap integral between the radial functions of the integral (multiplied by ulo_der) of a
  !>                                 local orbital and the flapw radial function with the same l
  !> @param[in]   duilon           : overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                                 derivative of the flapw radial function with the same l
  !> @param[in]   ulouilopn        : overlap integral between the radial functions of the integral of a local orbital and another
  !>                                 local orbital with the same l.
  !> @param[in]   gdp              : G-vectors of potentials and densities
  !> @param[in]   qpwcG            : Plane-wave expansion coefficients of the FFT of the gradient of the pseudo core-density ( Gauss
  !>                                 curve )
  !> @param[in]   mapKpq2K         : For a given k-point index and q-point index, this array gives the k-point set index of the
  !>                                 result k + q (already mapped back to Brillouin zone).
  !> @param[in]   ne               : Number of eigenvalues per k-point.
  !> @param[in]   eig              : Contains Kohn\--Sham eigenvalues.
  !> @param[in]   gbas             : G-basis vectors
  !> @param[in]   mapGbas          : For various k-points G-basis vectors occur more than once, thus they are only stored once in
  !>                                 juPhon. This pointer array contains the right index for GbasVec array to "unfold" G-basis
  !>                                 vectors again.
  !> @param[in]   z                : Kohn\--Sham eigenvectors.
  !> @param[in]   nv               : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]   nRadFun          : Number of radial functions per orbital quantum number l and atom type.
  !> @param[in]   iloTable         : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is
  !>                                 given.
  !> @param[in]   nobd             : Number of occupied bands per k-point and spin
  !> @param[in]   ilo2p            : mapping array giving the p value for given number of LO and itype
  !> @param[in]   gdp2Ind          : Stores the index of a potential and density G-Vector. The dimensions are the G-vector
  !>                                 components.
  !> @param[in]   gdp2iLim         : Stores the min and maxvals of gdp2Ind
  !> @param[in]   kpq2kPrVec       : Backfolding vector from k + q in 2nd Brillouin zone to k' in 1st Brillouin zone
  !> @param[in]   ylm              : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to
  !>                                 lmax created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   grRho0IR         : Plane-wave coefficients of the gradient of the interstitial unperturbed density
  !> @param[in]   grRho0MT         : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[in]   grVeff0MT_init   : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed external
  !>                                 potential (part of the effective) for initial Sternheimer iteration created in
  !>                                 m_jpGrVeff0::gengrveff0
  !> @param[in]   grVeff0MT_main   : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective
  !>                                 potential for a regular and the final Sternheimer iteration created in m_jpGrVeff0::gengrveff0
  !> @param[in]   dKernMTGPts      : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                                 respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   vxc1IRKern       : Plane-wave interstitial coefficients of the functional derivative of the xc-kernel with respect
  !>                                 to the density
  !> @param[in]   rho1MTCoreDispAt : Spherical-harmonic expansion coefficients of the gradient of the pseudo core-density ( Gauss
  !>                                 curve ) at the displaced atom.
  !> @param[in]   rho0IRpw         : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]   rho0MTsh         : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from
  !>                                 Fleur
  !> @param[in]   gausWts          : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   vEff0IRpwUw      : Plane-wave coefficients of the unperturbed and converged interstitial effective potential read
  !>                                 in from Fleur without Heaviside warping.
  !> @param[out]  rho1IRDS         : Plane-wave interstitial coefficients of the linear variation of the interstitial density
  !> @param[out]  rho1MT           : Spherical harmonic expansion coefficients of the full linear variation of the muffin-tin
  !>                                 density
  !> @param[out]  rho1MTDelta      : Spherical harmonic expansion coefficients of the full linear variation of the muffin-tin
  !>                                 density minus the gradient of the unperturbed density.
  !> @param[out]  vExt1IR_final    : Plane-wave coefficients of the interstitial converged external linear potential variation
  !> @param[out]  vHar1IR_final    : Plane-wave coefficients of the interstitial converged Hartree linear potential variation
  !> @param[out]  vExt1MT          : Spherical harmonic expansion coefficients of the complete converged external muffin-tin linear
  !>                                 potential variation
  !> @param[out]  vExt1MTDelta     : Spherical harmonic expansion coefficients of the difference of the linear variation of the
  !>                                 external potential between a finite q and q=0; only surface contribution no volume term.
  !> @param[out]  vExt1MTq0        : Spherical harmonic expansion coefficients of the full linear variation of the external
  !>                                 potential for q = 0
  !> @param[out]  vHar1MTDelta     : Spherical harmonic expansion coefficients of the difference of the linear variation of the
  !>                                 Hartree potential between a finite q and q=0; only surface contribution no volume term.
  !> @param[out]  vHar1MTq0        : Spherical harmonic expansion coefficients of the full linear variation of the Hartree
  !>                                 potential for q = 0
  !> @param[out]  vXc1MTDelta      : Spherical harmonic expansion coefficients of the of the linear variation of the xc potential
  !>                                 for finite q (if q/=0 otherwise it is set to zero).
  !> @param[out]  vXc1MTq0         : Spherical harmonic expansion coefficients of the full linear variation of the xc potential for
  !>                                 q = 0
  !> @param[out]  vEff1IR_final    : Plane-wave coefficients of the interstitial converged effective linear potential variation
  !> @param[out]  vHar1MT_final    : Spherical harmonic expansion coefficients of the complete converged Hartree muffin-tin linear
  !>                                 potential variation
  !> @param[out]  vEff1MT_final    : Spherical harmonic expansion coefficients of the complete converged effective muffin-tin linear
  !>                                 potential variation
  !> @param[out]  gpqdp            : Number of G-vectors for shifted G-set for q-point with index iqpt
  !>
  !> @details
  !> After the Sternheimer SCC has been initialized in m_jpSternheimer::initSternheimerSCC, for all displaced atoms and for a
  !> certain q, the Sternheimer self-consistency cycle is performed and the continuity(see log) of the resulting quantities checked.
  !> See also Section 7.4.5 (dissertation CRG).
  !> @note
  !> Handling the complex datatype for the zs. Can be removed when latest FLEUR convention is applied for z datatype.
  !>
  !> @todo
  !> write logFileOutput
  !> @todo
  !> add spin
  !> @note
  !> Spin is hard-coded
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine solveSternheimerSCC( fmpi,   atoms, sym, stars, lathar, cell, enpara, uds, input, kpts, qpts, results, usdus, logUnit,&
    & ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, gdp, mapKpq2K, ne, eig, gbas, mapGbas, z, nv, El, nRadFun, iloTable, &
    & nobd, ilo2p, gdp2Ind, gdp2iLim, kpq2kPrVec, qpwcG, iqpt, tdHS0, loosetd, ylm, grRho0IR, grRho0MT, grVeff0MT_init, grVeff0MT_main, &
    & dKernMTGPts, vxc1IRKern, rho1MTCoreDispAt, gausWts, rho1IRDS, rho1MT, vExt1MT, vEff1IR_final, vEff1MT_final, &
    & oneSternhCycle, ngpqdp, gpqdp, vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MTDelta, vExt1MTDelta, vExt1MTq0, &
    & vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, rho0IRpw, rho0MTsh, vEff0IRpwUw, noPts2chCont, vEff1MTnoVol, &
    & vExt1noqIR_final, rho1MTz0, vCoul1IRtempNoVol, vCoul1MTtempNoVol, vEff0MTsh )

    use m_types
    use m_jpDens1stVar, only : calcVarCTcorr, calcRho1IRValDS, calcKdepValRho1MT, multRadSolVzcRho1MT
    use m_juDFT_time, only : TimeStart, Timestop
    use m_jpSternhHF, only : tlmplm4V, uDHu, ctorsh
    use m_jpSternhPulaySurface, only : calcSfVeffFast, IRcoeffVeffUv
    use m_jpIOnMixing, only : UpdNCheckDens, loadDensity, storeZ1nG
    use m_juDFT_stop, only : juDFT_warn, juDFT_error
    USE m_npy

    implicit none

    ! Type parameters
    type(t_mpi),                  intent(in)  :: fmpi
     
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sym),                    intent(in)  :: sym
    type(t_stars),                  intent(in)  :: stars
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_cell),                   intent(in)  :: cell
    type(t_enpara),                 intent(in)  :: enpara
    type(t_usdus),                  intent(in)  :: uds
    type(t_input),                  intent(in)  :: input
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_results),                intent(in)  :: results
    type(t_usdus),                  intent(in)  :: usdus
    type(t_tlmplm),                 intent(in)  :: tdHS0

    ! Scalar parameters
    integer,                        intent(in)  :: logUnit
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: iqpt
    logical,                        intent(in)  :: oneSternhCycle
    integer,                        intent(in)  :: noPts2chCont

    ! Array parameters
    real,                           intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    integer,                        intent(in)  :: kveclo(:,:)
    real,                           intent(in)  :: uuilon(:, :)
    real,                           intent(in)  :: duilon(:, :)
    real,                           intent(in)  :: ulouilopn(:, :, :)
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: qpwcG(:, :)
    integer,                        intent(in)  :: mapKpq2K(:, :)
    integer,                        intent(in)  :: ne(:)
    real,                           intent(in)  :: eig(:,:,:)
    integer,                        intent(in)  :: gbas(:, :)
    integer,                        intent(in)  :: mapGbas(:, :, :)
    complex,                       intent(in)  :: z(:,:,:,:)
    integer,                        intent(in)  :: nv(:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: ilo2p(:, :)
    integer,                        intent(in)  :: gdp2Ind(:, :, :) ! will substitute shifted G-set
    integer,                        intent(in)  :: gdp2iLim(2, 3)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    complex,                        intent(in)  :: ylm(:, :)
    complex,                        intent(in)  :: grRho0IR(:, :)
    complex,                        intent(in)  :: grRho0MT(:, :, :, :)
    complex,                        intent(in)  :: grVeff0MT_init(:, :, :, :)
    complex,                        intent(in)  :: grVeff0MT_main(:, :, :, :)
    real,                           intent(in)  :: dKernMTGPts(:, :, :)
    complex,                        intent(in)  :: vxc1IRKern(:)
    complex,                        intent(in)  :: rho1MTCoreDispAt(:, :, :, :)
    complex,                        intent(in)  :: rho0IRpw(:, :)
    complex,                        intent(in)  :: rho0MTsh(:, :, :, :)
    real,                           intent(in)  :: gausWts(:) ! gaussian weights belonging to gausPts
    complex,                        intent(in)  :: vEff0IRpwUw(:, :)
    complex,           allocatable, intent(out) :: rho1IRDS(:, :, :)
    complex,           allocatable, intent(out) :: rho1MT(:, :, :, :, :)
    complex,           allocatable, intent(out) :: rho1MTDelta(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vExt1IR_final(:, :, :)
    complex,           allocatable, intent(out) :: vExt1noqIR_final(:, :, :)
    complex,           allocatable, intent(out) :: vHar1IR_final(:, :, :)
    complex,           allocatable, intent(out) :: vExt1MT(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vEff1MTnoVol(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vExt1MTDelta(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vExt1MTq0(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vHar1MTDelta(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vHar1MTq0(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vXc1MTDelta(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vXc1MTq0(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vEff1IR_final(:, :, :)
    complex,           allocatable, intent(out) :: vHar1MT_final(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vEff1MT_final(:, :, :, :, :)
    integer,           allocatable, intent(out) :: gpqdp(:, :) ! Deprecated, should be removed when G-set shift is removed
    complex,           allocatable, intent(out) :: rho1MTz0(:, :, :, :)
    complex,           allocatable, intent(out) :: vCoul1IRtempNoVol(:, :)
    complex,           allocatable, intent(out) :: vCoul1MTtempNoVol(:, :, :, :)
    complex,           optional,    intent(in)  :: vEff0MTsh(:, :, :, :)
    complex, intent(in) :: loosetd(:, :, :, :)

    ! Scalar variables
    logical                                     :: stern1stIt ! true, if 1st Sternheimer cycle iteration
    logical                                     :: stern2ndIt ! true, if 1st Sternheimer cycle iteration
    logical                                     :: sternRegIt ! true, if not 1st or final Sternheimer cycle iteration
    logical                                     :: sternFinIt
    integer                                     :: iatom
    integer                                     :: iDtype
    integer                                     :: ieqat
    integer                                     :: ispin
    integer                                     :: ikpt
    integer                                     :: idir
    integer                                     :: iDatom
    integer                                     :: itype
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: imesh
    integer                                     :: iG
    integer                                     :: iDeqat
    logical                                     :: densFfile
    logical                                     :: densFileOK
    integer                                     :: counter
    integer                                     :: coScale
    integer                                     :: ngpqdp ! Deprecated
    integer                                     :: ngpqdp2km ! Deprecated
    integer                                     :: iterations
    integer                                     :: maxlmp

    ! Type variables
    type(t_tlmplm)                              :: tdVx
    type(t_tlmplm)                              :: tdVy
    type(t_tlmplm)                              :: tdVz
    type(t_tlmplm)                              :: tdVx2
    type(t_tlmplm)                              :: tdVy2
    type(t_tlmplm)                              :: tdVz2

    ! Array variables
    complex,           allocatable              :: sumVMTs(:, :, :, :)
    complex,           allocatable              :: sumVMTs2(:, :, :, :)
    complex,           allocatable              :: rho1IRDSplus(:, :)
    complex,           allocatable              :: rho1MTplus(:, :, :, :)
    complex,           allocatable              :: rho1MTplus2(:, :, :, :)
    complex,           allocatable              :: z1nG(:, :, :)
    real,              allocatable              :: recEdiffME(:, :)
    real,              allocatable              :: cutContr(:, :) ! Deprecated !todo cutContr should be possibly an integer
    integer,           allocatable              :: nlo_atom(:)
    complex,           allocatable              :: rho1IRctC(:, :, :)
    complex,           allocatable              :: rho1MTctC(:, :, :, :, :)
    complex,           allocatable              :: uu(:,:,:)
    complex,           allocatable              :: du(:,:,:)
    complex,           allocatable              :: dd(:,:,:)
    complex,           allocatable              :: ud(:,:,:)
    complex,           allocatable              :: aclo(:,:,:)
    complex,           allocatable              :: bclo(:,:,:)
    complex,           allocatable              :: cclo(:,:,:,:)
    complex,           allocatable              :: uunmt(:,:,:,:,:)
    complex,           allocatable              :: udnmt(:,:,:,:,:)
    complex,           allocatable              :: dunmt(:,:,:,:,:)
    complex,           allocatable              :: ddnmt(:,:,:,:,:)
    complex,           allocatable              :: acnmt(:,:,:,:,:)
    complex,           allocatable              :: bcnmt(:,:,:,:,:)
    complex,           allocatable              :: ccnmt(:,:,:,:,:)
    complex,           allocatable              :: uu2(:,:,:)
    complex,           allocatable              :: du2(:,:,:)
    complex,           allocatable              :: dd2(:,:,:)
    complex,           allocatable              :: ud2(:,:,:)
    complex,           allocatable              :: aclo2(:,:,:)
    complex,           allocatable              :: bclo2(:,:,:)
    complex,           allocatable              :: cclo2(:,:,:,:)
    complex,           allocatable              :: uunmt2(:,:,:,:,:)
    complex,           allocatable              :: udnmt2(:,:,:,:,:)
    complex,           allocatable              :: dunmt2(:,:,:,:,:)
    complex,           allocatable              :: ddnmt2(:,:,:,:,:)
    complex,           allocatable              :: acnmt2(:,:,:,:,:)
    complex,           allocatable              :: bcnmt2(:,:,:,:,:)
    complex,           allocatable              :: ccnmt2(:,:,:,:,:)
    complex,           allocatable              :: veffUvIR(:, :)
    complex,           allocatable              :: lastDistance(:, :, :)
    complex,           allocatable              :: surfIntVFast(:, :, :)
    complex,           allocatable              :: vEff1IR(:, :)
    integer,           allocatable              :: gpqdp2Ind(:, :, :) ! Deprecated
    integer                                     :: gpqdp2iLim(2, 3) ! Deprecated
    INTEGER                                     :: killcont(9)
    logical                                     :: converged(3)
    integer,           allocatable              :: gdpIndex(:, :, :) ! Will substitute shifted G-set
    character(len=24)                           :: filename
    character(len=2)                            :: filename1

    complex                                     :: haa((atoms%lmaxd+1)**2,2,0:atoms%lmaxd,2,0:atoms%lmaxd), &
                                                  dhaa((atoms%lmaxd+1)**2,2,0:atoms%lmaxd,2,0:atoms%lmaxd)
    real                                        :: vEff0MTrsh(atoms%jmtd,(atoms%lmaxd+1)**2,atoms%nat,input%jspins)
    complex, allocatable :: mat_elH(:,:,:)
    complex, allocatable :: mat_elS(:,:)

    complex, allocatable :: loosetdx1(:, :, :, :), loosetdx2(:, :, :, :)
    complex, allocatable :: loosetdy1(:, :, :, :), loosetdy2(:, :, :, :)
    complex, allocatable :: loosetdz1(:, :, :, :), loosetdz2(:, :, :, :)

    ! Comments after activated quantities indicate the problems/differences to old juPhon.
    !CALL save_npy('rbas1.npy',rbas1)
    !CALL save_npy('rbas2.npy',rbas2)
    !CALL save_npy('El.npy',El) ! Size differs (certainly LO related). Max seems nicer, but the values match.
    !CALL save_npy('kveclo.npy',kveclo)
    !CALL save_npy('uuilon.npy',uuilon)
    !CALL save_npy('duilon.npy',duilon)
    !CALL save_npy('ulouilopn.npy',ulouilopn)
    !CALL save_npy('gdp.npy',gdp)
    !CALL save_npy('qpwcG.npy',qpwcG)
    !CALL save_npy('mapKpq2K.npy',mapKpq2K)
    !CALL save_npy('ne.npy',ne)
    !CALL save_npy('eig.npy',eig) ! Slightly different values; some high values uninitialized instead of fixed.
    !CALL save_npy('gbas.npy',gbas)
    !CALL save_npy('mapGbas.npy',mapGbas)
    !CALL save_npy('z0.npy',z) ! Not all eigenvalues sufficiently filled.
    !CALL save_npy('nv.npy',nv)
    !CALL save_npy('nRadFun.npy',nRadFun)
    !CALL save_npy('iloTable.npy',iloTable)
    !CALL save_npy('nobd.npy',nobd)
    !CALL save_npy('ilo2p.npy',ilo2p)
    !CALL save_npy('gdp2Ind.npy',gdp2Ind) ! will substitute shifted G-set
    !CALL save_npy('gdp2iLim.npy',gdp2iLim)
    !CALL save_npy('kpq2kPrVec.npy',kpq2kPrVec)
    !CALL save_npy('ylm.npy',ylm)
    !CALL save_npy('grRho0IR.npy',grRho0IR)
    !CALL save_npy('grRho0MT.npy',grRho0MT) ! Fixed.
    !CALL save_npy('grVeff0MT_init.npy',grVeff0MT_init)
    !CALL save_npy('grVeff0MT_main.npy',grVeff0MT_main)
    !CALL save_npy('dKernMTGPts.npy',dKernMTGPts) ! Fixed.
    !CALL save_npy('vxc1IRKern.npy',vxc1IRKern)
    !CALL save_npy('rho1MTCoreDispAt.npy',rho1MTCoreDispAt)
    !CALL save_npy('rho0IRpw.npy',rho0IRpw)
    !CALL save_npy('rho0MTsh.npy',rho0MTsh)
    !CALL save_npy('gausWts.npy',gausWts) ! gaussian weights belonging to gausPts
    !CALL save_npy('vEff0IRpwUw.npy', vEff0IRpwUw)
    !CALL save_npy('tuu.npy', loosetd(:, :, :, 1))
    !CALL save_npy('tud.npy', loosetd(:, :, :, 2))
    !CALL save_npy('tdu.npy', loosetd(:, :, :, 3))
    !CALL save_npy('tdd.npy', loosetd(:, :, :, 4))
    !stop
    ! Generate potential G-set which is shifted according to the q-vector so fulfills ||q + G|| < Gmax
    ! NOTE: This routine is deprecated, as it was decided to not shift the set of the G-vectors, therefore it is deactivated at the
    !       moment by just generating a G-set with no shift that is equivalent to the standard G-set for the density and the
    !       potential given by Fleur. After some convergance tests, this should be removed from all routines.
    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpts%bk(:, 1), gpqdp, gpqdp2Ind, gpqdp2iLim )
    ! Test that there is the same G-set as there would be for q = 0
    if ( ngpqdp /= ngdp ) then
      call juDFT_warn('ngpqdp /= ngdp', calledby='solveSternheimerSCC', hint='Do not shift the G-set!')
    else
      do iG = 1, ngdp
        if ( any( gdp(:, iG) /= gpqdp(:, iG) ) ) then
          call juDFT_warn('gpqdp /= gdp', calledby='solveSternheimerSCC', hint='Do not shift the G-set!')
        end if
      end do ! iG
    end if

    !todo check the dimensions of the arrays
    allocate( rho1IRDS(ngpqdp, 3, atoms%nat) )
    allocate( rho1IRDSplus(ngpqdp, 3) )
    allocate( rho1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat))
    allocate( rho1MTz0(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat))
    allocate( rho1MTDelta(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat))
    allocate( rho1MTplus(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3))
    allocate( rho1MTplus2(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3))
    allocate( sumVMTs(atoms%jmtd, (atoms%lmaxd + 1 )**2, 3, atoms%nat) )
    allocate( sumVMTs2(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, atoms%nat) )
    allocate( vExt1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vEff1MTnoVol(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vExt1MTDelta(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vExt1MTq0(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vHar1MTDelta(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vHar1MTq0(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vXc1MTDelta(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vXc1MTq0(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vEff1IR_final( ngpqdp, 3, atoms%nat ) )
    allocate( vHar1IR_final( ngpqdp, 3, atoms%nat ) )
    allocate( vExt1IR_final( ngpqdp, 3, atoms%nat ) )
    allocate( vExt1noqIR_final( ngdp, 3, atoms%nat ) )
    allocate( vEff1MT_final(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    allocate( vHar1MT_final(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat ))
    ! Note: The arrays for uu, du, dd and so on have spin index, which is now removed and should be added later
    allocate( uu( 0 : atoms%lmaxd, atoms%nat, 3 ), du( 0 : atoms%lmaxd, atoms%nat, 3 ), dd( 0 : atoms%lmaxd, atoms%nat, 3 ), &
                                                                                              & ud( 0 : atoms%lmaxd, atoms%nat, 3) )
    allocate( uunmt( 0 : atoms%lmaxd, 0: atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & ddnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & udnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & dunmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3) )
    allocate( acnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &bcnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &ccnmt(atoms%nlod, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
    allocate( aclo(atoms%nlod, atoms%nat, 3), bclo(atoms%nlod, atoms%nat, 3), cclo(atoms%nlod, atoms%nlod, atoms%nat, 3) )
    allocate( uu2( 0 : atoms%lmaxd, atoms%nat, 3 ), du2( 0 : atoms%lmaxd, atoms%nat, 3 ), dd2( 0 : atoms%lmaxd, atoms%nat, 3 ), &
                                                                                              & ud2( 0 : atoms%lmaxd, atoms%nat, 3) )
    allocate( uunmt2( 0 : atoms%lmaxd, 0: atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & ddnmt2( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & udnmt2( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3), &
      & dunmt2( 0 : atoms%lmaxd, 0 : atoms%lmaxd, (atoms%lmaxd + 1)**2, atoms%nat, 3) )
    allocate( acnmt2( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &bcnmt2( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &ccnmt2(atoms%nlod, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
    allocate( aclo2(atoms%nlod, atoms%nat, 3), bclo2(atoms%nlod, atoms%nat, 3), cclo2(atoms%nlod, atoms%nlod, atoms%nat, 3) )
    allocate( lastDistance(3, atoms%nat, qpts%nkpt) )
    allocate( rho1IRctC(ngpqdp, atoms%nat, 3), rho1MTctC(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat) )
    allocate( veffUvIR(3, (atoms%lmaxd + 1)**2) )
    !todo here a nobdMax could be introduced
    allocate( surfIntVFast(input%neig, maxval(nobd), 3) )
    allocate( z1nG(SIZE(z(:,1,1,1)), input%neig, 3) ) ! todo segfault in abcof if second argument not like this, this can be optimized
    allocate( recEdiffME(input%neig, maxval(nobd)) )
    allocate( cutContr(input%neig, maxval(nobd)) )


    rho1IRDS(:, :, :)            = cmplx(0., 0.)
    rho1IRDSplus(:, :)           = cmplx(0., 0.)
    rho1MT(:, :, :, :, :)        = cmplx(0., 0.)
    rho1MTz0(:, :, :, :)        = cmplx(0., 0.)
    rho1MTDelta(:, :, :, :, :)   = cmplx(0., 0.)
    rho1MTplus(:, :, :, :)       = cmplx(0., 0.)
    rho1MTplus2(:, :, :, :)       = cmplx(0., 0.)
    sumVMTs(:, :, :, :)          = cmplx(0., 0.)
    sumVMTs2(:, :, :, :)         = cmplx(0., 0.)
    vExt1MT(:, :, :, :, :)       = cmplx(0., 0.)
    vEff1MTnoVol(:, :, :, :, :)       = cmplx(0., 0.)
    vExt1MTDelta(:, :, :, :, :)  = cmplx(0., 0.)
    vExt1MTq0(:, :, :, :, :)     = cmplx(0., 0.)
    vHar1MTDelta(:, :, :, :, :)  = cmplx(0., 0.)
    vHar1MTq0(:, :, :, :, :)     = cmplx(0., 0.)
    vXc1MTDelta(:, :, :, :, :)   = cmplx(0., 0.)
    vXc1MTq0(:, :, :, :, :)      = cmplx(0., 0.)
    vEff1IR_final(:, :, :)       = cmplx(0., 0.)
    vHar1IR_final(:, :, :)       = cmplx(0., 0.)
    vExt1IR_final(:, :, :)       = cmplx(0., 0.)
    vExt1noqIR_final(:, :, :)       = cmplx(0., 0.)
    vEff1MT_final(:, :, :, :, :) = cmplx(0., 0.)
    vHar1MT_final(:, :, :, :, :) = cmplx(0., 0.)
    uu(:, :, :)                  = cmplx(0., 0.)
    du(:, :, :)                  = cmplx(0., 0.)
    dd(:, :, :)                  = cmplx(0., 0.)
    ud(:, :, :)                  = cmplx(0., 0.)
    uunmt(:, :, :, :, :)         = cmplx(0., 0.)
    ddnmt(:, :, :, :, :)         = cmplx(0., 0.)
    udnmt(:, :, :, :, :)         = cmplx(0., 0.)
    dunmt(:, :, :, :, :)         = cmplx(0., 0.)
    acnmt(:, :, :, :, :)         = cmplx(0., 0.)
    bcnmt(:, :, :, :, :)         = cmplx(0., 0.)
    ccnmt(:, :, :, :, :)         = cmplx(0., 0.)
    aclo(:, :, :)                = cmplx(0., 0.)
    bclo(:, :, :)                = cmplx(0., 0.)
    cclo(:, :, :, :)             = cmplx(0., 0.)
    uu2(:, :, :)                  = cmplx(0., 0.)
    du2(:, :, :)                  = cmplx(0., 0.)
    dd2(:, :, :)                  = cmplx(0., 0.)
    ud2(:, :, :)                  = cmplx(0., 0.)
    uunmt2(:, :, :, :, :)         = cmplx(0., 0.)
    ddnmt2(:, :, :, :, :)         = cmplx(0., 0.)
    udnmt2(:, :, :, :, :)         = cmplx(0., 0.)
    dunmt2(:, :, :, :, :)         = cmplx(0., 0.)
    acnmt2(:, :, :, :, :)         = cmplx(0., 0.)
    bcnmt2(:, :, :, :, :)         = cmplx(0., 0.)
    ccnmt2(:, :, :, :, :)         = cmplx(0., 0.)
    aclo2(:, :, :)                = cmplx(0., 0.)
    bclo2(:, :, :)                = cmplx(0., 0.)
    cclo2(:, :, :, :)             = cmplx(0., 0.)
    lastDistance(:, :, :)        = cmplx(0., 0.)
    rho1IRctC(:, :, :)           = cmplx(0., 0.)
    rho1MTctC(:, :, :, :, :)     = cmplx(0., 0.)
    z1nG(:, :, :)                = cmplx(0., 0.)
    recEdiffME(:, :)             = 0.
    cutContr(:, :)               = 0.

    ! Scaling cutoff parameter for surface integral in IR that is Rayleigh expanded, Aaron claims for the forces, it has to be
    ! 2 lmax?
    coScale = 1.

    ! This array stores which direction of displacement is already converged
    converged(:) = .false.

    ! Hard-code spin to index 1
    ispin = 1

    ! Counter for enumerating the file names for debugging during the iterations
    counter = 0

    ! If atom alpha is displaced also its coretail changes (in the IR and the MT of the remaining atoms). Thus, we have contribution
    ! contributions to the first variation of the density in the IR and the MTs. Since in FLEUR we also add this correction of the
    ! coretail to the atom inducing the coretail (which is wrong) we do it here also and subtract it later for the displaced atom
    ! alpha (as it is done in FLEUR).
    call TimeStart( 'Sh: calcCtCt' )
    call calcVarCTcorr( atoms, cell, ngpqdp, gpqdp, qpts%bk(:, iqpt), qpwcG, rho1IRctC, rho1MTctC )
    call Timestop( 'Sh: calcCtCt' )

    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        if (oneSternhCycle) then
          ! When a density is converged, there is the possibility to only perform the final run of Sternheimer filling all resulting
          ! arrays for post-processing. Then, no mixing is performed, so the distance file does not have to be opened. The flags for
          ! controlling the Sternheimer routiens are set according to that.
          write(*, '(a)') 'Final Sternheimer iteration mode activated!'
          stern1stIt = .false.
          stern2ndIt = .false.
          sternRegIt = .false.
          sternFinIt = .true.
        else
          ! If more than one Sternheimer cycle is set, it is assumed that the self-consistency cycle should start from the very
          ! beginning, i.e. create a starting density with only the external potential variation and so on and so are the flags for
          ! controlling the Sternheimer routines are set. As the density variation is to be mixed a distance file is opened tracking
          ! the evolution of the distance.
          !todo we have to change 789 to a name standing for the unitnumber
          write(filename, '(a10,i1,a4,i4)') 'distance', iDatom, 'qInd', iqpt
          open(789, file=filename, status='replace', action='write', form='formatted')
          stern1stIt = .true.
          stern2ndIt = .false.
          sternRegIt = .false.
          sternFinIt = .false.
        end if

        ! The arrays for storing the q+ linear variation of the density in the interstitial and the muffin-tin that result from the
        ! Sternheimer self-consistency cycle (without the gradient) are to be reset for every displaced atom
        rho1IRDSplus(:, :) = cmplx(0.0, 0.0)
        rho1MTplus(:, :, :, :) = cmplx(0.0, 0.0)
        rho1MTplus2(:, :, :, :) = cmplx(0.0, 0.0)

        ! If one wants to start from a certain iteration or a starting density read in from file, a (empty) file with the name
        ! start_cdn1 has to exist, it is then filled into tho rho1IRDSplus and rho1MTplus so that these arrays can be converged
        ! further
        ! Reset for new atom displaced
        inquire(file='start_cdn1', exist=densFfile)
        if (densFfile) then
          do idir = 1, 3
            write(filename, '(a10,i1,a4,i1,a4,i4)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
            inquire(file=filename, exist=densFileOK)
            if (densFileOK) then
              write(*, '(a,1x,i1,1x,a,1x,i4,1x,a,a)') 'Variation of density for displaced atom', iDatom, 'and q-point', iqpt, &
                                                                                                         & 'is read from ', filename
              call loadDensity(atoms, ngpqdp, filename, rho1IRDSplus(:, idir), rho1MTplus(:, :, :, idir))
            else
              write(*, '(a,1x,i1,1x,a,1x,i2,1x,a,a)') 'File for variation of density for displaced atom', iDatom, 'and q-point', &
                                                                                  & iqpt, 'with the name', filename, 'is not found.'
              CALL juDFT_error('Old stopping from t_juPhon',calledby ="solveSternheimerSCC")
            end if
          end do
          write(*, '(a)') 'Regular Sternheimer iteration mode activated!'
          !todo Check whether the density is stored only after the second run with the full effective potential activated
          ! If an initial iteration with only the external potential varaition activated was assumed, then after having read in the
          ! density variation from file a regular iteration should be performed and not any more an initial iteration.
          ! If, however, only the final iteration should be performed then stern1stIt was not activated so none of the controlling
          ! flags of the Sternheimer iteration is changed.
          if (stern1stIt) then
            stern1stIt = .false.
            sternRegIt = .true.
          end if ! stern1stIt
        end if ! densFile

        ! Reset counter for number of iterations (is to be used for log output later)
        iterations = 0

        ! Precalculation of potential for surface integral that is not k-dependent
        veffUvIR(:, :) = cmplx(0., 0.)
        !todo shift to sternheimerPUlay
        call IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, veffUvIR, vEff0IRpwUw )
        !CALL save_npy('veffUvIR.npy', veffUvIR) !Slightly different.

        ! The lmp index is comparable to the lm index but for every lm combination it includes also an index that runs over the
        ! matching coefficients for belonging to the u, the udot and the u_LO. Determine the max lmp for all atom types.
        maxlmp = maxval( (/ (sum( (/ ((2 * oqn_l + 1) * nRadFun(oqn_l,iDtype), oqn_l = 0, atoms%lmax(iDtype)) /) ), iDtype = 1,&
                                                                                                                  &atoms%ntype) /) )

        do ! self-consistency loop
          iterations = iterations + 1

          ! Knowing the index of the currently displaced atom and the phonon q-vector, we can now calculate the first variation of
          ! the effective potential and add it to the gradient of the unperturbed potential we have calculated outside the loop.
          ! This sum of the potentials is assigned to the tlmplm routine for every direction of displacement so that we can reuse
          ! the resulting integrals for the Sternheimer equation matrix element.

          call UpdPots( atoms, input, stars, cell, iDatom, iDtype, stern1stIt, ngdp, grRho0IR, grRho0MT, qpts%bk(1:3, iqpt), &
            & rho1IRDSplus, rho1MTplus, gdp, grVeff0MT_init, grVeff0MT_main, sumVMTs, ylm, gausWts, dKernMTGPts, vxc1IRKern, iqpt, &
            & sumVMTs2, sternFinIt, vExt1MT, vEff1MT_final, vExt1IR_final, vHar1IR_final, vHar1MT_final, vExt1MTDelta, vExt1MTq0, &
            & vHar1MTDelta, vHar1MTq0, vXc1MTDelta, vXc1MTq0, rho0IRpw, rho0MTsh, vEff1IR_final, ngpqdp, gpqdp, vEff1IR, &
            & noPts2chCont, logUnit, vEff1MTnoVol, vExt1noqIR_final, vCoul1IRtempNoVol, vCoul1MTtempNoVol )

            !CALL save_npy('sumVMTs.npy',sumVMTs)
            !CALL save_npy('sumVMTs2.npy',sumVMTs2)
            !CALL save_npy('vExt1IR_final.npy',vExt1IR_final)
            !CALL save_npy('vHar1IR_final.npy',vHar1IR_final)
            !CALL save_npy('vExt1MT.npy',vExt1MT)
            !CALL save_npy('vEff1MTnoVol.npy',vEff1MTnoVol)
            !CALL save_npy('vExt1MTDelta.npy',vExt1MTDelta)
            !CALL save_npy('vExt1MTq0.npy',vExt1MTq0)
            !CALL save_npy('vHar1MTDelta.npy',vHar1MTDelta)
            !CALL save_npy('vHar1MTq0.npy',vHar1MTq0)
            !CALL save_npy('vXc1MTDelta.npy',vXc1MTDelta)
            !CALL save_npy('vXc1MTq0.npy',vXc1MTq0)
            !CALL save_npy('vEff1IR_final.npy',vEff1IR_final)
            !CALL save_npy('vExt1noqIR_final.npy',vExt1noqIR_final)
            !CALL save_npy('vEff1MT_final.npy',vEff1MT_final)
            !CALL save_npy('vHar1MT_final.npy',vHar1MT_final)
            !CALL save_npy('vCoul1IRtempNoVol.npy',vCoul1IRtempNoVol)
            !CALL save_npy('vCoul1MTtempNoVol.npy',vCoul1MTtempNoVol)
            !CALL save_npy('vEff1IRsh.npy',vEff1IR)
            !stop

          ! resetting rho1IRDSplus for next iteration, because the current rho1IRDSplus is contained in the linear potential
          ! variations
          rho1IRDSplus(:, :) = cmplx(0.0, 0.0)
          ! Calculate the part of the MT potential matrix element which is not k-dependet. As the potential a complex part (Veff1)
          ! and we want, in a first step, obey to the Fleur routine structure of tlmplm and hssr_wu, we need a seperate treatment of
          ! the upper and lower half of this tlmplm matrices.
          ! We call the routine at the moment for every displacement direction seperately.
          !if (.false.) then
          !  open(111,file='000_dhmlrad',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
          !end if
          IF (ALLOCATED(loosetdx1)) THEN
              DEALLOCATE(loosetdx1, loosetdy1, loosetdz1, loosetdx2, loosetdy2, loosetdz2)
          END IF
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVx, loosetdx1, 1, 1, 1, sumVMTs, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )
          !if (.false.) then
          !  close(111)
          !end if
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVy, loosetdy1, 1, 1, 2, sumVMTs, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVz, loosetdz1, 1, 1, 3, sumVMTs, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )
          !if (.false.) then
          !  open(111,file='000_dhmlrad',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
          !end if
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVx2, loosetdx2, 1, 1, 1, sumVMTs2, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )
          !if (.false.) then
          !  close(111)
          !end if
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVy2, loosetdy2, 1, 1, 2, sumVMTs2, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )
          call tlmplm4V( atoms, lathar, enpara, uds, input, tdVz2, loosetdz2, 1, 1, 3, sumVMTs2, rbas1, rbas2, uuilon, duilon, &
            & ulouilopn, ilo2p, nlo_atom )

          if (.false..and..false.) then
            haa(:,:,:,:,:)  = cmplx(0.0,0.0)
            dhaa(:,:,:,:,:) = cmplx(0.0,0.0)

            call ctorsh(atoms,vEff0MTsh,vEff0MTrsh)
            call udHu(atoms,El,rbas1,rbas2,vEff0MTrsh,sumVMTs,haa,dhaa)
          end if

          uu(:, :, :)          = cmplx(0., 0.)
          du(:, :, :)          = cmplx(0., 0.)
          ud(:, :, :)          = cmplx(0., 0.)
          dd(:, :, :)          = cmplx(0., 0.)
          aclo(:, :, :)        = cmplx(0., 0.)
          bclo(:, :, :)        = cmplx(0., 0.)
          cclo(:, :, :, :)     = cmplx(0., 0.)
          uunmt(:, :, :, :, :) = cmplx(0., 0.)
          ddnmt(:, :, :, :, :) = cmplx(0., 0.)
          udnmt(:, :, :, :, :) = cmplx(0., 0.)
          dunmt(:, :, :, :, :) = cmplx(0., 0.)
          acnmt(:, :, :, :, :) = cmplx(0., 0.)
          bcnmt(:, :, :, :, :) = cmplx(0., 0.)
          ccnmt(:, :, :, :, :) = cmplx(0., 0.)

          uu2(:, :, :)          = cmplx(0., 0.)
          du2(:, :, :)          = cmplx(0., 0.)
          ud2(:, :, :)          = cmplx(0., 0.)
          dd2(:, :, :)          = cmplx(0., 0.)
          aclo2(:, :, :)        = cmplx(0., 0.)
          bclo2(:, :, :)        = cmplx(0., 0.)
          cclo2(:, :, :, :)     = cmplx(0., 0.)
          uunmt2(:, :, :, :, :) = cmplx(0., 0.)
          ddnmt2(:, :, :, :, :) = cmplx(0., 0.)
          udnmt2(:, :, :, :, :) = cmplx(0., 0.)
          dunmt2(:, :, :, :, :) = cmplx(0., 0.)
          acnmt2(:, :, :, :, :) = cmplx(0., 0.)
          bcnmt2(:, :, :, :, :) = cmplx(0., 0.)
          ccnmt2(:, :, :, :, :) = cmplx(0., 0.)

          z1nG(:, :, :) = cmplx(0., 0.)
          if (.false.) then
            open(109,file='000_coeffs',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
            close(109)
          end if
          do ikpt = 1 , kpts%nkpt ! additional k-points still have to be added

            !if (kpts%nkpt.eq.216) then
            !  !if (ikpt.gt.6) cycle
            !  if (ikpt.ne.165) cycle
            !else if (kpts%nkpt.eq.1728) then
            !  !if (ikpt.gt.11.or.mod(ikpt,2).ne.1) cycle
            !  if (ikpt.ne.1229) cycle
            !else
            !  cycle
            !end if

            !if (kpts%nkpt.eq.1728) then
            !  if (.not.((mod(kpts%bk(1,ikpt)*12.0,2.0).ne.0.0).or.(mod(kpts%bk(2,ikpt)*12.0,2.0).ne.0.0).or.(mod(kpts%bk(3,ikpt)*12.0,2.0).ne.0.0))) cycle
            !end if

            write (*, '("k-point index loop : ",i3," / "i3," = ",f6.2," %")') ikpt, kpts%nkpt, float(ikpt) / float(kpts%nkpt) * 100

            ! The energy differences are calculated for every atom which is redundant. But this is no bottleneck so keeping in
            ! memor is worse.
            recEdiffME(:, :) = 0.
            cutContr(:, :) = 0.

            call GenRecEdifME( ne(mapKpq2K(ikpt, iqpt)), nobd(ikpt, 1), eig(:, mapKpq2K(ikpt, iqpt), 1), eig(:, ikpt, 1), &
                                                                                                            & recEdiffME, cutContr )

            z1nG = cmplx(0.0, 0.0)
            ! Due to performance reasons we do not make a loop over idir (array of types, type of arrays)

            if (.false.) then
              allocate(mat_elH(nv(1,mapKpq2K(ikpt, iqpt)),nv(1,ikpt),3))
              allocate(mat_elS(nv(1,mapKpq2K(ikpt, iqpt)),nv(1,ikpt)))
              mat_elH(:, :, :) = cmplx(0., 0.)
              mat_elS(:, :) = cmplx(0., 0.)
            end if

            ! Calculate effective potential part of IR surface integral within Sternheimer equation
            surfIntVFast(:, :, :) = cmplx(0., 0.)
            if (.false.) then
              call calcSfVeffFast( atoms, input, kpts, qpts, cell, ikpt, mapKpq2K(ikpt, iqpt), iqpt, kpq2kPrVec, &
                & coScale, gbas, veffUvIR, nv, ne, nobd, &
                & mapGbas, z, iDtype, iDatom, surfIntVFast, mat_elH )
            else
               call calcSfVeffFast( atoms, input, kpts, qpts, cell, ikpt, mapKpq2K(ikpt, iqpt), iqpt, kpq2kPrVec, &
                & coScale, gbas, veffUvIR, nv, ne, nobd, &
                & mapGbas, z, iDtype, iDatom, surfIntVFast)
            end if

            !CALL save_npy("H0.npy",mat_elH)

            !CALL save_npy('surfIntVFast.npy',surfIntVFast)
            !stop

            killcont = [1,1,1,1,1,0,1,1,1]
            IF (.NOT.stern1stIt) killcont = [1,1,1,1,1,0,1,1,1]
            ! todo Due to performance reasons we do not make a loop over idir (array of types, type of arrays) is that good?
            if (.false.) then
            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVx, loosetdx1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 1, ngdp, nobd, z1nG(:, :, 1), nlo_atom, recEdiffME, kpq2kPrVec, tdVx2, loosetdx2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont, haa, dhaa, rbas1, mat_elH(:, :, 1), mat_elS)
              mat_elS(:, :) = cmplx(0., 0.)
            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVy, loosetdy1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 2, ngdp, nobd, z1nG(:, :, 2), nlo_atom, recEdiffME, kpq2kPrVec, tdVy2, loosetdy2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont, haa, dhaa, rbas1, mat_elH(:, :, 2), mat_elS)
              mat_elS(:, :) = cmplx(0., 0.)
            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVz, loosetdz1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 3, ngdp, nobd, z1nG(:, :, 3), nlo_atom, recEdiffME, kpq2kPrVec, tdVz2, loosetdz2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont, haa, dhaa, rbas1, mat_elH(:, :, 3), mat_elS)
            else
            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVx, loosetdx1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 1, ngdp, nobd, z1nG(:, :, 1), nlo_atom, recEdiffME, kpq2kPrVec, tdVx2, loosetdx2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont)

            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVy, loosetdy1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 2, ngdp, nobd, z1nG(:, :, 2), nlo_atom, recEdiffME, kpq2kPrVec, tdVy2, loosetdy2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont)

            call solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, uds, tdHS0, loosetd, tdVz, loosetdz1, stars, gdp, ne, nv, vEff1IR, &
              & eig(:, ikpt, ispin), eig(:, mapKpq2K(ikpt, iqpt), ispin), El, nRadFun, iloTable, mapGbas, gbas, &
              & z(:, :, ikpt, ispin), z(:, :, mapKpq2K(ikpt,iqpt), ispin), kveclo(:, :), iDtype, iDatom, ikpt, &
              & mapKpq2K(ikpt, iqpt), iqpt, 3, ngdp, nobd, z1nG(:, :, 3), nlo_atom, recEdiffME, kpq2kPrVec, tdVz2, loosetdz2, cutContr, &
              & surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont)
            end if

            !DEALLOCATE(loosetdx1, loosetdy1, loosetdz1, loosetdx2, loosetdy2, loosetdz2)

            if (.false.) then
              deallocate(mat_elH,mat_elS)
            end if

            if ( sternFinIt ) then
              call storeZ1nG( atoms, ikpt, iqpt, mapKpq2K, iDatom, nobd, nv, z1nG )
            end if

            ! Calculate k-dependent contributions of rho1 and perform a sum over k (outer loop within this very routine)
            do idir = 1, 3
              call calcRho1IRValDS( cell, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, mapKpq2K(ikpt,iqpt), idir, &
                & gbas(:, :), z(:, :, ikpt, ispin), z1nG, rho1IRDSplus(:, :), gpqdp2Ind, mapGbas, gpqdp2iLim, kpq2kPrVec)
            end do ! idir

            ! Calculate the k-dependent parts of the density. By looping over the k-points, we perform the sum over the k-points.
            if (sternFinIt) then
              call calcKdepValRho1MT( fmpi,   atoms, input, sym, cell, kpts, usdus, results, ikpt, mapKpq2K(ikpt, iqpt), iDatom, nv, &
                & mapGbas, gBas, nobd(:, 1), uu, du, ud, dd, aclo, bclo, cclo, uunmt, udnmt, dunmt, ddnmt, acnmt, bcnmt, ccnmt, &
                & z(:, :, ikpt, ispin), z1nG, kveclo, rbas1, rbas2, ilo2p, uu2, du2, ud2, dd2, aclo2, bclo2, cclo2, uunmt2, udnmt2, dunmt2, ddnmt2, acnmt2, bcnmt2, ccnmt2 )
            else
              call calcKdepValRho1MT( fmpi,   atoms, input, sym, cell, kpts, usdus, results, ikpt, mapKpq2K(ikpt, iqpt), iDatom, nv, &
                & mapGbas, gBas, nobd(:, 1), uu, du, ud, dd, aclo, bclo, cclo, uunmt, udnmt, dunmt, ddnmt, acnmt, bcnmt, ccnmt, &
                & z(:, :, ikpt, ispin), z1nG, kveclo, rbas1, rbas2, ilo2p )
            end if
          end do ! ikpt


          ! Reset rho1MTplus because it already has been used to create the potentials.
          rho1MTplus(:, :, :, :) = cmplx(0.0, 0.0)
          ! calculate the complete MT valence density
          call multRadSolVzcRho1MT( atoms, aclo, bclo, cclo, acnmt, bcnmt, ccnmt, rbas1, rbas2, uu, du, ud, dd, uunmt, udnmt, &
            & dunmt, ddnmt, ilo2p, rho1MTplus(:, :, :, :) )
          if (sternFinIt) then
            call multRadSolVzcRho1MT( atoms, aclo2, bclo2, cclo2, acnmt2, bcnmt2, ccnmt2, rbas1, rbas2, uu2, du2, ud2, dd2, uunmt2, udnmt2, &
              & dunmt2, ddnmt2, ilo2p, rho1MTplus2(:, :, :, :) )
          end if

          ! Add gradient of the full density for displaced atoms
          ! for the real first variation of the density because this term cancels away during the self-consistent calculation
          ! of the first-order density but we have to consider it when plotting the density in the end and using it for setting up
          ! the Dynamical Matrix in the end.
          ! Adding contributions to density variations not dependent on z1

          !if (.false.) then
            !open(110,file='000_rho1x_pw',form='FORMATTED',action='WRITE',status='replace')
            !open(109,file='000_rho1goodx_mt',form='FORMATTED',action='WRITE',status='replace')
            !open(111,file='000_rho1y_pw',form='FORMATTED',action='WRITE',status='replace')
            !open(112,file='000_rho1goody_mt',form='FORMATTED',action='WRITE',status='replace')
            !open(113,file='000_rho1z_pw',form='FORMATTED',action='WRITE',status='replace')
            !open(114,file='000_rho1goodz_mt',form='FORMATTED',action='WRITE',status='replace')
          !end if

          do idir = 1, 3
            do iG = 1, ngpqdp
              ! IR core contribution correction to the density variation
              ! We have a minus in the calculation of rho1IRctC
              !!!.FALSE.!!!
              rho1IRDSplus(iG, idir) = rho1IRDSplus(iG, idir)! + rho1IRctC(iG, iDatom, idir)
              !if (.false.) then
                !if (idir.eq.3) then
                 ! write(110,*) iG
                  !write(110,*) real(rho1IRDSplus(iG, 1)), aimag(rho1IRDSplus(iG, 1))
                  !write(111,*) iG
                  !write(111,*) real(rho1IRDSplus(iG, 2)), aimag(rho1IRDSplus(iG, 2))
                  !write(113,*) iG
                  !write(113,*) real(rho1IRDSplus(iG, 3)), aimag(rho1IRDSplus(iG, 3))
                !end if
              !end if
            end do ! iG

            iatom = 0
            do itype = 1, atoms%ntype
              do ieqat = 1, atoms%neq(itype)
                iatom = iatom + 1
                if ( iDatom == iatom ) then
                  do lm = 2, 4
                    do imesh = 1, atoms%jri(itype)
                      ! The correction of the core correction at the displaced atom. We need not to correct the core contribution of
                      ! the displaced atom as it is the reason for the core corrections in the non-displaced atoms
                      !!!.FALSE.!!!
                      rho1MTplus( imesh, lm, iatom, idir) = rho1MTplus( imesh, lm, iatom, idir) &
                                                                                      & - 0.0*rho1MTCoreDispAt( imesh, lm, iDatom, idir)
                    end do ! imesh
                  end do ! lm
                end if ! iDatom == iatom
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                    do imesh = 1, atoms%jri(itype)
                      ! MT core contribution correction to the density variation
                      ! we have a minus in the calculation of rho1MTctC
                      !!!.FALSE.!!!
                      rho1MTplus(imesh, lm, iatom, idir) = rho1MTplus(imesh, lm, iatom, idir) &
                                                                                      & + 0.0*rho1MTctC( imesh, lm, iatom, idir, iDatom)
                      !if (.false.) then
!                        if (idir.eq.3) then
!                          write(109,*) imesh, lm
!                          write(109,*) real(rho1MTplus(imesh, lm, iatom, 1)), aimag(rho1MTplus(imesh, lm, iatom, 1))
!                          write(109,*) real(rho1MTplus(imesh, lm, iatom, 1) - grRho0MT(imesh, lm, iatom, 1)), &
!                                    & aimag(rho1MTplus(imesh, lm, iatom, 1) - grRho0MT(imesh, lm, iatom, 1))
!                          write(112,*) imesh, lm
!                          write(112,*) real(rho1MTplus(imesh, lm, iatom, 2)), aimag(rho1MTplus(imesh, lm, iatom, 2))
!                          write(112,*) real(rho1MTplus(imesh, lm, iatom, 2) - grRho0MT(imesh, lm, iatom, 2)), &
!                                    & aimag(rho1MTplus(imesh, lm, iatom, 2) - grRho0MT(imesh, lm, iatom, 2))
!                          write(114,*) imesh, lm
!                          write(114,*) real(rho1MTplus(imesh, lm, iatom, 3)), aimag(rho1MTplus(imesh, lm, iatom, 3))
!                          write(114,*) real(rho1MTplus(imesh, lm, iatom, 3) - grRho0MT(imesh, lm, iatom, 3)), &
!                                    & aimag(rho1MTplus(imesh, lm, iatom, 3) - grRho0MT(imesh, lm, iatom, 3))
!                        end if
                      !end if
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! ieqat
            end do ! itype
          end do ! idir

          !if (.false.) then
            !close(109)
            !close(110)
            !close(111)
            !close(112)
            !close(113)
            !close(114)
          !end if

          ! After the Sternheimer equation is converged we jump out of the convergence loop having written out the effective
          ! potential variation, the change of the Kohn--Sham wavefunctions and the variation of the density to hard disk
          if ( sternFinIt ) then
            write(*, '(a)') 'Sternheimer self-consistency cycle terminated!'
            exit
          end if

          do idir = 1, 3
            ! If a direction might converge faster than the other, it should not be converged any further so that broyden does not
            ! move out of the minimum due to numerical fluctuations
            if (converged(idir)) then

              if (idir == 1) then
                write(789, '(i3,2x,i3,2(f20.5),1x)', advance='no') iDatom, iqpt, lastDistance(idir, iDatom, iqpt)
              else if (idir == 2) then
                write(789, '(2(f20.5),1x)', advance='no') lastDistance(idir, iDatom, iqpt)
              else if (idir == 3) then
                write(789, '(2(f20.5))') lastDistance(idir, iDatom, iqpt)
              end if

              cycle

            end if

            call UpdNCheckDens( atoms, stars, cell, input, ngpqdp, stern1stIt, stern2ndIt, sternRegIt, converged(idir), iDatom, &
              & iqpt, idir, gpqdp, lastDistance, rho1IRDSplus(:, idir), rho1MTplus(:, :, :, idir) )

          end do ! idir

          if (all(converged)) then
            write(*, '(a)') 'Converged linear density variation. Final run of Sternheimer self-consistency cycle!'
            sternFinIt = .true.
          end if

          ! Output of density variation for maintenance reasons
          if (.false.) then
            counter = counter + 1
            write(filename1, '(a1,i1)') 'i',counter
            open( 1000, file=filename1, status='replace', action='write', form='formatted')
            do idir = 1, 3
              do iG = 1,ngpqdp
                write(1000, '(2i8,2f15.8)') idir, iG, rho1IRDSplus(iG, idir)
              end do
            end do
            close(1000)
            write(filename1, '(a1,i1)') 'm',counter
            open( 1000, file=filename1, status='replace', action='write', form='formatted')
            do idir = 1, 3
              iatom = 0
              do itype = 1, atoms%ntype
                do ieqat = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  do oqn_l = 0, atoms%lmax(itype)
                    do mqn_m = -oqn_l, oqn_l
                      lm = oqn_l * (oqn_l + 1) + mqn_m + 1! todo does everythink start at 1?
                      do imesh = 1, atoms%jri(itype)
                        ! If we are at the displaced atom we have to add the gradient of the density
                        write(1000, '(4i8,2f15.8)') idir, iatom, lm, imesh, rho1MTplus(imesh, lm, iatom, idir)
                      end do
                    end do
                  end do
                end do
              end do
            end do
            close(1000)
          end if

        end do ! self-consistency-loop

        ! Close file giving out the distances
        if ( .not.oneSternhCycle ) then
          close(789)
        end if

        ! Having converged the Sternheimer equation for one displaced atom we assign this result to the collection variable.
        ! This might be one variable in the end for performance reasons, but for sake of security the first implementation is like
        ! this
        do idir = 1, 3
          do iG = 1, ngpqdp
            rho1IRDS(iG, idir, iDatom) = rho1IRDSplus(iG, idir)
          end do ! iG
        end do ! idir

        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              if (iatom == iDatom) then
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                    do imesh = 1, atoms%jri(itype)
                      ! If we are  the displaced atom we have to add the gradient of the density
                      rho1MTDelta(imesh, lm, iatom, idir, iDatom) = rho1MTplus(imesh, lm, iatom, idir)
                      rho1MT(imesh, lm, iatom, idir, iDatom) =  rho1MTplus(imesh, lm, iatom, idir) &
                        &- grRho0MT(imesh, lm, iDatom, idir)
                      rho1MTz0(imesh, lm, idir, iDatom) =  rho1MTplus2(imesh, lm, iatom, idir)
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              else
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                    do imesh = 1, atoms%jri(itype)
                      rho1MT(imesh, lm, iatom, idir, iDatom) = rho1MTplus(imesh, lm, iatom, idir)
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end if ! Displaced atom or not
            end do ! ieqat
          end do ! itype
        end do ! idir


        !Continuity of density variation
        call checkjuPhDens1( atoms, cell, ngpqdp, rho1IRDS(:, :, iDatom), rho1MT(:, :, :, :, iDatom), qpts%bk(1:3, iqpt), gpqdp, &
          & noPts2chCont, logUnit )

      end do ! iDeqat
    end do ! iDtype

  end subroutine solveSternheimerSCC

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Updates potential for next Sternheimer cycle
  !>
  !> @param[in]   atoms          : Atoms type, see types.f90.
  !> @param[in]   stars          : Stars type, see types.f90.
  !> @param[in]   cell           : Unit cell type, see types.f90.
  !> @param[in]   dimens         : Dimension type, see types.f90.
  !> @param[in]   iDatom         : Index of displaced atom
  !> @param[in]   iDtype         : Index of atom type to which displaced atom belongs
  !> @param[in]   first_run      : Switch indicating that the Sternheimer self-consistency cycle is in its first run resulting in
  !>                               generation of the adequate potentials (external and not effective) including a switch-off of the
  !>                               calculation of terms that cancel anyway.
  !> @param[in]   final_run      : Switch indicating that the Sternheimer self-consistency cycle is in its final run resulting in
  !>                               generation of the adequate potentials (external and not effective) including a switch-off of the
  !>                               calculation of terms that cancel anyway.
  !> @param[in]   ngdp           : Number of G-vectors for potentials and densities.
  !> @param[in]   iqpt           : Index of q-point in kpts file that is evaluated
  !> @param[in]   ngpqdp         : Number of G-vectors for shifted G-set for q-point with index iqpt.
  !> @param[in]   noPts2chCont   : Number of points for which continuity test should be performed.
  !> @param[in]   logUnit        : Unit number for juPhon.log.
  !> @param[in]   grRho0IR       : Plane-wave coefficients of the gradient of the interstitial unperturbed density
  !> @param[in]   grRho0MT       : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[in]   qpoint         : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]   rho1PW         : Plane-wave interstitial expansion coefficients of the linear density variation.
  !> @param[in]   rho1MT         : Spherical harmonic muffin-tin expansion coefficients of the linear density variation without the
  !>                               muffin-tin gradient of the density
  !> @param[in]   ylm            : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to
  !>                               lmax created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   gWghts         : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   dKernMTGPts    : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                               respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]   vxc1IRKern     : Plane-wave interstitial coefficients of the functional derivative of the xc-kernel with respect
  !>                               to the density
  !> @param[in]   gdp            : G-vectors of potentials and densities
  !> @param[in]   grVeff0MT_init : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed external potential
  !>                               (part of the effective) for initial Sternheimer iteration created in m_jpGrVeff0::gengrveff0
  !> @param[in]   grVeff0MT_main : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed effective
  !>                               potential for a regular and the final Sternheimer iteration created in m_jpGrVeff0::gengrveff0
  !> @param[in]   gpqdp          : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   rho0IRpw       : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]   rho0MTsh       : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[out]  sumVMTs        : Sum of the spherical harmonic expansion coefficients for the gradient of the unperturbed
  !>                               external(first iteration) or effective(regular/final iteration) potential and its linear
  !>                               variation. This potential is defined in 7.105b (dissertation CRG).
  !> @param[out]  sumVMTs2       : Sum of the spherical harmonic expansion coefficients for the gradient of the unperturbed
  !>                               external(first iteration) or effective(regular/final iteration) potential and its linear
  !>                               variation. This potential is defined in 7.105b (dissertation CRG) and decorated with a factor as
  !>                               pointed out around 7.108 (dissertation CRG).
  !> @param[out]  vExt1IR_final  : Plane-wave coefficients of the interstitial converged external linear potential variation
  !> @param[out]  vHar1IR_final  : Plane-wave coefficients of the interstitial converged Hartree linear potential variation
  !> @param[out]  vExt1MT        : Spherical harmonic expansion coefficients of the complete converged external muffin-tin linear
  !>                               potential variation
  !> @param[out]  vExt1MTDelta   : Spherical harmonic expansion coefficients of the difference of the linear variation of the
  !>                               external potential between a finite q and q=0; only surface contribution no volume term.
  !> @param[out]  vExt1MTq0      : Spherical harmonic expansion coefficients of the full linear variation of the external
  !>                               potential for q = 0
  !> @param[out]  vHar1MTDelta   : Spherical harmonic expansion coefficients of the difference of the linear variation of the
  !>                               Hartree potential between a finite q and q=0; only surface contribution no volume term.
  !> @param[out]  vHar1MTq0      : Spherical harmonic expansion coefficients of the full linear variation of the Hartree
  !>                               potential for q = 0
  !> @param[out]  vXc1MTDelta    : Spherical harmonic expansion coefficients of the difference of the linear variation of the
  !>                               external potential between a finite q and q=0; only surface contribution no volume term.
  !> @param[out]  vXc1MTq0       : Spherical harmonic expansion coefficients of the full linear variation of the external
  !>                               potential for q = 0
  !> @param[out]  vEff1IR_final  : Plane-wave coefficients of the interstitial converged effective linear potential variation
  !> @param[out]  vEff1MT_final  : Spherical harmonic expansion coefficients of the complete converged effective muffin-tin linear
  !>                               potential variation
  !> @param[out]  vHar1MT_final  : Spherical harmonic expansion coefficients of the complete converged Hartree muffin-tin linear
  !>                               potential variation
  !> @param[out]  vEff1IRsh      : Plane-wave coefficients of the interstitial effective linear variation of the (in first iteration
  !>                               the external) effective potential
  !>
  !> @note : For determining the final quantities the final_run switch must be set true.
  !> @details
  !> See also Section 7.4.5
  !>
  !> @ todo add spin
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine UpdPots( atoms, input, stars, cell, iDatom, iDtype, first_run, ngdp, grRho0IR, grRho0MT, qpoint, rho1PW, rho1MT, gdp,&
      & grVeff0MT_init, grVeff0MT_main, sumVMTs, ylm, gWghts, dKernMTGPts, vxc1IRKern, iqpt, sumVMTs2, final_run, vExt1MT, &
      & vEff1MT_final, vExt1IR_final, vHar1IR_final, vHar1MT_final, vExt1MTDelta, vExt1MTq0, vHar1MTDelta,  vHar1MTq0, vXc1MTDelta,&
      & vXc1MTq0, rho0IRpw, rho0MTsh, vEff1IR_final, ngpqdp, gpqdp, vEff1IRsh, noPts2chCont, logUnit, vEff1MTnoVol, &
      & vExt1noqIR_final, vCoul1IRtempNoVol, vCoul1MTtempNoVol )

    use m_types, only : t_atoms, t_input, t_stars, t_cell
    use m_jpVeff1, only : GenVeff1

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_stars),                  intent(in)  :: stars
    type(t_cell),                   intent(in)  :: cell

    ! Scalar parameters
    integer,                        intent(in)  :: iDatom
    integer,                        intent(in)  :: iDtype
    logical,                        intent(in)  :: first_run
    logical,                        intent(in)  :: final_run
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: iqpt
    integer,                        intent(in)  :: ngpqdp
    integer,                        intent(in)  :: noPts2chCont
    integer,                        intent(in)  :: logUnit

    ! Array parameters
    complex,                        intent(in)  :: grRho0IR(:, :)
    complex,                        intent(in)  :: grRho0MT(:, :, :, :)
    real,                           intent(in)  :: qpoint(:)
    complex,                        intent(in)  :: rho1PW(:, :)
    complex,                        intent(in)  :: rho1MT(:, :, :, :)
    complex,                        intent(in)  :: ylm(:, :)
    real,                           intent(in)  :: gWghts(:) ! gaussian weights belonging to gausPts
    real,                           intent(in)  :: dKernMTGPts(:, :, :)
    complex,                        intent(in)  :: vxc1IRKern(:)
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: grVeff0MT_init(:, :, :, :)
    complex,                        intent(in)  :: grVeff0MT_main(:, :, :, :)
    integer,                        intent(in)  :: gpqdp(:, :)
    complex,                        intent(in)  :: rho0IRpw(:, :)
    complex,                        intent(in)  :: rho0MTsh(:, :, :, :)
    complex,                        intent(out) :: sumVMTs(:, :, :, :)
    complex,                        intent(out) :: sumVMTs2(:, :, :, :)
    complex,                        intent(out) :: vExt1IR_final(:, :, :)
    complex,                        intent(out) :: vHar1IR_final(:, :, :)
    complex,                        intent(out) :: vExt1MT(:, :, :, :, :)
    complex,                        intent(out) :: vEff1MTnoVol(:, :, :, :, :)
    complex,                        intent(out) :: vExt1MTDelta(:, :, :, :, :)
    complex,                        intent(out) :: vExt1MTq0(:, :, :, :, :)
    complex,                        intent(out) :: vHar1MTDelta(:, :, :, :, :)
    complex,                        intent(out) :: vHar1MTq0(:, :, :, :, :)
    complex,                        intent(out) :: vXc1MTDelta(:, :, :, :, :)
    complex,                        intent(out) :: vXc1MTq0(:, :, :, :, :)
    complex,                        intent(out) :: vEff1IR_final(:, :, :)
    complex,                        intent(out) :: vExt1noqIR_final(:, :, :)
    complex,                        intent(out) :: vEff1MT_final(:, :, :, :, :)
    complex,                        intent(out) :: vHar1MT_final(:, :, :, :, :)
    complex,           allocatable, intent(out) :: vCoul1IRtempNoVol(:, :)
    complex,           allocatable, intent(out) :: vCoul1MTtempNoVol(:, :, :, :)
    complex,           allocatable, intent(out) :: vEff1IRsh(:, :)

    ! Scalar variables
    logical                                     :: harSw
    logical                                     :: extSw
    logical                                     :: xcSw
    logical                                     :: vExtFull
    logical                                     :: vHarNum
    integer                                     :: ieqat
    integer                                     :: iatom
    integer                                     :: idir
    integer                                     :: lm
    integer                                     :: itypeloc
    integer                                     :: iG
    integer                                     :: imesh
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: lm2
    integer                                     :: itype

    ! Array variables
    complex,           allocatable              :: vEff1IR(:, :)
    complex,           allocatable              :: vEff1MT(:, :, :, :)
    complex,           allocatable              :: vHar1IR(:, :)
    complex,           allocatable              :: vHar1MT(:, :, :, :)
    complex,           allocatable              :: vExt1IRtemp(:, :)
    complex,           allocatable              :: vExt1MTtemp(:, :, :, :)
    complex,           allocatable              :: vExt1IRtempNoVol(:, :)
    complex,           allocatable              :: vExt1MTtempNoVol(:, :, :, :)
    complex,           allocatable              :: vExt1IRtempNoVolnoq(:, :)
    complex,           allocatable              :: vExt1MTtempNoVolnoq(:, :, :, :)
    complex,           allocatable              :: vH1IRtempNoVol(:, :)
    complex,           allocatable              :: vH1MTtempNoVol(:, :, :, :)
    complex,           allocatable              :: rho1PWzero(:, :)
    complex,           allocatable              :: rho1MTzero(:, :, :, :)
    real                                        :: Gext(3)


    ! Init dummy variabels for intent(out)
    sumVMTs(:, :, :, :)          = cmplx(0., 0.)
    sumVMTs2(:, :, :, :)         = cmplx(0., 0.)
    vExt1MT(:, :, :, :, :)       = cmplx(0., 0.)
    vEff1MTnoVol(:, :, :, :, :)  = cmplx(0., 0.)
    vExt1MTDelta(:, :, :, :, :)  = cmplx(0., 0.)
    vExt1MTq0(:, :, :, :, :)     = cmplx(0., 0.)
    vHar1MTDelta(:, :, :, :, :)  = cmplx(0., 0.)
    vHar1MTq0(:, :, :, :, :)     = cmplx(0., 0.)
    vXc1MTDelta(:, :, :, :, :)   = cmplx(0., 0.)
    vXc1MTq0(:, :, :, :, :)      = cmplx(0., 0.)
    vExt1IR_final(:, :, :)       = cmplx(0., 0.)
    vExt1noqIR_final(:, :, :)       = cmplx(0., 0.)
    vHar1IR_final(:, :, :)       = cmplx(0., 0.)
    vEff1IR_final(:, :, :)       = cmplx(0., 0.)
    vEff1MT_final(:, :, :, :, :) = cmplx(0., 0.)
    vHar1MT_final(:, :, :, :, :) = cmplx(0., 0.)

    if ( first_run ) then

      write(*, '(a)') 'Setting up potential variations for first Sternheimer iteration...'

      ! Calculate the first-order external potential for initial self-consistency cyle
      harSw = .false.
      extSw = .true.
      xcSw = .false.
      !if (.false..or..FALSE.) then
        vExtFull = .true.
      !else
      !  vExtFull = .false.
      !end if
      vHarNum = .false.

      if (.false.) then
        open(109,file='000_V1in_mt',form='FORMATTED',action='WRITE',status='replace')
        open(110,file='000_V1in_pw',form='FORMATTED',action='WRITE',status='replace')
        open(111,file='000_Gqvec',form='FORMATTED',action='WRITE',status='replace')
      end if

      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vEff1IRsh, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )

      if (.false.) then
        close(109)
        close(110)
        close(111)
      end if

      write(logUnit, '(a)') 'Continuity of external potential variation for 1st iteration of Sternheimer SCC'
      write(logUnit, '(a)') '-------------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngpqdp, vEff1IRsh, vEff1MT, qpoint, gpqdp, noPts2chCont, logUnit  )

      ! Add the muffin-tin external linear potential variation with the gradient of the unperturbed external potential.
      ! As the tlmplm routines from FLEUR were recycled, that only fill (due to symmetry) one off-diagonal part of the Hamiltonian
      ! due to symmetry, for the non-symmetric potential, its complex expansion coefficient have to be adjusted afterwards so that
      ! the tlmplm routine can be used.
      sumVMTs = cmplx(0.0, 0.0)
      sumVMTs2 = cmplx(0.0, 0.0)
      iatom = 0

      if (.false.) then
        open(109,file='000_V1inGr_mt',form='FORMATTED',action='WRITE',status='replace')
      end if

      do itypeloc = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itypeloc)
          iatom = iatom + 1
          do idir = 1, 3
             do oqn_l = 0, atoms%lmax(itypeloc)
               do mqn_m = -oqn_l, oqn_l
                 lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                 lm2 = oqn_l * (oqn_l + 1) + 1 - mqn_m
                 do imesh = 1, atoms%jri(itypeloc)
                   sumVMTs(imesh, lm, idir, iatom) = vEff1MT(imesh, lm, iatom, idir) + grVeff0MT_init(imesh, lm, idir, iatom)

                   if (.false.) then
                     if (idir.eq.3) then
                       write(109,*) imesh, lm, iatom, 1
                       write(109,*) real(sumVMTs(imesh, lm, 1, iatom)), aimag(sumVMTs(imesh, lm, 1, iatom))
                       write(109,*) imesh, lm, iatom, 2
                       write(109,*) real(sumVMTs(imesh, lm, 2, iatom)), aimag(sumVMTs(imesh, lm, 2, iatom))
                       write(109,*) imesh, lm, iatom, 3
                       write(109,*) real(sumVMTs(imesh, lm, 3, iatom)), aimag(sumVMTs(imesh, lm, 3, iatom))
                     end if
                   end if
                   sumVMTs2(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTs(imesh, lm, idir, iatom))
                 end do
               end do
             end do
          end do
        end do
      end do

      if (.false.) then
        close(109)

      !do idir = 1, 3
      !  iatom = 0
      !  do itype = 1, atoms%ntype
      !    do ieqat = 1, atoms%neq(itype)
      !      iatom = iatom + 1
      !      do lm = 1, (atoms%lmax(itype) + 1)**2
      !        do imesh = 1, atoms%jri(itype)
      !          vEff1MT(imesh, lm, iatom, idir) = vEff1MT(imesh, lm, iatom, idir) + Vxc1MT(imesh, lm, iatom, idir)
      !          if (idir.eq.3) then
      !             write(109,*) imesh, lm, real(vEff1MT(imesh, lm, iatom, 1)), aimag(vEff1MT(imesh, lm, iatom, 1))
      !             write(109,*) imesh, lm, real(vEff1MT(imesh, lm, iatom, 2)), aimag(vEff1MT(imesh, lm, iatom, 2))
      !             write(109,*) imesh, lm, real(vEff1MT(imesh, lm, iatom, 3)), aimag(vEff1MT(imesh, lm, iatom, 3))
      !          end if
      !        end do
      !      end do
      !    end do
      !  end do
      !end do

      end if

    else if ( final_run ) then

      write(*, '(a)') 'Setting up potential variations for final run of Sternheimer iteration and dynamical matrix...'

      ! "Update" converged first-order effective potential for final self-consistency cycle
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      vExtFull = .false.!true for Alex
      vHarNum = .false.

      !if (.false.) then
        !open(109,file='000_V1fin_mt',form='FORMATTED',action='WRITE',status='replace')
        open(110,file='000_V1x_pw',form='FORMATTED',action='WRITE',status='replace')
        open(113,file='000_V1y_pw',form='FORMATTED',action='WRITE',status='replace')
        open(114,file='000_V1z_pw',form='FORMATTED',action='WRITE',status='replace')
      !end if
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW,  rho1MT, &
        & grRho0MT, gdp, vEff1IRsh, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )
      !if (.false.) then
        !close(109)
        close(110)
        close(113)
        close(114)
      !end if
      write(logUnit, '(a)') 'Continuity of effective potential variation for final iteration of Sternheimer SCC'
      write(logUnit, '(a)') '----------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngpqdp, vEff1IRsh, vEff1MT, qpoint, gpqdp, noPts2chCont, logUnit  )

      ! Add the muffin-tin effective linear potential variation with the gradient of the unperturbed effective potential.
      ! As the tlmplm routines from FLEUR were recycled, that only fill (due to symmetry) one off-diagonal part of the Hamiltonian
      ! due to symmetry, for the non-symmetric potential, its complex expansion coefficient have to be adjusted afterwards so that
      ! the tlmplm routine can be used.
      sumVMTs = cmplx(0.0, 0.0)
      sumVMTs2 = cmplx(0.0, 0.0)
      iatom = 0
      !if (.false.) then
        !open(109,file='000_V1goodx_mt',form='FORMATTED',action='WRITE',status='replace')
        !open(110,file='000_V1x_mt',form='FORMATTED',action='WRITE',status='replace')
        !open(111,file='000_V1goody_mt',form='FORMATTED',action='WRITE',status='replace')
        !open(112,file='000_V1y_mt',form='FORMATTED',action='WRITE',status='replace')
        !open(113,file='000_V1goodz_mt',form='FORMATTED',action='WRITE',status='replace')
        !open(114,file='000_V1z_mt',form='FORMATTED',action='WRITE',status='replace')
      !end if
      do itypeLoc = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itypeLoc)
          iatom = iatom + 1
          do idir = 1, 3
             do oqn_l = 0, atoms%lmax(itypeloc)
               do mqn_m = -oqn_l, oqn_l
                 lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                 lm2 = oqn_l * (oqn_l + 1) + 1 - mqn_m
                 do imesh = 1, atoms%jri(itypeloc)
                   sumVMTs(imesh, lm, idir, iatom) = vEff1MT(imesh, lm, iatom, idir) + grVeff0MT_main(imesh, lm, idir, iatom)
                   !if (.false.) then
                     !if (idir.eq.3) then
                      ! write(109,*) imesh, lm
                      ! write(109,*) real(sumVMTs(imesh, lm, 1, iatom)), aimag(sumVMTs(imesh, lm, 1, iatom))
                      ! write(109,*) real(vEff1MT(imesh, lm, iatom, 1)), aimag(vEff1MT(imesh, lm, iatom, 1))
                      ! write(111,*) imesh, lm
                      ! write(111,*) real(sumVMTs(imesh, lm, 2, iatom)), aimag(sumVMTs(imesh, lm, 2, iatom))
                      ! write(111,*) real(vEff1MT(imesh, lm, iatom, 2)), aimag(vEff1MT(imesh, lm, iatom, 2))
                      ! write(113,*) imesh, lm
                      ! write(113,*) real(sumVMTs(imesh, lm, 3, iatom)), aimag(sumVMTs(imesh, lm, 3, iatom))
                      ! write(113,*) real(vEff1MT(imesh, lm, iatom, 3)), aimag(vEff1MT(imesh, lm, iatom, 3))
                     !end if
                   !end if
                   sumVMTs2(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTs(imesh, lm, idir, iatom))
                   vEff1MTNoVol(imesh, lm, iatom, idir, iDatom) = vEff1MT(imesh, lm, iatom, idir) + grVeff0MT_main(imesh, lm, idir, iatom)
                 end do
              end do
            end do
          end do
        end do
      end do
      deallocate( vEff1MT )
      !if (.false.) then
        !close(109)
        !close(110)
        !close(111)
        !close(112)
        !close(113)
        !close(114)
        !NOstopNO!convergent NOstopNO
      !end if
      !NOstopNO!convergent NOstopNO

      ! Create Coulomb potential for Dynamical Matrix Surface  contribution.
      harSw = .true.
      extSw = .true.
      xcSw = .false.
      vExtFull = .false.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vCoul1IRtempNoVol, vCoul1MTtempNoVol, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )

      ! Create effective potential for Pulay Dynamical Matrix contribution.
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      vExtFull = .true.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vEff1IR, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )

      write(logUnit, '(a)') 'Continuity of effective potential variation for Pulay dynamical matrix'
      write(logUnit, '(a)') '----------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngpqdp, vEff1IR, vEff1MT, qpoint, gpqdp, noPts2chCont, logUnit  )

      ! Create Hartree potential
      harSw = .true.
      extSw = .false.
      xcSw = .false.
      vExtFull = .true.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vHar1IR, vHar1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )

      write(logUnit, '(a)') 'Continuity of Hartree potential variation for surface dynamical matrix'
      write(logUnit, '(a)') '----------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngpqdp, vHar1IR, vHar1MT, qpoint, gpqdp, noPts2chCont, logUnit  )

      ! Create external potential
      harSw = .false.
      extSw = .true.
      xcSw = .false.
      ! Actually we do not need it really
      allocate(rho1PWzero(ngpqdp, 3))
      allocate( rho1MTzero(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3))
      rho1PWzero = cmplx(0.0, 0.0)
      rho1MTzero = cmplx(0.0, 0.0)
      ! There is nothing canceling away as in Sternheimer therefore we need the full contribution of the external potential
      vExtFull = .true.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PWzero, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, &
        & gpqdp, vHarNum ) ! add spin

      ! Create external potential without volume term in MT
      harSw = .false.
      extSw = .true.
      xcSw = .false.
      vExtFull = .false.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PWzero, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtempNoVol, vExt1MTtempNoVol, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, &
        & iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

      !!!More terms needed for A.N. representation:

      harSw = .true.
      extSw = .false.
      xcSw = .false.
      vExtFull = .false.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vH1IRtempNoVol, vH1MTtempNoVol, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, &
        & iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

      harSw = .false.
      extSw = .true.
      xcSw = .false.
      vExtFull= .false.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, rho1PWzero, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtempNoVolnoq, vExt1MTtempNoVolnoq, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1, ngdp, gdp, &
        & vHarNum ) ! add spin


      write(logUnit, '(a)') 'Continuity of external potential variation for Hellmann-Feynman dynamical matrix'
      write(logUnit, '(a)') '--------------------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngpqdp, vExt1IRtemp, vExt1MTtemp, qpoint, gpqdp, noPts2chCont, logUnit  )

      do idir = 1, 3
        do iG = 1, ngpqdp
          vEff1IR_final(iG, idir, iDatom) = vEff1IR(iG, idir)
          ! We need the unwarped Hartree potential
          vHar1IR_final(iG, idir, iDatom) = vHar1IR(iG, idir)
          vExt1IR_final(iG, idir, iDatom) = vExt1IRtemp(iG, idir)
        end do
        do iG = 1, ngdp
          vExt1noqIR_final(iG, idir, iDatom) = vExt1IRtempNoVolnoq(iG, idir)
        end do
      end do

      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vHar1MT_final(imesh, lm, iatom, idir, iDatom) = vHar1MT(imesh, lm, iatom, idir)
                  vExt1MT(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                  vEff1MT_final(imesh, lm, iatom, idir, iDatom) = vEff1MT(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      deallocate(rho1PWzero)
      allocate(rho1PWzero(ngdp, 3))
      rho1PWzero = cmplx(0.0, 0.0)
      rho1MTzero = cmplx(0.0, 0.0)
      vExt1MTDelta(:, :, :, :, :) = cmplx(0., 0.)
      deallocate(vExt1IRtemp, vExt1MTtemp)
      harSw = .false.
      extSw = .true.
      xcSw = .false.
      vExtFull= .false.
      vHarNum = .false.
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PWzero, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, &
        & gpqdp, vHarNum )

      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vExt1MTDelta(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      harSw = .false.
      extSw = .true.
      xcSw = .false.
      vExtFull= .false.
      vHarNum = .false.
      deallocate(vExt1IRtemp, vExt1MTtemp)
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, rho1PWzero, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1, ngdp, gdp, &
        & vHarNum ) ! add spin

      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vExt1MTDelta(imesh, lm, iatom, idir, iDatom) = vExt1MTDelta(imesh, lm, iatom, idir, iDatom) &
                                                                                             & - vExt1MTtemp(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      vExt1MTq0 = cmplx(0., 0.)
      if (iqpt /= 1) then

        harSw = .false.
        extSw = .true.
        xcSw = .false.
        vExtFull= .true.
        vHarNum = .false.
        deallocate(vExt1IRtemp, vExt1MTtemp)
        call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, vExtFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, rho1PWzero,&
          & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1, ngdp, gdp,&
          & vHarNum ) ! add spin

        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    vExt1MTq0(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      vHar1MTDelta(:, :, :, :, :) = cmplx(0., 0.)

      ! Create Hartree potential snippet in MT contribution for q = 0
      harSw = .true.
      extSw = .false.
      xcSw = .false.
      vExtFull = .false.
      vHarNum = .true.
      deallocate(vExt1IRtemp, vExt1MTtemp)
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, -grRho0IR, &
        & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1,  ngdp, gdp,&
        & vHarNum )
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vHar1MTDelta(imesh, lm, iatom, idir, iDatom) = -vExt1MTtemp(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      ! Create rest of delta V_H
      harSw = .true.
      extSw = .false.
      xcSw = .false.
      vExtFull = .false.
      vHarNum = .false.
      deallocate(vExt1IRtemp, vExt1MTtemp)
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT,  gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, &
        & vHarNum )
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vHar1MTDelta(imesh, lm, iatom, idir, iDatom) = vHar1MTDelta(imesh, lm, iatom, idir, iDatom) &
                                                                                             & + vExt1MTtemp(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      vHar1MTq0(:, :, :, :, :) = cmplx(0., 0.)
      if (iqpt /= 1) then
        ! Create Hartree for q = 0
        harSw = .true.
        extSw = .false.
        xcSw = .false.
        vExtFull = .true.
        vHarNum = .false.
        deallocate(vExt1IRtemp, vExt1MTtemp)
        call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, -grRho0IR,&
          & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1, ngdp, &
          & gdp, vHarNum )


        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    vHar1MTq0(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      ! Create rest of delta V_xc
      vXc1MTDelta(:, :, :, :, :) = cmplx(0., 0.)
      harSw = .false.
      extSw = .false.
      xcSw = .true.
      vExtFull = .false.
      vHarNum = .false.
      deallocate(vExt1IRtemp, vExt1MTtemp)
      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, &
        & vHarNum )
      do idir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  vXc1MTDelta(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                end do
              end do
            end do
          end do
        end do
      end do

      vXc1MTq0(:, :, :, :, :) = cmplx(0., 0.)
      if (iqpt /= 1) then
        ! Create rest of delta V_xc
        harSw = .false.
        extSw = .false.
        xcSw = .true.
        vExtFull = .true.
        vHarNum = .false.
        deallocate(vExt1IRtemp, vExt1MTtemp)
        call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, [0., 0., 0.], rho0IRpw, rho0MTsh, rho1PWZero,&
          & rho1MTzero, grRho0MT, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, 1, ngdp, gdp,&
          & vHarNum )
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    vXc1MTq0(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir)
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

      !todo where is q=0 part of vxcdelta?

    else
      write(*, '(a)') 'Updating effective potential variation with new density variation for regular run.'
      ! Update first-order effective potential for regular self-consistency cycle to further converge
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      vExtFull = .false.!.true. for Alex
      vHarNum = .false.

      if (.false.) then
        open(109,file='000_V1out_mt',form='FORMATTED',action='WRITE',status='replace')
        open(112,file='000_V1out_pw',form='FORMATTED',action='WRITE',status='replace')
      end if

      call GenVeff1( input, stars, cell, atoms, harSw, extSw, xcSw, VextFull, ngdp, qpoint, rho0IRpw, rho0MTsh, rho1PW, rho1MT, &
        & grRho0MT, gdp, vEff1IRsh, vEff1MT, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum )

      if (.false.) then
        close(109)
        close(110)
      end if

      write(logUnit, '(a)') 'Continuity of effective potential variation for regular iteration of Sternheimer SCC'
      write(logUnit, '(a)') '--------------------------------------------------------------------------------'
      call checkjuPhDens1(atoms, cell, ngdp, vEff1IRsh, vEff1MT, qpoint, gdp, noPts2chCont, logUnit)

      sumVMTs = cmplx(0.0, 0.0)
      sumVMTs2 = cmplx(0.0, 0.0)
      iatom = 0
      if (.false.) then
        open(109,file='000_V1outGr_mt',form='FORMATTED',action='WRITE',status='replace')
      end if
      do itypeLoc = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itypeLoc)
          iatom = iatom + 1
          do idir = 1, 3
             do oqn_l = 0, atoms%lmax(itypeloc)
               do mqn_m = -oqn_l, oqn_l
                 lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                 lm2 = oqn_l * (oqn_l + 1) + 1 - mqn_m
                 do imesh = 1, atoms%jri(itypeloc)
                   sumVMTs(imesh, lm, idir, iatom) = vEff1MT(imesh, lm, iatom, idir) + grVeff0MT_main(imesh, lm, idir, iatom)
                   if (.false.) then
                     if (idir.eq.3) then
                       write(109,*) imesh, lm, iatom, 1
                       write(109,*) real(sumVMTs(imesh, lm, 1, iatom)), aimag(sumVMTs(imesh, lm, 1, iatom))
                       write(109,*) imesh, lm, iatom, 2
                       write(109,*) real(sumVMTs(imesh, lm, 2, iatom)), aimag(sumVMTs(imesh, lm, 2, iatom))
                       write(109,*) imesh, lm, iatom, 3
                       write(109,*) real(sumVMTs(imesh, lm, 3, iatom)), aimag(sumVMTs(imesh, lm, 3, iatom))
                     end if
                   end if
                   sumVMTs2(imesh, lm2, idir, iatom) = (-1)**(mqn_m) * conjg(sumVMTs(imesh, lm, idir, iatom))
                 end do
              end do
            end do
          end do
        end do
      end do
      !if (.false.) then
        close(109)
        !NOstopNO!1st iteration NOstopNO.
      !end if
      !NOstopNO!1st iteration NOstopNO

    end if

  end subroutine UpdPots

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Generates reciprocal energy difference for next Sternheimer cycle
  !>
  !> @details
  !>
  !> @param[in]  ne         : Number of eigenvalues per k-point.
  !> @param[in]  nobd       : Number of occupied bands per k-point and spin
  !> @param[in]  eigBra     : Contains Kohn\--Sham eigenvalues at k + q ( k')
  !> @param[in]  eigKet     : Contains Kohn\--Sham eigenvalues at k
  !> @param[out] recEdiffME : Reciprocal energy difference of eigBra and eigKet
  !> @param[out] cutContr   : Deprecated will be deleted
  !>
  !> @details
  !> See Section 7.4.1 and Equation 7.94 (dissertation CRG)
  !>
  !> @todo what is with cutContr?
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine GenRecEdifME( ne, nobd, eigBra, eigKet, recEdiffME, cutContr )

    implicit none

    integer,              intent(in)  :: ne
    integer,              intent(in)  :: nobd
    real,                 intent(in)  :: eigBra(:)
    real,                 intent(in)  :: eigKet(:)
    real,                 intent(out) :: recEdiffME(:, :)
    real,                 intent(out) :: cutContr(:, :)

    integer                           :: nband
    integer                           :: pband
    real                              :: eDiff

    recEdiffMe = 0.
    cutContr = 0.

    ! Generate inverse left side of Sternheimer equation to be multiplied to the right side. It is trivial because matrix element
    ! becomes diagonal in p and m so that we get matrix element for every p and n.
    ! todo find out best minimal value and maybe make this an input parameter!!!
    do nband = 1, nobd
      do pband = 1, ne
        eDiff = eigBra(pband) - eigKet(nband)
        if ( abs( eDiff ) < 1e-12 ) then !5e-3!vl. 10e-6 bei komplizierten Systemen
          ! in first order the matrix element vanishes for degenerated eigenvalues
          recEdiffMe(pband, nband) = 0
          cutContr(pband, nband) = 0
        else
          recEdiffMe(pband, nband) = 1. / eDiff
          cutContr(pband, nband) = 1.
        end if
      end do ! pband
    end do ! nband

  end subroutine GenRecEdifME

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calls routines to set-up Sternheimer equation and solves it for first-order wavefunction expansion coefficients.
  !>
  !> @param[in]  atoms        : Atoms type, see types.f90
  !> @param[in]  sym          : Symmetries type, see types.f90.
  !> @param[in]  cell         : Unit cell type, see types.f90.
  !> @param[in]  kpts         : K-points type, see types.f90.
  !> @param[in]  qpts         : Q-point set represented by a k-points type, see types.f90
  !> @param[in]  dimens       : Dimension type, see types.f90.
  !> @param[in]  usdus        : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]  td4HS0       : Tlmplm matrix type for Sternheimer Pulay matrix elements, see types.f90
  !> @param[in]  td4V         : Tlmplm matrix type for Sternheimer Hellmann-Feynman muffin-tin matrix element, see types.f90
  !> @param[in]  td4V2        : Tlmplm matrix type for Sternheimer Hellmann-Feynman muffin-tin matrix element, decorated with
  !>                            factor according to Equation 7.108, see types.f90
  !> @param[in]  stars        : Stars type, see types.f90.
  !> @param[in]  ikpt         : Index of k-point in k-point set
  !> @param[in]  ikpq         : Index of k + q (k' backfolded) in k-point set
  !> @param[in]  iqpt         : Index of q in q-point set
  !> @param[in]  idir         : Index of displacement direction
  !> @param[in]  ngdp         : Number of G-vectors for potentials and densities
  !> @param[in]  ngpqdp       : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  iDtype       : Index of atom type to which displaced atom belongs
  !> @param[in]  iDatom       : Index of displaced atom
  !> @param[in]  maxlmp       : Maximal size of an arry dimension containing a cumultative index for l m and the basis functions for
  !>                            this l (including LOs.)
  !> @param[in]  gdp          : G-vectors of potentials and densities
  !> @param[in]  ne           : Number of eigenvalues per k-point.
  !> @param[in]  nv           : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]  eigKet       : Contains Kohn\--Sham eigenvalues
  !> @param[in]  El           : Contains LAPW and LO energy parameters.
  !> @param[in]  nRadFun      : Number of radial functions per orbital quantum number l and atom type.
  !> @param[in]  iloTable     : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given
  !> @param[in]  ilst         : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !> @param[in]  GbasVec      : G-basis vectors
  !> @param[in]  zKet         : Kohn-Sham eigenvectors at k
  !> @param[in]  zBra         : Kohn-Sham eigenvectors at k+q
  !> @param[in]  kveclo       : Basis G-vectors of local orbitals.
  !> @param[in]  vEff1IR      : Plane-wave coefficients of the interstitial effective linear variation of the (in first iteration
  !>                            the external) effective potential
  !> @param[in]  nobd         : Number of occupied bands per k-point and spin
  !> @param[in]  nlo_atom     : Contains information about number of LOs at every atom
  !> @param[in]  recEdiffME   : Reciprocal energy difference of Kohn\--Sham eigenvalues at k + q (k') and k
  !> @param[in]  kpq2kPrVec   : Backfolding vector from k + q in 2nd Brillouin zone to k' in 1st Brillouin zone
  !> @param[in]  surfIntVFast : Surface integral part to the Sternheimer equation containing the effective potential
  !> @param[in]  gpqdp        : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  eigBra       : Deprecated: Should be removed after closing #45
  !> @param[in]  cutContr     : Deprecated: Should be removed after closing #45
  !> @param[out] z1G          : Linear variation of the Kohn\--Sham eigenvectors
  !>
  !> @details
  !> This routine calculates Equation 7.95 in dissertation of CRG. See Section 7.4.1 and Equation 7.94 (dissertation CRG)
  !>
  !> @note
  !> We have to consider all (occupied and unoccupied) bands p in the bra and only the occupied bands n in the kets.
  !> Using the OEP approach of Markus Betzinger, it should be possible to also only consider the occupied bands in the bras leading
  !> to significant runtime enhancements.
  !>
  !> @todo spin
  !> @todo cutcontr and eigBra?
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine solveSternheimerEq( fmpi,   atoms, input, sym, cell, kpts, qpts, usdus, td4HS0, loosetd, td4V, loosetd1, stars,  gdp, ne, nv, vEff1IR, eigKet,&
      & eigBra, El, nRadFun, iloTable, ilst, GbasVec, zKet, zBra, kveclo, iDtype, iDatom, ikpt, ikpq, iqpt, idir, ngdp, nobd, z1G, &
      & nlo_atom, recEdiffME, kpq2kPrVec, td4V2, loosetd2, cutContr, surfIntVFast, ngpqdp, gpqdp, maxlmp, killcont, haa, dhaa, rbas1, mat_elH, mat_elS)

    use m_types
     
    use m_abcof, only : abcof
    use m_abcof3
    use m_jpSternhHF, only : calcMEPotIR, calcVsumMT
    use m_jpSternhPulaySurface, only : CalcSintKinEps, CalcHS0MT
    use m_juDFT_stop, only : juDFT_warn
    !use m_intgr, only : intgr3LinIntp
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameter
    type(t_mpi),                  intent(in)  :: fmpi
     
    type(t_atoms),                  intent(in)  :: atoms
    type(t_input),                  intent(in)  :: input
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_usdus),                  intent(in)  :: usdus
    type(t_tlmplm),                 intent(in)  :: td4HS0
    type(t_tlmplm),                 intent(in)  :: td4V
    type(t_tlmplm),                 intent(in)  :: td4V2
    type(t_stars),                  intent(in)  :: stars

    ! Scalar parameter
    integer,                        intent(in)  :: ikpt
    integer,                        intent(in)  :: ikpq
    integer,                        intent(in)  :: iqpt
    integer,                        intent(in)  :: idir
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: ngpqdp ! todo : Deprecated, should be removed when finite q working
    integer,                        intent(in)  :: iDtype
    integer,                        intent(in)  :: iDatom
    integer,                        intent(in)  :: maxlmp

    ! Array parameter
    INTEGER, INTENT(IN) :: killcont(9)
    integer,                        intent(in)  :: gdp(:, :)
    integer,                        intent(in)  :: ne(:)
    integer,                        intent(in)  :: nv(:, :)
    real,                           intent(in)  :: eigKet(:)
    real,                           intent(in)  :: El(:, 0:, :, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    integer,                        intent(in)  :: iloTable(:, 0:, :)
    integer,                        intent(in)  :: ilst(:, :, :)
    integer,                        intent(in)  :: GbasVec(:, :)
    complex,                       intent(in)  :: zKet(:, :)
    complex,                       intent(in)  :: zBra(:, :)
    integer,                        intent(in)  :: kveclo(:,:)
    complex,                        intent(in)  :: vEff1IR(:, :)
    integer,                        intent(in)  :: nobd(:, :)
    integer,                        intent(in)  :: nlo_atom(:)
    real,                           intent(in)  :: recEdiffME(:, :)
    integer,                        intent(in)  :: kpq2kPrVec(:, :, :)
    complex,                        intent(in)  :: surfIntVFast(:, :, :)
    integer,                        intent(in)  :: gpqdp(:, :) ! todo : Deprecated, should be removed when finite q working
    real,                           intent(in)  :: eigBra(:) ! todo Deprecated: Should be removed after closing #45
    real,                           intent(in)  :: cutContr(:, :) ! todo Deprecated: Should be removed after closing #45
    complex,                        intent(out) :: z1G(:, :)
    complex, optional,              intent(in)  :: haa(:,:,0:,:,0:), dhaa(:,:,0:,:,0:)
    real,    optional,              intent(in)  :: rbas1(:, :, 0:, :, :)
    complex, optional,              intent(inout)  :: mat_elH(:,:), mat_elS(:,:)
    complex,                           intent(in)  :: loosetd1(:, :, :, :), loosetd2(:, :, :, :), loosetd(:, :, :, :)

    ! Type parameters
    type(t_noco)                                :: noco
    type(t_nococonv)                            :: nococonv
   
    TYPE(t_lapw) :: lapw
    TYPE (t_mat) :: zMatKet, zMatBra, zMatTilde, zMatBar

    ! Array parameters
    integer,           allocatable              :: ngoprI(:)
    complex,          allocatable              :: h0MTBv(:, :)
    complex,          allocatable              :: s0MTBv(:, :)
    complex,          allocatable              :: h0MTKv(:, :)
    complex,          allocatable              :: s0MTKv(:, :)
    complex,          allocatable              :: vEff1IRMat(:, :)
    complex,           allocatable              :: acofKet(:, :, :)
    complex,           allocatable              :: bcofKet(:, :, :)
    complex,           allocatable              :: ccofKet(:, :, :, :)
    complex,           allocatable              :: atestcofk2(:,:,:)
    complex,           allocatable              :: btestcofk2(:,:,:)
    complex,           allocatable              :: bascof_lo(:,:,:,:,:)
    complex,           allocatable              :: acofBra(:, :, :)
    complex,           allocatable              :: bcofBra(:, :, :)
    complex,           allocatable              :: ccofBra(:, :, :, :)
    complex,           allocatable              :: acofBar(:, :, :)
    complex,           allocatable              :: bcofBar(:, :, :)
    complex,           allocatable              :: ccofBar(:, :, :, :)
    complex                                     :: atestcofk((atoms%lmaxd+1)**2,nv(1,ikpt)), btestcofk((atoms%lmaxd+1)**2,nv(1,ikpt))
    complex                                     :: atestcofkq((atoms%lmaxd+1)**2,nv(1,ikpq)),btestcofkq((atoms%lmaxd+1)**2,nv(1,ikpq))
    complex                                     :: datestcofk((atoms%lmaxd+1)**2,nv(1,ikpt)), dbtestcofk((atoms%lmaxd+1)**2,nv(1,ikpt))
    complex                                     :: datestcofkq((atoms%lmaxd+1)**2,nv(1,ikpq)),dbtestcofkq((atoms%lmaxd+1)**2,nv(1,ikpq))
    complex                                     :: btemp(nv(1,ikpt))
    complex,           allocatable              :: acofTilde(:, :, :)
    complex,           allocatable              :: bcofTilde(:, :, :)
    complex,           allocatable              :: ccofTilde(:, :, :, :)
    complex,           allocatable              :: zBar(:, :)
    complex,           allocatable              :: zTilde(:, :)
    complex,           allocatable              :: mCoefK(:, :, :)
    complex,           allocatable              :: mCoefB(:, :, :)
    complex ,          allocatable              :: z1Band(:, :), hepss1band(:, :)
    complex ,          allocatable              :: vSumMT(:, :)
    complex,           allocatable              :: mCoefKv(:, :)
    complex,           allocatable              :: mCoefBv(:, :)
    complex,           allocatable              :: surfIntTeps(:, :)
    complex,           allocatable              :: surfInt(:, :)
    real                                        :: Gext(3)
    real                                        :: GpqCart(3)
    complex                                     :: a1(2*(atoms%lmaxd+1)**2,nv(1,ikpq)), b1(2*(atoms%lmaxd+1)**2,nv(1,ikpt))
    complex                                     :: a2(2*(atoms%lmaxd+1)**2,nv(1,ikpq)), b2(2*(atoms%lmaxd+1)**2,nv(1,ikpt))
    complex                                     :: b3(2*(atoms%lmaxd+1)**2,nv(1,ikpq))
    complex                                     :: vSumG(nv(1,ikpq),nv(1,ikpt)), vSumG2(nv(1,ikpq),nv(1,ikpt)), vSumG3(nv(1,ikpq),nv(1,ikpt))

    ! Scalar parameters
    integer                                     :: nmat
    integer                                     :: iBas
    integer                                     :: ieig
    integer                                     :: itype
    integer                                     :: iatom
    integer                                     :: ieqat
    integer                                     :: lmp
    integer                                     :: lm
    integer                                     :: oqn_l
    integer                                     :: mqn_m
    integer                                     :: p
    integer                                     :: ilo
    integer                                     :: nband
    integer                                     :: pband
    integer :: nk

    integer :: n1,n2, l, io, igp, lp, jo, lm2, k
    integer :: m, mp, ngp, ngpq, lmmax, i, iG, iGq
    logical :: didwe
    real    :: tempReal, tempImag, t2
    complex :: z1, z2, rgaunt
    complex :: hMatBand(ne(ikpq),nobd(ikpt,1)), sMatBand(ne(ikpq),nobd(ikpt,1))

    complex(8), external :: zdotc

    ! Initializiation
    z1G(:, :) = cmplx(0., 0.)

    ! Set noco stuff to false.
    noco%l_soc  = .FALSE.
    noco%l_noco = .FALSE.
    noco%l_ss   = .FALSE.
    ALLOCATE(noco%l_constrained(atoms%ntype))
    ALLOCATE(noco%l_unrestrictMT(atoms%ntype))
    ALLOCATE(noco%l_spinoffd_ldau(atoms%ntype))
    noco%l_constrained(atoms%ntype)   = .FALSE.
    noco%l_unrestrictMT(atoms%ntype)  = .FALSE.
    noco%l_spinoffd_ldau(atoms%ntype) = .FALSE.

    ! Calculate <ψ_{k'n'}^{(0)}|V^{(1)}_{eff}(q)|ψ_{kn}^{(0)}>_{IR}
    ! -------------------------------------------------------------
    allocate( vEff1IRMat(ne(ikpq), nobd(ikpt,1)) )
    vEff1IRMat(:, :) = cmplx(0., 0.)
    ! nmat has to be the dimension of the bras due to the concept of calcMEPotIR
    nmat = nv(1, ikpq) + atoms%nlotot
    if (present(mat_elH)) then
      call calcMEPotIR( stars, GbasVec(:, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(:, ilst(:nv(1, ikpt), ikpt, 1)), nv, &
        & vEff1IR(:,idir), zBra, zKet, gdp, nmat, ne(ikpq), nobd(ikpt,1), ikpt, iqpt, ikpq, ngdp, vEff1IRMat, kpq2kPrVec, idir, mat_elH)
    else
      call calcMEPotIR( stars, GbasVec(:, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(:, ilst(:nv(1, ikpt), ikpt, 1)), nv, &
        & vEff1IR(:,idir), zBra, zKet, gdp, nmat, ne(ikpq), nobd(ikpt,1), ikpt, iqpt, ikpq, ngdp, vEff1IRMat, kpq2kPrVec)
    end if

    ! Calculate <ψ_{k'n'}^{(0)}|V^{(1)}_{eff}(q)|ψ_{kn}^{(0)}>_{allMT}
    ! -------------------------------------------------------------

    ! The abcofs should not rotated in the muffin-tin but are all described by one global coordinate system. Therefore the rotation
    ! operator is just the unity operator (indexed by 1) for all atoms.
    allocate(ngoprI(atoms%nat))
    allocate(nococonv%alph(atoms%ntype), nococonv%beta(atoms%ntype))
    ngoprI(:) = 1

    ! Abcof for ψ_{kn}^{(0)}
    allocate(acofKet(nobd(ikpt,1), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bcofKet(nobd(ikpt,1), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      &ccofKet(-atoms%llod:atoms%llod, nobd(ikpt,1), atoms%nlod, atoms%nat))
    nmat = nv(1, ikpt) + atoms%nlotot
    acofKet(:, :, :) = cmplx(0., 0.)
    bcofKet(:, :, :) = cmplx(0., 0.)
    ccofKet(:, :, :, :) = cmplx(0., 0.)

    !if (.false.) then
    !  if (idir.eq.1) then
!        if (ikpt.eq.1) then
!          open(109,file='000_match',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
!        else
!          open(109,file='000_match',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
!        end if
!      end if

!      atestcofk(:,:) = cmplx(0., 0.)
!      btestcofk(:,:) = cmplx(0., 0.)
    !  allocate(atestcofk2(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), btestcofk2(MAXVAL(nv), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
    !  & bascof_lo(3,-atoms%llod:atoms%llod,4*atoms%llod+2,atoms%nlod,atoms%nat))

    !  atestcofk2(:,:,:) = cmplx(0., 0.)
    !  btestcofk2(:,:,:) = cmplx(0., 0.)
    !  bascof_lo(:,:,:,:,:) = cmplx(0.,0.)


    !  call abcof3(atoms%lmaxd,atoms%ntype,atoms%nat,sym%nop,MAXVAL(nv),input%jspins,1, &
    ! &                 atoms%lmaxd*(atoms%lmaxd+2),dimens%nbasfcn,atoms%llod,atoms%nlod,atoms%nlotot,sym%invtab, &
    ! &                 atoms%ntype,sym%mrot,ngoprI,atoms%taual,atoms%neq,atoms%lmax,atoms%rmt,cell%omtil, &
    ! &                 cell%bmat,cell%bbmat,kpts%bk(:,ikpt),GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
    ! &                 GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), &
    ! &                 nv(:, ikpt),  nmat, &
    ! &                 usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds(:, :, 1), usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, &
    ! & sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo,&
    ! & atoms%l_dulo, atoms%lapw_l, kveclo(:, ikpt),odi,ods,atestcofk2,btestcofk2,bascof_lo)

    ! if (idir.eq.1.and..false.) then
    ! open(109,file='000_match3',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
    !   do k=1, nv(1,ikpt)
    !     do l=0, atoms%lmaxd
    !       do m=-l, l
    !         lm=l*(l+1)+m
    !         write(109,*) lm+1, 1
    !         write(109,*) real(atestcofk2(k,lm,1)), aimag(atestcofk2(k,lm,1))
    !         write(109,*) lm+1, 2
    !         write(109,*) real(btestcofk2(k,lm,1)), aimag(btestcofk2(k,lm,1))
    !       end do
    !     end do
    !   end do
    !   close(109)
    ! end if

    ! deallocate(atestcofk2,btestcofk2,bascof_lo)
        !call abcof ( atoms%lmaxd, atoms%ntype, input%neig, nobd(ikpt,1), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
        !& SIZE(zKet(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
        !& atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
        !& GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt),  nmat, nobd(ikpt,1), &
        !& zKet(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, &
        !& sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo,&
        !& atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpt), odi, ods, &
        !& acofKet, bcofKet, ccofKet, atestcofk, btestcofk, idir, ikpt)

        nk=fmpi%k_list(ikpt)
        CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
        CALL zMatKet%init(.FALSE., nv(1, ikpt) + atoms%nlotot, nobd(ikpt, 1))
        zMatKet%data_c(:, :) = zKet(:nv(1, ikpt) + atoms%nlotot, :nobd(ikpt, 1))
        CALL abcof(input, atoms, sym, cell, lapw, nobd(ikpt, 1), usdus, noco, nococonv, 1,   &
                 & acofKet(:, 0:, :), bcofKet(:, 0:, :), &
                 & ccofKet(-atoms%llod:, :, :, :), zMatKet)
        !CALL save_npy("acofKet.npy", acofKet)

      !if (idir.eq.1) then
    !    close(109)
     ! end if
    !else
    !  call abcof ( atoms%lmaxd, atoms%ntype, input%neig, nobd(ikpt,1), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
    !    & SIZE(zKet(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
    !    & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
    !    & GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt),  nmat, nobd(ikpt,1), &
    !    & zKet(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, &
    !    & sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo,&
    !    & atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpt), odi, ods, &
    !    & acofKet, bcofKet, ccofKet)
    !end if

    ! Abcof for ψ_{k'n'}^{(0)}
    allocate(acofBra(ne(ikpq), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat))
    allocate(bcofBra(ne(ikpq), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat))
    allocate(ccofBra(-atoms%llod:atoms%llod, ne(ikpq), atoms%nlod, atoms%nat))

    acofBra(:, :, :) = cmplx(0., 0.)
    bcofBra(:, :, :) = cmplx(0., 0.)
    ccofBra(:, :, :, :) = cmplx(0., 0.)
    nmat = nv(1, ikpq) + atoms%nlotot

    !if (.false.) then
     ! atestcofkq(:,:) = cmplx(0., 0.)
      !btestcofkq(:,:) = cmplx(0., 0.)
      !if (idir.eq.1) then
!        if (ikpt.eq.1) then
!          open(109,file='000_matchq',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
!        else
!          open(109,file='000_matchq',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
!        end if
!      end if

!      call abcof ( atoms%lmaxd, atoms%ntype, input%neig, ne(ikpq), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
!        & SIZE(zBra(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
!        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), GbasVec(1, ilst(:nv(1, ikpq), ikpq, 1)), &
!        & GbasVec(2, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(3, ilst(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq),  nmat, ne(ikpq), zBra(:, :),&
!        & usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, sym%invsatnr, &
!        & usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1), usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
!        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpq), odi, ods, &
!        & acofBra, bcofBra, ccofBra, atestcofkq, btestcofkq, idir, ikpq)
!      if (idir.eq.1) then
!        close(109)
!      end if

!    else
!      call abcof ( atoms%lmaxd, atoms%ntype, input%neig, ne(ikpq), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
!        & SIZE(zBra(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
!        & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), GbasVec(1, ilst(:nv(1, ikpq), ikpq, 1)), &
!        & GbasVec(2, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(3, ilst(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq),  nmat, ne(ikpq), zBra(:, :),&
!        & usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, sym%invsatnr, &
!        & usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1), usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, atoms%l_dulo, &
!        & atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpq), odi, ods, &
!        & acofBra, bcofBra, ccofBra)
    nk=fmpi%k_list(ikpq)
    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpq, cell, fmpi)
    CALL zMatBra%init(.FALSE., nv(1, ikpq) + atoms%nlotot, ne(ikpq))
    zMatBra%data_c(:, :) = zBra(:nv(1, ikpq) + atoms%nlotot, :ne(ikpq))
    CALL abcof(input, atoms, sym, cell, lapw, ne(ikpq), usdus, noco, nococonv, 1,   &
            & acofBra(:, 0:, :), bcofBra(:, 0:, :), &
            & ccofBra(-atoms%llod:, :, :, :), zMatBra)
    !CALL save_npy("acofBra.npy", acofBra)
!    end if

    if (.false.) then
      ngp=nv(1, ikpt)
      ngpq=nv(1, ikpq)
      lmmax=(atoms%lmaxd+1)**2

      datestcofkq(:,:) = cmplx(0., 0.)
      dbtestcofkq(:,:) = cmplx(0., 0.)
      do iGq = 1, nv(1, ikpq)
        GpqCart(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(iGq, ikpq, 1)) + kpts%bk(1:3, ikpq))
        datestcofkq(:,iGq)=ImagUnit * GpqCart(idir) * atestcofkq(:,iGq)
        dbtestcofkq(:,iGq)=ImagUnit * GpqCart(idir) * btestcofkq(:,iGq)
      end do

      datestcofk(:,:) = cmplx(0., 0.)
      dbtestcofk(:,:) = cmplx(0., 0.)
      do iG = 1, nv(1, ikpt)
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(iG, ikpt, 1)) + kpts%bk(1:3, ikpt))
        datestcofk(:,iG)=ImagUnit * Gext(idir) * atestcofk(:,iG)
        dbtestcofk(:,iG)=ImagUnit * Gext(idir) * btestcofk(:,iG)
      end do

      if (.false.) then
      if (idir.eq.1) then
        if (ikpt.eq.1) then
          open(109,file='000_zs',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
          open(110,file='000_tlmplm0',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
          open(111,file='000_tlmplm1',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
        end if
        i=0
        do lp=0,atoms%lmaxd
          do mp=-lp,lp
            do io=1,2
              i=i+1
              b1(i,:)=0.0
              b2(i,:)=0.0
              b3(i,:)=0.0
              do l=0,atoms%lmaxd
                do m=-l,l
                  do jo=1,2
                    z1=0.0
                    z2=0.0
                    do oqn_l=0,atoms%lmaxd
                      if (mod(lp+oqn_l+l,2).eq.0) then
                        do mqn_m=-oqn_l,oqn_l
                          lm2=oqn_l*(oqn_l+1)+mqn_m+1
                          z1=z1+gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd)*dhaa(lm2,jo,l,io,lp)
                          if (mqn_m.gt.0) then
                            if (mod(mqn_m,2).eq.0) then
                              t2=(gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd)+gaunt1(lp,oqn_l,l,mp,-mqn_m,m,atoms%lmaxd))/sqrt(2.0)
                            else
                              t2=(gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd)-gaunt1(lp,oqn_l,l,mp,-mqn_m,m,atoms%lmaxd))/sqrt(2.0)
                            end if
                            rgaunt=cmplx(t2,0.0)
                          else if (mqn_m.lt.0) then
                            if (mod(mqn_m,2).eq.0) then
                              t2=(gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd)-gaunt1(lp,oqn_l,l,mp,-mqn_m,m,atoms%lmaxd))/sqrt(2.0)
                            else
                              t2=(gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd)+gaunt1(lp,oqn_l,l,mp,-mqn_m,m,atoms%lmaxd))/sqrt(2.0)
                            end if
                            rgaunt=cmplx(0.0,-t2)
                          else
                            rgaunt=cmplx(gaunt1(lp,oqn_l,l,mp,mqn_m,m,atoms%lmaxd),0.0)
                          end if
                          z2=z2+rgaunt*haa(lm2,jo,l,io,lp)
                        end do
                      end if
                    end do

                    if (ikpt.eq.1.and.io*jo.eq.1) then
                      write(109,*) l, m, lp, mp
                      write(109,*) jo, io, real(z1), aimag(z1)
                      write(110,*) l, m, lp, mp
                      write(110,*) jo, io, real(z2), aimag(z2)
                    end if
                    ! Kinetic surface term
                    if (l.eq.lp.and.m.eq.mp) then
                      if (io.eq.1) then
                        if (jo.eq.1) then
                          z2=z2+0.5*usdus%us(l,1,1)*usdus%dus(l,1,1)*atoms%rmt(1)**2
                        else
                          z2=z2+0.5*usdus%us(l,1,1)*usdus%duds(l,1,1)*atoms%rmt(1)**2
                        end if
                      else
                        if (jo.eq.1) then
                          z2=z2+0.5*usdus%uds(l,1,1)*usdus%dus(l,1,1)*atoms%rmt(1)**2
                        else
                          z2=z2+0.5*usdus%uds(l,1,1)*usdus%duds(l,1,1)*atoms%rmt(1)**2
                        end if
                      end if
                    end if
                    if (ikpt.eq.1.and.io*jo.eq.1) then
                      write(111,*) l, m, lp, mp
                      write(111,*) jo, io, real(z2), aimag(z2)
                    end if
                    lm=l*(l+1)+m+1
                    if (abs(dble(z1))+abs(aimag(z1)).gt.1.d-14) then
                      if (jo.eq.1) then
                        btemp=atestcofk(lm,1:ngp)
                        call zaxpy(ngp,z1*(ImagUnit**l),btemp,1,b1(i,1),lmmax)
                      else
                        btemp=btestcofk(lm,1:ngp)
                        call zaxpy(ngp,z1*(ImagUnit**l),btemp,1,b1(i,1),lmmax)
                      end if
                    end if
                    if (abs(dble(z2))+abs(aimag(z2)).gt.1.d-14) then
                      if (jo.eq.1) then
                        btemp=datestcofk(lm,1:ngp)
                        call zaxpy(ngp,z2*(ImagUnit**l),btemp,1,b3(i,1),lmmax)
                        btemp=atestcofk(lm,1:ngp)
                        call zaxpy(ngp,z2*(ImagUnit**l),btemp,1,b2(i,1),lmmax)
                      else
                        btemp=dbtestcofk(lm,1:ngp)
                        call zaxpy(ngp,z2*(ImagUnit**l),btemp,1,b3(i,1),lmmax)
                        btemp=btestcofk(lm,1:ngp)
                        call zaxpy(ngp,z2*(ImagUnit**l),btemp,1,b2(i,1),lmmax)
                      end if
                    end if
                  end do
                end do
              end do
              lmp=lp*(lp+1)+mp+1
              if (io.eq.1) then
                a1(i,1:ngpq)=atestcofkq(lmp,1:ngpq)*(ImagUnit**lp)
                a2(i,1:ngpq)=datestcofkq(lmp,1:ngpq)*(ImagUnit**lp)
              else
                a1(i,1:ngpq)=btestcofkq(lmp,1:ngpq)*(ImagUnit**lp)
                !a1(i,1:ngpq)=0
                a2(i,1:ngpq)=dbtestcofkq(lmp,1:ngpq)*(ImagUnit**lp)
                !a2(i,1:ngpq)=0
              end if
            end do
          end do
        end do
        if (ikpt.eq.1) then
          close(109)
          close(110)
          close(111)
        end if

        open(109,file='000_ME_V1_MT_newtest',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
        open(110,file='000_ME_H0_kGprq_MT_newtest',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
        open(111,file='000_ME_H0_kG_MT_newtest',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
        do iG=1, ngp
          do iGq=1, ngpq
            write(109,*) ikpt
            write(109,*) iGq, iG
            vSumG(iGq,iG)=zdotc(2*lmmax,a1(:,iGq),1,b1(:,iG),1)
            write(109,*) real(vSumG(iGq,iG)), aimag(vSumG(iGq,iG))
            vSumG2(iGq,iG)=zdotc(2*lmmax,a2(:,iGq),1,b2(:,iG),1)
            write(110,*) real(vSumG2(iGq,iG)), aimag(vSumG2(iGq,iG))
            vSumG3(iGq,iG)=zdotc(2*lmmax,a1(:,iGq),1,b3(:,iG),1)
            write(111,*) real(vSumG3(iGq,iG)), aimag(vSumG3(iGq,iG))
          end do
        end do
        close(109)
        close(110)
        close(111)
      end if
    end if
    end if
    ! Rearrange acofs, bcofs, ccofs ensuring an consistent handling of LO terms within the loop the structure compared to LAPWs.
    ! The lmp index is comparable to the lm index but for every lm combination it includes also an index that runs over the
    ! matching coefficients for belonging to the u, the udot and the u_LO
    allocate(mCoefK(nobd(ikpt,1), maxlmp, atoms%nat))
    allocate(mCoefB(ne(ikpq), maxlmp, atoms%nat))
    iatom  = 0
    mCoefK(:, :, :) =  cmplx(0., 0.)
    mCoefB(:, :, :) = cmplx(0., 0.)
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        lmp   =   0
        lm    =  -1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = lm + 1
            !p = 1
            lmp = lmp + 1
            mCoefK(:nobd(ikpt, 1), lmp, iatom) = acofKet(:nobd(ikpt, 1), lm, iatom)
            mCoefB(:ne(ikpq), lmp, iatom) = acofBra(:ne(ikpq), lm, iatom)
            !p = 2
            lmp = lmp + 1
            mCoefK(:nobd(ikpt, 1), lmp, iatom) = bcofKet(:nobd(ikpt, 1), lm, iatom)
            mCoefB(:ne(ikpq), lmp, iatom) = bcofBra(:ne(ikpq), lm, iatom)
            !LOs
            if (.false.) then
              do p = 3, nRadFun(oqn_l, itype)
                ilo = iloTable(p, oqn_l, itype)
                lmp = lmp + 1
                mCoefK(:nobd(ikpt, 1), lmp, iatom) = ccofKet(mqn_m, :nobd(ikpt, 1), ilo, iatom)
                mCoefB(:ne(ikpq), lmp, iatom) = ccofBra(mqn_m, :ne(ikpq), ilo, iatom)
              end do
            end if
          end do
        end do
      end do
    end do
    deallocate(acofKet, bcofKet, ccofKet, acofBra, bcofBra, ccofBra)
    ! Unite the matching coefficients with the given tlmplm integrals completing the calculation of the current matrix element.
    ! The sum over all muffin-tins takes place within CalcVsumMT
    allocate( vSumMT(ne(ikpq), nobd(ikpt, 1) ) )
    vSumMT(:, :) = cmplx(0., 0.)

    if (.false.) then
      inquire(file='000_tlmplm1x',exist=didwe)

      if (.not.didwe) then
        open(111,file='000_tlmplm1x',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
      end if

      if (present(mat_elH)) then
      call calcVsumMT( atoms, td4V, td4V2, loosetd1, loosetd2, ikpt, ikpq, ne, nobd, mCoefB, mCoefK, nRadFun, iloTable, nlo_atom, vSumMT, &
                       idir, nv(1, ikpt), nv(1, ikpq) , atestcofk, btestcofk, atestcofkq, btestcofkq, mat_elH)
      else
      call calcVsumMT( atoms, td4V, td4V2, loosetd1, loosetd2, ikpt, ikpq, ne, nobd, mCoefB, mCoefK, nRadFun, iloTable, nlo_atom, vSumMT, &
                       idir, nv(1, ikpt), nv(1, ikpq) , atestcofk, btestcofk, atestcofkq, btestcofkq)
      end if
      if (.not.didwe) then
        close(111)
      end if
    else
      call calcVsumMT( atoms, td4V, td4V2, loosetd1, loosetd2, ikpt, ikpq, ne, nobd, mCoefB, mCoefK, nRadFun, iloTable, nlo_atom, vSumMT )
    end if

    ! Calculate <\tilde{ψ}_{k'n'}^{(0)}|H|ψ_{kn}^{(0)}>_{MT} and <ψ_{k'n'}^{(0)}|H|\tilde{ψ}_{kn}^{(0)}>_{MT} and additionally the
    ! overlaps <\tilde{ψ}_{k'n'}^{(0)}|ψ_{kn}^{(0)}>_{MT} and <ψ_{k'n'}^{(0)}|\tilde{ψ}_{kn}^{(0)}>_{MT}
    ! ----------------------------------------------------------------------------------------------------------------------------

    ! In order to set up the matching coefficients for the bra tilde wave function the wave function coefficients are decorated with
    ! a factor for every G.
    allocate(zTilde(SIZE(zBra(:,1)), ne(ikpq)))
    zTilde(:, :) = cmplx(0., 0.)
    do ieig = 1, ne(ikpq)
      do iBas = 1, nv(1, ikpq)
        GpqCart(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(iBas, ikpq, 1)) + kpts%bk(1:3, ikpq))
        zTilde(iBas, ieig) = ImagUnit * GpqCart(idir) * zBra(iBas, ieig)
      end do
      if (.false.) then
        ! todo LO: LOs mit berücksichtigen, auch mit Gs?, + or - q?
        do iBas = nv(1, ikpq) + 1, nv(1, ikpq) + atoms%nlotot
          GpqCart(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(kveclo(iBas - nv(1, ikpq), ikpq), ikpq, 1)) &
                                                                                                             & + kpts%bk(1:3, ikpq))
          zTilde(iBas, ieig) = ImagUnit * Gpqcart(idir) * zBra(iBas, ieig)
        end do
      end if
    end do

    ! Calculate matching coefficients for \tilde{ψ}_{k'n'}^{(0)}
    allocate(acofTilde(ne(ikpq), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bcofTilde(ne(ikpq), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      &ccofTilde(-atoms%llod:atoms%llod, ne(ikpq), atoms%nlod, atoms%nat))

    acofTilde(:, :, :) = cmplx(0., 0.)
    bcofTilde(:, :, :) = cmplx(0., 0.)
    ccofTilde(:, :, :, :) = cmplx(0., 0.)
    nmat = nv(1, ikpq) + atoms%nlotot

    !call abcof ( atoms%lmaxd, atoms%ntype, ne(ikpq), ne(ikpq), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
     ! & SIZE(zTilde(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
     ! & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpq), GbasVec(1, ilst(:nv(1, ikpq), ikpq, 1)), &
     ! & GbasVec(2, ilst(:nv(1, ikpq), ikpq, 1)), GbasVec(3, ilst(:nv(1, ikpq), ikpq, 1)), nv(:, ikpq),  nmat, ne(ikpq), &
     ! & zTilde(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, &
     ! & sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1), usdus%dulos(:, :, 1), atoms%llo, atoms%nlo, &
     ! & atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpq), odi, ods, &
     ! & acofTilde, bcofTilde, ccofTilde)

      nk=fmpi%k_list(ikpq)
      CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpq, cell, fmpi)
      CALL zMatTilde%init(.FALSE., nv(1, ikpq) + atoms%nlotot, ne(ikpq))
      zMatTilde%data_c(:, :) = zTilde(:nv(1, ikpq) + atoms%nlotot, :ne(ikpq))
      CALL abcof(input, atoms, sym, cell, lapw, ne(ikpq), usdus, noco, nococonv, 1,   &
               & acofTilde(:, 0:, :), bcofTilde(:, 0:, :), &
               & ccofTilde(-atoms%llod:, :, :, :), zMatTilde)
      !CALL save_npy("acofTilde.npy", acofTilde)
    deallocate(zTilde)

    ! In order to set up the matching coefficients for the ket tilde wave function the wave function coefficients are decorated with
    ! a factor for every G.
    allocate(zBar(SIZE(zKet(:,1)), nobd(ikpt,1)))
    zBar(:, :) = cmplx(0., 0.)
    ! is that correct that it oes to ne(ikpt and not ne(ikpq)
    do ieig = 1, nobd(ikpt, 1)
      do iBas = 1, nv(1, ikpt)
        Gext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(iBas, ikpt, 1)) + kpts%bk(1:3, ikpt))
        zBar(iBas, ieig) = ImagUnit * Gext(idir) * zKet(iBas, ieig)
      end do
      if (.false.) then
        do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
          Gext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) &
                                                                                                             & + kpts%bk(1:3, ikpt))
          zBar(iBas, ieig) = ImagUnit * Gext(idir) * zKet(iBas, ieig)
        end do
      end if
    end do

    ! Calculate matching coefficients for \tilde{ψ}_{kn}^{(0)}
    allocate(acofBar(nobd(ikpt,1), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), bcofBar(nobd(ikpt,1), 0:atoms%lmaxd*(atoms%lmaxd+2), atoms%nat), &
      &ccofBar(-atoms%llod:atoms%llod, nobd(ikpt,1), atoms%nlod, atoms%nat))
    acofBar(:, :, :) = cmplx(0., 0.)
    bcofBar(:, :, :) = cmplx(0., 0.)
    ccofBar(:, :, :, :) = cmplx(0., 0.)

    nmat = nv(1, ikpt) + atoms%nlotot
    !call abcof ( atoms%lmaxd, atoms%ntype, nobd(ikpt,1), nobd(ikpt,1), atoms%nat, sym%nop, MAXVAL(nv), input%jspins, atoms%lmaxd*(atoms%lmaxd+2), &
    !  & SIZE(zBar(:,1)), atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, &
    !  & atoms%lmax, atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
    !  & GbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), GbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt),  nmat, nobd(ikpt,1), &
    !  & zBar(:, :), usdus%us(:, :, 1), usdus%dus(:, :, 1), usdus%uds, usdus%duds(:, :, 1), usdus%ddn(:, :, 1), sym%invsat, &
    !  & sym%invsatnr, usdus%ulos(:, :, 1), usdus%uulon(:, :, 1), usdus%dulon(:, :, 1),  usdus%dulos(:, :, 1), atoms%llo, atoms%nlo,&
    !  & atoms%l_dulo, atoms%lapw_l, noco%l_noco, noco%l_ss, 1, nococonv%alph, nococonv%beta, nococonv%qss, kveclo(:, ikpt), odi, ods, &
    !  & acofBar, bcofBar, ccofBar)
    nk=fmpi%k_list(ikpt)
    CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell, fmpi)
    CALL zMatBar%init(.FALSE., nv(1, ikpt) + atoms%nlotot, nobd(ikpt, 1))
    zMatBar%data_c(:, :) = zBar(:nv(1, ikpt) + atoms%nlotot, :nobd(ikpt, 1))
    CALL abcof(input, atoms, sym, cell, lapw, nobd(ikpt, 1), usdus, noco, nococonv, 1,   &
             & acofBar(:, 0:, :), bcofBar(:, 0:, :), &
             & ccofBar(-atoms%llod:, :, :, :), zMatBar)
    !CALL save_npy("acofBar.npy", acofBar)
    !STOP
    deallocate (zBar)

    ! Rearrange Hamiltonian acofs, bcofs, ccofs ensuring a consistent handling of LOs within the loop structure compared to LAPWs
    allocate( mCoefKv(nobd(ikpt,1), maxlmp), mCoefBv(ne(ikpq), maxlmp) )
    mCoefKv(:, :) = cmplx(0., 0.)
    mCoefBv(:, :) = cmplx(0., 0.)
    lmp = 0
    lm = -1
    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        lm = lm + 1
        !p = 1
        lmp = lmp + 1
        mCoefKv(:nobd(ikpt,1), lmp) = acofBar(:nobd(ikpt,1), lm, iDatom)
        mCoefBv(:ne(ikpq), lmp) = acofTilde(:ne(ikpq), lm, iDatom)
        !p = 2
        lmp = lmp + 1
        mCoefKv(:nobd(ikpt,1), lmp) = bcofBar(:nobd(ikpt,1), lm, iDatom)
        mCoefBv(:ne(ikpq), lmp) = bcofTilde(:ne(ikpq), lm, iDatom)
        !LOs
        if (.false.) then
          do p = 3, nRadFun(oqn_l, iDtype)
            ilo = iloTable(p, oqn_l, iDtype)
            lmp = lmp + 1
            mCoefKv(:nobd(ikpt,1), lmp) = ccofBar(mqn_m, :nobd(ikpt,1), ilo, iDatom)
            mCoefBv(:ne(ikpq), lmp) = ccofTilde(mqn_m, :ne(ikpq), ilo, iDatom)
          end do
        end if
      end do
    end do
    deallocate(acofTilde, bcofTilde, ccofTilde, acofBar, bcofBar, ccofBar)

    allocate( h0MTBv(ne(ikpq), nobd(ikpt,1)), s0MTBv(ne(ikpq), nobd(ikpt,1)) )
    allocate( h0MTKv(ne(ikpq), nobd(ikpt,1)), s0MTKv(ne(ikpq), nobd(ikpt,1)) )
    h0MTBv(:, :) = cmplx(0., 0.)
    s0MTBv(:, :) = cmplx(0., 0.)
    h0MTKv(:, :) = cmplx(0., 0.)
    s0MTKv(:, :) = cmplx(0., 0.)
    ! Calculate <\tilde{ψ}_{k'n'}^{(0)}|H|ψ_{kn}^{(0)}>_{MT} given the tlmplm integrals of the non-spherical potential as well as
    ! the overlap <\tilde{ψ}_{k'n'}^{(0)}|ψ_{kn}^{(0)}>_{MT}
    if (.false.) then
      inquire(file='000_tlmplm0',exist=didwe)

      if (.not.didwe) then
        open(111,file='000_tlmplm0',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
      end if

      if (present(mat_elH)) then
      call calcHS0MT( atoms, usdus, td4HS0, loosetd, ikpt, ikpq, iDtype, iDatom, ne, nobd, El, mCoefBv, mCoefK(:, :, iDatom), nRadFun, &
        & iloTable, nlo_atom, s0MTBv, h0MTBv, &
        & idir,nv(1,ikpt),nv(1,ikpq),atestcofk,btestcofk,atestcofkq,btestcofkq,datestcofk,dbtestcofk,datestcofkq,dbtestcofkq,&
        & mat_elH, mat_elS)
      else
      call calcHS0MT( atoms, usdus, td4HS0, loosetd, ikpt, ikpq, iDtype, iDatom, ne, nobd, El, mCoefBv, mCoefK(:, :, iDatom), nRadFun, &
        & iloTable, nlo_atom, s0MTBv, h0MTBv, &
        & idir,nv(1,ikpt),nv(1,ikpq),atestcofk,btestcofk,atestcofkq,btestcofkq,datestcofk,dbtestcofk,datestcofkq,dbtestcofkq)
      end if

      if (.not.didwe) then
        close(111)
      end if
    else
      call calcHS0MT( atoms, usdus, td4HS0, loosetd, ikpt, ikpq, iDtype, iDatom, ne, nobd, El, mCoefBv, mCoefK(:, :, iDatom), nRadFun, &
        & iloTable, nlo_atom, s0MTBv, h0MTBv )
    end if

    ! Calculate <ψ_{k'n'}^{(0)}|H|\tilde{ψ}_{kn}^{(0)}>_{MT} given the tlmplm integrals of the non-spherical potential as well as
    ! the overlap <ψ_{k'n'}^{(0)}|\tilde{ψ}_{kn}^{(0)}>_{MT}
    call calcHS0MT( atoms, usdus, td4HS0, loosetd, ikpt, ikpq, iDtype, iDatom, ne, nobd, El, mCoefB(:, :, iDatom), mCoefKv, nRadFun, &
      & iloTable, nlo_atom, s0MTKv, h0MTKv )
    deallocate( mCoefBv, mCoefKv )
    deallocate( mCoefK, mCoefB )


    ! Calculate the surface integrals ∮dS ψ_{k'n'}^{*IR(0)} T ψ_{k'n'}^{IR(0)} and ∮dS ψ_{k'n'}^{*IR(0)} ψ_{k'n'}^{IR(0)}
    allocate( surfIntTeps(ne(ikpq), nobd(ikpt,1)) )
    allocate( surfInt(ne(ikpq), nobd(ikpt,1)) )
    surfIntTeps(:, :) = cmplx(0., 0.)
    surfInt(:, :) = cmplx(0., 0.)

    if (present(mat_elH)) then
    call calcSintKinEps(atoms, cell, kpts, qpts, iDtype, iDatom, nobd(ikpt,1), ne(ikpq), ikpt, ikpq, iqpt, idir, nv(1, :), GbasVec,&
      & ilst, zBra, zKet, surfIntTeps, surfInt, kpq2kPrVec, mat_elH, mat_elS)
    else
    call calcSintKinEps(atoms, cell, kpts, qpts, iDtype, iDatom, nobd(ikpt,1), ne(ikpq), ikpt, ikpq, iqpt, idir, nv(1, :), GbasVec,&
      & ilst, zBra, zKet, surfIntTeps, surfInt, kpq2kPrVec)
    end if

    if (present(mat_elH)) then
      hMatBand(:,:) = cmplx(0.0,0.0)
      sMatBand(:,:) = cmplx(0.0,0.0)
      do iGq=1,nv(1,ikpq)
        do iG=1,nv(1,ikpt)
          do nBand = 1, nobd(ikpt, 1)
          hMatBand(:ne(ikpq), nBand) = hMatBand(:ne(ikpq), nBand) + &
                               & conjg(zBra(iGq,:ne(ikpq)))*mat_elH(iGq,iG)*zKet(iG, nBand)
          sMatBand(:ne(ikpq), nBand) = sMatBand(:ne(ikpq), nBand) + &
                               & conjg(zBra(iGq,:ne(ikpq)))*mat_elS(iGq,iG)*zKet(iG, nBand)
          end do
        end do
      end do
    end if

    ! Sum up right side multiplied with inverted left side in Sternheimer equation for given direction and k-point, k+q-point.
    ! As we do not support polyatomic metals currently, epsilon1 is not added here, but tested to be zero in a seperate test.
    allocate( z1Band( ne(ikpq), nobd(ikpt,1) ) )
    ALLOCATE(hepss1band(ne(ikpq),nobd(ikpt,1)))
    z1Band(:, :) = cmplx(0., 0.)
    hepss1band(:, :) = cmplx(0., 0.)
    do nBand = 1, nobd(ikpt, 1)
      do pBand = 1, ne(ikpq) ! pBand >> nBand
         hepss1band(pBand, nBand) = 1.0 * ( &
                                  &   killcont(1) * Veff1IRMat(pBand, nBand) &
                                  & + killcont(2) * vSumMT(pBand, nBand) &
                                  & + killcont(3) * h0MTBv(pBand, nBand) &
                                  & + killcont(4) * h0MTKv(pBand, nBand) &
                                  & + killcont(5) * surfIntTeps(pBand, nBand) &
                                  & + killcont(6) * surfIntVFast(pBand, nBand, idir) &
                                  & - killcont(7) * eigKet(nBand) * s0MTBv(pBand, nBand) &
                                  & - killcont(8) * eigKet(nBand) * s0MTKv(pBand, nBand) &
                                  & - killcont(9) * eigKet(nBand) * surfInt(pBand, nBand) &
                                  & )
         z1Band(pBand, nBand) = -recEdiffME(pBand, nBand) * ( &
                              &   killcont(1) * Veff1IRMat(pBand, nBand) &
                              & + killcont(2) * vSumMT(pBand, nBand) &
                              & + killcont(3) * h0MTBv(pBand, nBand) &
                              & + killcont(4) * h0MTKv(pBand, nBand) &
                              & + killcont(5) * surfIntTeps(pBand, nBand) &
                              & + killcont(6) * surfIntVFast(pBand, nBand, idir) &
                              & - killcont(7) * eigKet(nBand) * s0MTBv(pBand, nBand) &
                              & - killcont(8) * eigKet(nBand) * s0MTKv(pBand, nBand) &
                              & - killcont(9) * eigKet(nBand) * surfInt(pBand, nBand) &
                              & )

        !if (present(mat_elH).and..false.) then
        !  z1Band(pBand, nBand) = -recEdiffME(pBand, nBand) * &
        !                       & (hMatBand(pBand,nBand) - eigKet(nBand)*sMatBand(pBand,nBand))
        !end if
      end do
    end do

    if (iqpt == 1) then ! q is zero, Goldstone mode.
      do nBand = 1, nobd(ikpt, 1)
        if ( abs(z1Band(nBand, nBand)) > 1e-7 ) then
          call juDFT_warn( 'Variation of wave function expansion coefficient do not cancel for occupied occupied combinations!', &
            &calledby='solveSternheimerEq', hint='Debug.' )
        end if
      end do ! nBand
    end if

    ! Transformation from basis-function to wave-function space
    z1G(:nv(1, ikpq) + atoms%nlotot, :nobd(ikpt,1)) = matmul( zBra(:nv(1, ikpq) + atoms%nlotot, :ne(ikpq)), &
                                                                                                & z1Band(:ne(ikpq), :nobd(ikpt,1)) )

    !if (.false.) then
!      if (ikpt.eq.1.and.idir.eq.1) then
!        open(109,file='000_z1',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
!        open(110,file='000_z00',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
!      else
!        open(109,file='000_z1',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
!        if (idir.eq.1) then
!          open(110,file='000_z00',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
!        end if
!      end if
!      do n2=1, nobd(ikpt,1)
!        do n1=1, nv(1, ikpq) + atoms%nlotot
!          Gext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(n1, ikpq, 1)) + kpts%bk(1:3, ikpq))
!          write(109,*) ikpq, n1, n2,  idir
!          write(109,*) Gext(1), Gext(2), Gext(3)
!          write(109,*) real(z1G(n1,n2)), aimag(z1G(n1,n2))
!        end do
!        if (idir.eq.1) then
!          do n1=1, nv(1, ikpt) + atoms%nlotot
!            Gext(1:3) = matmul( cell%bmat(1:3, 1:3), GbasVec(1:3, ilst(n1, ikpt, 1)) + kpts%bk(1:3, ikpt))
!            write(110,*) ikpt, n1, n2
!            write(110,*) Gext(1), Gext(2), Gext(3)
!            write(110,*) real(zKet(n1,n2)), aimag(zKet(n1,n2))
!          end do
!        end if
!      end do
!      close(109)
!      if (idir.eq.1) then
!        close(110)
!      end if
    !end if

    !if (.FALSE.) then
      deallocate(z1Band)
    !end if

  end subroutine solveSternheimerEq

  subroutine genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )

    use m_types

    implicit none

    ! Type parameters
    type(t_stars),          intent(in)   :: stars
    type(t_cell),           intent(in)   :: cell
    type(t_input),          intent(in)   :: input

    ! Scalar parameters
    integer,                intent(out)  :: ngpqdp
    integer,                intent(out)  :: ngpqdp2km

    ! Array parameters
    real,                   intent(in)   :: qpoint(:)
    integer,  allocatable,  intent(out)  :: gpqdp(:, :)
    integer,  allocatable,  intent(out)  :: gpqdp2Ind(:, :, :)
    integer,                intent(out)  :: gpqdp2iLim(2, 3)

    ! Scalar variables
    integer                              :: ngrest
    integer                              :: iGx
    integer                              :: iGy
    integer                              :: iGz
    integer                              :: iG

    ! Array variables
    integer,  allocatable                :: gpqdptemp2kmax(:, :)
    integer,  allocatable                :: gpqdptemprest(:, :)
    integer                              :: Gint(3)
    real                                 :: Gpqext(3)

    allocate( gpqdptemp2kmax(3, (2 * stars%mx1 + 1) * (2 * stars%mx2 + 1) * (2 * stars%mx3 +  1)), &
            & gpqdptemprest(3, (2 * stars%mx1 + 1) * (2 * stars%mx2 + 1) * (2 * stars%mx3 +  1)) )

    ngpqdp = 0
    ngpqdp2km = 0
    ngrest = 0
    gpqdptemp2kmax(:, :) = cmplx(0., 0.)
    gpqdptemprest(:, :) = cmplx(0., 0.)
    ! From all possible G-vectors in a box, only these are accepted which are element of a sphere with radius gmax which is shifted.
    ! We need a little bit more than k*d because they are thought for a Gmax ball that is not shifted by a q, i.e. |G+q|<Gmax
    do iGx = -(stars%mx1 + 3), (stars%mx1 + 3)
      do iGy = -(stars%mx2 + 3), (stars%mx2 + 3)
        do iGz = -(stars%mx3 + 3), (stars%mx3 + 3)
          Gint = [iGx, iGy, iGz]
          Gpqext =  matmul(cell%bmat, real(Gint(1:3) + qpoint(1:3))) !transform from internal to external coordinates
          if (norm2(Gpqext) <= input%gmax) then
            ngpqdp = ngpqdp + 1
            ! Sort G-vectors
            if ( norm2(Gpqext) <= 2 * input%rkmax ) then
              ngpqdp2km = ngpqdp2km + 1
              gpqdptemp2kmax(1:3, ngpqdp2km) = Gint(1:3)
            else
              ngrest = ngrest + 1
              gpqdptemprest(1:3, ngrest) = Gint(1:3)
            end if
          end if
        end do !iGz
      end do !iGy
    end do !iGx
    allocate(gpqdp(3, ngpqdp))
    ! Mapping array from G-vector to G-vector index
    gpqdp(:, :) = 0
    gpqdp(1:3, 1:ngpqdp2km) = gpqdptemp2kmax(1:3, 1:ngpqdp2km)
    gpqdp(1:3, ngpqdp2km + 1 : ngpqdp) = gpqdptemprest(1:3, 1:ngrest)

    gpqdp2iLim(1, 1) = minval(gpqdp(1, :))
    gpqdp2iLim(2, 1) = maxval(gpqdp(1, :))
    gpqdp2iLim(1, 2) = minval(gpqdp(2, :))
    gpqdp2iLim(2, 2) = maxval(gpqdp(2, :))
    gpqdp2iLim(1, 3) = minval(gpqdp(3, :))
    gpqdp2iLim(2, 3) = maxval(gpqdp(3, :))

    allocate(gpqdp2Ind(gpqdp2iLim( 1, 1) : gpqdp2iLim( 2, 1), gpqdp2iLim( 1, 2) : gpqdp2iLim( 2, 2), gpqdp2iLim( 1, 3) : &
                                                                                                               & gpqdp2iLim( 2, 3)))
    gpqdp2Ind = 0
    do iG = 1, ngpqdp
      gpqdp2Ind(gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG)) = iG
    end do

  end subroutine genPertPotDensGvecs

  subroutine checkjuPhDens1(atomsT, cellT, ngdp, fIR, fMufT, qpt, gdp, noPts2chCont, logFileUnit)

      use m_types
      use m_ylm_old
      USE m_sphpts

      implicit none

      type(t_atoms),        intent(in) :: atomsT
      type(t_cell),         intent(in) :: cellT
      integer,              intent(in) :: noPts2chCont !number of points to check continuity
      integer,              intent(in) :: ngdp
      integer,              intent(in) :: gdp(:, :)
      integer,              intent(in)            :: logFileUnit ! why not as intent(in)?
      complex,             intent(in) :: fIR(:, :)
      complex,             intent(in) :: fMufT(:, :, :, :)
      real, intent(in) :: qpt(:)
  !    real,                 intent(in)  :: rmt(:) !MT-radii of atom type itype
  !    integer,              intent(in)  :: ntype ! number of atom types
  !    integer,              intent(in)  :: nop !number of symmetry operations
  !    integer,              intent(in)  :: nq3 ! number of stars
  !    logical,              intent(in)  :: symor !are symmetry operations symorphic or noct
  !    integer,              intent(in)  :: kv3(:, :) !reciprocal g-vector of star
  !    integer,              intent(in)  :: mrot(:, :, :) !symmetry operations, rotation matrices
  !    real,              intent(in)  :: tau(:, :) !symmetry operations, translation vectors
  !    integer,              intent(in)  :: invtab(:) !lists which symmetry operations are inverse to symmetry operations with indesx isym
  !    real,                 intent(in)  :: bmat(:, :) !reciprocal Bravais matrix transposed
  !    real,                 intent(in)  :: amat(:, :) !Bravais matrix
  !    real,                 intent(in)  :: pos(:, :) !Postions of atoms in Brillouin zone
  !    integer,              intent(in)  :: nmem(:, :) !number of lattice harmonics
  !    integer,              intent(in)  :: llh(:, :) !ls of  lattice harmonics
  !    complex,              intent(in)  :: clnu(:, :, :) !phasefactors in linear combination of lattice harmonics
  !    logical,              intent(in)  :: density ! is function density?
  !    complex,              intent(in)  :: fpw(:, :)
  !    real,                 intent(in)  :: fr(:, :, :, :)
  !    integer,              intent(in)  :: nstr(:)
  !    integer,              intent(in)  :: nat
  !    integer,              intent(in)  :: ntypsy(:)
  !    integer,              intent(in)  :: ngopr(:)
  !    integer,              intent(in)  :: lmax(:)
  !    integer,              intent(in)  :: nlh(:)
  !    integer,              intent(in)  :: mlh(:, :, :)
  !    integer,              intent(in)  :: jri(:)
  !    integer,              intent(in)  :: lmaxd
  !    integer,              intent(in)  :: neq(:)
  !
  !
      integer                           :: oqn_l
      integer                           :: mqn_m
      integer                           :: lm_lcontr
      integer                           :: lm
      integer                           :: irandPt
      integer                           :: ieq !loop variable
      real                              :: euclidNorm !variable to store eculidian L2 norm of vector
      real,   allocatable               :: randPtsCart(:, :, :) !cartesian dimension, number of Pts, itype, random points
      real,   allocatable               :: randPtsGrid(:, :, :) !random points on grid, coordinates, number of points
      real,   allocatable               :: randPtsMTLoc(:, :, :) !random points in local atomic coordinate (coordinates, number of points, ntype)
      integer                           :: itype !certain atom type
  !    complex                           :: starAtRanPt(nq3) !value of star at certain point ! ng3 = nq3
      complex,   allocatable               :: sumfuncvalI(:, :, :) !sum of function values coming from interstitial
      integer                         :: iGvec
      integer                          :: idirec
  !    real                              :: scaling !scaling factor which makes difference between density and potential or wavefunction
  !    integer                           :: icoord !loop variable
      integer                           :: imatmulc, imatmulr !loop variable
      real                              :: rvec(3) !randoom point
      real                              :: rvec_int(3)
      real  :: rotatedVector(3) !rotated vector after symmetry operation
  !    integer                           :: istar !loop variable
      integer                           :: iatom
  !    integer                           :: symAt_temp
      integer                           :: symOpr_temp
      complex                           :: ylm((atomsT%lmaxd + 1)**2) !TODO eigentlich lmax(itype)
      complex,    allocatable              :: sumOfMTFuncVal(:, :,:)
  !    integer                           :: ilath
  !    real                              :: linCombBas
  !    integer                           :: l_temp
  !    integer                           :: imem
  !    integer                           :: lm_temp
      real                              :: av(3), dms(3), rms(3)
  !    integer                           :: ispin

      allocate(randPtsCart(3, noPts2chCont, atomsT%nat), randPtsGrid(3, noPts2chCont, atomsT%ntype), randPtsMTLoc(3, noPts2chCont, atomsT%ntype) )

      allocate(sumfuncvalI(3, noPts2chCont,atomsT%nat), sumOfMTFuncVal(3, noPts2chCont,atomsT%nat))
     !Generate noPts2chCont random points on the MT sphere
      randPtsCart(:, :, :) = cmplx(0., 0.)
      iatom = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          iatom = iatom + 1
          call sphpts(randPtsCart(:, :, iatom), noPts2chCont, atomsT%rmt(itype), atomsT%pos(:, iatom))
  !      do irandPt = 1, noPts2chCont
  !        call random_number(randPtsCart(:, irandPt, itype))
  !        euclidNorm = norm2(randPtsCart(:, irandPt, itype))
  !        randPtsCart(:, irandPt, itype) = randPtsCart(:, irandPt, itype) / euclidNorm * atomsT%rmt(itype)
        end do
      end do


      !Evaluate Stars for this random points at each MT sphere boundary in the unit cell
      sumfuncValI = cmplx(0.0) !todo this is probably unneccessary, because it is already set to 0 in the loop
      iatom       = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          iatom = iatom + 1
          do irandPt = 1, noPts2chCont
  !          rvec = atomsT%pos(:,iatom) + randPtsCart(:,irandPt,itype)
            rvec(:) = randPtsCart(:, irandPt, iatom)
            ! transform rvec in internal, real-space coordinates --> rvec_int
            !call cotra1(rvec, rvec_int, cellT%bmat) !todo do this inline from now on!!!
            rvec_int = matmul(cellT%bmat, rvec) / tpi_const

            ! Evaluate stars at random point
            sumfuncValI(:, irandPt,iatom) = cmplx(0., 0.) !todo somehow this does not make sense
            do iGvec = 1, ngdp !nq3 is number of stars
              do idirec = 1, 3
                sumfuncValI(idirec, irandPt,iatom) = sumfuncvalI(idirec, irandPt,iatom) + fIR(iGvec, idirec) * exp(ImagUnit * tpi_const * dot_product(gdp(:, iGvec) + qpt(:), rvec_int(:))) !jsp noch einfügen, warum Multiplikation with G-vectors
              end do ! idirec
            end do  !iGvec
          end do  !irandPt
        end do  !ieq
      end do  !itype

      if (all(qpt(:) < 1e-6)) then
        if (any(aimag(sumfuncValI(:, :, :)) > 1e-8)) then
          write(*, *) "Imaginary part of observable potential should be 0"
          !NOstopNO
        end if
      end if


      iatom = 0
      sumOfMTFuncVal = cmplx(0., 0.)
      do itype = 1, atomsT%ntype
        !Evaluate MT sphere function at same random points
        randPtsGrid = 0.0
        do ieq = 1, atomsT%neq(itype)
          iatom       = iatom + 1
          !symOpr_temp = 1!atomsT%ngopr(iatom) ! symmetry operation mapping local to global coordinate system
          do irandPt = 1, noPts2chCont
            !rvec = randPtsCart(:,irandPt,itype)
            ! We have to subtract the atom position to get into a local coordinate system of atom iatom
            rvec = randPtsCart(:,irandPt,iatom) - atomsT%pos(:, iatom)
           ! if (symOpr_temp /= 1) then
           !   !Transform into internal real space coordinates
           !   call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
           !   do imatmulr = 1, 3 !todo optimize this
           !     rotatedVector(imatmulr) = 0.0
           !     do imatmulc = 1, 3
           !       rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) *rvec_int(imatmulc)
           !     end do
           !   end do
           !   call cotra0(rotatedVector, rvec, cellT%amat)
           ! end if
  !          call ylmnorm_init(atomsT%lmaxd + 1)
            !call ylm4(atomsT%lmax(itype) + 1, rvec, ylm)
            call ylm4(atomsT%lmax(itype) , rvec, ylm)
  !          call ylmnorm_init(atomsT%lmaxd)

            do oqn_l = 0, atomsT%lmax(itype)! + 1
              lm_lcontr = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_lcontr + mqn_m
                do idirec = 1, 3
                  sumOfMTFuncVal(idirec, irandPt, iatom) = sumOfMTFuncVal(idirec, irandPt, iatom) + fMufT(atomsT%jri(itype), lm, iatom, idirec) * ylm(lm)
                end do ! idir
              end do
            end do
          end do
        end do
      end do

      if (all(qpt(:) < 1e-6)) then
        if (any(aimag(sumOfMTFuncVal(:, :, :)) > 1e-8)) then
          write(*, *) "Imaginary part of observable potential should be 0"
          !NOstopNO
        end if
      end if

      !write ( logFileUnit, '(a)' ) "Continuity of the unperturbed potential's gradient:"
      !write(logFileUnit,*)         '---------------------------------------------------'
      iatom = 0
      do itype = 1,atomsT%ntype
        do ieq = 1,atomsT%neq(itype)
          iatom = iatom + 1
          write(logFileUnit,*) '  Atom:',iatom
          write(logFileUnit,'(a)') '  coordinate(cartesian)      IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
          do irandPt = 1,noPts2chCont
              write(logFileUnit,'(3f8.4,3x,3(4f12.8, 3x))') randPtsCart(:,irandPt,iatom),sumfuncValI(1, irandPt,iatom),sumOfMTFuncVal(1, irandPt,iatom), sumfuncValI(2, irandPt,iatom),sumOfMTFuncVal(2, irandPt,iatom), sumfuncValI(3, irandPt,iatom),sumOfMTFuncVal(3, irandPt,iatom)
          end do
          do idirec = 1, 3
            call fitchkNeg(real(sumfuncvalI(idirec, :,iatom)), real(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
          end do
          write(logFileUnit, *)
        !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
          write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'real average absolute value interstitial :',av(1), av(2), av(3)
          write(logFileUnit,'(2x,A,6f15.8,A/)')               'real rms(x, y, z), real dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
          av(:) = 0.
          rms(:) = 0.
          dms(:) = 0.
          do idirec = 1, 3
            call fitchkNeg(aimag(sumfuncvalI(idirec, :,iatom)), aimag(sumOfMTFuncVal(idirec, :,iatom)), noPts2chCont, av(idirec), rms(idirec), dms(idirec))
          end do
          write(logFileUnit, *)
        !  write(logFileUnit,'(a)') '  average deviation          IR x-comp.  MT x-comp.     IR y-comp.  MT y-comp.     IR z-comp.  MT z-comp.'
          write(logFileUnit,'(2x, A,11x, 3(f15.8, 15x))'   )  'aimag average absolute value interstitial :',av(1), av(2), av(3)
          write(logFileUnit,'(2x,A,6f15.8,A/)')               'aimag rms(x, y, z), dms(x, y, z):          ',rms(1), rms(2), rms(3), dms(1), dms(2), dms(3),' in %'
        end do
      end do
    end subroutine checkjuPhDens1

    subroutine fitchkNeg(f1, f2, n, av, rms, dmx)

      implicit none

      integer, intent(in) :: n

      real,    intent(in) :: f1(n)
      real,    intent(in) :: f2(n)

      real,    intent(out) :: av
      real,    intent(out) :: dmx
      real,    intent(out) :: rms


      real                 :: d
      integer              :: i

      av = 0.
      rms = 0.
      dmx = 0.
      do i = 1,n
         av = av + abs(f1(i))
         d = (f1(i)-f2(i))**2
         dmx = max(d,dmx)
         rms = rms + d
      end do
      av = av/n
      if (abs(av).LT.1.e-30) then
         rms = 0.
         dmx = 0.
         return
      end if
      rms = sqrt(rms/n)/av*100.
      dmx = sqrt(dmx)/av*100.

    end subroutine fitchkNeg

end module m_jpSternheimer
