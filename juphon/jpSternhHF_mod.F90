!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Helping module containing all routines for Hellmann-Feynman contributions to the Sternheimer SCC.
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> Contains routines implementing the Hellmann-Feynman contributions to the Sternheimer equation and is submodule of
!> m_jpsternheimer
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpDens1stVar.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpSternhHF

#include "cppmacro.h"

  implicit none

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Sets up interstitial part of Optimized Hellmann-Feynman Contribution to the Sternheiemr equation according to Chapter 7.4.2 in
  !> PhD thesis CRG.
  !>
  !> @note
  !> A pseudocode is given in the dissertation of CRG.
  !> @note
  !> This routine is based on a routine, which Gregor Michaliczek has used for his dissertation.
  !>
  !> @todo spin
  !> @todo LOs
  !>
  !> @details
  !> See 7.101 and 7.102 (dissertation CRG)
  !>
  !> @param[in]   stars      : Stars type, see types.f90.
  !> @param[in]   dimens     : Dimension type, see types.f90.
  !> @param[in]   nmatBra    : Number of basis functions at k + q (k') including LOs
  !> @param[in]   nrBraBands : Number of all bands at k + q (k')
  !> @param[in]   nrKetBands : Number of occupied bands at k
  !> @param[in]   ikpt       : Index of k-point in k-point set
  !> @param[in]   ikpq       : Index of k + q (k' backfolded) in k-point set
  !> @param[in]   ngpqdp     : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   iqpt       : Index of q in q-point set
  !> @param[in]   GbasBra    : G-basis vectors at k + q (k')
  !> @param[in]   GbasKet    : G-basis vectors at k
  !> @param[in]   nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]   vEff1IR    : Plane-wave coefficients of the interstitial effective linear variation of the (in first iteration
  !>                           the external) effective potential
  !> @param[in]   zBra       : Kohn-Sham eigenvectors at k+q
  !> @param[in]   zKet       : Kohn-Sham eigenvectors at k
  !> @param[in]   gpqdp      : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]   kpq2kPrVec : Backfolding vector from k + q in 2nd Brillouin zone to k' in 1st Brillouin zone
  !> @param[out]  MeVeff1IR  : Interstitial part of Optimized Hellmann-Feynman Contribution to the Sternheiemr equation according to
  !>                           Chapter 7.4.2 in PhD thesis CRG.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcMEPotIR( stars, dimens, GbasBra, GbasKet, nv, vEff1IR, zBra, zKet, gpqdp, nmatBra, nrBraBands, nrKetBands, ikpt,&
      & iqpt, ikpq, ngpqdp, MeVeff1IR, kpq2kPrVec ,idirx, mat_elH)

#include "recycledRoutines/cpp_arch.h"
#include "recycledRoutines/cpp_double.h"

      use m_types, only : t_stars, t_dimension
      use m_jpConstants, only : compPhon

      implicit none

      include "fftw3.f" ! include FFTW constants

      ! Type parameter
      type(t_stars),                 intent(in)  :: stars
      type(t_dimension),             intent(in)  :: dimens

      ! Scalar parameter
      integer,                       intent(in)  :: nmatBra
      integer,                       intent(in)  :: nrBraBands
      integer,                       intent(in)  :: nrKetBands
      integer,                       intent(in)  :: ikpt
      integer,                       intent(in)  :: ikpq
      integer,                       intent(in)  :: ngpqdp
      integer,                       intent(in)  :: iqpt
      integer, optional,             intent(in)  :: idirx

      ! Array parameter
      integer,                       intent(in)  :: GbasBra(:, :)
      integer,                       intent(in)  :: GbasKet(:, :)
      integer,                       intent(in)  :: nv(:, :)
      complex,                       intent(in)  :: vEff1IR(:)
      MCOMPLEX,                      intent(in)  :: zBra(:,:)
      MCOMPLEX,                      intent(in)  :: zKet(:,:)
      integer,                       intent(in)  :: gpqdp(:, :)
      integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
      MCOMPLEX,                      intent(out) :: MeVeff1IR(:,:)
      complex, optional,             intent(inout) :: mat_elH(:,:)

      ! Scalar variables
      integer                                    :: i
      integer                                    :: j
      integer                                    :: iv
      integer                                    :: il, im, in
      integer                                    :: idir
      integer                                    :: ifftds
      integer                                    :: iG1
      integer*8                                     backwardPlan
      integer*8                                     forwardPlan


      ! Array variables
      integer,           allocatable             :: iv1d(:,:)
      integer,           allocatable             :: iv1dBra(:,:)
      integer,           allocatable             :: igfft(:)
      complex,           allocatable             :: vEff1(:)
      complex,           allocatable             :: tempGrid(:)
      complex,           allocatable             :: zFFTBox(:)
      complex,           allocatable             :: vThetaZ(:), V1T0(:), MEVeff1IRG(:,:)
      integer                                    :: gMirr(3)
      integer                                    :: nfft(3)

      ! For BLAS
      complex                                       CPP_BLAS_cdotc
      external                                      CPP_BLAS_cdotc

      ! For FFTW
      external                                      dfftw_plan_dft_3d
      external                                      dfftw_execute_dft
      external                                      dfftw_destroy_plan

      ! Initialization
      allocate( iv1d(dimens%nvd, 1) )
      allocate( iv1dBra(dimens%nvd, 1) )
      allocate( vThetaZ(dimens%nbasfcn) )
      allocate( v1T0(nv(1, ikpq)) )
      allocate( MEVeff1IRG(nv(1, ikpq),nv(1, ikpt)) )
      iv1d(:, :) = 0
      iv1dBra(:, :) = 0
      vThetaZ(:) = cmplx(0., 0.)
      MeVeff1IR(:, :) = cmplx(0., 0.)

      ! Standard size of FFT mesh in Fleur
      nfft(1) = 3 * stars%k1d
      nfft(2) = 3 * stars%k2d
      nfft(3) = 3 * stars%k3d

      ifftds = product(nfft)

      ! Generate mapping array for putting potential variation to FFT mesh.
      allocate( igfft(ngpqdp) )
      gMirr = 0
      igfft = 0
      do iG1 = 1, ngpqdp
        do idir = 1, 3
          gMirr(idir) = gpqdp(idir, iG1)
          if ( gMirr(idir) < 0 ) then
            gMirr(idir) = gMirr(idir) + nfft(idir)
          end if
        end do
        igfft(iG1) = gMirr(1) + nfft(1) * gMirr(2) + nfft(1) * nfft(2) * gMirr(3)
        gMirr = 0
      end do

      ! Generate FFT plan for potential variation, this is required when using the fftw library.
      allocate (vEff1(0:ifftds-1))
      vEff1 = (0.0,0.0)
      call Dfftw_plan_dft_3d(backwardPlan, nfft(1),nfft(2),nfft(3), vEff1, vEff1, FFTW_BACKWARD, FFTW_ESTIMATE)
      vEff1 = (0.0,0.0)

      ! Map effective potential to FFT array
      do iG1 = 1, ngpqdp
         vEff1(igfft(iG1)) = vEff1IR(iG1)
      end do

      ! Perform 3D FFT on "vEff1" from reciprocal space to real space and destroy the plan
      call Dfftw_execute_dft(backwardPlan, vEff1, vEff1)
      call Dfftw_destroy_plan(backwardPlan)

      ! Generate the mapping array between FFT mesh and plane-wave expansion for the ket wave function
      do iv = 1, nv(1, ikpt)
         il = GbasKet(1, iv)
         im = GbasKet(2, iv)
         in = GbasKet(3, iv)
         if (il < 0) then
           il = il + nfft(1)
         end if
         if (im < 0) then
            im = im + nfft(2)
         end if
         if (in < 0) then
            in = in + nfft(3)
         end if
         iv1d(iv, 1) = il + nfft(1) * im + nfft(1) * nfft(2) * in
      end do

      ! Generate the mapping array between FFT mesh and plane-wave expansion for the bra wave function. As we are using k', we have
      ! to shift the G-vectors here, so that the correct G-vectors are selected.
      do iv = 1, nv(1, ikpq)
         il = GbasBra(1, iv) + kpq2kPrVec(1, ikpt, iqpt)
         im = GbasBra(2, iv) + kpq2kPrVec(2, ikpt, iqpt)
         in = GbasBra(3, iv) + kpq2kPrVec(3, ikpt, iqpt)
         if (il < 0) then
            il = il + nfft(1)
         end if
         if (im < 0) then
            im = im + nfft(2)
         end if
         if (in < 0) then
            in = in + nfft(3)
         end iF
         iv1dBra(iv, 1) = il + nfft(1) * im + nfft(1) * nfft(2) * in
      end do


      ! Generate FFT plan for cumulative variable containing step function, linear potential variation and ket wave function,
      ! this is required when using the fftw library.
      allocate (zFFTBox(0:ifftds - 1))
      zFFTBox = (0.0,0.0)
      CALL Dfftw_plan_dft_3d(backwardPlan, nfft(1), nfft(2), nfft(3), zFFTBox, zFFTBox, FFTW_BACKWARD, FFTW_MEASURE)
      zFFTBox = (0.0,0.0)
      CALL Dfftw_plan_dft_3d(forwardPlan, nfft(1), nfft(2), nfft(3), zFFTBox, zFFTBox, FFTW_FORWARD, FFTW_MEASURE)

      if (compPhon) then
        do i = 1, nv(1, ikpt)
          ! Map ket basis function expansion coefficient from plane-wave expansion representation to FFT mesh of cumultative variable
          zFFTBox = (0., 0.)
          zFFTBox(iv1d(i, 1)) = (1., 0.)

          ! Perform 3D FFT on cumultative variable from reciprocal space to real space and destroy the plan
          call Dfftw_execute_dft(backwardPlan, zFFTBox, zFFTBox)

          ! Multiply to the ket wave function expansion coefficients the linear variation of the effective potential and the step
          ! function in real space exploiting the convolution theorem
          do j = 0, ifftds - 1
            zFFTBox(j) = zFFTBox(j) * vEff1(j) * stars%ufft(j)
          end do

          call dfftw_execute_dft(forwardPlan, zFFTBox, zFFTBox)
          v1T0 = cmplx(0.0,0.0)
          do iv = 1, nv(1, ikpq)
            v1T0(iv) = zFFTBox(iv1dBra(iv, 1))
            v1T0(iv) = v1T0(iv) / ifftds
          end do
          do j = 1, nv(1,ikpq)
            MeVeff1IRG(j,i) = v1T0(j)
          end do
        end do
        !if (ikpt.eq.1.and.idirx.eq.1) then
        !  open(111,file='000_ME_V1Thet0_IR',form='FORMATTED',position='append',action='WRITE',status='REPLACE')
        !else
        !  open(111,file='000_ME_V1Thet0_IR',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
        !end if
        !do j = 1, nv(1,ikpt)
        !  do i = 1, nv(1,ikpq)
        !    write(111,*) i, j, idirx
        !    write(111,*) real(MeVeff1IRG(i,j)), aimag(MeVeff1IRG(i,j))
        !  end do
        !end do
        !close(111)

        if (present(mat_elH)) then
          mat_elH(:,:) = mat_elH(:,:) + MeVeff1IRG(:,:)
        end if
      end if
      do i = 1, nrKetBands
         ! Map ket basis function expansion coefficient from plane-wave expansion representation to FFT mesh of cumultative variable
         zFFTBox = (0., 0.)
         do iv = 1, nv(1, ikpt)
            zFFTBox(iv1d(iv, 1)) = zKet(iv, i)
         end do

         ! Perform 3D FFT on cumultative variable from reciprocal space to real space and destroy the plan
         call Dfftw_execute_dft(backwardPlan, zFFTBox, zFFTBox)

         ! Multiply to the ket wave function expansion coefficients the linear variation of the effective potential and the step
         ! function in real space exploiting the convolution theorem
         do j = 0, ifftds - 1
            zFFTBox(j) = zFFTBox(j) * vEff1(j) * stars%ufft(j)
         end do

         ! Perform 3D FFT on cumultative variable from real space to reciprocal space
         call dfftw_execute_dft(forwardPlan, zFFTBox, zFFTBox)

         ! Map the cumultative variable back to the expansion of the plane waves and fix a factor ifftds that has emerged due to the
         ! FFT. The cumultative variable corresponds to 7.103 (dissertation CRG).
         vThetaZ = cmplx(0.0,0.0)
         do iv = 1, nv(1, ikpq)
            vThetaZ(iv) = zFFTBox(iv1dBra(iv, 1))
            vThetaZ(iv) = vThetaZ(iv) / ifftds
         end do

         ! Project onto bra wave function corresponding to 7.102i (dissertation CRG)
         do j = 1, nrBraBands
         ! The zs are multiplied by 1 / omtil because the overlap matrix gives 1.
         ! We convoluted the Psiket with the warped potential so have in principle the multiplication of two quantities in the
         ! interstitial which have to be multiplied by an omtil still.
         ! conjugation of cdotc zBra is being conjugated automatically
            MeVeff1IR(j,i) = MeVeff1IR(j, i) + CPP_BLAS_cdotc(nmatBra, zBra(1, j),1, vThetaZ, 1)!nmatBra von bra!
         end do

      end do

      ! Destroy plans so that there is no garbage
      call dfftw_destroy_plan(forwardPlan)
      call dfftw_destroy_plan(backwardPlan)

  end subroutine calcMEPotIR


  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calculates the tlmplm integrals for the matrix element <Ψ|Veff1 + gradVeff0|Ψ>_α (recycled subroutine tlmplm from Fleur).
  !>
  !> @param[in]  atoms     : Atoms type, see types.f90
  !> @param[in]  lathar    : Lattice harmonics type, see types.f90.
  !> @param[in]  dimens    : Dimension type, see types.f90.
  !> @param[in]  enpara    : Energy parameter type, see types.f90.
  !> @param[in]  usdus     : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]  input     : Input type, see types.f90.
  !> @param[out] td        : Tlmplm matrix type for Sternheimer Hellmann-Feynman muffin-tin matrix element, see types.f90
  !> @param[in]  jspin     : Physical spin
  !> @param[in]  jsp       : Spin index for data
  !> @param[in]  idir      : Index of displacement direction
  !> @param[in]  vr        : Sum of the spherical harmonic expansion coefficients for the gradient of the unperturbed external
  !>                         (first iteration) or effective(regular/final iteration) potential and its linear variation. This
  !>                         potential is defined in 7.105b (dissertation CRG).
  !> @param[in]  rbas1     : Large components of radial solution, its energy derivative and u_LO
  !> @param[in]  rbas2     : Small components of radial solution, its energy derivative and u_LO
  !> @param[in]  uuilon    : Overlap integral between the radial functions of the integral (multiplied by ulo_der) of a local
  !>                         orbital and the flapw radial function with the same l
  !> @param[in]  duilon    : Overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                         derivative of the flapw radial function with the same l
  !> @param[in]  ulouilopn : Overlap integral between the radial functions of the integral of a local orbital and another local
  !>                         orbital with the same l.
  !> @param[in]  ilo2p     : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[out] nlo_atom  : Contains information about number of LOs at every atom
  !>
  !> @details
  !> From tlmplm: Sets up the t(l'm',lm) matrices for each atom. These matrices are k-point independent quantities.
  !> m. weinert 1986. See also Equation 7.106, Figure 7.6, 7.107 and 7.108, Algorithms 10 - 13 (dissertation CRG)
  !>
  !> @note
  !> In order to set up the total matrix also the potential with 7.105b also has to be decorated with a factor as in 7.108
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine tlmplm4V( atoms, lathar, dimens, enpara, usdus, input, td, jspin, jsp, idir, vr, rbas1, rbas2, uuilon, duilon, &
      & ulouilopn, ilo2p, nlo_atom )

!    use m_intgr, only : intgr3
    use m_intgr, only : intgr3LinIntp
    use m_gaunt, only : gaunt1
    use m_types, only : t_atoms, t_sphhar, t_enpara, t_usdus, t_input, t_tlmplm, t_dimension
    use m_juDFT_NOstopNO, only : juDFT_error
    use m_jpConstants, only : compPhon

    implicit none

    ! Type parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_sphhar),              intent(in)  :: lathar
    type(t_dimension),           intent(in)  :: dimens
    type(t_enpara),              intent(in)  :: enpara
    type(t_usdus),               intent(in)  :: usdus
    type(t_input),               intent(in)  :: input
    type(t_tlmplm),              intent(out) :: td

    ! Scalar Parameters
    integer,                     intent(in)  :: jspin !physical spin
    integer,                     intent(in)  :: jsp ! spin index for data
    integer,                     intent(in)  :: idir

    ! Array Arguments
    complex,                     intent(in)  :: vr(:,:,:, :)
    real,                        intent(in)  :: rbas1(:,:,0:,:,:)
    real,                        intent(in)  :: rbas2(:,:,0:,:,:)
    real,                        intent(in)  :: uuilon(:,:)
    real,                        intent(in)  :: duilon(:,:)
    real,                        intent(in)  :: ulouilopn(:,:,:)
    integer,                     intent(in)  :: ilo2p(:, :)
    integer,        allocatable, intent(out) :: nlo_atom(:)

    ! Local Scalars
    complex                                  :: cil
    real                                     :: tempReal
    real                                     :: tempImag
    real                                     :: tempRbas
    integer                                  :: i
    integer                                  :: l
    integer                                  :: l2
    integer                                  :: lamda
    integer                                  :: lm
    integer                                  :: lmin
    integer                                  :: lmin0
    integer                                  :: lmp
    integer                                  :: lmpl
    integer                                  :: lmplm
    integer                                  :: lmx
    integer                                  :: lmxx
    integer                                  :: lp
    integer                                  :: lp1
    integer                                  :: lpl
    integer                                  :: mp
    integer                                  :: mu
    integer                                  :: n
    integer                                  :: na
    integer                                  :: m
    integer                                  :: err
    integer                                  :: mlotot
    integer                                  :: mlolotot
    integer                                  :: mlot_d
    integer                                  :: mlolot_d
    integer                                  :: ieqat
    integer                                  :: lmsph
    !     ..
    !     .. Local Arrays ..
    integer,        allocatable              :: indt(:)
    complex,        allocatable              :: dvd(:, :)
    complex,        allocatable              :: dvu(:, :)
    complex,        allocatable              :: uvd(:, :)
    complex,        allocatable              :: uvu(:, : )
    real,           allocatable              :: xReal(:)
    real,           allocatable              :: xImag(:)

    ! Initialization which is made in eigen, only relevant for the LOs
    mlotot = 0
    mlolotot = 0
    do n = 1, atoms%ntype
      do na = 1, atoms%neq(n)
       mlotot = mlotot + atoms%nlo(n)
       mlolotot = mlolotot + atoms%nlo(n) * (atoms%nlo(n) + 1 ) / 2
      end do
    end do
    mlot_d = max(mlotot, 1)
    mlolot_d = max(mlolotot, 1)

    err = 0
    allocate( indt(0:dimens%lmplmd) )
    allocate( dvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( dvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( uvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( uvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( xReal(atoms%jmtd), xImag(atoms%jmtd) )
    allocate( td%tuu(0:dimens%lmplmd, atoms%nat, 1), stat=err )
    allocate( td%tud(0:dimens%lmplmd, atoms%nat, 1), stat=err )
    allocate( td%tdd(0:dimens%lmplmd, atoms%nat, 1), stat=err )
    allocate( td%tdu(0:dimens%lmplmd, atoms%nat, 1), stat=err )
    allocate( td%tdulo(0:dimens%lmd, -atoms%llod:atoms%llod, mlot_d, 1), stat=err )
    allocate( td%tuulo(0:dimens%lmd, -atoms%llod:atoms%llod, mlot_d, 1), stat=err )
    allocate( td%tuloulo(-atoms%llod:atoms%llod, -atoms%llod:atoms%llod, mlolot_d, 1), stat=err )
    allocate( td%ind(0:dimens%lmd, 0:dimens%lmd, atoms%nat, 1), stat=err )
    if (err.ne.0) then
       write (*,*) 'eigen: an error occured during allocation of'
       write (*,*) 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
       CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="eigen")
    end if

    indt(:) = 0
    dvd(:, :) = cmplx(0., 0.)
    dvu(:, :) = cmplx(0., 0.)
    uvd(:, :) = cmplx(0., 0.)
    uvu(:, :) = cmplx(0., 0.)
    xReal(:) = 0.
    xImag(:) = 0.
    td%tuu(:, :, :) = cmplx(0.,0.)
    td%tud(:, :, :) = cmplx(0.,0.)
    td%tdd(:, :, :) = cmplx(0.,0.)
    td%tdu(:, :, :) = cmplx(0.,0.)
    td%tdulo(:, :, :, :) = cmplx(0.,0.)
    td%tuulo(:, :, :, :) = cmplx(0.,0.)
    td%tuloulo(:, :, :, :) = cmplx(0., 0.)
    td%ind(:, :, :, :) = -9999

    na = 0
    do n = 1, atoms%ntype
      do ieqat = 1, atoms%neq(n)
        na = na + 1
        do lp = 0, atoms%lmax(n)
          ! Generation of triangular numbers. NOTE: lp1 is integer
          lp1 = (lp * (lp + 1)) / 2
          ! Loop over ket but only in a triangular mode so l <= l'
          do l = 0, lp
            lpl = lp1 + l
            ! Loop over components of the potential
            do lamda = 0, atoms%lmax(n)
              lmin = lp - l
              lmx = lp + l
              ! We only have a contribution according to the Gaunt selection rules. The triangular condition can be derived from an
              ! inequalities system so that we can set up the triangular conditions for every index of the Gaunt coefficients.
              ! Furthermore l'+l+lamda should be even
              if ((mod(lamda + lmx, 2) .eq. 1) .or. (lamda.lt.lmin) .or. (lamda.gt.lmx)) then
                do m = -lamda, lamda
                  lm = lamda * (lamda + 1) + m + 1
                  uvu(lpl, lm) = cmplx(0.0, 0.0)
                  dvd(lpl, lm) = cmplx(0.0, 0.0)
                  uvd(lpl, lm) = cmplx(0.0, 0.0)
                  dvu(lpl, lm) = cmplx(0.0, 0.0)

                  !if (compPhon) then
                  !  write(111,*) lm, lp, l, n , na
                  !  write(111,*) 1, 1, real(uvu(lpl,lm)), aimag(uvu(lpl,lm))
                  !  write(111,*) 2, 1, real(dvu(lpl,lm)), aimag(dvu(lpl,lm))
                  !  write(111,*) 1, 2, real(uvd(lpl,lm)), aimag(uvd(lpl,lm))
                  !  write(111,*) 2, 2, real(dvd(lpl,lm)), aimag(dvd(lpl,lm))
                  !end if
                end do ! m
              else
                xReal = 0
                xImag = 0
                do m = -lamda, lamda
                  lm = lamda * (lamda + 1) + m + 1

                  !if (compPhon) then
                  !  write(111,*) lm, lp, l, n , na
                  !end if

                  ! Set up integrals Sigma 7.107d/e (dissertation CRG). Note, that the effective potential variation is complex in
                  ! general, therefore we need a two-fold integration using a real integration routine.

                  ! Calculate the integral <u|V|u>
                  do i = 1,atoms%jri(n)
                    tempRbas = ( rbas1(i, 1, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                    xReal(i) =  tempRbas * real(vr(i,lm,idir,na))
                    xImag(i) =  tempRbas * aimag(vr(i,lm,idir,na))
                  end do ! i
                  call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                  call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                  uvu(lpl,lm) = cmplx(tempReal, tempImag)

                  ! Calculate the integral <uDot|V|u>
                  do i = 1,atoms%jri(n)
                    tempRbas = ( rbas1(i, 2, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                    xReal(i) = tempRbas * real(vr(i,lm,idir,na))
                    xImag(i) = tempRbas * aimag(vr(i,lm,idir,na))
                  end do ! i
                  call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                  call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                  dvu(lpl,lm) = cmplx(tempReal, tempImag) ! changed

                  ! Calculate the integral <u|V|uDot>
                  do i = 1,atoms%jri(n)
                    tempRbas = (rbas1(i, 1, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                    xReal(i) = tempRbas * real(vr(i, lm, idir,na))
                    xImag(i) = tempRbas * aimag(vr(i, lm, idir,na))
                  end do ! i
                  call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                  call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                  uvd(lpl,lm) = cmplx(tempReal, tempImag)

                  ! Calculte the integral <uDot|V|uDot>
                  do i = 1,atoms%jri(n)
                    tempRbas = (rbas1(i, 2, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                    xReal(i) = tempRbas * real(vr(i, lm, idir, na))
                    xImag(i) = tempRbas * aimag(vr(i, lm, idir, na))
                  end do ! i
                  call intgr3LinIntp(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal, 1)
                  call intgr3LinIntp(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag, 1)
                  dvd(lpl, lm) = cmplx(tempReal, tempImag)

                  !if (compPhon) then
                  !  write(111,*) 1, 1, real(uvu(lpl,lm)), aimag(uvu(lpl,lm))
                  !  write(111,*) 2, 1, real(dvu(lpl,lm)), aimag(dvu(lpl,lm))
                  !  write(111,*) 1, 2, real(uvd(lpl,lm)), aimag(uvd(lpl,lm))
                  !  write(111,*) 2, 2, real(dvd(lpl,lm)), aimag(dvd(lpl,lm))
                  !end if

                end do ! m
              end if ! Gaunt selection rules
            end do ! lambda
          end do ! l
        end do ! lp
        ! Generate the various t(l'm',lm) matrices for l'm'.ge.lm
        indt(:) = 0
        ! Loop over l'm'
        do lp = 0, atoms%lmax(n)
          ! NOTE: lp1 is integer
          lp1 = (lp * (lp + 1)) / 2
          do mp = -lp, lp
            lmp = lp * (lp + 1) + mp
            lmpl = (lmp * (lmp + 1)) / 2
            do lamda = 0, atoms%lmax(n)
              lmin0 = abs(lp - lamda)
              ! l' - lambda > l' and -l' + lambda > l' leads to 0 <= lamda <= 2 l'
              if (lmin0 .GT. lp) cycle
              ! lmxx is in principle lp. To ensure, that the oddness of lambda is always compensated (here strictly speaking only at
              ! lmxx) by l, the modulo operation is subtracted. It is subtracted to not exceed lmxx. The lmxx boundary was given by
              ! Gaunt selection rules, actually reading l' + lambda but as we have l <= l' lmax is only l'. We group lambda and the
              ! modulo operation in lmxx + l' + lambda, so that we only have to ensure l + l' is even.
              lmxx = lp - mod(lamda, 2)
              do mu = -lamda, lamda
                ! Collection index of lamda and mu!
                lmsph = lamda * (lamda + 1) + 1 + mu
                ! selection rule mp = m + mu!
                m = mp - mu
                ! lmin = max(|l' - lamda|, |m' - mu|) to ensure that |m = m' - mu| <= l while fulfilling m = m' - mu.
                ! Imagine for example l' = 1, lamda = 1, m' = -1,mu = 1. Independent of the selection rule, l >= |m|
                lmin = max(lmin0, abs(m))
                ! Serves only for ensuring l + l' + lambda = even, if either lmxx is odd or lmin odd, with l2, lmxx gives its
                ! oddness or eveness to lmin, so as lmxx was l' and l is interated in steps of size 2 l + l' is always even, because
                ! either l and l' are even or l and l' are odd, therefore l + l' + lambda is even because the eveness is ensured by
                ! the modulo operation mod(lambda, 2).
                l2 = abs(lmxx - lmin)
                ! Corrects lmin so that mod(lmin, 2) = mod(lmxx, 2), lmin is corrected upwards to not exceed the limits given by
                ! either the Gaunt selection rules (soft) or l >= |m| (hard) condition.
                lmin = lmin + mod(l2, 2)
                do l = lmin, lmxx, 2
                  ! collection index of l and m
                  lm = l * (l + 1) + m
                  ! lm is always <= lmp!
                  if (lm .GT. lmp) cycle
                  ! index for uvu, dvd, uvd, dvu to access only those relevant after the Gaunt selection rules
                  lpl = lp1 + l
                  ! index similiar to lpl but with ms in it (after the selection Gaunt selection rules)
                  lmplm = lmpl + lm
                  cil = gaunt1(lp,lamda,l,mp,mu,m,atoms%lmaxd)
                  td%tuu(lmplm,na,jsp) = td%tuu(lmplm,na,jsp) + cil*uvu(lpl,lmsph)
                  td%tdd(lmplm,na,jsp) = td%tdd(lmplm,na,jsp) + cil*dvd(lpl,lmsph)
                  td%tud(lmplm,na,jsp) = td%tud(lmplm,na,jsp) + cil*uvd(lpl,lmsph)
                  td%tdu(lmplm,na,jsp) = td%tdu(lmplm,na,jsp) + cil*dvu(lpl,lmsph)
                  ! Logical matrix where there are non-vanishing matrix entries in the tlmplm matrices
                  indt(lmplm) = 1
                end do ! l
              end do ! mem
            end do ! lh
          end do ! mp
        end do ! lp

        ! set up mapping array
        do lp = 0,atoms%lmax(n)
          do mp = -lp,lp
            lmp = lp* (lp + 1) + mp
            do l = 0, atoms%lmax(n)
              do m = -l, l
                lm = l * (l + 1) + m
                if (lmp.ge.lm) then
                  lmplm = (lmp * (lmp + 1)) / 2 + lm
                  if (indt(lmplm).NE.0) then
                    td%ind(lmp,lm,na,jsp) = lmplm
                  else
                    td%ind(lmp,lm,na,jsp) = -9999
                  end if
                else
                  ! As we have lm > lmp here (this was not calculated within the routine), we have to transpose, i.e. interchanging
                  ! lmp and lm for the packed storage index.
                  lmplm = (lm * (lm + 1)) / 2 + lmp
                  if (indt(lmplm).NE.0) then
                    td%ind(lmp,lm,na,jsp) = -lmplm
                  else
                    td%ind(lmp,lm,na,jsp) = -9999
                  end if
                end if
              end do ! m
            end do ! l
          end do ! mp
        end do ! lp

        ! set up the t-matrices for the local orbitals, if there are any
        if (.false.) then
        !if (atoms%nlo(n).ge.1) then ! This is to be commented in, when LOs are implemented again
          call tlo4V(atoms, enpara, lathar, usdus, input, td, jspin, jsp, n, na, vr(:, :, idir, na), rbas1, rbas2, uuilon, duilon, &
                                                                                                       & ulouilopn, ilo2p, nlo_atom)
        else
          allocate( nlo_atom(atoms%nat) )
          nlo_atom = 0
        end if
      end do ! ieqat
    end do ! n

  end subroutine tlmplm4V


  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calculates the tlo integrals for the matrix element <Ψ|Veff1 + gradVeff0|Ψ>_α (recycled subroutine tlmplm from Fleur).
  !>
  !> @details
  !> sets up the extra t-matrix elements due to the local orbitals.
  !> only non=zero elements are calculated
  !> ************** ABBREVIATIONS *****************************************
  !> tuulo      : t-matrix element of the lo and the apw radial fuction
  !> tdulo      : t-matrix element of the lo and the energy derivativ of
  !>              the apw radial fuction
  !> tuloulo    : t-matrix element of two los
  !> ***********************************************************************
  !>
  !> p.kurz jul. 1996
  !>
  !> @todo review documentation
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine tlo4V(atoms, enpara, lathar, usdus, input, tlmplm, jspin, jsp, ntyp, na, vr, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, nlo_atom)
    !use m_intgr, only : intgr3
    use m_intgr, only : intgr3LinIntp
    use m_gaunt, only : gaunt1
    use m_types
    use m_juDFT_NOstopNO, only : juDFT_error
    use m_jpConstants, only : iu
    implicit none

    type(t_atoms),       intent(in)    :: atoms
    type(t_enpara),      intent(in)    :: enpara
    type(t_sphhar),      intent(in)    :: lathar
    type(t_usdus),       intent(in)    :: usdus
    type(t_input),       intent(in)    :: input
    type(t_tlmplm),      intent(inout) :: tlmplm
    !     ..
    !     .. Scalar Arguments ..
    integer,             intent (in)   :: jspin,jsp,ntyp ,na
    !     ..
    !     .. Array Arguments ..
    integer,             intent(in)    :: ilo2p(:, :)
    complex,             intent (in)   :: vr(atoms%jmtd,(atoms%lmaxd + 1)**2 ) ! changed
    real,                intent(in)    :: rbas1(:,:,0:,:,:)
    real,                intent(in)    :: rbas2(:,:,0:,:,:)
    real,                intent (in)   :: uuilon(atoms%nlod,atoms%ntypd),duilon(atoms%nlod,atoms%ntypd)
    real,                intent (in)   :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntypd)
    integer, allocatable, intent(out)  :: nlo_atom(:)
    !     ..
    !     .. Local Scalars ..
    complex                            :: cil
    integer                            :: lmpp ! changed
    integer                            :: i,l,lh,lm ,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,lpmin,lpmin0,lpp ,mem,mp,mpp,m,lmx,mlo,mlolo !changed
    real                               :: tempReal, tempImag  !changed
    integer                            :: mytype
    integer                            :: myatom
    integer                            :: myeqat
    !     ..
    !     .. Local Arrays ..
    real                               :: xReal(atoms%jmtd),xImag(atoms%jmtd) ! changed
    complex                            :: ulovulo(atoms%nlod*(atoms%nlod+1)/2,(atoms%lmaxd + 1)**2) ! changed
    complex                            :: uvulo(atoms%nlod,0:atoms%lmaxd,(atoms%lmaxd + 1)**2),dvulo(atoms%nlod,0:atoms%lmaxd,(atoms%lmaxd + 1)**2) ! changed
    !     ..

    !todo this can be shifted to tlo4H because this is not called for every q-point
    allocate(nlo_atom(atoms%nat))
    myatom = 0
    do mytype = 1, atoms%ntype
      do myeqat = 1, atoms%neq(mytype)
        myatom = myatom + 1
        nlo_atom(myatom) = atoms%nlo(mytype)
      end do
    end do


    do lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       do lp = 0,atoms%lmax(ntyp)
          lmin = ABS(lp-l)
          !               lmin = lp - l
          lmx = lp + l
          do lpp = 0, atoms%lmax(ntyp) ! changed ! todo to l + 1????
             if ((mod(l+lp+lpp,2).eq.1) .or. (lpp.lt.lmin) .or.&
                  (lpp.GT.lmx)) then
                do mpp = -lpp, lpp ! changed
                lm = lpp * (lpp + 1) + 1 + mpp ! changed
                uvulo(lo,lp,lm) = 0.0 ! changed
                dvulo(lo,lp,lm) = 0.0 ! changed
                end do ! changed
             else
                do mpp = -lpp, lpp ! changed
                lmpp = lpp * (lpp + 1) + 1 + mpp ! changed
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,1,lp,ntyp,1)*rbas1(i, ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+ rbas2(i, 1, lp, ntyp, 1) * rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp)) ! changed
                   xImag(i) = (rbas1(i,1,lp,ntyp,1)*rbas1(i, ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+ rbas2(i, 1, lp, ntyp, 1) * rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp)) ! changed
                end do
                call intgr3LinIntp(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal, 1) ! changed
                call intgr3LinIntp(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag, 1) ! changed
                uvulo(lo,lp,lmpp) = cmplx(tempReal, tempImag) ! changed
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,2,lp,ntyp,1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp),ntyp,1)+ rbas2(i,2,lp,ntyp,1)*rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp)) ! changed
                   xImag(i) = (rbas1(i,2,lp,ntyp,1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp),ntyp,1)+ rbas2(i,2,lp,ntyp,1)*rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp)) ! changed
                end do !changed
                call intgr3LinIntp(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal, 1) ! changed
                call intgr3LinIntp(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag, 1) ! changed
                dvulo(lo,lp,lmpp) = cmplx(tempReal, tempImag) ! changed
                end do !changed
             end if
          end do
       end do
    end do
    loplo = 0
    do lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       do lo = 1,lop
          l = atoms%llo(lo,ntyp)
          loplo = loplo + 1
          IF (loplo>size(ulovulo,1))  CALL juDFT_error("loplo too large!!!" ,calledby ="tlo")
          do lpp = 0, atoms%lmax(ntyp) ! changed
             lmin = ABS(lp - l)
             lmx = lp + l
             if ((mod(l+lp+lpp,2).eq.1).or.(lpp.lt.lmin).or.(lpp.gt.lmx)) then
               do mpp = -lpp, lpp ! changed
               lmpp = lpp * (lpp + 1) + 1 + mpp ! changed
                ulovulo(loplo,lmpp) = 0.0 ! changed
               end do ! changed
             else
               do mpp = -lpp, lpp ! changed
               lmpp = lpp * (lpp + 1) + 1 + mpp ! changed
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+rbas2(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas2(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp)) ! changed
                   xImag(i) = (rbas1(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+rbas2(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas2(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp)) ! changed
                end do ! changed
                call intgr3LinIntp(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal, 1) ! changed
                call intgr3LinIntp(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag, 1) ! changed
                ulovulo(loplo,lmpp) = cmplx(tempReal, tempImag) !changed
                end do ! changed
             end if
          end do
       end do
    end do
    !---> generate the different t matrices
    !---> but first initialize them ( done in eigen )
    !
    !---> generate the t-matrices. for optimal performance consider only
    !---> those combinations of l,l',l'',m,m',m'' that satisfy the three
    !---> conditions for non-zero gaunt-coeff. i.e.
    !---> |l - l''| <= l' <= l + l'' (triangular condition)
    !---> m' = m + m'' and l + l' + l'' even
    !---> loop over the local orbitals
    mlo=sum(nlo_atom(:na-1))
    do lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       do m = -l,l
          !--->       loop over the lattice harmonics
          do lpp = 0, atoms%lmax(ntyp) ! changed
             lpmin0 = ABS(l-lpp)
             lpmax0 = l + lpp
             !--->          check that lpmax is smaller than the max l of the
             !--->          wavefunction expansion at this atom
             lpmax = MIN(lpmax0,atoms%lmax(ntyp))
             !--->          make sure that l + l'' + lpmax is even
             lpmax = lpmax - MOD(l+lpp+lpmax,2)
             do mpp = -lpp, lpp ! changed
            !    mpp = lathar%mlh(mem,lh,atoms%ntypsy(na)) ! changed
            lmpp = lpp * (lpp + 1) + 1 + mpp ! changed
                mp = m + mpp
                lpmin = MAX(lpmin0,ABS(mp))
                !--->             make sure that l + l'' + lpmin is even
                lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                !--->             loop over l'
                do lp = lpmin,lpmax,2
                   lmp = lp* (lp+1) + mp
              !     cil = gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd + 1) ! changed
                   cil = gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd) ! changed
                   tlmplm%tuulo(lmp,m,lo+mlo,jsp) = &
                        tlmplm%tuulo(lmp,m,lo+mlo,jsp) + cil*uvulo(lo,lp,lmpp) ! changed
                   tlmplm%tdulo(lmp,m,lo+mlo,jsp) = &
                        tlmplm%tdulo(lmp,m,lo+mlo,jsp) + cil*dvulo(lo,lp,lmpp) ! changed
                end do
             end do
          end do
       end do
    end do
    !---> generate the t-matrix including two local orbitals for lo' >= lo
    !---> loop over lo'
    !mlolo=dot_product(atoms%nlo(:ntyp-1),atoms%nlo(:ntyp-1)+1)/2
    mlolo=dot_product(nlo_atom(:na-1),nlo_atom(:na-1)+1)/2
    do lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       do mp = -lp,lp
          !--->       loop over the lattice harmonics
          do lpp = 0, atoms%lmax(ntyp)  ! changed
             do mpp = -lpp, lpp ! changed
             lmpp = lpp * (lpp + 1) + 1 + mpp ! changed
                m = mp - mpp
                !--->             loop over lo
                do lo = 1,lop
                   l = atoms%llo(lo,ntyp)
                   loplo = ((lop-1)*lop)/2 + lo
                   if ((abs(l-lpp).le.lp) .and. (lp.le. (l+lpp)) .and.&
                        (mod(l+lp+lpp,2).eq.0) .and. (abs(m).le.l)) then
                      cil = gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd) ! changed
                      tlmplm%tuloulo(mp,m,loplo+mlolo,jsp) = tlmplm%tuloulo(mp,m,loplo+mlolo,jsp) + cil*ulovulo(loplo,lmpp) ! changed
                   end if
                end do
             end do
          end do
       end do
    end do
  end subroutine tlo4V


  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the matrix element Psi_k+q Veff1 Psi_k in the muffin-tin, as part of the Sternheimer equation Hellmann--Feynman
  !> contribution, see 7.110 (Dissertation CRG).
  !>
  !> @details
  !> See also Algorithm 14 (dissertation CRG)
  !>
  !>  ! Type Parameters
  !>  @param[in]  atoms    : Atoms type, see types.f90
  !>  @param[in]  td4V1    : Tlmplm matrix type for Sternheimer Hellmann-Feynman muffin-tin matrix element, see types.f90
  !>  @param[in]  td4V2    : Tlmplm matrix type for Sternheimer Hellmann-Feynman muffin-tin matrix element, decorated with
  !>                         factor according to Equation 7.108, see types.f90
  !>  @param[in]  ikpt     : Index of k-point in k-point set
  !>  @param[in]  ikpq     : Index of k + q (k' backfolded) in k-point set
  !>  @param[in]  ne       : Number of eigenvalues per k-point
  !>  @param[in]  nobd     : Number of occupied bands per k-point and spin
  !>  @param[in]  mCoefBp  : Matching coefficient at k+q as in Psi^(0)
  !>  @param[in]  mCoefKb  : Matching coefficient at k as in Psi^(0)
  !>  @param[in]  nRadFun  : Number of radial functions per orbital quantum number l and atom type.
  !>  @param[in]  iloTable : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !>  @param[in]  nlo_atom : Contains information about number of LOs at every atom
  !>  @param[out] vSum     : Matrix element Psi_k+q Veff1 Psi_k in the muffin-tin, as part of the Sternheimer equation
  !>                         Hellmann--Feynman contribution, see 7.106 (Dissertation CRG).
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcVsumMT( atoms, td4V1, td4V2, ikpt, ikpq, ne, nobd, mCoefBp, mCoefKb, nRadFun, iloTable, nlo_atom, vSum , &
                       & idir, nvk, nvkq, almkg, blmkg, almkgq, blmkgq, mat_elH)

    use m_types
    use m_jpConstants, only : iu, compPhon

    implicit none

    ! Type Parameters
    type(t_atoms),              intent(in)  :: atoms
    TYPE(t_tlmplm),             intent(in)  :: td4V1
    TYPE(t_tlmplm),             intent(in)  :: td4V2

    ! Scalar Parameters
    integer,                    intent(in)  :: ikpt
    integer,                    intent(in)  :: ikpq

    ! Array Parameters
    integer,                    intent(in)  :: ne(:)
    integer,                    intent(in)  :: nobd(:, :)
    complex,                    intent(in)  :: mCoefBp(:, :, :)
    complex,                    intent(in)  :: mCoefKb(:, :, :)
    integer,                    intent(in)  :: nRadFun(0:, :)
    integer,                    intent(in)  :: iloTable(:, 0:, :)
    integer,                    intent(in)  :: nlo_atom(:)
    MCOMPLEX,                   intent(out) :: vSum(:,:)
    integer, optional,          intent(in)  :: idir
    integer, optional,          intent(in)  :: nvk, nvkq
    complex, optional,          intent(in)  :: almkg(:,:), blmkg(:,:), almkgq(:,:), blmkgq(:,:)
    complex, optional,          intent(inout) :: mat_elH(:,:)

    ! Local scalars
    integer                                 :: mlo
    integer                                 :: mlolo
    integer                                 :: nband
    integer                                 :: pband
    integer                                 :: coBsh
    integer                                 :: loBra
    integer                                 :: coKsh
    integer                                 :: loKet
    integer                                 :: lmloB
    integer                                 :: lmB
    integer                                 :: mB
    integer                                 :: lmloK
    integer                                 :: lmK
    integer                                 :: lK
    integer                                 :: mK
    integer                                 :: ind, iG, iGq
    integer                                 :: indN
    integer                                 :: itype
    integer                                 :: iatom
    integer                                 :: ieqat
    integer                                 :: lB
    integer                                 :: loBraKet
    complex                                 :: utu
    complex                                 :: dtu
    complex                                 :: utd
    complex                                 :: dtd
    complex                                 :: utulo
    complex                                 :: dtulo
    complex                                 :: ulotu
    complex                                 :: ulotd
    complex                                 :: ulotulo

    ! Local Arrays
    complex,        allocatable             :: ax(:)
    complex,        allocatable             :: bx(:)
    complex,        allocatable             :: cx(:, :)
    complex,        allocatable             :: ax2(:)
    complex,        allocatable             :: bx2(:)
    complex,        allocatable             :: vSumG(:,:)
    vSum(:, :) = cmplx(0., 0.)

    allocate(ax(nobd(ikpt, 1)), bx(nobd(ikpt, 1)))
    ax(:) = cmplx(0., 0.)
    bx(:) = cmplx(0., 0.)

    if (compPhon) then
      allocate(ax2(nvk), bx2(nvk))
      ax2(:) = cmplx(0., 0.)
      bx2(:) = cmplx(0., 0.)
      allocate(vSumG(nvkq,nvk))
      vSumG(:,:) = cmplx(0., 0.)
    end if

    iatom = 0
    do itype = 1, atoms%ntype
      if (.false.) then
        ! Determine shift of LO-type tlmplm arrays for current atom type
        if ( atoms%nlo(itype) > 0 ) then
          allocate(cx(nobd(ikpt, 1), 2: maxval(nRadFun) - 1))
        else
          allocate(cx(nobd(ikpt, 1), 2: 3))
        end if
      end if
      ! Determine shift of LO-type td4V1 arrays for current atom type
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        ! Determine shift of LO-type tlmplm arrays for current atom type
        mlo = 0
        mlolo = 0
        if(atoms%nlo(itype) > 0) then
          mlo = sum( nlo_atom(:iatom - 1) )
          mlolo = dot_product( nlo_atom(:iatom - 1), nlo_atom(:iatom - 1) + 1 ) / 2
        end if

        ! The lmloK(B) are cumultative indices incorporating not only the l and m but also the radial u and udot (later also LOs).
        lmloB = 1
        lmB = -1
        do lB = 0, atoms%lmax(itype)
          do mB = -lB,lB
            lmB = lmB + 1

            ax = cmplx(0.0, 0.0)
            bx = cmplx(0.0, 0.0)

            if (compPhon) then
              ax2 = cmplx(0.0, 0.0)
              bx2 = cmplx(0.0, 0.0)
            end if

            if (.false.) then
              cx = cmplx(0.0, 0.0)
            end if

            lmloK = 1
            lmK = -1
            do lK = 0, atoms%lmax(itype)
              do mK = -lK, lK
                lmK = lmK + 1

                ! As in tlmplm4V there is a check whether a respective matrix entry is set, the type of the index is consistent with
                ! the type of tuu tdu, tud and tdd. The diagonal is only in once and taken from the regular triangular matrix l'>=l.
                if (lmB .ge. lmK) then
                  ind = td4V1%ind(lmB, lmK, iatom, 1)
                else
                  ind = td4V2%ind(lmB, lmK, iatom, 1)
                end if
                if (ind /= -9999) THEN

                  ! For calculating the full matrix that is interated by the matching coefficients one has to unpack the tlmplm
                  ! integrals for the potential only calculated for the triangular part.
                  if (ind >= 0) THEN
                     ! The i^l and i^l' stem from the matching coefficients
                     utu = td4V1%tuu(ind, iatom, 1)
                     dtu = td4V1%tdu(ind, iatom, 1)
                     utd = td4V1%tud(ind, iatom, 1)
                     dtd = td4V1%tdd(ind, iatom, 1)

                     !if (compPhon) then
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 1
                     !  write(111,*) real(utu), aimag(utu)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 1
                     !  write(111,*) real(dtu), aimag(dtu)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 2
                     !  write(111,*) real(utd), aimag(utd)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 2
                     !  write(111,*) real(dtd), aimag(dtd)
                     !end if

                     utu=utu*(iu**(lK - lB))
                     dtu=dtu*(iu**(lK - lB))
                     utd=utd*(iu**(lK - lB))
                     dtd=dtd*(iu**(lK - lB))
                  else
                     ! the sign of the index is only a hint in which triangle we are
                     ! The i^l and i^l' stem from the matching coefficients, they are not only interated in a triangular but the
                     ! full matrix, therefore, they are not complex conjugated.
                     indN = -ind
                     utu = conjg(td4V2%tuu(indn, iatom, 1))
                     dtd = conjg(td4V2%tdd(indn, iatom, 1))
                     utd = conjg(td4V2%tdu(indn, iatom, 1))
                     dtu = conjg(td4V2%tud(indn, iatom, 1))

                     !if (compPhon) then
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 1
                     !  write(111,*) real(utu), aimag(utu)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 1
                     !  write(111,*) real(dtu), aimag(dtu)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 2
                     !  write(111,*) real(utd), aimag(utd)
                     !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 2
                     !  write(111,*) real(dtd), aimag(dtd)
                     !end if

                     utu=utu*(iu**(lK - lB))
                     dtu=dtu*(iu**(lK - lB))
                     utd=utd*(iu**(lK - lB))
                     dtd=dtd*(iu**(lK - lB))
                  end if

                  do nBand = 1, nobd(ikpt, 1)
                    ! Sum of contributions to be multiplied with an acof.
                    ax(nBand) = ax(nBand) + utu * mCoefKb(nBand, lmloK, iatom)  + utd *  mCoefKb(nBand, lmloK + 1, iatom)
                    ! Sum of contributions to be multiplied with an bcof.
                    bx(nBand) = bx(nBand) + dtu * mCoefKb(nBand, lmloK, iatom)  + dtd *  mCoefKb(nBand, lmloK + 1, iatom)

                    coKsh = 1
                    if (.false.) then
                      do loKet = 3, nRadFun(lK, itype)
                        coKsh = coKsh + 1
                        utulo = (iu**(lK - lB)) * conjg(td4V1%tuulo(lmB, mK, iloTable(loKet, lK, itype) + mlo, 1))
                        dtulo = (iu**(lK - lB)) * conjg(td4V1%tdulo(lmB, mK, iloTable(loKet, lK, itype) + mlo, 1))
                        ax(nBand) = ax(nBand) + utulo * mCoefKb(nBand, lmloK + coKsh, iatom)
                        bx(nBand) = bx(nBand) + dtulo * mCoefKb(nBand, lmloK + coKsh, iatom)
                      end do
                    end if ! LO false

                    coBsh = 1
                    if (.false.) then
                      do loBra = 3, nRadFun(lB, itype)
                        coBsh = coBsh + 1
                        !todo LO: the tuu tdu tud tdd have the inverse conjugation in hlomat
                        ! todo LO: if the LOs are sorted in a special way here even an exit statement could be possible
                        ! what is with the ls they either have to be equal or not
                        ! we only calculate a triangular matrix aren't we missing entries?
                        ! don't forget the shift mlo and mlolo
                        ! indices have to be vice versa to utulo
                        ulotu =  (iu**(lK - lB)) * td4V1%tuulo(lmK, mB, iloTable(loBra, lB, itype) + mlo, 1)
                        ulotd =  (iu**(lK - lB)) * td4V1%tdulo(lmK, mB, iloTable(loBra, lB, itype) + mlo, 1)
                        cx(nBand, coBsh) = cx(nBand, coBsh) + ulotu *  mCoefKb(nBand, lmloK, iatom) &
                                                                                        & + ulotd * mCoefKb(nBand, lmloK + 1, iatom)

                        coKsh = 1
                        do loKet = 3, nRadFun(lK, itype)
                          coKsh = coKsh + 1
                        ! This construct ist due to the fact that for LOs only the triangular matrix values are stored in tuloulo
                        if ( iloTable(loBra, lB, itype) < iloTable(loKet, lK, itype ) ) then
                          loBraKet = ( ( iloTable(loKet, lK, itype) - 1 ) * iloTable(loKet, lK, itype) ) / 2 &
                                                                                                      & + iloTable(loBra, lB, itype)
                          ulotulo = (iu**(lK - lB)) * td4V1%tuloulo(mK, mB, loBraKet + mlolo, 1)
                        else
                          loBraKet = ( ( iloTable(loBra, lB, itype) - 1 ) * iloTable(loBra, lB, itype) ) / 2 &
                                                                                                      & + iloTable(loKet, lK, itype)
                          ulotulo = (iu**(lK - lB)) * conjg(td4V1%tuloulo(mB, mK, loBraKet + mlolo, 1))
                        end if
                        cx(nBand, coBsh) = cx(nBand, coBsh) + ulotulo *  mCoefKb(nBand, lmloK + coKsh, iatom)
                        end do ! loKet
                      end do ! loBra
                    end if ! LO false

                  end do ! nBand

                  if (compPhon) then
                    do iG=1, nvk
                      ax2(iG) = ax2(iG) + utu * almkg(lmK+1, iG) + utd *  blmkg(lmK+1, iG)
                      bx2(iG) = bx2(iG) + dtu * almkg(lmK+1, iG) + dtd *  blmkg(lmK+1, iG)
                    end do
                  end if

                else
                  !if (compPhon) then
                  !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 1
                  !  write(111,*) 0, 0
                  !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 1
                  !  write(111,*) 0, 0
                  !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 2
                  !  write(111,*) 0, 0
                  !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 2
                  !  write(111,*) 0, 0
                  !end if
                end if ! ind /= -9999
                lmloK = lmloK + nRadFun(lK, itype)
              end do ! mp
            end do ! lp

            do nBand = 1, nobd(ikpt, 1)
              do pBand = 1, ne(ikpq)
                vSum(pBand, nBand) = vSum(pBand, nBand) + conjg(mCoefBp(pBand, lmloB, iatom)) * ax(nBand)&
                  & + conjg(mCoefBp(pBand, lmloB + 1, iatom)) * bx(nBand)
                if (.false.) then
                  coBsh = 1
                  do loBra = 3, nRadFun(lB, itype)
                    coBsh = coBsh + 1
                    vSum(pBand, nBand) = vSum(pBand, nBand) + cx(nBand, coBsh) * conjg(mCoefBp(pBand, lmloB + coBsh, iatom))
                  end do ! loBra
                end if ! LO false
              end do ! pband
            end do ! nBand

            if (compPhon) then
              do iGq=1, nvkq
                do iG=1, nvk
                  vSumG(iGq,iG)=vSumG(iGq,iG)+conjg(almkgq(lmB+1, iGq))*ax2(iG)+conjg(blmkgq(lmB+1, iGq))*bx2(iG)
                end do
              end do
            end if

            lmloB = lmloB + nRadFun(lB, itype)
          end do ! mB
        end do ! lB
      end do ! ieqat
      if ( allocated(cx) ) deallocate (cx)
    end do ! itype

    !if (compPhon) then
    !    if (ikpt.eq.1.and.idir.eq.1) then
    !      open(109,file='000_ME_V1_MT',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
    !    else
    !      open(109,file='000_ME_V1_MT',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')  
    !    end if
    !    do iG=1, nvk
    !      do iGq=1, nvkq
    !        write(109,*) iGq, iG, idir
    !        write(109,*) real(vSumG(iGq,iG)), aimag(vSumG(iGq,iG))
    !      end do
    !    end do
    !    close(109)
    !end if

    if (present(mat_elH)) then
      mat_elH(:,:) = mat_elH + vSumG(:,:)
    end if

  end subroutine calcVsumMT

  subroutine udHu(atoms,El,rbas1,rbas2,vEff0MT,sumVMTs,haa,dhaa)

    use m_intgr, only : intgr3LinIntp
    use m_JPConstants
    use m_types, only : t_atoms

    implicit none

    type(t_atoms),  intent(in) :: atoms

    real,    intent(in)  :: El(:, 0:, :, :), rbas1(:,:,0:,:,:), rbas2(:,:,0:,:,:)
    real,    intent(in)  :: vEff0MT(:,:,:,:)
    complex, intent(in)  :: sumVMTs(:,:,:,:)
    complex, intent(out) :: haa(:,:,0:,:,0:), dhaa(:,:,0:,:,0:)

    integer :: l, io, lp, jo, oqn_l, mqn_m, lm
    real    :: xReal(atoms%jri(1)), xImag(atoms%jri(1))
    complex :: t1(atoms%jri(1))
    real    :: tempReal, tempImag

    haa(:,:,:,:,:)  = cmplx(0., 0.)
    dhaa(:,:,:,:,:) = cmplx(0., 0.)

    !open(111,file='000_testthis0',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
    !open(112,file='000_testthis1',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
    !do lm=1, (atoms%lmaxd+1)**2
    !  do io=1, atoms%jri(1)
    !    write(111,*) vEff0MT(io,lm,1,1)
    !    write(112,*) real(sumVMTs(io,lm,1,1)), aimag(sumVMTs(io,lm,1,1))
    !  end do
    !end do
    !close(111)
    !close(112)
    open(111,file='000_dhaa',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')
    open(112,file='000_haa',form='FORMATTED',position='append',action='WRITE',status='UNKNOWN')

    do l=0,atoms%lmaxd
      do io=1,2
        do lp=0,atoms%lmaxd
          do jo=1,2

            if (lp.eq.l) then
              t1=rbas1(:,io,l,1,1)*rbas1(:,jo,lp,1,1)
              xReal=real(t1)
              xImag=aimag(t1)

              call intgr3LinIntp(xReal, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempReal, 1)
              call intgr3LinIntp(xImag, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempImag, 1)

              haa(1,jo,lp,io,l)=El(1, l, 1, 1)*CMPLX(tempReal,tempImag)*sqrt(fpi)
              write(112,*) 1, l, l
              write(112,*) jo, io, real(haa(1,jo,lp,io,l)), aimag(haa(1,jo,lp,io,l))
            else
              haa(1,jo,lp,io,l)=0.0
              write(112,*) 1, l, l
              write(112,*) jo, io, 0.0, 0.0
            end if

            lm=0
            do oqn_l=0,atoms%lmaxd
              do mqn_m=-oqn_l,oqn_l

                lm=lm+1
                write(111,*) lm, lp, l

                t1=rbas1(:,io,l,1,1)*rbas1(:,jo,lp,1,1)*sumVMTs(:,lm,1,1)
                xReal=real(t1)
                xImag=aimag(t1)

                call intgr3LinIntp(xReal, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempReal, 1)
                call intgr3LinIntp(xImag, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempImag, 1)

                dhaa(lm,jo,lp,io,l)=CMPLX(tempReal,tempImag)

                write(111,*) jo, io, real(dhaa(lm,jo,lp,io,l)), aimag(dhaa(lm,jo,lp,io,l))

                if (l.ge.lp.and.lm.ge.2) then
                  t1=rbas1(:,io,l,1,1)*rbas1(:,jo,lp,1,1)*vEff0MT(:,lm,1,1)
                  xReal=real(t1)
                  xImag=aimag(t1)

                  call intgr3LinIntp(xReal, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempReal, 1)
                  call intgr3LinIntp(xImag, atoms%rmsh(:,1), atoms%dx(1), atoms%jri(1), tempImag, 1)

                  haa(lm,jo,lp,io,l)=CMPLX(tempReal,tempImag)
                  haa(lm,io,l,jo,lp)=haa(lm,jo,lp,io,l)

                  write(112,*) lm, lp, l
                  write(112,*) jo, io, real(haa(lm,jo,lp,io,l)), aimag(haa(lm,jo,lp,io,l))

                  if (l.ne.lp) then
                    write(112,*) lm, l, lp
                    write(112,*) io, jo, real(haa(lm,jo,lp,io,l)), aimag(haa(lm,jo,lp,io,l))
                  end if
                end if

              end do
            end do

          end do
        end do
      end do
    end do

    close(111)
    close(112)

  end subroutine udHu

  subroutine ctorsh(atoms,mtc,mtr)

    use m_types, only : t_atoms

    implicit none

    type(t_atoms), intent(in) :: atoms
    complex, intent(in) :: mtc(:,:,:,:)
    real, intent(out) :: mtr(:,:,:,:)

    integer l,m,lm1,lm2

    lm1=0
    do l=0,atoms%lmaxd
      lm2=lm1+2*(l+1)
      do m=-l,-1
        lm1=lm1+1
        lm2=lm2-1
        if (mod(m,2).ne.0) then
          mtr(:,lm1,:,:)=-(aimag(mtc(:,lm1,:,:))+aimag(mtc(:,lm2,:,:)))/sqrt(2.0)
        else
          mtr(:,lm1,:,:)=(aimag(mtc(:,lm2,:,:))-aimag(mtc(:,lm1,:,:)))/sqrt(2.0)
        end if
      end do
      lm1=lm1+1
      lm2=lm2-1
      mtr(:,lm1,:,:)=real(mtc(:,lm1,:,:))
      do m=1,l
        lm1=lm1+1
        lm2=lm2-1
        if (mod(m,2).ne.0) then
          mtr(:,lm1,:,:)=(real(mtc(:,lm1,:,:))-real(mtc(:,lm2,:,:)))/sqrt(2.0)
        else
          mtr(:,lm1,:,:)=(real(mtc(:,lm1,:,:))+real(mtc(:,lm2,:,:)))/sqrt(2.0)
        end if
      end do
    end do
  end subroutine ctorsh

end module m_jpSternhHF
