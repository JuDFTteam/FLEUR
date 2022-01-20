!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Helping module containing all routines for Pulay and Surface contributions to the Sternheimer SCC.
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> Contains routines implementing the Pulay and surface contributions to the Sternheimer equation and is submodule of
!> m_jpsternheimer
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpDens1stVar.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpSternhPulaySurface

    USE m_constants

  implicit none

  contains

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calculates the tlmplm integrals for the matrix element <Ψ|Η_0 - ε|Ψ>_α (recycled subroutine tlmplm from Fleur).
  !>
  !> @details
  !> From tlmplm: Sets up the t(l'm',lm) matrices (written by Weinert 1986) for the displaced muffin-tin α. These matrices are
  !> k-point independent quantities.
  !> See also 7.117 (dissertation CRG)
  !>
  !>  @param[in]  atoms     : Atoms type, see types.f90
  !>  @param[in]  dimens    : Dimension type, see types.f90.
  !>  @param[in]  enpara    : Energy parameter type, see types.f90.
  !>  @param[in]  usdus     : Type containing quantities consisting of the radial solutions, see types.f90.
  !>  @param[in]  input     : Input type, see types.f90.
  !>  @param[out] td        : Tlmplm matrix type for Sternheimer Pulay matrix elements, see types.f90
  !>  @param[in]  jsp       : Spin
  !>  @param[in]  logUnit   : Unit number for juPhon.log.
  !>  @param[in]  rbas1     : Large components of radial solution, its energy derivative and u_LO
  !>  @param[in]  rbas2     : Small components of radial solution, its energy derivative and u_LO
  !>  @param[in]  uuilon    : overlap integral between the radial functions of the integral (multiplied by ulo_der) of a local
  !>                          orbital and the flapw radial function with the same l
  !>  @param[in]  duilon    : overlap integral between the radial functions of the integral of a local orbital and the energy
  !>                          derivative of the flapw radial function with the same l
  !>  @param[in]  ulouilopn : overlap integral between the radial functions of the integral of a local orbital and another local
  !>                          orbital with the same l.
  !>  @param[in]  ilo2p     : mapping array giving the p value for given number of LO and itype
  !>  @param[in]  vr0Sph    : Spherical harmonic coefficients of unperturbed and converged muffin-tin effective potential
  !>                             parsed from Fleur
  !>
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine tlmplm4H0( atoms, enpara, usdus, input, td, loosetdout, jsp, logUnit, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, &
                                                                                                                          & vr0Sph )

    use m_types, only : t_atoms, t_enpara, t_usdus, t_input, t_tlmplm
    use m_intgr, only : intgr3!LinIntp ! TODO: Is this ok?
    use m_gaunt, only : gaunt1
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type Parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_enpara),              intent(in)  :: enpara
    type(t_usdus),               intent(in)  :: usdus
    type(t_input),               intent(in)  :: input
    type(t_tlmplm),              intent(out) :: td

    ! Scalar Parameters
    integer,                     intent(in)  :: jsp
    integer,                     intent(in)  :: logUnit

    ! Array Parameters
    real,                        intent(in)  :: rbas1(:,:,0:,:,:)
    real,                        intent(in)  :: rbas2(:,:,0:,:,:)
    real,                        intent(in)  :: uuilon(:,:)
    real,                        intent(in)  :: duilon(:,:)
    real,                        intent(in)  :: ulouilopn(:,:,:)
    integer,                     intent(in)  :: ilo2p(:, :)
    complex,                     intent(in)  :: vr0Sph(:, :, :)
    COMPLEX, ALLOCATABLE,        INTENT(OUT) :: loosetdout(:, :, :, :)

    ! Scalar Variables
    complex                                  :: cil
    real                                     :: tempReal
    real                                     :: tempImag
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
    real                                     :: tempRbas

    ! Array Variables
    integer,        allocatable              :: indt(:)
    complex,        allocatable              :: dvd(:, :)
    complex,        allocatable              :: dvu(:, :)
    complex,        allocatable              :: uvd(:, :)
    complex,        allocatable              :: uvu(:, : )
    real,           allocatable              :: xReal(:)
    real,           allocatable              :: xImag(:)

    ! Initialization which is taken from eigen.F90 in Fleur
    mlotot = 0
    mlolotot = 0
    DO n = 1, atoms%ntype
      do na = 1, atoms%neq(n)
       mlotot = mlotot + atoms%nlo(n)
       mlolotot = mlolotot + atoms%nlo(n) * (atoms%nlo(n) + 1) / 2
      end do
    ENDDO
    mlot_d = max(mlotot,1)
    mlolot_d = max(mlolotot,1)

    err = 0
    allocate( indt(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2) )
    allocate( dvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( dvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( uvd(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( uvu(0:atoms%lmaxd * (atoms%lmaxd + 3) / 2,(atoms%lmaxd + 1)**2 ) )
    allocate( xReal(atoms%jmtd), xImag(atoms%jmtd) )
    !allocate(td%tuu(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2,atoms%nat,1),stat=err)
    !allocate(td%tud(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2,atoms%nat,1),stat=err)
    !allocate(td%tdd(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2,atoms%nat,1),stat=err)
    !allocate(td%tdu(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2,atoms%nat,1),stat=err)
    allocate(loosetdout(0:(atoms%lmaxd*(atoms%lmaxd+2)* (atoms%lmaxd*(atoms%lmaxd+2)+3))/2,atoms%nat,1,4),stat=err)
    allocate(td%tdulo(0:atoms%lmaxd*(atoms%lmaxd+2),-atoms%llod:atoms%llod,mlot_d,1,1),stat=err) ! TODO: These needed a second spin index.
    allocate(td%tuulo(0:atoms%lmaxd*(atoms%lmaxd+2),-atoms%llod:atoms%llod,mlot_d,1,1),stat=err)
    allocate(td%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,mlolot_d,1,1), stat=err)
    allocate(td%ind(0:atoms%lmaxd*(atoms%lmaxd+2),0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,1),stat=err )
    if (err.ne.0) then
       write (logUnit,'(a)') 'eigen: an error occured during allocation of'
       write (logUnit,'(a)') 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
       CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="eigen")
    end if

    indt(:) = 0
    dvd(:, :) = cmplx(0., 0.)
    dvu(:, :) = cmplx(0., 0.)
    uvd(:, :) = cmplx(0., 0.)
    uvu(:, :) = cmplx(0., 0.)
    xReal(:) = 0.
    xImag(:) = 0.
    !td%tuu(:, :, :) = cmplx(0., 0.)
    !td%tdd(:, :, :) = cmplx(0., 0.)
    !td%tud(:, :, :) = cmplx(0., 0.)
    !td%tdu(:, :, :) = cmplx(0., 0.)
    loosetdout(:, :, :, :) = cmplx(0., 0.)
    td%tdulo(:, :, :, :, 1) = cmplx(0., 0.) ! TODO: These needed a second spin index.
    td%tuulo(:, :, :, :, 1) = cmplx(0., 0.)
    td%tuloulo(:, :, :, :, 1) = cmplx(0., 0.)
    td%ind(:, :, :, :) = -9999

    na = 0
    do n = 1, atoms%ntype
      do ieqat = 1, atoms%neq(n)
      na = na + 1
       !generate the irreducible integrals (u(l'):v(lamda,nu):u(l)) for l' >= l
       do lp = 0, atoms%lmax(n)
          ! NOTE: lp1 is integer
          ! Generation of triangular numbers. NOTE: lp1 is integer
          lp1 = (lp * (lp + 1)) / 2
          do l = 0, lp
             lpl = lp1 + l
             ! loop over non-spherical components of the potential the spherical part (l = 0) will be mixed with the kinetic energy.
            do lamda = 1, atoms%lmax(n)
              lmin = lp - l
              lmx = lp + l
              ! We only have a contribution according to the Gaunt selection rules. The triangular condition can be derived from an
              ! inequalities system so that we can set up the triangular conditions for every index of the Gaunt coefficients.
              ! Furthermore l'+l+lamda should be even
              if ((mod(lamda + lmx, 2) .eq. 1) .or. (lamda.lt.lmin) .or. (lamda.gt.lmx)) then
                ! No contribution
                do m = -lamda, lamda ! changed
                  lm = lamda * (lamda + 1) + m + 1 ! changed
                  uvu(lpl, lm) = 0.0
                  dvd(lpl, lm) = 0.0
                  uvd(lpl, lm) = 0.0
                  dvu(lpl, lm) = 0.0
                end do
              else
                ! Gaunt coefficient is not zero for this branch. As we have complex expansion coefficients, we have to calculate
                ! 2 real integrals. This was done for security reasons and to have a similiar method to tlmplm4V which can be
                ! tested.
                ! todo In a second step we can use the lattice harmonic coefficients to speed up the routine. In principle we still
                ! have the Fleur symmetry for the unperturbed quantities but should benchmark the alternative routine against this
                ! algorithm.
                xReal = 0.
                xImag = 0.
                do m = -lamda, lamda
                  lm = lamda * (lamda + 1) + m + 1
                  ! Calculate the integral <u|V|u>
                  do i = 1, atoms%jri(n)
                    tempRbas = ( rbas1(i, 1, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                    xReal(i) =  tempRbas * real(vr0SpH(i, lm, na))
                    xImag(i) =  tempRbas * aimag(vr0SpH(i, lm, na))
                  end do
                  call intgr3(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal) ! TODO: Is this ok?
                  call intgr3(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag) ! TODO: Is this ok?
                  uvu(lpl, lm) = cmplx(tempReal, tempImag)
                  ! Calculate the integral <uDot|V|u>
                  do i = 1, atoms%jri(n)
                    tempRbas = ( rbas1(i, 2, lp, n, 1) * rbas1(i, 1, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 1, l, n, 1) )
                    xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                    xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                  end do
                  call intgr3(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal) ! TODO: Is this ok?
                  call intgr3(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag) ! TODO: Is this ok?
                  dvu(lpl, lm) = cmplx(tempReal, tempImag)
                  ! Calculate the integral <u|V|uDot>
                  do i = 1,atoms%jri(n)
                    tempRbas = (rbas1(i, 1, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 1, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                    xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                    xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                  end do
                  call intgr3(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal) ! TODO: Is this ok?
                  call intgr3(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag) ! TODO: Is this ok?
                  uvd(lpl, lm) = cmplx(tempReal, tempImag)
                  ! Calculte the integral <uDot|V|uDot>
                  do i = 1,atoms%jri(n)
                    tempRbas = (rbas1(i, 2, lp, n, 1) * rbas1(i, 2, l, n, 1) + rbas2(i, 2, lp, n, 1) * rbas2(i, 2, l, n, 1) )
                    xReal(i) = tempRbas * real(vr0SpH(i, lm, na))
                    xImag(i) = tempRbas * aimag(vr0SpH(i, lm, na))
                  end do
                  call intgr3(xReal, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempReal) ! TODO: Is this ok?
                  call intgr3(xImag, atoms%rmsh(1,n), atoms%dx(n), atoms%jri(n), tempImag) ! TODO: Is this ok?
                  dvd(lpl, lm) = cmplx(tempReal, tempImag)
                end do ! m
              end if ! Gaunt selection rules
            end do ! lh
          end do ! l
        end do ! lp
        ! generate the various t(l'm',lm) matrices for l'm'.ge.lm
        indt(:) = 0
        ! loop over l'm'
        do lp = 0, atoms%lmax(n)
          ! NOTE: lp1 is integer
          lp1 = (lp * (lp + 1)) / 2
          do mp = -lp, lp
            lmp = lp * (lp + 1) + mp
            ! NOTE: lmpl is integer
            lmpl = (lmp * (lmp + 1)) / 2
            ! Loop over the non-spherical (l \= 0) components of the potential
            do lamda = 1, atoms%lmax(n)
              lmin0 = abs(lp - lamda)
              ! l' - lambda > l' and -l' + lambda > l' leads to 0 <= lamda <= 2 l'
              IF (lmin0.GT.lp) cycle
              ! lmxx is in principle lp. To ensure, that the oddness of lambda is always compensated (here strictly speaking only at
              ! lmxx) by l, the modulo operation is subtracted. It is subtracted to not exceed lmxx. The lmxx boundary was given by
              ! Gaunt selection rules, actually reading l' + lambda but as we have l <= l' lmax is only l'. We group lambda and the
              ! modulo operation in lmxx + l' + lambda, so that we only have to ensure l + l' is even.
              lmxx = lp - mod(lamda, 2)
              do mu = -lamda, lamda
                !collection index of lamda and mu!
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
                  lm = l * (l+1) + m
                  ! lm is always <= lmp!
                  IF (lm.GT.lmp) CYCLE
                  ! index for uvu, dvd, uvd, dvu to access only those relevant after the Gaunt selection rules
                  lpl = lp1 + l
                  ! index similiar to lpl but with ms in it (after the selection Gaunt selection rules)
                  lmplm = lmpl + lm
                  ! After having found out that the Gaunt coefficient is not zero we multiply it.
                  cil = gaunt1(lp, lamda, l, mp, mu, m, atoms%lmaxd)

                  !td%tuu(lmplm, na, jsp) = td%tuu(lmplm, na, jsp) + cil * uvu(lpl, lmsph)
                  !td%tdd(lmplm, na, jsp) = td%tdd(lmplm, na, jsp) + cil * dvd(lpl, lmsph)
                  !td%tud(lmplm, na, jsp) = td%tud(lmplm, na, jsp) + cil * uvd(lpl, lmsph)
                  !td%tdu(lmplm, na, jsp) = td%tdu(lmplm, na, jsp) + cil * dvu(lpl, lmsph)
                  loosetdout(lmplm, na, jsp, 1) = loosetdout(lmplm, na, jsp, 1) + cil * uvu(lpl, lmsph)
                  loosetdout(lmplm, na, jsp, 2) = loosetdout(lmplm, na, jsp, 2) + cil * uvd(lpl, lmsph)
                  loosetdout(lmplm, na, jsp, 3) = loosetdout(lmplm, na, jsp, 3) + cil * dvu(lpl, lmsph)
                  loosetdout(lmplm, na, jsp, 4) = loosetdout(lmplm, na, jsp, 4) + cil * dvd(lpl, lmsph)
                  ! Logical matrix where there are non-vanishing matrix entries in the tlmplm matrices
                  indt(lmplm) = 1
                end do ! l
              end do ! mem
            end do ! lh
          end do ! mp
        end do ! lp

        ! set up mapping array
        do lp = 0, atoms%lmax(n)
          do mp = -lp, lp
            lmp = lp * (lp + 1) + mp
            do l = 0, atoms%lmax(n)
              do m = -l, l
                lm = l * (l + 1) + m
                if (lmp.ge.lm) then
                  lmplm = (lmp * (lmp + 1)) / 2 + lm
                  if (indt(lmplm).NE.0) THEN
                    td%ind(lmp, lm, na, jsp) = lmplm
                  else
                    td%ind(lmp, lm, na, jsp) = -9999
                  end if
                else
                  ! As we have lm > lmp here (this was not calculated within the routine), we have to transpose, i.e. interchanging
                  ! lmp and lm for the packed storage index.
                  lmplm = (lm* (lm + 1)) / 2 + lmp
                  if (indt(lmplm).NE.0) THEN
                    td%ind(lmp, lm, na, jsp) = -lmplm
                  else
                    td%ind(lmp, lm, na, jsp) = -9999
                  end if
                end if
              end do ! m
            end do ! l
          end do ! mp
        end do ! lp

        ! set up the t-matrices for the local orbitals, if there are any
        if (.false.) then
          if (atoms%nlo(n).ge.1) then
            call tlo4HS0(atoms, enpara, usdus, input, td, input%jspins, jsp, n, na, vr0SpH(:, :, na), rbas1, rbas2, uuilon, &
                                                                                                         & duilon, ulouilopn, ilo2p)
          end if
        end if

      end do !ieqat
    end do ! n

  end subroutine tlmplm4H0

  !>--------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calculates the tlo integrals for the matrix element <Ψ|H - ε|Ψ>_α (recycled subroutine tlo from Fleur).
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
  !> @todo build in El
  !>--------------------------------------------------------------------------------------------------------------------------------
  subroutine tlo4HS0(atoms, enpara, usdus, input, tlmplm, jspin, jsp, ntyp, na, vr, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p )

    !use m_intgr, only : intgr3
    use m_intgr, only : intgr3!LinIntp ! TODO: Is this ok?
    use m_gaunt, only : gaunt1
    use m_types
    use m_juDFT_stop, only : juDFT_error
    implicit none

    type(t_atoms),         intent(in)    :: atoms
    type(t_enpara),        intent(in)    :: enpara
    type(t_usdus),         intent(in)    :: usdus
    type(t_input),         intent(in)    :: input
    type(t_tlmplm),        intent(inout) :: tlmplm
    !     ..
    !     .. Scalar Arguments ..
    integer,               intent (in)   :: jspin,jsp,ntyp ,na
    !     ..
    !     .. Array Arguments ..
    complex,               intent (in)   :: vr(atoms%jmtd,(atoms%lmaxd + 1)**2)
    real,                  intent(in)    :: rbas1(:,:,0:,:,:)
    real,                  intent(in)    :: rbas2(:,:,0:,:,:)
    real,                  intent (in)   :: uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
    real,                  intent (in)   :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
    integer,               intent (in)   :: ilo2p(:, :)
    !     ..
    !     .. Local Scalars ..
    real                                 :: tempReal, tempImag  !changed
    complex                              :: cil
    integer                              :: lmpp
    integer                              :: i,l,lh,lm ,lmin,lmp,lo,lop,loplo,lp,lpmax,lpmax0,lpmin,lpmin0,lpp ,mem,mp,mpp,m,lmx,mlo,mlolo
    integer                              :: mytype
    integer                              :: myatom
    integer                              :: myeqat
    !     ..
    !     .. Local Arrays ..
    real                                 :: xReal(atoms%jmtd), xImag(atoms%jmtd)
    complex                              :: ulovulo(atoms%nlod*(atoms%nlod+1)/2,(atoms%lmaxd + 1)**2)
    complex                              :: uvulo(atoms%nlod,0:atoms%lmaxd,(atoms%lmaxd + 1)**2),dvulo(atoms%nlod,0:atoms%lmaxd,(atoms%lmaxd + 1)**2)
    integer, allocatable                 :: nlo_atom(:)

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
          do lpp = 1, atoms%lmax(ntyp)
             if ((mod(l+lp+lpp,2).eq.1) .or. (lpp.lt.lmin) .or.&
                  (lpp.GT.lmx)) then
                  do mpp = - lpp, lpp
                    lm = lpp * (lpp + 1) + 1 + mpp ! changed
                    uvulo(lo,lp,lm) = 0.0
                    dvulo(lo,lp,lm) = 0.0
                  end do
             else
               do mpp = - lpp, lpp
               lmpp = lpp * (lpp + 1) + 1 + mpp
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,1,lp,ntyp,1)*rbas1(i, ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+ rbas2(i, 1, lp, ntyp, 1) * rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp))
                   xImag(i) = (rbas1(i,1,lp,ntyp,1)*rbas1(i, ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+ rbas2(i, 1, lp, ntyp, 1) * rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp))
                end do
                call intgr3(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal) ! TODO: Is this ok?
                call intgr3(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag) ! TODO: Is this ok?
                uvulo(lo,lp,lmpp) = cmplx(tempReal, tempImag)
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,2,lp,ntyp,1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp),ntyp,1)+ rbas2(i,2,lp,ntyp,1)*rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp))
                   xImag(i) = (rbas1(i,2,lp,ntyp,1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp),ntyp,1)+ rbas2(i,2,lp,ntyp,1)*rbas2(i,ilo2p(lo, ntyp),atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp))
                end do
                call intgr3(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal) ! TODO: Is this ok?
                call intgr3(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag) ! TODO: Is this ok?
                dvulo(lo,lp,lmpp) = cmplx(tempReal, tempImag)
               end do
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
          do lpp = 1, atoms%lmax(ntyp)
             lmin = ABS(lp - l)
             lmx = lp + l
             if ((mod(l+lp+lpp,2).eq.1).or.(lpp.lt.lmin).or.(lpp.gt.lmx)) then
               do mpp = -lpp, lpp
                lmpp = lpp * (lpp + 1) + 1 + mpp
                ulovulo(loplo,lmpp) = 0.0
                end do
             else
               do mpp = -lpp, lpp
               lmpp = lpp * (lpp + 1) + 1 + mpp
                do i = 1,atoms%jri(ntyp)
                   xReal(i) = (rbas1(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+rbas2(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas2(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1))*real(vr(i,lmpp))
                   xImag(i) = (rbas1(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas1(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1)+rbas2(i,ilo2p(lop, ntyp), atoms%llo(lop, ntyp), ntyp, 1)*rbas2(i,ilo2p(lo, ntyp), atoms%llo(lo, ntyp), ntyp, 1))*aimag(vr(i,lmpp))
                end do
                call intgr3(xReal,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempReal) ! TODO: Is this ok?
                call intgr3(xImag,atoms%rmsh(:,ntyp),atoms%dx(ntyp),atoms%jri(ntyp),tempImag) ! TODO: Is this ok?
                ulovulo(loplo, lmpp) = cmplx(tempReal, tempImag)
                end do
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
          do lpp = 1, atoms%lmax(ntyp)
             lpmin0 = ABS(l-lpp)
             lpmax0 = l + lpp
             !--->          check that lpmax is smaller than the max l of the
             !--->          wavefunction expansion at this atom
             lpmax = MIN(lpmax0,atoms%lmax(ntyp))
             !--->          make sure that l + l'' + lpmax is even
             lpmax = lpmax - MOD(l+lpp+lpmax,2)
             do mpp = -lpp, lpp
             lmpp = lpp * (lpp + 1) + 1 + mpp
                mp = m + mpp
                lpmin = MAX(lpmin0,ABS(mp))
                !--->             make sure that l + l'' + lpmin is even
                lpmin = lpmin + MOD(ABS(lpmax-lpmin),2)
                !--->             loop over l'
                do lp = lpmin,lpmax,2
                   lmp = lp* (lp+1) + mp
                   cil = gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                   tlmplm%tuulo(lmp,m,lo+mlo,jsp,1) = & ! TODO: These needed a second spin index.
                        tlmplm%tuulo(lmp,m,lo+mlo,jsp,1) + cil*uvulo(lo,lp,lmpp)
                   tlmplm%tdulo(lmp,m,lo+mlo,jsp,1) = &
                        tlmplm%tdulo(lmp,m,lo+mlo,jsp,1) + cil*dvulo(lo,lp,lmpp)
                end do
             end do
          end do
       end do
    end do
    !---> generate the t-matrix including two local orbitals for lo' >= lo
    !---> loop over lo'
    mlolo=dot_product(nlo_atom(:na-1),nlo_atom(:na-1)+1)/2
    do lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       do mp = -lp,lp
          !--->       loop over the lattice harmonics
          do lpp = 1,atoms%lmax(ntyp)
             do mpp = - lpp, lpp
               lmpp = lpp * (lpp + 1) + 1 + mpp
                m = mp - mpp
                !--->             loop over lo
                do lo = 1,lop
                   l = atoms%llo(lo,ntyp)
                   loplo = ((lop-1)*lop)/2 + lo
                   if ((abs(l-lpp).le.lp) .and. (lp.le. (l+lpp)) .and.&
                        (mod(l+lp+lpp,2).eq.0) .and. (abs(m).le.l)) then
                      cil = gaunt1(lp,lpp,l,mp,mpp,m,atoms%lmaxd)
                      tlmplm%tuloulo(mp,m,loplo+mlolo,jsp,1) = tlmplm%tuloulo(mp,m,loplo+mlolo,jsp,1) + cil*ulovulo(loplo,lmpp) ! TODO: These needed a second spin index.
                   end if
                end do
             end do
          end do
       end do
    end do
    !---> add the diagonal terms from the muffin-tin hamiltonian. these
    !---> terms have to be made hermitian. if second variation is switched
    !---> on, the t-matrices contain only the contributions from the
    !---> non-spherical hamiltonian.
    !todo here ello0 is already used, in the future switch to El!
    do lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       do m = -l,l
          lm = l* (l+1) + m
          tlmplm%tuulo(lm,m,lo+mlo,jsp,1) = tlmplm%tuulo(lm,m,lo+mlo,jsp,1) + 0.5 * usdus%uulon(lo,ntyp,jspin) *& ! TODO: These needed a second spin index.
               ( enpara%el0(l,ntyp,jspin)+enpara%ello0(lo,ntyp,jspin) )
          tlmplm%tdulo(lm,m,lo+mlo,jsp,1) = tlmplm%tdulo(lm,m,lo+mlo,jsp,1) + 0.5 * usdus%dulon(lo,ntyp,jspin) *& ! TODO: These needed a second spin index.
               ( enpara%el0(l,ntyp,jspin)+enpara%ello0(lo,ntyp,jspin) ) + 0.5 * usdus%uulon(lo,ntyp,jspin)
          if (atoms%ulo_der(lo,ntyp).GE.1) THEN
             tlmplm%tuulo(lm,m,lo+mlo,jsp,1) = tlmplm%tuulo(lm,m,lo+mlo,jsp,1) + 0.5 * uuilon(lo,ntyp) ! TODO: These needed a second spin index.
             tlmplm%tdulo(lm,m,lo+mlo,jsp,1) = tlmplm%tdulo(lm,m,lo+mlo,jsp,1) + 0.5 * duilon(lo,ntyp)
          endif
          !+apw_lo
          if (atoms%l_dulo(lo,ntyp)) THEN
             tlmplm%tuulo(lm,m,lo+mlo,jsp,1) = tlmplm%tuulo(lm,m,lo+mlo,jsp,1) + 0.5 ! TODO: These needed a second spin index.
             tlmplm%tdulo(lm,m,lo+mlo,jsp,1) = 0.0
          endif
          !+apw_lo
        end do
    end do
    do lop = 1,atoms%nlo(ntyp)
       lp = atoms%llo(lop,ntyp)
       do lo = atoms%lo1l(lp,ntyp),lop
          loplo = ((lop-1)*lop)/2 + lo
          do m = -lp,lp
             tlmplm%tuloulo(m,m,loplo+mlolo,jsp,1) = tlmplm%tuloulo(m,m,loplo+mlolo,jsp,1) + 0.5* (enpara%ello0(lop,ntyp,jspin)+& ! TODO: These needed a second spin index.
                  enpara%ello0(lo,ntyp,jspin))* usdus%uloulopn(lop,lo,ntyp,jspin) + 0.5* (ulouilopn(lop,lo,ntyp) +&
                  ulouilopn(lo,lop,ntyp))
          end do
       end do
    end do

  end subroutine tlo4HS0

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the Pulay incomplete basis correction in Sternheimer, namely <Psi_k+q | (H - eps) | Psi_k >
  !>
  !> @details
  !> See 7.115 - 7.122 (dissertation CRG)
  !>
  !> @details
  !> @param[in]   atoms    : Atoms type, see types.f90
  !> @param[in]   usdus    : Type containing quantities consisting of the radial solutions, see types.f90.
  !> @param[in]   tlmplm   : Tlmplm matrix type for Sternheimer Pulay matrix elements, see types.f90
  !> @param[in]   ikpt     : Index of k-point in k-point set
  !> @param[in]   ikpq     : Index of k + q (k' backfolded) in k-point set
  !> @param[in]   itype    : Index of atom type to which displaced atom belongs
  !> @param[in]   iatom    : Index of displaced atom
  !> @param[in]   ne       : Number of eigenvalues per k-point.
  !> @param[in]   nobd     : Number of occupied bands per k-point and spin
  !> @param[in]   El       : Contains LAPW and LO energy parameters.
  !> @param[in]   mCoefBp  : Matching coefficient as defined in Equation 7.112a and 7.112b (dissertation CRG) or canonical matching
  !>                         coefficient of common muffin-tin wave function at k
  !> @param[in]   mCoefKb  : Matching coefficient as defined in Equation 7.91c (dissertation CRG) or canonical matching coefficient
  !>                         of common muffin-tin wave function at k
  !> @param[in]   nRadFun  : Number of radial functions per orbital quantum number l and atom type.
  !> @param[in]   iloTable : Number of local orbital if orbital quantum number l, atom type and index p > 2 from nRadFun is given.
  !> @param[in]   nlo_atom : Contains information about number of LOs at every atom
  !> @param[out]  s        : Overlap of the Pulay incomplete basis correction in Sternheimer
  !> @param[out]  h        : Hamiltonian of the Pulay incomplete basis correction in Sternheimer
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcHS0MT( atoms, usdus, tlmplm, loosetdin, ikpt, ikpq, itype, iatom, ne, nobd, El, mCoefBp, mCoefKb, nRadFun, iloTable, &
      & nlo_atom, s, h, idir, nvk, nvkq, almkg, blmkg, almkgq, blmkgq, dalmkg, dblmkg, dalmkgq, dblmkgq, mat_elH, mat_elS )

    use m_types

    implicit none

    ! Type Parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_usdus),               intent(in)  :: usdus
    type(t_tlmplm),              intent(in)  :: tlmplm

    ! Scalar Parameters
    integer,                     intent(in)  :: ikpt
    integer,                     intent(in)  :: ikpq
    integer,                     intent(in)  :: itype
    integer,                     intent(in)  :: iatom

    ! Array Parameters
    integer,                     intent(in)  :: ne(:)
    integer,                     intent(in)  :: nobd(:, :)
    real,                        intent(in)  :: El(:, 0:, :, :)
    complex,                     intent(in)  :: mCoefBp(:, :)
    complex,                     intent(in)  :: mCoefKb(:, :)
    integer,                     intent(in)  :: nRadFun(0:, :)
    integer,                     intent(in)  :: iloTable(:, 0:, :)
    integer,                     intent(in)  :: nlo_atom(:)
    complex,                    intent(out) :: s(:, :)
    complex,                    intent(out) :: h(:,:)
    integer, optional,           intent(in)  :: idir
    integer, optional,           intent(in)  :: nvk, nvkq
    complex, optional,           intent(in)  :: almkg(:,:), blmkg(:,:), almkgq(:,:), blmkgq(:,:)
    complex, optional,           intent(in)  :: dalmkg(:,:), dblmkg(:,:), dalmkgq(:,:), dblmkgq(:,:)
    complex, optional,           intent(inout) :: mat_elH(:,:), mat_elS(:,:)
    COMPLEX, INTENT(IN) :: loosetdin(:, :, :, :)

    ! Local Variables
    integer                                  :: mlo
    integer                                  :: mlolo
    integer                                  :: lmlo, lm
    integer                                  :: oqn_l
    integer                                  :: mqn_m
    integer                                  :: nband
    integer                                  :: pband
    integer                                  :: coBsh
    integer                                  :: loBra
    integer                                  :: coKsh
    integer                                  :: loKet
    integer                                  :: lmloB
    integer                                  :: lmB
    integer                                  :: mB
    integer                                  :: lmloK
    integer                                  :: lmK
    integer                                  :: lK
    integer                                  :: mK
    integer                                  :: ind, iG, iGq
    integer                                  :: indN
    integer                                  :: loBraKet
    integer                                  :: lB
    complex                                 :: hpnEl, hEl, hElq
    complex                                 :: spnEl, sEl, sElq
    complex                                 :: spnElo
    complex                                 :: hpnAdd, hAdd, hAddq
    complex                                  :: utu
    complex                                  :: dtu
    complex                                  :: utd
    complex                                  :: dtd, utu2, utd2, dtu2, dtd2
    complex                                  :: utulo
    complex                                  :: dtulo
    complex                                  :: ulotu
    complex                                  :: ulotd
    complex                                  :: ulotulo
    complex                                  :: spnAdd

    complex,        allocatable              :: ax(:)
    complex,        allocatable              :: bx(:)
    complex,        allocatable              :: cx(:, :)

    complex,        allocatable              :: ax2(:), ax3(:), ax4(:)
    complex,        allocatable              :: bx2(:), bx3(:), bx4(:)
    complex,        allocatable              :: hSumG(:,:), hSumGq(:,:), sSumG(:,:), sSumGq(:,:), hsurf(:,:)

    allocate(ax(nobd(ikpt, 1)), bx(nobd(ikpt, 1)))
    ax = cmplx(0., 0.)
    bx = cmplx(0., 0.)
    !allocate(ax4(nobd(ikpt, 1)), bx4(nobd(ikpt, 1)))
    !ax4 = cmplx(0., 0.)
    !bx4 = cmplx(0., 0.)

    s(:, :) = cmplx(0., 0.)
    h(:, :) = cmplx(0., 0.)

    !allocate(hsurf,mold=h)
    !hsurf(:, :) = cmplx(0., 0.)

    if (.FALSE..and.present(almkg)) then
      allocate(ax2(nvk), bx2(nvk),ax3(nvk), bx3(nvk))
      ax2(:) = cmplx(0., 0.)
      bx2(:) = cmplx(0., 0.)
      ax3(:) = cmplx(0., 0.)
      bx3(:) = cmplx(0., 0.)
      allocate(hSumG(nvkq,nvk),hSumGq(nvkq,nvk),sSumG(nvkq,nvk),sSumGq(nvkq,nvk))
      hSumG(:,:) = cmplx(0., 0.)
      hSumGq(:,:) = cmplx(0., 0.)
      sSumG(:,:) = cmplx(0., 0.)
      sSumGq(:,:) = cmplx(0., 0.)
    end if

    lmlo = 0
    lm = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        ! go to acof of next lm
        lmlo = lmlo + 1
        lm = lm +1

        if (oqn_l.eq.0) then
          utu2=0.5*usdus%us(0,1,1)*usdus%dus(0,1,1)*atoms%rmt(1)**2
          utd2=0.5*usdus%us(0,1,1)*usdus%duds(0,1,1)*atoms%rmt(1)**2
          dtu2=0.5*usdus%uds(0,1,1)*usdus%dus(0,1,1)*atoms%rmt(1)**2
          dtd2=0.5*usdus%uds(0,1,1)*usdus%duds(0,1,1)*atoms%rmt(1)**2
        end if

        !Attention: the LO spherical part of the Hamiltonian matrix element is taken care of within tlo routine. Possibly we can
        ! improve performance here if we add the spherical part of the Hamiltonian to the tlmplm integrals to the tlmplm integrals
        ! variable, which is, at the moment, only done temporarily within hlomat routine of Kurz and then subracted at the end to
        ! be consisten with FLeur routines. One could try whether it is sufficient to only calculate the overlapp matrix here and
        ! do the rest outside because the tlmplm routines are only called once this type of routines mor than once.

        do nBand = 1, nobd(ikpt, 1)
          do pBand = 1, ne(ikpq)

            spnEl = 0
            hpnAdd = 0

            ! If index lmlo then acof is meant, if lmlo + 1 then bcof is meant
            spnEl = conjg(mCoefBp(pBand, lmlo)) * mCoefKb(nBand, lmlo) + conjg(mCoefBp(pBand, lmlo + 1))&
              & * mCoefKb(nBand, lmlo + 1) * usdus%ddn(oqn_l, itype, 1)

            hpnAdd = conjg(mCoefBp(pband, lmlo)) * mCoefKb(nBand, lmlo + 1)
            !CRGfix
            !hpnAdd = 0.5 * (conjg(mCoefBp(pband, lmlo)) * mCoefKb(nBand, lmlo + 1) + &
            !              & conjg(mCoefBp(pband, lmlo + 1)) * mCoefKb(nBand, lmlo ))

            spnElo = 0
            ! Initialize shift so that next increment shifts to ccof, with coKsh=1 we are at bcof, incrementing gives first ccof
            coKsh = 1
            if (.false.) then
              do loKet = 3, nRadFun(oqn_l, itype)
              coKsh = coKsh + 1
                spnElo = spnElo + ( conjg(mCoefBp(pBand, lmlo))  * usdus%uulon(iloTable(loKet, oqn_l, itype), itype, 1) + &
                  & conjg(mCoefBp(pBand, lmlo + 1))  * usdus%dulon(iloTable(loKet, oqn_l, itype), itype, 1) ) * &
                  & mCoefKb(nBand, lmlo + coKsh)
              end do
            end if

            spnAdd = 0
            ! Initialize shift so that next increment shifts to ccof, with coBsh=1 we are at bcof, incrementing gives first ccof
            coBsh = 1
            if (.false.) then
              do loBra = 3, nRadFun(oqn_l, itype)
                coBsh = coBsh + 1
                !todo LO: where is hpnElo?
                spnElo = spnElo + conjg(mCoefBp(pBand, lmlo + coBsh)) * ( mCoefKb(nBand, lmlo)  &
                & * usdus%uulon(iloTable(loBra, oqn_l, itype), itype, 1) + mCoefKb(nBand, lmlo + 1) &
                & * usdus%dulon(iloTable(loBra, oqn_l, itype), itype, 1) )
                coKsh = 1
                do loKet = 3, nRadFun(oqn_l, itype)
                  coKsh = coKsh + 1
                  ! todo LO: we have to adopt the loop structure to the structure below, otherwise not all acofs are multiplied to
                  ! the ccofs
                  spnAdd = spnAdd + conjg(mCoefBp(pBand, lmlo + coBsh))  * &
                    &usdus%uloulopn(iloTable(loBra, oqn_l, itype), iloTable(loKet, oqn_l, itype), itype, 1) * &
                    &mCoefKb(nBand, lmlo + coKsh)
                end do ! loKet
              end do ! loBra
            end if

            hpnEl = El(1, oqn_l, itype, 1) * spnEl
            s(pBand, nBand) = s(pBand, nBand) + spnEl + spnElo + spnAdd
            h(pBand, nBand) = h(pBand, nBand) + hpnEl + hpnAdd

            if (oqn_l.eq.0) then
              !hsurf(pBand, nBand) = hsurf(pBand, nBand) + conjg(mCoefBp(pBand, lmlo))   * utu2 * mCoefKb(nBand, lmlo) &
              !                                        & + conjg(mCoefBp(pBand, lmlo))   * utd2 * mCoefKb(nBand, lmlo+1) &
              !                                        & + conjg(mCoefBp(pBand, lmlo+1)) * dtu2 * mCoefKb(nBand, lmlo) &
              !                                        & + conjg(mCoefBp(pBand, lmlo+1)) * dtd2 * mCoefKb(nBand, lmlo+1)
              if (.FALSE.) then
                h(pBand, nBand) = h(pBand, nBand) + conjg(mCoefBp(pBand, lmlo))   * utu2 * mCoefKb(nBand, lmlo) &
                                                & + conjg(mCoefBp(pBand, lmlo))   * utd2 * mCoefKb(nBand, lmlo+1) &
                                                & + conjg(mCoefBp(pBand, lmlo+1)) * dtu2 * mCoefKb(nBand, lmlo) &
                                                & + conjg(mCoefBp(pBand, lmlo+1)) * dtd2 * mCoefKb(nBand, lmlo+1)
              end if
            end if

          end do ! pband
        end do ! nband

        if (.FALSE..and.present(almkg)) then
          do iG=1, nvk
            do iGq=1, nvkq
              sEl = 0
              hAdd = 0
              sElq = 0
              hAddq = 0

              sEl = conjg(almkgq(lm, iGq)) * dalmkg(lm, iG) + conjg(blmkgq(lm, iGq)) * dblmkg(lm, iG) * usdus%ddn(oqn_l, itype, 1)
              hAdd = conjg(almkgq(lm, iGq)) * dblmkg(lm, iG)

              sElq = conjg(dalmkgq(lm, iGq)) * almkg(lm, iG) + conjg(dblmkgq(lm, iGq)) * blmkg(lm, iG) * usdus%ddn(oqn_l, itype, 1)
              hAddq = conjg(dalmkgq(lm, iGq)) * blmkg(lm, iG)

              hEl = El(1, oqn_l, itype, 1) * sEl
              hElq = El(1, oqn_l, itype, 1) * sElq
              sSumG(iGq, iG) = sSumG(iGq, iG) + sEl
              hSumG(iGq, iG) = hSumG(iGq, iG) + hEl + hAdd

              sSumGq(iGq, iG) = sSumGq(iGq, iG) + sElq
              hSumGq(iGq, iG) = hSumGq(iGq, iG) + hElq + hAddq

              if (oqn_l.eq.0) then
                hSumG(iGq, iG) = hSumG(iGq, iG) + conjg(almkgq(lm, iGq)) * utu2 * dalmkg(lm, iG) &
                                              & + conjg(almkgq(lm, iGq)) * utd2 * dblmkg(lm, iG) &
                                              & + conjg(blmkgq(lm, iGq)) * dtu2 * dalmkg(lm, iG) &
                                              & + conjg(blmkgq(lm, iGq)) * dtd2 * dblmkg(lm, iG)
                hSumGq(iGq, iG) = hSumGq(iGq, iG) + conjg(dalmkgq(lm, iGq)) * utu2 * almkg(lm, iG) &
                                                & + conjg(dalmkgq(lm, iGq)) * utd2 * blmkg(lm, iG) &
                                                & + conjg(dblmkgq(lm, iGq)) * dtu2 * almkg(lm, iG) &
                                                & + conjg(dblmkgq(lm, iGq)) * dtd2 * blmkg(lm, iG)
              end if
            end do
          end do
        end if

        ! Add bcof plus (if any) the radial LO functions for this l, so that we are 1 before the last coBsh
        lmlo = lmlo + coBsh
      end do ! mqn_m
    end do ! oqn_l


    ! Determine shift of LO-type tlmplm arrays for current atom type
    mlo = 0
    mlolo = 0
    if ( atoms%nlo(itype) > 0 ) then
      allocate(cx(nobd(ikpt, 1), 2: maxval(nRadFun) - 1))
      mlo = sum( nlo_atom(:iatom - 1) )
      mlolo = dot_product( nlo_atom(:iatom - 1), nlo_atom(:iatom - 1) + 1 ) / 2
    end if

    ! Calculate non-spherical part of the Hamiltonian
    lmloB = 1
    lmB = -1
    do lB = 0, atoms%lmax(itype)
      do mB = -lB,lB
        lmB = lmB + 1
        ax = cmplx(0., 0.)
        bx = cmplx(0., 0.)
        cx = cmplx(0., 0.)
        !ax4 = cmplx(0., 0.)
        !bx4 = cmplx(0., 0.)

        if (.FALSE.) then
          ax2 = cmplx(0., 0.)
          bx2 = cmplx(0., 0.)
          ax3 = cmplx(0., 0.)
          bx3 = cmplx(0., 0.)
        end if

        lmloK = 1
        lmK = -1
        do lK = 0,atoms%lmax(itype)
          do mK = -lK,lK
            lmK = lmK + 1

            ind = tlmplm%ind(lmB, lmK, iatom, 1) !todo LO:we should allow contributions with same l in ind sothat los work correctly
            if (ind /= -9999) THEN

              !todo REVIEW ALGORITHM
              if (ind >= 0) THEN !todo these are not conjg in kurz hlomat
                !utu = iu**(lK - lB) * tlmplm%tuu(ind, iatom, 1)
                !dtu = iu**(lK - lB) * tlmplm%tdu(ind, iatom, 1)
                !utd = iu**(lK - lB) * tlmplm%tud(ind, iatom, 1)
                !dtd = iu**(lK - lB) * tlmplm%tdd(ind, iatom, 1)
                !utu = tlmplm%tuu(ind, iatom, 1)
                !dtu = tlmplm%tdu(ind, iatom, 1)
                !utd = tlmplm%tud(ind, iatom, 1)
                !dtd = tlmplm%tdd(ind, iatom, 1)
                utu = loosetdin(ind, iatom, 1, 1)
                utd = loosetdin(ind, iatom, 1, 2)
                dtu = loosetdin(ind, iatom, 1, 3)
                dtd = loosetdin(ind, iatom, 1, 4)

                utu2=0.5*usdus%us(lB,1,1)*usdus%dus(lB,1,1)*atoms%rmt(1)**2
                utd2=0.5*usdus%us(lB,1,1)*usdus%duds(lB,1,1)*atoms%rmt(1)**2
                dtu2=0.5*usdus%uds(lB,1,1)*usdus%dus(lB,1,1)*atoms%rmt(1)**2
                dtd2=0.5*usdus%uds(lB,1,1)*usdus%duds(lB,1,1)*atoms%rmt(1)**2

                if (.FALSE..and.(lB.eq.lK.and.mB.eq.mK)) then
                  utu=utu+utu2
                  utd=utd+utd2
                  dtu=dtu+dtu2
                  dtd=dtd+dtd2
                end if

                !if (.FALSE.) then
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 1
                !  write(111,*) real(utu), aimag(utu)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 1
                !  write(111,*) real(dtu), aimag(dtu)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 2
                !  write(111,*) real(utd), aimag(utd)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 2
                !  write(111,*) real(dtd), aimag(dtd)
                !end if

                utu=utu*(ImagUnit**(lK - lB))
                dtu=dtu*(ImagUnit**(lK - lB))
                utd=utd*(ImagUnit**(lK - lB))
                dtd=dtd*(ImagUnit**(lK - lB))
                utu2=utu2*(ImagUnit**(lK - lB))
                dtu2=dtu2*(ImagUnit**(lK - lB))
                utd2=utd2*(ImagUnit**(lK - lB))
                dtd2=dtd2*(ImagUnit**(lK - lB))
              else !todo these are conjg in kurz hlomat
                indN = -ind
                !utu = iu**(lK - lB) * conjg(tlmplm%tuu(indN, iatom, 1))
                !dtd = iu**(lK - lB) * conjg(tlmplm%tdd(indN, iatom, 1))
                !utd = iu**(lK - lB) * conjg(tlmplm%tdu(indN, iatom, 1))
                !dtu = iu**(lK - lB) * conjg(tlmplm%tud(indN, iatom, 1))
                !utu = conjg(tlmplm%tuu(indN, iatom, 1))
                !utd = conjg(tlmplm%tdu(indN, iatom, 1))
                !dtu = conjg(tlmplm%tud(indN, iatom, 1))
                !dtd = conjg(tlmplm%tdd(indN, iatom, 1))
                utu = conjg(loosetdin(ind, iatom, 1, 1))
                dtu = conjg(loosetdin(ind, iatom, 1, 2))
                utd = conjg(loosetdin(ind, iatom, 1, 3))
                dtd = conjg(loosetdin(ind, iatom, 1, 4))

                !if (.FALSE..and.(lB.eq.lK.and.mB.eq.mK)) then
                !  write(4009,*) 'Huh. We got here.'
                !  utu=utu+0.5*usdus%us(lB,1,1)*usdus%dus(lB,1,1)*atoms%rmt(1)**2
                !  utd=utd+0.5*usdus%us(lB,1,1)*usdus%duds(lB,1,1)*atoms%rmt(1)**2
                !  dtu=dtu+0.5*usdus%uds(lB,1,1)*usdus%dus(lB,1,1)*atoms%rmt(1)**2
                !  dtd=dtd+0.5*usdus%uds(lB,1,1)*usdus%duds(lB,1,1)*atoms%rmt(1)**2
                !end if

                !if (.FALSE.) then
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 1
                !  write(111,*) real(utu), aimag(utu)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 1
                !  write(111,*) real(dtu), aimag(dtu)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 1, 2
                !  write(111,*) real(utd), aimag(utd)
                !  write(111,*) lB*(lB+1)+1+mB, lK*(lK+1)+1+mK, 2, 2
                !  write(111,*) real(dtd), aimag(dtd)
                !end if

                utu=utu*(ImagUnit**(lK - lB))
                dtu=dtu*(ImagUnit**(lK - lB))
                utd=utd*(ImagUnit**(lK - lB))
                dtd=dtd*(ImagUnit**(lK - lB))
              end if ! ind >= or < 0

              do nBand = 1, nobd(ikpt, 1)
                ax(nBand) = ax(nBand) + utu * mCoefKb(nBand, lmloK) + utd * mCoefKb(nBand, lmloK + 1)
                bx(nBand) = bx(nBand) + dtu * mCoefKb(nBand, lmloK) + dtd * mCoefKb(nBand, lmloK + 1)
                !ax4(nBand) = ax4(nBand) + utu2 * mCoefKb(nBand, lmloK) + utd2 * mCoefKb(nBand, lmloK + 1)
                !bx4(nBand) = bx4(nBand) + dtu2 * mCoefKb(nBand, lmloK) + dtd2 * mCoefKb(nBand, lmloK + 1)

                coKsh = 1
                ! LOs deactivated
                if (.false.) then
                  do loKet = 3, nRadFun(lK, itype)
                    coKsh = coKsh + 1
                    ! cannot be considered seperately!!!
                    ! todo LO: utulo bzw. dtulo ausgeben mit entsprechend lmB ud mK indices
                    utulo = ImagUnit**(lK - lB) * conjg(tlmplm%tuulo(lmB, mK, iloTable(loKet, lK, itype) + mlo, 1, 1)) ! TODO: These needed a second spin index.
                    dtulo = ImagUnit**(lK - lB) * conjg(tlmplm%tdulo(lmB, mK, iloTable(loKet, lK, itype) + mlo, 1, 1))
                    ax(nBand) = ax(nBand) + utulo * mCoefKb(nBand, lmloK + coKsh)

                    bx(nBand) = bx(nBand) + dtulo * mCoefKb(nBand, lmloK + coKsh)  ! todo LO: should ax and bx really be an array
                  end do
                end if ! LOs false

                coBsh = 1
                if (.false.) then
                  do loBra = 3, nRadFun(lB, itype)
                    coBsh = coBsh + 1
                    ! todo LO: if the LOs are sorted in a special way here even an exit statement could be possible
                    ! what is with the ls they either have to be equal or not
                    ! we only calculate a triangular matrix aren't we missing entries?
                    ! don't forget the shift mlo and mlolo
                    ! indices have to be vice versa to utulo
                    ! todo LO: ulotu utulo vertauschen einzeln ausgeben für eine lm lmp kombi, cx(nBand?)
                    ulotu = ImagUnit**(lK - lB) * tlmplm%tuulo(lmK, mB, iloTable(loBra, lB, itype) + mlo, 1, 1) ! TODO: These needed a second spin index.
                    ulotd = ImagUnit**(lK - lB) * tlmplm%tdulo(lmK, mB, iloTable(loBra, lB, itype) + mlo, 1, 1)
                    ! give out ulotu, ulotd, are all lmloK acessed here?
                    ! the second index gives the shift of the LOs and determins which ccof has to be used
                    cx(nBand, coBsh) = cx(nBand, coBsh) + ulotu * mCoefKb(nBand, lmloK) + ulotd * mCoefKb(nBand, lmloK + 1)
                    !todo LO: this  has to be a bit more beneath otherwise more than one LO is not working!
                    coKsh = 1 !when debugging LOs check whether coKsh is set correctly
                    do loKet = 3, nRadFun(lK, itype)
                      coKsh = coKsh + 1
                      ! This construct ist due to the fact that for LOs only the triangular matrix values are stored in tuloulo
                      if ( iloTable(loBra, lB, itype) < iloTable(loKet, lK, itype ) ) then
                        loBraKet = ( ( iloTable(loKet, lK, itype) - 1 ) * iloTable(loKet, lK, itype) ) / 2 &
                                                                                                      & + iloTable(loBra, lB, itype)
                        ulotulo = ImagUnit**(lK - lB) * tlmplm%tuloulo(mK, mB, loBraKet + mlolo, 1, 1) ! TODO: These needed a second spin index.
                      else
                        loBraKet = ( ( iloTable(loBra, lB, itype) - 1 ) * iloTable(loBra, lB, itype) ) / 2 &
                                                                                                      & + iloTable(loKet, lK, itype)
                        ulotulo = ImagUnit**(lK - lB) * conjg(tlmplm%tuloulo(mB, mK, loBraKet + mlolo, 1, 1)) ! TODO: These needed a second spin index.
                      end if
                      cx(nBand, coBsh) =  cx(nBand, coBsh) + ulotulo * mCoefKb(nBand, lmloK + coKsh)
                    end do ! loKet
                  end do ! loBra
                end if

              end do ! nBand
              if (.FALSE..and.present(almkg)) then
                do iG=1, nvk
                  ax2(iG) = ax2(iG) + utu*almkg(lmK+1, iG) + utd*blmkg(lmK+1, iG)
                  bx2(iG) = bx2(iG) + dtu*almkg(lmK+1, iG) + dtd*blmkg(lmK+1, iG)
                  ax3(iG) = ax3(iG) + utu*dalmkg(lmK+1, iG) + utd*dblmkg(lmK+1, iG)
                  bx3(iG) = bx3(iG) + dtu*dalmkg(lmK+1, iG) + dtd*dblmkg(lmK+1, iG)
                end do
              end if
            else
              !if (.FALSE.) then
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
          end do ! m
        end do ! l
        do nBand = 1, nobd(ikpt, 1)
          do pBand = 1, ne(ikpq)
          ! Complex conjugation because due to the recycling of tlmplm and tlo we have to calculate here h* instead of h
            h(pBand, nBand) = h(pBand, nBand) +  conjg(mCoefBp(pBand, lmloB)) * ax(nBand) &
                                                                                    & + conjg(mCoefBp(pBand, lmloB + 1)) * bx(nBand)
            !hsurf(pBand, nBand) = hsurf(pBand, nBand) +  conjg(mCoefBp(pBand, lmloB)) * ax4(nBand) &
                                                                                    !& + conjg(mCoefBp(pBand, lmloB + 1)) * bx4(nBand)
            !write(112,*) pband, nband
            !write(112,*) hsurf(pBand, nBand)

            if (.false.) then
              coBsh = 1
              do loBra = 3, nRadFun(lB, itype)
                coBsh = coBsh + 1
                ! coBsh = coefficient Bra Shift
                h(pBand, nBand) = h(pBand, nBand) +  cx(nBand, coBsh) * conjg(mCoefBp(pBand, lmloB + coBsh))
              end do
            end if
          end do !pBand
        end do ! nBand
        if (.FALSE..and.present(almkg)) then
          do iGq=1, nvkq
            do iG=1, nvk
              hSumG(iGq,iG)=hSumG(iGq,iG)+conjg(almkgq(lmB+1, iGq))*ax3(iG)+ &
                                         & conjg(blmkgq(lmB+1, iGq))*bx3(iG)
              hSumGq(iGq,iG)=hSumGq(iGq,iG)+conjg(dalmkgq(lmB+1, iGq))*ax2(iG)+ &
                                         & conjg(dblmkgq(lmB+1, iGq))*bx2(iG)
            end do
          end do
        end if

        lmloB = lmloB + nRadFun(lB, itype)
      end do ! mp
    end do ! lp

    if (.false..and..FALSE..and.present(almkg)) then
        if (ikpt.eq.1.and.idir.eq.1) then
          open(109,file='000_ME_H0_kG_MT',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
          open(110,file='000_ME_H0_kGprq_MT',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
          open(113,file='000_ME_S0_kG_MT',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
          open(112,file='000_ME_S0_kGprq_MT',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
        else
          open(109,file='000_ME_H0_kG_MT',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
          open(110,file='000_ME_H0_kGprq_MT',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
          open(113,file='000_ME_S0_kG_MT',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
          open(112,file='000_ME_S0_kGprq_MT',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
        end if
        do iG=1, nvk
          do iGq=1, nvkq
            write(109,*) iGq, iG, idir
            write(109,*) real(hSumG(iGq,iG)), aimag(hSumG(iGq,iG))
            write(110,*) iGq, iG, idir
            write(110,*) real(hSumGq(iGq,iG)), aimag(hSumGq(iGq,iG))
            write(113,*) iGq, iG, idir
            write(113,*) real(sSumG(iGq,iG)), aimag(sSumG(iGq,iG))
            write(112,*) iGq, iG, idir
            write(112,*) real(sSumGq(iGq,iG)), aimag(sSumGq(iGq,iG))
          end do
        end do
        close(109)
        close(110)
        close(113)
        close(112)

    end if

    if (present(mat_elH)) then
      mat_elH(:,:) = mat_elH(:,:) + hSumG(:,:) + hSumGq(:,:)
      mat_elS(:,:) = mat_elS(:,:) + sSumG(:,:) + sSumGq(:,:)
    end if

    !deallocate(ax,ax2,ax3,bx,bx2,bx3)
    !deallocate(ax4,bx4)

  end subroutine calcHS0MT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the kinetic energy and overlap part of the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !>
  !> @details
  !> See Equations 7.125 and 7.126 and Algorithm 15 (dissertation CRG)
  !>
  !> @param[in]  atoms      : Atoms type, see types.f90
  !> @param[in]  cell       : Unit cell type, see types.f90.
  !> @param[in]  kpts       : K-points type, see types.f90.
  !> @param[in]  qpts       : Q-point set represented by a k-points type, see types.f90
  !> @param[in]  iDatom     : Index of displaced atom
  !> @param[in]  iDtype     : Index of atom type to which displaced atom belongs
  !> @param[in]  nrKetBands : Number of occupied bands at k
  !> @param[in]  nrBraBands : Number of all bands at k + q
  !> @param[in]  ikpt       : Index of k-point in k-point set
  !> @param[in]  ikpq       : Index of k + q (k' backfolded) in k-point set
  !> @param[in]  iqpt       : Index of q in q-point set
  !> @param[in]  iDdir      : Index of displacement direction
  !> @param[in]  nv         : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]  Gbas       : G-basis vectors
  !> @param[in]  ilst       : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !> @param[in]  kpq2kPrVec : Backfolding vector from k + q in 2nd Brillouin zone to k' in 1st Brillouin zone
  !> @param[in]  zBra       : Kohn-Sham eigenvectors at k+q
  !> @param[in]  zKet       : Kohn-Sham eigenvectors at k
  !> @param[out] sIntTeps   : Kinetic energy part of the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !> @param[out] sInt       : Overlap part of the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcSintKinEps( atoms, cell, kpts, qpts, iDtype, iDatom, nrKetBands, nrBraBands, ikpt, ikpq, iqpt, iDdir, nv, &
      Gbas, ilst, zBra, zKet, sIntTeps, sInt, kpq2kPrVec, mat_elT, mat_elS )

    use m_types, only : t_atoms, t_cell, t_kpts
    use m_ylm
    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_cell),              intent(in)  :: cell
    type(t_kpts),              intent(in)  :: kpts
    type(t_kpts),              intent(in)  :: qpts

    ! Scalar parameters
    integer,                   intent(in)  :: iDatom
    integer,                   intent(in)  :: iDtype
    integer,                   intent(in)  :: nrKetBands
    integer,                   intent(in)  :: nrBraBands
    integer,                   intent(in)  :: ikpt
    integer,                   intent(in)  :: ikpq
    integer,                   intent(in)  :: iqpt
    integer,                   intent(in)  :: iDdir

    ! Array parameters
    integer,                   intent(in)  :: nv(:)
    integer,                   intent(in)  :: Gbas(:, :)
    integer,                   intent(in)  :: ilst(:, :, :)
    integer,                   intent(in)  :: kpq2kPrVec(:, :, :)
    complex,                  intent(in)  :: zBra(:, :)
    complex,                  intent(in)  :: zKet(:, :)
    complex,                   intent(out) :: sIntTeps(:, :)
    complex,                   intent(out) :: sInt(:, :)
    complex, optional,         intent(inout) :: mat_elT(:,:)
    complex, optional,         intent(inout) :: mat_elS(:,:)
    ! Local scalar variable
    ! constPreFacs : constant prefactors
    ! sumM         : sum of m-dependent factors
    ! i* : loop indices
    ! GGqNorm : holds norm of variables of arguments of exponential function being Raileigh expanded
    ! lm: index for l and m
    ! tempNoMdep : temporary complex prefactor which is not dependent on m
    complex                                :: constPreFacs
    complex                                :: sumM
    integer                                :: iBband
    integer                                :: iKband
    integer                                :: iKbas
    integer                                :: iBbas
    real                                   :: GGqNorm
    integer                                :: lm, iG, iGq
    complex                                :: tempNoMdep

    ! Local array variable
    ! basisTemp : Temparray containing the basis vectors
    ! kinEffect : Temparray containing the effect of the kinetic operator
    complex,       allocatable             :: basisTemp(:, :), surfIntTG(:, :), surfIntSG(:, :)
    real,          allocatable             :: kinEffect(:), kinEffect2(:,:)
    real                                   :: GGqInt(3)
    real                                   :: GGqCart(3)
    complex                                :: ylm(4)
    real                                   :: sphBesJ(0 : 1)
    real                                   :: GpkKCart(3), gpkqBCart(3)

    ! Initialization
    allocate(basisTemp(nv(ikpq), nv(ikpt)))
    allocate(kinEffect(nv(ikpt)))
    allocate(kinEffect2(nv(ikpq),nv(ikpt)))

    basisTemp(:, :) = cmplx(0., 0.)
    kinEffect(:) = 0.
    sIntTeps(:, :) = cmplx(0., 0.)
    sInt(:, :) = cmplx(0., 0.)

    if (.FALSE.) then
      allocate(surfIntTG(nv(ikpq), nv(ikpt)))
      allocate(surfIntSG(nv(ikpq), nv(ikpt)))
      surfIntTG(:, :) = cmplx(0., 0.)
      surfIntSG(:, :) = cmplx(0., 0.)
    end if

    ! The minus comes from 7.121c
    constPreFacs = -fpi_const * ImagUnit / cell%omtil * atoms%rmt(iDtype)**2

    ! For sake of performance, we first calculate the G-vector dependent quantities.
    do iKbas = 1, nv(ikpt)

      ! Effect of kinetic energy operator
      gpkKCart(1:3) = matmul( cell%bmat(1:3, 1:3), Gbas(1:3, ilst(iKbas, ikpt, 1)) + kpts%bk(1:3, ikpt) )
      kinEffect(iKbas)= (norm2(gpkKCart))**2

      do iBbas = 1, nv(ikpq)

        GGqInt(1:3) = Gbas(1:3, ilst(iKbas, ikpt, 1)) - Gbas(1:3, ilst(iBbas, ikpq, 1)) - qpts%bk(1:3, iqpt) &
                                                                                                     & - kpq2kPrVec(1:3, ikpt, iqpt)
        GGqCart(1:3) = matmul(cell%bmat(1:3, 1:3), GGqInt(1:3))
        gpkqBCart(1:3) = matmul( cell%bmat(1:3, 1:3), Gbas(1:3, ilst(iBbas, ikpq, 1)) + kpts%bk(1:3, ikpq) )
        !CRGfix
        !gpkqBCart(1:3) = matmul( cell%bmat(1:3, 1:3), Gbas(1:3, ilst(iBbas, ikpq, 1)) + kpts%bk(1:3, ikpq) &
        !                       & + qpts%bk(1:3, iqpt) + kpq2kPrVec(1:3, ikpt, iqpt))
        kinEffect2(iBbas,iKbas)=0.5*(gpkqBCart(1)*gpkKCart(1)+gpkqBCart(2)*gpkKCart(2)+gpkqBCart(3)*gpkKCart(3))
        !CRGfix
        !kinEffect2(iBbas,iKbas)=(norm2(gpkqBCart))**2
        GGqNorm = norm2(GGqCart(:))
        ! We only have l = 1 because we have used the orthogonality relation of the spherical
        ! harmonic between the Y_lm from the natural coordinates and the Y_lm from the Rayleigh expansion
        ylm(:) = cmplx(0., 0.)
        call ylm4(1, GGqCart, ylm)
        sphBesJ(:) = 0.
        call sphbes(1, GGqNorm * atoms%rmt(iDtype), sphBesJ)
        tempNoMdep =  sphBesJ(1) * exp(ImagUnit * tpi_const * dot_product(real(GGqInt(1:3)), atoms%taual(1:3, iDatom)))
        sumM = (0., 0.)
        ! For l = 1, we only have -1 < m < 1, so lm = 2, 3, 4
        do lm = 2, 4
          ! We do not have the conjugated natural coordinate expansion coefficients because only the Y_lm were conjugated
          ! by multiplying them with a factor (-1)^m canceling away with the same factor stemming from ylm not being
          ! conjugated, although it is claimed in the Rayleigh expansion
          sumM = sumM + conjg(c_im(iDdir, lm - 1)) * conjg(ylm(lm))
        end do ! t
        basisTemp(iBbas, iKbas) = tempNoMdep * sumM
      end do ! iBbas
    end do ! iKbas

    ! Then we multiply the band dependent quantities
    do iKband = 1, nrKetBands
      do iBband = 1, nrBraBands
        do iKbas = 1, nv(ikpt)
          do iBbas = 1, nv(ikpq)
            if (.FALSE.) then
              sIntTeps(iBband, iKband) = sIntTeps(iBband, iKband) + constPreFacs  &
                & * conjg(zBra(iBbas, iBband)) * (kinEffect2(iBbas,iKbas)) * basisTemp(iBbas, iKbas) * zKet(iKbas, iKband)
              !CRGfix
              !sIntTeps(iBband, iKband) = sIntTeps(iBband, iKband) + constPreFacs  &
              !  & * conjg(zBra(iBbas, iBband)) * 0.25*(kinEffect2(iBbas,iKbas)+kinEffect(iKbas)) &
              !  & * basisTemp(iBbas, iKbas) * zKet(iKbas, iKband)
            else
              sIntTeps(iBband, iKband) = sIntTeps(iBband, iKband) + constPreFacs  &
                & * conjg(zBra(iBbas, iBband)) * (0.5 * kinEffect(iKbas)) * basisTemp(iBbas, iKbas) * zKet(iKbas, iKband)
            end if
            sInt(iBband, iKband) = sInt(iBband, iKband) + constPreFacs  &
              & * conjg(zBra(iBbas, iBband)) * basisTemp(iBbas, iKbas) * zKet(iKbas, iKband)
          end do ! iBbas
        end do ! iKbas
      end do !iKband
    end do ! iBband

    if (.FALSE.) then
      if (ikpt.eq.1.and.iDdir.eq.1) then
        open(109,file='000_ME_eThet1_IR',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
        open(110,file='000_ME_Thet1_IR',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
      else
        open(109,file='000_ME_eThet1_IR',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
        open(110,file='000_ME_Thet1_IR',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
      end if
      do iG=1, nv(ikpt)
        do iGq=1, nv(ikpq)
          surfIntTG(iGq,iG)=constPreFacs*kinEffect2(iGq,iG)*basisTemp(iGq, iG)
          surfIntSG(iGq,iG)=constPreFacs*basisTemp(iGq, iG)
          !write(109,*) ikpt, iGq, iG, iDdir
          !write(109,*) real(surfIntTG(iGq,iG)), aimag(surfIntTG(iGq,iG))
          !write(110,*) iGq, iG, iDdir
          !write(110,*) real(surfIntSG(iGq,iG)), aimag(surfIntSG(iGq,iG))
        end do
      end do
      close(109)
      close(110)

    end if

    if (present(mat_elS)) then
      mat_elS(:,:) = mat_elS(:,:) + surfIntSG(:,:)
      mat_elT(:,:) = mat_elT(:,:) + surfIntTG(:,:)
    end if

  end subroutine calcSintKinEps


  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the lm coefficients for the unification of normal vector and potential for the potential part of the the surface
  !> integral Psi_k+q^IR (H - eps) Psi_k^IR
  !>
  !> @details
  !> See 7.127 (dissertation CRG)
  !>
  !> @param[in]   atoms       : Atoms type, see types.f90
  !> @param[in]   stars       : Stars type, see types.f90.
  !> @param[in]   cell        : Unit cell type, see types.f90.
  !> @param[in]   iDtype      : Index of atom type to which displaced atom belongs
  !> @param[in]   iDatom      : Index of displaced atom
  !> @param[in]   ngdp        : Number of G-vectors for potentials and densities
  !> @param[in]   coScale     : Scaling factor for lmax cutoff of Rayleigh expansion
  !> @param[in]   gdp         : G-vectors of potentials and densities
  !> @param[in]   vEff0IRpwUw : Spherical harmonic coefficients of unperturbed and converged muffin-tin effective potential
  !>                            parsed from Fleur
  !> @param[out]  veffUvIR    : lm coefficients for the unification of normal vector and potential for the potential part of the
  !>                            the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !>
  !> @details
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine IRcoeffVeffUv( atoms, stars, cell, iDtype, iDatom, ngdp, coScale, gdp, veffUvIR, vEff0IRpwUw )

    use m_gaunt, only : gaunt1
    use m_sphbes
    use m_ylm
    use m_types, only : t_atoms, t_stars, t_cell
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameter
    type(t_atoms),             intent(in)  :: atoms
    type(t_stars),             intent(in)  :: stars
    type(t_cell),              intent(in)  :: cell

    ! Scalar parameter
    integer,                   intent(in)  :: iDtype
    integer,                   intent(in)  :: iDatom
    integer,                   intent(in)  :: ngdp
    integer,                   intent(in)  :: coScale

    ! Array parameter
    integer,                   intent(in)  :: gdp(:, :)
    complex,                   intent(in)  :: vEff0IRpwUw(:, :)
    complex,                   intent(out) :: veffUvIR(:, :)

    ! Scalar variables
    integer                                :: iG
    complex                                :: phaseFac
    complex                                :: factG
    integer                                :: oqn_l
    integer                                :: mqn_m
    integer                                :: oqn_l2p
    integer                                :: mqn_m1p
    complex                                :: factL
    integer                                :: lm
    complex                                :: factLM
    integer                                :: mqn_m2p
    integer                                :: lm2p
    real                                   :: gauntCoeff
    integer                                :: idir
    integer                                :: lmaxScale

    ! Array variables
    real,          allocatable             :: sbes(:)
    complex,       allocatable             :: ylm(:)
    real                                   :: gCart(3)

    !todo discuss cutoffs!
    lmaxScale = coScale * atoms%lmax(iDtype)

    if ( coScale /= 1 ) CALL juDFT_error('Overthink norm of ylm',calledby ="IRcoeffVeffUv")

    allocate( sbes(0:lmaxScale) )
    allocate( ylm((lmaxScale + 1)**2) )

    sbes(:) = 0.
    ylm(:) = cmplx(0., 0.)

    veffUvIR(:, :) = cmplx(0., 0.)

    ! Unite the unit vector with the effective potential determing a coefficient dependent on direction and lm
    do iG  = 1, ngdp

      phaseFac = exp( ImagUnit * tpi_const * dot_product(gdp(1:3, iG), atoms%taual(1:3, iDatom)) )
      factG = fpi_const * vEff0IRpwUw(iG, 1) * phaseFac
      !if (norm2(real(gdp(1:3, iG))).lt.10e-6) then
      !!  factG = 1.53*factG
      !  do idir = 1, 3
      !    do lm = 2, 4
      !      veffUvIR(idir, lm) = veffUvIR(idir, lm) + c_im(idir, lm - 1) * vEff0IRpwUw(iG, 1)
      !    end do
      !  end do ! idir
      !  cycle
      !end if
      gCart(1:3) = matmul( cell%bmat(1:3, 1:3), gdp(1:3, iG) )
      sbes(:) = 0.
      call sphbes(lmaxScale, norm2(gCart(1:3)) * atoms%rmt(iDtype), sbes)

      ylm(:) = cmplx(0., 0.)
      call ylm4( lmaxScale, gCart, ylm )

      do oqn_l = 0, lmaxScale
        factL = factG * ImagUnit**oqn_l * sbes(oqn_l)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          factLM = factL * conjg(ylm(lm))
          do mqn_m1p = -1, 1
            do oqn_l2p = abs(oqn_l - 1), min(oqn_l + 1, atoms%lmax(iDtype))
              if (mod(oqn_l+oqn_l2p,2).eq.0) cycle ! quite superfluous
              mqn_m2p = mqn_m + mqn_m1p
              ! Despite Gaunt selection rule |m''|< l''
              if ( abs(mqn_m2p) > oqn_l2p ) cycle
              lm2p = oqn_l2p * (oqn_l2p + 1) + 1 + mqn_m2p
              gauntCoeff = gaunt1( oqn_l2p, 1, oqn_l, mqn_m2p, mqn_m1p, mqn_m, lmaxScale )
              do idir = 1, 3
                veffUvIR(idir, lm2p) = veffUvIR(idir, lm2p) + factLM * c_im(idir, mqn_m1p + 2) * gauntCoeff
              end do ! idir
            end do ! oqn_l2p
          end do ! mqn_m1p
        end do ! mqn_m
      end do ! oqn_l
    end do ! iG

  end subroutine  IRcoeffVeffUv

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the potential part of the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !>
  !> @details
  !> See 7.128 (dissertation CRG)
  !>
  !> @param[in]  atoms   : Atoms type, see types.f90
  !> @param[in]  kpts    : K-points type, see types.f90.
  !> @param[in]  cell    : Unit cell type, see types.f90.
  !> @param[in]  dimens  : Dimension type, see types.f90.
  !> @param[in]  ikpt    : Index of k-point in k-point set
  !> @param[in]  ikpq    : Index of k + q (k' backfolded) in k-point set
  !> @param[in]  coScale : Scaling factor for lmax cutoff of Rayleigh expansion
  !> @param[in]  iDtype  : Index of atom type to which displaced atom belongs
  !> @param[in]  iDatom  : Index of displaced atom
  !> @param[in]  gbas    : G-basis vectors
  !> @param[in]  zeta    : lm coefficients for the unification of normal vector and potential for the potential part of the
  !>                       the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !> @param[in]  nv      : Number of LAPW G-basis vectors for given k-point.
  !> @param[in]  ne      : Number of eigenvalues per k-point.
  !> @param[in]  nobd    : Number of occupied bands per k-point and spin
  !> @param[in]  ilst    : For various k-points G-basis vectors occur more than once, thus they are only stored once in juPhon.
  !> @param[in]  z       : Kohn-Sham eigenvectors
  !> @param[out) surfInt : Calculates the potential part of the surface integral Psi_k+q^IR (H - eps) Psi_k^IR
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcSfVeffFast( atoms, input, kpts, qpts, cell, ikpt, ikpq, iqpt, kpq2kPrVec, coScale, gbas, zeta, nv, ne, nobd, ilst, z, iDtype, iDatom, &
      & surfInt, mat_el )

    use m_types, only : t_atoms, t_input, t_kpts, t_cell
    use m_gaunt, only : gaunt1
    use m_ylm
    use m_sphbes

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in)  :: atoms
    type(t_input),             intent(in)  :: input
    type(t_kpts),              intent(in)  :: kpts
    type(t_kpts),              intent(in)  :: qpts
    type(t_cell),              intent(in)  :: cell

    ! Scalar parameters
    integer,                   intent(in)  :: ikpt
    integer,                   intent(in)  :: ikpq
    integer,                   intent(in)  :: iqpt
    integer,                   intent(in)  :: kpq2kPrVec(:, :, :)
    integer,                   intent(in)  :: coScale
    integer,                   intent(in)  :: iDtype
    integer,                   intent(in)  :: iDatom

    ! Array parameters
    integer,                   intent(in)  :: gbas(:, :)
    complex,                   intent(in)  :: zeta(:, :)
    integer,                   intent(in)  :: nv(:, :)
    integer,                   intent(in)  :: ne(:)
    integer,                   intent(in)  :: nobd(:, :)
    integer,                   intent(in)  :: ilst(:, :, :)
    complex,                  intent(in)  :: z(:, :, :, :)
    complex,                   intent(out) :: surfInt(:, :, :)
    complex, optional,         intent(inout) :: mat_el(:,:,:)

    ! Scalar variables
    real                                   :: pref, pref2
    complex                                :: phaseFac
    integer                                :: ibandK
    integer                                :: ibandB
    integer                                :: oqn_l
    integer                                :: lm_pre
    complex                                :: factL
    integer                                :: mqn_m
    integer                                :: lm
    integer                                :: iG, iGq
    integer                                :: idir
    integer                                :: oqn_l2p
    integer                                :: lm2pPre
    integer                                :: lm2p
    integer                                :: oqn_l1p
    integer                                :: lm1pPre
    integer                                :: mqn_m1p
    integer                                :: lm1p
    complex                                :: gauntCoeff
    integer                                :: mqn_m2p
    integer                                :: lmaxScaled
    real                                   :: GGqNorm
    complex                                :: sumM(3), sumL(3)

    ! Array variables
    complex,       allocatable             :: psiKCoeff(:, :)
    complex,       allocatable             :: psiBCoeff(:, :)
    complex,       allocatable             :: basCoeff(:, :)
    complex,       allocatable             :: coeffkg(:, :), coeffkgq(:, :), surfIntG(:,:,:)
    complex,       allocatable             :: ylm(:)
    real,          allocatable             :: sbes(:)
    real,          allocatable             :: gk(:,:), gkq(:,:)
    real                                   :: gpk(3),gpkq(3),GGqInt(3)
    real                                   :: gpkCart(3),GGqCart(3)


    ! Aaron suggests that for the forces he had expanded the Rayleigh expansion up to 2 lmax
    lmaxScaled = coScale * atoms%lmax(iDtype)

    allocate( psiKCoeff(nobd(ikpt, 1), (lmaxScaled + 1)**2) )
    allocate( psiBCoeff(input%neig, (lmaxScaled + 1)**2))
    allocate( basCoeff(MAXVAL(nv), (lmaxScaled + 1)**2) )
    allocate( ylm((lmaxScaled + 1)**2) )
    allocate( sbes(0: lmaxScaled) )

    psiKCoeff(:, :) = cmplx(0., 0.)
    psiBCoeff(:, :) = cmplx(0., 0.)
    basCoeff(:, :) = cmplx(0., 0.)
    ylm(:) = cmplx(0., 0.)
    sbes(:) = 0.

    surfInt(:, :, :) = cmplx(0., 0.)
    surfIntG(:, :, :) = cmplx(0., 0.)

    if (.FALSE.) then
      allocate( coeffkg(nv(1, ikpt), (lmaxScaled + 1)**2) )
      allocate( coeffkgq(nv(1, ikpq), (lmaxScaled + 1)**2) )
      allocate( gk(3,nv(1, ikpt)))
      allocate( gkq(3,nv(1, ikpq)))
      allocate(surfIntG(nv(1, ikpq),nv(1, ikpt),3))
      coeffkg(:, :) = cmplx(0., 0.)
      coeffkgq(:, :) = cmplx(0., 0.)
    end if

    if ( coScale /= 1 ) CALL juDFT_error('ylmNorm is not correctly set in calcSfVeffFast',calledby ="calcSfVeffFast")

    pref = -fpi_const**2 / cell%omtil * atoms%rmt(iDtype)**2
    ! Calculate wavefunction at k, start with basis function part
    do iG = 1, nv(1, ikpt)
      gpk(1:3) = gbas(1:3, ilst(iG, ikpt, 1)) + kpts%bk(1:3, ikpt)
      phaseFac = exp( ImagUnit * tpi_const * dot_product(gpk(:), atoms%taual(1:3, iDatom)) )
      gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gpk(1:3) )
      if (.FALSE.) then
        gk(1:3,iG) = gpkCart(1:3)
      end if
      ylm(:) = cmplx(0., 0.)
      call ylm4( lmaxScaled, gpkCart, ylm )

      sbes(:) = 0.
      call sphbes( lmaxScaled, norm2(gpkCart(1:3)) * atoms%rmt(iDtype), sbes)
      do oqn_l = 0, lmaxScaled
        lm_pre = oqn_l * (oqn_l + 1) + 1
        factL = ImagUnit**oqn_l * sbes(oqn_l) * phaseFac
        do mqn_m = -oqn_l, oqn_l
          lm = lm_pre + mqn_m
          basCoeff(iG, lm) = basCoeff(iG, lm) + factL * conjg(ylm(lm))
        end do ! mqn_m
      end do ! oqn_l
    end do ! iG

    if (.FALSE.) then
      coeffkg(:,:)=basCoeff(:nv(1, ikpt),:)
    end if

    ! Multiply the wave function expansion coefficients
    do oqn_l = 0, lmaxScaled
      lm_pre = oqn_l * (oqn_l + 1) + 1
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do ibandK = 1, nobd(ikpt, 1)
          psiKCoeff(ibandK, lm) = dot_product(conjg(basCoeff(1:nv(1, ikpt), lm)), z(1:nv(1, ikpt), ibandK, ikpt, 1))
        end do ! ibandK
      end do ! mqn_m
    end do ! oqn_l


    ! Calculate conjugated wavefunction at k', start with basis function part
    basCoeff(:, :) = cmplx(0., 0.)
    do iG = 1, nv(1, ikpq)
      gpk(1:3) = gbas(1:3, ilst(iG, ikpq, 1)) + kpts%bk(1:3, ikpq)
      phaseFac = exp( -ImagUnit * tpi_const * dot_product(gpk(:), atoms%taual(1:3, iDatom)) )

      gpkCart(1:3) = matmul( cell%bmat(1:3, 1:3), gpk(1:3) )

      if (.FALSE.) then
        gkq(1:3,iG) = gpkCart(1:3)
      end if

      ylm(:) = cmplx(0., 0.)
      call ylm4( lmaxScaled, gpkCart, ylm )

      sbes(:) = 0.
      call sphbes(lmaxScaled, norm2(gpkCart(1:3)) * atoms%rmt(iDtype), sbes)
      do oqn_l = 0, lmaxScaled
        lm_pre = oqn_l * (oqn_l + 1) + 1
        factL = conjg(ImagUnit)**oqn_l * sbes(oqn_l) * phaseFac
        do mqn_m = -oqn_l, oqn_l
          lm = lm_pre + mqn_m
          basCoeff(iG, lm) = basCoeff(iG, lm) + factL * ylm(lm)
        end do ! mqn_m
      end do ! oqn_l
    end do ! iG

    if (.FALSE.) then
      coeffkgq(:,:)=basCoeff(:nv(1, ikpq),:)
    end if

    ! Multiply the wave function expansion coefficients
    do oqn_l = 0, lmaxScaled
      lm_pre = oqn_l * (oqn_l + 1) + 1
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do ibandB = 1, ne(ikpq)
          psiBCoeff(ibandB, lm) = dot_product(z(1:nv(1, ikpq), ibandB, ikpq, 1), basCoeff(1:nv(1, ikpq), lm))
        end do ! ibandK
      end do ! mqn_m
    end do ! oqn_l

    ! Unite with precalculated product of potential and unit vector of surface integral. Gaunt selection rules are used.
    do oqn_l2p = 0, lmaxScaled
      lm2pPre = oqn_l2p * (oqn_l2p + 1) + 1
      do mqn_m2p = -oqn_l2p, oqn_l2p
        lm2p = lm2pPre + mqn_m2p
        do oqn_l = 0, lmaxScaled
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do mqn_m = -oqn_l, oqn_l
            lm = lm_pre + mqn_m
            ! Gaunt selection rules, but cutoff at lmax
            do oqn_l1p = abs( oqn_l2p - oqn_l ), min(oqn_l2p + oqn_l, atoms%lmax(iDtype))
              lm1pPre = oqn_l1p * (oqn_l1p + 1) + 1
              ! Gaunt selection rules
              mqn_m1p = mqn_m2p + mqn_m
              if ( abs(mqn_m1p) > oqn_l1p ) cycle
              lm1p = lm1pPre + mqn_m1p
              gauntCoeff = gaunt1(oqn_l1p, oqn_l2p, oqn_l, mqn_m1p, mqn_m2p, mqn_m, lmaxScaled)
              do idir = 1, 3
                do ibandK = 1, nobd(ikpt, 1)
                  do ibandB = 1, ne(ikpq)
                    surfInt(ibandB, ibandK, idir) = surfInt(ibandB, ibandK, idir) + pref * zeta(idir, lm2p) * gauntCoeff &
                                                                                 & * psiBCoeff(ibandB, lm1p) * psiKCoeff(ibandK, lm)
                  end do ! ibandB
                end do ! ibandK

                !if (.FALSE.) then
                !  do iGq=1, nv(1, ikpq)
                !    do iG=1, nv(1, ikpt)
                !      surfIntG(iGq,iG,idir)=surfIntG(iGq,iG,idir)+pref*zeta(idir,lm2p)*gauntCoeff*coeffkgq(iGq,lm1p)*coeffkg(iG,lm)
                !    end do
                !  end do
                !end if

              end do ! idir
            end do ! oqn_l1p
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2p
    end do ! oqn_l2p

    if (.FALSE.) then
      pref2 = -fpi_const / cell%omtil * atoms%rmt(iDtype)**2
      do iG = 1, nv(1, ikpt)
        do iGq = 1, nv(1, ikpq)
          GGqInt(1:3) = Gbas(1:3, ilst(iG, ikpt, 1)) - Gbas(1:3, ilst(iGq, ikpq, 1)) &
                    & - qpts%bk(1:3, iqpt) - kpq2kPrVec(1:3, ikpt, iqpt)
          GGqCart(1:3) = matmul(cell%bmat(1:3, 1:3), GGqInt(1:3))

          GGqNorm = norm2(GGqCart(:))

          ylm(:) = cmplx(0., 0.)
          call ylm4(lmaxScaled, GGqCart, ylm)

          sbes(:) = 0.
          call sphbes(lmaxScaled, GGqNorm * atoms%rmt(iDtype), sbes)

          sumL = (0., 0.)
          do oqn_l=0, lmaxScaled
            lm_pre = oqn_l * (oqn_l + 1) + 1
            factL = ImagUnit**oqn_l * sbes(oqn_l) * exp(ImagUnit * tpi_const * dot_product(real(GGqInt(1:3)), atoms%taual(1:3, iDatom)))
            sumM = (0., 0.)
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              sumM(:) = sumM(:) + conjg(zeta(:,lm)) * conjg(ylm(lm))
            end do
            sumL(:) = sumL(:) + factL*sumM(:)
          end do
          surfIntG(iGq,iG,:) = pref2 * sumL(:)
        end do
      end do
      !if (ikpt.eq.1) then
      !  open(109,file='000_ME_kGprqkG',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
      !  open(110,file='000_ME_V0Thet1_IR',form='FORMATTED',action='WRITE',position='append',status='REPLACE')
      !else
      !  open(109,file='000_ME_kGprqkG',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
      !  open(110,file='000_ME_V0Thet1_IR',form='FORMATTED',action='WRITE',position='append',status='UNKNOWN')
      !end if
      !do iG=1, nv(1, ikpt)
      !  do iGq=1, nv(1, ikpq)
      !    write(109,*) iGq, iG
      !    write(109,*) gkq(1,iGq), gkq(2,iGq), gkq(3,iGq)
      !    write(109,*) gk(1,iG), gk(2,iG), gk(3,iG)
      !    write(109,*) kpts%bk(1, ikpt), kpts%bk(2, ikpt), kpts%bk(3, ikpt)
      !    write(110,*) iGq, iG, 1
      !    write(110,*) real(surfIntG(iGq,iG,1)), aimag(surfIntG(iGq,iG,1))
      !    write(110,*) iGq, iG, 2
      !    write(110,*) real(surfIntG(iGq,iG,2)), aimag(surfIntG(iGq,iG,2))
      !    write(110,*) iGq, iG, 3
      !    write(110,*) real(surfIntG(iGq,iG,3)), aimag(surfIntG(iGq,iG,3))
      !  end do
      !end do
      !close(109)
      !close(110)
    end if
    if (present(mat_el)) then
      mat_el(:,:,:) = mat_el(:,:,:) + surfIntG(:,:,:)
    end if
  end subroutine calcSfVeffFast

end module m_jpSternhPulaySurface
