module m_jp2ndOrdQuant

  use m_types

  implicit none

  contains


  ! should be correct, has been reviewed
  subroutine GenPsDens2ndOrd(atoms, cell, ngpqdp, G0index, gpqdp, qpt, psDens2ndOrd, testMode)

    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell

    ! Scalar parameter
    integer,                        intent(in) :: ngpqdp
    logical,                        intent(in) :: testMode
    integer,                        intent(out):: G0index

    ! Array parameter
    integer,                        intent(in) :: gpqdp(:, :)
    real,                           intent(in) :: qpt(:)
    complex,           allocatable, intent(out):: psDens2ndOrd(:, :, :, :)

    ! Scalar variable
    integer                                    :: itype
    integer                                    :: ii
    integer                                    :: jj
    integer                                    :: iG
    integer                                    :: iatom
    integer                                    :: iatomTemp
    integer                                    :: ieqat

    ! Array variable
    real,              allocatable             :: sbes(:, :, :)
    real,              allocatable             :: psDensMat(:, :, :) ! Matrix part of pseudo density
    real,              allocatable             :: Gpqext(:, :)
    real,              allocatable             :: Gpq(:, :)
    complex,           allocatable             :: phaseFactor(:, :)
    real,              allocatable             :: GpqextRmt(:, :)
    real,              allocatable             :: prefactor(:)

    ! Initializiation of local arrays
    allocate(psDensMat(ngpqdp, 3, 3))
    allocate(psDens2ndOrd(ngpqdp, 3, atoms%nat, 3))
    allocate(Gpqext(3, ngpqdp))
    allocate(Gpq(3, ngpqdp))
    allocate( sbes( 0 : MAXVAL(atoms%ncv) + 1, ngpqdp, atoms%ntype) )
    allocate( phaseFactor(ngpqdp, atoms%nat) )
    allocate( GpqextRmt(ngpqdp, atoms%ntype) )
    allocate( prefactor(atoms%ntype) )

    sbes(:, :, :) = 0.
    GpqextRmt(:, :) = 0.
    psDensMat(:, :, :) = 0.
    psDens2ndOrd(:, :, :, :) = cmplx(0.0, 0.0)
    Gpqext(:, :) = 0.
    Gpq(:, :) = 0.
    G0index = -1
    prefactor(:) = 1.
    phaseFactor(:, :) = cmplx(0., 0.)

    ! Precalculate the matrix-like part of the pseudodensity
    ! todo If we optimize for the dynamical matrix in the end, we need one non-linear run through memory. As we have to evaluate every matrix
    !      element of the dynamical matrix seperately, it makes sense to shift the indices of the 3x3 matrix after iG to avoid redundant
    !      operations later.
    do iG = 1, ngpqdp
      ! If denominator gets zero skip G-vector, that happens only for q = 0, basically.
      if ( norm2( gpqdp(1:3, iG) + qpt(1:3) ) <= 1e-9 )  then
        G0index = iG
        cycle
      end if

      ! (G+q)(G+q)^T and ((G+q)(G+q)^T - (G+q)^2/3)
      Gpq(1:3, iG) = real(gpqdp(1:3, iG) + qpt(1:3))
      Gpqext(1:3, iG) = matmul(Gpq(1:3, iG), cell%bmat(1:3, 1:3))
      if ( testMode ) then
        psDensMat(iG, 1:3, 1:3) = outerProduct(Gpqext(1:3, iG), Gpqext(1:3, iG))
      else
        psDensMat(iG, 1:3, 1:3) = outerProduct(Gpqext(1:3, iG), Gpqext(1:3, iG)) - (norm2(Gpqext(1:3, iG))**2 * id3x3(1:3, 1:3) / 3.)
      end if
    end do ! iG

    iatom = 0
    do itype = 1, atoms%ntype

      ! Calculate double factorial (2 N + 7)!! in 7.78 (PhD thesis Klüppelberg). Note, atoms%ncv(itype) is already the value taken from
      ! Table 1 in the Weinert paper. For the orbital quantum number l = 2, this leads to
      ! (2 * atoms%ncv(itype) + 3)!! = (2 * N + 2 * l + 3)!! = (2 * N + 7)!!.
      ! The right hand side of the recent equation is consistent with 7.78 (PhD thesis Klüppelberg), where N is the Weinert parameter.
      do ii = 1, 2 * atoms%ncv(itype) + 3, 2
        prefactor(itype) = prefactor(itype) * ii
      end do ! ii

      ! Complete prefactor.
      prefactor(itype) = prefactor(itype) * atoms%zatom(itype) / cell%omtil

      iatomTemp = iatom
      do iG = 1, ngpqdp

        ! No pseudo density contribution for G + q = 0
        if (iG == G0index) cycle

        GpqextRmt(iG, itype) = norm2(Gpqext(1:3, iG)) * atoms%rmt(itype)
        ! sbes is initialized within sphbes with first parameter of sphbes
        call sphbes(atoms%ncv(itype) + 1, GpqextRmt(iG, itype), sbes(0:atoms%ncv(itype) + 1, iG, itype))
        iatom = iatomTemp
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          phaseFactor(iG, iatom) = exp(-tpi_const * ImagUnit * dot_product(Gpq(1:3, iG), atoms%taual(1:3, iatom)))
        end do ! ieqat
      end do ! iG
    end do ! itype

    do jj = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do ii = 1, 3
            do iG = 1, ngpqdp
              if (iG == G0index) cycle
              ! Note that atoms%ncv = N + l, so for l = 2, we need only an exponent or the spherical Bessel function, respectively,
              ! of ncv + 1.
              psDens2ndOrd(iG, ii, iatom, jj) = prefactor(itype) / GpqextRmt(iG, itype)**(atoms%ncv(itype) + 1) &
                &                           * sbes(atoms%ncv(itype) + 1, iG, itype) * psDensMat(iG, ii, jj) * phaseFactor(iG, iatom)
            end do ! iG
          end do ! jj
        end do ! ieqat
      end do ! itype
    end do ! jj

  end subroutine GenPsDens2ndOrd

  subroutine CalcIIEnerg2MatElem( atoms, cell, qpt, ngpqdp, gpqdp, E2ndOrdII )

    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell

    ! Scalar parameter
    integer,                        intent(in) :: ngpqdp

    ! Array parameter
    integer,                        intent(in) :: gpqdp(:, :)
    real,                           intent(in) :: qpt(:)
    complex,                        intent(out):: E2ndOrdII(:, :)

    ! Scalar variable
    integer                                    :: G0index
    integer                                    :: iG
    integer                                    :: iAdir
    integer                                    :: iBdir
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: t
    logical                                    :: testMode
    integer                                    :: iAtype
    integer                                    :: iBtype
    integer                                    :: iAatom
    integer                                    :: iBatom
    integer                                    :: iAeqat
    integer                                    :: iBeqat


    ! Array variables
    complex,           allocatable             :: psDens2ndOrd(:, :, :, :)
    complex,           allocatable             :: expAlpha(:)
    real,              allocatable             :: Gpqext(:, :)
    real,              allocatable             :: sbes(:, :)

    allocate( expAlpha(ngpqdp) )
    allocate( Gpqext(3, ngpqdp))
    allocate( sbes(0:atoms%lmaxd, ngpqdp) )

    E2ndOrdII(:, :) = cmplx(0.0, 0.0)
    expAlpha(:) = cmplx(0.0, 0.0)
    Gpqext(:, :) = 0.0
    sbes(:, :) = 0.0

    ! Generates pseudodensity as it is given in 7.95 Aaron Phd thesis, or for q = 0 in 7.78 Aaron Phd thesis
    ! We leave the test mode on, leading to the fact that the trace is not subtracted. Therefore, we have a non-vanishing diagonal.
    testMode = .true.
    call GenPsDens2ndOrd(atoms, cell, ngpqdp, G0index, gpqdp, qpt, psDens2ndOrd, testMode)

    ! Leave it here so it needs not be calculated 3N x 3N times.
    do iG = 1, ngpqdp
      Gpqext(1:3, iG) = matmul(real(gpqdp(1:3, iG) + qpt(1:3)), cell%bmat(1:3, 1:3))
    end do ! iG

    iAatom = 0
    do iAtype = 1, atoms%ntype
      sbes(:, :) = 0.0
      do iG = 1, ngpqdp
        ! If we precalcalculate the scalar factor we have 7 multiplications less
        !todo we can also only calculate it to 0 not to lmaxd for performance reasons
        call Sphbes( atoms%lmaxd, norm2(Gpqext(1:3, iG)) * atoms%rmt(iAtype), sbes(0:atoms%lmaxd, iG) )
      end do
      do iAeqat = 1, atoms%neq(iAtype)
        iAatom = iAatom + 1
        expAlpha(:) = cmplx(0.0, 0.0)
        do iG = 1, ngpqdp
          expAlpha(iG) = exp(tpi_const * ImagUnit * dot_product(gpqdp(1:3, iG) + qpt(1:3), atoms%taual(1:3, iAatom)))
        end do
        do iAdir = 1, 3
          iBatom = 0
          do iBtype = 1, atoms%ntype
            do iBeqat = 1, atoms%neq(iBtype)
              iBatom = iBatom + 1
              if (iBatom /= iAatom) then
                do iBdir = 1, 3
                  ! Add contribution 7.99
                  do iG = 1, ngpqdp
                    if (iG == G0index) cycle
                    E2ndOrdII( iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1) ) = &
                      & E2ndOrdII( iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1) ) &
                      + atoms%zatom(iAtype) * fpi_const / &
                                                & norm2(Gpqext(1:3, iG))**2 * psDens2ndOrd(iG, iBdir, iBatom, iAdir) * expAlpha(iG)
                  end do ! iG
                end do ! iBdir
              else ! iBatom = iAatom
                oqn_l = 0
                mqn_m = 0
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do iBdir = 1, 3
                  do iG = 1, ngpqdp
                    if (iG == G0index) cycle
                    E2ndOrdII( iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1) ) = &
                      E2ndOrdII( iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1) ) + atoms%zatom(iAtype) * fpi_const *               &
                      & psDens2ndOrd(iG, iBdir, iBatom, iAdir) * sbes(0, iG) / norm2(Gpqext(:, iG))**2 * expAlpha(iG)
                  end do
                  ! Constant term is switched off because it is subtracted away anyway in Equation 7.89, as q-independent
                  if (.false.) then
                    do t = -1, 1
                      E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1)) = &
                        & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom -1)) &
                                 - atoms%zatom(iAtype) * atoms%zatom(iBtype) / atoms%rmt(iBtype)**3 * &
                                            & ( 3 / fpi_const * c_im(iBdir, t + 2) * c_im(iAdir, 2 - t) * (-1)**t  )
                    end do ! t
                  end if ! Constant term switched off
                end do ! iBdir
              end if ! iBatom = iAatom ?
            end do ! iBeqat
          end do ! iBtype

          ! Constant term is switched off because it is subtracted away anyway in Equation 7.89 because it is not dependent on q
          if (.false.) then
            ! Factor 3 because it is within the t-sum
            E2ndOrdII(iAdir + 3 * (iAatom - 1), iAdir + 3 * (iAatom -1)) = E2ndOrdII(iAdir, iAdir) + 3 * atoms%zatom(iAtype) * atoms%zatom(iAtype) / atoms%rmt(iAtype)**3
          end if ! Constant term switched off

        end do ! iAdir
      end do ! iAeqat
    end do ! iAtype

  end subroutine CalcIIEnerg2MatElem

  SUBROUTINE getConstTerm(atoms, cell, constTerm)
     TYPE(t_atoms), INTENT(IN) :: atoms
     TYPE(t_cell),  INTENT(IN) :: cell

     REAL, INTENT(OUT) :: constTerm(:, :)

     INTEGER :: iAlpha, iBeta

     REAL :: tauExtAlpha(3), tauExtBeta(3), vecR(3)

     constTerm = 0.0

     DO iBeta = 1, atoms%nat
        tauExtBeta = MATMUL(cell%amat, atoms%taual(:, iBeta))
        DO iAlpha = 1, atoms%nat
           IF (iBeta==iAlpha) CYCLE
           tauExtAlpha = MATMUL(cell%amat, atoms%taual(:, iAlpha))
           vecR = tauExtAlpha - tauExtBeta
           constTerm(3*(iBeta-1) + 1:3*(iBeta-1) + 3, 3*(iAlpha-1) + 1:3*(iAlpha-1) + 3) &
           & = atoms%zatom(iBeta) * atoms%zatom(iBeta) * (3*outerProduct(vecR,vecR)-id3x3*norm2(vecR)**2) / norm2(vecR)**5
       END DO
     END DO
  END SUBROUTINE

  subroutine CalcIIEnerg2(atoms, cell, qpts, stars, input, iqpt, ngdp, gdp, E2ndOrdII)

    use m_sphbes

    implicit none


    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: qpts
    type(t_stars),                  intent(in) :: stars
    type(t_input),                  intent(in) :: input

    ! Scalar parameter
    integer,                        intent(in) :: iqpt
    integer,                        intent(in) :: ngdp

    ! Array parameter
    integer,                        intent(in) :: gdp(:, :)
    complex,           allocatable, intent(out):: E2ndOrdII(:, :)

    ! Scalar variables
    integer                                    :: iAtype
    integer                                    :: iBtype, iCtype
    integer                                    :: iAeqat
    integer                                    :: iBeqat
    integer                                    :: iAatom
    integer                                    :: iBatom
    integer                                    :: iAdir
    integer                                    :: iBdir, iCdir
    integer                                    :: ngpqdp2km
    integer                                    :: ngpqdp

    LOGICAL :: oldStuff

    ! Array variables
    complex,           allocatable             :: E2ndOrdIIatFinQ(:, :)
    complex,           allocatable             :: E2ndOrdIIatQ0(:, :)
    real,              allocatable             :: constTerm(:, :)
    integer,           allocatable             :: gpqdp(:, :)

    oldStuff = .FALSE.

    ! We get the same results for -q and q, probably because Eii2 is a real quantity
    ! todo only generate G-vectors once in the beginning #56, leave the -q version here so that we can test it to be the same
    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpts%bk(1:3, iqpt), gpqdp )

#ifdef DEBUG_MODE
    if (.false.) then
      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpts%bk(1:3, iqpt), gpqdp )
    end if
#endif

    ! Create final Eii2 matrix and temporary array for passing to the subroutine
    allocate( E2ndOrdII(3 * atoms%nat, 3 * atoms%nat) )
    allocate( E2ndOrdIIatFinQ(3 * atoms%nat, 3 * atoms%nat) )
    allocate( E2ndOrdIIatQ0(3 * atoms%nat, 3 * atoms%nat) )
    allocate( constTerm(3 * atoms%nat, 3 * atoms%nat) )

    E2ndOrdII = cmplx(0.0, 0.0)
    E2ndOrdIIatFinQ = cmplx(0.0, 0.0)
    E2ndOrdIIatQ0 = cmplx(0.0, 0.0)

    ! Call the routine for q = 0
    call CalcIIEnerg2MatElem(atoms, cell, [0.0,0.0,0.0], ngdp, gdp, E2ndOrdIIatQ0)

    ! Call the routine for finite q
    call CalcIIEnerg2MatElem(atoms, cell, -qpts%bk(1:3, iqpt), ngpqdp, gpqdp, E2ndOrdIIatFinQ)

    CALL getConstTerm(atoms, cell, constTerm)

    !write(4543,*) E2ndOrdIIatQ0
    !write(4544,*) E2ndOrdIIatFinQ

#ifdef DEBUG_MODE
    if (.false.) then
      ! We get the same results for -q and q, probably because Eii2 is a real quantity
      call CalcIIEnerg2MatElem(atoms, cell, qpts%bk(1:3, iqpt), ngpqdp, gpqdp, E2ndOrdIIatFinQ)
    end if
#endif
    iAatom = 0
    do iAtype = 1, atoms%ntype
      do iAeqat = 1, atoms%neq(iAtype)
        iAatom = iAatom + 1
        do iAdir = 1, 3
          iBatom = 0
          do iBtype = 1, atoms%ntype
            do iBeqat = 1, atoms%neq(iBtype)
              iBatom = iBatom + 1
              do iBdir = 1, 3
                IF (oldStuff) THEN
                E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                   & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                   & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                   & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                ELSE
                   E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) = &
                 & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) + &
                 & E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                   IF (iBatom/=iAatom) THEN
                  !    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) = &
                  !  & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) - &
                  !  & constTerm(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                   ELSE
                      DO iCtype = 1, atoms%ntype
                        !DO iCdir = 1, 3
                           E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) = &
                         & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) - &
                         & E2ndOrdIIatQ0(iBdir + 3 * (iCtype - 1), iAdir + 3 * (iAatom - 1))
                        !END DO
                        IF (iCtype==iBtype) CYCLE
                        DO iCdir = 1, 3
                        !   E2ndOrdII(iCdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) = &
                        ! & E2ndOrdII(iCdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) + &
                        ! & constTerm(iCdir + 3 * (iCtype - 1), iAdir + 3 * (iAatom - 1))
                        END DO
                      END DO
                   END IF
                END IF
              end do ! iBdir
            end do ! iBeqat
          end do ! iBtype
        end do ! iAdir
      end do ! iAeqat
    end do ! iAtype

  end subroutine CalcIIEnerg2

  ! Creates interstitial and muffin-tin coefficients of the second-variation external potential required for the Dynamical Matrix using
  ! the Weinert method.
  ! Main subroutine for creating the second-order V_ext interstitial coefficients and muffin-tin coefficients using the Weinert method.
  ! Here the output is a matrix!!! So the factor of 2 is multiplied in the dynamical matrix not before multipliying the Qs.
  subroutine GenVext2(atoms, cell, ngdp, gdp, vExt2IR, vExt2MT, testMode)

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell

    ! Scalar parameter
    integer,                        intent(in) :: ngdp
    logical,                        intent(in) :: testMode

    ! Array parameter
    integer,                        intent(in) :: gdp(:, :)
    complex,           allocatable, intent(out):: vExt2IR(:, :, :, :)
    complex,           allocatable, intent(out):: vExt2MT(:, :, :, :)

    ! Scalar variables
    integer                                    :: G0index
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: ii
    integer                                    :: jj
    integer                                    :: iG

    ! Array variables
    real,              allocatable             :: Gext(:, :)
    complex,           allocatable             :: psDens2ndOrd(:, :, :, :)

    allocate(vExt2IR(ngdp, 3, 3, atoms%nat))
    allocate(Gext(3, ngdp))
    vExt2IR = cmplx(0.0, 0.0)

    call GenPsDens2ndOrd(atoms, cell, ngdp, G0index, gdp, [0., 0., 0.], psDens2ndOrd, testMode)

    do iG = 1, ngdp
      Gext(:, iG) = matmul(gdp(:, iG), cell%bmat)
    end do

    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
          do jj = 1, 3
            do ii = 1, 3
              do iG = 1, ngdp
                if ( iG /= G0index ) then
                  vExt2IR(iG, ii, jj, iatom) = fpi_const * psDens2ndOrd(iG, ii, iatom, jj ) / norm2(Gext(:, iG))**2
                end if
              end do
            end do
          end do
        end do
    end do ! iG

    call GenVext2MT(atoms, cell, ngdp, gdp, testMode, vExt2IR, vExt2MT)

  end subroutine GenVext2

  ! Generates second-order V_ext coefficients for the muffin-tin region.
  subroutine GenVext2MT(atoms, cell, ngdp, gdp, testMode, vExt2IR, vExt2MT)

    use m_gaunt, only : gaunt1
    use m_sphbes
    use m_jpGrVeff0, only : Phasy1nSym

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell

    ! Scalar parameters
    integer,                    intent(in)  :: ngdp

    ! Array parameters
    integer,                    intent(in)  :: gdp(:, :)
    logical,                    intent(in)  :: testMode
    complex,                    intent(in)  :: vExt2IR(:, :, :, :)
    complex,       allocatable, intent(out) :: vExt2MT(:, :, :, :)
    ! Scalar variables
    integer                                 :: itype
    integer                                 :: ii
    integer                                 :: jj
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: t
    integer                                 :: tPrime
    integer                                 :: lm
    integer                                 :: imesh
    integer                                 :: iatom
    integer                                 :: ieqat
    integer                                 :: iG
    integer                                 :: iDtype
    integer                                 :: iDeqat
    integer                                 :: iDatom

    !complex(kind=16) :: kindTest ?? This is not portable 

    ! Array variables
    ! lmax is 2, so the upper limit of lm is (2 + 1)**3 = 9
    real                                    :: Gext(3)
    real                                    :: sbes(0:atoms%lmaxd)
    real,          allocatable              :: prfMesh(:, :)
    complex,       allocatable              :: pylm(:, :)
    complex,       allocatable              :: surfIntMat(:, :, :, :, :)
    complex,       allocatable              :: scalFac(:)
    real,          allocatable              :: recMesh(:, :)
    complex                                 :: volIntMat(9, 3, 3)
    complex                                 :: volIntMatTest(9, 3, 3)

    complex :: testMat(3, 3)

    allocate( recMesh(atoms%jmtd, atoms%ntype) )
    recMesh(:, :) = 0.

    ! Calculate the lm dependent matrix part of the 2nd order external potential
    volIntMat = cmplx(0.0, 0.0)
    volIntMatTest = cmplx(0.0, 0.0)

    ! l = 0, m = 0
    !if (testMode) then
    if (.false.) then
      volIntMat(1, 1, 1) = cmplx(sqrt(fpi_const), 0.)
      volIntMat(1, 2, 2) = cmplx(sqrt(fpi_const), 0.)
      volIntMat(1, 3, 3) = cmplx(sqrt(fpi_const), 0.)
    end if

    ! l = 1 has vanishing Gaunt coefficients.

    ! l = 2,  m = -2
    volIntMat(5, 1, 1) = cmplx(sqrt(3 * tpi_const / 5), 0.)
    volIntMat(5, 2, 2) = cmplx(-sqrt(3 * tpi_const / 5), 0.)
    volIntMat(5, 1, 2) = cmplx(0., sqrt(3 * tpi_const / 5))
    volIntMat(5, 2, 1) = cmplx(0., sqrt(3 * tpi_const / 5))

    ! l = 2, m = -1
    volIntMat(6, 1, 3) = cmplx(sqrt(3 * tpi_const / 5), 0.)
    volIntMat(6, 3, 1) = cmplx(sqrt(3 * tpi_const / 5), 0.)
    volIntMat(6, 2, 3) = cmplx(0., sqrt(3 * tpi_const / 5))
    volIntMat(6, 3, 2) = cmplx(0., sqrt(3 * tpi_const / 5))

    ! l = 2, m = 0
    volIntMat(7, 1, 1) = cmplx(-sqrt(fpi_const / 5), 0.)
    volIntMat(7, 2, 2) = cmplx(-sqrt(fpi_const / 5), 0.)
    volIntMat(7, 3, 3) = cmplx(4 * sqrt(pi_const / 5), 0.)

    ! l = 2,  m = 1
    volIntMat(8, 1, 3) = cmplx(-sqrt(3 * tpi_const / 5), 0.)
    volIntMat(8, 3, 1) = cmplx(-sqrt(3 * tpi_const / 5), 0.)
    volIntMat(8, 3, 2) = cmplx(0., sqrt(3 * tpi_const / 5))
    volIntMat(8, 2, 3) = cmplx(0., sqrt(3 * tpi_const / 5))

    ! l = 2, m = 2
    volIntMat(9, 1, 1) = cmplx(sqrt(3 * tpi_const / 5), 0.)
    volIntMat(9, 2, 2) = cmplx(-sqrt(3 * tpi_const / 5), 0.)
    volIntMat(9, 1, 2) = cmplx(0., -sqrt(3 * tpi_const / 5))
    volIntMat(9, 2, 1) = cmplx(0., -sqrt(3 * tpi_const / 5))

    if (.false.) then
    do jj = 1, 3
      do ii = 1, 3
        do oqn_l = 0, 2, 2
          do t = -1, 1
            do tPrime = -1, 1
              ! For l = 0, we only have m = 0
              if (abs(t + tPrime) > oqn_l) cycle
              lm = oqn_l * (oqn_l + 1) + 1 + t + tPrime
              ! + 2 because c_im array starts at 1 and not at - 1
              volIntMatTest(lm, ii, jj) = volIntMatTest(lm, ii, jj) + 3 * c_im(ii, t + 2) * c_im(jj, tPrime + 2)                   &
              &  * gaunt1(oqn_l, 1, 1, t + tPrime, t, tPrime, 2)
            end do
          end do
        end do ! oqn_l
      end do ! ii
    end do ! jj
   ! if (.not.testMode) then
   !   volIntMat(1, :, :) = cmplx(0.0, 0.0)
   ! end if

    do lm = 1, 9
      do jj = 1, 3
        do ii = 1, 3
          write(1104, '(3(i8),2(f15.8))') jj, ii, lm, volIntMatTest(lm, ii, jj)
          write(1105, '(3(i8),2(es15.8))') jj, ii, lm, volIntMatTest(lm, ii, jj) - volIntMat(lm, ii, jj)
        end do ! lm
      end do ! ii
      write(1104, *)
      write(1105, *)
    end do ! jj

    do jj = 1, 3
      write(*, '(i8,1x,2(es15.8))') jj, volIntMat(1, jj, jj)
    end do ! lm
    write(*, *) 'shack'

  end if


!    testMat = 0
!    do lm = 1, 10
!      testMat(:, :) = testMat(:, :) + volIntMat(lm, :, :)
!    end do
!    do ii = 1, 3
!      do jj = 1, 3
!        write(*, *) testMat(jj, ii)
!      end do
!    end do
!
!    NOstopNO
!    do jj = 1, 9
!      write(*, *) jj
!      do ii = 1, 3
!        write(*, '(3(2(f15.8),2x))') volIntMat(jj, ii, 1), volIntMat(jj, ii, 2), volIntMat(jj, ii, 3)
!      end do
!      write(*, *)
!    end do
!   ! do ii = 1, 3
!   !   do jj = 1, 3
!   !     write(*, *) jj, ii, c_im(jj, ii)
!   !   end do
!   ! end do
!   ! NOstopNO
!!    write(*, *) 'gaunt'
!    write(*, *) gaunt1(0, 1, 1, 0, 0, 0, 2)
!    write(*, *) gaunt1(0, 1, 1, 0, -1, 1, 2)
!    write(*, *) gaunt1(0, 1, 1, 0, 1, -1, 2)
!    NOstopNO
   ! NOstopNO
!    write(*, *) gaunt1(0, 1, 1, 0, -1, 1, 2)
!    NOstopNO
!    write(*, *) 1. / 4. / 3.14 !    NOstopNO
!    do t = 1, 3
!      do tPrime = 1, 3
!        write(*, *) tPrime, t, c_im(tPrime, t)
!      end do
!    end do
!    NOstopNO


    ! Calculate the mesh dependent prefactor
    allocate( prfMesh(atoms%jmtd, atoms%ntype) )
    prfMesh = cmplx(0.0, 0.0)
    do itype = 1, atoms%ntype
      ! Numerical optimization. At the MT boundary, we set the factor zero because it is analytically zero to prevent error propagation.
      do imesh = 1, atoms%jri(itype) - 1
        prfMesh(imesh, itype) = atoms%zatom(itype) * (1 - (atoms%rmsh(imesh, itype) / atoms%rmt(itype))**5)
        recMesh(imesh, itype) = atoms%rmsh(imesh, itype)**(-3)
      end do
    end do

    ! Calculate the surface integral part of the Dirichelet boundary problem
    allocate( surfIntMat((atoms%lmaxd + 1)**2, 3, 3, atoms%nat, atoms%nat) )
    allocate( pylm(( atoms%lmaxd + 1 )**2, atoms%nat) )
    allocate( scalFac(( atoms%lmaxd + 1 )**2) )
    surfIntMat(:, :, :, :, :) = cmplx(0.0, 0.0)
    pylm = cmplx(0.0, 0.0)
    scalFac = cmplx(0.0, 0.0)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do iG = 1, ngdp
          if ( all( gdp(:, iG)  == 0 ) ) then
            !todo we should determine the G=0 vector once or put it to the beginning at all
            cycle
          end if
          Gext(1:3) = matmul( gdp(1:3, iG), cell%bmat(1:3, 1:3) )
          pylm(:, :) = cmplx(0.0, 0.0)
          scalFac(:) = cmplx(0.0, 0.0)
          sbes(:) = cmplx(0., 0.)
          call phasy1nSym( atoms, cell, gdp(:, iG), [0., 0., 0.], pylm )
          call Sphbes( atoms%lmax(itype), norm2( Gext ) * atoms%rmt(itype), sbes )
          ! If we precalcalculate the scalar factor we have 7 multiplications less
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + mqn_m + 1
              scalFac(lm) = sbes(oqn_l) * pylm(lm, iatom)
            end do
          end do
          iDatom = 0
          do iDtype = 1, atoms%ntype
            do iDeqat = 1, atoms%neq(iDtype)
              iDatom = iDatom + 1
              do jj = 1, 3
                do ii = 1, 3
                  do oqn_l = 0, atoms%lmax(itype)
                    do mqn_m = -oqn_l, oqn_l
                      lm = oqn_l * (oqn_l + 1) + mqn_m + 1
                      surfIntMat(lm, ii, jj, iDatom, iatom) = surfIntMat(lm, ii, jj, iDatom, iatom) + vExt2IR(iG, ii, jj, iDatom)  &
                                                                                                                     & * scalFac(lm)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(pylm, scalFac)

    ! Add the volume and surface integral contribution of the Dirichelet boundary problem for determing the MT coefficients of the 2nd
    ! order external potential
    ! We have chosen this loop structure sothat the if-statement is as outermost as possible
    allocate(vExt2MT(atoms%jmtd, (atoms%lmaxd + 1)**2, 3 * atoms%nat, 3 * atoms%nat ))
    vExt2MT = cmplx(0.0, 0.0)
    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        do jj = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              if ( iatom == iDatom ) then
                do ii = 1, 3
                ! if trace is subtracted there is no l = 0 component
                  !do oqn_l = 0, 2, 2
                  do oqn_l = 2, 2
                  !do oqn_l = 0, 2
                    do mqn_m = -oqn_l, oqn_l
                      lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                      do imesh = 1, atoms%jri(itype)
                        vExt2MT(imesh, lm, ii + (iatom - 1) * 3, jj + (iDatom - 1) * 3) =                                          &
                                                                & vExt2MT(imesh, lm, ii + (iatom - 1) * 3,  jj + (iDatom - 1) * 3) &
                                                              ! We calculate prfMesh * volIntMat first due to numerical stability
                                                         & - (prfMesh(imesh, itype) * volIntMat(lm, ii, jj)) * recMesh(imesh, itype)
                      end do ! imesh
                    end do ! mqn_m
                  end do ! oqn_l
                end do ! ii
              end if ! within displaced muffin-tin

              ! Add surface integral contribution.
              do ii = 1, 3
                do oqn_l = 0, atoms%lmax(itype)
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itype)
                      vExt2MT(imesh, lm, ii + (iatom - 1) * 3, jj + (iDatom - 1) * 3) =                                            &
                                                                 & vExt2MT(imesh, lm, ii + (iatom - 1) * 3, jj + (iDatom - 1) * 3) &
                                    & + (atoms%rmsh(imesh, itype) / atoms%rmt(itype))**oqn_l * surfIntMat(lm, ii, jj, iDatom, iatom)
                    end do ! imesh
                  end do ! mqn_m
                end do ! oqn_l
              end do ! ii

            end do ! ieqat
          end do ! itype
        end do ! ii
      end do ! iDeqat
    end do ! iDtype

  end subroutine GenVext2MT

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Auxiliary method which is similiar to FLEUR's phasy routine but does not implement symmetry.
  !>
  !>
  !> @details
  !> This method calculates for a given G-vector and q-vector: 4 π i^l exp((Gvec + qpoint) * taual) Y*_lm(Gvec + qpoint)
  !> Opposite to the original phasy1 routine, here, stars and any type of symmetry is not considered here.
  !>
  !> @param[in]  atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in]  cellT      : Contains unit-cell related quantities; definition of its members in types.F90 file.
  !> @param[in]  Gvec       : Current G-vector
  !> @param[in]  qpoint     : Current q-point
  !> @param[out] pylm       : Result of the routine, is atom resolved due to lack of symmetry in first variation of the potentials.
  !--------------------------------------------------------------------------------------------------------------------------------------
  ! Deprecated
  subroutine phasy1lp2nSym(atomsT, cellT, Gvec, qptn, pylm)

    use m_ylm
    use m_types

    implicit none

    ! Scalar Type Arguments
    type(t_atoms),  intent(in)  ::  atomsT
    type(t_cell),   intent(in)  ::  cellT

    ! Array Arguments
    integer,        intent(in)  ::  Gvec(:)
    real,           intent(in)  ::  qptn(:)
    complex,        intent(out) ::  pylm((atomsT%lmaxd + 3)**2, atomsT%nat)

    !---------------------------------------------------------------------------------------------------------------------------------
    ! Local Scalar Variables
    ! iatom : runs over all atoms
    ! itype : runs over all types
    ! ineq  : runs over all equivalent atoms of one atom type
    ! lm    : encodes oqn_l and mqn_m
    ! sf    : stores exponential function
    ! csf   : stores exponential function times 4 π i^l
    ! x     : stores argument of exponential function
    ! mqn_m : magnetic quantum number m
    ! oqn_l : orbital quantum number l
    ! ll1   : auxiliary variable to calculate lm
    !---------------------------------------------------------------------------------------------------------------------------------
    integer                     ::  iatom
    integer                     ::  itype
    integer                     ::  ineq
    integer                     ::  lm
    complex                     ::  sf
    complex                     ::  csf
    real                        ::  x
    integer                     ::  mqn_m
    integer                     ::  oqn_l
    integer                     ::  ll1

    !---------------------------------------------------------------------------------------------------------------------------------
    ! Local Array Variables
    ! fpiul: stores 4 π i^l
    ! Gqext: stores G + q in external coordinates
    ! ylm  : stores Y_lm
    !---------------------------------------------------------------------------------------------------------------------------------
    complex                     ::  fpiul(0:atomsT%lmaxd + 2)
    real                        ::  Gqext(3)
    complex                     ::  ylm((atomsT%lmaxd + 3)**2)

    ! calculates 4 π i^l resolved for every l, not divided by nop because no loop over symmetry operations
    do oqn_l = 0, atomsT%lmaxd + 2
       fpiul(oqn_l) = fpi_const * ImagUnit**oqn_l
    enddo


    ! calculates Y*_lm(\vec{G} + \vec{q}) for every l and m. The argument Gqext must be in external coordinates.
    Gqext = matmul(Gvec + qptn, cellT%bmat)
    call Ylmnorm_init( atomsT%lmaxd + 2 )
    call ylm4(atomsT%lmaxd + 2, Gqext, ylm)
    call Ylmnorm_init( atomsT%lmaxd)
    ylm = conjg(ylm)


    ! calculates first exp(i (G + q) tau)  and multiplies recent factors before storing the final result to pylm
    iatom = 1
    pylm = cmplx(0.,0.)
    do itype = 1, atomsT%ntype
       do ineq = 1, atomsT%neq(itype)
          x = tpi_const * dot_product(Gvec + qptn, atomsT%taual(:, iatom))
          sf = exp(ImagUnit *  x)
          do oqn_l = 0, atomsT%lmax(itype) + 2
             ll1 = oqn_l * (oqn_l + 1) + 1
             csf = fpiul(oqn_l) * sf
             do mqn_m = -oqn_l, oqn_l
                lm = ll1 + mqn_m
                pylm(lm, iatom) = csf * ylm(lm)
             enddo ! mqn_m
          enddo ! oqn_l
          iatom = iatom + 1
       enddo ! ineq
    enddo ! itype

  end subroutine phasy1lp2nSym

  function outerProduct(a, b)

    implicit none

    real, intent(in) :: a(:)
    real, intent(in) :: b(:)

    real             :: outerProduct(size(a), size(b))

    outerProduct(:, :) = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))

  end function outerProduct

  function outerProductME(a, b, i, j)

    use m_juDFT_stop, only : juDFT_error

    implicit none

    complex,    intent(in)  :: a(:)
    complex,    intent(in)  :: b(:)
    integer,    intent(in)  :: i
    integer,    intent(in)  :: j

    complex                 :: outerProductME

    if (i > size(a) .or. j > size(b)) then
      call juDFT_error( 'Wished matrix element is out of scope of outer product matrix.', calledby='outerProductME', &
        & hint='Choose smaller indices.')
    end if

    outerProductME = a(i) * b(j)

  end function outerProductME

  subroutine genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp )

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
    gpqdptemp2kmax(:, :) = 0
    gpqdptemprest(:, :) = 0
    ! From all possible G-vectors in a box, only these are accepted which are element of a sphere with radius gmax which is shifted.
    ! We need a little bit more than k*d because they are thought for a Gmax ball that is not shifted by a q, i.e. |G+q|<Gmax
    do iGx = -(stars%mx1 + 3), (stars%mx1 + 3)
      do iGy = -(stars%mx2 + 3), (stars%mx2 + 3)
        do iGz = -(stars%mx3 + 3), (stars%mx3 + 3)
          Gint = [iGx, iGy, iGz]
          Gpqext =  matmul(real(Gint(1:3) + qpoint(1:3)),cell%bmat) !transform from internal to external coordinates
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

  end subroutine genPertPotDensGvecs

end module m_jp2ndOrdQuant
