module m_jpSetupDynMatHelper

  implicit none

  contains
  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Read-in of z1 from Sternheimer cycle for a certain k- and q-point.
  !>
  !> @details
  !> The z1 from the Sternheimer loop are written out as they depend on the k-points. They are read in here again for a certain k- and q-
  !> point.
  !>
  !> @note
  !> The filling for the z1 is not optimal in memory access, but minimizes array sizes, and ensures a linear run through the array later.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine readInz1( atoms, ikpt, iqpt, ikpq, nobd, nv, z1nG )

    use m_types, only : t_atoms
    use mod_juPhonUtils, only : fopen, fclose
    use m_juDFT_NOstopNO, only : juDFT_warn

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms

    ! Scalar parameter
    integer,                       intent(in)  :: ikpt
    integer,                       intent(in)  :: iqpt
    integer,                       intent(in)  :: ikpq

    ! Array parameter
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: nv(:, :)
    complex,                       intent(out) :: z1nG(:, :, :, :)

    ! Scalar variables
    integer                                    :: iDtype
    integer                                    :: iDeqat
    integer                                    :: iDatom
    integer                                    :: iband
    integer                                    :: idir
    integer                                    :: iBas
    logical                                    :: st

    ! Array variables
    character(len=:), allocatable              :: filename
    character(len=15)                          :: filenameTemp


    iDatom = 0
    z1nG(:, :, :, :) = cmplx(0.0, 0.0)
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        if (ikpt < 10 .and. iqpt >= 10 ) then
          write(filenameTemp, '(a,i1,a,i0,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else if( ikpt >= 10 .and. iqpt < 10) then
          write(filenameTemp, '(a,i1,a,i2.2,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else if( ikpt < 10 .and. iqpt < 10) then
          write(filenameTemp, '(a,i1,a,i2.2,a,i2.2)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        else
          write(filenameTemp, '(a,i1,a,i0,a,i0)') 'z1a', iDatom, 'q', iqpt, 'k', ikpt
        end if
        filename = trim(filenameTemp)
        inquire (file=filename, exist=st)
        if ( st ) then
          call fopen( 1000, name=filename, status='old', action='read', form='unformatted')
          rewind(1000)
          ! ATTENTION: Do not change order of loops, otherwise wrong input!
          do iband = 1, nobd(ikpt, 1)
            do idir = 1, 3
              do iBas = 1, nv(1, ikpq) + atoms%nlotot
                read(1000) z1nG(iBas, idir, iDatom, iband)
              end do ! iBas
            end do ! idir
          end do ! iband
          call fclose(1000)
        else
          call juDFT_warn( filename//' not found! Its z1 is set to zero.', calledby='readInz1', hint='Perform calculation to gather&
            & required z1.' )
        end if
      end do ! iDeqat
    end do ! iDtype

  end subroutine readInz1

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the action of a given effective potential expansion coefficients onto basis function
  !>
  !> @details
  !> 
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcFnsphVarphi(atoms, itype, oqn_lVmin, nRadFun, varphi1, varphi2, fSh, fShNsphVarphi)

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    type(t_atoms), intent(in)  :: atoms

    integer,       intent(in)  :: itype
    integer,       intent(in)  :: oqn_lVmin

    integer,       intent(in)  :: nRadFun(0:, :)
    real,          intent(in)  :: varphi1(:, :, 0:)
    real,          intent(in)  :: varphi2(:, :, 0:)
    complex,       intent(in)  :: fSh(:, :)
    complex,       intent(out) :: fShNsphVarphi(:, :, 0:, :)

    integer                    :: lmp
    integer                    :: oqn_l
    integer                    :: mqn_m
    integer                    :: iradf
    integer                    :: oqn_lV
    integer                    :: lmV_pre
    integer                    :: mqn_mV
    integer                    :: lmV
    integer                    :: mqn_m1Pr
    integer                    :: oqn_l1Pr
    integer                    :: lm1Pr
    integer                    :: imesh
    real                       :: gauntFactor
    complex                    :: vEff0Gaunt

    ! It might be slightly faster to interchange the oqn_lV and oqn_l loops, but then then bugs are easier to produce.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        ! We have due to the p loop a factor 2 without LOs and with LOs 2 + number of LOs but otherwise we would not run linearly through
        ! the arrays.
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do oqn_lV = oqn_lVmin, atoms%lmax(itype)
            lmV_pre = oqn_lV * (oqn_lV + 1)
            do mqn_mV = -oqn_lV, oqn_lV
              lmV = lmV_pre + mqn_mV + 1
              ! The ket is only given until lmax, therefore we need no extension of rbas1/rbas2.
              ! Gaunt selection rule for m
              mqn_m1Pr = mqn_m + mqn_mV
              ! Gaunt selection rule for the l of the bra and
              ! in the non-spherical part of this routine, we have to consider that the bra (due to twofold gradient) has contributions up
              ! to lmax + 2 and the non-spherical potential which has a cut-off of lmax and is multiplied with the ket having a cutoff of
              ! lmax. The product is a quantity expanded until 2 * lmax (+ 1). Since we multiply it with the bra the product only has to
              ! calculated until lmax + 2. Combined with the Gaunt selection rules, we calculate either until a l' from where the 
              ! Gaunt coefficients are zero or if this is larger than lmax + 2 we do not need them so we NOstopNOhere. On the other
              ! hand the minimal l' should be given by the lower border of the Gaunt coefficients or by m' which is also determined
              ! according to a Gaunt selection rule. l' should always be larger or equals m'. Otherwise we get lm indices which are
              ! negative.
       !       do oqn_l1Pr = max(abs(oqn_lV - oqn_l), abs(mqn_m1Pr)), min(oqn_lV + oqn_l, atoms%lmax(itype) + 2)
              do oqn_l1Pr = max(abs(oqn_lV - oqn_l), abs(mqn_m1Pr)), min(oqn_lV + oqn_l, atoms%lmax(itype))
                lm1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + mqn_m1Pr
       !         gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l, mqn_m1Pr, mqn_mV, mqn_m, atoms%lmaxd + 2)
                gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l, mqn_m1Pr, mqn_mV, mqn_m, atoms%lmaxd)
                do imesh = 1, atoms%jri(itype)
                  vEff0Gaunt = fSh(imesh, lmV) * gauntFactor
                  fShNsphVarphi(1, imesh, lm1Pr, lmp) = fShNsphVarphi(1, imesh, lm1Pr, lmp) + vEff0Gaunt * varphi1(imesh, iradf, oqn_l)
                  fShNsphVarphi(2, imesh, lm1Pr, lmp) = fShNsphVarphi(2, imesh, lm1Pr, lmp) + vEff0Gaunt * varphi2(imesh, iradf, oqn_l)
                end do ! imesh
              end do ! oqn_l
            end do ! mqn_mV
          end do ! oqn_lV
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

  end subroutine calcFnsphVarphi

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the action of a given effective potential expansion coefficients onto basis function
  !>
  !> @details
  !> 
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcFnsphGrVarphi(atoms, itype, oqn_lVmin, mqn_m2Pr, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, fSh,&
                                                                                                                    & fShNsphVarphi)

    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    type(t_atoms), intent(in)  :: atoms

    integer,       intent(in)  :: itype
    integer,       intent(in)  :: oqn_lVmin
    integer,       intent(in)  :: mqn_m2Pr

    integer,       intent(in)  :: nRadFun(0:, :)
    real,          intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,          intent(in)  :: grVarphiCh2(:, :, :, -1:)
    integer,       intent(in)  :: grVarphiChLout(:, 0:)
    integer,       intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    complex,       intent(in)  :: fSh(:, :)
    complex,       intent(out) :: fShNsphVarphi(:, :, 0:, :)

    integer                    :: lmp
    integer                    :: oqn_l
    integer                    :: mqn_m
    integer                    :: iradf
    integer                    :: oqn_lV
    integer                    :: lmV_pre
    integer                    :: mqn_mV
    integer                    :: lmV
    integer                    :: mqn_m1Pr
    integer                    :: oqn_l1Pr
    integer                    :: lm1Pr
    integer                    :: imesh
    real                       :: gauntFactor
    complex                    :: vEff0Gaunt
    integer                    :: mqn_m3Pr
    integer                    :: ichan
    integer                    :: oqn_l3Pr

    fShNsphVarphi = cmplx(0., 0.)
    ! It might be slightly faster to interchange the oqn_lV and oqn_l loops, but then then bugs are easier to produce.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        mqn_m3Pr = grVarphiChMout(mqn_m, mqn_m2Pr)
        ! We have due to the p loop a factor 2 without LOs and with LOs 2 + number of LOs but otherwise we would not run linearly through
        ! the arrays.
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do ichan = 1, 2
            oqn_l3Pr = grVarphiChLout(ichan, oqn_l)
            if (oqn_l3Pr < 0 .or. abs(mqn_m3Pr) > oqn_l3Pr .or. oqn_l3Pr > atoms%lmax(itype)) cycle
            do oqn_lV = oqn_lVmin, atoms%lmax(itype)
              lmV_pre = oqn_lV * (oqn_lV + 1)
              do mqn_mV = -oqn_lV, oqn_lV
                lmV = lmV_pre + mqn_mV + 1
                ! The ket is only given until lmax, therefore we need no extension of rbas1/rbas2.
                ! Gaunt selection rule for m
                mqn_m1Pr = mqn_m3Pr + mqn_mV
                ! Gaunt selection rule for the l of the bra and
                ! in the non-spherical part of this routine, we have to consider that the bra (due to twofold gradient) has contributions up
                ! to lmax + 2 and the non-spherical potential which has a cut-off of lmax and is multiplied with the ket having a cutoff of
                ! lmax. The product is a quantity expanded until 2 * lmax (+ 1). Since we multiply it with the bra the product only has to
                ! calculated until lmax + 2. Combined with the Gaunt selection rules, we calculate either until a l' from where the 
                ! Gaunt coefficients are zero or if this is larger than lmax + 2 we do not need them so we NOstopNOhere. On the other
                ! hand the minimal l' should be given by the lower border of the Gaunt coefficients or by m' which is also determined
                ! according to a Gaunt selection rule. l' should always be larger or equals m'. Otherwise we get lm indices which are
                ! negative.
!                do oqn_l1Pr = max(abs(oqn_lV - oqn_l3Pr), abs(mqn_m1Pr)), min(oqn_lV + oqn_l3Pr, atoms%lmax(itype) + 2)
                do oqn_l1Pr = max(abs(oqn_lV - oqn_l3Pr), abs(mqn_m1Pr)), min(oqn_lV + oqn_l3Pr, atoms%lmax(itype))
                  lm1Pr = oqn_l1Pr * (oqn_l1Pr + 1) + mqn_m1Pr
!                  gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l3Pr, mqn_m1Pr, mqn_mV, mqn_m3Pr, atoms%lmaxd + 2)
                  gauntFactor = gaunt1( oqn_l1Pr, oqn_lV, oqn_l3Pr, mqn_m1Pr, mqn_mV, mqn_m3Pr, atoms%lmaxd)
                  do imesh = 1, atoms%jri(itype)
                    vEff0Gaunt = fSh(imesh, lmV) * gauntFactor
                    fShNsphVarphi(1, imesh, lm1Pr, lmp) = fShNsphVarphi(1, imesh, lm1Pr, lmp) + vEff0Gaunt &
                                                                                        & * grVarphiCh1(imesh, ichan, lmp, mqn_m2Pr)
                    fShNsphVarphi(2, imesh, lm1Pr, lmp) = fShNsphVarphi(2, imesh, lm1Pr, lmp) + vEff0Gaunt &
                                                                                        & * grVarphiCh2(imesh, ichan, lmp, mqn_m2Pr)
                  end do ! imesh
                end do ! oqn_l
              end do ! mqn_mV
            end do ! oqn_lV
          end do ! ichan
        end do ! iradf
      end do ! mqn_m
    end do ! oqn_l

  end subroutine calcFnsphGrVarphi

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Calculates the action of a Hamiltonian, of the non-spherical part of the effective unperturbed potential and of the gradient of
  !> the non-spherical effective potential onto basis functions.
  !>
  !> @details
  !> The action of the potentials onto the basis functions are needed for the calculation of the last term within A.53 PhD thesis
  !> Aaron Klueppelberg (PhDAK). The action of the Hamiltonian onto the basisfunctions is needed severaltimes within A.50, A.51 and
  !> A.53
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
      & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi, r2grVeff0Varphi )

    use m_types, only : t_atoms, t_sphhar
    use m_jpConstants, only : c_im, c_mi, iu
    use mod_juPhonUtils, only : calcGrR2FinLH

!    use mod_juPhonUtils, only : fopen, fclose
    implicit none

    ! Type parameter
    type(t_atoms),        intent(in)  :: atoms
    type(t_sphhar),       intent(in)  :: lathar

    ! Scalar parameter
    integer,              intent(in)  :: itype
    integer,              intent(in)  :: iatom
    integer,              intent(in)  :: lmpMax

    ! Array parameter
    real,                 intent(in)  :: El(:, 0:, :, :)
    real,                 intent(in)  :: varphi1(:,:,0:)
    real,                 intent(in)  :: varphi2(:,:,0:)
    integer,              intent(in)  :: nRadFun(0:, :)
    complex,              intent(in)  :: vEff0MtSpH(:, :)
    real,                 intent(in)  :: vEff0MtLh(:, 0:, :)
    complex,              intent(in)  :: clnu_atom(:, 0:, :)
    integer,              intent(in)  :: nmem_atom(0:, :)
    integer,              intent(in)  :: mlh_atom(:, 0:, :)
    real,                 intent(in)  :: grVarphiCh1(:, :, :, -1:)
    real,                 intent(in)  :: grVarphiCh2(:, :, :, -1:)
    integer,              intent(in)  :: grVarphiChLout(:, 0:)
    integer,              intent(in)  :: grVarphiChMout(-atoms%lmaxd:, -1:)
    complex,              intent(out) :: hVarphi(:, :, 0:, :)
    complex,              intent(out) :: vEff0NsphGrVarphi(:, :, :, :, -1:)
    complex,              intent(out) :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,              intent(out) :: r2grVeff0Varphi(:, :, :, :, :)



    ! Scalar variables
    integer                           :: lmp
    integer                           :: oqn_l
    integer                           :: mqn_m
    integer                           :: oqn_l1Pr
    integer                           :: mqn_m1Pr
    integer                           :: imesh
    integer                           :: iradf
    real                              :: rInv
    integer                           :: idir
    integer                           :: lm
    integer                           :: lm_pre
    integer                           :: lm1Pr
    integer                           :: lm1Pr_pre
    integer                           :: irel
    integer                           :: mqn_m2PrC
    integer                           :: ilh

    ! Array variables
    real,           allocatable       :: hSph(:, :, :)
    real,           allocatable       :: r2Veff0SphLh(:, :, :)
    real,           allocatable       :: r2Veff0Lh(:, :, :)
    complex,        allocatable       :: r2grVeff0SphSh( :, :, :, : )
    complex,        allocatable       :: r2grVeff0Sh( :, :, :, : )
    real,           allocatable       :: recMesh(:)
    complex,        allocatable       :: grVeff0ShRobust(:, :, :)
    complex,        allocatable       :: vnsphEff0Varphi(:, :, :, :)

!    real,           allocatable       :: vxc0mt(:, :, :, :)
!
    allocate( hSph(2, atoms%jmtd, lmpMax), vnsphEff0Varphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )

    hSph(:, :, :) = 0.
    vnsphEff0Varphi = cmplx(0., 0.)

    ! Spherical part only contributes until lmax and not lmax + 2 as we have only the ket here expanded until lmax and in the overlap we
    ! cut there at lmax.
    ! It should be the best compromise to close the loops here over oqn_l, mqn_m and p otherwise we have redundant code or if-clauses in
    ! the inner loops. With this solution we only count the loops of oqn_l, mqn_m and p twice.
    ! If we would have run the p loop from 1, this had led to if-clauses in the inner loops.
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      do mqn_m = -oqn_l, oqn_l
        lmp = lmp + 1
        do imesh = 1, atoms%jri(itype)
          hSph(1, imesh, lmp) = El(1, oqn_l, itype, 1) * varphi1(imesh, 1, oqn_l)
          hSph(2, imesh, lmp) = El(1, oqn_l, itype, 1) * varphi2(imesh, 1, oqn_l)
        end do ! imesh
        lmp = lmp + 1
        do imesh = 1, atoms%jri(itype)
          hSph(1, imesh, lmp) = varphi1(imesh, 1, oqn_l) + El(1, oqn_l, itype, 1) * varphi1(imesh, 2, oqn_l)
          hSph(2, imesh, lmp) = varphi2(imesh, 1, oqn_l) + El(1, oqn_l, itype, 1) * varphi2(imesh, 2, oqn_l)
        end do ! imesh
        do iradf = 3, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          do imesh = 1, atoms%jri(itype)
            hSph(1, imesh, lmp) = El(iradf - 1, oqn_l, itype, 1) * varphi1(imesh, iradf, oqn_l)
            hSph(2, imesh, lmp) = El(iradf - 1, oqn_l, itype, 1) * varphi2(imesh, iradf, oqn_l)
          end do ! imesh
        end do ! p
      end do ! mqn_m
    end do ! oqn_l

    ! Action of normal non-spherical potential onto MT basis functions without basis matching coefficients
!    call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, vEff0MtSpH, vnsphEff0Varphi )
    call CalcFnsphVarphi( atoms, itype, 1, nRadFun, varphi1, varphi2, vEff0MtSpH, vnsphEff0Varphi )

    ! Action of complete Hamiltonian onto MT basis functions without basis matching coefficients
    lmp = 0
    do oqn_l = 0, atoms%lmax(itype)
      lm_pre = oqn_l * (oqn_l + 1)
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
        do iradf = 1, nRadFun(oqn_l, itype)
          lmp = lmp + 1
          ! Add action of spherical Hamiltonian to the diagonal of lmp and l'm'p'.
          do imesh = 1, atoms%jri(itype)
            do irel = 1, 2
              hVarphi(irel, imesh, lm, lmp) = hSph(irel, imesh, lmp)
            end do ! irel
          end do ! imesh
          ! We store the bra index until lmax + 2 because the double gradient of the radial solution is expanded maximally until
          ! lmax + 2
          !do oqn_l1Pr = 0, atoms%lmax(itype) + 2
          do oqn_l1Pr = 0, atoms%lmax(itype) !+ 2
            lm1Pr_pre = oqn_l1Pr * (oqn_l1Pr + 1)
            do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
              lm1Pr = lm1Pr_pre + mqn_m1Pr
              do imesh = 1, atoms%jri(itype)
                do irel = 1, 2
                  ! The factors i^l actually assigned to the matching coefficients in theoretical equations are excluded from abcof
                  ! and abcof3. Therefore, we multiply them here, so that, we do not need to care for them later.
                  hVarphi(irel, imesh, lm1Pr, lmp) = hVarphi(irel, imesh, lm1Pr, lmp) + vnsphEff0Varphi(irel, imesh, lm1Pr, lmp)
!                  hVarphi(irel, imesh, lm1Pr, lmp) = vnsphEff0Varphi(irel, imesh, lm1Pr, lmp)
                end do ! irel
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! iradf
      end do ! mqn_m1Pr
    end do ! oqn_l1Pr
    deallocate(hSph)

    ! Construction of the special non-spherical potential required for the calculation of <gradVarphi|H|gradVarphi>

    ! Remove numerically critical singularity at the core
    ! !todo test this!
    ! allocate( recMesh(atoms%jmtd) )
    ! recMesh = 0.
    ! do imesh = 1, atoms%jri(itype)
    !   rInv = 1. / atoms%rmsh(imesh, itype)
    !   recMesh(imesh) =  rInv / atoms%rmsh(imesh, itype)
    !   r2VeffSpLh(imesh, 0) = r2VeffSpLh(imesh, 0) - atoms%zatom(itype) * rInv
    ! end do ! imesh
    !!call CalcGrFLhNatAt(atoms, lathar, itype, iatom, clnu_atom, nmem_atom, mlh_atom, vEff0MtLh(:, :, itype), grVeff0ShRobust)

    ! ! The analytical gradient 1 / r of the l = 0 channel is expanded in spherical harmonics and added to the robust numerically derived
    ! ! gradient part in the l = 1 channel
    ! ! todo c_mi is the transposed vesion of c_im!
    ! do idir = -1, 1
    !   do mqn_m = -1, 1
    !     ! l * (l + 1) = 2
    !     lm = 2 + mqn_m
    !     do imesh = 1, atoms%jri(itype)
    !       grVeff0ShRobust(imesh, lm, idir) = grVeff0ShRobust(imesh, lm, idir) - recMesh(imesh) * c_mi(mqn_m + 2, idir + 2)
    !     end do ! imesh
    !   end do ! mqn_m
    ! end do ! idir
    ! deallocate(recMesh)

    allocate(r2Veff0SphLh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype))
    r2Veff0SphLh(:, :, :) = cmplx(0.0, 0.0)

!    allocate( vXC0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
!    vxc0mt(:, :, :, :) = cmplx(0., 0.)
!    ! XC potential in the MT read in from FLEUR
!    call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
!    read(1000) vXC0MT(:, :, :, :)
!    call fclose(1000)
!    do ilh = 0, lathar%nlhd
    do imesh = 1, atoms%jri(itype)
      r2Veff0SphLh(imesh, 0, itype) = vEff0Mtlh(imesh, 0, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      !r2Veff0SphLh(imesh, ilh, itype) = vEff0Mtlh(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = vxc0mt(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = (vEff0Mtlh(imesh, ilh, itype) - vxc0mt(imesh, ilh, itype, 1)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
    end do
!    end do

    allocate( r2grVeff0SphSh( atoms%jmtd, (atoms%lmax(itype) + 2)**2, atoms%nat, 3 ) )
    r2grVeff0SphSh = cmplx(0.0, 0.0)

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0SphLh, r2grVeff0SphSh )

    allocate(r2Veff0Lh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype))
    r2Veff0Lh(:, :, :) = cmplx(0.0, 0.0)
    do ilh = 0, lathar%nlhd
      do imesh = 1, atoms%jri(itype)
        r2Veff0Lh(imesh, ilh, itype) = vEff0Mtlh(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh
    end do ! ilh

    allocate( r2grVeff0Sh( atoms%jmtd, (atoms%lmax(itype) + 2)**2, atoms%nat, 3 ) )
    r2grVeff0Sh = cmplx(0.0, 0.0)
    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0Lh, r2grVeff0Sh )


  !  do idir = 1, 3
  !    do oqn_l = 0, atoms%lmax(itype) + 1
  !      do mqn_m = -oqn_l, oqn_l
  !        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
  !        do imesh = 1, atoms%jmtd
  !          write(4000, '(4i8,2f15.8)') idir, oqn_l, mqn_m, imesh, r2grVeff0SphSh(imesh, lm, 1, idir)
  !        end do ! imesh
  !      end do ! mqn_m
  !    end do ! oqn_l
  !  end do ! idir


    !todo before commnit: is it critical that lmPr is running until lmax + 2
    ! NOTE: The main difference comes from the Coulomb potential, because here we have a difference in Weinert or numerical, in the
    ! xc potential, we both do a numerical gradient so do not feel a big difference between chain rule to the kernel or direct 
    ! numerical gradient. In the Weinert method, for Neon we only have a contribuiton of 1e-6 from the non-spherical contributions of
    ! the potential. The main difference comes from the spherical component. Here, we only have a difference of 1e-4 between using
    ! Weinert or not.
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2grVeff0SphSh(:, :, iatom, idir), &
                                                                                            & r2grVeff0SphVarphi(:, :, :, :, idir) )
    end do

    r2grVeff0Varphi(:, :, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2grVeff0Sh(:, :, iatom, idir), &
                                                                                            & r2grVeff0Varphi(:, :, :, :, idir) )
    end do

    vEff0NsphGrVarphi = cmplx(0., 0.)
    do mqn_m2PrC = -1, 1
      call calcFnsphGrVarphi( atoms, itype, 1, mqn_m2PrC, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, vEff0MtSpH, &
                                                                                        & vEff0NsphGrVarphi(:, :, :, :, mqn_m2PrC) )
!      call calcFnsphGrVarphi( atoms, itype, 0, mqn_m2PrC, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, vEff0MtSpH, &
!                                                                                        & vEff0NsphGrVarphi(:, :, :, :, mqn_m2PrC) )
    end do ! mqn_m2Pr

    !todo should be up to lmax + 1
    !todo put in the generation of the spherical potential in here
  end subroutine CalcHnGrV0Varphi

  ! For q = 0 the sum of the IR Hellmann-Feynman integrals contributing is this surface integral
  subroutine CalcSurfIntIRDynMat( atoms, cell, ngdp1, ngdp2, gdp1, gdp2, rho0IRpw, grVext0IR, qpoint, surfInt )

    use m_types
    use m_jPConstants, only : iu, tpi, fpi, c_im
    use m_ylm_old
    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),        intent(in)  :: atoms
    type(t_cell),         intent(in)  :: cell

    ! Scalar parameter
    integer,              intent(in)  :: ngdp1
    integer,              intent(in)  :: ngdp2

    ! Array parameter
    integer,              intent(in)  :: gdp1(:, :)
    integer,              intent(in)  :: gdp2(:, :)
    complex,              intent(in)  :: rho0IRpw(:)
    complex,              intent(in)  :: grVext0IR(:, :)
    real,                 intent(in)  :: qpoint(:)
    complex,              intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                           :: idirC
    integer                           :: idirR
    integer                           :: iatom
    integer                           :: itype
    integer                           :: ieqat
    integer                           :: iG
    integer                           :: it
    integer                           :: iGp
    complex                           :: phaseFac
    complex                           :: tSummedCitY1t
    complex                           :: pref
    complex                           :: surfIntNonMat

    ! Array variables
    real                              :: gSum(3)
    real                              :: gSumCart(3)
    complex                           :: ylm(4)
    real                              :: sbes(0:1)

    surfInt(:, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      pref = fpi * iu * atoms%rmt(itype) * atoms%rmt(itype)
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do iG = 1, ngdp1
          do iGp = 1, ngdp2

            gSum(1:3) = gdp1(1:3, iG) + gdp2(1:3, iGp) + qpoint(1:3)
            gSumCart(1:3) = matmul( cell%bmat(1:3, 1:3), gSum(1:3) )

            ylm(:) = cmplx(0., 0.)
            call ylm4( 1, gSumCart, ylm )

            sbes(:) = 0
            call sphbes(1, norm2(gSumCart) * atoms%rmt(itype), sbes)

            phaseFac = exp( iu * tpi * dot_product(gSum(1:3), atoms%taual(1:3, iatom)))

            surfIntNonMat = pref * phaseFac * sbes(1) * rho0IRpw(iG)

            do idirC = 1, 3
              do idirR = 1, 3
                ! Corresponds to the magnetic quantum number m for the orbital quantum number l = 1
                tSummedCitY1t = cmplx(0., 0.)
                do it = -1, 1
                  tSummedCitY1t = tSummedCitY1t + c_im(idirR, it + 2) * ylm(it + 3)
                end do ! it

                surfInt(idirR, idirC) = surfInt(idirR, idirC) - surfIntNonMat * tSummedCitY1t * grVext0IR(iGp, idirC)
              end do ! idirR
            end do ! idirC
          end do ! iGp
        end do ! iG
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntIRDynMat

  subroutine CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MT, surfInt)

    use m_types, only : t_atoms, t_sphhar
    use m_gaunt, only : Gaunt1
    use m_jpConstants, only : c_im

    implicit none

    ! Type parameters
    type(t_atoms),                     intent(in)  :: atoms
    type(t_sphhar),                    intent(in)  :: lathar

    ! Array parameters
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    real,                              intent(in)  :: rho0MT(:, 0:, :)
    complex,                           intent(in)  :: grVext0MT(:, :, :, :)
    complex,                           intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                                        :: iatom
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: ptsym
    integer                                        :: ilh
    integer                                        :: oqn_l
    integer                                        :: imem
    integer                                        :: mqn_m
    integer                                        :: mqn_mp
    integer                                        :: oqn_lpp
    integer                                        :: mqn_mpp
    integer                                        :: lmpp
    real                                           :: gauntFactor
    integer                                        :: idirC
    integer                                        :: idirR

    surfInt(:, :) = 0
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            do mqn_mp = -1, 1
              ! lmax + 1 is okay for grVext1
             ! do oqn_lpp = abs(oqn_l - 1), oqn_l + 1
              do oqn_lpp = 0, atoms%lmax(itype)
                do mqn_mpp = -oqn_lpp, oqn_lpp
             !   mqn_mpp = mqn_m + mqn_mp
                lmpp = oqn_lpp * (oqn_lpp + 1) + 1 + mqn_mpp
!                gauntFactor = Gaunt1( oqn_lpp, oqn_l, 1, mqn_mpp, mqn_m, mqn_mp, atoms%lmax(itype) + 1)
                gauntFactor = Gaunt1( oqn_lpp, oqn_l, 1, mqn_mpp, mqn_m, mqn_mp, atoms%lmax(itype))
                do idirC = 1, 3
                  do idirR = 1, 3
                  !todo check whether rho0MT or and c_im should be conjugated or grVext0MT should e conjugatd. This current solution gives the best results.
                    surfInt(idirR, idirC) = surfInt(idirR, idirC) + c_im(idirR, mqn_mp + 2)  * atoms%rmt(itype)**2         &
                      & * rho0MT(atoms%jri(itype), ilh, itype) * clnu_atom(imem, ilh, iatom) *                           &
                      & conjg(grVext0MT(atoms%jri(itype), lmpp, idirC, iatom)) * gauntFactor
                  end do ! idirC
                end do ! idirR
                end do
              end do ! oqn_lpp
            end do ! mqn_mp
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntMTDynMat


end module m_jpSetupDynMatHelper
