module m_jpSetupDynMatDeprecated

  implicit none

  !>-------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst
  !>
  !> @brief
  !> Calculates the core Pulay contribution to the dynamical matrix given in equation 7.115 PhD thesis Aaron Klueppelberg (PhDAK)
  !>
  !> @details
  !> This routine reads the core density from cdnc file written by fleur and multiplies it with the non-spherical part of the
  !> unperturbed potential according to 7.115 PhDAK. This contribution only matters if the effective potential has non-vanishing
  !> coefficients in the l = 2 channel.
  !>
  !> @note
  !> Note: One could make this routine faster by pulling out the Gaunt coefficients from the gradient routine and using it after the
  !> calculation of the radial integrals. In this way one would only have to calculate the radial integrals for every l and then
  !> multiply the Gaunt coefficients. This would be more efficient but we choose the traditional way to minimize the error
  !> probability. After this programm is running one has a good benchmark.
  !>
  !> @attention
  !> cndc file required!
  !>
  !> @todo
  !> Account for the spin.
  !>-------------------------------------------------------------------------------------------------------------------------------------
  subroutine EvalPuHepsBraKetsCore( atoms, lathar, input, clnu_atom, nmem_atom, mlh_atom, vEff0MT, dynMatPu )

    use m_types
    use mod_juPhonUtils, only : fopen, fclose
    use m_intgr, only : intgr0, intgr3NoIntp
    use m_jpPotDensHelper, only : calcGrR2FinLH
    use m_JPConstants, only : fpi

    implicit none

    ! Type parameter
    type(t_atoms),               intent(in)    :: atoms
    type(t_sphhar),              intent(in)    :: lathar
    type(t_input),               intent(in)    :: input

    ! Array parameter
    complex,                     intent(in)    :: clnu_atom(:, 0:, :)
    integer,                     intent(in)    :: nmem_atom(0:, :)
    integer,                     intent(in)    :: mlh_atom(:, 0:, :)
    real,                        intent(in)    :: vEff0MT(:, 0:, :)
    complex,                     intent(inout) :: dynMatPu(:, :)


    ! Local scalars
    integer                                    :: ispin
    integer                                    :: itype
    integer                                    :: iatom
    integer                                    :: imesh
    real                                       :: qcore
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm_pre
    integer                                    :: lm
    integer                                    :: ilh
    integer                                    :: idir
    real                                       :: intgrReal
    real                                       :: intgrImag

    ! Local arrays
    real,           allocatable                :: rhoCoreRead( :, :, : )
    real,           allocatable                :: r2rho0Core( :, :, :, : )
    complex,        allocatable                :: grVeff0MT(:, :, :, :)
    complex,        allocatable                :: r2grRho0Core(:, :, :, :)
    complex,        allocatable                :: grRho0Core(:, :, :, :)
    real,           allocatable                :: r2(:, :)
    real,           allocatable                :: fReal(:)
    real,           allocatable                :: fImag(:)
!    complex,        allocatable                :: vEff0MTnSph(:, :, :)
    real,           allocatable                :: r2Veff0MT(:, :, :, :)
    complex,        allocatable                :: r2GrVeff0MT(:, :, :, :)


    allocate( fReal(atoms%jmtd) )
    allocate( fImag(atoms%jmtd) )
    allocate( r2rho0Core( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins ) )
    allocate( rhoCoreRead( atoms%jmtd, atoms%ntype, input%jspins ) )
    allocate( r2(atoms%jmtd, atoms%ntype))

    rhoCoreRead(:, :, :) = 0.
    ! Read in core density from cdnc
    call fopen(1000, name='cdnc', status='old', action='read', form='unformatted')
    do ispin = 1, input%jspins
      if ( ispin == 1 ) then
        rewind( 1000 )
      end if
      do itype = 1, atoms%ntype
        ! read in core density
        read( 1000 ) (rhoCoreRead( imesh, itype, ispin ), imesh=1, atoms%jri( itype ))

        ! core density is currently multiplied by r**2*sqrt(4*pi) such that intrg0 gives correct core charge
        qcore = 0.0
        call intgr0(rhoCoreRead, atoms%rmsh(1,itype), atoms%dx(itype),atoms%jri(itype), qcore)
        write (*, *) 'Core charge: ', qcore

        ! skip kinetic enrgy of the core
        read( 1000 )

      end do
      read( 1000 )
    end do
    call fclose(1000)
    !todo cored wie man richtig die core dichte ausließt

    !todo we should think about the r within the potential here!
    ! Construct full density
    r2Rho0Core(:, :, :, :) = 0.
    do ispin = 1, input%jspins
      do itype = 1, atoms%ntype
        do imesh = 1, atoms%jri(itype)
          r2(imesh, itype) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          r2Rho0Core(imesh, 0, itype, ispin) = r2Rho0Core(imesh, 0, itype, ispin) + rhoCoreRead(imesh, itype, ispin) / sqrt(fpi)
        end do
      end do
    end do

    ! Generate gradient of unperturbed effective potential
    ! Here to
    ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
    ! subtraction of small numbers close to the core.
    allocate( r2Veff0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Veff0MT(:, :, :, :) = 0.
    grVeff0MT(:, :, :, :) = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          r2Veff0MT(imesh, ilh, itype, 1) = veff0MT(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MT(:, :, :, 1), r2GrVeff0MT )
    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2rho0Core(:, :, :, 1), r2GrRho0Core)

    grVeff0MT(:, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                grVeff0MT(imesh, lm, iatom, idir) = r2GrVeff0MT(imesh, lm, iatom, idir) / r2(imesh, itype)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    allocate( grRho0Core(atoms%jmtd, ( atoms%lmaxd + 2)**2, atoms%nat, 3) )
    grRho0Core(:, :, :, :) = cmplx(0., 0.)
    do idir = 1, 3
      do ispin = 1, input%jspins
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do imesh = 1, atoms%jri(itype)
              grRho0Core(imesh, 2:4, iatom, idir) = r2GrRho0Core(imesh, 2:4, iatom, idir) / r2(imesh, itype)
            end do ! imesh
          end do ! ieqat
        end do ! itype
      end do ! ispin
    end do ! idir

!    allocate( vEff0MTnSph(atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
!    vEff0MTnSph(:, :, :) = cmplx(0., 0.)
!    vEff0MTnSph(:atoms%jmtd, 0:lathar%nlhd, :atoms%ntype) = vEff0MT(:atoms%jmtd, 0:lathar%nlhd, :atoms%ntype )


    ! Calculate the gradient of the full density
    !todo gradient before or after summation over bands and k-points because cndc is already summed?
!    write(*, *) 'Use the routine which calculates r2the gradient'
!    call calcGrFinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0Core(:, :, :, 1), grRho0Core )


!    write(*, *) 'Use the routine which calculates r2the gradient'
    ! Calculate the gradient of potential
!    call calcGrFinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, vEff0MTnSph, grVeff0MT )

    !NOTE: The integral only has a contribution if the non-spherical potential has a l = 2 compoenent, whose gradient scatters into
    ! the l = 1 channel, so that the multiplication with the core density has a contribution as this has only a l = 0 channel
    ! scattering into the l = 1 channel.

    ! Integration given in the last line of 7.115
    do idirC = 1, 3
      do idirR = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                fReal(:) = 0.0
                fImag(:) = 0.0
                do imesh = 1, atoms%jri(itype)
                  fReal(imesh) = real ( conjg(grRho0Core(imesh, lm, iatom, idirR)) * grVeff0MT(imesh, lm, iatom, idirC) * r2(imesh, itype) )
                  fImag(imesh) = aimag( conjg(grRho0Core(imesh, lm, iatom, idirR)) * grVeff0MT(imesh, lm, iatom, idirC) * r2(imesh, itype) )
                end do ! imesh
                call intgr3NoIntp(fReal, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrReal)
                call intgr3NoIntp(fImag, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), intgrImag)
                dynMatPu(idirR + (iatom - 1) * 3, idirC + (iatom - 1) * 3) =  &
                  & dynMatPu(idirR + (iatom - 1) * 3, idirC + (iatom - 1) * 3) - cmplx(intgrReal, intgrImag)
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idirR
    end do ! idirC

  end subroutine EvalPuHepsBraKetsCore

end module m_jpSetupDynMatDeprecated


