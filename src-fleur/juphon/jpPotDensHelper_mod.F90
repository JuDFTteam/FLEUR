!-----------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: JuPhon Potential and Density Helper Module
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> This module contains all auxiliary routines related to the calculation of potential and density.
!>
!> @details
!> There are routines which are needed for the calculation of both the (un)perturbed density and the (un)perturbed potential. This
!> module is intended as a helper module which contains such routines.
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpPotDensHelper.pdf'>document</a>.
!-----------------------------------------------------------------------------------------------------------------------------------
module m_jpPotDensHelper

  implicit none

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Converts expansion coefficient of star expansion into expansion coefficients of plane-wave coefficients for interstitial
  !> unperturbed density or potential.
  !>
  !> @details
  !> See als Equation 7.4 (dissertation CRG)
  !>
  !> @param[in]  stars       : Stars type, see types.f90
  !> @param[in]  ngdp        : Number of G-vectors for potentials and densities
  !> @param[in]  stExpandQ   : Expansion coefficients of quantity expanded in stars
  !> @param[in]  gdp         : G-vectors of potentials and densities
  !> @param[out] gVecExpandQ : Expansion coefficients of quantity expanded in plane waves
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine convertStar2G(stExpandQ, gVecExpandQ, stars, ngdp, gdp, cell)

    use m_types
    use m_jpConstants, only : compPhon

    implicit none

    ! Scalar Arguments
    type(t_stars),  intent(in)   :: stars
    integer,        intent(in)   :: ngdp

    ! Array Arguments
    complex,        intent(in)   :: stExpandQ(:)
    integer,        intent(in)   :: gdp(:, :)
    complex,        intent(out)  :: gVecExpandQ(ngdp)
    type(t_cell), optional, intent(in) :: cell

    ! Local Scalar Variables
    integer                      :: iGvec
    real                         :: Gext(3)
    gVecExpandQ = 0
    do iGvec = 1, ngdp
      gVecExpandQ(iGvec) = stExpandQ(stars%ig(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec))) * stars%rgphs(gdp(1, iGvec), gdp(2, iGvec), gdp(3, iGvec))
      if (compPhon) then
        write(108,*) iGvec       
        if (present(cell)) then
          Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iGvec))
          write(108,*) Gext(1), Gext(2), Gext(3)
        end if
        write(109,*) iGvec
        write(109,*) real(gVecExpandQ(iGvec)), aimag(gVecExpandQ(iGvec))
      end if
    end do
  end subroutine convertStar2G

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
  !>
  !> @details
  !> See Equation 7.5 (dissertation CRG)
  !>
  !> @param[in]  atoms      : Atoms type, see types.f90
  !> @param[in]  lathar     : Lattice harmonics type, see types.f90.
  !> @param[in]  mlh_atom   : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[in]  nmem_atom  : Number of lattice harmonic members for every atom.
  !> @param[in]  clnu_atom  : Phase mediating between stars and plane waves.
  !> @param[out] spHarECoef : Expansion coefficients of quantity expanded in spherical harmonics
  !>
  !> @todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
  !        maybe construct a pointer and run only over them to make it memory efficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine convertLH2SphHarm(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, lHarECoef, spHarECoef)

    use m_types
    use m_jpConstants, only : compPhon

    implicit none

    ! Scalar Argument
    type(t_atoms),               intent(in)  :: atoms
    type(t_sphhar),              intent(in)  :: lathar

    ! Array Arguments
    integer,                     intent(in)  :: mlh_atom(:,0:,:)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    complex,                     intent(in)  :: clnu_atom(:,0:,:)
    real,                        intent(in)  :: lHarECoef(:, 0:, :, :)
    complex,                     intent(out) :: spHarECoef(:, :, :, :)

    ! Local Scalar Variables
    integer                                  :: itype
    integer                                  :: ieqat
    integer                                  :: iatom
    integer                                  :: ptsym
    integer                                  :: ilh
    integer                                  :: oqn_l
    integer                                  :: lm_pre
    integer                                  :: imem
    integer                                  :: mqn_m
    integer                                  :: lm
    integer                                  :: imesh

    spHarECoef(:, :, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            do imesh = 1, atoms%jri(itype)
              spHarECoef(imesh, lm, iatom, 1) = spHarECoef(imesh, lm, iatom, 1) + lHarECoef(imesh, ilh, itype, 1) &
                                                                                                 & * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

    if (compPhon) then
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do lm = 1, atoms%lmaxd*(atoms%lmaxd+2) + 1
            do imesh = 1, atoms%jri(itype)
              write(109,*) imesh, lm, iatom
              write(109,*) real(spHarECoef(imesh, lm, iatom, 1)), aimag(spHarECoef(imesh, lm, iatom, 1))
            end do ! imesh
          end do ! lm
        end do ! ieqat
      end do ! itype
    end if

  end subroutine convertLH2SphHarm

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine generates the G-vectors which are used for constructing the potential and the density.
  !>
  !> @details
  !> They are fulfilling |G| < Gmax
  !>
  !> @param[in]  starsT   : Contains stars-related quantities; definition of its members in types.F90 file.
  !> @param[in]  cellT    : Contains unit-cell related quantities; definition of its members in types.F90 file.
  !> @param[out] input    : Input type, see types.f90.
  !> @param[out] ngpd     : Number of G-vectors used for density and potential.
  !> @param[out] ngdp2km  : Number of G-vectors for potentials and densities which are smaller than 2 kmax.
  !> @param[out] gdp      : Contains G-vectors in internal coordinates used for density and potential, see for dimensions.
  !> @param[out] gdp2Ind  : Stores the index of a potential and density G-Vector. The dimensions are the G-vector components.
  !> @param[out] gdp2iLim : Stores the min and maxvals of gdp2Ind
  !> @param[in]  testMode : Indexarray stores all G-vectors until Gmax instead of 2kmax if set .true.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine genPotDensGvecs(starsT, cellT, inputT, ngdp, ngdp2km, gdp, gdp2Ind, gdp2iLim, testMode )

    use m_types
    use m_juDFT

    implicit none

    ! Variables:
    ! Gx      : loop variable to run over all possible x-components of an internal G-vector
    ! Gy      : loop variable to run over all possible y-components of an internal G-vector
    ! Gz      : loop variable to run over all possible z-components of an internal G-vector
    ! gdptemp : auxiliary array for intermediate storage of accepted G-vectors
    ! Gint    : temporary variable to store current candidate for a G-vector in internal coordinates
    ! Gext    : temporary variable to store current candidate for a G-vector in external coordinates

    ! Scalar parameters
    type(t_stars),          intent(in)   :: starsT
    type(t_cell),           intent(in)   :: cellT
    type(t_input),          intent(in)   :: inputT

    integer,                intent(out)  :: ngdp
    integer,                intent(out)  :: ngdp2km
    logical,                intent(in)   :: testMode

    ! Array parameter
    integer,  allocatable,  intent(out)  :: gdp(:, :)
    integer,  allocatable,  intent(out)  :: gdp2Ind(:, :, :)
    integer,                intent(out)  :: gdp2iLim(2, 3)

    ! Local scalar variables
    integer                              :: Gx, Gy, Gz, iG
    integer                              :: ngrest

    ! Local array variables
    integer                              :: gdptemp2kmax(3, (2 * starsT%k1d + 1) * (2 * starsT%k2d + 1) * (2 * starsT%k3d +  1))
    integer                              :: gdptemprest(3, (2 * starsT%k1d + 1) * (2 * starsT%k2d + 1) * (2 * starsT%k3d +  1))
    integer                              :: Gint(3)
    real                                 :: Gext(3)

    ! From all possible G-vectors in a box, only these are accepted which are element of a sphere with radius gmax. As benchmark,
    ! it is checked, whether the current G-vector candidate is contained within any star
    ngdp = 0
    gdptemp2kmax = 0
    gdptemprest = 0
    ngdp2km = 0
    ngrest = 0
    ! The idea to sort the G-vectors until 2kmax before all other Gvectors stems from M. Betzinger
    do Gx = -starsT%k1d, starsT%k1d
      do Gy = -starsT%k2d, starsT%k2d
        do Gz = -starsT%k3d, starsT%k3d
          Gint = [Gx, Gy, Gz]
          Gext =  matmul(cellT%bmat, Gint) !transform from internal to external coordinates
          if (norm2(Gext) <= starsT%gmaxInit) then
#ifdef DEBUG_MODE
            if (starsT%ig(Gx, Gy, Gz) <= 0) then
              call juDFT_error('Inconsistency in determination of G-vectors for potential or density', calledby='genPotDensGvecs', &
                 & hint='Check whether selection methods correct.', file='jpPotDens.F90', line=52)
            endif
#endif
            ngdp = ngdp + 1
            ! Sort G-vectors
            if ( norm2(Gext) <= 2 * inputT%rkmax ) then
              ngdp2km = ngdp2km + 1
              gdptemp2kmax(:, ngdp2km) = Gint(:)
            else
              ngrest = ngrest + 1
              gdptemprest(:, ngrest) = Gint(:)
            end if
          endif
        enddo !Gz
      enddo !Gy
    enddo !Gx
    allocate(gdp(3, ngdp))
    gdp(:, :ngdp2km) = gdptemp2kmax(:, :ngdp2km)
    gdp(:, ngdp2km + 1 : ngdp) = gdptemprest(:, :ngrest)


    ! Create mapping array from G-vector to G-vector index up to Gmax
    gdp2iLim(1, 1) = minval(gdp(1, :))
    gdp2iLim(2, 1) = maxval(gdp(1, :))
    gdp2iLim(1, 2) = minval(gdp(2, :))
    gdp2iLim(2, 2) = maxval(gdp(2, :))
    gdp2iLim(1, 3) = minval(gdp(3, :))
    gdp2iLim(2, 3) = maxval(gdp(3, :))
    allocate(gdp2Ind(gdp2iLim( 1, 1) : gdp2iLim( 2, 1), gdp2iLim( 1, 2) : gdp2iLim( 2, 2), gdp2iLim( 1, 3) : gdp2iLim( 2, 3)))
    gdp2Ind = 0
    if (testMode) then
      do iG = 1, ngdp2km
        gdp2Ind(gdptemp2kmax(1, iG), gdptemp2kmax(2, iG), gdptemp2kmax(3, iG)) = iG
      end do
    else
      do iG = 1, ngdp
        gdp2Ind(gdptemp2kmax(1, iG), gdptemp2kmax(2, iG), gdptemp2kmax(3, iG)) = iG
      end do
    end if

  end subroutine genPotDensGvecs

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine generates the G-vectors which are used for constructing the linear variation of the density and potential in
  !> the interstitial.
  !>
  !> @details
  !> The G-vectors are fulfilling |G + q| < Gmax
  !>
  !> @param[in]  stars      : Stars type, see types.f90
  !> @param[in]  cell       : Unit cell type, see types.f90
  !> @param[in]  input      : Input type, see types.f90.
  !> @param[in]  ngpqdp     : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  ngpqdp2km  : Number of G-vectors for shifted G-set for q-point with index iqpt which are smaller than 2 kmax.
  !> @param[in]  qpoint     : Qpoint in internal coordinates as in qpts%bk(1:3, iqpt)
  !> @param[in]  gpqdp      : G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  gpqdp2Ind  : Stores the index of a shifted G-Vector for q-point with index iqpt. The dimensions are the G-vector
  !>                          components.
  !> @param[out] gpqdp2iLim : Stores the min and maxvals of gpqdp2Ind
  !>
  !> @note this is deprecated
  !---------------------------------------------------------------------------------------------------------------------------------
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

    allocate( gpqdptemp2kmax(3, (2 * stars%k1d + 1) * (2 * stars%k2d + 1) * (2 * stars%k3d +  1)), &
            & gpqdptemprest(3, (2 * stars%k1d + 1) * (2 * stars%k2d + 1) * (2 * stars%k3d +  1)) )

    ngpqdp = 0
    ngpqdp2km = 0
    ngrest = 0
    gpqdptemp2kmax(:, :) = cmplx(0., 0.)
    gpqdptemprest(:, :) = cmplx(0., 0.)
    ! From all possible G-vectors in a box, only these are accepted which are element of a sphere with radius gmax which is shifted.
    ! We need a little bit more than k*d because they are thought for a Gmax ball that is not shifted by a q, i.e. |G+q|<Gmax
    do iGx = -(stars%k1d + 3), (stars%k1d + 3)
      do iGy = -(stars%k2d + 3), (stars%k2d + 3)
        do iGz = -(stars%k3d + 3), (stars%k3d + 3)
          Gint = [iGx, iGy, iGz]
          Gpqext =  matmul(cell%bmat, real(Gint(1:3) + qpoint(1:3))) !transform from internal to external coordinates
          if (norm2(Gpqext) <= stars%gmaxInit) then
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

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This routine enables the usage of a global coordinate system instead of local coordinate systems which are rotated for a
  !> non-representativ atom to fulfill the symmetry requirements.
  !>
  !> @details
  !>  If symmetry is considered  the non-representative atoms of an atom type are treated with a rotated coordinate system. Since
  !>  for moving atoms governed by a phonon mode either almost every or every symmetry is broken. So there is no need to use a
  !>  rotated coordinate system for symmetry reasons. Instead we use a global coordinate system. But quantities like the unperturbed
  !>  density coming from FLEUR are stored in these local coordinate systems. Hence this has to be transformed. The following method
  !>  calculates the new lattice-harmonics expansion coefficients and new properties of the rotated lattice harmonics.
  !>
  !>  @param[in]   atoms     : Atoms type, see types.f90
  !>  @param[in]   sym       : Symmetries type, see types.f90.
  !>  @param[in]   cell      : Unit cell type, see types.f90.
  !>  @param[in]   lathar    : Lattice harmonics type, see types.f90.
  !>  @param[out]  memd_atom : Maximal number of members in all lattice harmonics for all atoms.
  !>  @param[out]  clnu_atom : Phase mediating between stars and plane waves.
  !>  @param[out]  nmem_atom : Number of lattice harmonic members for every atom.
  !>  @param[out]  mlh_atom  : Magnetic quantum number m of lattice harmonic members for every atom.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine genRotLh(atomsT, symT, cellT, latharT, clnu_atom, memd_atom, nmem_atom, mlh_atom)

    use m_types
    use m_dwigner

    implicit none

    ! Variables:
    ! atomsT     : atoms type, see m_types for documentation
    ! symT       : sym type, see m_types for documentation
    ! cellT      : cell type, see m_types for documentation
    ! latharT    : lattice harmonic type, see m_types for documentation
    ! memd_atom  : maximal number of members of spherical harmonics in global coordinates
    ! clnu_atom  : phase factors in lattice harmonic expansion in global coordinates for every single atom
    ! nmem_atom  : number of members in lattice harmonic for every single atom
    ! mlh_atom   : ms of members of a lattice harmoinc for every single atom
    ! d_wgnJp    : Wigner matrices to rotate from local coordinate systems into global coordinate system.

    ! Scalar variables
    type(t_atoms),                 intent(in)  :: atomsT
    type(t_sym),                   intent(in)  :: symT
    type(t_cell),                  intent(in)  :: cellT
    type(t_sphhar),                intent(in)  :: latharT
    integer,                       intent(out) :: memd_atom

    !Array variables
    complex,          allocatable, intent(out) :: clnu_atom(:, :, :)
    integer,          allocatable, intent(out) :: nmem_atom(:, :)
    integer,          allocatable, intent(out) :: mlh_atom(:, :, :)

    ! Local Array variabel
    complex                                    :: d_wgnJp(-atomsT%lmaxd:atomsT%lmaxd, -atomsT%lmaxd:atomsT%lmaxd, 0:atomsT%lmaxd, &
                                                                                                                         & symT%nop)

    ! Generate Wigner matrices which rotate the lattice harmonics and communicate between the l and ms and the l and m's of the
    ! rotated and back-rotated system.

    call d_wigner(symT%nop, symT%mrot, cellT%bmat, atomsT%lmaxd, d_wgnJp(:, :, 1:, :symT%nop))
    d_wgnJp(:, :, 0, :) = 1

    memd_atom=1
    ! Generate new parameters * _atom for the rotated global coordinate system
    call rotate_clnu(latharT, symT, atomsT, d_wgnJp, clnu_atom, memd_atom, nmem_atom, mlh_atom)

  end subroutine genRotLh

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This method creates quantities of the lattice harmonics which alter when rotating back from the local coordinate systems to the
  !> global coordinate system.
  !>
  !> @details
  !>
  !>  @param[in]   latharT   : Lattice harmonics type, see types.f90.
  !>  @param[in]   symT      : Symmetries type, see types.f90.
  !>  @param[in]   atomsT    : Atoms type, see types.f90.
  !>  @param[out]  memd_atom : Maximal number of members in all lattice harmonics for all atoms.
  !>  @param[in]   d_wgn     : Wigner matrix mediating between rotated and not-rotated coordinate system
  !>  @param[out]  mlh_atom  : Magnetic quantum number m of lattice harmonic members for every atom.
  !>  @param[out]  nmem_atom : Number of lattice harmonic members for every atom.
  !>  @param[out]  clnu_atom : Phase mediating between stars and plane waves.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine rotate_clnu(latharT, symT, atomsT, d_wgn, clnu_atom, memd_atom, nmem_atom, mlh_atom)

    use m_JPConstants
    use m_ylm_old
    use m_types
    use m_juDFT

    implicit none

    ! Variables:
    !
    ! latharT    : lattice harmonic type, see m_types for documentation
    ! symT       : sym_Type, see m_types for documentation
    ! atomsT     : atoms_Type, see m_types for documentation
    ! memd_atom  : maximal number of members of lattice harmonics for every atom
    ! d_wgn      : Wigner matrices
    ! mlh_atom   : ms of members of the lattice harmonics for every atom
    ! nmem_atom  : number of members of the lattice harmonics for every atom
    ! clnu_atom  : phase factor for the lattice harmonic expansion for every atom
    ! iop        : number of symmetry operation which rotates back
    ! itype      : runs over atom types
    ! ieq        : runs over equivalent atoms of an atom type
    ! iatom      : runs over atoms
    ! ilh        : runs over lattice harmonics
    ! imem       : runs over members of a lattice harmonic
    ! itypsym    : point symmetry of iatom
    ! l          : orbital quantum number l
    ! m          : magnetic quantum number m
    ! ok         : stores error code after allocation of a new variable otherwise it is 0
    ! carr1      : stores phase factors clnu in local coordinate systems
    ! carr2      : stores phase factors clnu in global coordinate system

    ! - scalars -
    type(t_sphhar),       intent(in)  ::  latharT
    type(t_sym),          intent(in)  ::  symT
    type(t_atoms),        intent(in)  ::  atomsT
    integer,              intent(out) ::  memd_atom

    ! - arrays -
    complex,              intent(in)  ::  d_wgn(-atomsT%lmaxd:, -atomsT%lmaxd:, 0:, :)
    integer, allocatable, intent(out) ::  mlh_atom(:,:,:)
    integer, allocatable, intent(out) ::  nmem_atom(:, :)
    complex, allocatable, intent(out) ::  clnu_atom(:,:,:)

    ! - local scalars -
    integer                           ::  iop
    integer                           ::  itype
    integer                           ::  ieq
    integer                           ::  iatom
    integer                           ::  ilh
    integer                           ::  imem
    integer                           ::  itypsym
    integer                           ::  l
    integer                           ::  m
    integer                           ::  ok
    integer                           :: imqn_ms, imqn_mz ! delete it

    ! - local arrays -
    complex                           ::  carr1(-atomsT%lmaxd:atomsT%lmaxd)
    complex                           ::  carr2(-atomsT%lmaxd:atomsT%lmaxd)

    ! This block calculates maximum number of members in a lattice harmonic
    allocate(nmem_atom(0:latharT%nlhd, atomsT%nat))
    iatom = 0
    memd_atom = 0
    do itype = 1, atomsT%ntype
      do ieq = 1, atomsT%neq(itype)
        iatom   = iatom + 1
        itypsym = atomsT%ntypsy(iatom)
        iop = symT%invtab(atomsT%ngopr(iatom)) ! todo is this the correct operation??
!        iop = atomsT%ngopr(iatom) ! todo is this the correct operation??
        do ilh = 0, latharT%nlh(itypsym)
          l     = latharT%llh(ilh, itypsym)
          carr1 = 0
          carr2 = 0
          do imem = 1, latharT%nmem(ilh, itypsym)
            m        = latharT%mlh (imem, ilh, itypsym)
            carr1(m) = latharT%clnu(imem, ilh, itypsym)
          end do

          carr2(-l:l) = matmul(transpose(d_wgn(-l:l,-l:l,l,iop)), carr1(-l:l))
!          carr2(-l:l) = matmul(conjg(d_wgn(-l:l,-l:l,l,iop)), carr1(-l:l))
!         do m = -l,l
!            WRITE(200,'(2i4,4f15.10)') l,m,carr1(m),carr2(m)
!         end do
          imem = 0
          do m = -l, l
            if (abs(carr2(m)) > 1E-8 ) then
              imem = imem + 1
            end if
          end do
          memd_atom = max(memd_atom, imem)
        end do
      end do
    end do

    ! now the intent(out) quantities can be calculated
    allocate(clnu_atom(memd_atom,0:latharT%nlhd,atomsT%nat), mlh_atom(memd_atom,0:latharT%nlhd,atomsT%nat), stat = ok  )
    if( ok .ne. 0 ) then
      call juDFT_error('error during allocation of clnu_atom/mlh_atom', calledby='rotate_clnu', hint='storage full?', file='jpPotDens.F90', line=299)
    endif

    mlh_atom = 0
    clnu_atom = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atomsT%ntype
       do ieq = 1, atomsT%neq(itype)
        iatom = iatom + 1
        itypsym = atomsT%ntypsy(iatom)
        iop = symT%invtab(atomsT%ngopr(iatom)) ! todo is this the correct operation??
     !   write(*, *) iatom, iop
!        iop = atomsT%ngopr(iatom)
        do ilh = 0, latharT%nlh(itypsym)
          l = latharT%llh(ilh, itypsym)
          carr1 = 0
          carr2 = 0
          do imem = 1, latharT%nmem(ilh, itypsym)
             m = latharT%mlh (imem, ilh, itypsym)
             carr1(m) = latharT%clnu(imem, ilh, itypsym)
          end do
          carr2(-l:l) = matmul(transpose(d_wgn(-l:l, -l:l, l, iop)), carr1(-l:l))
     !     write(*, *) iatom, ilh
     !     do imem = 1, latharT%nmem(ilh, itypsym)
     !       write (*, *)
     !       write(*, *) carr1(latharT%mlh(imem, ilh, itypsym))
     !       write(*, *) latharT%clnu(imem, ilh, itypsym)
     !       write(*, *) carr2(latharT%mlh(imem, ilh, itypsym))
     !       write (*, *)
     !     end do
!          do imqn_mz = -l, l
!            do imqn_ms = -l, l
!              write(*, *) imqn_mz, imqn_ms, d_wgn(imqn_mz, imqn_ms, l, iop)
!            end do
!          end do
          ! todo think about transposing, and conjugating and using the right iop
          imem = 0
          if ( abs(carr2(0)) > 1e-8) then
            imem = imem + 1
            clnu_atom(imem, ilh, iatom) = carr2(0)
            mlh_atom(imem, ilh, iatom) = 0
          end if
          do m = 1, l
            if( (abs(carr2(m)) .gt. 1E-8) ) then
              imem                      = imem + 1
              clnu_atom(imem, ilh, iatom) = carr2(m)
              mlh_atom (imem, ilh, iatom) = m
            end if
            if( (abs(carr2(-m)) .gt. 1E-8) ) then
              imem                      = imem + 1
              clnu_atom(imem, ilh, iatom) = carr2(-m)
              mlh_atom (imem, ilh, iatom) = -m
            end if
          end do
          nmem_atom(ilh,iatom) = imem

          ! this loop is only for testing reason, carr should be completely 0 after it otherwise lattice harmonics have been rotated falsely
#ifdef DEBUG_MODE
            do imem = 1, nmem_atom(ilh, iatom)
               m = mlh_atom(imem, ilh, iatom)
               carr2(m) = 0
            end do
!            write (*, *)
!            do m = -l, l
!              write (*, *) m, carr2(m)
!            end do
!            write (*, *) carr2
!            write (*, *)
            if (any(abs(carr2) .gt. 1E-8)) then
              call juDFT_error('rotation error', calledby='rotate_clnu', hint='Check algorithm rotating the lattice harmonics.', file='jpPotDens.F90', line=307)
            endif
#endif
          end do
        end do
      end do
!      write (*, *) 'clnu'
!      iatom = 1
!      do itype = 1, atomsT%ntype
!        do ieq = 1, atomsT%neq(itype)
!          itypsym = atomsT%ntypsy(iatom)
!          do ilh = 0, latharT%nlh(itypsym)
!            do imem = 1, latharT%nmem(ilh, itypsym)
!              write (*, *) iatom, itypsym, ilh, imem
!              write(*, *) latharT%clnu(imem, ilh, itypsym), clnu_atom(imem, ilh, iatom)
!              write(*, *)
!            end do
!          end do
!          iatom = iatom + 1
!        end do
!      end do

!      write(*, *) 'mlh'
!      iatom = 1
!      do itype = 1, atomsT%ntype
!        do ieq = 1, atomsT%neq(itype)
!          itypsym = atomsT%ntypsy(iatom)
!          do ilh = 0, latharT%nlh(itypsym)
!            do imem = 1, latharT%nmem(ilh, itypsym)
!              write (*, *) iatom, itypsym, ilh, imem
!              write (*, *) latharT%mlh(imem, ilh, itypsym), mlh_atom(imem, ilh, iatom)
!              write (*, *)
!            end do
!          end do
!          iatom = iatom + 1
!        end do
!      end do
!      NOstopNO
   end subroutine rotate_clnu

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculating the plane-wave expansion coefficients of the functional derivative of any interstitial kernel with respect to the
  !> unperturbed density.
  !>
  !> @note
  !> This routine is tested as we compared Fleur and juPhon potentials
  !> @note
  !> At the moment, we use the x-alpha potential, but in principle in this routine the call to libxc can be integrated
  !> @note
  !> This routine is similiar to the routine that calculates the interstitial x-alpha kernel in Fleur.
  !>
  !> @param[in] stars : Stars type, see types.f90
  !> @param[in] gdp   : G-vectors of potentials and densities
  !> @param[in] ngdp  : Number of G-vectors for potentials and densities
  !> @param[in] qpw   : Star coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in] fg3G  : Plane wave interstitial coefficients of the functional derivative of the xc-kernel with respect to the
  !>                    density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcIRdVxcKern(stars, gdp, ngdp, qpw, fg3G)

    use m_fft3d
    use m_types, only : t_stars

    implicit none

    ! Type parameter
    type(t_stars),              intent(in) :: stars

    ! Scalar parameter
    integer,                    intent(in) :: ngdp

    ! Array parameter
    integer,                    intent(in) :: gdp(:, :)
    complex,                    intent(in) :: qpw(:)
    complex,       allocatable, intent(out):: fg3G(:)

    ! Scalar variables
    integer                                :: ifftd
    integer                                :: iG

    ! Array variables
    real,          allocatable             :: af3(:), bf3(:)
    real,          allocatable             :: VxcIRKern(:)
    complex,       allocatable             :: fg3(:)

    ! Init and allocate output variable
    allocate(fg3G(ngdp))
    fg3G(:) = cmplx(0., 0.)

    ! Size of FFT mesh
    ifftd = 27 * stars%k1d * stars%k2d * stars%k3d
    allocate( af3(0:ifftd - 1), bf3(0:ifftd - 1), VxcIRKern(ifftd), fg3(stars%n3d))
    af3(:) = 0.
    bf3(:) = 0.
    VxcIRKern(:) = 0.
    fg3(:) = cmplx(0., 0.)

    ! FFT of qpw from reciprocal to direct space
    ! We can use the fft which is desiged for star expanded quantities because here we have the qpw which is given in stars.
    call fft3d(af3, bf3, qpw, stars, +1)

    if ( any( abs(bf3) >= 1e-8 ) ) then
      write(*, *) 'Warning: FFT in calcIRdVxcKern has complex contributions!'
    end if
    bf3 = 0

    ! Calculate functional derivative of kernel based on qpw in direct space because we know the representation in direct space.
    call calcKernDerOnGrid(ifftd, 1., af3, VxcIRKern)

    ! Back-FFT to direct space. We only have made a functional derivative which does not break the symmetry so we can still use this
    ! kind of fft. We use the fft of the new fleur version, if we have no optional parameter the function is the same as in old
    ! fleur.
    CALL fft3d(VxcIRKern, bf3, fg3, stars, -1)

    ! Transform from star representation to plane wave representation
    do iG = 1, ngdp
      fg3G(iG) = fg3G(iG) + fg3(stars%ig(gdp(1, iG), gdp(2, iG), gdp(3, iG))) * stars%rgphs(gdp(1, iG), gdp(2, iG), gdp(3, iG))
    end do

  end subroutine calcIRdVxcKern

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Convolutes kernel with density for interstitial suitable for gradient of unperturbed potential and for
  !> first variation of the potential
  !>
  !> @note
  !> This routine is tested as we compared Fleur and juPhon potentials
  !>
  !> @param[in]  stars     : Stars type, see types.f90
  !> @param[in]  ngdp      : Number of G-vectors for potentials and densities
  !> @param[in]  ngpqdp    : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  iqpt      : Index of q-point in kpts file that is evaluated
  !> @param[in]  f1IR      : Plane-wave interstitial coefficients of first quantity relating to pdG2FouMv
  !> @param[in]  f2IR      : Plane-wave interstitial coefficients of second quantity relating to pdG2FouM
  !> @param[in]  pdG2FouM  : Mapping array that maps f1IR plane-wave interstitial coefficiens to unique point on FFT mesh
  !> @param[in]  pdG2FouMv : Mapping array that maps f2IR plane-wave interstitial coefficiens to unique point on FFT mesh
  !> @param[out] f3IR      : Plane-wave interstitial coefficients of convoluted quantity of f1IR and f2IR.
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine convolGrRhoKern(stars, ngdp, ngpqdp, f1IR, f2IR, pdG2FouM, pdG2FouMv, f3IR, iqpt)

    use m_cfft
    use m_types

    implicit none

    ! Type parameter
    ! -------------
    type(t_stars),        intent(in)  :: stars

    ! Scalar parameter
    integer,              intent(in)  :: ngdp
    integer,              intent(in)  :: ngpqdp
    integer,              intent(in)  :: iqpt

    ! Array parameter
    !----------------
    ! 2 quantities to convolute
    complex,              intent(in)  :: f1IR(:)
    complex,              intent(in)  :: f2IR(:)
    ! maps from a G to respective mesh entry of FFT mesh
    integer,              intent(in)  :: pdG2FouM(:)
    integer,              intent(in)  :: pdG2FouMv(:)
    ! Convoluted quantity
    complex,              intent(out) :: f3IR(:)

    ! Maximal length of FFT mesh
    integer                           :: ifftd
    ! Loop index
    integer                           :: iG
    ! Scaling factor for canceling FFT artefacts.
    real                              :: scaling
    integer                           :: ii

    ! Real and imaginary parts of quantities to be convoluted (1, 2) and convoluted quantity 3
    real,    allocatable              :: rf1(:)
    real,    allocatable              :: if1(:)
    real,    allocatable              :: rf2(:)
    real,    allocatable              :: if2(:)
    real,    allocatable              :: rf3(:)
    real,    allocatable              :: if3(:)


    ! Length of FFT mesh. The cartesian dimensions are stored sequentially. The G's expand maximally from -k_i to k_i (from this
    ! cube a ball of radius gmax is cut off) so 3 * k_i should be in principle large enough for all given G's components. And we
    ! have some free space to avoid aliasing. However, this factor of 3 is only working well for symmetrized quantities which is
    ! why we have to expand the k1d, k2d, k3d if we expand our quantities in plane waves.
    ifftd = 27 * stars%k1d * stars%k2d * stars%k3d

    ! All quantities have the size of the FFT mesh.
    allocate( rf1(0:ifftd - 1), if1(0:ifftd - 1), rf2(0:ifftd - 1), if2(0:ifftd - 1), rf3(0:ifftd - 1), if3(0:ifftd - 1) )

    rf1(:) = 0
    rf2(:) = 0
    rf3(:) = 0
    if1(:) = 0
    if2(:) = 0
    if3(:) = 0

    ! Extract the real and imaginary part of the expansion coefficients of the quantities to be convoluted for every G-vector and
    ! map them with pdG2FouM to its respective mesh point.
    do iG = 1, ngpqdp
      rf1(pdG2FouMv(iG)) = real(f1IR(iG))
      if1(pdG2FouMv(iG)) = aimag(f1IR(iG))
    end do

    do iG = 1, ngdp ! is it kimax - 1 or kimax?
      rf2(pdG2FouM(iG)) = real(f2IR(iG))
      if2(pdG2FouM(iG)) = aimag(f2IR(iG))
    end do

      ! Complex FFT of 1st quantity into real space, this is done as it is done in FLEUR fft3d
      call Cfft(rf1, if1, ifftd, 3 * stars%k1d, 3 * stars%k1d, 1)
      call Cfft(rf1, if1, ifftd, 3 * stars%k2d, 9 * stars%k1d * stars%k2d, 1)
      call Cfft(rf1, if1, ifftd, 3 * stars%k3d, ifftd, 1)

      ! Complex FFT of 2nd quantity into real space
      call Cfft(rf2, if2, ifftd, 3 * stars%k1d, 3 * stars%k1d, 1)
      call Cfft(rf2, if2, ifftd, 3 * stars%k2d, 9 * stars%k1d * stars%k2d, 1)
      call Cfft(rf2, if2, ifftd, 3 * stars%k3d, ifftd, 1)

      ! Exploiting the convolution theorem
      do ii = 0, ifftd - 1
        rf3(ii) = real( cmplx(rf1(ii),if1(ii)) * cmplx(rf2(ii), if2(ii)))
        if3(ii) = aimag( cmplx(rf1(ii),if1(ii)) * cmplx(rf2(ii), if2(ii)))
      end do

# ifdef DEBUG_MODE
      ! It is sufficient to test it only for iqpt == 1 because the gradient of rho is calculated outside the q-loop and the density
      ! variation has only to be real for q = 0
      if(iqpt == 1) then
        if ( any(abs(if3) > 1e-7 )) then
          write(*, *) 'Convolution FFT has imaginary components'
          !NOstopNO'Convolution FFT has imaginary components'
        end if
        if3 = 0
      end if
#endif

      ! Complex FFT of convoluted quantity into reciprocal space
      call Cfft(rf3, if3, ifftd, 3 * stars%k1d, 3 * stars%k1d, -1)
      call Cfft(rf3, if3, ifftd, 3 * stars%k2d, 9 * stars%k1d * stars%k2d, -1)
      call Cfft(rf3, if3, ifftd, 3 * stars%k3d, ifftd, -1)

      ! We have to care for the artefacts of this FFT
      scaling = 1. / ifftd
      f3IR = cmplx( 0.0, 0.0 )
      ! Map convoluted quantity from FFT mesh to plane-wave expansion coefficient representation.
      do iG = 1, ngpqdp !kimax is max G-vector
        f3IR(iG) = f3IR(iG) + cmplx(rf3(pdG2FouMv(iG)), if3(pdG2FouMv(iG))) * scaling
      end do

  end subroutine convolGrRhoKern

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculating the functional derivative of any MT kernel with respect to the unperturbed density
  !>
  !> @param[in] atoms       : Atoms type, see types.f90
  !> @param[in] dimens      : Dimension type, see types.f90
  !> @param[in] sphhar      : Lattice harmonics type, see types.f90.
  !> @param[in] rho0MT      : Lattice harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[in] nmem_atom   : Number of lattice harmonic members for every atom.
  !> @param[in] clnu_atom   : Phase mediating between stars and plane waves.
  !> @param[in] mlh_atom    : Magnetic quantum number m of lattice harmonic members for every atom.
  !> @param[in] gWghts      : Gauss weights for Gauss quadrature
  !> @param[in] ylm         : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
  !> @param[in] dKernMTGPts : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                         respect to the density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcMTdVxcKern(atoms, dimens, sphhar, rho0MT, nmem_atom, clnu_atom, mlh_atom, gWghts, ylm, dKernMTGPts)

    use m_gaussp
    use m_types, only : t_atoms, t_dimension, t_sphhar
    use m_ylm_old

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_dimension),              intent(in)  :: dimens
    type(t_sphhar),                 intent(in)  :: sphhar

    ! Array parameters
    real,                           intent(in)  :: rho0MT(:, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:, 0:, :)
    integer,                        intent(in)  :: mlh_atom(:, 0:, :)
    real,              allocatable, intent(out) :: gWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable, intent(out) :: ylm(:, :)
    real,              allocatable, intent(out) :: dKernMTGPts(:, :, :)

    ! Scalar local variables
    integer                                     :: igmesh
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ieqat
    integer                                     :: ptSym
    integer                                     :: irmesh
    integer                                     :: ilh
    integer                                     :: oqn_l
    integer                                     :: imem
    integer                                     :: mqn_m
    integer                                     :: lm

    ! Local Array Variables
    real,              allocatable              :: gPts(:, :) ! gaussian points to exactly integrate spherial harmonics
    complex,           allocatable              :: rhoMTGpts(:, :)

    allocate( gWghts(dimens%nspd), ylm(dimens%nspd, ( atoms%lmaxd + 1 )**2), &
      & dKernMTGPts(dimens%nspd, atoms%jmtd, atoms%nat), gPts(3, dimens%nspd), rhoMTGpts(dimens%nspd, atoms%jmtd) )

    gWghts(:) = 0.
    ylm(:, :) = cmplx(0., 0.)
    dKernMTGPts(:, :, :) = 0.
    gPts(:, :) = 0.
    rhoMTGpts(:, :) = cmplx(0., 0.)

    ! generates dimension%nspd points on a spherical shell with radious 1.0. The angular mesh is equidistant in phi, theta are the
    ! zeros of the legendre polynomials.?
    ! Gaussian points to exactly integrate product of spherical harmonics up to lmax
    ! todo understand the gaussp routine ! Gauss quadrature on wikipedia!
    call gaussp( atoms%lmaxd, gPts, gWghts )

    ! The ylms are note filled linearily as a compromise to perform a linear run through all the other arrays
    do igmesh = 1, dimens%nspd
      call Ylm4( atoms%lmaxd, gPts(:, igmesh), ylm(igmesh, :) )
    end do

    dKernMTGPts(:, :, :) = 0
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ! Evaluate the unperturbed density on the Gauss mesh then calculate the functional derivative of the density
        ptSym = atoms%ntypsy(iatom)
        rhoMTGpts(:, :) = cmplx(0., 0.)
        do ilh = 0, sphhar%nlh(ptsym)
          oqn_l = sphhar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = oqn_l * ( oqn_l + 1 ) + mqn_m + 1
            do irmesh = 1, atoms%jri(itype)
              do igmesh = 1, dimens%nspd
                rhoMTGpts(igmesh, irmesh) = rhoMTGpts(igmesh, irmesh) + rho0MT(irmesh, ilh, itype) * clnu_atom(imem, ilh, iatom) &
                                          & * ylm(igmesh, lm)
              end do ! igmesh
            end do ! irmesh
          end do ! imem
        end do ! ilh
        if ( any( abs(aimag( rhoMTGpts )) > 1e-7 ) ) then
          write(*, *) 'Warning rhoMTGpts has imaginary components.'
        end if
        ! Calculate functional derivative of the density on the Gauss mesh for every atom
        do irmesh = 1, atoms%jri(itype)
          call calcKernDerOnGrid(dimens%nspd, 1., real(rhoMTGpts(:, irmesh)), dKernMTGPts(:, irmesh, iatom))
        end do ! irmesh
      end do ! ieqat
    end do ! itype

  end subroutine calcMTdVxcKern

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Convolute gradient of unperturbed density with functional derivative of xc - kernel
  !>
  !> @param[in]  atoms       : Atoms type, see types.f90
  !> @param[in]  dimens      : Dimension type, see types.f90
  !> @param[in]  grRho0MT    : Spherical harmonic coefficients of the gradient of the muffin-tin unperturbed density
  !> @param[in]  gWghts      : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  dKernMTGpts : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                           respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  ylm         : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
  !>                           created in m_jppotdenshelper::calcmtdvxckern
  !> @param[out] grVxc0MT    : Spherical harmonic muffin-tin coefficients of the gradient of the unperturbed xc potential
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine convolMTgrVeff0dKern(atoms, dimens, grRho0MT, dKernMTGPts, gWghts, ylm, grVxc0MT)

    use m_types, only : t_atoms, t_dimension

    implicit none

    ! Type parameter
    type(t_atoms),     intent(in)  :: atoms
    type(t_dimension), intent(in)  :: dimens

    ! Array parameter
    complex,           intent(in)  :: grRho0MT(:, :, :, :)
    complex,           intent(in)  :: ylm(:, :)
    real,              intent(in)  :: gWghts(:) ! gaussian weights belonging to gausPts
    real,              intent(in)  :: dKernMTGPts(:, :, :)
    complex,           intent(out) :: grVxc0MT(:, :, :, :)

    ! Local scalar variables
    integer                        :: idir
    integer                        :: iatom
    integer                        :: itype
    integer                        :: ieqat
    integer                        :: igmesh ! Loop variable over sampling points of spherical Gauss mesh
    integer                        :: irmesh ! Loop variable over radial MT mesh
    integer                        :: oqn_l
    integer                        :: lm_lonly !reduce multiplication when calculating lm
    integer                        :: mqn_m
    integer                        :: lm
    complex                        :: grVxcMTKernAdd

    ! Local allocatable variables
    complex, allocatable           :: grRhoMTGpts(:, :) !grRhoMT on Gauss mesh
    real, allocatable              :: grVxcMTKernGPts(:, :)

    grVxc0MT(:, :, :, :) = cmplx(0., 0.)

    allocate( grRhoMTGpts(dimens%nspd, atoms%jmtd), grVxcMTKernGPts(dimens%nspd, atoms%jmtd) )
    grRhoMTGpts(:, :) = cmplx(0., 0.)
    grVxcMTKernGPts(:, :) = 0.

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          grRhoMTGpts(:, :) = cmplx(0., 0.)
          grVxcMTKernGPts(:, :) = 0.
          ! Density's gradient has l + 1 entries
          do oqn_l = 0, atoms%lmax(itype)
            lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lOnly + mqn_m
              ! Evaluate grRho on spherical Gauss mesh in order to apply Gauss quadrature.
              do irmesh = 1, atoms%jri(itype)
                do igmesh = 1, dimens%nspd
                  grRhoMTGpts(igmesh, irmesh) = grRhoMTGpts(igmesh, irmesh) + grRho0MT(irmesh, lm, iatom, idir) * ylm(igmesh, lm)
                end do ! igmesh
              end do ! irmesh
            end do ! mqn_m
          end do ! oqn_l
#if DEBUG_MODE
          if ( any( aimag( grRhoMTGpts ) > 1e-7 ) ) then
            write(*, *) 'Warning rhoMTGpts has imaginary components.'
          end if
#endif
          ! On the spherical Gauss mesh the integral reduces to a weighted (gWghts) sum (over all sampling points on the Gauss mesh)
          ! of the MT exchange-correlation kernel and either the density's gradient or the first variation of the gradient.
          do irmesh = 1, atoms%jri(itype)
            do igmesh = 1, dimens%nspd
              grVxcMTKernGPts(igmesh, irmesh) = grVxcMTKernGpts(igmesh, irmesh) + real(grRhoMTGpts(igmesh, irmesh)) &
                                              & * dKernMTGPts(igmesh, irmesh, iatom) * gWghts(igmesh)
            end do
          end do
          do oqn_l = 0, atoms%lmax(itype)
            lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lOnly + mqn_m
              do irmesh = 1, atoms%jri(itype)
              ! Back-transformation of the MT coefficients. Now they are expansion coefficients of the MT grid.
                grVxcMTKernAdd = dot_product( grVxcMTKernGPts(:dimens%nspd, irmesh), conjg(ylm(: dimens%nspd, lm)) )
              ! Add this contribution to MT exchange-correlation contribution to the potential
                grVxc0MT(irmesh, lm, iatom, idir) = grVxc0MT(irmesh, lm, iatom, idir) + grVxcMTKernAdd
              end do ! irmesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! itype

  end subroutine convolMTgrVeff0dKern

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Convolute the first variation of density with the functional derivative of the xc - kernel with respect to the unperturbed
  !> density.
  !>
  !> @param[in]  atoms       : Atoms type, see types.f90
  !> @param[in]  dimens      : Dimension type, see types.f90
  !> @param[in]  iqpt        : Index of q-point in kpts file that is evaluated
  !> @param[in]  rho1MT      : Spherical harmonic muffin-tin expansion coefficients of the linear density variation
  !> @param[in]  gWghts      : Gauss weights for Gauss quadrature created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  dKernMTGpts : Spherical harmonic muffin-tin coefficients of the functional derivative of the xc-kernel with
  !>                           respect to the density created in m_jppotdenshelper::calcmtdvxckern
  !> @param[in]  ylm         : Set of spherical harmonics whose arguments are the unit vectors of the Gauss mesh points up to lmax
  !>                           created in m_jppotdenshelper::calcmtdvxckern
  !> @param[out] grVeff0MT   : Spherical harmonic muffin-tin coefficients of the first-order variation of the xc potential
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine convolMTVeff1dKern(atoms, dimens, iqpt, rho1MT, gWghts, dKernMTGPts, ylm, vxc1MT)
    use m_types

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms
    type(t_dimension),             intent(in)  :: dimens

    ! Scalar parameter
    integer,                       intent(in)   :: iqpt
    ! Array parameter
    complex,                       intent(in)  :: rho1MT(:, :, :, :)
    complex,                       intent(in)  :: ylm(:, :)
    real,                          intent(in)  :: gWghts(:) ! gaussian weights belonging to gausPts
    real,                          intent(in)  :: dKernMTGPts(:, :, :)
    complex,                       intent(out) :: vxc1MT(:, :, :, :)

    ! Local scalar variables
    integer                                    :: idir
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: igmesh ! Loop variable over sampling points of spherical Gauss mesh
    integer                                    :: irmesh ! Loop variable over radial MT mesh
    integer                                    :: oqn_l
    integer                                    :: lm_lonly !reduce multiplication when calculating lm
    integer                                    :: mqn_m
    integer                                    :: lm
    complex                                    :: vxcMT1KernAdd

    ! Local allocatable variables
    complex,           allocatable             :: rhoMT1Gpts(:, :) !grRhoMT on Gauss mesh
    complex,           allocatable             :: vxcMT1KernGPts(:, :)

    vxc1MT(:, :, :, :) = cmplx(0., 0.)

    allocate( rhoMT1Gpts(dimens%nspd, atoms%jmtd), vxcMT1KernGPts(dimens%nspd, atoms%jmtd) )
    rhoMT1Gpts(:, :) = cmplx(0., 0.)
    vxcMT1KernGPts(:, :) = cmplx(0., 0.)

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          rhoMT1Gpts(:, :) = 0
          vxcMT1KernGPts(:, :) = 0
          do oqn_l = 0, atoms%lmax(itype)
            lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lOnly + mqn_m
              ! Evaluate grRho on spherical Gauss mesh in order to apply Gauss quadrature.
              do irmesh = 1, atoms%jri(itype)
                do igmesh = 1, dimens%nspd
                  rhoMT1Gpts(igmesh, irmesh) = rhoMT1Gpts(igmesh, irmesh) + rho1MT(irmesh, lm, iatom, idir) * ylm(igmesh, lm)
                end do ! igmesh
              end do ! irmesh
            end do ! mqn_m
          end do ! oqn_l
#if DEBUG_MODE
          if (iqpt == 1) then
            if ( any( aimag( rhoMT1Gpts ) > 1e-7 ) ) then
              write(*, *) 'Warning rhoMTGpts has imaginary components.'
            end if
          end if
#endif
          ! On the spherical Gauss mesh the integral reduces to a weighted (gWghts) sum (over all sampling points on the Gauss mesh)
          ! of the MT exchange-correlation kernel and either the density's gradient or the first variation of the gradient.
          do irmesh = 1, atoms%jri(itype)
            do igmesh = 1, dimens%nspd
              vxcMT1KernGPts(igmesh, irmesh) = vxcMT1KernGPts(igmesh, irmesh) + rhoMT1Gpts(igmesh, irmesh) &
                                              & * dKernMTGPts(igmesh, irmesh, iatom) * gWghts(igmesh)
            end do
          end do
          do oqn_l = 0, atoms%lmax(itype)
            lm_lonly = oqn_l * ( oqn_l + 1 ) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_lOnly + mqn_m
              do irmesh = 1, atoms%jri(itype)
              ! Back-transformation of the MT coefficients. Now they are expansion coefficients of the MT grid.
                vxcMT1KernAdd = dot_product( ylm(: dimens%nspd, lm), vxcMT1KernGPts(:dimens%nspd, irmesh) )
              ! Add this contribution to MT exchange-correlation contribution to the potential
                vxc1MT(irmesh, lm, iatom, idir) = vxc1MT(irmesh, lm, iatom, idir) + vxcMT1KernAdd
              end do ! irmesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! itype

  end subroutine convolMTVeff1dKern

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the x-alpha functional derivative of the xc-kernel with respect to the unperturbed density in real space.
  !>
  !> @param[in] nGridPts        : Number of grid points of the real space mesh of rhoMTGpts
  !> @param[in] alpha           : Parameter alpha in the x-alpha kernel.
  !> @param[in] rhoMTGpts       : Unperturbed density on a real-space mesh
  !> @param[in] grVxcMTKernGPts : Functional derivative of the x-alpha kernel with respect to the density.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcKernDerOnGrid(nGridPts, alpha, rhoMTGpts, grVxcMTKernGPts)

    use m_jPConstants, only : pi
    use m_types

    implicit none

    integer,           intent(in)  :: nGridPts
    real,              intent(in)  :: alpha
    real,              intent(in)  :: rhoMTGpts(:)
    real,              intent(out) :: grVxcMTKernGPts(:)

    real                           :: prfac
    integer                        :: imesh
    real                           :: rhoCapped

    grVxcMTKernGPts(:) = 0.
    prfac = 1. / 9. / pi
    do imesh = 1, nGridPts
      rhoCapped = max( 1e-15, rhoMTGpts(imesh) )
      ! We merged all factors from FLEUR together converted it in hartree units (division by 2) and multiplied it with the 1 / 3
      ! from the derivative which is the prefactor used here so we have 1 / 3 * (3 / pi)**(1/3) * rho**(-2 / 3)
      grVxcMTKernGPts(imesh) = - prfac**(1. / 3.) * alpha * rhoCapped**(-2. / 3.)
    end do

  end subroutine calcKernDerOnGrid

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Warps interstitial potential, i.e. convolutes intersitial potential with stepfunction.
  !>
  !> @note
  !> This routine is recycled and adapted from Fleur warping
  !>
  !> @param[in]  stars      : Stars type, see types.f90.
  !> @param[in]  ngpqdp     : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  idir       : Index of displacement direction
  !> @param[in]  gpqdp      : Number of G-vectors for shifted G-set for q-point with index iqpt
  !> @param[in]  pot        : Interstitial plane-wave expanded quantitiy to be warped
  !> @param[out] pot_warped : Warped interstitial plane-wave expanded quantitiy
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine WarpIRPot( stars, ngpqdp, idir, gpqdp, pot, pot_warped )

    use m_types
    use m_cfft

    implicit none

    ! Type parameters
    type(t_stars),                                intent(in)  :: stars

    ! Scalar parameters
    integer,                                      intent(in)  :: ngpqdp
    integer,                                      intent(in)  :: idir

    ! Array parameters
    integer,                                      intent(in)  :: gpqdp(:, :)
    complex,                                      intent(in)  :: pot(:, :)
    complex,                                      intent(out) :: pot_warped(:)

    ! Array variables
    integer,                         allocatable              :: igfft(:)
    real, allocatable :: mygfft(:, :)
    real,                            allocatable              :: gfft(:, :)
    integer                                                  :: nfft(3)
    integer                                                   :: gabs(3)


    ! Scalar variables
    integer                                                   :: ifftd
    integer                                                   :: iG
    real                                                      :: scaling
    integer                                                   :: imesh
    integer :: idir1

    REAL                :: pgfftF(0:(2*stars%k1d+1)*(2*stars%k2d+1)*(2*stars%k3d+1)-1)
    INTEGER             :: igfftF(0:(2*stars%k1d+1)*(2*stars%k2d+1)*(2*stars%k3d+1)-1,2)

    nfft = [3 * stars%k1d, 3 * stars%k2d, 3 * stars%k3d]
    ifftd = product(nfft)

    allocate( igfft(ngpqdp), gfft(0: ifftd - 1, 2) )
    allocate( mygfft(0: ifftd - 1, 2) )
    gabs = 0
    igfft = 0
    do iG = 1, ngpqdp
      do idir1 = 1, 3
        if ( gpqdp(idir1, iG) >= 0 ) then
          gabs(idir1) = gpqdp(idir1, iG)
        else
          gabs(idir1) = gpqdp(idir1, iG) + nfft(idir1)
        end if
      end do
      igfft(iG) = gabs(1) + gabs(2) * nfft(1) + gabs(3) * nfft(1) * nfft(2)
    end do

!    open(1000, file='igfftjp', form='formatted')
!    do iG = 1, ngdp
!    !  write (1000, '(i10,2x,i3,2x,i3,2x,i3,2x,i8)') ( gdp(1, iG) + stars%k1d ) + 100 * ( gdp(2, iG) + stars%k2d ) + 10000 * ( gdp(3, iG) + stars%k3d ), &
!    !    & gdp(1, iG), gdp(2, iG), gdp(3, iG), igfft(iG)
!      write (1000, '(i10,2x,i3,2x,i3,2x,i3,2x,i8)') iG, &
!        & gdp(1, iG), gdp(2, iG), gdp(3, iG), igfft(iG)
!    end do
!    close(1000)

    gfft = 0
    pot_warped=0

      do iG = 1, ngpqdp
        gfft(igfft(iG), 1) = real(pot(iG, idir))
        gfft(igfft(iG), 2) = aimag(pot(iG, idir))
      end do
   !   mygfft = 0
   !   mygfft = gfft

!      open(1000, file='foo', form='formatted')
!      do imesh = 0, ubound( gfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)


      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(1), nfft(1), 1)
      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(2), nfft(1) * nfft(2), 1)
      call cfft(gfft(0, 1), gfft(0, 2), ifftd, nfft(3), ifftd, 1)

!      open(1000, file='foo1', form='formatted')
!      do imesh = 0, ubound( gfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)

   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(1), nfft(1), 1)
   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(2), nfft(1) * nfft(2), 1)
   !   call cfft(mygfft(:, 1), mygfft(:, 2), ifftd, nfft(3), ifftd, 1)


!      open(1000, file='foo2', form='formatted')
!      do imesh = 0, ubound( mygfft, dim=1)
!        write (1000, '(i8,f20.13)') imesh, mygfft(imesh, 1)
!        write (1000, '(i8,f20.13)') imesh, mygfft(imesh, 2)
!      end do
!      close(1000)

!    write (29999, *) stars%ufft
      !gfft = 0
      do imesh = 0, ifftd-1
        gfft(imesh, :) = gfft(imesh, :) * stars%ufft(imesh) ! todo is ufft correctlz initialized
       ! gfft(imesh, :) = stars%ufft(imesh) ! todo is ufft correctlz initialized
      end do

      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(1), nfft(1), -1)
      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(2), nfft(1) * nfft(2), -1)
      call cfft(gfft(:, 1), gfft(:, 2), ifftd, nfft(3), ifftd, -1)

!      open(1000, file='foo5', form='formatted')
!      do imesh = 0, ubound( mygfft, dim=1)
!        write (1000, '(i8,f25.15)') imesh, gfft(imesh, 1)
!        write (1000, '(i8,f25.15)') imesh, gfft(imesh, 2)
!      end do
!      close(1000)

!      igfftF = 0
!      pgfftF = 0
!      open(1000, file='mapArray', form='unformatted')
!      read(1000) igfftF, pgfftF
!      close(1000)

      !write (*, *) 'kimax und ngdp', stars%kimax, ngdp

      scaling = 1. / real(ifftd)
      do iG=1, ngpqdp
        pot_warped(iG) = pot_warped(iG) +  cmplx( gfft(igfft(iG), 1), gfft(igfft(iG), 2) ) * scaling
      end do
  end subroutine warpIRPot

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Auxiliary method which is similiar to FLEUR's phasy routine but does not implement symmetry.
  !>
  !> @details
  !> This method calculates for a given G-vector and q-vector: 4 π i^l exp((Gvec + qpoint) * taual) Y*_lm(Gvec + qpoint)
  !> Opposite to the original phasy1 routine, here, stars and any type of symmetry is not considered here.
  !>
  !> @param[in]  atoms  : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in]  cell   : Contains unit-cell related quantities; definition of its members in types.F90 file.
  !> @param[in]  Gvec   : G-vector to be evaluated
  !> @param[in]  qpoint : q-vector to be evaluated
  !> @param[out] pylm   : Result described in details, atom resolved due to lack of symmetry in first variation of the potentials.
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine phasy1nSym(atoms, cell, Gvec, qptn, pylm)

    use m_jPConstants, only : tpi, fpi, iu
    use m_ylm_old
    use m_types
    use m_cotra

    implicit none

    ! Scalar Type Arguments
    type(t_atoms),  intent(in)  ::  atoms
    type(t_cell),   intent(in)  ::  cell

    ! Array Arguments
    integer,        intent(in)  ::  Gvec(:)
    real,           intent(in)  ::  qptn(:)
    complex,        intent(out) ::  pylm(:, :)

    !-------------------------------------------------------------------------------------------------------------------------------
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
    !-------------------------------------------------------------------------------------------------------------------------------
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

    !-------------------------------------------------------------------------------------------------------------------------------
    ! Local Array Variables
    ! fpiul: stores 4 π i^l
    ! Gqext: stores G + q in external coordinates
    ! ylm  : stores Y_lm
    !-------------------------------------------------------------------------------------------------------------------------------
    complex, allocatable        ::  fpiul(:)
    complex, allocatable        ::  ylm(:)
    real                        ::  Gqext(3)

    allocate(fpiul(0:atoms%lmaxd))
    fpiul(:) = cmplx(0., 0.)
    ! calculates 4 π i^l resolved for every l, not divided by nop because no loop over symmetry operations
    do oqn_l = 0, atoms%lmaxd
       fpiul(oqn_l) = fpi * iu**oqn_l
    enddo


    ! calculates Y*_lm(\vec{G} + \vec{q}) for every l and m. The argument Gqext must be in external coordinates.
    allocate(ylm((atoms%lmaxd + 1)**2))
    ylm(:) = cmplx(0., 0.)
    Gqext(1:3) = matmul(cell%bmat(1:3, 1:3), real(Gvec(1:3) + qptn(1:3)))
    call ylm4(atoms%lmaxd, Gqext(1:3), ylm)
    ylm = conjg(ylm)


    ! calculates first exp(i (G + q) tau)  and multiplies recent factors before storing the final result to pylm
    iatom = 1
    pylm(:, :) = cmplx(0.,0.)
    do itype = 1, atoms%ntype
       do ineq = 1, atoms%neq(itype)
          x = tpi * dot_product(real(Gvec(1:3) + qptn(1:3)), atoms%taual(1:3, iatom))
          sf = exp(iu *  x)
          do oqn_l = 0, atoms%lmax(itype)
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

  end subroutine phasy1nSym

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates expansion coefficients of analytical gradient of function expanded in plane-waves
  !> @param[in]  stars       : Stars type, see types.f90
  !> @param[in]  ngdp        : Number of G-vectors for potentials and densities
  !> @param[in]  stExpandQ   : Expansion coefficients of quantity expanded in stars
  !> @param[in]  gdp         : G-vectors of potentials and densities
  !> @param[out] gVecExpandQ : Expansion coefficients of quantity expanded in plane waves
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcAnalytGradIRMat(cell, ngpqdp, gpqdp, qpoint, grExpCoeff)

    use m_types, only : t_cell
    use m_jPConstants, only : iu

    implicit none

    ! Type parameters
    type(t_cell), intent(in)  :: cell

    ! Scalar parameters
    integer                   :: ngpqdp

    ! Array parameters
    integer,      intent(in)  :: gpqdp(:, :)
    real,         intent(in)  :: qpoint(:)
    complex,      intent(out) :: grExpCoeff(:, :)

    ! Scalar variables
    integer                   :: iG
    integer                   :: idir

    ! Array variables
    real                      :: Gpqext(3)
    complex, allocatable      :: expCoeff(:)

    allocate(expCoeff(ngpqdp))
    expCoeff(:) = cmplx(0., 0.)

    do iG = 1, ngpqdp
      Gpqext(1:3) = matmul(cell%bmat(1:3, 1:3), gpqdp(1:3, iG) + qpoint(1:3))
      do idir = 1, 3
        grExpCoeff(iG, idir) = iu * Gpqext(idir) * expCoeff(iG)
      end do ! idir
    end do ! iG

  end subroutine CalcAnalytGradIRMat

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculate the interstitial part of the multipole moments for the gradient or the linear variation of the Hartree potential
  !> that result from a discontinuity of the unperturbed density.
  !>
  !> @param[in]  atoms         : Atoms type, see types.f90
  !> @param[in]  cell          : Unit cell type, see types.f90
  !> @param[in]  iDtype        : Index of atom type to which displaced atom belongs
  !> @param[in]  iDatom        : Index of atom that is displaced
  !> @param[in]  ngdp          : Number of G-vectors for potentials and densities
  !> @param[in]  gdp           : G-vectors of potentials and densities
  !> @param[in]  rho0IRpw      : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[out] qlmHartSurfIR : Interstitial part of the multipole moments for the gradient or the linear variation of the Hartree
  !>                             potential that result from a discontinuity of the unperturbed density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcQlmHarSurfIR( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw, qlmHartSurfIR )

    use m_types, only : t_atoms, t_cell
    use m_jPConstants, only : fpi, c_im, iu, tpi
    use m_ylm_old, only : ylm4
    use m_sphbes, only : sphbes
    use m_gaunt, only : Gaunt1

    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms
    type(t_cell),               intent(in)  :: cell

    ! Scalar parameter
    integer,                    intent(in)  :: ngdp
    integer,                    intent(in)  :: iDtype
    integer,                    intent(in)  :: iDatom

    ! Array parameter
    integer,                    intent(in)  :: gdp(:, :)
    complex,                    intent(in)  :: rho0IRpw(:)
    complex,                    intent(out) :: qlmHartSurfIR(:, :)

    ! Scalar variables
    integer                                 :: idir
    integer                                 :: iG
    integer                                 :: oqn_l
    integer                                 :: mqn_m
    integer                                 :: oqn_l1p
    integer                                 :: mqn_m1p
    integer                                 :: mqn_m2p
    integer                                 :: lm
    integer                                 :: lm1p
    complex                                 :: phaseFac
    complex                                 :: pref
    complex                                 :: temp1
    complex                                 :: temp2
    complex                                 :: temp3
    real                                    :: gauntFactor

    ! Array variables
    complex,       allocatable              :: ylm(:)
    real,          allocatable              :: sbes(:)
    real                                    :: gExt(3)

    ! Init assumed-shape dummy array
    qlmHartSurfIR(:, :) = cmplx(0., 0.)

    allocate( ylm( (atoms%lmax(iDtype) + 1)**2 ) )
    allocate( sbes(0:atoms%lmax(iDtype)) )
    ylm(:)   = cmplx(0., 0.)
    sbes(:) = 0.

    pref = fpi * atoms%rmt(iDtype) * atoms%rmt(iDtype)
    do iG = 1, ngdp
      gExt(1:3) = matmul( cell%bmat(1:3, 1:3), real(gdp(1:3, iG)) )

      ylm(:) = cmplx(0., 0.)
      call ylm4( atoms%lmax(iDtype), gExt(1:3), ylm )

      sbes(:) = 0
      call sphbes(atoms%lmax(iDtype), norm2(gExt(1:3)) * atoms%rmt(iDtype), sbes)

      phaseFac = exp( iu * tpi * dot_product(gdp(1:3, iG), atoms%taual(1:3, iDatom)))

      do oqn_l = 0, atoms%lmax(iDtype)
        temp1 = pref * phaseFac * atoms%rmt(iDtype)**oqn_l * rho0IRpw(iG)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do oqn_l1p = 0, atoms%lmax(iDtype)
            temp2 = temp1 * sbes(oqn_l1p) * iu**oqn_l1p
            do mqn_m1p = -oqn_l1p, oqn_l1p
              lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
              temp3 = temp2 * conjg(ylm(lm1p))
              do mqn_m2p = -1, 1
                gauntFactor = Gaunt1( oqn_l, oqn_l1p, 1, mqn_m, mqn_m1p, mqn_m2p, atoms%lmax(iDtype))
                do idir = 1, 3
                  qlmHartSurfIR(idir, lm) = qlmHartSurfIR(idir, lm) + c_im(idir, mqn_m2p + 2) * gauntFactor * temp3
                end do ! idir
              end do ! mqn_mpp
            end do ! mqn_mp
          end do ! oqn_lp
        end do ! mqn_m
      end do ! oqn_l
    end do ! iG

  end subroutine CalcQlmHarSurfIR

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculate muffin-tin part of the multipole moments for the gradient or the linear variation of the Hartree potential that
  !> result from a discontinuity of the unperturbed density.
  !>
  !> @param[in]  atoms         : Atoms type, see types.f90
  !> @param[in]  iDtype        : Index of atom type to which displaced atom belongs
  !> @param[in]  iDatom        : Index of atom that is displaced
  !> @param[in]  rho0MTsh      : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[out] qlmHartSurfMT : Muffin-tin part of the multipole moments for the gradient or the linear variation of the Hartree
  !>                             potential that result from a discontinuity of the unperturbed density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcQlmHarSurMT(atoms, iDtype, iDatom, rho0MTsh, qlmHartSurfMT)

    use m_types, only : t_atoms
    use m_gaunt, only : Gaunt1
    use m_jpConstants, only : c_im

    implicit none

    ! Type parameters
    type(t_atoms), intent(in)  :: atoms

    ! Scalar parameter
    integer,       intent(in)  :: iDtype
    integer,       intent(in)  :: iDatom

    ! Array parameters
    complex,       intent(in)  :: rho0MTsh(:, :, :, :)
    complex,       intent(out) :: qlmHartSurfMT(:, :)

    ! Scalar variables
    integer                    :: oqn_l
    integer                    :: oqn_l1p
    integer                    :: mqn_m1p
    integer                    :: mqn_m2p
    integer                    :: mqn_m
    integer                    :: lm
    integer                    :: lm1p
    integer                    :: idir
    real                       :: gauntFactor

    qlmHartSurfMT(:, :) = cmplx(0., 0.)

    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
        do oqn_l1p = 0, atoms%lmax(iDtype)
          do mqn_m1p = -oqn_l1p, oqn_l1p
            lm1p = oqn_l1p * (oqn_l1p + 1) + 1 + mqn_m1p
            do mqn_m2p = -1, 1
              gauntFactor = Gaunt1( oqn_l, oqn_l1p, 1, mqn_m, mqn_m1p, mqn_m2p, atoms%lmax(iDtype))
              do idir = 1, 3
                qlmHartSurfMT(idir, lm) = qlmHartSurfMT(idir, lm) + c_im(idir, mqn_m2p + 2) * atoms%rmt(iDtype)**(oqn_l + 2)     &
                  & * rho0MTsh(atoms%jri(iDtype), lm1p, iDatom, 1) * gauntFactor
              end do ! idir
            end do ! mqn_m2p
          end do ! mqn_m1p
        end do ! oqn_l1p
      end do ! mqn_m
    end do ! oqn_l

  end subroutine CalcQlmHarSurMT

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculate multipole moments for the gradient or the linear variation of the Hartree potential that result from a discontinuity
  !> of the unperturbed density.
  !>
  !> @details
  !> Main routine that calls routines for calculating the interstitial and the muffin-tin contribution
  !>
  !> @param[in]  atoms       : Atoms type, see types.f90
  !> @param[in]  cell        : Unit cell type, see types.f90
  !> @param[in]  iDtype      : Index of atom type to which displaced atom belongs
  !> @param[in]  iDatom      : Index of atom that is displaced
  !> @param[in]  ngdp        : Number of G-vectors for potentials and densities
  !> @param[in]  gdp         : G-vectors of potentials and densities
  !> @param[in]  rho0IRpw    : Plane-wave coefficients of the unperturbed and converged interstitial density parsed from Fleur
  !> @param[in]  rho0MTsh    : Spherical harmonic coefficients of unperturbed and converged muffin-tin density parsed from Fleur
  !> @param[out] qlmHartSurf : Multipole moments for the gradient or the linear variation of the Hartree potential that result from
  !>                           a discontinuity of the unperturbed density
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine CalcQlmHarSurf( atoms, cell, iDtype, iDatom, ngdp, gdp, rho0IRpw, rho0MTsh, qlmHartSurf )

    use m_types, only : t_atoms, t_cell

    implicit none


    ! Type parameters
    type(t_atoms),             intent(in)  :: atoms
    type(t_cell),              intent(in)  :: cell

    ! Scalar parameter
    integer,                   intent(in)  :: iDtype
    integer,                   intent(in)  :: iDatom
    integer,                   intent(in)  :: ngdp

    ! Array parameters
    integer,                   intent(in)  :: gdp(:, :)
    complex,                   intent(in)  :: rho0IRpw(:, :)
    complex,                   intent(in)  :: rho0MTsh(:, :, :, :)
    complex,                   intent(out) :: qlmHartSurf(:, :)

    ! Scalar variables
    integer                                :: oqn_l
    integer                                :: mqn_m
    integer                                :: lm
    integer                                :: idir

    ! Array variables
    complex,       allocatable             :: qlmHartSurfIR(:, :)
    complex,       allocatable             :: qlmHartSurfMT(:, :)

    qlmHartSurf(:, :) = cmplx(0., 0.)

    allocate( qlmHartSurfIR(3, (atoms%lmax(iDtype) + 1)**2) )
    allocate( qlmHartSurfMT(3, (atoms%lmax(iDtype) + 1)**2) )
    qlmHartSurfIR(:, :) = cmplx(0., 0.)
    qlmHartSurfMT(:, :) = cmplx(0., 0.)

    call CalcQlmHarSurfIR( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw(:, 1), qlmHartSurfIR )
    call CalcQlmHarSurMT(atoms, iDtype, iDatom, rho0MTsh, qlmHartSurfMT)

    do oqn_l = 0, atoms%lmax(iDtype)
      do mqn_m = -oqn_l, oqn_l
        lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
        do idir = 1, 3
          qlmHartSurf(lm, idir) = qlmHartSurfMT(idir, lm)
          qlmHartSurf(lm, idir) = qlmHartSurf(lm, idir) - qlmHartSurfIR(idir, lm)
          if (.false.) then
            write(2401, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurfIR(idir, lm)
            write(2402, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurfMT(idir, lm)
            write(2403, '(3i8,2f15.8)') oqn_l, mqn_m, idir, qlmHartSurf(lm, idir)
          end if
        end do ! idir
      end do ! mqn_m
    end do ! oqn_l

  end subroutine CalcQlmHarSurf

end module m_jpPotDensHelper
