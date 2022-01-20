! This is module is thought for routines to plot calculated observables for production
module m_jpPlotObservables

  implicit none

  contains

    subroutine PlotGradV0( atoms, cell, lathar, stars, dimens, input, sym, memd_atom, ngdp, logUnit, plot3DgrVhar0Sw, &
        & plot3DgrVext0Sw, plot3DgrVxc0Sw, plot3DgrVeff0Sw, gdp, rho0IRst, rho0MTlh, nmem_atom, mlh_atom, clnu_atom, grRho0IR,&
        & grRho0MT, gaussWghts, ylm, dKernMTGPts, grVxcIRKern )

      use m_types
      use m_jpPotDensHelper, only : convertStar2G

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_stars),                  intent(in)  :: stars
      type(t_dimension),              intent(in)  :: dimens
      type(t_input),                  intent(in)  :: input
      type(t_sym),                    intent(in)  :: sym

      ! Scalar parameters
      integer,                        intent(in)  :: memd_atom
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: logUnit
      logical,                        intent(in)  :: plot3DgrVhar0Sw
      logical,                        intent(in)  :: plot3DgrVext0Sw
      logical,                        intent(in)  :: plot3DgrVxc0Sw
      logical,                        intent(in)  :: plot3DgrVeff0Sw


      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      complex,                        intent(in)  :: rho0IRst(:, :) !n3d !todo  with spin
      real,                           intent(in)  :: rho0MTlh(:, 0:, :, :) ! todo with spin
      integer,                        intent(in)  :: nmem_atom(0:, :)
      integer,                        intent(in)  :: mlh_atom(:,0:,:)
      complex,                        intent(in)  :: clnu_atom(:,0:,:)
      complex,                        intent(in)  :: grRho0IR(:, :)
      complex,                        intent(in)  :: grRho0MT(:, :, :, :)
      real,                           intent(in)  :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in)  :: ylm(:, : )
      real,                           intent(in)  :: dKernMTGPts(:, :, :)
      complex,                        intent(in)  :: grVxcIRKern(:)

      ! Scalar variables
      integer                                     :: iatom
      integer                                     :: itype
      integer                                     :: ieqat
      integer                                     :: ptsym
      integer                                     :: ilh
      integer                                     :: oqn_l
      integer                                     :: imem
      integer                                     :: mqn_m
      integer                                     :: lm_pre
      integer                                     :: lm
      integer                                     :: imesh

      ! Array variables
      complex,           allocatable              :: rho0IRpw(:, :)
      complex,           allocatable              :: rho0MTsh(:, :, :, :)

      allocate(  rho0IRpw(ngdp, 1) )
      rho0IRpw(:, :) = cmplx(0., 0.)
      call convertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

      allocate( rho0MTsh( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      rho0MTsh(:, :, :, :) = cmplx(0., 0.)
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = lm_pre + mqn_m
              !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
              ! maybe construct a pointer and run only over them to make it memory efficient.
              do imesh = 1, atoms%jri(itype)
                rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MTlh(imesh, ilh, itype, 1) &
                                                                                                      & * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
        end do ! ieqat
      end do ! itype

      if ( plot3DgrVhar0Sw ) then
      write(*, *) 'plot3DgrVhar0'
        call Plot3DgrVhar0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, rho0IRpw, rho0MTsh, grRho0IR, &
          & grRho0MT )
      end if ! plot3DgrVhar0

      if ( plot3DgrVext0Sw ) then
      write(*, *) 'plot3DgrVext0'
        call Plot3DgrVext0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp )
      end if ! plot3DgrVext0

      if ( plot3DgrVxc0Sw ) then
      write(*, *) 'plot3DgrVxc0'
        call Plot3DgrVxc0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, grRho0IR, grRho0MT, gaussWghts, ylm,&
          & dKernMTGpts, grVxcIRKern )
      end if ! plot3DgrVxc0

      if ( plot3DgrVeff0Sw ) then
      write(*, *) 'plot3DgrVeff0'
        call Plot3DgrVeff0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, rho0IRpw, rho0MTsh, grRho0IR, &
          & grRho0MT, gaussWghts, ylm, dKernMTGpts, grVxcIRKern )
      end if ! plot3DgrVeff0

    end subroutine PlotGradV0

    subroutine plotConvergedQuantities( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, memd_atom, logUnit, &
        & plotVhar1fccXZSw, plotVext1fccXZSw, plotVxc1fccXZSw, plotVeff1fccXZSw, gdp, rho0IRst, rho0MTlh, nmem_atom, mlh_atom, &
        & clnu_atom, rho1IR, rho1MTfull, gaussWghts, ylm, dKernMTGPts, vxc1IRKern, rho1MTDelta, plotRho1fccXZSw )

      use m_types
      use m_jpPotDensHelper, only : ConvertStar2G, genPertPotDensGvecs

      implicit none

      type(t_atoms),                  intent(in) :: atoms
      type(t_stars),                  intent(in) :: stars
      type(t_cell),                   intent(in) :: cell
      type(t_sphhar),                 intent(in) :: lathar
      type(t_dimension),              intent(in) :: dimens
      type(t_sym),                    intent(in) :: sym
      type(t_input),                  intent(in) :: input
      type(t_kpts),                   intent(in) :: qpts

      ! Scalar paramters
      integer,                        intent(in) :: ngdp
      integer,                        intent(in) :: iqpt
      integer,                        intent(in) :: memd_atom
      integer,                        intent(in) :: logUnit
      logical,                        intent(in) :: plotVhar1fccXZSw
      logical,                        intent(in) :: plotVext1fccXZSw
      logical,                        intent(in) :: plotVxc1fccXZSw
      logical,                        intent(in) :: plotVeff1fccXZSw
      logical,                        intent(in) :: plotRho1fccXZSw

      ! Array parameters
      integer,                        intent(in) :: gdp(:, :)
      complex,                        intent(in) :: rho0IRst(:, :) !n3d !todo  with spin
      real,                           intent(in) :: rho0MTlh(:, 0:, :, :) ! todo with spin
      integer,                        intent(in) :: nmem_atom(0:, :)
      integer,                        intent(in) :: mlh_atom(:,0:,:)
      complex,                        intent(in) :: clnu_atom(:,0:,:)
      complex,                        intent(in) :: rho1IR(:, :, :)
      complex,                        intent(in) :: rho1MTfull(:, :, :, :, :)
      real,                           intent(in) :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in) :: ylm(:, : )
      real,                           intent(in) :: dKernMTGPts(:, :, :)
      complex,                        intent(in) :: vxc1IRKern(:)
      complex,                        intent(in) :: rho1MTDelta(:, :, :, :, :)

      integer                                    :: ngpqdp
      integer                                    :: ngpqdp2km
      integer,           allocatable             :: gpqdp2Ind(:, :, :)
      integer,           allocatable             :: gpqdp(:, :)
      complex,           allocatable             :: rho0IRpw(:, :)
      complex,           allocatable             :: rho0MTsh(:, :, :, :)
      integer                                    :: itype
      integer                                    :: ieqat
      integer                                    :: iatom
      integer                                     :: ptsym
      integer                                     :: ilh
      integer                                     :: oqn_l
      integer                                     :: imem
      integer                                     :: mqn_m
      integer                                     :: lm_pre
      integer                                     :: lm
      integer                                     :: imesh
      integer                                    :: gpqdp2iLim(2, 3)

      allocate(  rho0IRpw(ngdp, 1) )
      rho0IRpw(:, :) = cmplx(0., 0.)
      call convertStar2G( rho0IRst(:, 1), rho0IRpw(:, 1), stars, ngdp, gdp )

      allocate( rho0MTsh( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      rho0MTsh(:, :, :, :) = cmplx(0., 0.)
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = lm_pre + mqn_m
              !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
              ! maybe construct a pointer and run only over them to make it memory efficient.
              do imesh = 1, atoms%jri(itype)
                rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MTlh(imesh, ilh, itype, 1) &
                                                                                                      & * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
        end do ! ieqat
      end do ! itype
      !call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpts%bk(1:3, iqpt), gpqdp, gpqdp2Ind, gpqdp2iLim )
      call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpts%bk(1:3, 1), gpqdp, gpqdp2Ind, gpqdp2iLim )
      deallocate(gpqdp2Ind)

      if (plotRho1fccXZSw) then
        call plotRho1Final( atoms, dimens, stars, lathar, input, sym, cell, qpts, iqpt, ngdp, gdp, rho1IR, rho1MTfull, rho1MTDelta )
      end if

      if (plotVhar1fccXZSw) then
      write(*, *) 'Plot Vhar1'
        call PlotVhar1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
          & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho0IRpw, rho0MTsh, rho1IR, rho1MTFull )
      end if

      if (plotVext1fccXZSw) then
      write(*, *) 'Plot Vext1'
        call PlotVext1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
          & gpqdp, nmem_atom, mlh_atom, clnu_atom )
      end if

      if (plotVxc1fccXZSw) then
      write(*, *) 'Plot Vxc1'
        call PlotVxc1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
          & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho1IR, rho1MTFull, gaussWghts, ylm, dKernMTGPts, vxc1IRKern )
      end if

      if (plotVeff1fccXZSw) then
      write(*, *) 'Plot Veff1'
        call PlotVeff1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
          & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho0IRpw, rho0MTsh, rho1IR, rho1MTFull, gaussWghts, ylm, dKernMTGPts, vxc1IRKern)
      end if

    end subroutine plotConvergedQuantities

    subroutine plotRho1Final( atoms, dimens, stars, lathar, input, sym, cell, qpts, iqpt, ngdp, gdp, rho1IRDS, rho1MTfull, rho1MTDelta )

      use m_types

      implicit none

      type(t_atoms),                  intent(in) :: atoms
      type(t_dimension),              intent(in) :: dimens
      type(t_stars),                  intent(in) :: stars
      type(t_sphhar),                 intent(in) :: lathar
      type(t_input),                  intent(in) :: input
      type(t_sym),                    intent(in) :: sym
      type(t_cell),                   intent(in) :: cell
      type(t_kpts),                   intent(in) :: qpts

      integer,                        intent(in) :: iqpt
      integer,                        intent(in) :: ngdp

      integer,                        intent(in) :: gdp(:, :)
      complex,                        intent(in) :: rho1IRDS(:, :, :)
      complex,                        intent(in) :: rho1MTDelta(:, :, :, :, :)
      complex,                        intent(in) :: rho1MTFull(:, :, :, :, :)

      integer                                    :: iDtype
      integer                                    :: iDeqat
      integer                                    :: iDatom
      integer                                    :: idir
      character(len=25)                          :: plotname

      write(*, *) 'Plot full charge density variation'
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do idir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotRho1fullVesta', idir, 'at',iDatom, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, qpts%bk(:, iqpt), &
                                    & rho1IRDS(:, idir, iDatom), rho1MTFull(:, :, :, idir, iDatom), plotname, iqpt, .FALSE., "plot_inp")
          end do ! idir
        end do ! iDeqat
      end do ! iDtype

      write(*, *) 'Plot charge density variation without gradient'
      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          do idir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotRho1NoGrVesta', idir, 'at',iDatom, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, qpts%bk(:, iqpt), &
                                    & rho1IRDS(:, idir, iDatom), rho1MTDelta(:, :, :, idir, iDatom), plotname, iqpt, .FALSE., "plot_inp")
          end do ! idir
        end do ! iDeqat
      end do ! iDtype

    end subroutine plotRho1Final

    subroutine PlotVhar1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
        & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho0IRpw, rho0MTsh, rho1IR, rho1MTfull )

      use m_types
      use m_jpVeff1, only : GenVeff1

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_stars),                  intent(in)  :: stars
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_dimension),              intent(in)  :: dimens
      type(t_sym),                    intent(in)  :: sym
      type(t_input),                  intent(in)  :: input
      type(t_kpts),                   intent(in)  :: qpts

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: iqpt
      integer,                        intent(in)  :: ngpqdp
      integer,                        intent(in)  :: memd_atom
      integer,                        intent(in)  :: logUnit

      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      integer,                        intent(in)  :: gpqdp(:, :)
      integer,                        intent(in)  :: nmem_atom(0:, :)
      integer,                        intent(in)  :: mlh_atom(:,0:,:)
      complex,                        intent(in)  :: clnu_atom(:,0:,:)
      complex,                        intent(in)  :: rho0IRpw(:, :)
      complex,                        intent(in)  :: rho0MTsh(:, :, :, :)
      complex,                        intent(in)  :: rho1IR(:, :, :)
      complex,                        intent(in)  :: rho1MTfull(:, :, :, :, :)

      ! Scalar variables
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      logical                                     :: vExtFull
      logical                                     :: vHarNum
      logical                                     :: expndTillLmaxP1
      integer                                     :: iDdir
      integer                                     :: iDtype
      integer                                     :: iDeqat
      integer                                     :: iDatom

      ! Array variables
      complex,           allocatable              :: vHar1IR(:, :)
      complex,           allocatable              :: vHar1MT(:, :, :, :)
      complex,           allocatable              :: rho0IRDummy(:, :)
      complex,           allocatable              :: rho0MTDummy(:, :, :, :)
      complex,           allocatable              :: ylmDummy(:, :)
      complex,           allocatable              :: vxc1IRKernDummy(:)
      real,              allocatable              :: dKernMTGPtsDummy(:, :, :)
      real,              allocatable              :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,           allocatable              :: grRho0MTZero(:, :, :, :)
      character(len=29)                           :: plotname


      ! Create external potential
      harSw = .true.
      extSw = .false.
      xcSw = .false.
      ! There is nothing canceling away as in Sternheimer, therefore we need the full contribution of the external potential
      vExtFull = .true.
!todo we don't need rho0IR and rho0MT in the argument list
      vHarNum = .false.
      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      allocate( grRho0MTZero(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      allocate( gaussWghts(dimens%nspd) )
      allocate( ylmDummy(dimens%nspd, ( atoms%lmaxd + 1 )**2) )
      allocate( vxc1IRKernDummy(ngdp) )
      allocate( dKernMTGPtsDummy(dimens%nspd, atoms%jmtd, atoms%nat) )
     allocate(  rho0IRDummy(ngdp, 1) )
     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      rho0IRDummy(:, :) = cmplx(0., 0.)
      rho0MTDummy(:, :, :, :) = 0.
      ! rho1MTfull contains already the gradient added to the density, as this routine is actually developed for the Sternheimer SCC
      grRho0MTZero(:, :, :, :) = cmplx(0., 0.)

      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpts%bk(1:3, iqpt), rho0IRpw, rho0MTsh, &
            & rho1IR(:, :, iDatom), rho1MTfull(:, :, :, :, iDatom), grRho0MTZero, gdp, vHar1IR, vHar1MT, vxc1IRKernDummy, ylmDummy,&
            & dKernMTGPtsDummy, gaussWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

          do iDdir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotVhar1fccXZ2DdAt', iDatom, 'dDir', iDdir, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gpqdp, ngpqdp, qpts%bk(1:3, iqpt), &
                                          & vHar1IR(:, iDdir), vHar1MT(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp")
          end do
        end do ! iDeqat
      end do ! iDtype

    end subroutine PlotVhar1fccXZ

    subroutine PlotVext1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, &
        & gdp, gpqdp, nmem_atom, mlh_atom, clnu_atom )

      use m_types
      use m_jpVeff1, only : GenVeff1

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_stars),                  intent(in)  :: stars
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_dimension),              intent(in)  :: dimens
      type(t_sym),                    intent(in)  :: sym
      type(t_input),                  intent(in)  :: input
      type(t_kpts),                   intent(in)  :: qpts

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: iqpt
      integer,                        intent(in)  :: ngpqdp
      integer,                        intent(in)  :: memd_atom
      integer,                        intent(in) :: logUnit

      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      integer,                        intent(in)  :: gpqdp(:, :)
      integer,                        intent(in)  :: nmem_atom(0:, :)
      integer,                        intent(in)  :: mlh_atom(:,0:,:)
      complex,                        intent(in)  :: clnu_atom(:,0:,:)

      ! Scalar variables
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      logical                                     :: vExtFull
      logical                                     :: vHarNum
      logical                                     :: expndTillLmaxP1
      integer                                     :: iDdir
      integer                                     :: iDtype
      integer                                     :: iDeqat
      integer                                     :: iDatom

      ! Array variables
      complex,           allocatable              :: vExt1IR(:, :)
      complex,           allocatable              :: vExt1MT(:, :, :, :)
      complex,           allocatable              :: rho1PWzero(:, :)
      complex,           allocatable              :: rho1MTzero(:, :, :, :)
      complex,           allocatable              :: rho0IRDummy(:, :)
      complex,           allocatable              :: rho0MTDummy(:, :, :, :)
      complex,           allocatable              :: ylmDummy(:, :)
      complex,           allocatable              :: vxc1IRKernDummy(:)
      real,              allocatable              :: dKernMTGPtsDummy(:, :, :)
      real,              allocatable              :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,           allocatable              :: grRho0MTDummy(:, :, :, :)
      character(len=29)                           :: plotname


      ! Create external potential
      harSw = .false.
      extSw = .true.
      xcSw = .false.
      ! There is nothing canceling away as in Sternheimer, therefore we need the full contribution of the external potential
      vExtFull = .true.
!todo we don't need rho0IR and rho0MT in the argument list
      vHarNum = .false.
      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      ! Actually we do not need it really
      allocate(rho1PWzero(ngpqdp, 3))
      allocate( rho1MTzero(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3))
      rho1PWzero = cmplx(0.0, 0.0)
      rho1MTzero = cmplx(0.0, 0.0)

     allocate(  rho0IRDummy(ngdp, 1) )

     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      allocate( grRho0MTDummy(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      allocate( gaussWghts(dimens%nspd) )
      allocate( ylmDummy(dimens%nspd, ( atoms%lmaxd + 1 )**2) )
      allocate( vxc1IRKernDummy(ngdp) )
      allocate( dKernMTGPtsDummy(dimens%nspd, atoms%jmtd, atoms%nat) )
      grRho0MTDummy(:, :, :, :) = cmplx(0., 0.)
     rho0MTDummy(:, :, :, :) = cmplx(0., 0.)
     rho0IRDummy(:, :) = cmplx(0., 0.)
     gaussWghts(:) = 0.
     ylmDummy(:, :) = cmplx(0., 0.)
     vxc1IRKernDummy(:) = cmplx(0.)
     dKernMTGPtsDummy(:, :, :) = cmplx(0., 0.)

      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpts%bk(1:3, iqpt), rho0IRDummy, &
            & rho0MTDummy, rho1PWzero, rho1MTzero(:, :, :, :), grRho0MTDummy, gdp, vExt1IR, vExt1MT, vxc1IRKernDummy, ylmDummy, &
            & dKernMTGPtsDummy, gaussWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

          do iDdir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotVext1fccXZ2DdAt', iDatom, 'dDir', iDdir, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gpqdp, ngpqdp, qpts%bk(1:3, iqpt), &
                                          & vExt1IR(:, iDdir), vExt1MT(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp")
          end do
        end do ! iDeqat
      end do ! iDtype

    end subroutine PlotVext1fccXZ

    subroutine PlotVxc1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp,   &
        & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho1IR, rho1MTfull, gaussWghts, ylm, dKernMTGPts, vxc1IRKern )

      use m_types
      use m_jpVeff1, only : GenVeff1

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_stars),                  intent(in)  :: stars
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_dimension),              intent(in)  :: dimens
      type(t_sym),                    intent(in)  :: sym
      type(t_input),                  intent(in)  :: input
      type(t_kpts),                   intent(in)  :: qpts

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: iqpt
      integer,                        intent(in)  :: ngpqdp
      integer,                        intent(in)  :: memd_atom
      integer,                        intent(in) :: logUnit

      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      integer,                        intent(in)  :: gpqdp(:, :)
      integer,                        intent(in)  :: nmem_atom(0:, :)
      integer,                        intent(in)  :: mlh_atom(:,0:,:)
      complex,                        intent(in)  :: clnu_atom(:,0:,:)
      complex,                        intent(in)  :: rho1IR(:, :, :)
      complex,                        intent(in)  :: rho1MTfull(:, :, :, :, :)
      real,                           intent(in)  :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in)  :: ylm(:, : )
      real,                           intent(in)  :: dKernMTGPts(:, :, :)
      complex,                        intent(in)  :: vxc1IRKern(:)

      ! Scalar variables
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      logical                                     :: vExtFull
      logical                                     :: vHarNum
      logical                                     :: expndTillLmaxP1
      integer                                     :: iDdir
      integer                                     :: iDtype
      integer                                     :: iDeqat
      integer                                     :: iDatom
      integer :: idir
      integer :: iG
      integer :: oqn_l
      integer :: mqn_m
      integer :: imesh
      integer :: lm

      ! Array variables
      complex,           allocatable              :: vxc1IR(:, :)
      complex,           allocatable              :: vxc1MT(:, :, :, :)
      complex,           allocatable              :: rho0IRDummy(:, :)
      complex,           allocatable              :: rho0MTDummy(:, :, :, :)
      complex,           allocatable              :: grRho0MTZero(:, :, :, :)
      character(len=29)                           :: plotname


      ! Create external potential
      harSw = .false.
      extSw = .false.
      xcSw = .true.
      ! There is nothing canceling away as in Sternheimer, therefore we need the full contribution of the external potential
      ! vExtFull is actually not required for the xc potential, provided one uses the charge density variation with the gradient addeed
      vExtFull = .true.
!todo we don't need rho0IR and rho0MT in the argument list
      vHarNum = .false.
      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      allocate(  rho0IRDummy(ngdp, 1) )
      allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      allocate( grRho0MTZero(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      rho0MTDummy(:, :, :, :) = cmplx(0., 0.)
      rho0IRDummy(:, :) = cmplx(0., 0.)
      ! rho1MTfull contains already the gradient added to the density, as this routine is actually developed for the Sternheimer SCC
      grRho0MTZero(:, :, :, :) = cmplx(0., 0.)

      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpts%bk(1:3, iqpt), rho0IRDummy, rho0MTDummy, &
            & rho1IR(:, :, iDatom), rho1MTfull(:, :, :, :, iDatom), grRho0MTZero, gdp, vxc1IR, vxc1MT, vxc1IRKern, ylm, &
            & dKernMTGPts, gaussWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

!
!      write(*, *) qpts%bk(1:3, iqpt)
!      write(*, *) ngdp
!      write(*, *) ngpqdp
!      write(5423, '(f15.8)') vxc1IRKern
!      do idir = 1, 3
!        do iG = 1, ngdp
!          write(5421, '(2i8, 2f15.8)') idir, iG, vxc1IR(iG, idir)
!          write(5424, '(2i8, 2f15.8)') idir, iG, rho1IR(iG, idir, 1)
!          write(5425, '(2i8, 2f15.8)') idir, iG, gdp(idir, iG)
!          write(5426, '(2i8, 2f15.8)') idir, iG, gpqdp(idir, iG)
!        end do ! iG
!      end do ! idir
!          do idir = 1, 3
!             do oqn_l = 0, atoms%lmax(iDtype)
!               do mqn_m = -oqn_l, oqn_l
!                 lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                 do imesh = 1, atoms%jri(iDtype)
!                   write(5422, '(3i8,2f15.8)') idir, lm, imesh, vxc1MT(imesh, lm, iDatom, idir)
!                 end do
!               end do
!             end do
!          end do
!
          do iDdir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotVxc1fccXZ2DdAt', iDatom, 'dDir', iDdir, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gpqdp, ngpqdp, qpts%bk(1:3, iqpt), &
                                          & vxc1IR(:, iDdir), vxc1MT(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp")
          end do
        end do ! iDeqat
      end do ! iDtype

    end subroutine PlotVxc1fccXZ

    subroutine PlotVeff1fccXZ( atoms, stars, cell, lathar, dimens, sym, input, qpts, ngdp, iqpt, ngpqdp, memd_atom, logUnit, gdp, &
        & gpqdp, nmem_atom, mlh_atom, clnu_atom, rho0IRpw, rho0MTsh, rho1IR, rho1MTfull, gaussWghts, ylm, dKernMTGPts, vxc1IRKern )

      use m_types
      use m_jpVeff1, only : GenVeff1

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_stars),                  intent(in)  :: stars
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_dimension),              intent(in)  :: dimens
      type(t_sym),                    intent(in)  :: sym
      type(t_input),                  intent(in)  :: input
      type(t_kpts),                   intent(in)  :: qpts

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: iqpt
      integer,                        intent(in)  :: ngpqdp
      integer,                        intent(in)  :: memd_atom
      integer,                        intent(in) :: logUnit

      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      integer,                        intent(in)  :: gpqdp(:, :)
      integer,                        intent(in)  :: nmem_atom(0:, :)
      integer,                        intent(in)  :: mlh_atom(:,0:,:)
      complex,                        intent(in)  :: clnu_atom(:,0:,:)
      complex,                        intent(in)  :: rho0IRpw(:, :)
      complex,                        intent(in)  :: rho0MTsh(:, :, :, :)
      complex,                        intent(in)  :: rho1IR(:, :, :)
      complex,                        intent(in)  :: rho1MTfull(:, :, :, :, :)
      real,                           intent(in)  :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in)  :: ylm(:, : )
      real,                           intent(in)  :: dKernMTGPts(:, :, :)
      complex,                        intent(in)  :: vxc1IRKern(:)

      ! Scalar variables
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      logical                                     :: vExtFull
      logical                                     :: vHarNum
      logical                                     :: expndTillLmaxP1
      integer                                     :: iDdir
      integer                                     :: iDtype
      integer                                     :: iDeqat
      integer                                     :: iDatom

      ! Array variables
      complex,           allocatable              :: vxc1IR(:, :)
      complex,           allocatable              :: vxc1MT(:, :, :, :)
      complex,           allocatable              :: rho0IRDummy(:, :)
      complex,           allocatable              :: rho0MTDummy(:, :, :, :)
      complex,           allocatable              :: grRho0MTZero(:, :, :, :)
      character(len=29)                           :: plotname


      ! Create external potential
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      ! There is nothing canceling away as in Sternheimer, therefore we need the full contribution of the external potential
      vExtFull = .true.
!todo we don't need rho0IR and rho0MT in the argument list
      vHarNum = .false.
      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      allocate( grRho0MTZero(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
     allocate(  rho0IRDummy(ngdp, 1) )
     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      rho0IRDummy(:, :) = cmplx(0., 0.)
      rho0MTDummy(:, :, :, :) = 0.
      ! rho1MTfull contains already the gradient added to the density, as this routine is actually developed for the Sternheimer SCC
      grRho0MTZero(:, :, :, :) = cmplx(0., 0.)

      iDatom = 0
      do iDtype = 1, atoms%ntype
        do iDeqat = 1, atoms%neq(iDtype)
          iDatom = iDatom + 1
          call GenVeff1( stars, cell, atoms, dimens, harSw, extSw, xcSw, vExtFull, ngdp, qpts%bk(1:3, iqpt), rho0IRDummy, &
            & rho0MTDummy, rho1IR(:, :, iDatom), rho1MTfull(:, :, :, :, iDatom), grRho0MTZero, gdp, vxc1IR, vxc1MT, vxc1IRKern, &
            & ylm, dKernMTGPts, gaussWghts, iDatom, iDtype, iqpt, ngpqdp, gpqdp, vHarNum ) ! add spin

          do iDdir = 1, 3
            plotname = ''
            write(plotname, '(a,i1,a,i1,a)') 'plotVeff1fccXZ2DdAt', iDatom, 'dDir', iDdir, '.xsf'
            ! Create file of density variation for Vesta
            call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gpqdp, ngpqdp, qpts%bk(1:3, iqpt), &
                                          & vxc1IR(:, iDdir), vxc1MT(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp")
          end do
        end do ! iDeqat
      end do ! iDtype

    end subroutine PlotVeff1fccXZ


    ! Note: If one wants to plot the gradient as it is given by Aaron, one has to construct a new quantitiy for polyatomic systems.
    !       This routine just always gives the gradient term, i.e. the volume term in the muffin-tins. As actually, it is not used
    !       anywhere and does not make sense for the gradient to displace atom and one not, we always plot all atoms displaced.
    subroutine Plot3DgrVext0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp )

      use m_types
      use m_jpGrVeff0, only : GenGrVeff0

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_stars),                  intent(in)  :: stars
      type(t_dimension),              intent(in)  :: dimens
      type(t_input),                  intent(in)  :: input
      type(t_sym),                    intent(in)  :: sym

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: logUnit


      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)

      ! Scalar variables
      logical                                     :: testGoldstein
      logical                                     :: gradMtTermSw
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      integer                                     :: iDatom
      integer                                     :: iDdir
      integer                                     :: iqpt
      logical                                     :: expndTillLmaxP1
      integer                                     :: itype
      integer                                     :: iatom
      integer                                     :: ieqat
      integer                                     :: lm
      integer                                     :: oqn_l
      integer                                     :: mqn_m
      integer                                     :: imesh

      ! Array variables
      complex,           allocatable              :: grVext0IR(:, :)
      complex,           allocatable              :: grVext0MT(:, :, :, :)
      complex,           allocatable              :: grRho0IRDummy(:, :)
      complex,           allocatable              :: grRho0MTDummy(:, :, :, :)
      real,              allocatable              :: gaussWghtsDummy(:) ! gaussian weights belonging to gausPts
      complex,           allocatable              :: ylmDummy(:, : )
      real,              allocatable              :: dKernMTGPtsDummy(:, :, :)
      complex,           allocatable              :: grVxcIRKernDummy(:)
      complex,           allocatable              :: grVext0MTplot(:, :, :, :)
      complex,           allocatable              :: rho0IRpwDummy(:, :)
      complex,           allocatable              :: rho0MTshDummy(:, :, :, :)
      character(len=18)                           :: plotname

      allocate( grRho0IRDummy(ngdp, 3), grRho0MTDummy( atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      allocate( gaussWghtsDummy(dimens%nspd) )
      allocate( ylmDummy(dimens%nspd, ( atoms%lmaxd + 1 )**2) )
      allocate( dKernMTGPtsDummy(dimens%nspd, atoms%jmtd, atoms%nat) )
      allocate( grVxcIRKernDummy(ngdp) )
      allocate( grVext0MTplot(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
      allocate( rho0IRpwDummy(ngdp, 1) )
      allocate( rho0MTshDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )

      grRho0IRDummy(:, :) = cmplx(0., 0.)
      grRho0MTDummy(:, :, :, :) = cmplx(0., 0.)
      gaussWghtsDummy(:) = 0.
      ylmDummy(:, :) = cmplx(0., 0.)
      dKernMTGPtsDummy(:, :, :) = 0.
      grVxcIRKernDummy(:) = 0.
      grVext0MTplot(:, :, :, :) = cmplx(0., 0.)
      rho0IRpwDummy(:, :) = cmplx(0., 0.)
      rho0MTshDummy(:, :, :, :) = cmplx(0., 0.)

      harSw = .false.
      extSw = .true.
      xcSw = .false.
      testGoldstein = .false.
      gradMtTermSw = .true.
      call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpwDummy, rho0MTshDummy, grRho0IRDummy, &
        & grRho0MTDummy, gaussWghtsDummy, ylmDummy, dKernMTGPtsDummy, grVxcIRKernDummy, testGoldstein, gradMtTermSw, grVext0IR, &
        & grVext0MT ) ! add spin

      ! Reorder grVextMT
      do iDdir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVext0MTplot(imesh, lm, iatom, iDdir) = grVext0MT(imesh, lm, iDdir, iatom)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! iDdir

      ! We deal with q = (0, 0, 0)
      iqpt = 1

      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      do iDdir = 1, 3
        plotname = ''
        write(plotname, '(a,i1,a)') 'plotgrVext03D', iDdir, '.xsf'
        ! Create file of density variation for Vesta
        call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, [0., 0., 0.], &
                                 & grVext0IR(:, iDdir), grVext0MTplot(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp3D")
      end do

    end subroutine Plot3DgrVext0

    ! Note: If one wants to plot the gradient as it is given by Aaron, one has to construct a new quantitiy for polyatomic systems.
    !       This routine just always gives the gradient term, i.e. the volume term in the muffin-tins. As actually, it is not used
    !       anywhere and does not make sense for the gradient to displace atom and one not, we always plot all atoms displaced.
    subroutine Plot3DgrVhar0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, rho0IRpw, rho0MTsh, grRho0IR, &
        & grRho0MT )

      use m_types
      use m_jpGrVeff0, only : GenGrVeff0

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_stars),                  intent(in)  :: stars
      type(t_dimension),              intent(in)  :: dimens
      type(t_input),                  intent(in)  :: input
      type(t_sym),                    intent(in)  :: sym

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: logUnit


      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      complex,                        intent(in)  :: rho0IRpw(:, :)
      complex,                        intent(in)  :: rho0MTsh(:,  :, :, :)
      complex,                        intent(in)  :: grRho0IR(:, :)
      complex,                        intent(in)  :: grRho0MT(:, :, :, :)

      ! Scalar variables
      logical                                     :: testGoldstein
      logical                                     :: gradMtTermSw
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      integer                                     :: iDatom
      integer                                     :: iDdir
      integer                                     :: iqpt
      logical                                     :: expndTillLmaxP1
      integer                                     :: itype
      integer                                     :: iatom
      integer                                     :: ieqat
      integer                                     :: lm
      integer                                     :: oqn_l
      integer                                     :: mqn_m
      integer                                     :: imesh

      ! Array variables
      complex,           allocatable              :: grVhar0IR(:, :)
      complex,           allocatable              :: grVhar0MT(:, :, :, :)
      real,              allocatable              :: gaussWghtsDummy(:) ! gaussian weights belonging to gausPts
      complex,           allocatable              :: ylmDummy(:, : )
      real,              allocatable              :: dKernMTGPtsDummy(:, :, :)
      complex,           allocatable              :: grVxcIRKernDummy(:)
      complex,           allocatable              :: grVhar0MTplot(:, :, :, :)
      character(len=18)                           :: plotname

      allocate( gaussWghtsDummy(dimens%nspd) )
      allocate( ylmDummy(dimens%nspd, ( atoms%lmaxd + 1 )**2) )
      allocate( dKernMTGPtsDummy(dimens%nspd, atoms%jmtd, atoms%nat) )
      allocate( grVxcIRKernDummy(ngdp) )
      allocate( grVhar0MTplot(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )

      gaussWghtsDummy(:) = 0.
      ylmDummy(:, :) = cmplx(0., 0.)
      dKernMTGPtsDummy(:, :, :) = 0.
      grVxcIRKernDummy(:) = 0.
      grVhar0MTplot(:, :, :, :) = cmplx(0., 0.)

      harSw = .true.
      extSw = .false.
      xcSw = .false.
      testGoldstein = .false.
      gradMtTermSw = .true.
      ! todo shift allocation of grVext0IR, grVext0MT outwards
      call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, &
        & gaussWghtsDummy, ylmDummy, dKernMTGPtsDummy, grVxcIRKernDummy, testGoldstein, gradMtTermSw, grVhar0IR, grVhar0MT ) ! add spin

      ! Reorder grVharMT
      do iDdir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVhar0MTplot(imesh, lm, iatom, iDdir) = grVhar0MT(imesh, lm, iDdir, iatom)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! iDdir

      ! We deal with q = (0, 0, 0)
      iqpt = 1

      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      do iDdir = 1, 3
        plotname = ''
        write(plotname, '(a,i1,a)') 'plotgrVhar03D', iDdir, '.xsf'
        ! Create file of density variation for Vesta
        call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, [0., 0., 0.], &
                              & grVhar0IR(:, iDdir), grVhar0MTplot(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp3D")
      end do

    end subroutine Plot3DgrVhar0

    ! Note: If one wants to plot the gradient as it is given by Aaron, one has to construct a new quantitiy for polyatomic systems.
    !       This routine just always gives the gradient term, i.e. the volume term in the muffin-tins. As actually, it is not used
    !       anywhere and does not make sense for the gradient to displace atom and one not, we always plot all atoms displaced.
    subroutine Plot3DgrVxc0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, grRho0IR, grRho0MT, gaussWghts, &
        & ylm, dKernMTGpts, grVxcIRKern )

      use m_types
      use m_jpGrVeff0, only : GenGrVeff0

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_stars),                  intent(in)  :: stars
      type(t_dimension),              intent(in)  :: dimens
      type(t_input),                  intent(in)  :: input
      type(t_sym),                    intent(in)  :: sym

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: logUnit


      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      complex,                        intent(in)  :: grRho0IR(:, :)
      complex,                        intent(in)  :: grRho0MT(:, :, :, :)
      real,                           intent(in)  :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in)  :: ylm(:, : )
      real,                           intent(in)  :: dKernMTGPts(:, :, :)
      complex,                        intent(in)  :: grVxcIRKern(:)

      ! Scalar variables
      logical                                     :: testGoldstein
      logical                                     :: gradMtTermSw
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      integer                                     :: iDatom
      integer                                     :: iDdir
      integer                                     :: itype
      integer                                     :: iqpt
      logical                                     :: expndTillLmaxP1
      integer                                     :: iatom
      integer                                     :: ieqat
      integer                                     :: lm
      integer                                     :: oqn_l
      integer                                     :: mqn_m
      integer                                     :: imesh

      ! Array variables
      complex,           allocatable              :: grVxc0IR(:, :)
      complex,           allocatable              :: grVxc0MT(:, :, :, :)
      complex,           allocatable              :: grVxc0MTplot(:, :, :, :)
      complex,           allocatable              :: rho0IRpwDummy(:, :)
      complex,           allocatable              :: rho0MTshDummy(:, :, :, :)
      character(len=17)                           :: plotname

      allocate( rho0IRpwDummy(ngdp, 1) )
      allocate( rho0MTshDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
      allocate( grVxc0MTplot(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )

      grVxc0MTplot(:, :, :, :) = cmplx(0., 0.)
      rho0IRpwDummy(:, :) = cmplx(0., 0.)
      rho0MTshDummy(:, :, :, :) = cmplx(0., 0.)

      harSw = .false.
      extSw = .false.
      xcSw = .true.
      testGoldstein = .false.
      gradMtTermSw = .true.
      call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpwDummy, rho0MTshDummy, grRho0IR, grRho0MT, &
        & gaussWghts, ylm, dKernMTGPts, grVxcIRKern, testGoldstein, gradMtTermSw, grVxc0IR, grVxc0MT ) ! add spin

      ! Reorder grVharMT
      do iDdir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVxc0MTplot(imesh, lm, iatom, iDdir) = grVxc0MT(imesh, lm, iDdir, iatom)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! iDdir

      ! We deal with q = (0, 0, 0)
      iqpt = 1

      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      do iDdir = 1, 3
        plotname = ''
        write(plotname, '(a,i1,a)') 'plotgrVxc03D', iDdir, '.xsf'
        ! Create file of density variation for Vesta
        call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, [0., 0., 0.], &
                                & grVxc0IR(:, iDdir), grVxc0MTplot(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp3D")
      end do

    end subroutine Plot3DgrVxc0

    ! Note: If one wants to plot the gradient as it is given by Aaron, one has to construct a new quantitiy for polyatomic systems.
    !       This routine just always gives the gradient term, i.e. the volume term in the muffin-tins. As actually, it is not used
    !       anywhere and does not make sense for the gradient to displace atom and one not, we always plot all atoms displaced.
    subroutine Plot3DgrVeff0( atoms, cell, lathar, stars, dimens, input, sym, ngdp, logUnit, gdp, rho0IRpw, rho0MTsh, &
        & grRho0IR, grRho0MT, gaussWghts, ylm, dKernMTGpts, grVxcIRKern )

      use m_types
      use m_jpGrVeff0, only : GenGrVeff0

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in)  :: atoms
      type(t_cell),                   intent(in)  :: cell
      type(t_sphhar),                 intent(in)  :: lathar
      type(t_stars),                  intent(in)  :: stars
      type(t_dimension),              intent(in)  :: dimens
      type(t_input),                  intent(in)  :: input
      type(t_sym),                    intent(in)  :: sym

      ! Scalar parameters
      integer,                        intent(in)  :: ngdp
      integer,                        intent(in)  :: logUnit


      ! Array parameters
      integer,                        intent(in)  :: gdp(:, :)
      complex,                        intent(in)  :: rho0IRpw(:, :) !n3d !todo  with spin
      complex,                        intent(in)  :: rho0MTsh(:, :, :, :) ! todo with spin
      complex,                        intent(in)  :: grRho0IR(:, :) ! todo with spin
      complex,                        intent(in)  :: grRho0MT(:, :, :, :)
      real,                           intent(in)  :: gaussWghts(:) ! gaussian weights belonging to gausPts
      complex,                        intent(in)  :: ylm(:, : )
      real,                           intent(in)  :: dKernMTGPts(:, :, :)
      complex,                        intent(in)  :: grVxcIRKern(:)

      ! Scalar variables
      logical                                     :: testGoldstein
      logical                                     :: gradMtTermSw
      logical                                     :: harSw
      logical                                     :: extSw
      logical                                     :: xcSw
      integer                                     :: iDatom
      integer                                     :: iDdir
      integer                                     :: iqpt
      logical                                     :: expndTillLmaxP1
      integer                                     :: iatom
      integer                                     :: ieqat
      integer                                     :: itype
      integer                                     :: lm
      integer                                     :: oqn_l
      integer                                     :: mqn_m
      integer                                     :: imesh

      ! Array variables
      complex,           allocatable              :: grVeff0IR(:, :)
      complex,           allocatable              :: grVeff0MT(:, :, :, :)
      complex,           allocatable              :: grVeff0MTplot(:, :, :, :)
      character(len=18)                           :: plotname

      allocate( grVeff0MTplot(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )

      grVeff0MTplot(:, :, :, :) = cmplx(0., 0.)

      harSw = .true.
      extSw = .true.
      xcSw = .true.
      testGoldstein = .false.
      gradMtTermSw = .true.
      call GenGrVeff0( atoms, cell, dimens, stars, ngdp, harSw, extSw, xcSw, gdp, rho0IRpw, rho0MTsh, grRho0IR, grRho0MT, &
        & gaussWghts, ylm, dKernMTGPts, grVxcIRKern, testGoldstein, gradMtTermSw, grVeff0IR, grVeff0MT ) ! add spin

      ! Reorder grVharMT
      do iDdir = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVeff0MTplot(imesh, lm, iatom, iDdir) = grVeff0MT(imesh, lm, iDdir, iatom)
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! iDdir

      ! We deal with q = (0, 0, 0)
      iqpt = 1

      ! Quantities expanded until lmax and incorporating the MT gradient are expanded until lmax + 1
      expndTillLmaxP1 = .false.

      do iDdir = 1, 3
        plotname = ''
        write(plotname, '(a,i1,a)') 'plotgrVeff03D', iDdir, '.xsf'
        ! Create file of density variation for Vesta
        call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, [0., 0., 0.], &
                                  & grVeff0IR(:, iDdir), grVeff0MTplot(:, :, :, iDdir), plotname, iqpt, expndTillLmaxP1, "plot_inp3D")
      end do

    end subroutine Plot3DgrVeff0

    subroutine plotVestPotDens0Wrap( atoms, stars, sym, cell, dimens, lathar, input, sliceplot, rho0IR, rho0MT )

      use m_types

      implicit none

      ! Type parameter
      type(t_atoms),      intent(in) :: atoms
      type(t_stars),      intent(in) :: stars
      type(t_sym),        intent(in) :: sym
      type(t_cell),       intent(in) :: cell
      type(t_dimension), intent(in) :: dimens
      type(t_sphhar),     intent(in) :: lathar
      type(t_input),      intent(in) :: input
      type(t_sliceplot),  intent(in) :: sliceplot

      complex,            intent(in) :: rho0IR(:, :)
      real,               intent(in) :: rho0MT(:, :, :, :)

      type(t_oneD)                   :: oneD
      type(t_vacuum)                 :: vacuum

      oneD%odi%d1 = .false.
      call plotVestPotDens0(oneD, dimens, stars, vacuum, lathar, atoms, input, sym, cell, sliceplot, rho0IR, rho0MT)

    end subroutine plotVestPotDens0Wrap
      ! if 3 plots are chosen vesta can plot it, we still have to find out whether slices are possible in vesta
      !be careful this routine needs
      SUBROUTINE plotVestPotDens0(oneD,dimension,&
     &     stars,vacuum,sphhar,atoms,&
     &     input,sym,cell,sliceplot,&
     &      qpw, rho)
!    *****************************************************
      USE m_outcdn
      USE m_loddop
      USE m_xsf_io
      use m_juDFT_NOstopNO, only : juDFT_error
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_oneD),INTENT(IN)     :: oneD
      TYPE(t_dimension),INTENT(IN):: dimension
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_vacuum),INTENT(IN)   :: vacuum
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_cell),INTENT(IN)     :: cell
      TYPE(t_sliceplot),INTENT(IN):: sliceplot
      complex, intent(in) :: qpw(:, :)
      real, intent(in) :: rho(:, 0:, :, :)

!     .. Local Scalars ..
      REAL          :: s,tec,qint,xdnout,fermiEnergyTemp
      INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,jsp,na,nplo
      INTEGER       :: nplot,nq,nt,jm,jspin,iter,ii1,ii2
      LOGICAL       :: twodim,newform,l_qfix
!     ..
!     .. Local Arrays ..
      COMPLEX :: qpwCurrent(stars%ng3,input%jspins),rhtxy(vacuum%nmzxyd,stars%ng2-1,2,input%jspins)
      REAL    :: rhoCurrent(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),rht(vacuum%nmzd,2,input%jspins)
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      COMPLEX :: cdom(1),cdomvz(1,1),cdomvxy(1,1,1)
      INTEGER :: grid(3)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(atoms%jmtd)
      REAL    :: point(3)
      CHARACTER (len=30) :: filename
      CHARACTER (len=7)  :: append
      integer :: ispin
      integer :: itype
      integer :: ilh
      integer :: imesh
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename


      INQUIRE(file ="plot_inp",exist= newform)
      IF (.NOT.newform) THEN !no input file exists, create a template and
                            !exit
         OPEN(20,file ="plot_inp")
         WRITE(20,'(i2,a5,l1)') 2,",xsf=",.true.
         WRITE(20,*) "&PLOT twodim=t,cartesian=t"
         WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
         WRITE(20,*) "  filename='plot1' /"
         WRITE(20,*) "&PLOT twodim=f,cartesian=f"
         WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
         WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
         WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
         WRITE(20,*) "  grid(1)=30  grid(2)=30  grid(3)=30  "
         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.0 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No plot_inp file found. Created a template"
         CALL juDFT_error("Missing input for plot; modify plot_inp"&
     &        ,calledby ="plotdop")
      ENDIF

      rhoCurrent = 0.0
      qpwCurrent = 0.0
      qpwCurrent = qpw


      ! todo remove this redundancy
      ! rhoCurrent = rho ! uncomment when redundancy removed
      ! Multiply factor r^2 beacuse it is divided out somewhere else.
      do ispin = 1, input%jspins
        do itype = 1, atoms%ntype
          do ilh = 0,sphhar%nlhd
            do imesh = 1, atoms%jri(itype)
              rhoCurrent(imesh, ilh, itype, ispin) = rho(imesh, ilh, itype, ispin) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
            end do
          end do
        end do
      end do

      !<--Perhaps only the core charge should be plotted
      IF (input%score) THEN
         OPEN (17,file='cdnc',form='unformatted',status='old')
         REWIND 17
         DO jspin = 1,input%jspins
            DO nt = 1,atoms%ntype
               jm = atoms%jri(nt)
               READ (17) (rhocc(i),i=1,jm)
               DO i = 1,atoms%jri(nt)
                  rhoCurrent(i,0,nt,jspin) = rhoCurrent(i,0,nt,jspin) - rhocc(i)/2.0&
     &                 /SQRT( pi_const )
               ENDDO
               READ (17) tec
            ENDDO
            READ (17) qint
            qpwCurrent(1,jspin) = qpwCurrent(1,jspin) - qint/cell%volint
         ENDDO
         CLOSE (17)
      END IF
      !>
      !>
      !<-- Open the plot_inp file for input
      OPEN (18,file='plot_inp')
      READ(18,'(i2,5x,l1)')    nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (xsf) THEN
         OPEN(55,file="plot0.xsf")
         CALL xsf_WRITE_atoms(&
     &                        55,atoms,input%film,oneD%odi%d1,cell%amat)
      ENDIF
      !write(*, *) cell%amat
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (.NOT.twodim.AND.ANY(grid<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(cell%amat,vec1)
            vec2=matmul(cell%amat,vec2)
            vec3=matmul(cell%amat,vec3)
            zero=matmul(cell%amat,zero)
         ENDIF
         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
         IF (xsf) THEN
            CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero&
     &           ,grid)
         ELSE
           OPEN (55,file = filename,form='formatted')
         ENDIF
         !loop over spins
         DO jsp = 1,input%jspins
            !loop over all points
            DO iz = 0,grid(3)-1
               DO iy = 0,grid(2)-1
                  xloop:DO ix = 0,grid(1)-1
                    point = zero + vec1*REAL(ix)/(grid(1)-1) +&
     &                             vec2*REAL(iy)/(grid(2)-1)
                    IF (.NOT.twodim) point = point +&
     &                             vec3*REAL(iz)/(grid(3)-1)
                    !Check if the point is in MT-sphere
                    ii1 = 3
                    ii2 = 3
                    ii3 = 3
                    DO  i1 = -ii1,ii1
                       DO  i2 = -ii2,ii2
                          DO  i3 = -ii3,ii3
                             pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
                             na = 0
                             DO nt = 1,atoms%ntype
                                DO nq = 1,atoms%neq(nt)
                                   na   = na + 1
                                   s  = SQRT(dot_PRODUCT(atoms%pos(:,na)&
     &                                  -pt,atoms%pos(:,na)-pt))
                                   IF (s<atoms%rmsh(atoms%jri(nt),nt)) THEN
                                      CALL outcdn(&
     &                                     pt,nt,na,0,1,jsp,sliceplot,stars,&
     &                                     vacuum,sphhar,atoms,sym,cell,oneD, &
     &                                     qpwCurrent,rhtxy,rhoCurrent,rht,xdnout)
                                      IF (xsf) THEN
                                         write(55,*) xdnout
                                      ELSE
                                            WRITE(55,'(4e15.7)') point ,xdnout
                                      ENDIF
                                      CYCLE xloop
                                   ENDIF
                                ENDDO
                             ENDDO !nt
                          ENDDO
                       ENDDO
                    ENDDO !i1
                    CALL outcdn(&
     &                   point,0,0,0,2,jsp,sliceplot,stars,&
     &                   vacuum,sphhar,atoms,sym,cell,oneD,&
     &                   qpwCurrent,rhtxy,rhoCurrent,rht,xdnout)
                    IF (xsf) THEN
                       WRITE(55,*) xdnout
                    ELSE
                       WRITE(55,'(4e15.7)') point(:),xdnout
                    ENDIF
                 ENDDO xloop
              ENDDO
           ENDDO !z-loop
           IF (xsf.AND.jsp /= input%jspins) CALL xsf_WRITE_newblock(55,twodim&
     &          ,vec1,vec2,vec3,zero,grid)
        ENDDO !Spin-loop
        IF (xsf) THEN
           CALL xsf_WRITE_endblock(55,twodim)
        ELSE
           CLOSE(55)
        ENDIF
      ENDDO   !nplot
      CLOSE(18)
      IF (xsf) CLOSE(55)
      RETURN
      END SUBROUTINE plotVestPotDens0

      SUBROUTINE plotVestPotDens1(dimension,&
     &     stars,sphhar,atoms,&
     &     input,sym,cell, gdp, ngdp, qpt,&
     &      qpw, rho, plotname, iqpt, dens, plot_inpFileName)
!    *****************************************************
      USE m_outcdn
      USE m_loddop
      USE m_xsf_io
      use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn
      use m_types
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_dimension),INTENT(IN):: dimension
      TYPE(t_stars),INTENT(IN)    :: stars
      TYPE(t_sphhar),INTENT(IN)   :: sphhar
      TYPE(t_atoms),INTENT(IN)    :: atoms
      TYPE(t_input),INTENT(IN)    :: input
      TYPE(t_sym),INTENT(IN)      :: sym
      TYPE(t_cell),INTENT(IN)     :: cell
      complex, intent(in) :: qpw(:)
      complex, intent(in) :: rho(:, :, :)
      integer, intent(in) :: gdp(:, :)
      real, intent(in) :: qpt(:)
      integer, intent(in) :: ngdp
      integer, intent(in) :: iqpt
      logical, intent(in) :: dens

!     .. Local Scalars ..
      REAL          :: s,tec,qint,xdnout,fermiEnergyTemp
      INTEGER       :: i,i1,i2,i3,ii3,ix,iy,iz,jsp,na,nplo
      INTEGER       :: nplot,nq,nt,jm,jspin,iter,ii1,ii2
      LOGICAL       :: twodim,newform,l_qfix
!     ..
!     .. Local Arrays ..
      REAL    :: pt(3),vec1(3),vec2(3),vec3(3),zero(3)
      COMPLEX :: cdom(1),cdomvz(1,1),cdomvxy(1,1,1)
      INTEGER :: grid(3)
      LOGICAL :: cartesian,xsf
      REAL    :: rhocc(atoms%jmtd)
      REAL    :: point(3)
      CHARACTER (len=30) :: filename
      CHARACTER (len=7)  :: append
      character(len=*),  intent(in) ::plotname
      character(len=*),  intent(in) ::plot_inpFileName
      NAMELIST /plot/twodim,cartesian,vec1,vec2,vec3,grid,zero,filename


      INQUIRE(file =plot_inpFileName,exist= newform)
      IF (.NOT.newform) THEN !no input file exists, create a template and
                            !exit
         OPEN(20,file =plot_inpFileName)
         WRITE(20,'(i2,a5,l1)') 2,",xsf=",.true.
         WRITE(20,*) "&PLOT twodim=t,cartesian=t"
         WRITE(20,*) "  vec1(1)=10.0 vec2(2)=10.0"
         WRITE(20,*) "  filename='plot1' /"
         WRITE(20,*) "&PLOT twodim=f,cartesian=f"
         WRITE(20,*) "  vec1(1)=1.0 vec1(2)=0.0 vec1(3)=0.0 "
         WRITE(20,*) "  vec2(1)=0.0 vec2(2)=1.0 vec2(3)=0.0 "
         WRITE(20,*) "  vec3(1)=0.0 vec3(2)=0.0 vec3(3)=1.0 "
         WRITE(20,*) "  grid(1)=30  grid(2)=30  grid(3)=30  "
         WRITE(20,*) "  zero(1)=0.0 zero(2)=0.0 zero(3)=0.0 "
         WRITE(20,*) "  filename ='plot2' /"
         CLOSE(20)
         WRITE(*,*) "No " // plot_inpFileName // "file found. Created a template"
         CALL juDFT_error("Missing input for plot; modify" // plot_inpFileName&
     &        ,calledby ="plotdop")
      ENDIF

!      ! Initialize data
!      rho = 0.0
!      rho = rho
!      qpw = cmplx(0.0, 0.0)
!      qpw = qpw
!
!      !<--Perhaps only the core charge should be plotted
!      IF (input%score) THEN
!         OPEN (17,file='cdnc',form='unformatted',status='old')
!         REWIND 17
!         DO jspin = 1,input%jspins
!            DO nt = 1,atoms%ntype
!               jm = atoms%jri(nt)
!               READ (17) (rhocc(i),i=1,jm)
!               DO i = 1,atoms%jri(nt)
!                  rho(i,0,nt,jspin) = rho(i,0,nt,jspin) - rhocc(i)/2.0&
!     &                 /SQRT( pi_const )
!               ENDDO
!               READ (17) tec
!            ENDDO
!            READ (17) qint
!            qpw(1,jspin) = qpw(1,jspin) - qint/cell%volint
!         ENDDO
!         CLOSE (17)
!      END IF
      !>
      !>
      !<-- Open the plot_inp file for input
      OPEN (18,file=plot_inpFileName)
      READ(18,'(i2,5x,l1)')    nplot,xsf
      ! If xsf is specified we create an input file for xcrysden
      IF (xsf) THEN
         OPEN(55,file=plotname)
         CALL xsf_WRITE_atoms(&
     &                        55,atoms,input%film,.false.,cell%amat)
      ENDIF
      !<-- Loop over all plots
      DO nplo=1,nplot
         ! the defaults
         twodim = .TRUE.;cartesian=.TRUE.;grid=(/100,100,100/)
         vec1 = (/0.,0.,0./);vec2=(/0.,0.,0./);vec3=(/0.,0.,0./)
         zero = (/0.,0.,0./);filename="default"
         ! Read namelist plot defined above
         READ(18,plot)
         IF (twodim.AND.ANY(grid(1:2)<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (.NOT.twodim.AND.ANY(grid<1)) &
     &        CALL juDFT_error("Illegal grid size in plot",calledby&
     &        ="plotdop")
         IF (twodim) grid(3) = 1
         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(cell%amat,vec1)
            vec2=matmul(cell%amat,vec2)
            vec3=matmul(cell%amat,vec3)
            zero=matmul(cell%amat,zero)
         ENDIF
         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo
         IF (xsf) THEN
            CALL xsf_WRITE_header(55,twodim,filename,vec1,vec2,vec3,zero&
     &           ,grid)
         ELSE
           OPEN (55,file = filename,form='formatted')
         ENDIF
         !loop over spins
         DO jsp = 1,input%jspins
            !loop over all points
            DO iz = 0,grid(3)-1
               DO iy = 0,grid(2)-1
                  xloop:DO ix = 0,grid(1)-1
                    point = zero + vec1*REAL(ix)/(grid(1)-1) +&
     &                             vec2*REAL(iy)/(grid(2)-1)
                    IF (.NOT.twodim) point = point +&
     &                             vec3*REAL(iz)/(grid(3)-1)
                    !Check if the point is in MT-sphere
                    ii1 = 4
                    ii2 = 4
                    ii3 = 4
                    DO  i1 = -ii1,ii1
                       DO  i2 = -ii2,ii2
                          DO  i3 = -ii3,ii3
                             pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
                             na = 0
                             DO nt = 1,atoms%ntype
                                DO nq = 1,atoms%neq(nt)
                                   na   = na + 1
                                   ! Is the distance between the atoms and the point under consideration smaller or larger than R_MT
                                   s  = SQRT(dot_PRODUCT(atoms%pos(:,na)&
     &                                  -pt,atoms%pos(:,na)-pt))
                                   IF (s<atoms%rmsh(atoms%jri(nt),nt)) THEN
!                                      CALL outcdn(&
!     &                                     pt,nt,na,0,1,jsp,sliceplot,stars,&
!     &                                     vacuum,sphhar,atoms,sym,cell,oneD, &
!     &                                     qpw,rhtxy,rho,rht,xdnout)

                     if (dens) then
                     call evalPotDensAtPtlp1(pt, nt, na, 1, jsp, stars, sphhar, atoms, sym, cell, qpw, rho, ngdp, gdp, qpt, xdnout, iqpt,[-i1, -i2, -i3])
                   else
                     call evalPotDensAtPt(pt, nt, na, 1, jsp, stars, sphhar, atoms, sym, cell, qpw, rho, ngdp, gdp, qpt, xdnout, iqpt,[-i1, -i2, -i3])
                   end if
                                      IF (xsf) THEN
                                         write(55,*) xdnout
                                      ELSE
                                            WRITE(55,'(4e15.7)') point ,xdnout
                                      ENDIF
                                      CYCLE xloop
                                   ENDIF
                                ENDDO
                             ENDDO !nt
                          ENDDO
                       ENDDO
                    ENDDO !i1
                    if (dens) then
                    call evalPotDensAtPtlp1(point, 0, 0, 2, jsp, stars, sphhar, atoms, sym, cell, qpw, rho, ngdp, gdp, qpt, xdnout, iqpt, [0, 0, 0])
                  else
                    call evalPotDensAtPt(point, 0, 0, 2, jsp, stars, sphhar, atoms, sym, cell, qpw, rho, ngdp, gdp, qpt, xdnout, iqpt, [0, 0, 0])
                  end if
!                    CALL outcdn(&
!     &                   point,0,0,0,2,jsp,sliceplot,stars,&
!     &                   vacuum,sphhar,atoms,sym,cell,oneD,&
!     &                   qpw,rhtxy,rho,rht,xdnout)
                    IF (xsf) THEN
                       WRITE(55,*) xdnout
                    ELSE
                       WRITE(55,'(4e15.7)') point(:),xdnout
                    ENDIF
                 ENDDO xloop
              ENDDO
           ENDDO !z-loop
           IF (xsf.AND.jsp /= input%jspins) CALL xsf_WRITE_newblock(55,twodim&
     &          ,vec1,vec2,vec3,zero,grid)
        ENDDO !Spin-loop
        IF (xsf) THEN
           CALL xsf_WRITE_endblock(55,twodim)
        ELSE
           CLOSE(55)
        ENDIF
      ENDDO   !nplot
      CLOSE(18)
      IF (xsf) CLOSE(55)
      RETURN
    END SUBROUTINE plotVestPotDens1

    SUBROUTINE evalPotDensAtPt(&
     &                  p,n,na,iflag,jsp,stars,&
     &                  sphhar,atoms,sym,cell,&
     &                  qpw,rho,ngdp, gdp,qpt,&
     &                  xdnout, iqpt, latVec)
!
      USE m_ylm_old
      use m_types
      use m_constants

      use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn

      IMPLICIT NONE
!
!      TYPE(t_sliceplot),INTENT(IN) :: sliceplot
      TYPE(t_stars),INTENT(IN)     :: stars
   !   TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
!      TYPE(t_oneD),INTENT(IN)      :: oneD


!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iflag,jsp,n,na!,iv
      REAL,    INTENT (OUT) :: xdnout
      integer, intent(in) :: ngdp
      integer, intent(in) :: iqpt
!-odim
!+odim
!     ..
!     .. Array Arguments ..
      integer, intent(in) :: gdp(:, :)
      real, intent(in) :: qpt(:)
      integer, intent(in) :: latVec(3)
      COMPLEX, INTENT (IN) :: qpw(:) !(stars%ng3,dimension%jspd)
 !     COMPLEX, INTENT (IN) :: rhtxy(:,:,:,:) !(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd)
      complex,    INTENT (IN) :: rho(:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
 !     REAL,    INTENT (IN) :: rht(:,:,:) !(vacuum%nmzd,2,dimension%jspd)
      REAL,    INTENT (INOUT) :: p(3)
!     ..
!     .. Local Scalars ..
      REAL delta,s,sx,xx1,xx2,rrr,phi
      INTEGER i,j,jp3,jr,k,lh,mem,nd,nopa,ivac,ll1,lm ,gzi,m
      COMPLEX ci
      complex xd1, xd2
!     ..
!     .. Local Arrays ..
      COMPLEX sf2(stars%ng2),sf3(stars%ng3),ylm((atoms%lmaxd+1)**2)
      REAL rcc(3),x(3)
      complex, allocatable :: pwb(:) !plane -wave basis
      complex :: xdnoutTest
      integer :: iG, oqn_l, mqn_m


!     ..
!     ..
      ci = cmplx(0.,1.)
!      ivac=iv


      if (iflag.ne.1) then
      if (iflag.ne.0) THEN
!     ---> interstitial part
      !CALL cotra1(p(1),rcc,cell%bmat)
      rcc=matmul(cell%bmat,p)/tpi_const
      allocate(pwb(ngdp))
      pwb(:) = 0
      do iG = 1, ngdp
        pwb(iG) = exp( tpi_const *cmplx(0., 1.) * dot_product(gdp(:, iG) + qpt(:), rcc(:)) )
      end do

      xdnoutTest = cmplx(0.0, 0.0)
      ! if we use dot_product we have to conjugate the first entry!
      if (iqpt == 1) then
        xdnoutTest = dot_product(conjg(qpw(:ngdp)), pwb(:ngdp))
      else
        xdnoutTest = 2 * real(dot_product(conjg(qpw(:ngdp)), pwb(:ngdp)))
      end if
      !xdnoutTest =sum(pwb(:))

      if ( aimag(xdnoutTest) > 1e-8 ) then
        call juDFT_error( "Complex point in density or potential within interstitial region.", calledby ="evalPotDensAtPt", hint='Fix bug.')
        write(*, *) 'xdnout =', xdnoutTest
        NOstopNO
      end if

      xdnout = real(xdnoutTest)

      RETURN
      ENDIF
      endif
     ! xdnout = 0.
!-odim
!     ----> m.t. part

      !nd = atoms%ntypsy(na)
      !nopa = atoms%ngopr(na)
      !IF (oneD%odi%d1) nopa = oneD%ods%ngopr(na)
      sx = 0.0
      DO  i = 1,3
         x(i) = p(i) - atoms%pos(i,na)
         sx = sx + x(i)*x(i)
      enddo
      sx = sqrt(sx)
     ! IF (nopa.NE.1) THEN
!... !switch to internal units
     !    !CALL cotra1(x,rcc,cell%bmat)
     !    rcc=matmul(cell%bmat,x)/tpi_const
!... !rotate into representative
     !    DO  i = 1,3
     !       p(i) = 0.
     !       DO  j = 1,3
     !         IF (.NOT.oneD%odi%d1) THEN
     !          p(i) = p(i) + sym%mrot(i,j,nopa)*rcc(j)
     !         ELSE
     !          p(i) = p(i) + oneD%ods%mrot(i,j,nopa)*rcc(j)
     !         END IF
     !    enddo
     !    enddo
!... !switch back to cartesian units
     !    !CALL cotra0(p,x,cell%amat)
     !    x=matmul(cell%amat,p)
     ! END IF
      DO j = atoms%jri(n),2,-1
         IF (sx.GE.atoms%rmsh(j,n)) EXIT
      ENDDO
      jr = j
   !   call ylmnorm_init(atoms%lmax(n) + 1)
      !call ylm4(atoms%lmax(n) + 1, x, ylm)
      call ylm4(atoms%lmax(n), x, ylm)
   !   call ylmnorm_init(atoms%lmax(n))
!      CALL ylm4(&
!     &          atoms%lmax(n),x,&
!     &          ylm)
      xd1 = cmplx(0.0, 0.0)
      xd2 = cmplx(0.0, 0.0)
      !DO  lh = 0, sphhar%nlh(nd)
      do oqn_l = 0, atoms%lmax(n)! + 1
        ll1 = oqn_l * ( oqn_l + 1 ) + 1
        !s = 0.0
        !DO mem = 1,sphhar%nmem(lh,nd)
        DO mqn_m = - oqn_l, oqn_l
          lm = ll1 + mqn_m
      !    s = s + real( sphhar%clnu(mem,lh,nd)*ylm(lm) )
      !  ENDDO
          !IF (sliceplot%plpot) THEN
            if (iqpt == 1) then
              xd1 = xd1 + rho(jr,lm,na) * ylm(lm)
            else
              xd1 = xd1 + 2 * real(rho(jr,lm,na) * ylm(lm) * exp(ci * tpi_const * dot_product(qpt(:), latVec(:))))
            end if

          !ELSE
          !  xd1 = xd1 + rho(jr,lm,n,jsp) * ylm(lm) / (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
          !END IF
          IF (jr.EQ.atoms%jri(n)) CYCLE
          !IF (sliceplot%plpot) THEN
          if (iqpt == 1) then
             xd2 = xd2 + rho(jr+1,lm,na)*ylm(lm)
          else
             xd2 = xd2 + 2 * real(rho(jr+1,lm,na)*ylm(lm) * exp(ci * tpi_const * dot_product(qpt(:), latVec(:))))
          end if

          !ELSE
          !   xd2 = xd2 + rho(jr+1,lm,n,jsp)*ylm(lm)/ (atoms%rmsh(jr+1,n)*atoms%rmsh(jr+1,n))
          !END IF
        ENDDO
      end do
      xdnoutTest = cmplx(0.0, 0.0)
      IF (jr.EQ.atoms%jri(n)) THEN
         xdnoutTest = xd1
      ELSE
         xdnoutTest = xd1 + (xd2-xd1) *&
     &                  (sx-atoms%rmsh(jr,n)) / (atoms%rmsh(jr+1,n)-atoms%rmsh(jr,n))
      END IF
      if ( aimag(xdnoutTest) > 1e-8 ) then
        !write(1003, '(2(f15.8))') xdnoutTest
        call juDFT_error( "Complex point in density or potential within muffin-tin region.", calledby ="evalPotDensAtPt", hint='Fix bug.')
      end if

      xdnout = real(xdnoutTest)
!
      RETURN
    END SUBROUTINE evalPotDensAtPt

    SUBROUTINE evalPotDensAtPtlp1(&
     &                  p,n,na,iflag,jsp,stars,&
     &                  sphhar,atoms,sym,cell,&
     &                  qpw,rho,ngdp, gdp,qpt,&
     &                  xdnout, iqpt, latVec)
!
      USE m_ylm_old
      use m_types
      use m_constants

      use m_juDFT_NOstopNO, only : juDFT_error, juDFT_warn

      IMPLICIT NONE
!
!      TYPE(t_sliceplot),INTENT(IN) :: sliceplot
      TYPE(t_stars),INTENT(IN)     :: stars
   !   TYPE(t_vacuum),INTENT(IN)    :: vacuum
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
!      TYPE(t_oneD),INTENT(IN)      :: oneD


!     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: iflag,jsp,n,na!,iv
      REAL,    INTENT (OUT) :: xdnout
      integer, intent(in) :: ngdp
      integer, intent(in) :: iqpt
!-odim
!+odim
!     ..
!     .. Array Arguments ..
      integer, intent(in) :: gdp(:, :)
      real, intent(in) :: qpt(:)
      integer, intent(in) :: latVec(3)
      COMPLEX, INTENT (IN) :: qpw(:) !(stars%ng3,dimension%jspd)
 !     COMPLEX, INTENT (IN) :: rhtxy(:,:,:,:) !(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd)
      complex,    INTENT (IN) :: rho(:,:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd)
 !     REAL,    INTENT (IN) :: rht(:,:,:) !(vacuum%nmzd,2,dimension%jspd)
      REAL,    INTENT (INOUT) :: p(3)
!     ..
!     .. Local Scalars ..
      REAL delta,s,sx,xx1,xx2,rrr,phi
      INTEGER i,j,jp3,jr,k,lh,mem,nd,nopa,ivac,ll1,lm ,gzi,m
      COMPLEX ci
      complex xd1, xd2
!     ..
!     .. Local Arrays ..
      COMPLEX sf2(stars%ng2),sf3(stars%ng3),ylm((atoms%lmaxd+1)**2)
      REAL rcc(3),x(3)
      complex, allocatable :: pwb(:) !plane -wave basis
      complex :: xdnoutTest
      integer :: iG, oqn_l, mqn_m


!     ..
!     ..
      ci = cmplx(0.,1.)
!      ivac=iv


      if (iflag.ne.1) then
      if (iflag.ne.0) THEN
!     ---> interstitial part
      !CALL cotra1(p(1),rcc,cell%bmat)
      rcc=matmul(cell%bmat,p)/tpi_const
      allocate(pwb(ngdp))
      pwb(:) = 0
      do iG = 1, ngdp
        pwb(iG) = exp( tpi_const *cmplx(0., 1.) * dot_product(gdp(:, iG) + qpt(:), rcc(:)) )
      end do

      xdnoutTest = cmplx(0.0, 0.0)
      ! if we use dot_product we have to conjugate the first entry!
      if (iqpt == 1) then
        xdnoutTest = dot_product(conjg(qpw(:ngdp)), pwb(:ngdp))
      else
        xdnoutTest = 2 * real(dot_product(conjg(qpw(:ngdp)), pwb(:ngdp)))
      end if
      !xdnoutTest =sum(pwb(:))

      if ( aimag(xdnoutTest) > 1e-8 ) then
        call juDFT_error( "Complex point in density or potential within interstitial region.", calledby ="evalPotDensAtPt", hint='Fix bug.')
        write(*, *) 'xdnout =', xdnoutTest
        NOstopNO
      end if

      xdnout = real(xdnoutTest)

      RETURN
      ENDIF
      endif
     ! xdnout = 0.
!-odim
!     ----> m.t. part

      !nd = atoms%ntypsy(na)
      !nopa = atoms%ngopr(na)
      !IF (oneD%odi%d1) nopa = oneD%ods%ngopr(na)
      sx = 0.0
      DO  i = 1,3
         x(i) = p(i) - atoms%pos(i,na)
         sx = sx + x(i)*x(i)
      enddo
      sx = sqrt(sx)
     ! IF (nopa.NE.1) THEN
!... !switch to internal units
     !    !CALL cotra1(x,rcc,cell%bmat)
     !    rcc=matmul(cell%bmat,x)/tpi_const
!... !rotate into representative
     !    DO  i = 1,3
     !       p(i) = 0.
     !       DO  j = 1,3
     !         IF (.NOT.oneD%odi%d1) THEN
     !          p(i) = p(i) + sym%mrot(i,j,nopa)*rcc(j)
     !         ELSE
     !          p(i) = p(i) + oneD%ods%mrot(i,j,nopa)*rcc(j)
     !         END IF
     !    enddo
     !    enddo
!... !switch back to cartesian units
     !    !CALL cotra0(p,x,cell%amat)
     !    x=matmul(cell%amat,p)
     ! END IF
      DO j = atoms%jri(n),2,-1
         IF (sx.GE.atoms%rmsh(j,n)) EXIT
      ENDDO
      jr = j
      call ylmnorm_init(atoms%lmax(n) + 1)
      !call ylm4(atoms%lmax(n) + 1, x, ylm)
      call ylm4(atoms%lmax(n) + 1, x, ylm)
      call ylmnorm_init(atoms%lmax(n))
!      CALL ylm4(&
!     &          atoms%lmax(n),x,&
!     &          ylm)
      xd1 = cmplx(0.0, 0.0)
      xd2 = cmplx(0.0, 0.0)
      !DO  lh = 0, sphhar%nlh(nd)
      do oqn_l = 0, atoms%lmax(n) + 1
        ll1 = oqn_l * ( oqn_l + 1 ) + 1
        !s = 0.0
        !DO mem = 1,sphhar%nmem(lh,nd)
        DO mqn_m = - oqn_l, oqn_l
          lm = ll1 + mqn_m
      !    s = s + real( sphhar%clnu(mem,lh,nd)*ylm(lm) )
      !  ENDDO
          !IF (sliceplot%plpot) THEN
            if (iqpt == 1) then
              xd1 = xd1 + rho(jr,lm,na) * ylm(lm)
            else
              xd1 = xd1 + 2 * real(rho(jr,lm,na) * ylm(lm) * exp(ci * tpi_const * dot_product(qpt(:), latVec(:))))
            end if

          !ELSE
          !  xd1 = xd1 + rho(jr,lm,n,jsp) * ylm(lm) / (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
          !END IF
          IF (jr.EQ.atoms%jri(n)) CYCLE
          !IF (sliceplot%plpot) THEN
          if (iqpt == 1) then
             xd2 = xd2 + rho(jr+1,lm,na)*ylm(lm)
          else
             xd2 = xd2 + 2 * real(rho(jr+1,lm,na)*ylm(lm) * exp(ci * tpi_const * dot_product(qpt(:), latVec(:))))
          end if

          !ELSE
          !   xd2 = xd2 + rho(jr+1,lm,n,jsp)*ylm(lm)/ (atoms%rmsh(jr+1,n)*atoms%rmsh(jr+1,n))
          !END IF
        ENDDO
      end do
      xdnoutTest = cmplx(0.0, 0.0)
      IF (jr.EQ.atoms%jri(n)) THEN
         xdnoutTest = xd1
      ELSE
         xdnoutTest = xd1 + (xd2-xd1) *&
     &                  (sx-atoms%rmsh(jr,n)) / (atoms%rmsh(jr+1,n)-atoms%rmsh(jr,n))
      END IF
      if ( aimag(xdnoutTest) > 1e-8 ) then
        !write(1003, '(2(f15.8))') xdnoutTest
        call juDFT_error( "Complex point in density or potential within muffin-tin region.", calledby ="evalPotDensAtPtlp1", hint='Fix bug.')
      end if

      xdnout = real(xdnoutTest)
!
      RETURN
    END SUBROUTINE evalPotDensAtPtlp1

  !Plot potential routine extracted from checkGradPot0s
  subroutine plotPotPathUC(cellT, atomsT, vGradCoul0IR, ngdp, gdp,  gradVrjuPhon, paPoXCo, paPoYCo, paPoZCo, qpoint)

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_juDFT_NOstopNO
    use mod_juPhonUtils, only : fopen, fclose

    implicit none

!    type(t_sym),        intent(in) :: symT
    type(t_cell),       intent(in) :: cellT
!    type(t_sphhar),     intent(in) :: sphharT
!    type(t_stars),      intent(in) :: starsT
    type(t_atoms),      intent(in) :: atomsT
!    type(t_input),      intent(in) :: inputT
    complex,            intent(in) :: vGradCoul0IR(:, :)
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: gdp(:, :)
!    integer,            intent(in) :: mlh_atom(:, 0:, :)
!    integer,            intent(in) :: nmem_atom(0:, :)
!    complex,            intent(in) :: clnu_atom(:, 0:, :)
!    real,               intent(in) :: rho0MT(:, 0:, :, :)
    real,               intent(in) :: paPoXCo
    real,               intent(in) :: paPoYCo
    real,               intent(in) :: paPoZCo
  complex,             intent(in) :: gradVrjuPhon(:, :, :, :)
  real, intent(in) :: qpoint(:)
!    logical,            intent(in) :: harSw
!    logical,            intent(in) :: extSw
!    logical,            intent(in) :: xcSw
!
    integer                        :: iGvec
!    real                           :: Gext(3)
!    complex                        :: vCoulBenchmark(ngdp, 3)
!    complex                        :: gradVcMT(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%nat, 3)
!    real                           :: sqr4pi3
!    real                           :: Vc0nSym(atomsT%jmtd, (atomsT%lmaxd + 1 +  1)**2, atomsT%nat) ! todo is this the correct dimension
    integer                        :: iatom
    integer                        :: itype
!    integer                        :: ineq
!    integer                        :: symType
!    integer                        :: lh
    integer                        :: oqn_l
!    integer                        :: lm_temp
!    integer                        :: mem
!    integer                        :: mems
    integer                        :: mqn_m
    integer                        :: mqn_mpp
    integer                        :: lm
!    integer                        :: lmMinus
!    integer                        :: lmPlus
!    real                           :: Vc0MTDerived(atomsT%jmtd)
!    real                           :: tempGaunt1, tempGaunt2
    integer                        :: imesh
!    complex                        :: vpw_G(ngdp)
!    complex                        :: vpw_G_coul(ngdp)
!    complex                        :: vpw_G_hart(ngdp)
!    complex                        :: vpw_G_xc(ngdp)
!    complex              :: vpwCoulUwStars1(starsT%n3d, 1)
!    complex              :: vpwCoulUwStars2(starsT%n3d, 1)
!    complex              :: vpwStar(starsT%n3d, 1)
!    complex              :: vpwStarDiff(starsT%n3d, 1)
!    complex                        :: vpwRfleur(3)
    complex                        :: vpwRjuPhon(3)
!    complex                        :: exponentialjuPhon
!    complex                        :: exponentialFLEUR
    complex                        :: exponential
    integer           :: isIR
    integer           :: wasIR
    integer          :: c1, c2, c3
    real                           :: ucpath(3)
    real                              :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                           :: ucpathExt(3)
    real                      :: ucpathExtc(3), ucpathExta(3)
    real                           :: dx
    integer                        :: ii
    real                           :: x
!    complex                        :: newvpwFLEUR(starsT%n3d, 1) ! make this declaration more general!
!    complex                        :: newvpw_G(ngdp)
!    complex                        :: onlyPotFLEUR
!    complex                        :: onlypsqFLEUR
    integer                        :: ieq
    integer                        :: ext2brav(3)
    integer                        :: ext2bravLock(3)
    integer                        :: atomLock
    integer                        :: itypeLock
!    logical                        :: pseudoPot = .true.
    logical                        :: realPot = .true.
    integer                   :: iterations
!  real   :: vrFleur(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
!  real   :: vrFleurDiff(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
!  real   :: vrFleurxc(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
!  real   :: vrTemp(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
!  complex,      allocatable :: gradVrFleur(:, :, :, :)
  real    :: direc(3)
  real :: direcExt(3)
  complex                                     :: ylm((atomsT%lmaxd + 1)**2)
!  complex                                     :: vMTfleur(3, atomsT%jmtd)
  complex                                     :: vMTjuPhon(3, atomsT%jmtd)
  integer                          :: idirec
!  real                              :: av, dmx, rms
!  real                              :: IRvalue(1), MTvalue(1)
!  integer                             :: jspin
!  integer                               :: ilh
!  real        :: potFleurMT(atomsT%jmtd)
!  integer    :: symAt_temp, symOpr_temp
!  real      :: rvec(3), rvec_int(3)
!  integer    :: imatmulr, imatmulc
!  real :: rotatedVector(3)
!  integer :: ilath
!  real :: linCombBas
!  integer :: l_temp
!  integer :: imem
!
!logical :: fBenchSw

    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
  !  fBenchSw = .true.
!    vpwStar = 0
  !  vpwStarDiff = 0
  !  vrFLEUR = 0
  !  vrFLEURDiff = 0

!    !todo do this with pottot and potcoul
!    if ( harSw .and. .not.extSw .and. .not.xcSw ) then ! Hartree potential
!      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( .not.harSw .and. extSw .and. .not.xcSw ) then ! external potential
!      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen( 1000, name='vpw_hart', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStarDiff
!      call fclose( 1000 )
!
!      vpwStar = vpwStar - vpwStarDiff
!
!
!      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!      call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
!        read(1000) vrFLEURDiff
!      call fclose(1000)
!
!      vrFLEUR = vrFLEUR - vrFLEURDiff
!
!    else if ( .not.harSw .and. .not.extSw .and. xcSw ) then ! exchange-correlation potential
!      call fopen( 1000, name='vpw_xc', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( harSw .and. extSw .and. .not.xcSw ) then ! Coulomb potential
!      call fopen( 1000, name='vpw_coul', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!
!      call fopen(1000, name='v0MTFLEUR_coul', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else if ( harSw .and. extSw .and. xcSw ) then ! effective potential
!      call fopen( 1000, name='vpw_eff', status='old', action='read', form='unformatted' )
!        read( 1000 ) vpwStar
!      call fclose( 1000 )
!      call fopen(1000, name='v0MTFLEUR_eff', status='old', action='read', form='unformatted')
!        read(1000) vrFLEUR
!      call fclose(1000)
!
!    else
!      call juDFT_warn('FLEUR benchmark for this potential combination not available. Setting benchmark curve equals zero.', &
!        & calledby='checkGradPot0s', hint='Choose Hartree, external, exchange-correlation, Coulomb or effective potential', &
!        & file='jpTestPotential_mod.f90')
!      fBenchSw = .false.
!
!    end if

!    call ConvertStar2G( vpwStar(:, 1), vpw_G, starsT, ngdp, gdp )

!    do iGvec = 1, ngdp
!      Gext(:) = matmul(cellT%bmat, gdp(:, iGvec))
!      vCoulBenchmark(iGvec, :) = iu * Gext(:) * vpw_G(iGvec)
!    enddo
!    call calcGrFinLH(atomsT, sphharT, clnu_atom, nmem_atom, mlh_atom, vrFleur(:, :, :, 1), gradVrFleur)


!    call fopen(1111, name='pathPseudoPot', status='replace', action='write', form='formatted')
!    call fopen(1222, name='potFLEUR', status='replace', action='write', form='formatted')
    call fopen(1333, name='pathPot', status='replace', action='write', form='formatted')
!    call fopen(1444, name='pathPotCont', status='replace', action='write', form='formatted')


    dx = 1. / 500.
    x = 0 - dx
    wasIR = -1
    do ii = 0, 500
      x = x + dx
      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
      iatom = 1
      isIR = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          do c1 = -1, 1
            do c2 = -1, 1
              do c3 = -1, 1
                ext2brav = [c1, c2, c3]
                ucpathExt= matmul(cellT%amat, ucpath - atomsT%taual(:, iatom) - ext2brav)
                if ( norm2 ( ucpathExt ) <= atomsT%rmt(itype) ) then
                  isIR = iatom ! todo atomLock is redundant!
                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                  ext2bravLock = ext2brav
                  itypeLock = itype
                end if
              end do
            end do
          end do
          iatom = iatom + 1
        end do
      end do

!      if ( pseudoPot ) then
!
!        vpwRfleur = cmplx(0.,0.)
!        vpwRjuPhon = cmplx(0.,0.)
!
!        do iGvec = 1, ngdp
!          exponential = exp(iu * tpi * dot_product(gdp(:, iGvec), ucpath))
!          vpwRfleur(:)  = vpwRfleur(:)  + vCoulBenchmark(iGvec, :) * exponential
!          vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :)   * exponential
!        end do
!
!        write (1111, '(13(es16.8E3, 2x),i2, 3(2x, i2))') x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), &
!          & aimag(vpwRfleur(2)), real(vpwRfleur(3)), aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)),     &
!          & real(vpwRjuPhon(2)), aimag(vpwRjuPhon(2)), real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), &
!          & ext2bravLock(2), ext2bravLock(3)
!
!      end if

      if ( realPot ) then ! if no pseudoPotMode
        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

!            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
!            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpoint, ucpath))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
!              onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponentialFLEUR ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

!            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            !onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpoint, ucpath))
!              vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              !write (*, *) iGvec, gdp(:, iGvec), vpwRjuPhon
              !write (*, *)
              !read(*, *)
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do

            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR > 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000 !todo does this loop have to go to 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atomsT%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            !vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            !onlyPotFLEUR = cmplx(0., 0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpoint, ucpathC))
              !vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
!            call ylmnorm_init(atomsT%lmaxd + 1)
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
!            call ylmnorm_init(atomsT%lmaxd)

            !vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock)! + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
            !        vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
            vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm) * exp(tpi * iu * dot_product(qpoint, ext2bravLock))
                  end do
                end do
              end do
            end do

            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), 0, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            do imesh = atomsT%jri(itypeLock), 1, -1
              write(1333, '(7(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atomsT%rmsh(imesh, itypeLock), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

           ! do idirec = 1, 3
           !   if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
           !     write(1444, '(8(es16.8E3, 2x))') &
           !       &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
           !       &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
           !   else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
           !     write(1444, '(8(es16.8E3, 2x))') &
           !       &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
           !       &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
           !   else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
           !     write(1444, '(8(es16.8E3, 2x))') &
           !       &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
           !       &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
           !   else
           !     write(1444, '(8(es16.8E3, 2x))') &
           !       &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
           !       &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
           !       &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
           !       &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
           !   end if
           ! end do

          else if (wasIR > 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atomsT%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            direcExt = direcExt / norm2(direcExt)
       !     call ylmnorm_init(atomsT%lmaxd + 1)
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
       !     call ylmnorm_init(atomsT%lmaxd)

            !vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atomsT%lmax(itypeLock) !+ 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atomsT%jri(itypeLock)
            !        vMTfLEUR(idirec, imesh) = vMTFLEUR(idirec, imesh) + gradVrFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, idirec, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end do
            do imesh = 1, atomsT%jri(itypeLock)
              write(1333, '(7(es20.10E3, 2x),i2, 3(2x, i2))')&
                &atomsT%rmsh(imesh, itypeLock), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                &aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            !todo seems to be the real ptoential here
          !  iatom = 0
           ! potFleurMT = 0
           !   !Evaluate MT sphere function at same random points
           ! symAt_temp  = atomsT%ntypsy(atomLock) ! symmetry of atom iatom
           ! symOpr_temp = atomsT%ngopr (atomLock) ! symmetry operation mapping local to global coordinate system
           ! do imesh = 1, atomsT%jri(itypeLock)
           !   rvec = direcExt
           !   if (symOpr_temp /= 1) then
           !   !Transform into internal real space coordinates
           !   !todo might be wrong bmat
           !     call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
           !     do imatmulr = 1, 3
           !       rotatedVector(imatmulr) = 0.0
           !       do imatmulc = 1, 3
           !         rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) * rvec_int(imatmulc)
           !       end do
           !     end do
           !     call cotra0(rotatedVector, rvec, cellT%amat)
           !   end if
           !   call ylm4(atomsT%lmax(itypeLock), rvec, ylm)

           !   do ilath = 0, sphharT%nlh(symAt_temp)
           !     linCombBas = 0.0
           !     l_temp = sphharT%llh(ilath, symAt_temp) * (sphharT%llh(ilath, symAt_temp) + 1) + 1 !understand this
           !     do imem = 1, sphharT%nmem(ilath, symAt_temp)
           !       lm_temp = l_temp + sphharT%mlh(imem, ilath, symAt_temp)
           !       linCombBas = linCombBas + real(clnu_atom(imem, ilath, symAt_temp) * ylm(lm_temp))
           !     end do
           !       potFleurMT(imesh) = potFleurMT(imesh) + vrFleur(imesh, ilath, itypeLock, 1) * linCombBas
           !   end do
           ! end do
!            call fopen(1020, name='potMTfleur', status='replace', action='write', form='formatted')
!            do imesh = 1, atomsT%jri(itypeLock)
!              write(1020, '(3(es16.8E3, 2x))') atomsT%rmsh(imesh, itypeLock), potFleurMT(imesh), real(vMTfLEUR(1, imesh))
!            end do
!            call fclose(1020)


          !  vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            !onlyPotFLEUR = cmplx(0., 0.)
            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpoint, ucpathC))
          !    vpwRfleur(:) = vpwRfleur(:) + vCoulBenchmark(iGvec, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
!            write (1222, '(13(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR) ! todo the lines down really needed?
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point ! todo this is proably the wrong unit number, and was changed from 1009 to 1333
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), 0, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

!            do idirec = 1, 3
!              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
!              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              else
!                write(1444, '(8(es16.8E3, 2x))') &
!                  &real(vMTfleur(idirec, atomsT%jri(itypeLock))), real(vpwRfleur(idirec)), &
!                  &abs(real(vMTfleur(idirec, atomsT%jri(itypeLock))) - real(vpwRfleur(idirec))), -9e9,&
!                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
!                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
!              end if
!            end do

          end if! if still in MT then wait for crossing

        end if !is not first point on path
      end if ! is not in pseudopot mode
      wasIR = isIR
    end do
!    call fclose(1111)
!    call fclose(1222)
    call fclose(1333)
!    call fclose(1444)
  end subroutine plotPotPathUC

  ! we have to read in from disk and plot what we want to have
  subroutine plot1stDensVar(symT, cellT, inputT, sphharT, atomsT, vGradCoul0IR, ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, &
      &rho0MT, qpt, gradVrjuPhon, paPoXCo, paPoYCo, paPoZCo, iDatom, counter, iqpt )

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use mod_juPhonUtils, only : fopen, fclose
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_jpPotDensHelper, only : convertStar2G
    use m_juDFT_NOstopNO

    implicit none

    type(t_sym),        intent(in) :: symT
    type(t_cell),       intent(in) :: cellT
    type(t_sphhar),     intent(in) :: sphharT
    type(t_atoms),      intent(in) :: atomsT
    type(t_input),      intent(in) :: inputT
    complex,            intent(in) :: vGradCoul0IR(:, :)
    integer,            intent(in) :: ngdp
    integer,            intent(in) :: gdp(:, :)
    integer,            intent(in) :: mlh_atom(:, 0:, :)
    integer,            intent(in) :: nmem_atom(0:, :)
    complex,            intent(in) :: clnu_atom(:, 0:, :)
    real,               intent(in) :: rho0MT(:, 0:, :, :)
    real,               intent(in) :: qpt(:)
    real,               intent(in) :: paPoXCo
    real,               intent(in) :: paPoYCo
    real,               intent(in) :: paPoZCo
    integer,            intent(in) :: iDatom
    integer,            intent(in) :: counter
    integer,            intent(in) :: iqpt

    integer                        :: iGvec
    real                           :: Gext(3)
    complex                        :: vCoulBenchmark(ngdp, 3)
    complex                        :: gradVcMT(atomsT%jmtd, (atomsT%lmaxd + 1)**2, atomsT%nat, 3)
    real                           :: sqr4pi3
    real                           :: Vc0nSym(atomsT%jmtd, (atomsT%lmaxd + 1 )**2, atomsT%nat) ! todo is this the correct dimension
    integer                        :: iatom
    integer                        :: itype
    integer                        :: ineq
    integer                        :: symType
    integer                        :: lh
    integer                        :: oqn_l
    integer                        :: lm_temp
    integer                        :: mem
    integer                        :: mems
    integer                        :: mqn_m
    integer                        :: mqn_mpp
    integer                        :: lm
    integer                        :: lmMinus
    integer                        :: lmPlus
    real                           :: Vc0MTDerived(atomsT%jmtd)
    real                           :: tempGaunt1, tempGaunt2
    integer                        :: imesh
    complex                        :: vpw_G(ngdp)
    complex                        :: vpw_G_coul(ngdp)
    complex                        :: vpw_G_hart(ngdp)
    complex                        :: vpw_G_xc(ngdp)
    complex                        :: vpwRfleur(3)
    complex                        :: vpwRjuPhon(3)
    complex                        :: exponentialjuPhon
    complex                        :: exponentialFLEUR
    complex                        :: exponential
    integer           :: isIR
    integer           :: wasIR
    integer          :: c1, c2, c3
    real                           :: ucpath(3)
    real                              :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                           :: ucpathExt(3)
    real                      :: ucpathExtc(3), ucpathExta(3)
    real                           :: dx
    integer                        :: ii
    real                           :: x
    complex                        :: newvpw_G(ngdp)
    complex                        :: onlyPotFLEUR
    complex                        :: onlypsqFLEUR
    integer                        :: ieq
    integer                        :: ext2brav(3)
    integer                        :: ext2bravLock(3)
    integer                        :: atomLock
    integer                        :: itypeLock
    logical                        :: pseudoPot = .true.
    logical                        :: realPot = .true.
    integer                   :: iterations
  real   :: vrFleur(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurDiff(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrFleurxc(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  real   :: vrTemp(atomsT%jmtd, 0:sphharT%nlhd, atomsT%ntypd, inputT%jspins)
  complex,             intent(in) :: gradVrjuPhon(:, :, :, :)
  complex,      allocatable :: gradVrFleur(:, :, :, :)
  real    :: direc(3)
  real :: direcExt(3)
  complex                                     :: ylm((atomsT%lmaxd + 1)**2)
  complex                                     :: vMTfleur(3, atomsT%jmtd)
  complex                                     :: vMTjuPhon(3, atomsT%jmtd)
  integer                          :: idirec
  real                              :: av, dmx, rms
  real                              :: IRvalue(1), MTvalue(1)
  integer                             :: jspin
  integer                               :: ilh
  real        :: potFleurMT(atomsT%jmtd)
  integer    :: symAt_temp, symOpr_temp
  real      :: rvec(3), rvec_int(3)
  integer    :: imatmulr, imatmulc
  real :: rotatedVector(3)
  integer :: ilath
  real :: linCombBas
  integer :: l_temp
  integer :: imem
  character(len=11) :: filename

logical :: fBenchSw

!TODO Fix bug that the interstitial part is recognized as MT part because isIR should be 1 if additionally IR point is added at MT surface
! after bisection method
    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
    fBenchSw = .true.
    vrFLEUR = 0
    vrFLEURDiff = 0


    !call fopen(1222, name='densFLEUR', status='replace', action='write', form='formatted')
    write(filename, '(a,i1,i2)') 'pathDens', iDatom, counter
    call fopen(1333, name=filename, status='replace', action='write', form='formatted')
    write(filename, '(a,i1,i1)') 'pathDenC', iDatom,counter
    call fopen(1444, name=filename, status='replace', action='write', form='formatted')

    dx = 1. / 600.
    x = 0 - dx
    wasIR = -1
    do ii = 0, 600
      x = x + dx
      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
      iatom = 1
      isIR = 0
      do itype = 1, atomsT%ntype
        do ieq = 1, atomsT%neq(itype)
          do c1 = -1, 1
            do c2 = -1, 1
              do c3 = -1, 1
                ext2brav = [c1, c2, c3]
                ucpathExt= matmul(cellT%amat, ucpath - atomsT%taual(:, iatom) - ext2brav)
                if ( norm2 ( ucpathExt ) < atomsT%rmt(itype) ) then
                  isIR = iatom ! todo atomLock is redundant!
                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                  ext2bravLock = ext2brav
                  itypeLock = itype
                end if
              end do
            end do
          end do
          iatom = iatom + 1
        end do
      end do

      if ( realPot ) then ! if no pseudoPotMode
        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product((gdp(:, iGvec) + qpt(:)), ucpath))
              if (iqpt /= 1) then
                vpwRjuPhon(:) = vpwRjuPhon(:) + 2 * real(vGradCoul0IR(iGvec, :) * exponential)! todo has been changed, seems to be added so that starting point ist correct!
              else
                vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
              end if
              !todo floating is inexact here for some reasons....
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponentialFLEUR ! pure potential
            end do

            !write (1222, '(3(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0.,0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpt(:), ucpath))
              if (iqpt /= 1) then
                vpwRjuPhon(:) = vpwRjuPhon(:) + 2 * real(vGradCoul0IR(iGvec, :) * exponential)
              else
                vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              end if
              !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do

            !write (1222, '(3(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR > 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000 !todo does this loop have to go to 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atomsT%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)

            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpt(:), ucpathC))
              if (iqpt /= 1) then
                vpwRjuPhon(:) = vpwRjuPhon(:) + 2 * real(vGradCoul0IR(iGvec, :) * exponential)
              else
                vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              end if
            !onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            !write (1222, '(3(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR)
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            !direcExt = direcExt / norm2(direcExt)
            !call ylmnorm_init(atomsT%lmaxd + 1)
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
            !call ylmnorm_init(atomsT%lmaxd)

            vMTjuPhon = 0
            if (iqpt /= 1) then
              do idirec = 1, 3
                do oqn_l = 0, atomsT%lmax(itypeLock) !+ 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atomsT%jri(itypeLock)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + 2 * real(gradVrjuPhon(imesh, lm, atomLock, idirec) * ylm(lm) * exp(iu * tpi * dot_product(qpt(:), ext2bravLock(:))))
                    end do
                  end do
                end do
              end do
            else
              do idirec = 1, 3
                do oqn_l = 0, atomsT%lmax(itypeLock) !+ 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atomsT%jri(itypeLock)
                    vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, atomLock, idirec) * ylm(lm)
                    end do
                  end do
                end do
              end do
            end if

            do imesh = atomsT%jri(itypeLock), 1, -1
              write(1333, '(7(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atomsT%rmsh(imesh, itypeLock), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            do idirec = 1, 3
              if ( ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(4(es16.8E3, 2x))') &
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else
                write(1444, '(4(es16.8E3, 2x))') &
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          else if (wasIR > 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cellT%amat, ucpathC - atomsT%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atomsT%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atomsT%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atomsT%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cellT%amat, direc)
            !direcExt = direcExt / norm2(direcExt)
            !call ylmnorm_init(atomsT%lmaxd + 1)
            call ylm4(atomsT%lmaxd, direcExt, ylm) ! todo really until l + 1
            !call ylmnorm_init(atomsT%lmaxd)

            vMTjuPhon = 0
            if (iqpt /= 1) then
              do idirec = 1, 3
                do oqn_l = 0, atomsT%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atomsT%jri(itypeLock)
                      vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + 2 * real(gradVrjuPhon(imesh, lm, atomLock, idirec) * ylm(lm))
                    end do
                  end do
                end do
              end do
            else
              do idirec = 1, 3
                do oqn_l = 0, atomsT%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atomsT%jri(itypeLock)
                      vMTjuPhon(idirec, imesh) = vMTjuPhon(idirec, imesh) + gradVrjuPhon(imesh, lm, atomLock, idirec) * ylm(lm)
                    end do
                  end do
                end do
              end do
            end if
            do imesh = 1, atomsT%jri(itypeLock)
              write(1333, '(7(es20.10E3, 2x),i2, 3(2x, i2))')&
                &atomsT%rmsh(imesh, itypeLock), real(vMTjuPhon(1, imesh)) ,aimag(vMTjuPhon(1, imesh)),&
                &real(vMTjuPhon(2, imesh)), aimag(vMTjuPhon(2, imesh)), real(vMTjuPhon(3, imesh)), &
                &aimag(vMTjuPhon(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            !todo seems to be the real ptoential here
            iatom = 0
            potFleurMT = 0
              !Evaluate MT sphere function at same random points
            symAt_temp  = atomsT%ntypsy(atomLock) ! symmetry of atom iatom
            symOpr_temp = atomsT%ngopr (atomLock) ! symmetry operation mapping local to global coordinate system
            do imesh = 1, atomsT%jri(itypeLock)
              rvec = direcExt
              if (symOpr_temp /= 1) then
              !Transform into internal real space coordinates
              !todo might be wrong bmat
                call cotra1(rvec, rvec_int, cellT%bmat) !allocate and deallocate randPtsGrid
                do imatmulr = 1, 3
                  rotatedVector(imatmulr) = 0.0
                  do imatmulc = 1, 3
                    rotatedVector(imatmulr) = rotatedVector(imatmulr) + symT%mrot(imatmulr, imatmulc, symOpr_temp) * rvec_int(imatmulc)
                  end do
                end do
                call cotra0(rotatedVector, rvec, cellT%amat)
              end if
              call ylm4(atomsT%lmax(itypeLock), rvec, ylm)

              do ilath = 0, sphharT%nlh(symAt_temp)
                linCombBas = 0.0
                l_temp = sphharT%llh(ilath, symAt_temp) * (sphharT%llh(ilath, symAt_temp) + 1) + 1 !understand this
                do imem = 1, sphharT%nmem(ilath, symAt_temp)
                  lm_temp = l_temp + sphharT%mlh(imem, ilath, symAt_temp)
                  linCombBas = linCombBas + real(clnu_atom(imem, ilath, symAt_temp) * ylm(lm_temp))
                end do
                  potFleurMT(imesh) = 0!potFleurMT(imesh) + rho0MT(imesh, ilath, itypeLock, 1) * linCombBas
           !       potFleurMT(imesh) = potFleurMT(imesh) + rho0MT(imesh, ilath, itypeLock, 1) * linCombBas
              end do
            end do
           ! call fopen(1020, name='potMTfleur', status='replace', action='write', form='formatted')
           ! do imesh = 1, atomsT%jri(itypeLock)
           !   write(1020, '(3(es16.8E3, 2x))') atomsT%rmsh(imesh, itypeLock), potFleurMT(imesh), real(vMTfLEUR(1, imesh))
           ! end do
           ! call fclose(1020)

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            onlyPotFLEUR = cmplx(0., 0.)
            do iGvec = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iGvec) + qpt(:), ucpathC))
              if (iqpt /= 1) then
                vpwRjuPhon(:) = vpwRjuPhon(:) + 2 * real(vGradCoul0IR(iGvec, :) * exponential)
              else
                vpwRjuPhon(:) = vpwRjuPhon(:) + vGradCoul0IR(iGvec, :) * exponential
              end if
            !  onlyPotFLEUR = onlyPotFLEUR + vpw_G(iGvec) * exponential ! pure potential
            end do
            write (1333, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            !write (1222, '(3(es16.8E3, 2x))') x, real(onlyPotFLEUR), aimag(onlyPotFLEUR) ! todo the lines down really needed?
            write (1009, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            do idirec = 1, 3
              if ( vpwRjuPhon(idirec) /= 0 ) then
                write(1444, '(4(es16.8E3, 2x))') &
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else
                write(1444, '(4(es16.8E3, 2x))') &
                  &real(vMTjuPhon(idirec, atomsT%jri(itypeLock))), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(idirec, atomsT%jri(itypeLock))) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          end if! if still in MT then wait for crossing

        end if !is not first point on path
      end if ! is not in pseudopot mode
      wasIR = isIR
    end do
    call fclose(1333)
    call fclose(1444)
    end subroutine plot1stDensVar

  subroutine TestPlotVext2( atoms, lathar, input, stars, cell, dimens, sym, harSw, extSw, xcSw, ngdp, paPoXCo, paPoYCo, paPoZCo,   &
      & memd_atom, gdp, mlh_atom, nmem_atom, clnu_atom, vExt2IR, vExt2MT, logUnit )

    use m_juDFT_NOstopNO
    use mod_juPhonUtils, only : calcGrFinLH, calcGrFinSH, outerProduct
    use m_jpVeff1, only : GenVeff1
    use m_jpConstants, only : iu
    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),               intent(in) :: atoms
    type(t_sphhar),              intent(in) :: lathar
    type(t_input),               intent(in) :: input
    type(t_stars),               intent(in) :: stars
    type(t_cell),                intent(in) :: cell
    type(t_dimension),           intent(in) :: dimens
    type(t_sym),                 intent(in) :: sym

    ! Scalar parameters
    logical,                     intent(in) :: harSw
    logical,                     intent(in) :: extSw
    logical,                     intent(in) :: xcSw
    integer,                     intent(in) :: ngdp
    real,                        intent(in) :: paPoXCo
    real,                        intent(in) :: paPoYCo
    real,                        intent(in) :: paPoZCo
    integer,                     intent(in) :: memd_atom
    integer,  intent(in) :: logUnit

    ! Array parameters
    integer,                     intent(in) :: gdp(:, :)
    integer,                     intent(in) :: mlh_atom(:, 0:, :)
    integer,                     intent(in) :: nmem_atom(0:, :)
    complex,                     intent(in) :: clnu_atom(:, 0:, :)
    complex,                     intent(in) :: vExt2IR(:, :, :, :)
    complex,                     intent(in) :: vExt2MT(:, :, :, :)

    ! Scalar variables
    logical                           :: benchAvailSw
    integer                           :: iG
    integer                           :: ii
    integer                           :: jj
    real                              :: resolution
    integer                           :: iDatom
    integer                           :: iDtype
    integer                           :: iDeqat
    logical                           :: harSwPot
    logical                           :: extSwPot
    logical                           :: xcSwPot
    logical                           :: vExtFull
    logical                           :: vHarNum

    ! Array variables
    complex,        allocatable :: vIRstFleur(:, :)
    complex,        allocatable :: vDiffIRstFleur(:, :)
    real,           allocatable :: vMTlhFleur(:, :, :, :)
    real,           allocatable :: vDiffMTlhFleur(:, :, :, :)
    complex,        allocatable :: vIRGsFleur(:)
    complex,        allocatable :: gr2VpwFleur(:, :)
    complex,        allocatable :: grVMTshFleur(:, :, :, :)
    complex,        allocatable :: gr2VMTcolTemp(:, :, :, :)
    complex,        allocatable :: gr2VMTshFleur(:, :, :, :)
    real                        :: Gext(3)
    complex,        allocatable :: rho0IRDummy(:, :)
    complex,        allocatable :: rho0MTDummy(:, :, :, :)
    complex,        allocatable :: grRho0MTDummy(:, :, :, :)
    complex,        allocatable :: rho1PWDummy(:, :)
    complex,        allocatable :: rho1MTDummy(:, :, :, :)
    complex,        allocatable :: vxc1IRKernDummy(:)
    complex,        allocatable :: ylmDummy(:, :)
    real,           allocatable :: dKernMTGPtsDummy(:, :, :)
    real,           allocatable :: gWghtsDummy(:)
    complex,        allocatable :: vEff1IR(:, :)
    complex,        allocatable :: vEff1MT(:, :, :, :)


    allocate( vIRstFleur(stars%n3d, 1) )
    allocate( vDiffIRstFleur(stars%n3d, 1) )
    allocate( vMTlhFleur(atoms%jmtd, 0:lathar%nlhd, atoms%ntypd, input%jspins) )
    allocate( vDiffMTlhFleur(atoms%jmtd, 0:lathar%nlhd, atoms%ntypd, input%jspins) )
    allocate( vIRGsFleur(ngdp) )
    allocate( gr2VpwFleur(ngdp, 3))

    allocate(vxc1IRKernDummy(1))
    allocate(ylmDummy(1, 1))
    allocate(dKernMTGPtsDummy(1, 1, 1))
    allocate(gWghtsDummy(1))

    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw
    benchAvailSw = .true.
    vIRstFleur = 0
    vDiffIRstFleur = cmplx(0., 0.)
    vMTlhFleur = 0.
    vDiffMTlhFleur = 0.
    gr2VpwFleur = cmplx(0.0, 0.0)

    harSwPot = .false.
    extSwPot = .true.
    xcSwPot = .false.
    vHarNum = .false.

    allocate ( rho1PWDummy(ngdp, 3), rho1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3),                                         &
      & grRho0MTDummy(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
     allocate(  rho0IRDummy(ngdp, 1) )

     allocate( rho0MTDummy( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
    rho0IRDummy = cmplx(0.0, 0.0)
    rho0MTDummy = cmplx(0.0, 0.0)
    rho1PWDummy = cmplx(0.0, 0.0)
    rho1MTDummy = cmplx(0.0, 0.0)

    resolution = 600
    vExtFull = .true.
    ! todo Think about what happens in polyatomic systems
    write(*, *) 'testplotvext2 might not work for polyatomic systems'
    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        ! We choose q = 0 because we only want to have Vext2 for q = 0, because of 7.8 one can extract the q dependency so Vext1 is
        ! basically only the gradVext0, then we apply again the gradient and have Vext2. So for q = 0 we do not need to shift the
        ! the G basis set by q
        call GenVeff1( stars, cell, atoms, dimens, harSwPot, extSwPot, xcSwPot, vExtFull, ngdp, [0., 0., 0.], rho0IRDummy, rho0MTDummy, rho1PWDummy, &
          & rho1MTDummy(:, :, :, :), grRho0MTDummy, gdp, vEff1IR, vEff1MT, vxc1IRKernDummy, ylmDummy, dKernMTGPtsDummy, gWghtsDummy,&
          & iDatom, iDtype, 1, ngdp, gdp, vHarNum ) ! add spin
        do jj = 1, 3
          do iG = 1, ngdp
            Gext(:) = matmul(cell%bmat, gdp(:, iG))
            gr2VpwFleur(iG, :) = (-1) * iu * Gext(:) * vEff1IR(iG, jj)
          enddo
          call calcGrFinSH(atoms, -vEff1MT(:, :, :, jj), gr2VMTshFleur)
          call plotCompVnVbench( atoms, cell, paPoXCo, paPoYCo, paPoZCo, jj, ngdp, resolution, gdp, gr2VpwFleur(:, :),             &
            & vExt2IR(:, :, jj, iDatom), gr2VmtshFleur, vExt2MT(:, :, :, jj + (iDatom - 1) * 3) )
          deallocate(gr2VmtshFleur)
        end do
        deallocate(Veff1IR, vEff1MT)
      end do
    end do

    write(logUnit, '(a)') '---------------------------------------------------------------------'
    write(logUnit, '(a)') 'Test Vext2 routines'
    write(logUnit, *)
    write(logUnit, '(a)') '  --> Plot data Vext2Plot Vext2Cont using the plotVext2Test.py script'
    write(logUnit, '(a)') '---------------------------------------------------------------------'

!    do jj = 1, 3
!      do ii = 1, 3
!        do iG = 1, ngdp
!          write(1378, '(6(i8),2(es15.8,1x))') ii, jj, iG, gdp(1, iG), gdp(2, iG), gdp(3, iG), Vext2IR(iG, ii, jj, 1)
!        end do
!      end do
!    end do

  end subroutine TestPlotVext2

  ! Scans Eii with q through the unit cell
  subroutine ScanEii2q( atoms, cell, stars, dimens, input, ngdp, gdp )

    use m_jp2ndOrdQuant, only : CalcIIEnerg2MatElem
    use m_jpPotDensHelper, only : genPertPotDensGvecs
    use m_types

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_stars),                  intent(in) :: stars
    type(t_dimension),              intent(in) :: dimens
    type(t_input),                  intent(in) :: input

    ! Scalar parameter
    integer,                        intent(in) :: ngdp

    ! Array parameter
    integer,                        intent(in) :: gdp(:, :)

    ! Scalar variable
    integer                                    :: iAtype
    integer                                    :: iBtype
    integer                                    :: iAatom
    integer                                    :: iBatom
    integer                                    :: iAeqat
    integer                                    :: iBeqat
    integer                                    :: iAdir
    integer                                    :: iBdir
    real                                       :: x
    real                                       :: dx
    integer                                    :: ii
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: ngpqdp2km
    integer                                    :: ngpqdp
    logical                                    :: fcc = .true.
    logical                                    :: bcc = .false.

    ! Array variable
    real                                       :: qpoint(3)
    real                                       :: qpointStore(3)
    integer                                    :: gpqdp2iLim(2, 3)
    integer,           allocatable             :: gpqdp(:, :)
    integer,           allocatable             :: gpqdp2Ind(:, :, :)
    complex,           allocatable             :: E2ndOrdIIat(:, :)
    complex,           allocatable             :: E2ndOrdIIatFinQ(:, :)
    complex,           allocatable             :: E2ndOrdIIatQ0(:, :)
    complex,           allocatable             :: E2ndOrdII(:, :)

    allocate( E2ndOrdII(3 * atoms%nat, 3 * atoms%nat) )
    allocate( E2ndOrdIIat(3 * atoms%nat, 3 * atoms%nat) )
    allocate( E2ndOrdIIatFinQ(3 * atoms%nat, 3 * atoms%nat) )
    allocate( E2ndOrdIIatQ0(3 * atoms%nat, 3 * atoms%nat) )
    E2ndOrdIIat(:, :) = cmplx(0.0, 0.0)
    E2ndOrdIIatFinQ(:, :) = cmplx(0.0, 0.0)
    E2ndOrdIIatQ0(:, :) = cmplx(0.0, 0.0)

    if (.false.) then
      dx = 1. / 100.
      x = -dx
      do ii = 0, 100
        x = x + dx
        qpoint = [x, x, x]
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIatFinQ = cmplx(0.0, 0.0)
        E2ndOrdIIatQ0 = cmplx(0.0, 0.0)

        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIatQ0 )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIatFinQ )

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
                    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                      & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                      & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                      & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                  end do ! iBdir
                end do ! iBeqat
              end do ! iBtype
            end do ! iAdir
          end do ! iAeqat
        end do ! iAtype

        if (atoms%nat == 1) then
          write(2132, '(f15.8,1x,9(f15.8, 1x, f15.8))') x, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3),&
                                                         & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3),&
                                                         & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
        else if (atoms%nat == 2) then
          write(2132, '(f15.8,1x,36(f15.8, 1x, f15.8))') x, E2ndOrdII(1, 1),  E2ndOrdII(1, 2), E2ndOrdII(1, 3)&
                                                         &, E2ndOrdII(2, 1),  E2ndOrdII(2, 2), E2ndOrdII(2, 3)&
                                                         &, E2ndOrdII(3, 1),  E2ndOrdII(3, 2), E2ndOrdII(3, 3)&
                                                         &, E2ndOrdII(4, 1),  E2ndOrdII(4, 2), E2ndOrdII(4, 3)&
                                                         &, E2ndOrdII(5, 1),  E2ndOrdII(5, 2), E2ndOrdII(5, 3)&
                                                         &, E2ndOrdII(6, 1),  E2ndOrdII(6, 2), E2ndOrdII(6, 3)&
                                                         &, E2ndOrdII(1, 4),  E2ndOrdII(1, 5), E2ndOrdII(1, 6)&
                                                         &, E2ndOrdII(2, 4),  E2ndOrdII(2, 5), E2ndOrdII(2, 6)&
                                                         &, E2ndOrdII(3, 4),  E2ndOrdII(3, 5), E2ndOrdII(3, 6)&
                                                         &, E2ndOrdII(4, 4),  E2ndOrdII(4, 5), E2ndOrdII(4, 6)&
                                                         &, E2ndOrdII(5, 4),  E2ndOrdII(5, 5), E2ndOrdII(5, 6)&
                                                         &, E2ndOrdII(6, 4),  E2ndOrdII(6, 5), E2ndOrdII(6, 6)
        end if ! number of atoms
      end do ! ii
    return
    end if
    if (fcc) then
      write(*, *) 'Gamma to X'
      ! Gamma nach X
      dx = 1. / 500.
      x = 0 - dx
      do ii = 0, 500
        x = x + dx
        qpoint(1) = 0.
        qpoint(2) = x * 0.5
        qpoint(3) = x * 0.5
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIatFinQ = cmplx(0., 0.)
        E2ndOrdIIatQ0 = cmplx(0.0, 0.0)

        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIatQ0 )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIatFinQ )

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
                    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                      & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                      & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                      & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                  end do ! iBdir
                end do ! iBeqat
              end do ! iBtype
            end do ! iAdir
          end do ! iAeqat
        end do ! iAtype
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR

        if (atoms%nat == 1) then
          write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                        &  E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                        &  E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
        else if (atoms%nat == 2) then
          write(2133, '(3(f15.8,2x),1x,36(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1),  E2ndOrdII(1, 2), E2ndOrdII(1, 3)&
                                                         &, E2ndOrdII(2, 1),  E2ndOrdII(2, 2), E2ndOrdII(2, 3)&
                                                         &, E2ndOrdII(3, 1),  E2ndOrdII(3, 2), E2ndOrdII(3, 3)&
                                                         &, E2ndOrdII(4, 1),  E2ndOrdII(4, 2), E2ndOrdII(4, 3)&
                                                         &, E2ndOrdII(5, 1),  E2ndOrdII(5, 2), E2ndOrdII(5, 3)&
                                                         &, E2ndOrdII(6, 1),  E2ndOrdII(6, 2), E2ndOrdII(6, 3)&
                                                         &, E2ndOrdII(1, 4),  E2ndOrdII(1, 5), E2ndOrdII(1, 6)&
                                                         &, E2ndOrdII(2, 4),  E2ndOrdII(2, 5), E2ndOrdII(2, 6)&
                                                         &, E2ndOrdII(3, 4),  E2ndOrdII(3, 5), E2ndOrdII(3, 6)&
                                                         &, E2ndOrdII(4, 4),  E2ndOrdII(4, 5), E2ndOrdII(4, 6)&
                                                         &, E2ndOrdII(5, 4),  E2ndOrdII(5, 5), E2ndOrdII(5, 6)&
                                                         &, E2ndOrdII(6, 4),  E2ndOrdII(6, 5), E2ndOrdII(6, 6)
        end if ! number of atoms

      end do ! ii

      ! X nach K
      write(*, *) 'X to K'
      dx = 1. / 375.
      x = 0
      do ii = 1, 375
        x = x + dx
        qpoint(1) = 0.  + x * 0.375
        qpoint(2) = 0.5 + x * 0.25
        qpoint(3) = 0.5 - x * 0.125
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIatFinQ = cmplx(0., 0.)
        E2ndOrdIIatQ0 = cmplx(0.0, 0.0)
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIatQ0 )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIatFinQ )

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
                    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                      & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                      & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                      & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                  end do ! iBdir
                end do ! iBeqat
              end do ! iBtype
            end do ! iAdir
          end do ! iAeqat
        end do ! iAtype
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR

        if (atoms%nat == 1) then
          write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
        else if (atoms%nat == 2) then
          write(2133, '(3(f15.8,2x),1x,36(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1),  E2ndOrdII(1, 2), E2ndOrdII(1, 3)&
                                                         &, E2ndOrdII(2, 1),  E2ndOrdII(2, 2), E2ndOrdII(2, 3)&
                                                         &, E2ndOrdII(3, 1),  E2ndOrdII(3, 2), E2ndOrdII(3, 3)&
                                                         &, E2ndOrdII(4, 1),  E2ndOrdII(4, 2), E2ndOrdII(4, 3)&
                                                         &, E2ndOrdII(5, 1),  E2ndOrdII(5, 2), E2ndOrdII(5, 3)&
                                                         &, E2ndOrdII(6, 1),  E2ndOrdII(6, 2), E2ndOrdII(6, 3)&
                                                         &, E2ndOrdII(1, 4),  E2ndOrdII(1, 5), E2ndOrdII(1, 6)&
                                                         &, E2ndOrdII(2, 4),  E2ndOrdII(2, 5), E2ndOrdII(2, 6)&
                                                         &, E2ndOrdII(3, 4),  E2ndOrdII(3, 5), E2ndOrdII(3, 6)&
                                                         &, E2ndOrdII(4, 4),  E2ndOrdII(4, 5), E2ndOrdII(4, 6)&
                                                         &, E2ndOrdII(5, 4),  E2ndOrdII(5, 5), E2ndOrdII(5, 6)&
                                                         &, E2ndOrdII(6, 4),  E2ndOrdII(6, 5), E2ndOrdII(6, 6)
        end if ! number of atoms

      end do ! ii

      ! K nach Gamma
      write(*, *) 'K to Gamma'
      dx = 1. / 750.
      x = 0
      do ii = 1, 750
        x = x + dx
        qpoint(1) = 0.375 - x * 0.375
        qpoint(2) = 0.75 - x * 0.75
        qpoint(3) = 0.375 -x * 0.375
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIatFinQ = cmplx(0., 0.)
        E2ndOrdIIatQ0 = cmplx(0.0, 0.0)
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIatQ0 )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIatFinQ )

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
                    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                      & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                      & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                      & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                  end do ! iBdir
                end do ! iBeqat
              end do ! iBtype
            end do ! iAdir
          end do ! iAeqat
        end do ! iAtype
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR

        if (atoms%nat ==1) then
          write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
        else if (atoms%nat == 2) then
          write(2133, '(3(f15.8,2x),1x,36(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1),  E2ndOrdII(1, 2), E2ndOrdII(1, 3)&
                                                         &, E2ndOrdII(2, 1),  E2ndOrdII(2, 2), E2ndOrdII(2, 3)&
                                                         &, E2ndOrdII(3, 1),  E2ndOrdII(3, 2), E2ndOrdII(3, 3)&
                                                         &, E2ndOrdII(4, 1),  E2ndOrdII(4, 2), E2ndOrdII(4, 3)&
                                                         &, E2ndOrdII(5, 1),  E2ndOrdII(5, 2), E2ndOrdII(5, 3)&
                                                         &, E2ndOrdII(6, 1),  E2ndOrdII(6, 2), E2ndOrdII(6, 3)&
                                                         &, E2ndOrdII(1, 4),  E2ndOrdII(1, 5), E2ndOrdII(1, 6)&
                                                         &, E2ndOrdII(2, 4),  E2ndOrdII(2, 5), E2ndOrdII(2, 6)&
                                                         &, E2ndOrdII(3, 4),  E2ndOrdII(3, 5), E2ndOrdII(3, 6)&
                                                         &, E2ndOrdII(4, 4),  E2ndOrdII(4, 5), E2ndOrdII(4, 6)&
                                                         &, E2ndOrdII(5, 4),  E2ndOrdII(5, 5), E2ndOrdII(5, 6)&
                                                         &, E2ndOrdII(6, 4),  E2ndOrdII(6, 5), E2ndOrdII(6, 6)
        end if ! number of atoms

      end do ! ii

      ! Gamma nach L
      write(*, *) 'Gamma to L'
      dx = 1. / 500.
      x = 0
      do ii = 1, 500
        x = x + dx
        qpoint(1) = 0. + x * 0.5
        qpoint(2) = 0. + x * 0.5
        qpoint(3) = 0. + x * 0.5
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        E2ndOrdIIatFinQ = cmplx(0., 0.)
        E2ndOrdIIatQ0 = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIatQ0 )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIatFinQ )

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
                    E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) =          &
                      & E2ndOrdII(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))      &
                      & - E2ndOrdIIatQ0(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1)) &
                      & + E2ndOrdIIatFinQ(iBdir + 3 * (iBatom - 1), iAdir + 3 * (iAatom - 1))
                  end do ! iBdir
                end do ! iBeqat
              end do ! iBtype
            end do ! iAdir
          end do ! iAeqat
        end do ! iAtype
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR
        !E2ndOrdIIat = cmplx(0.0, 0.0)
        !call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        !do idirR = 1, 3
        !  do idirC = 1, 3
        !    E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
        !  end do ! idirC
        !end do ! idirR

        if (atoms%nat == 1) then
          write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
        else if (atoms%nat == 2) then
          write(2133, '(3(f15.8,2x),1x,36(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1),  E2ndOrdII(1, 2), E2ndOrdII(1, 3)&
                                                         &, E2ndOrdII(2, 1),  E2ndOrdII(2, 2), E2ndOrdII(2, 3)&
                                                         &, E2ndOrdII(3, 1),  E2ndOrdII(3, 2), E2ndOrdII(3, 3)&
                                                         &, E2ndOrdII(4, 1),  E2ndOrdII(4, 2), E2ndOrdII(4, 3)&
                                                         &, E2ndOrdII(5, 1),  E2ndOrdII(5, 2), E2ndOrdII(5, 3)&
                                                         &, E2ndOrdII(6, 1),  E2ndOrdII(6, 2), E2ndOrdII(6, 3)&
                                                         &, E2ndOrdII(1, 4),  E2ndOrdII(1, 5), E2ndOrdII(1, 6)&
                                                         &, E2ndOrdII(2, 4),  E2ndOrdII(2, 5), E2ndOrdII(2, 6)&
                                                         &, E2ndOrdII(3, 4),  E2ndOrdII(3, 5), E2ndOrdII(3, 6)&
                                                         &, E2ndOrdII(4, 4),  E2ndOrdII(4, 5), E2ndOrdII(4, 6)&
                                                         &, E2ndOrdII(5, 4),  E2ndOrdII(5, 5), E2ndOrdII(5, 6)&
                                                         &, E2ndOrdII(6, 4),  E2ndOrdII(6, 5), E2ndOrdII(6, 6)
        end if ! number of atoms
      end do ! ii
    else if (bcc) then
      write(*, *) 'Gamma to H'
      ! Gamma nach H
      dx = 1. / 500.
      x = 0 - dx
      do ii = 0, 500
        x = x + dx
        qpoint(1) = 0. - x * 0.5
        qpoint(2) = x * 0.5
        qpoint(3) = x * 0.5
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR

        write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                      &  E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                      &  E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)

      end do ! ii

      ! H nach P
      write(*, *) 'H to P'
      dx = 1. / 750.
      x = 0
      do ii = 1, 750
        x = x + dx
        qpoint(1) = -0.5 + x * 0.75
        qpoint(2) =  0.5 - x * 0.25
        qpoint(3) =  0.5 - x * 0.25
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR

        write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)

      end do ! ii

      ! P nach Gamma
      write(*, *) 'P to Gamma'
      dx = 1. / 250.
      x = 0
      do ii = 1, 250
        x = x + dx
        qpoint(1) = 0.25 - x * 0.25
        qpoint(2) = 0.25 - x * 0.25
        qpoint(3) = 0.25 - x * 0.25
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR

        write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)

      end do ! ii

      ! Gamma nach N
      write(*, *) 'Gamma to N'
      dx = 1. / 500.
      x = 0
      do ii = 1, 500
        x = x + dx
        qpoint(1) = 0.
        qpoint(2) = 0. + x * 0.5
        qpoint(3) = 0.
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR

        write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
      end do ! ii

      ! N nach H
      write(*, *) 'N to H'
      dx = 1. / 500.
      x = 0
      do ii = 1, 500
        x = x + dx
        qpoint(1) = 0. - x * 0.5
        qpoint(2) = 0.5
        qpoint(3) = 0. + x * 0.5
        do idirR = 1, 3
          qpoint(idirR) = mod( qpoint(idirR), 1.)
        end do ! idirR
        if ( all( qpoint(:) < 1e-12 ) ) qpoint(:) = 0.
        qpointStore(1:3) = qpoint(1:3)

        E2ndOrdII = cmplx(0.0, 0.0)
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, -qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
        call CalcIIEnerg2MatElem( atoms, cell, dimens, [0., 0., 0.], ngdp, gdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = -E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR
        E2ndOrdIIat = cmplx(0.0, 0.0)
        call CalcIIEnerg2MatElem( atoms, cell, dimens, -qpoint, ngpqdp, gpqdp, E2ndOrdIIat )
        do idirR = 1, 3
          do idirC = 1, 3
            E2ndOrdII(idirR, idirC) = E2ndOrdII(idirR, idirC) + E2ndOrdIIat(idirR, idirC)
          end do ! idirC
        end do ! idirR

        write(2133, '(3(f15.8,2x),1x,9(f15.8, 1x, f15.8))') qpointStore, E2ndOrdII(1, 1), E2ndOrdII(1, 2), E2ndOrdII(1, 3), &
                                                                       & E2ndOrdII(2, 1), E2ndOrdII(2, 2), E2ndOrdII(2, 3), &
                                                                       & E2ndOrdII(3, 1), E2ndOrdII(3, 2), E2ndOrdII(3, 3)
      end do ! ii
    end if

  end subroutine ScanEii2q

  ! we have to read in from disk and plot what we want to have
  subroutine plotPathUCScalar(cell, atoms, fIRexpCoeff, ngdp, gdp, qpt, fMTexpCoeff, paPoXCo, paPoYCo, paPoZCo, counter, iqpt, lmaxp1, nrDX )

    use m_types
    use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
    use mod_juPhonUtils, only : fopen, fclose
    use m_gaunt
    use m_ylm_old
    use m_cotra
    use m_jpPotDensHelper, only : convertStar2G
    use m_juDFT_NOstopNO

    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in) :: atoms
    type(t_cell),                  intent(in) :: cell

    ! Scalar parameters
    integer,                       intent(in) :: ngdp
    real,                          intent(in) :: paPoXCo
    real,                          intent(in) :: paPoYCo
    real,                          intent(in) :: paPoZCo
    integer,                       intent(in) :: counter
    integer,                       intent(in) :: iqpt
    logical,                       intent(in) :: lmaxp1
    integer,                       intent(in) :: nrDX

    ! Array parameters
    integer,                       intent(in) :: gdp(:, :)
    real,                          intent(in) :: qpt(:)
    complex,                       intent(in) :: fIRexpCoeff(:)
    complex,                       intent(in) :: fMTexpCoeff(:, :, :)

    ! Scalar variables
    integer                                   :: iG
    integer                                   :: iatom
    integer                                   :: itype
    integer                                   :: oqn_l
    integer                                   :: mqn_m
    integer                                   :: lm
    integer                                   :: imesh
    complex                                   :: exponential
    integer                                   :: isIR
    integer                                   :: wasIR
    integer                                   :: c1, c2, c3
    real                                      :: dx
    integer                                   :: ii
    real                                      :: x
    integer                                   :: ieqat
    integer                                   :: atomLock
    integer                                   :: itypeLock
    integer                                   :: iterations
    integer                                   :: idir
    complex                                   :: fIRonMesh

    !todo really all of ucpaths required?

    real                                      :: ucpath(3)
    real                                      :: ucpathA(3), ucpathB(3), ucpathC(3)
    real                                      :: ucpathExt(3)
    real                                      :: ucpathExtc(3), ucpathExta(3)
    integer                                   :: ext2brav(3)
    integer                                   :: ext2bravLock(3)
    logical                                   :: realPot = .true.
    real                                      :: direc(3)
    real                                      :: direcExt(3)
    character(len=11)                         :: filename
    complex,          allocatable             :: ylm(:)
    complex,          allocatable             :: fMTonMesh(:)

    if (lmaxp1) then
      allocate(ylm((atoms%lmaxd + 2)**2))
    else
      allocate(ylm((atoms%lmaxd + 1)**2))
    end if
    allocate(fMTonMesh(atoms%jmtd))

    ylm(:) = cmplx(0., 0.)
    fMTonMesh(:) = cmplx(0., 0.)

!TODO Fix bug that the interstitial part is recognized as MT part because isIR should be 1 if additionally IR point is added at MT surface
! after bisection method
    ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw

    call fopen(1000, name='pathPot', status='replace', action='write', form='formatted')

    ! Calculate vector or point on path, respectively.
    dx = 1. / real(nrDX)
    x = 0 - dx
    wasIR = -1
    do ii = 0, nrDX
      x = x + dx
      ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
      iatom = 1
      isIR = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          do c1 = -1, 1
            do c2 = -1, 1
              do c3 = -1, 1
                ext2brav = [c1, c2, c3]
                ucpathExt= matmul(cell%amat, ucpath - atoms%taual(:, iatom) - ext2brav)
                if ( norm2 ( ucpathExt ) < atoms%rmt(itype) ) then
                  isIR = iatom ! todo atomLock is redundant!
                  atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                  ext2bravLock = ext2brav
                  itypeLock = itype
                end if
              end do
            end do
          end do
          iatom = iatom + 1
        end do
      end do

      if ( wasIR == -1 ) then ! initial point of path
        if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

          fIRonMesh = cmplx(0.,0.)

          do iG = 1, ngdp ! calculate Interstitial Potential
            exponential = exp(iu * tpi * dot_product((gdp(:, iG) + qpt(:)), ucpath))
            if (iqpt /= 1) then
              fIRonMesh = fIRonMesh + 2 * real(fIRexpCoeff(iG) * exponential)! todo has been changed, seems to be added so that starting point ist correct!
            else
              fIRonMesh = fIRonMesh + fIRexpCoeff(iG) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
            end if
          end do

          write (1000, '(3(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
            &x, real(fIRonMesh), aimag(fIRonMesh), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

        end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
      else ! wasIR is available
        if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

          fIRonMesh = cmplx(0.,0.)

          do iG = 1, ngdp ! calculate Interstitial Potential
            exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpath))
            if (iqpt /= 1) then
              fIRonMesh = fIRonMesh + 2 * real(fIRexpCoeff(iG) * exponential)
            else
              fIRonMesh = fIRonMesh + fIRexpCoeff(iG) * exponential
            end if
          end do

          write (1000, '(3(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
            &x, real(fIRonMesh), aimag(fIRonMesh), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

        else if (wasIR == 0 .and. isIR > 0 ) then ! crosses MT surface into the atom

          ucpathA = ucpath ! ucpathA is in atom
          ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
          do iterations = 1, 1000 !todo does this loop have to go to 1000
            ucpathC = (ucpathA + ucpathB) / 2
            ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
            if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
            !  ! mtsphere found
              exit
            else
              if ( ( norm2(ucpathExtc) - atoms%rmt( itypeLock ) ) < 0) then
                ucpathA = ucpathC
              else
                ucpathB = ucpathC
              end if
            end if
          end do
          if (iterations > 1000) then
            write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
          end if

          fIRonMesh = cmplx(0.,0.)

          do iG = 1, ngdp ! calculate Interstitial Potential
            exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpathC))
            if (iqpt /= 1) then
              fIRonMesh = fIRonMesh + 2 * real(fIRexpCoeff(iG) * exponential)
            else
              fIRonMesh = fIRonMesh + fIRexpCoeff(iG) * exponential
            end if
          end do
          write (1000, '(3(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
            &x, real(fIRonMesh), aimag(fIRonMesh), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
          direcExt = matmul(cell%amat, direc)
          if (lmaxp1) then
            call ylmnorm_init(atoms%lmaxd + 1)
            call ylm4(atoms%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylmnorm_init(atoms%lmaxd)
          else
            call ylm4(atoms%lmaxd, direcExt, ylm) ! todo really until l + 1
          end if


          fMTonMesh = 0
          if (lmaxp1) then
            if (iqpt /= 1) then
              do oqn_l = 0, atoms%lmax(itypeLock) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                  fMTonMesh(imesh) = fMTonMesh(imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock) * ylm(lm) * exp(iu * tpi * dot_product(qpt(:), ext2bravLock(:))))
                  end do
                end do
              end do
            else
              do oqn_l = 0, atoms%lmax(itypeLock) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                  fMTonMesh(imesh) = fMTonMesh(imesh) + fMTexpCoeff(imesh, lm, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end if
          else
            if (iqpt /= 1) then
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                  fMTonMesh(imesh) = fMTonMesh(imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock) * ylm(lm) * exp(iu * tpi * dot_product(qpt(:), ext2bravLock(:))))
                  end do
                end do
              end do
            else
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                  fMTonMesh(imesh) = fMTonMesh(imesh) + fMTexpCoeff(imesh, lm, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end if
          end if

          do imesh = atoms%jri(itypeLock), 1, -1
            write(1000, '(3(es20.10E3, 2x), i2, 3(2x, i2))')&
              &-atoms%rmsh(imesh, itypeLock), real(fMTonMesh(imesh)) ,aimag(fMTonMesh(imesh)),&
              &atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
          end do

        else if (wasIR > 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

          ucpathA = ucpath ! is out of atom
          ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
          do iterations = 1, 1000 ! really so much iterations needed?
            ucpathC = (ucpathA + ucpathB) / 2
            ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
            if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              ! mtsphere found
              exit
            else
              if ( ( norm2( ucpathExtc ) - atoms%rmt( itypeLock ) ) > 0) then
                ucpathA = ucpathC
              else
                ucpathB = ucpathC
              end if
            end if
          end do
          if (iterations >= 1000) then
            write(*, *) 'Warning: Iteration not converged'
          end if

          direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
          direcExt = matmul(cell%amat, direc)

          if (lmaxp1) then
            call ylmnorm_init(atoms%lmaxd + 1)
            call ylm4(atoms%lmaxd + 1, direcExt, ylm) ! todo really until l + 1
            call ylmnorm_init(atoms%lmaxd)
          else
            call ylm4(atoms%lmaxd, direcExt, ylm) ! todo really until l + 1
          end if
!todo why is here no exponential function???
          fMTonMesh = 0
          if (lmaxp1) then
            if (iqpt /= 1) then
              do oqn_l = 0, atoms%lmax(itypeLock) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(imesh) = fMTonMesh(imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock) * ylm(lm))
                  end do
                end do
              end do
            else
              do oqn_l = 0, atoms%lmax(itypeLock) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(imesh) = fMTonMesh(imesh) + fMTexpCoeff(imesh, lm, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end if
          else
            if (iqpt /= 1) then
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(imesh) = fMTonMesh(imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock) * ylm(lm))
                  end do
                end do
              end do
            else
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(imesh) = fMTonMesh(imesh) + fMTexpCoeff(imesh, lm, atomLock) * ylm(lm)
                  end do
                end do
              end do
            end if
          end if
          do imesh = 1, atoms%jri(itypeLock)
            write(1000, '(3(es20.10E3, 2x),i2, 3(2x, i2))')&
              &atoms%rmsh(imesh, itypeLock), real(fMTonMesh(imesh)) ,aimag(fMTonMesh(imesh)),&
              &atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
          end do

          fIRonMesh = cmplx(0.,0.)
          do iG = 1, ngdp ! calculate Interstitial Potential
            exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpathC))
            if (iqpt /= 1) then
              fIRonMesh = fIRonMesh + 2 * real(fIRexpCoeff(iG) * exponential)
            else
              fIRonMesh = fIRonMesh + fIRexpCoeff(iG) * exponential
            end if
          end do
          write (1000, '(3(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
            &x, real(fIRonMesh), aimag(fIRonMesh), &
            &isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

        end if! if still in MT then wait for crossing

      end if !is not first point on path
      wasIR = isIR
    end do
    call fclose(1000)

    end subroutine plotPathUCScalar

    ! we have to read in from disk and plot what we want to have
    subroutine plotPathUCVector(cell, atoms, fIRexpCoeff, ngdp, gdp, qpt, fMTexpCoeff, paPoXCo, paPoYCo, paPoZCo, counter, iqpt, lmaxp1, nrDX )

      use m_types
      use m_jPConstants, only : iu, pi, tpi, fpi, Tmatrix
      use mod_juPhonUtils, only : fopen, fclose
      use m_gaunt
      use m_ylm_old
      use m_cotra
      use m_jpPotDensHelper, only : convertStar2G
      use m_juDFT_NOstopNO

      implicit none

      ! Type parameters
      type(t_atoms),                 intent(in) :: atoms
      type(t_cell),                  intent(in) :: cell

      ! Scalar parameters
      integer,                       intent(in) :: ngdp
      real,                          intent(in) :: paPoXCo
      real,                          intent(in) :: paPoYCo
      real,                          intent(in) :: paPoZCo
      integer,                       intent(in) :: counter
      integer,                       intent(in) :: iqpt
      logical,                       intent(in) :: lmaxp1
      integer,                       intent(in) :: nrDX

      ! Array parameters
      integer,                       intent(in) :: gdp(:, :)
      real,                          intent(in) :: qpt(:)
      complex,                       intent(in) :: fIRexpCoeff(:, :)
      complex,                       intent(in) :: fMTexpCoeff(:, :, :, :)

      ! Scalar variables
      integer                                   :: iG
      integer                                   :: iatom
      integer                                   :: itype
      integer                                   :: oqn_l
      integer                                   :: mqn_m
      integer                                   :: lm
      integer                                   :: imesh
      complex                                   :: exponential
      integer                                   :: isIR
      integer                                   :: wasIR
      integer                                   :: c1, c2, c3
      real                                      :: dx
      integer                                   :: ii
      real                                      :: x
      integer                                   :: ieqat
      integer                                   :: atomLock
      integer                                   :: itypeLock
      integer                                   :: iterations
      integer                                   :: idir

      !todo really all of ucpaths required?

      complex                                   :: fIRonMesh(3)
      real                                      :: ucpath(3)
      real                                      :: ucpathA(3), ucpathB(3), ucpathC(3)
      real                                      :: ucpathExt(3)
      real                                      :: ucpathExtc(3), ucpathExta(3)
      integer                                   :: ext2brav(3)
      integer                                   :: ext2bravLock(3)
      logical                                   :: realPot = .true.
      real                                      :: direc(3)
      real                                      :: direcExt(3)
      character(len=11)                         :: filename
      complex,          allocatable             :: ylm(:)
      complex,          allocatable             :: fMTonMesh(:, :)

      if (lmaxp1) then
        allocate(ylm((atoms%lmaxd + 2)**2))
      else
        allocate(ylm((atoms%lmaxd + 1)**2))
      end if
      allocate(fMTonMesh(3, atoms%jmtd))

      ylm(:) = cmplx(0., 0.)
      fMTonMesh(:, :) = cmplx(0., 0.)

!T  ODO Fix bug that the interstitial part is recognized as MT part because isIR should be 1 if additionally IR point is added at MT surface
!   after bisection method
      ! in this test the gradient of the Coulomb potential in the interstitial region calculated with Aaron's method must be equal to iG qpw

      call fopen(1000, name='pathPot', status='replace', action='write', form='formatted')

      ! Calculate vector or point on path, respectively.
      dx = 1. / real(nrDX)
      x = 0 - dx
      wasIR = -1
      do ii = 0, nrDX
        x = x + dx
        ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
        iatom = 1
        isIR = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            do c1 = -1, 1
              do c2 = -1, 1
                do c3 = -1, 1
                  ext2brav = [c1, c2, c3]
                  ucpathExt= matmul(cell%amat, ucpath - atoms%taual(:, iatom) - ext2brav)
                  if ( norm2 ( ucpathExt ) < atoms%rmt(itype) ) then
                    isIR = iatom ! todo atomLock is redundant!
                    atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                    ext2bravLock = ext2brav
                    itypeLock = itype
                  end if
                end do
              end do
            end do
            iatom = iatom + 1
          end do
        end do

        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

            fIRonMesh = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product((gdp(:, iG) + qpt(:)), ucpath))
              if (iqpt /= 1) then
                fIRonMesh(:) = fIRonMesh(:) + 2 * real(fIRexpCoeff(iG, :) * exponential)! todo has been changed, seems to be added so that starting point ist correct!
              else
                fIRonMesh(:) = fIRonMesh(:) + fIRexpCoeff(iG, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
              end if
            end do

            write (1000, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(fIRonMesh(1)), aimag(fIRonMesh(1)), real(fIRonMesh(2)),  aimag(fIRonMesh(2)), &
              &real(fIRonMesh(3)), aimag(fIRonMesh(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

            fIRonMesh = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpath))
              if (iqpt /= 1) then
                fIRonMesh(:) = fIRonMesh(:) + 2 * real(fIRexpCoeff(iG, :) * exponential)
              else
                fIRonMesh(:) = fIRonMesh(:) + fIRexpCoeff(iG, :) * exponential
              end if
            end do

            write (1000, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(fIRonMesh(1)), aimag(fIRonMesh(1)), real(fIRonMesh(2)),  aimag(fIRonMesh(2)), &
              &real(fIRonMesh(3)), aimag(fIRonMesh(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR > 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000 !todo does this loop have to go to 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atoms%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            fIRonMesh = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpathC))
              if (iqpt /= 1) then
                fIRonMesh(:) = fIRonMesh(:) + 2 * real(fIRexpCoeff(iG, :) * exponential)
              else
                fIRonMesh(:) = fIRonMesh(:) + fIRexpCoeff(iG, :) * exponential
              end if
            end do
            write (1000, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(fIRonMesh(1)), aimag(fIRonMesh(1)), real(fIRonMesh(2)),  aimag(fIRonMesh(2)), &
              &real(fIRonMesh(3)), aimag(fIRonMesh(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cell%amat, direc)
            !call ylmnorm_init(atoms%lmaxd + 1)
            call ylm4(atoms%lmaxd, direcExt, ylm) ! todo really until l + 1
            !call ylmnorm_init(atoms%lmaxd)

            fMTonMesh = 0
            if (iqpt /= 1) then
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(idir, imesh) = fMTonMesh(idir, imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock, idir) * ylm(lm) * exp(iu * tpi * dot_product(qpt(:), ext2bravLock(:))))
                    end do
                  end do
                end do
              end do
            else
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itypeLock)
                    fMTonMesh(idir, imesh) = fMTonMesh(idir, imesh) + fMTexpCoeff(imesh, lm, atomLock, idir) * ylm(lm)
                    end do
                  end do
                end do
              end do
            end if

            do imesh = atoms%jri(itypeLock), 1, -1
              write(1000, '(7(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atoms%rmsh(imesh, itypeLock), real(fMTonMesh(1, imesh)) ,aimag(fMTonMesh(1, imesh)),&
                &real(fMTonMesh(2, imesh)), aimag(fMTonMesh(2, imesh)), real(fMTonMesh(3, imesh)), &
                aimag(fMTonMesh(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

          else if (wasIR > 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atoms%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cell%amat, direc)
            !call ylmnorm_init(atoms%lmaxd + 1)
            call ylm4(atoms%lmaxd, direcExt, ylm) ! todo really until l + 1
            !call ylmnorm_init(atoms%lmaxd)

            fMTonMesh = 0
            if (iqpt /= 1) then
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itypeLock)
                      fMTonMesh(idir, imesh) = fMTonMesh(idir, imesh) + 2 * real(fMTexpCoeff(imesh, lm, atomLock, idir) * ylm(lm))
                    end do
                  end do
                end do
              end do
            else
              do idir = 1, 3
                do oqn_l = 0, atoms%lmax(itypeLock)! + 1
                  do mqn_m = -oqn_l, oqn_l
                    lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                    do imesh = 1, atoms%jri(itypeLock)
                      fMTonMesh(idir, imesh) = fMTonMesh(idir, imesh) + fMTexpCoeff(imesh, lm, atomLock, idir) * ylm(lm)
                    end do
                  end do
                end do
              end do
            end if
            do imesh = 1, atoms%jri(itypeLock)
              write(1000, '(7(es20.10E3, 2x),i2, 3(2x, i2))')&
                &atoms%rmsh(imesh, itypeLock), real(fMTonMesh(1, imesh)) ,aimag(fMTonMesh(1, imesh)),&
                &real(fMTonMesh(2, imesh)), aimag(fMTonMesh(2, imesh)), real(fMTonMesh(3, imesh)), &
                &aimag(fMTonMesh(3, imesh)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            fIRonMesh = cmplx(0.,0.)
            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG) + qpt(:), ucpathC))
              if (iqpt /= 1) then
                fIRonMesh(:) = fIRonMesh(:) + 2 * real(fIRexpCoeff(iG, :) * exponential)
              else
                fIRonMesh(:) = fIRonMesh(:) + fIRexpCoeff(iG, :) * exponential
              end if
            end do
            write (1000, '(7(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(fIRonMesh(1)), aimag(fIRonMesh(1)), real(fIRonMesh(2)),  aimag(fIRonMesh(2)), &
              &real(fIRonMesh(3)), aimag(fIRonMesh(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if! if still in MT then wait for crossing

        end if !is not first point on path
        wasIR = isIR
      end do
      call fclose(1000)
      end subroutine plotPathUCVector

  ! todo Later this routine should be called by the gradVeff test as well!
  subroutine plotCompVnVbench(atoms, cell, paPoXCo, paPoYCo, paPoZCo, colInd, ngdp, resolution, gdp, gr2VpwFleur, gr2VIRjuPhon, gr2VmtshFleur, &
      & gr2VMTjuPhon )

    use m_types
    use m_jPConstants, only : iu, tpi!, pi, tpi, fpi, Tmatrix
    use m_ylm_old
    use mod_juPhonUtils, only : fopen, fclose
!
    implicit none

   ! Type parameters
   type(t_atoms),      intent(in) :: atoms
   type(t_cell),       intent(in) :: cell

   ! Scalar parameters
   real,               intent(in) :: paPoXCo
   real,               intent(in) :: paPoYCo
   real,               intent(in) :: paPoZCo
   integer,            intent(in) :: colInd
   integer,            intent(in) :: ngdp
   real,               intent(in) :: resolution

   ! Array parameters
   integer,            intent(in) :: gdp(:, :)
   complex,            intent(in) :: gr2VpwFleur(:, :)
   complex,            intent(in) :: gr2VIRjuPhon(:, :)
   complex,            intent(in) :: gr2VmtshFleur(:, :, :, :)
   complex,            intent(in) :: gr2VMTjuPhon(:, :, :)


   ! Scalar variables
   logical                        :: realPot = .true.
   real                           :: dx
   character(len=14)              :: filename1
   character(len=14)              :: filename2
   integer                        :: wasIR
   integer                        :: isIR
   real                           :: x
   integer                        :: ii
   integer                        :: iatom
   integer                        :: itype
   integer                        :: ieq
   integer                        :: c1, c2, c3
   integer                        :: atomLock
   integer                        :: itypeLock
   integer                        :: iG
   complex                        :: exponential
   integer                        :: iterations
   integer                        :: idirec
   integer                        :: oqn_l
   integer                        :: mqn_m
   integer                        :: lm
   integer                        :: imesh

   ! Array variables
   real                           :: ucpath(3)
   integer                        :: ext2brav(3)
   real                           :: ucpathExt(3)
   integer                        :: ext2bravLock(3)
   complex                        :: vpwRfleur(3)
   complex                        :: vpwRjuPhon(3)
   real                           :: ucpathA(3), ucpathB(3), ucpathC(3)
   real                           :: ucpathExtc(3)!, ucpathExta(3)
   real                           :: direc(3)
   real                           :: direcExt(3)
   complex, allocatable           :: ylm(:)
   complex, allocatable           :: vMTfleur(:, :)
   complex, allocatable           :: vMTjuPhon(:, :)
   complex, allocatable           :: vFleurMT(:)

   allocate( ylm( (atoms%lmaxd + 1)**2 ) )
   allocate( vMTfleur( atoms%jmtd, 3) )
   allocate( vMTjuPhon( atoms%jmtd, 3) )
   allocate( vFleurMT( atoms%jmtd ) )



   write(filename1, '(a12,i1)') 'Vexteff2Plot', colInd
   write(filename2, '(a12,i1)') 'Vexteff2Cont', colInd

   call fopen(1333, name=filename1, status='replace', action='write', form='formatted')
   call fopen(1444, name=filename2, status='replace', action='write', form='formatted')

   dx = 1. / real(resolution)
   x = 0 - dx
   wasIR = -1
   do ii = 0, resolution
     x = x + dx
     ucpath = [paPoXCo * x, paPoYCo * x, paPoZCo * x] !CHANGE
     iatom = 0
     isIR = 0
     do itype = 1, atoms%ntype
       do ieq = 1, atoms%neq(itype)
         iatom = iatom + 1
         do c1 = -1, 1
           do c2 = -1, 1
             do c3 = -1, 1
               ext2brav = [c1, c2, c3]
               ucpathExt= matmul(cell%amat, ucpath - atoms%taual(:, iatom) - ext2brav)
               if ( norm2 ( ucpathExt ) < atoms%rmt(itype) ) then
                 isIR = iatom ! todo atomLock is redundant!
                 atomLock = iatom ! maybe also possible with two exit statements if MT found then lock variables are unneccssary
                 ext2bravLock = ext2brav
                 itypeLock = itype
               end if
             end do
           end do
         end do
       end do
     end do

      if ( realPot ) then ! if no pseudoPotMode
        if ( wasIR == -1 ) then ! initial point of path
          if ( isIR == 0) then ! is in interstitial, so point can be simply plotted

            vpwRfleur(:) = cmplx(0.,0.)
            vpwRjuPhon(:) = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + gr2VpwFleur(iG, :) * exponential
              vpwRjuPhon = vpwRjuPhon + gr2VIRjuPhon(iG, :) * exponential ! todo has been changed, seems to be added so that starting point ist correct!
            end do

            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          end if ! if not in interstitial we wait till reaching the MT surface enabling to solve the Dirichelet boundary problem
        else ! wasIR is available
          if ( ( wasIR == 0 .and. isIR == 0 )) then ! is between atoms

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG), ucpath))
              vpwRfleur(:) = vpwRfleur(:) + gr2VpwFleur(iG, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + gr2VIRjuPhon(iG, :) * exponential
            end do

            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

          else if (wasIR == 0 .and. isIR > 0 ) then ! crosses MT surface into the atom

            ucpathA = ucpath ! ucpathA is in atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] !todo make these dx vector general! ucpathB is outside atom CHANGE
            do iterations = 1, 1000
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) < 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
              !  ! mtsphere found
                exit
              else
                if ( ( norm2(ucpathExtc) - atoms%rmt( itypeLock ) ) < 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations > 1000) then
              write(*, *) 'Warning: Iteration not converged' !todo juDFT warning
            end if

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)

            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + gr2VpwFleur(iG, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + gr2VIRjuPhon(iG, :) * exponential
            end do
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), wasIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cell%amat, direc)
            direcExt = direcExt / norm2(direcExt)
  !          call ylmnorm_init(atoms%lmaxd + 2)
            call ylm4(atoms%lmaxd, direcExt, ylm)
  !          call ylmnorm_init(atoms%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    vMTfLEUR(imesh, idirec) = vMTFLEUR(imesh, idirec) + gr2VmtshFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(imesh, idirec) = vMTjuPhon(imesh, idirec) + gr2VMTjuPhon(imesh, lm, idirec + (atomLock - 1) * 3) * ylm(lm)
                  end do
                end do
              end do
            end do

            do imesh = atoms%jri(itypeLock), 1, -1
              write(1333, '(13(es20.10E3, 2x), i2, 3(2x, i2))')&
                &-atoms%rmsh(imesh, itypeLock), real(vMTfleur(imesh, 1)), aimag(vMTfleur(imesh, 1)), &
                &real(vMTfleur(imesh, 2)), aimag(vMTfleur(imesh, 2)), real(vMTfleur(imesh, 3)), &
                &aimag(vMTfleur(imesh, 3)), real(vMTjuPhon(imesh, 1)) ,aimag(vMTjuPhon(imesh, 1)),&
                &real(vMTjuPhon(imesh, 2)), aimag(vMTjuPhon(imesh, 2)), real(vMTjuPhon(imesh, 3)), &
                aimag(vMTjuPhon(imesh, 3)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec )), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          else if (wasIR > 0 .and. isIR == 0 ) then ! crosses MT surface out from atom

            ucpathA = ucpath ! is out of atom
            ucpathB = ucpath - [paPoXCo * dx, paPoYCo * dx, paPoZCo * dx] ! is in atom, todo make these dx vector general!
            do iterations = 1, 1000 ! really so much iterations needed?
              ucpathC = (ucpathA + ucpathB) / 2
              ucpathExtc = matmul(cell%amat, ucpathC - atoms%taual(:, atomLock) - ext2bravLock)
              if ( ( abs( norm2 ( ucpathExtc ) - atoms%rmt(itypeLock) ) <= 1e-7 ) .or. norm2( ucpathB - ucpathA ) < 1e-7 ) then
                ! mtsphere found
                exit
              else
                if ( ( norm2( ucpathExtc ) - atoms%rmt( itypeLock ) ) > 0) then
                  ucpathA = ucpathC
                else
                  ucpathB = ucpathC
                end if
              end if
            end do
            if (iterations >= 1000) then
              write(*, *) 'Warning: Iteration not converged'
            end if

            direc = ucpathC - atoms%taual(:, atomLock) - ext2bravLock
            direcExt = matmul(cell%amat, direc)
            direcExt = direcExt / norm2(direcExt)
!            call ylmnorm_init(atoms%lmaxd + 2)
            call ylm4(atoms%lmaxd, direcExt, ylm) ! todo really until l + 1
!            call ylmnorm_init(atoms%lmaxd)

            vMTfleur = 0
            vMTjuPhon = 0
            do idirec = 1, 3
              do oqn_l = 0, atoms%lmax(itypeLock)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itypeLock)
                    vMTfLEUR(imesh, idirec) = vMTFLEUR(imesh, idirec) + gr2VmtshFleur(imesh, lm, atomLock, idirec) * ylm(lm)
                    vMTjuPhon(imesh, idirec) = vMTjuPhon(imesh, idirec) + gr2VMTjuPhon(imesh, lm, idirec + (atomLock - 1) * 3) * ylm(lm)
                  end do
                end do
              end do
            end do
            do imesh = 1, atoms%jri(itypeLock)
              write(1333, '(13(es20.10E3, 2x),i2, 3(2x, i2))')&
                &atoms%rmsh(imesh, itypeLock), real(vMTfleur(imesh, 1)), aimag(vMTfleur(imesh, 1)), &
                &real(vMTfleur(imesh, 2)), aimag(vMTfleur(imesh, 2)), real(vMTfleur(imesh, 3)), &
                &aimag(vMTfleur(imesh, 3)), real(vMTjuPhon(imesh, 1)) ,aimag(vMTjuPhon(imesh, 1)),&
                &real(vMTjuPhon(imesh, 2)), aimag(vMTjuPhon(imesh, 2)), real(vMTjuPhon(imesh, 3)), &
                &aimag(vMTjuPhon(imesh, 3)), atomLock, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)
            end do

            vpwRfleur = cmplx(0.,0.)
            vpwRjuPhon = cmplx(0.,0.)
            do iG = 1, ngdp ! calculate Interstitial Potential
              exponential = exp(iu * tpi * dot_product(gdp(:, iG), ucpathC))
              vpwRfleur(:) = vpwRfleur(:) + gr2VpwFleur(iG, :) * exponential
              vpwRjuPhon(:) = vpwRjuPhon(:) + gr2VIRjuPhon(iG, :) * exponential
            end do
            write (1333, '(13(es16.8E3, 2x),i2, 3(2x, i2))') & ! write out calculated point
              &x, real(vpwRfleur(1)), aimag(vpwRfleur(1)), real(vpwRfleur(2)), aimag(vpwRfleur(2)),  real(vpwRfleur(3)), &
              &aimag(vpwRfleur(3)), real(vpwRjuPhon(1)), aimag(vpwRjuPhon(1)), real(vpwRjuPhon(2)),  aimag(vpwRjuPhon(2)), &
              &real(vpwRjuPhon(3)), aimag(vpwRjuPhon(3)), isIR, ext2bravLock(1), ext2bravLock(2), ext2bravLock(3)

            do idirec = 1, 3
              if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) == 0 ) .and. ( vpwRjuPhon(idirec) /= 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))) / abs(real(vpwRjuPhon(idirec)))
              else if ( ( vpwRfleur(idirec) /= 0 ) .and. ( vpwRjuPhon(idirec) == 0 ) ) then
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))) / abs(real(vpwRfleur(idirec))),&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), -9e9
              else
                write(1444, '(8(es16.8E3, 2x))') &
                  &real(vMTfleur(atoms%jri(itypeLock), idirec)), real(vpwRfleur(idirec)), &
                  &abs(real(vMTfleur(atoms%jri(itypeLock), idirec)) - real(vpwRfleur(idirec))), -9e9,&
                  &real(vMTjuPhon(atoms%jri(itypeLock), idirec)), real(vpwRjuPhon(idirec)), &
                  &abs(real(vMTjuPhon(atoms%jri(itypeLock), idirec)) - real(vpwRjuPhon(idirec))), -9e9
              end if
            end do

          end if! if still in MT then wait for crossing

        end if !is not first point on path
      end if ! is not in pseudopot mode
      wasIR = isIR
    end do
    call fclose(1333)
    call fclose(1444)
  end subroutine plotCompVnVbench

end module m_jpPlotObservables
