module m_jpIOnMixing

  implicit none

  contains

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Gregor Michaliczek
  !>
  !> @brief
  !> Updates and mixes density variation for next Sternheimer cycle
  !>
  !> @details
  !>
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine UpdNCheckDens( atoms, stars, cell, input, ngpqdp, stern1stIt, stern2ndIt, sternRegIt, sternFinIt, iDatom, iqpt, idir,    &
      & gpqdp, lastDistance, rho1IRout, rho1MTout)

    use m_types
    !use m_jpDens1stVar, only : loadDensity, noSymMetric, packMixVector, storeDensity, unpackMixVector, broydenNsym
    use m_stmix_old

    implicit none

    ! Type parameters
    type(t_atoms),                 intent(in)    :: atoms
    type(t_stars),                 intent(in)    :: stars
    type(t_cell),                  intent(in)    :: cell
    type(t_input),                 intent(in)    :: input

    ! Scalar parameters
    logical,                       intent(inout) :: stern1stIt
    logical,                       intent(inout) :: stern2ndIt
    logical,                       intent(inout) :: sternRegIt
    logical,                       intent(inout) :: sternFinIt
    integer,                       intent(in)    :: iDatom
    integer,                       intent(in)    :: iqpt
    integer,                       intent(in)    :: ngpqdp
    integer,                       intent(in)    :: idir

    ! Array parameters
    integer,                       intent(in)    :: gpqdp(:, :)
    complex,                       intent(inout) :: lastDistance(:, :, :)
    complex,                       intent(inout) :: rho1IRout(:)
    complex,                       intent(inout) :: rho1MTout(:, :, :)

    ! Local type
    type(t_noco)                                 :: noco

    ! Local scalars
    integer                                      :: mmap
    integer                                      :: nmap, nmapMT
    real                                         :: distR
    real                                         :: distI
 !   real                                        :: distIR
 !   real                                        :: distRI

    ! Local arrays
    real,              allocatable               :: sm(:), fsm(:)!, fsmMet(:), smMet(:)
    character(len=24)                            :: filename
    complex,           allocatable               :: rho1IRin(:)
    complex,           allocatable               :: rho1MTin(:, :, :)

#include "cpp_double.h"
    !External functions
    real CPP_BLAS_sdot
    external CPP_BLAS_sdot

    ! In first iteration we have a first guess for the linear density variation (only external contribution) which we don't mix but
    ! use immediately for the next iteration after writing it to disc. As the first mixing is between the second and third iteration, we
    ! return from this routine also if we are in the second iteration
    if ( stern1stIt .or. stern2ndIt ) then
      write(filename, '(a10,i1,a4,i1,a4,i4)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
      call storeDensity(atoms, ngpqdp, filename, rho1IRout, rho1MTout)
      if ( stern1stIt .and. idir == 3) then
        stern1stIt = .false.
        stern2ndIt = .true.
        sternFinIt = .false.
        sternRegIt = .true.
        write(*, '(a)') 'Density variation of 1st Sternheimer SCC iteration stored! No mixing!'
      else if (stern2ndIt .and. idir == 3) then
        stern2ndIt = .false.
        sternRegIt = .true.
        write(*, '(a)') 'Density variation of 2nd Sternheimer SCC iteration stored! No mixing!'
      end if
      return
    end if

    ! If this point is reached we are at a regular iteration. Load the linear density variation which has been stored to one of the *cdn1* files
    allocate( rho1IRin(ngpqdp), rho1MTin(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat) )

    rho1IRin(:) = cmplx(0.0, 0.0)
    rho1MTin(:, :, :) = cmplx(0.0, 0.0)
    write(filename, '(a10,i1,a4,i1,a4,i4)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
    call loadDensity(atoms, ngpqdp, filename, rho1IRin, rho1MTin)

    ! Both the interstitial and the muffin-tin density are complex, so we get a factor of 2
    mmap = 2 * (ngpqdp  + atoms%jmtd * (atoms%lmaxd + 1)**2 * atoms%nat)

    allocate(sm(mmap), fsm(mmap))
    sm(:) = 0.
    fsm(:) = 0.
!    allocate(fsmMet(mmap), smMet(mmap))
    ! Here the metric could also be initiated

    ! pack input and output density into mixing vector
    call packMixVector( atoms, nmap, nmapMT, ngpqdp, rho1IRin, rho1MTin, sm )
    call packMixVector( atoms, nmap, nmapMT, ngpqdp, rho1IRout, rho1MTout, fsm )
!    call packMixVector( atoms, nmap, nmapMT, ngpqdp, conjg(rho1IRin), conjg(rho1MTin), smMet )
!    call packMixVector( atoms, nmap, nmapMT, ngpqdp, conjg(rho1IRout), conjg(rho1MTout), fsmMet )

    ! Calculate difference fsm - sm into fsm
    fsm(:nmap) = fsm(:nmap) - sm(:nmap)
!    fsmMet(:nmap) = fsmMet(:nmap) - smMet(:nmap)

    ! Mix with chosen mixing method
    noco%l_noco = .false.
    if (input%imix.EQ.0) then
      write(*, *)'straight mixing'
      CALL stmix_old(atoms,input,noco, nmap,nmap,fsm, sm)
    else if (input%imix == 7) then
      write(*, '(a)') 'General Anderson Mixing of displacement direction: '
      call broydenNsym(cell, stars, atoms, input, mmap, nmap, nmapMT, nmap, fsm, sm, idir, ngpqdp, gpqdp, iqpt, iDatom)
    end if

    ! Reset new input density and write result into new input density
    rho1IRout(:) = cmplx(0.0, 0.0)
    rho1MTout(:, :, :) = cmplx(0.0, 0.0)
    call unPackMixVector(atoms, ngpqdp, sm, rho1IRout, rho1MTout)

    ! Create metric for broyden and for distance
    call noSymMetric(atoms, stars, cell, idir, ngpqdp, nmap, nmapMT, gpqdp, fsm, sm)

    ! Calculate distance and update distance array
    !todo performance of sdot?
    distR = 0.0
    distI= 0.0
 !   distIR = 0.0
 !   distRI = 0.0
    distR = CPP_BLAS_sdot(nmap / 2, fsm(1), 2, sm(1), 2)
    distI = CPP_BLAS_sdot(nmap / 2, fsm(2), 2, sm(2), 2)
!    distR = CPP_BLAS_sdot(nmap / 2, fsmMet(1), 2, sm(1), 2)
!    distI = CPP_BLAS_sdot(nmap / 2, fsmMet(2), 2, sm(2), 2)
!    distIR = CPP_BLAS_sdot(nmap / 2, fsmMet(1), 2, sm(2), 2)
!    distRI = CPP_BLAS_sdot(nmap / 2, fsmMet(2), 2, sm(1), 2)
    lastDistance(idir, iDatom, iqpt) = cmplx(1000 * sqrt( abs( (distR / cell%vol ) )), 1000 * sqrt( abs( distI / cell%vol ) ))
!    lastDistance(idir, iDatom, iqpt) = cmplx(1000 * sqrt( abs( (distR - distI) / cell%vol ) ), 1000 * sqrt( abs( (distIR + distRI) / cell%vol ) ))
!    distR = CPP_BLAS_sdot(nmap, fsmMet(1), 1, sm(1), 1)
!    lastDistance(idir, iDatom, iqpt) = cmplx(1000 * sqrt( abs( distR / cell%vol ) ), 1000 * sqrt( abs( distR / cell%vol ) ))


    ! If we only want to have one iteration
  !  if( idir == 3 ) then
  !    sternFinIt = .true.
  !  end if

    if (idir == 1) then
      write(789, '(i3,2x,i3,2(f20.5),1x)', advance='no') iDatom, iqpt, lastDistance(idir, iDatom, iqpt)
    else if (idir == 2) then
      write(789, '(2(f20.5),1x)', advance='no') lastDistance(idir, iDatom, iqpt)
    else if (idir == 3) then
      write(789, '(2(f20.5))') lastDistance(idir, iDatom, iqpt)
    end if
    !write(1030, *) idir, iDatom, lastDistance(idir, iDatom, iqpt)
    ! If distance below a tolerance the iteration has been converged
    if ( abs(lastDistance(idir, iDatom, iqpt)) < 4e-5 ) then
      sternFinIt = .true.
      return
    else
      sternFinIt = .false.
    end if

    write(filename, '(a10,i1,a4,i1,a4,i4)') 'JPcdn1_Dat', iDatom, 'Ddir', idir, 'qInd', iqpt
    call storeDensity(atoms, ngpqdp, filename, rho1IRout, rho1MTout)
!    if (idir == 3) then
!      sternFinIt = .true.
!    end if

  end subroutine UpdNCheckDens

  ! Initialize Mixing vector for mixing routines and fill it with density variations
  subroutine packMixVector(atoms, nmap, iMVecMT, ngdp, rho1IR, rho1MT, mixVec)

    use m_types

    implicit none

    ! Type parameter
    type(t_atoms)        :: atoms

    ! Scalar parameter
    integer, intent(in)  :: ngdp
    integer, intent(out)  :: nmap
    integer, intent(out) :: iMVecMT

    ! Array parameters
    complex, intent(in)  :: rho1IR(:)
    complex, intent(in)  :: rho1MT(:, :, :)
    real,    intent(out) :: mixVec(:)

    ! Local scalars
    integer              :: iMVec
    integer              :: ispin
    integer              :: iG
    integer              :: iatom
    integer              :: itype
    integer              :: ieqat
    integer              :: oqn_l
    integer              :: mqn_m
    integer              :: lm
    integer              :: imesh

    ! Initialize mixing vector
    do iMVec = 1, nmap
      mixVec(iMVec) = 0.0
    end do

    ! Pack mixing vector with linear density variation arrays
    iMVec = 0
    iMVecMT = 0
    do iG = 1, ngdp
      iMVec = iMVec + 1
      mixVec(iMVec) = real(rho1IR(iG))
      iMVec = iMVec + 1
      mixVec(iMVec) = aimag(rho1IR(iG))
    end do ! iG
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              iMVecMT = iMVecMT + 1
              iMVec = iMVec + 1
              mixVec(iMVec) = real(rho1MT(imesh, lm, iatom))
              iMVecMT = iMVecMT + 1
              iMVec = iMVec + 1
              mixVec(iMVec) = aimag(rho1MT(imesh, lm, iatom))
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

    nmap = iMVec

  end subroutine packMixVector

  ! Reset arrays for linear density variation and unpack mixed mixing vector into them
  subroutine unPackMixVector(atoms, ngdp, mixVec, rho1IR, rho1MT)

    use m_types

    implicit none

    ! Type parameter
    type(t_atoms)        :: atoms

    ! Scalar parameter
    integer              :: ngdp

    ! Array parameter
    real,    intent(in)  :: mixVec(:)
    complex, intent(out) :: rho1IR(:)
    complex, intent(out) :: rho1MT(:, :, :)

    integer              :: iMVec
    integer              :: ispin
    integer              :: iG
    integer              :: iatom
    integer              :: itype
    integer              :: ieqat
    integer              :: oqn_l
    integer              :: mqn_m
    integer              :: lm
    integer              :: imesh


    ! Reset linear density variation arrays
    rho1IR(:) = cmplx(0.0, 0.0)

    rho1MT(:, :, :) = cmplx(0.0, 0.0)

    ! Unpack mixed vector into linear density variation arrays
    iMVec = 0
    do iG = 1, ngdp
      iMVec = iMVec + 1
      rho1IR(iG) = cmplx(mixVec(iMVec), mixVec(iMVec + 1))
      iMVec = iMVec + 1
    end do ! iG
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              iMVec = iMVec + 1
              rho1MT(imesh, lm, iatom) = cmplx(mixVec(iMVec), mixVec(iMVec + 1))
              iMVec = iMVec + 1
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

  end subroutine unPackMixVector

  subroutine storeDensity(atoms, ngpqdp, filename, rho1IR, rho1MT)

    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),    intent(in) :: atoms

    ! Scalar parameters
    integer,          intent(in) :: ngpqdp

    ! Array parameters
    character(len=*), intent(in) :: filename
    complex,          intent(in) :: rho1IR(:)
    complex,          intent(in) :: rho1MT(:, :, :)

    ! Scalar variables
    integer                      :: idir
    integer                      :: iG
    integer                      :: itype
    integer                      :: ieqat
    integer                      :: iatom
    integer                      :: oqn_l
    integer                      :: mqn_m
    integer                      :: lm
    integer                      :: imesh

    open(1000, file=filename, status='replace', action='write', form='unformatted')
    ! Dimensions of IR array
    write(1000) ngpqdp, 3
    ! Write IR array
    do iG = 1, ngpqdp
      write(1000) rho1IR(iG)
    end do ! iG

    write(1000) atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              write(1000) rho1MT(imesh, lm, iatom)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype
    close(1000)

  end subroutine storeDensity

  subroutine loadDensity(atoms, ngpqdp, filename, rho1IR, rho1MT)

    use m_types, only : t_atoms
    use m_juDFT_stop, only : juDFT_error

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms

    ! Scalar parameter
    integer,                       intent(in)  :: ngpqdp

    ! Array parameters
    character(len=*),              intent(in)  :: filename
    complex,                       intent(out) :: rho1IR(:)
    complex,                       intent(out) :: rho1MT(:, :, :)

    ! Local scalars
    integer                                    :: ngpqdp_proof
    integer                                    :: dir_proof
    integer                                    :: iG
    integer                                    :: jmtd_proof
    integer                                    :: lm_proof
    integer                                    :: nat_proof
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: imesh

    ! Initialize dummy arrays (This is the actual performance-optimal way of initializing an array
    rho1IR(:) = cmplx(0.0, 0.0)

   rho1MT(:, :, :) = cmplx(0.0, 0.0)

    open(1000, file=filename, status='old', action='read', form='unformatted')
    ! Dimensions of IR array
    read(1000) ngpqdp_proof, dir_proof
    if ((ngpqdp_proof /= ngpqdp) .or. (dir_proof /= 3)) then
      call juDFT_error( 'Density file JPcdn1 inconsistent', calledby='loadDensity',                                                 &
        & hint='Dimensions of the linear IR density variation are not correct' )
    end if
    ! read IR array
    do iG = 1, ngpqdp
      read(1000) rho1IR(iG)
    end do ! iG

    read(1000) jmtd_proof, lm_proof, nat_proof, dir_proof
    write(109,*) jmtd_proof, lm_proof, nat_proof, dir_proof
    write(109,*) atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3
    if ( ( jmtd_proof /= atoms%jmtd ) .or. ( lm_proof /= (atoms%lmaxd + 1)**2 ) .or. (nat_proof /= atoms%nat) .or.                  &
      & (dir_proof /= 3) )then
      call juDFT_error( 'Density file JPcdn1 inconsistent', calledby='loadDensity', &
                 & hint='Dimensions of the linear MT density variation are not correct' )
    end if
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + mqn_m + 1
            do imesh = 1, atoms%jri(itype)
              read(1000) rho1MT(imesh, lm, iatom)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype
    close(1000)

  end subroutine loadDensity

!  subroutine archiveDensity(atoms, filename, rho1IR, rho1MT)
!
!    use m_types
!    use mod_juPhonUtils, only : fopen, fclose
!
!    implicit none
!
!    type(t_atoms),                 intent(in) :: atoms
!    character(len=*),              intent(in) :: filename
!    complex,                       intent(in) :: rho1IR(:)
!    complex,                       intent(in) :: rho1MT(:, :, :)
!
!   ! call fopen(1000, name=filename, status='replace', action='write', form='unformatted')
!   ! ! Dimensions of IR array
!   ! write(1000) ngdp, 3
!   ! ! Write IR array
!   !   do iG = 1, ngdp
!   !     write(1000) rho1IR(iG)
!   !   end do ! iG
!
!   ! write(1000) atoms%jmtd, (atoms%lmax + 2)**2, atoms%nat, 3
!   !   iatom = 0
!   !   do itype = 1, atoms%ntype
!   !     do ieqat = 1, atoms%neq(itype)
!   !       iatom = iatom + 1
!   !       do oqn_l = 0, atoms%lmax(itype) + 1
!   !         do mqn_m = -oqn_l, oqn_l
!   !           lm = oqn_l * (oqn_l + 1) + mqn_m + 1
!   !           do imesh = 1, atoms%jri(itype)
!   !             write(1000) rho1MT(imesh, lm, iatom)
!   !           end do ! imesh
!   !         end do ! mqn_m
!   !       end do ! oqn_l
!   !     end do ! ieqat
!   !   end do ! itype
!   ! call fclose(1000)
!
!  end subroutine archiveDensity
!
!  subroutine unpackArchivedDens(atoms, filename, rho1IR, rho1MT)
!
!    use m_types
!    use mod_juPhonUtils, only : fopen, fclose
!
!    implicit none
!
!    type(t_atoms),                 intent(in) :: atoms
!    character(len=*),              intent(in) :: filename
!    complex,                       intent(in) :: rho1IR(:, :)
!    complex,                       intent(in) :: rho1MT(:, :, :, :)
!
!   ! call fopen(1000, name=filename, status='replace', action='write', form='unformatted')
!   ! ! Dimensions of IR array
!   ! write(1000) ngdp, 3
!   ! ! Write IR array
!   !   do iG = 1, ngdp
!   !     write(1000) rho1IR(iG)
!   !   end do ! iG
!
!   ! write(1000) atoms%jmtd, (atoms%lmax + 2)**2, atoms%nat, 3
!   !   iatom = 0
!   !   do itype = 1, atoms%ntype
!   !     do ieqat = 1, atoms%neq(itype)
!   !       iatom = iatom + 1
!   !       do oqn_l = 0, atoms%lmax(itype) + 1
!   !         do mqn_m = -oqn_l, oqn_l
!   !           lm = oqn_l * (oqn_l + 1) + mqn_m + 1
!   !           do imesh = 1, atoms%jri(itype)
!   !             write(1000) rho1MT(imesh, lm, iatom)
!   !           end do ! imesh
!   !         end do ! mqn_m
!   !       end do ! oqn_l
!   !     end do ! ieqat
!   !   end do ! itype
!   ! call fclose(1000)
!
!  end subroutine unpackArchivedDens

  subroutine noSymMetric(atoms, stars, cell, idir, ngdp, nmap, nmapMT, gdp, mixVec, metMixVec)

    use m_types
    use m_jpSetupDynMat, only : warpIRPot

    implicit none

    ! Type parameter
    type(t_atoms), intent(in)  :: atoms
    type(t_stars), intent(in)  :: stars
    type(t_cell),  intent(in)  :: cell

    ! Scalar parameter
    integer,       intent(in)  :: idir
    integer,       intent(in)  :: ngdp
    integer,       intent(in)  :: nmap
    integer,       intent(in)  :: nmapMT

    ! Array parameters
    integer,       intent(in)  :: gdp(:, :)
    real,          intent(in)  :: mixVec(:)
    real,          intent(out) :: metMixVec(:)

    ! Local scalars
    integer                    :: itype
    real                       :: dxn
    real                       :: dxn2
    real                       :: dxn4
    integer                    :: ieqat
    integer                    :: oqn_l
    integer                    :: mqn_m
    integer                    :: iMVec
    integer                    :: iMVecMT
    integer                    :: iG
    integer                    :: imesh

    ! Local arrays
    ! Metric for MT
    real,    allocatable       :: g(:)
    complex, allocatable       :: rho1IRtemp(:, :)
    complex, allocatable       :: w_rho1IRtemp(:)

    ! metric for MT is r^2 dr/di di = r^3 dx
    ! simpson integration used, weights for first and last point: 1, weights forthe rest alternating: 2 or 4
    allocate (g(nmapMT), rho1IRtemp(ngdp, 3), w_rho1IRtemp(ngdp))
    g = 0.0
    metMixVec(:) = 0.0
    iMVec = 0
    do itype = 1, atoms%ntype
      dxn = atoms%dx(itype) / 3.0e0
      dxn2 =2.0e0 * dxn
      dxn4 =4.0e0 * dxn
      ! We run other all atoms but we don't need the index of the individual atom
      do ieqat = 1, atoms%neq(itype)
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            iMVec = iMVec + 1
            ! 1st-point weight: 1
            g(iMVec) = dxn * atoms%rmsh(1,itype)*3
            do imesh = 2, atoms%jri(itype) - 1, 2
              iMVec = iMVec + 2
              ! periodically changing weight: 4 and 2
              g(iMVec - 1) = dxn4 * atoms%rmsh(imesh,     itype)**3
              g(iMVec)     = dxn2 * atoms%rmsh(imesh + 1, itype)**3
            end do
            ! Last-point weight: 1. Special indexing due to loop header of previous loop
            iMVec = iMVec + 1 - MOD(atoms%jri(itype),2)
            g(iMVec) = dxn * atoms%rmsh(atoms%jri(itype),itype)**3
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

    ! Apply metric to interstitial region (metric here = step function)
    ! + multiply metric g with s_in for MT (store in sout)
    ! map s_in on a complex help array ag3
    iMVec = 0
    do iG = 1, ngdp
      iMVec = iMVec + 1
      rho1IRtemp(iG, idir) = cmplx(mixVec(iMVec), mixVec(iMVec + 1))
      iMVec = iMVec + 1
    end do

    ! todo move out direction from routine
    call warpIRPot(stars, ngdp, idir, gdp, rho1IRtemp, w_rho1IRtemp(:))

    ! interstitial
    iMVec = 0
    do iG = 1, ngdp
      iMVec = iMVec + 1
      MetMixVec(iMVec) = cell%omtil * real(w_rho1IRtemp(iG))
      iMVec = iMVec + 1
      MetMixVec(iMVec) = cell%omtil * aimag(w_rho1IRtemp(iG))
     ! write(*, *) 'nicht 1/omtil'
    end do
    ! MT
    do iMVecMT = iMVec + 1, nmap
      MetMixVec(iMVecMT) = g(iMVecMT - iMVec) * mixVec(iMVecMT)
    end do

  end subroutine noSymMetric


  !################################################################
  !     IMIX = 3 : BROYDEN'S FIRST METHOD
  !     IMIX = 5 : BROYDEN'S SECOND METHOD
  !     IMIX = 7 : GENERALIZED ANDERSEN METHOD
  !     sm   : input charge density of iteration m
  !            afterwards update rho(m+1)
  !     fm   : output minus input charge density of iteration m
  !     sm1  : input charge density of iteration m-1
  !     fm1   : output minus inputcharge density of iteration m-1
  !################################################################
  ! Taken from workstation release branch fleur at mix/broyden folder
  SUBROUTINE broydenNsym(cell, stars, atoms, input, mmap, nmaph, nmapmt, nmap, fm, sm, idir, ngdp, gdp, iqpt, iDatom)

    USE m_types
    USE m_broyd_io
    use m_juDFT_stop, only : juDFT_error

    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_atoms),INTENT(IN)   :: atoms

    TYPE(t_hybdat)             :: hybdat
    ! Scalar Arguments
    INTEGER, INTENT (IN)        :: mmap,nmap
    INTEGER, INTENT (IN)        :: nmapmt
    integer, intent(in) :: iqpt
    integer, intent(in) :: iDatom
    integer, intent(in) :: idir
    integer,                       intent(in)    :: ngdp

    ! Array Arguments
    REAL,    INTENT (IN)    :: fm(nmap)
    REAL,    INTENT (INOUT) :: sm(nmap)
    integer,                       intent(in)    :: gdp(:, :)

    ! Local Scalars
    INTEGER         :: i,it,k,nit,npos,iread,nmaph, mit
    REAL            :: bm,dfivi,fmvm,smnorm,vmnorm,alphan
    LOGICAL         :: l_exist
    REAL, PARAMETER :: one=1.0, zero=0.0

    ! Local Arrays
    REAL, ALLOCATABLE :: am(:)
    REAL, ALLOCATABLE :: fm1(:),sm1(:),ui(:),um(:),vi(:),vm(:)

#include "cpp_double.h"
    ! External Functions
    REAL CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot

    ! External Subroutines
    EXTERNAL CPP_BLAS_saxpy,CPP_BLAS_sscal

    dfivi = zero


    ALLOCATE (fm1(mmap),sm1(mmap),ui(mmap),um(mmap),vi(mmap),vm(mmap))
    ALLOCATE (am(input%maxiter+1))

    fm1 = 0.0
    sm1 = 0.0
    ui  = 0.0
    um  = 0.0
    vi  = 0.0
    vm  = 0.0
    am  = 0.0

    mit = 0
    hybdat%l_calhf = .false.
    l_exist = initBroydenHistory(input,hybdat,nmap, idir, iqpt, iDatom) ! returns true if there already exists a Broyden history
    IF(.NOT.l_exist) mit = 1

    IF (mit.NE.1) THEN
       ! load input charge density (sm1) and difference of
       ! in and out charge densities (fm1) from previous iteration (m-1)

       CALL readLastIterInAndDiffDen(hybdat,nmap,mit,alphan,sm1(:nmap),fm1(:nmap), idir, iqpt, iDatom)
       IF (ABS(input%alpha-alphan) > 0.0001) THEN
          WRITE (6,*) 'mixing parameter has been changed; reset'
          WRITE (6,*) 'broyden algorithm or set alpha to',alphan
          CALL juDFT_error("mixing parameter (input) changed", calledby ="broyden")
       END IF

       ! generate F_m   - F_(m-1)  ... sm1
       !      and rho_m - rho_(m-1) .. fm1
       sm1(1:nmap) = sm(1:nmap) - sm1(1:nmap)
       fm1(1:nmap) = fm(1:nmap) - fm1(1:nmap)
    END IF

    ! save F_m and rho_m for next iteration
    nit = mit +1
    IF (nit > input%maxiter+1) nit = 1
    CALL writeLastIterInAndDiffDen(hybdat,nmap,nit,input%alpha,sm,fm, idir, iqpt, iDatom)

    IF (mit.EQ.1) THEN
       !     update for rho for mit=1 is straight mixing
       !     sm = sm + alpha*fm
       CALL CPP_BLAS_saxpy(nmap,input%alpha,fm,1,sm,1)
    ELSE
       !     |vi> = w|vi>
       !     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
       um(:nmap) = input%alpha * fm1(:nmap) + sm1(:nmap)
       iread = MIN(mit-1,input%maxiter+1)
       DO it = 2, iread
          CALL readUVec(input,hybdat,nmap,it-mit,mit,ui, idir, iqpt, iDatom)
          CALL readVVec(input,hybdat,nmap,it-mit,mit,dfivi,vi, idir, iqpt, iDatom)

          am(it) = CPP_BLAS_sdot(nmap,vi,1,fm1,1)
          ! calculate um(:) = -am(it)*ui(:) + um
          CALL CPP_BLAS_saxpy(nmap,-am(it),ui,1,um,1)
          WRITE(6,FMT='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it)
       END DO

       IF (input%imix.EQ.3) THEN
          !****************************************
          !        broyden's first method
          !****************************************

          ! convolute drho(m) with the metric: |fm1> = w|sm1>
         ! CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
         !             mmap,nmaph,mapmt,mapvac2,sm1,fm1,l_pot)

          ! calculate the norm of sm1 : <sm1|w|sm1>
          smnorm = CPP_BLAS_sdot(nmap,sm1,1,fm1,1)

          ! generate vm = alpha*sm1  - \sum <ui|w|sm1> vi
          vm(:) = input%alpha * fm1(:)
          DO it = 2,iread
             CALL readUVec(input,hybdat,nmap,it-mit,mit,ui, idir, iqpt, iDatom)
             CALL readVVec(input,hybdat,nmap,it-mit,mit,dfivi,vi, idir, iqpt, iDatom)
             bm = CPP_BLAS_sdot(nmap,ui,1,fm1,1)
             ! calculate vm(:) = -bm*vi(:) + vm
             CALL CPP_BLAS_saxpy(nmap,-bm,vi,1,vm,1)
             !write(6,FMT='(5x,"<ui|w|Fm> for it",i2,5x,f10.6)') it, bm
          END DO

          ! complete evaluation of vm
          ! vmnorm = <um|w|sm1>-<sm1|w|sm1>
          vmnorm = CPP_BLAS_sdot(nmap,fm1,1,um,1) - smnorm
          ! if (vmnorm.lt.tol_10) NOstopNO

          CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)

       ELSE IF (input%imix.EQ.5) THEN
          !****************************************
          !      broyden's second method
          !****************************************

          ! multiply fm1 with metric matrix and store in vm:  w |fm1>
         ! CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
         !             mmap,nmaph,mapmt,mapvac2,fm1,vm,l_pot)

          ! calculate the norm of fm1 and normalize vm it: vm = wfm1 / <fm1|w|fm1>
          vmnorm = one / CPP_BLAS_sdot(nmap,fm1,1,vm,1)
          CALL CPP_BLAS_sscal(nmap,vmnorm,vm,1)

       ELSE IF (input%imix.EQ.7) THEN
          !****************************************
          !      generalized anderson method
          !****************************************

          ! calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
          ! convolute fm1 with the metrik and store in vm
         ! CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
         !             mmap,nmaph,mapmt,mapvac2,fm1,vm,l_pot)
           call noSymMetric(atoms, stars, cell, idir, ngdp, nmap, nmapMT, gdp, fm1, vm)

          DO it = 2,iread
             CALL readVVec(input,hybdat,nmap,it-mit,mit,dfivi,vi, idir, iqpt, iDatom)
             ! calculate vm(:) = -am(it)*dfivi*vi(:) + vm
             CALL CPP_BLAS_saxpy(nmap,-am(it)*dfivi,vi,1,vm,1)
          END DO

          vmnorm = CPP_BLAS_sdot(nmap,fm1,1,vm,1)
          ! if (vmnorm.lt.tol_10) NOstopNO

          ! calculate vm(:) = (1.0/vmnorm)*vm(:)
          CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)

          ! save dfivi(mit) for next iteration
          dfivi = vmnorm

       END IF

       ! write um,vm and dfivi on file broyd.?

       npos=mit-1
       IF (mit.GT.input%maxiter+1) npos = MOD(mit-2,input%maxiter)+1

       CALL writeUVec(input,hybdat,nmap,mit,um, idir, iqpt, iDatom)
       CALL writeVVec(input,hybdat,nmap,mit,dfivi,vm, idir, iqpt, iDatom)

       ! update rho(m+1)
       ! calculate <fm|w|vm>
       fmvm = CPP_BLAS_sdot(nmap,vm,1,fm,1)
       ! calculate sm(:) = (1.0-fmvm)*ui(:) + sm
       CALL CPP_BLAS_saxpy(nmap,one-fmvm,um,1,sm,1)
    END IF

    DEALLOCATE (fm1,sm1,ui,um,vi,vm,am)

  END SUBROUTINE broydenNsym

  subroutine storeZ1nG( atoms, ikpt, iqpt, mapKpq2K, iDatom, nobd, nv, z1nG )

    use m_types, only : t_atoms
    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in) :: atoms

    ! Scalar parameter
    integer,                       intent(in) :: ikpt
    integer,                       intent(in) :: iqpt
    integer,                       intent(in) :: iDatom

    ! Array parameter
    integer,                       intent(in) :: mapKpq2K(:, :) ! todo possible error here
    complex,                       intent(in) :: z1nG(:, :, :)
    integer,                       intent(in) :: nobd(:, :)
    integer,                       intent(in) :: nv(:, :)

    ! Scalar variables
    integer                                   :: iG
    integer                                   :: iband
    integer                                   :: idir

    ! Array variables
    character(len=:), allocatable             :: filename
    character(len=15)                         :: filenameTemp

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

    open( 1000, file=filename, status='replace', action='write', form='unformatted')
    rewind(1000)

    ! We have a loop over the occupied bands at k, not k + q, as the kets in the Sternheimer equation feature a k-point k and not k + q.
    do iband = 1, nobd(ikpt, 1)
      do idir = 1, 3
        ! The G-basis vectors of the z1nG come from a z0 evaluated at k + q
        do iG = 1, nv(1, mapKpq2K(ikpt, iqpt)) + atoms%nlotot
          write(1000) z1nG(iG, iband, idir)
        end do ! iG
      end do ! iband
    end do ! idir

    close(1000)

  end subroutine storeZ1nG

end module m_jpIOnMixing
