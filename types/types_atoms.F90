!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifndef CPP_MANAGED
#define CPP_MANAGED 
#endif
MODULE m_types_atoms
  use m_types_econfig

  TYPE t_utype
      SEQUENCE
      REAL :: u, j         ! the actual U and J parameters
      REAL :: theta,phi   !the rotation angles by which the density metrics is rotated
      INTEGER :: l        ! the l quantum number to which this U parameter belongs
      INTEGER :: atomType ! The atom type to which this U parameter belongs
      LOGICAL :: l_amf ! logical switch to choose the "around mean field" LDA+U limit
   END TYPE t_utype


   TYPE t_atoms
      !<no of types
      INTEGER :: ntype
      !<total-no of atoms
      INTEGER :: nat
      !<dimensions of LO's
      INTEGER ::nlod
      INTEGER ::llod
      INTEGER ::nlotot
      !lmaxd=maxval(lmax)
      INTEGER:: lmaxd
      ! no of lda+us
      INTEGER ::n_u
      ! dimensions
      INTEGER :: jmtd
      !No of element
      INTEGER, ALLOCATABLE ::nz(:)
      !atoms per type
      INTEGER, ALLOCATABLE::neq(:)
      !radial grid points
      INTEGER, ALLOCATABLE::jri(:)
      !core states
      TYPE(t_econfig),ALLOCATABLE::econf(:)
      !lmax
      INTEGER, ALLOCATABLE::lmax(:)
      !lmax non-spherical
      INTEGER, ALLOCATABLE::lnonsph(:)
      !expansion of pseudo-charge
      INTEGER, ALLOCATABLE::ncv(:)
      !no of LO
      INTEGER, ALLOCATABLE::nlo(:)
      !l of LO (nlo,ntype)
      INTEGER, ALLOCATABLE::llo(:, :)
      !lmax for lapw (ntype)
      INTEGER, ALLOCATABLE::lapw_l(:)
      !first LO with a given l (max(nlo
      INTEGER, ALLOCATABLE::lo1l(:, :)
      !??
      INTEGER, ALLOCATABLE::ulo_der(:, :)
      !no of LOs per l (max(nlo1),ntype
      INTEGER, ALLOCATABLE::nlol(:, :)
      !true if LO is formed by \dot u (
      LOGICAL, ALLOCATABLE::l_dulo(:, :)
      !no of op that maps atom into
      INTEGER, ALLOCATABLE::ngopr(:)
      !symetry of atom (nat)
      INTEGER, ALLOCATABLE::ntypsy(:)
      !no of sphhar for atom type(ntype
      INTEGER, ALLOCATABLE ::nlhtyp(:)
      !atom mapped to by inversion (nat
      INTEGER, ALLOCATABLE ::invsat(:)
      !Calaculate forces for this atom?
      LOGICAL, ALLOCATABLE :: l_geo(:)
      !MT-Radius (ntype)
      REAL, ALLOCATABLE CPP_MANAGED::rmt(:)
      !log increment(ntype)
      REAL, ALLOCATABLE::dx(:)
      !vol of MT(ntype)
      REAL, ALLOCATABLE::volmts(:)
      !radial grid points(max(jri),ntyp
      REAL, ALLOCATABLE::rmsh(:, :)
      !charge of nucleus(ntype)
      REAL, ALLOCATABLE::zatom(:)
      !initial mag moment(ntype)
      REAL, ALLOCATABLE::bmu(:)
      !pos of atom (absol) (3,nat)
      REAL, ALLOCATABLE::pos(:, :)
      !pos of atom (relat)(3,nat)
      REAL, ALLOCATABLE CPP_MANAGED::taual(:, :)
      !labels
      CHARACTER(LEN=20), ALLOCATABLE :: label(:)
      CHARACTER(len=20), ALLOCATABLE :: speciesName(:)
      !name and other data of explicitely provided xc functional
      CHARACTER(len=4), ALLOCATABLE :: namex(:)
      INTEGER, ALLOCATABLE :: icorr(:)
      INTEGER, ALLOCATABLE :: igrd(:)
      INTEGER, ALLOCATABLE :: krla(:)
      LOGICAL, ALLOCATABLE :: relcor(:)
      !lda_u information(ntype)
      TYPE(t_utype), ALLOCATABLE::lda_u(:)
      INTEGER, ALLOCATABLE :: relax(:, :) !<(3,ntype)
      INTEGER, ALLOCATABLE :: nflip(:) !<flip magnetisation of this atom
  ! CONTAINS
  !    procedure :: nsp => calc_nsp_atom
    contains
      PROCEDURE :: same_species
   END TYPE t_atoms

 contains
   LOGICAL function same_species(atoms,n,nn)
     use m_judft
     implicit none
     class(t_atoms),INTENT(IN)::atoms
     integer,intent(in)::n,nn

     if (n>atoms%ntype.or.nn>atoms%ntype) call judft_error("Same species checked for non-existing atom")

     same_species=atoms%nz(n)==atoms%nz(nn)
     same_species=same_species.and.atoms%jri(n)==atoms%jri(nn)
     same_species=same_species.and.atoms%dx(n)==atoms%dx(nn)
     same_species=same_species.and.atoms%rmt(n)==atoms%rmt(nn)
     same_species=same_species.and.atoms%lmax(n)==atoms%lmax(nn)
     same_species=same_species.and.atoms%lnonsph(n)==atoms%lnonsph(nn)
     same_species=same_species.and.atoms%nlo(n)==atoms%nlo(nn)
     if (atoms%nlo(n)==atoms%nlo(nn)) same_species=same_species.and.all(atoms%llo(:,n)==atoms%llo(:,nn))
     same_species=same_species.and.atoms%lapw_l(n)==atoms%lapw_l(nn)
     same_species=same_species.and.atoms%l_geo(n)==atoms%l_geo(nn)
     same_species=same_species.and.trim(atoms%econf(n)%coreconfig)==trim(atoms%econf(nn)%coreconfig)
     same_species=same_species.and.trim(atoms%econf(n)%valenceconfig)==trim(atoms%econf(nn)%valenceconfig)
     same_species=same_species.and.trim(atoms%econf(n)%valenceconfig)==trim(atoms%econf(nn)%valenceconfig)
  end function

 END MODULE m_types_atoms
