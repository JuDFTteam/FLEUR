!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

MODULE m_inpgen_make_inp
   IMPLICIT NONE
   private
   PUBLIC :: make_inp_xml
CONTAINS

   

   SUBROUTINE make_inp_xml(simple_inp_file_id, filename,  l_addPath, &
                           profileName, l_include, l_explicit,l_onlyIdentitySym)
      USE m_juDFT
      USE m_read_inpgen_input
      USE m_make_crystal
      USE m_make_atomic_defaults
      USE m_make_defaults
      USE m_make_kpoints
      USE m_make_magnetism
      USE m_winpxml

      USE m_types_input
      USE m_types_atoms
      USE m_types_cell
      USE m_types_sym
      USE m_types_noco
      USE m_types_vacuum
      USE m_types_banddos
      USE m_types_hybinp
      USE m_types_xcpot_inbuild_nofunction
      USE m_types_forcetheo
      USE m_types_kpts
      USE m_types_gfinp
      USE m_types_hub1inp
      USE m_types_enpara
      USE m_types_sliceplot
      USE m_types_stars
      USE m_types_dfpt
      USE m_types_profile
      USE m_constants
      USE m_types_enparaXML
      USE m_types_mpinp

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN) :: simple_inp_file_id
      CHARACTER(LEN=*), INTENT(IN) :: filename,  profileName
      LOGICAL, INTENT(IN) :: l_addPath, l_include(4), l_explicit(4),l_onlyIdentitySym

      ! Type declarations
      TYPE(t_input)            :: input
      TYPE(t_atoms)            :: atoms
      TYPE(t_cell)             :: cell
      TYPE(t_sym)              :: sym
      TYPE(t_noco)             :: noco
      TYPE(t_vacuum)           :: vacuum
      TYPE(t_banddos)          :: banddos
      TYPE(t_mpinp)            :: mpinp
      TYPE(t_hybinp)           :: hybinp
      TYPE(t_xcpot_inbuild_nf) :: xcpot
      TYPE(t_enpara)           :: enpara
      TYPE(t_forcetheo)        :: forcetheo
      TYPE(t_kpts)             :: kpts(2)
      TYPE(t_sliceplot)        :: sliceplot
      TYPE(t_stars)            :: stars
      TYPE(t_gfinp)            :: gfinp
      TYPE(t_hub1inp)          :: hub1inp
      TYPE(t_enparaXML)        :: enparaxml
      TYPE(t_dfpt)           :: dfpt
      TYPE(t_profile)          :: profile

      ! Local variables
      INTEGER                  :: iKpts, numKpts
      CHARACTER(LEN=40)        :: kptsSelection
      CHARACTER(LEN=40)        :: kpts_str(2)
      CHARACTER(LEN=40)        :: kptsName(2)
      CHARACTER(LEN=500)       :: kptsPath(2)
      INTEGER                  :: kptsBZintegration(2)
      LOGICAL                  :: kptsGamma(2), l_check

      REAL,    ALLOCATABLE :: atompos(:, :),atomid(:),mag_mom(:,:)
      CHARACTER(len=20), ALLOCATABLE :: atomLabel(:)
      CHARACTER(len=40) :: k_set_selected

      ! Initialize k-point arrays
      kpts_str(:) = ""
      kptsPath(:) = ""
      kptsName(:) = ""
      kptsBZintegration = BZINT_METHOD_HIST
      kptsGamma = .FALSE.

      ! Initialize profile
      CALL profile%init(profileName)

      ! Read input parameters
      CALL read_inpgen_input(profile, atompos, atomid, mag_mom, atomlabel, &
                            kpts_str(1), kptsName(1), kptsPath(1), kptsBZintegration(1), &
                            kptsGamma(1), input, sym, noco, vacuum, stars, xcpot, &
                            cell, hybinp, simple_inp_file_id)

      ! Determine the space group and map atoms to groups
      CALL make_crystal(input%film, l_onlyIdentitySym, atomid, atompos, mag_mom, atomlabel, &
                       vacuum%dvac, noco, cell, sym, atoms)

      ! Generate magnetic settings
      CALL make_magnetism(input, noco, atoms, mag_mom)

      ! Set all atom-related parameters
      CALL make_atomic_defaults(input, vacuum, profile, cell, atoms, enpara)

      ! Set all defaults that have not been specified before
      CALL make_defaults(atoms, sym, cell, vacuum, input, stars, xcpot, profile, &
                        noco, banddos, mpinp, hybinp, sliceplot)

      ! Setup k-points
      numKpts = 2
      IF (l_addPath) THEN
         l_check = .TRUE.
         CALL add_special_points_default(kpts(numKpts), input%film, cell, l_check)
         IF (l_check) THEN
            kpts_str(2) = 'band=240'
            kptsPath(2) = 'default'
            WRITE(kptsName(2), '(a,i0)') "path-", 2
         ELSE
            numKpts = 1
         END IF
      END IF

      DO iKpts = 1, numKpts
         CALL make_kpoints(kpts(iKpts), cell, sym,  input%film, &
                          noco%l_ss .OR. noco%l_soc, kptsBZintegration(iKpts), &
                          kptsGamma(ikpts).or. hybinp%l_hybrid, kpts_str(iKpts), kptsName(iKpts), &
                          kptsPath(iKpts),.false.)
                        
      END DO

      input%bz_integration = kptsBZintegration(1)
      k_set_selected=kpts(1)%kptsName
      ! Write inp.xml file
      CALL w_inpxml(atoms, vacuum, input, stars, sliceplot, forcetheo, banddos, &
                   dfpt, cell, sym, xcpot, noco, mpinp, hybinp, kpts(:numkpts), &
                   k_set_selected, enpara, gfinp, hub1inp, l_explicit, l_include, &
                   filename)


   END SUBROUTINE make_inp_xml

END MODULE m_inpgen_make_inp