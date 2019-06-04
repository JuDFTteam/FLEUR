!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wann_optional
  USE m_juDFT
CONTAINS
  SUBROUTINE wann_optional(input,kpts,atoms,sym,cell,oneD,noco,wann)
    !**************************************************
    !     Make preparations for the calculation of 
    !     Wannier functions.
    !     Frank Freimuth
    !**************************************************
    USE m_types
    USE m_wann_read_inp
    USE m_wann_projgen
    USE m_wann_kpointgen
    USE m_wann_w90kpointgen
    USE m_wann_kptsreduc
    USE m_wann_kptsreduc2
    USE m_wann_wan90prep
    USE m_wann_dipole3
    USE m_wann_dipole

    IMPLICIT NONE

    TYPE(t_input),     INTENT(IN)    :: input
    TYPE(t_kpts),      INTENT(IN)    :: kpts
    TYPE(t_atoms),     INTENT(IN)    :: atoms
    TYPE(t_sym),       INTENT(IN)    :: sym
    TYPE(t_cell),      INTENT(IN)    :: cell
    TYPE(t_oneD),      INTENT(IN)    :: oneD
    TYPE(t_noco),      INTENT(IN)    :: noco
    TYPE(t_wann),      INTENT(INOUT) :: wann

    INTEGER       :: num_wann(2)
    LOGICAL       :: l_nocosoc

    l_nocosoc=noco%l_noco.OR.noco%l_soc

    !-----read the input file to determine what to do
    CALL wann_read_inp(input,.TRUE.,wann)

    !-----generate projection-definition-file
    IF(wann%l_projgen) THEN
       CALL wann_projgen(atoms%ntype,atoms%neq,atoms%nat,atoms%zatom,l_nocosoc,wann)
       wann%l_stopopt=.TRUE.
    ENDIF

    !-----generate k-point-files
    IF(wann%l_kpointgen) THEN
       CALL wann_kpointgen()
       wann%l_stopopt=.TRUE.
    ENDIF
    IF(wann%l_w90kpointgen) THEN
       CALL wann_w90kpointgen()
       wann%l_stopopt=.TRUE.
    ENDIF

    !-----find Wannier-irreducible part of BZ
    IF(wann%l_kptsreduc)THEN
       CALL wann_kptsreduc(sym%nop,sym%mrot,cell%bmat,sym%tau,input%film, oneD%odi%d1,l_nocosoc)
       wann%l_stopopt=.TRUE.
    ENDIF

    !-----find Wannier-irreducible part of BZ
    IF(wann%l_kptsreduc2)THEN
       CALL wann_kptsreduc2(wann%mhp, sym%nop,sym%mrot,cell%bmat,sym%tau,input%film, oneD%odi%d1,l_nocosoc)
       wann%l_stopopt=.TRUE.
    ENDIF

    !-----generate WF1.win and bkpts
    IF(wann%l_prepwan90)THEN
       CALL wann_wan90prep(input,kpts, input%jspins,cell%amat,cell%bmat, atoms%nat,atoms%taual,&
            atoms%zatom,atoms%ntype, atoms%ntype,atoms%neq,wann%l_bzsym,input%film, oneD%odi%d1,&
            wann%l_ms,wann%l_sgwf,wann%l_socgwf, wann%aux_latt_const,wann%param_file,wann%l_dim, &
            wann%wan90version)
    ENDIF

    !-----calculate polarization, if not wannierize
    !-----if wannierize, then calculate polarization later (after wannierize)
    IF(wann%l_dipole3.AND..NOT.wann%l_wannierize)THEN
       num_wann(1)=wann%band_max(1)-wann%band_min(1)+1
       num_wann(2)=wann%band_max(2)-wann%band_min(2)+1
       CALL wann_dipole3(input%jspins,cell%omtil,atoms%nat,atoms%pos, cell%amat,cell%bmat,atoms%taual,&
            num_wann, atoms%ntype,atoms%neq,atoms%zatom,l_nocosoc)
       wann%l_stopopt=.TRUE.
    ENDIF

    !-----calculate polarization, if not wannierize
    !-----if wannierize, then calculate polarization later (after wannierize)
    IF(wann%l_dipole.AND..NOT.wann%l_wannierize)THEN
       CALL wann_dipole(input%jspins,cell%omtil,atoms%nat,atoms%pos, cell%amat,atoms%ntype,&
            atoms%neq,atoms%zatom)
       wann%l_stopopt=.TRUE.
    ENDIF

    IF(wann%l_stopopt)  CALL juDFT_end("wann_optional done",1) ! The 1 is temporarily. Should be mpi%irank.

  END SUBROUTINE wann_optional
END MODULE m_wann_optional
