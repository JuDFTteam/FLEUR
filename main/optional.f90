!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_optional
  USE m_juDFT
CONTAINS
  SUBROUTINE OPTIONAL(mpi, atoms,sphhar,vacuum,DIMENSION,&
       stars,input,sym, cell, sliceplot,obsolete, xcpot, noco, oneD)
    !
    !----------------------------------------
    ! this routine is called by: fleur.F90
    !
    ! optional --+-- plot -+- loddop
    !            |         +- outcdn -+- cotra0
    !            |                    +- cotra1
    !            |                    +- starf2 -- spgrot
    !            |                    +- starf3
    !            |                    +- ylm3
    !            +-- stden -+- atom2 -+- setcor
    !            |          |         +- stpot1
    !            |          |         +- differ -+- inwint
    !            |          |         |          +- outint
    !            |          |         +- vxcall (-> see vgen.F) or:
    !            |          |         +- potl0 -+- grdchlh
    !            |          |         |         +- mkgl0
    !            |          |         |         +- vxcallg (-> see vgen.F)
    !            |          |         +- intgr1
    !            |          +- cdnovlp -+- spgrot
    !            |          |           +- rcerf --wofz
    !            |          |           +- diflgr
    !            |          |           +- qpw_to_nmt -+- phasy1 -+- spgrot
    !            |          |                          |          +- ylm3
    !            |          |                          +- sphbes
    !            |          +- qfix -- cdntot -+- intgr3
    !            |          |                  +- qsf
    !            |          |                  +- pwint -- spgrot
    !            |          +- wrtdop
    !            |          +- points -- qranf
    !            |          +- sphpts -- qranf
    !            |          +- checkdop -+- starf3
    !            |                       +- cotra0
    !            |                       +- starf2 -- spgrot
    !            |                       +- fitchk
    !            |                       +- cotra1
    !            |                       +- ylm3
    !            +-- cdnsp -+- readDensity
    !            |          +- writeDensity
    !            |          +- intgr3
    !            +-- flipcdn -+- readDensity
    !            |            +- writeDensity
    !            +-- f2u -- wrtdop
    !            +-- u2f -- loddop
    !            +-- bmt -+- readDensity
    !                     +- wrtdop
    !----------------------------------------
    USE m_bmt
    USE m_plotdop
    USE m_stden
    USE m_cdnsp
    USE m_flipcdn
    USE m_cdn_io
    USE m_f2u
    USE m_u2f
    USE m_types
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..

    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_dimension),INTENT(IN):: DIMENSION
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_obsolete),INTENT(IN) :: obsolete
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(INOUT) :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_xcpot),INTENT(IN)    :: xcpot
    TYPE(t_sliceplot),INTENT(IN):: sliceplot
    !     ..
    !     .. Local Scalars ..
    INTEGER :: it, archiveType
    CHARACTER*10 :: cdnfname
    LOGICAL :: strho
    !     ..
    it = 1

    IF (mpi%irank == 0) THEN
10     IF (sliceplot%plpot) input%score = .FALSE.
       IF (sliceplot%iplot) THEN
          CALL timestart("Plotting")
          IF (input%strho)  CALL juDFT_error("strho = T and iplot=T",calledby&
               &        ="optional")
          IF (noco%l_noco) THEN
             cdnfname = 'cdn'
             CALL plotdop(&
                  &           oneD,stars,vacuum,sphhar,atoms,&
                  &           input,sym,cell,sliceplot,&
                  &           noco%l_noco,cdnfname)
             cdnfname = 'mdnx'
             CALL plotdop(&
                  &           oneD,stars,vacuum,sphhar,atoms,&
                  &           input,sym,cell,sliceplot,&
                  &           noco%l_noco,cdnfname)
             cdnfname = 'mdny'
             CALL plotdop(&
                  &           oneD,stars,vacuum,sphhar,atoms,&
                  &           input,sym,cell,sliceplot,&
                  &           noco%l_noco,cdnfname)
             cdnfname = 'mdnz'
             CALL plotdop(&
                  &           oneD,stars,vacuum,sphhar,atoms,&
                  &           input,sym,cell,sliceplot,&
                  &           noco%l_noco,cdnfname)
          ELSE
             IF (sliceplot%slice) THEN
                cdnfname = 'cdn_slice'
             ELSE
                cdnfname = 'cdn1'
             ENDIF
             CALL plotdop(&
                  &           oneD,stars,vacuum,sphhar,atoms,&
                  &           input,sym,cell,sliceplot,&
                  &           noco%l_noco,cdnfname)
          ENDIF
          CALL timestop("Plotting")
       END IF
    ENDIF ! mpi%irank == 0
    !
    !     --->generate starting charge density
    !
    strho=input%strho
    IF (.NOT.(strho.OR.obsolete%l_f2u.OR.obsolete%l_u2f.OR.sliceplot%iplot)) THEN
       archiveType = CDN_ARCHIVE_TYPE_CDN1_const
       IF (noco%l_noco) THEN
          archiveType = CDN_ARCHIVE_TYPE_NOCO_const
       END IF
       strho = .NOT.isDensityFilePresent(archiveType)
    ENDIF

    IF (strho) THEN
       strho=input%total 
       input%total = .FALSE.
       !
       CALL timestart("generation of start-density")
       CALL stden(mpi,&
            &              sphhar,stars,atoms,sym,&
            &              DIMENSION,vacuum,&
            &              input,&
            &              cell,&
            &              xcpot,&
            &              obsolete,&
            &              oneD)
       !
       input%total=strho
       CALL timestop("generation of start-density")
    END IF
    IF (mpi%irank == 0) THEN
       !
       !     --->generate spin polarized charge density
       !
       IF (input%swsp) THEN
          CALL timestart("optional: spin polarized density")
          CALL cdnsp(&
               &              atoms,input,vacuum,sphhar,&
               &              stars,sym,oneD,cell,dimension)
          !
          CALL timestop("optional: spin polarized density")
       END IF
       !
       !     --->flip magnetic moments
       !
       IF (input%lflip) THEN

          CALL timestart('optional: flip magnetic moments')
          CALL flipcdn(&
               &                atoms,input,vacuum,sphhar,&
               &                stars,sym,oneD,cell,noco%l_noco)
          !
          CALL timestop('optional: flip magnetic moments')
       END IF

       IF (obsolete%l_u2f) THEN

          CALL timestart('optional: conversion to formatted')
          CALL u2f(&
               &           stars,input,atoms,sphhar,vacuum,&
               &           cell,sym,noco%l_noco)
          !
          CALL timestop('optional: conversion to formatted')
       ENDIF

       IF (obsolete%l_f2u) THEN

          CALL timestart('optional: conversion to unformatted')
          CALL f2u(&
               &           stars,input,atoms,sphhar,vacuum,&
               &           cell,sym,noco%l_noco)
          !
          CALL timestop('optional: conversion to unformatted')
       ENDIF

       IF (input%l_bmt) THEN
          CALL bmt(&
               &           stars,input,noco,atoms,sphhar,vacuum,&
               &           cell,sym,oneD)
       ENDIF

    ENDIF ! mpi%irank == 0

  END SUBROUTINE OPTIONAL
END MODULE m_optional
