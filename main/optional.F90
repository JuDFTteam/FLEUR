!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_optional
  USE m_juDFT
#ifdef CPP_MPI 
  use mpi 
#endif
CONTAINS
  SUBROUTINE OPTIONAL(fmpi, atoms,sphhar,vacuum,&
       stars,input,sym, cell, sliceplot, xcpot, noco)
    !
    !----------------------------------------
    ! this routine is called by: fleur.F90
    !
    ! optional  stden -+- atom2 -+- setcor
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
    !            |                       +- starf2 -- spgrot
    !            |                       +- fitchk
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
    USE m_stden
    USE m_cdnsp
    USE m_flipcdn
    USE m_cdn_io
    USE m_types


    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..

    TYPE(t_mpi),INTENT(IN)      :: fmpi
    TYPE(t_atoms),INTENT(IN)    :: atoms

    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_stars),INTENT(IN)    :: stars
     
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_cell),INTENT(IN)     :: cell
    CLASS(t_xcpot),INTENT(IN)   :: xcpot
    TYPE(t_sliceplot),INTENT(IN):: sliceplot
    !     ..
    !     .. Local Scalars ..
    INTEGER :: it, archiveType, atomsCounter
    CHARACTER*10 :: cdnfname
    LOGICAL :: strho
    LOGICAL :: stateCheck=.TRUE.
#ifdef CPP_MPI
    INTEGER :: ierr
#endif
    !     ..
    it = 1


    !
    !     --->generate starting charge density
    !
    strho=input%strho
    IF (.NOT.(strho.OR.(sliceplot%iplot.NE.0))) THEN
       archiveType = CDN_ARCHIVE_TYPE_CDN1_const
       IF (any(noco%l_unrestrictMT)) THEN
          archiveType = CDN_ARCHIVE_TYPE_FFN_const
       ELSE IF (noco%l_noco) THEN
          archiveType = CDN_ARCHIVE_TYPE_NOCO_const
       END IF
       IF (fmpi%irank == 0) THEN
          strho = .NOT.isDensityFilePresent(archiveType)
       END IF
#ifdef CPP_MPI
       CALL MPI_BCAST(strho,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
#endif
    ENDIF
    IF (strho) THEN
       !strho=input%total
       !input%total = .FALSE.
       !
       CALL timestart("generation of start-density")
       IF (input%jspins.EQ.2) THEN
          DO atomsCounter=1, atoms%ntype
             IF(.NOT.MAXVAL(ABS(atoms%econf(atomsCounter)%Occupation(:,1)-atoms%econf(atomsCounter)%Occupation(:,2))).EQ.0)stateCheck=.FALSE.
          END DO
       END IF
       IF (stateCheck.AND.(input%jspins.EQ.2)) CALL juDFT_warn("You're setting up a spin-polarized calculation (jspins=2) without any actual polarization given in the systems occupation. You're sure you want that?", calledby = "optional")
       CALL stden(fmpi,sphhar,stars,atoms,sym,vacuum,&
                  input,cell,xcpot,noco )
       !
       !input%total=strho
       CALL timestop("generation of start-density")
    END IF
    IF (fmpi%irank == 0) THEN
       !
       !     --->generate spin polarized charge density
       !
       IF (input%swsp) THEN
          CALL timestart("optional: spin polarized density")
          CALL cdnsp(atoms,input,vacuum,sphhar,stars,sym,noco ,cell)
          !
          CALL timestop("optional: spin polarized density")
       END IF
       !
       !     --->flip magnetic moments
       !
       IF (input%lflip) THEN

          CALL timestart('optional: flip magnetic moments')
          CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco ,cell,toGlobal=.true.)
          print *,"TODO,toGlobal in optional should be removed"
          CALL timestop('optional: flip magnetic moments')
       END IF



       IF (input%l_bmt) THEN
          CALL bmt(stars,input,noco,atoms,sphhar,vacuum,cell,sym )
       ENDIF

    ENDIF ! fmpi%irank == 0

    IF (input%strho)          CALL juDFT_end("starting density generated",fmpi%irank)
    IF (input%swsp)           CALL juDFT_end("spin polarised density generated",fmpi%irank)
    IF (input%lflip)          CALL juDFT_end("magnetic moments flipped",fmpi%irank)
    IF (input%l_bmt)          CALL juDFT_end('"cdnbmt" written',fmpi%irank)

  END SUBROUTINE OPTIONAL
END MODULE m_optional
