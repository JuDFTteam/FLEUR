!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_optional
  USE m_juDFT
CONTAINS
  SUBROUTINE OPTIONAL(mpi, atoms,sphhar,vacuum,&
       stars,input,sym, cell, sliceplot, xcpot, noco, oneD)
    !
    !----------------------------------------
    ! this routine is called by: fleur.F90
    !
    ! optional --+-- plot -+- loddop
    !            |         +- outcdn -+- starf2 -- spgrot
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

    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_atoms),INTENT(IN)    :: atoms

    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(INOUT) :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_cell),INTENT(IN)     :: cell
    CLASS(t_xcpot),INTENT(IN)   :: xcpot
    TYPE(t_sliceplot),INTENT(IN):: sliceplot
    !     ..
    !     .. Local Scalars ..
    INTEGER :: it, archiveType
    CHARACTER*10 :: cdnfname
    LOGICAL :: strho
#ifdef CPP_MPI
    include 'mpif.h'
    INTEGER :: ierr(2)
#endif
    !     ..
    it = 1

 !   IF ((sliceplot%iplot.NE.0 ).AND. (mpi%irank==0) ) THEN
 !      IF (noco%l_noco) THEN
 !         CALL pldngen(mpi,sym,stars,atoms,sphhar,vacuum,&
 !              cell,input,noco,oneD,sliceplot)
 !      ENDIF
 !   ENDIF


    !
    !     --->generate starting charge density
    !
    strho=input%strho
    IF (.NOT.(strho.OR.(sliceplot%iplot.NE.0))) THEN
       archiveType = CDN_ARCHIVE_TYPE_CDN1_const
       IF (noco%l_noco) THEN
          archiveType = CDN_ARCHIVE_TYPE_NOCO_const
       END IF
       IF (mpi%irank == 0) THEN
          strho = .NOT.isDensityFilePresent(archiveType)
       END IF
#ifdef CPP_MPI
       CALL MPI_BCAST(strho,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
    ENDIF
    IF (strho) THEN
       strho=input%total
       input%total = .FALSE.
       !
       CALL timestart("generation of start-density")
       CALL stden(mpi,sphhar,stars,atoms,sym,vacuum,&
                  input,cell,xcpot,noco,oneD)
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
          CALL cdnsp(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell)
          !
          CALL timestop("optional: spin polarized density")
       END IF
       !
       !     --->flip magnetic moments
       !
       IF (input%lflip) THEN

          CALL timestart('optional: flip magnetic moments')
          CALL flipcdn(atoms,input,vacuum,sphhar,stars,sym,noco,oneD,cell)
          !
          CALL timestop('optional: flip magnetic moments')
       END IF



       IF (input%l_bmt) THEN
          CALL bmt(stars,input,noco,atoms,sphhar,vacuum,cell,sym,oneD)
       ENDIF

    ENDIF ! mpi%irank == 0

    IF (input%strho)          CALL juDFT_end("starting density generated",mpi%irank)
    IF (input%swsp)           CALL juDFT_end("spin polarised density generated",mpi%irank)
    IF (input%lflip)          CALL juDFT_end("magnetic moments flipped",mpi%irank)
    IF (input%l_bmt)          CALL juDFT_end('"cdnbmt" written',mpi%irank)

  END SUBROUTINE OPTIONAL
END MODULE m_optional
