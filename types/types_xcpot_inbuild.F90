!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_xcpot_inbuild
   !This module contains the xcpot-type used for the in-build xc-implementations
   USE m_types_xcpot_data
   USE m_types_xcpot

   USE m_judft
   IMPLICIT NONE
   PRIVATE
   REAL, PARAMETER, PRIVATE :: hrtr_half = 0.5
   CHARACTER(len=4),PARAMETER:: xc_names(20)=[&
                                'l91 ','x-a ','wign','mjw ','hl  ','bh  ','vwn ','pz  ', &
                                'pw91','pbe ','rpbe','Rpbe','wc  ','PBEs', &
                                'pbe0','hse ','vhse','lhse','exx ','hf  ']

   LOGICAL,PARAMETER:: priv_gga(20)=[&
                       .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE.,.TRUE.]

   LOGICAL,PARAMETER:: priv_hybrid(20)=[&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.]

   REAL, PARAMETER       ::  amix_pbe0 = 0.25
   REAL, PARAMETER       ::  amix_hse  = 0.25
   REAL, PARAMETER       ::  amix_hf   = 1.00

   TYPE, EXTENDS(t_xcpot):: t_xcpot_inbuild
#ifdef CPP_MPI
      INTEGER             :: icorr=0 !not private to allow bcasting it around
#else
      INTEGER,PRIVATE     :: icorr=0
#endif

      TYPE(t_xcpot_data)   :: DATA

      LOGICAL,ALLOCATABLE :: lda_atom(:)

   CONTAINS
      !overloading t_xcpot:
      PROCEDURE        :: is_gga=>xcpot_is_gga
      PROCEDURE        :: is_hybrid=>xcpot_is_hybrid
      PROCEDURE        :: get_exchange_weight=>xcpot_get_exchange_weight
      PROCEDURE        :: get_vxc=>xcpot_get_vxc
      PROCEDURE        :: get_exc=>xcpot_get_exc
      !Overloading t_fleursetup
      PROCEDURE,PASS :: broadcast=>broadcast_xcpot
      PROCEDURE,PASS :: WRITE=>WRITE_xcpot
      PROCEDURE,PASS :: READ=>READ_xcpot
      PROCEDURE,PASS :: read_xml=>read_xml_xcpot
      !not overloaded
      PROCEDURE        :: get_name=>xcpot_get_name
      PROCEDURE        :: is_name=>xcpot_is_name
      PROCEDURE        :: init=>xcpot_init
   END TYPE t_xcpot_inbuild
   PUBLIC t_xcpot_inbuild
 CONTAINS
   SUBROUTINE broadcast_xcpot(tt,mpi_comm,origin)
#ifdef CPP_MPI
     USE m_bc_tool
#endif
     IMPLICIT NONE
     CLASS(t_xcpot_inbuild),INTENT(INOUT):: tt
    INTEGER,INTENT(IN)               :: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL      :: origin 
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER :: pe,ierr
    
    IF (PRESENT(origin)) THEN
       pe=origin
    ELSE
       pe=0
    ENDIF

    
   CALL MPI_BCAST(tt%icorr,1,MPI_INTEGER,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%krla,1,MPI_INTEGER,pe,mpi_comm,ierr)
   CALL mpi_bc(tt%lda_atom,pe,mpi_comm)
   CALL MPI_BCAST(tt%DATA%is_rpbe ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_wc,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_hse ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_pbes ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_pbe0 ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_bh ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%is_mjw ,1,MPI_LOGICAL,pe,mpi_comm,ierr)
   

   CALL MPI_BCAST(tt%DATA%uk,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%um,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
   CALL MPI_BCAST(tt%DATA%exchange_weight,1,MPI_DOUBLE_PRECISION,pe,mpi_comm,ierr)
#endif
 END SUBROUTINE broadcast_xcpot

 SUBROUTINE write_xcpot(tt, unit, iotype, v_list, iostat, iomsg)
   IMPLICIT NONE
   CLASS(t_xcpot_inbuild),INTENT(IN):: tt
  INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    WRITE(unit,*,IOSTAT=iostat) '"xcpot_inbuild":{'

    CALL json_print(unit,"icorr",tt%icorr)
    CALL json_print(unit,"DATA%krla",tt%DATA%krla)
    CALL json_print(unit,"lda_atom",tt%lda_atom)
    CALL json_print(unit,"DATA%is_rpbe ",tt%DATA%is_rpbe )
    CALL json_print(unit,"DATA%is_wc",tt%DATA%is_wc)
    CALL json_print(unit,"DATA%is_hse ",tt%DATA%is_hse )
    CALL json_print(unit,"DATA%is_pbes ",tt%DATA%is_pbes )
    CALL json_print(unit,"DATA%is_pbe0 ",tt%DATA%is_pbe0 )
    CALL json_print(unit,"DATA%is_bh ",tt%DATA%is_bh )
    CALL json_print(unit,"DATA%is_mjw ",tt%DATA%is_mjw )
    
    
    CALL json_print(unit,"DATA%uk",tt%DATA%uk)
    CALL json_print(unit,"DATA%um",tt%DATA%um)
    CALL json_print(unit,"DATA%exchange_weight",tt%DATA%exchange_weight)
    
    
    WRITE(unit,*,IOSTAT=iostat) '}'
    
  END SUBROUTINE write_xcpot

  SUBROUTINE read_xcpot(tt, unit, iotype, v_list, iostat, iomsg)
    IMPLICIT NONE
    CLASS(t_xcpot_inbuild),INTENT(INOUT):: tt
    INTEGER, INTENT(IN)             :: unit
    CHARACTER(*), INTENT(IN)        :: iotype
    INTEGER, INTENT(IN)             :: v_list(:)
    INTEGER, INTENT(OUT)            :: iostat
    CHARACTER(*), INTENT(INOUT)     :: iomsg

    CALL json_open_class("xcpot_inbuild",unit,iostat)
    IF (iostat.NE.0)   RETURN

    CALL json_read(unit,"icorr",tt%icorr)
    CALL json_read(unit,"DATA%krla",tt%DATA%krla)
    CALL json_read(unit,"lda_atom",tt%lda_atom)
    CALL json_read(unit,"DATA%is_rpbe ",tt%DATA%is_rpbe )
    CALL json_read(unit,"DATA%is_wc",tt%DATA%is_wc)
    CALL json_read(unit,"DATA%is_hse ",tt%DATA%is_hse )
    CALL json_read(unit,"DATA%is_pbes ",tt%DATA%is_pbes )
    CALL json_read(unit,"DATA%is_pbe0 ",tt%DATA%is_pbe0 )
    CALL json_read(unit,"DATA%is_bh ",tt%DATA%is_bh )
    CALL json_read(unit,"DATA%is_mjw ",tt%DATA%is_mjw )
    
    
    CALL json_read(unit,"DATA%uk",tt%DATA%uk)
    CALL json_read(unit,"DATA%um",tt%DATA%um)
    CALL json_read(unit,"DATA%exchange_weight",tt%DATA%exchange_weight)
    
    CALL json_close_class(unit,iostat)
    
  END SUBROUTINE read_xcpot

  SUBROUTINE read_xml_xcpot(tt)
    USE m_xmlIntWrapFort
    USE m_calculator
    USE m_inp_xml
    IMPLICIT NONE
    CLASS(t_xcpot_inbuild),INTENT(OUT):: tt


    CHARACTER(len=10):: namex
    LOGICAL          :: relcor
    CHARACTER(len=30):: xpath
    INTEGER          :: ntype,i
    IF (xmlGetNumberOfNodes('/fleurInput/xcFunctional/LibXCID')   == 0 &
         .AND. xmlGetNumberOfNodes('/fleurInput/xcFunctional/LibXCName') == 0) THEN
       IF (xmlGetNumberOfNodes('/fleurInput/calculationSetup/cutoffs/@GmaxXC')==1) THEN
          tt%gmaxxc= evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@GmaxXC'))
       ELSE
          tt%gmaxxc=  evaluateFirstOnly(xmlGetAttributeValue('/fleurInput/calculationSetup/cutoffs/@Gmax'))
       ENDIF
       relcor = evaluateFirstBoolOnly(xmlGetAttributeValue('/fleurInput/xcFunctional/@relativisticCorrections'))
       namex=TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL('/fleurInput/xcFunctional/@name')))))
       ntype=xmlGetNumberOfNodes('/fleurInput/atomGroups/atomGroup')
       CALL tt%init(namex,relcor,ntype)
       DO i=1,ntype
          xpath=inp_xml_speciesxpath_for_group(i)
          IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xpath))//'/special')==1) THEN
             tt%lda_atom(i)=evaluateFirstBoolOnly(TRIM(ADJUSTL(xmlGetAttributeValue(TRIM(ADJUSTL(xPath))//'/special/@lda'))))
          ELSE
             tt%lda_atom(i)=.FALSE.
          ENDIF
       END DO
       
    ENDIF
  END SUBROUTINE read_xml_xcpot

    
   CHARACTER(len=4) FUNCTION xcpot_get_name(xcpot)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild),INTENT(IN)    :: xcpot
      IF (xcpot%icorr==0) CALL judft_error("xc-potential not initialized",calledby="types_xcpot.F90")
      xcpot_get_name=xc_names(xcpot%icorr)
   END FUNCTION xcpot_get_name

   SUBROUTINE xcpot_init(xcpot,namex,relcor,ntype)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild),INTENT(INOUT)    :: xcpot
      CHARACTER(len=*),INTENT(IN)  :: namex
      LOGICAL,INTENT(IN)           :: relcor
      INTEGER,INTENT(IN)           :: ntype
      INTEGER:: n
      !Determine icorr from name

      ALLOCATE(xcpot%lda_atom(ntype))
      xcpot%lda_atom=.FALSE.
      xcpot%icorr=0
      DO n=1,SIZE(xc_names)
         IF (TRIM(ADJUSTL(namex))==TRIM(xc_names(n))) THEN
            xcpot%icorr=n
         ENDIF
      ENDDO
      if (xcpot%icorr==0) CALL judft_error("Unkown xc-potential:"//namex,calledby="types_xcpot.F90")
      xcpot%data%krla=MERGE(1,0,relcor)

      !Code from exchpbe to speed up determination of constants
      IF (xcpot%is_name("rpbe")) THEN
         xcpot%data%uk=1.2450
      ELSE
         xcpot%data%uk=0.8040
      ENDIF
      IF (xcpot%is_name("PBEs")) THEN     ! pbe_sol
         xcpot%data%um=0.123456790123456d0
      ELSE
         xcpot%data%um=0.2195149727645171e0
      ENDIF
      xcpot%data%is_hse=xcpot%is_name("hse").OR.xcpot%is_name("lhse").OR.xcpot%is_name("vhse")
      xcpot%data%is_rpbe=xcpot%is_name("Rpbe") !Rpbe
      xcpot%data%is_wc=xcpot%is_name("wc")
      xcpot%data%is_pbes=xcpot%is_name("PBEs")
      xcpot%data%is_pbe0=xcpot%is_name("pbe0")
      xcpot%data%is_mjw=xcpot%is_name("mjw")
      xcpot%data%is_bh=xcpot%is_name("bh")
      xcpot%DATA%exchange_weight=xcpot%get_exchange_weight()

   END SUBROUTINE xcpot_init

   LOGICAL FUNCTION xcpot_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild),INTENT(IN):: xcpot
      xcpot_is_gga=priv_gga(xcpot%icorr)
   END FUNCTION xcpot_is_gga

   LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild),INTENT(IN):: xcpot
      xcpot_is_hybrid=priv_hybrid(xcpot%icorr)
   END FUNCTION xcpot_is_hybrid

   FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_inbuild),INTENT(IN):: xcpot

      REAL:: a_ex

      a_ex=-1
      IF (xcpot%is_name("pbe0")) a_ex=amix_pbe0
      IF (xcpot%is_name("hf")) a_ex=amix_hf
      IF (xcpot%is_name("hse")) a_ex=amix_hse
      IF (xcpot%is_name("vhse")) a_ex=amix_hse
   END FUNCTION xcpot_get_exchange_weight


   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad)
!
      USE m_xcxal, ONLY : vxcxal
      USE m_xcwgn, ONLY : vxcwgn
      USE m_xcbh,  ONLY : vxcbh
      USE m_xcvwn, ONLY : vxcvwn
      USE m_xcpz,  ONLY : vxcpz
      USE m_vxcl91
      USE m_vxcwb91
      USE m_vxcpw91
      USE m_vxcepbe
      IMPLICIT NONE
!c
!c---> running mode parameters
!c
      CLASS(t_xcpot_inbuild),INTENT(IN) :: xcpot
      INTEGER, INTENT (IN)     :: jspins
!c
!c---> charge density
!c
      REAL,INTENT (IN) :: rh(:,:)
!c
!c---> xc potential
!c
      REAL, INTENT (OUT) :: vx (:,:)
      REAL, INTENT (OUT) :: vxc(:,:)

      ! optional arguments for GGA
      TYPE(t_gradients),INTENT(INOUT),OPTIONAL::grad
!c
!c ---> local scalars
      INTEGER :: ngrid
      REAL, PARAMETER :: hrtr_half = 0.5

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0,iofile=6
      REAL,PARAMETER      :: sprsv=0.0
      LOGICAL,PARAMETER   :: lwbc=.false. ! l-white-bird-current (ta)
!c
!c.....------------------------------------------------------------------
!c
!c-----> determine exchange correlation potential
!c
      vx (:,:) = 0.0
      vxc(:,:) = 0.0
      ngrid=SIZE(rh,1)

      IF (xcpot%is_gga()) THEN
         IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
         IF (xcpot%is_name("l91")) THEN    ! local pw91
            CALL vxcl91(jspins,ngrid,ngrid,rh,grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid), grad%g2rt(:ngrid),&
                 grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),grad%gggrd(:ngrid),&
                 grad%gzgr(:ngrid), vx(:ngrid,:),vxc(:ngrid,:), isprsv,sprsv)
         ELSEIF (xcpot%is_name("pw91")) THEN  ! pw91
            IF (lwbc) THEN
               CALL vxcwb91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid),&
                 grad%g2rt(:ngrid),grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),&
                 grad%gggrd(:ngrid),grad%gzgr(:ngrid), vx(:ngrid,:),vxc(:ngrid,:), idsprs,isprsv,sprsv)
            ELSE

               CALL vxcpw91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt(:ngrid),grad%agru(:ngrid),grad%agrd(:ngrid),&
                 grad%g2rt(:ngrid),grad%g2ru(:ngrid),grad%g2rd(:ngrid),grad%gggrt(:ngrid),grad%gggru(:ngrid),&
                 grad%gggrd,grad%gzgr, vx(:ngrid,:),vxc(:ngrid,:), idsprs,isprsv,sprsv)

            ENDIF
         ELSE  ! pbe or similar
            CALL vxcepbe(xcpot%DATA,jspins,ngrid,ngrid,rh(:ngrid,:), grad%agrt,grad%agru,grad%agrd,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd, vx(:ngrid,:),vxc(:ngrid,:))
         ENDIF
      ELSE  !LDA potentials
         IF (xcpot%is_name("x-a"))  THEN   ! X-alpha method
            CALL vxcxal(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("wign")) THEN    ! Wigner interpolation formula
            CALL vxcwgn(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("mjw").OR.xcpot%is_name("bh")) THEN ! von Barth,Hedin correlation
            CALL vxcbh(iofile,xcpot%data,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))

         ELSEIF (xcpot%is_name("vwn")) THEN     ! Vosko,Wilk,Nusair correlation
            CALL vxcvwn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("pz")) THEN     ! Perdew,Zunger correlation
            CALL vxcpz(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
         ELSEIF (xcpot%is_name("hf")) THEN
            ! Hartree-Fock  calculation: X-alpha potential is added to generate a rational local potential,
            !                            later it is subtracted again
            CALL juDFT_error('HF should now be treated as a GGA functional', calledby='xcpot_get_vxc')
            CALL vxcxal(xcpot%data%krla,jspins, ngrid,ngrid,rh(:ngrid,:), vx(:ngrid,:),vxc(:ngrid,:))
            !         vxc=0
         ELSEIF (xcpot%is_name("exx")) THEN
            ! if exact exchange calculation do nothing
            vxc = 0
         ELSE
            CALL juDFT_error("Unkown LDA potential",calledby="type xcpot")
         ENDIF
      ENDIF
!
!-----> hartree units
!
      vx=hrtr_half*vx
      vxc=hrtr_half*vxc

   END SUBROUTINE xcpot_get_vxc

!***********************************************************************
   SUBROUTINE xcpot_get_exc(xcpot,jspins,rh, exc,grad)
!***********************************************************************
      USE m_xcxal, ONLY : excxal
      USE m_xcwgn, ONLY : excwgn
      USE m_xcbh,  ONLY : excbh
      USE m_xcvwn, ONLY : excvwn
      USE m_xcpz,  ONLY : excpz
      USE m_excl91
      USE m_excwb91
      USE m_excpw91
      USE m_excepbe
      IMPLICIT NONE
!c
!c---> running mode parameters
!c
      CLASS(t_xcpot_inbuild),INTENT(IN) :: xcpot
      INTEGER, INTENT (IN)     :: jspins
!c
!c---> charge density
!c
      REAL,INTENT (IN) :: rh(:,:)
!c
!c---> xc energy density
!c
      REAL, INTENT (OUT) :: exc(:)

      ! optional arguments for GGA
      TYPE(t_gradients),OPTIONAL,INTENT(IN)::grad
!c
!c ---> local scalars
      INTEGER :: ngrid
      REAL, PARAMETER :: hrtr_half = 0.5

      !used to be dummy arguments for testing
      INTEGER,PARAMETER   :: idsprs=0,isprsv=0,iofile=6
      REAL,PARAMETER      :: sprsv=0.0
      LOGICAL,PARAMETER   :: lwbc=.false. ! l-white-bird-current (ta)
!c
!c-----> determine exchange correlation energy density
!c
      exc(:) = 0.0
      ngrid=SIZE(rh,1)
      IF (xcpot%is_gga()) THEN
         IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_exc for a GGA potential without providing derivatives")
         IF (xcpot%is_name("l91")) THEN  ! local pw91
            CALL excl91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd,grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, isprsv,sprsv)
         ELSEIF (xcpot%is_name("pw91")) THEN     ! pw91
            IF (lwbc) THEN
               CALL excwb91(ngrid,ngrid,rh(:ngrid,1),rh(:ngrid,2),grad%agrt,grad%agru,grad%agrd, grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, idsprs,isprsv,sprsv)
            ELSE
               CALL excpw91(jspins,ngrid,ngrid,rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd, grad%g2rt,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd,grad%gzgr, exc, idsprs,isprsv,sprsv)
            ENDIF
         ELSE
            CALL excepbe(xcpot%data,jspins,ngrid,ngrid, rh(:ngrid,:),grad%agrt,grad%agru,grad%agrd,grad%g2ru,grad%g2rd,grad%gggrt,grad%gggru,grad%gggrd, exc)
         ENDIF
      ELSE !LDA
         IF (xcpot%is_name("x-a"))  THEN   ! X-alpha method
            CALL excxal(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("wign")) THEN    ! Wigner interpolation formula
            CALL excwgn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("mjw").OR.xcpot%is_name("bh")) THEN ! von Barth,Hedin correlation
            CALL excbh(iofile,xcpot%data,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("vwn")) THEN     ! Vosko,Wilk,Nusair correlation
            CALL excvwn(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("pz")) THEN     ! Perdew,Zunger correlation
            CALL excpz(iofile,xcpot%data%krla,jspins, ngrid,ngrid,rh, exc)
         ELSEIF (xcpot%is_name("hf") .OR. xcpot%is_name("exx")) THEN
            CALL juDFT_error('HF should now be treated as a GGA functional', calledby='xcpot_get_exc')
            exc=0
         ELSE
            CALL juDFT_error("Unkown LDA potential",calledby="type xcpot")
         ENDIF
      ENDIF
!c-----> hartree units
      exc= hrtr_half*exc

   END SUBROUTINE xcpot_get_exc

   LOGICAL FUNCTION xcpot_is_name(xcpot,name)
      CLASS(t_xcpot_inbuild),INTENT(IN):: xcpot
      CHARACTER(len=*),INTENT(IN)  :: name
      xcpot_is_name=(TRIM(xc_names(xcpot%icorr))==TRIM((name)))
   END FUNCTION xcpot_is_name

END MODULE m_types_xcpot_inbuild
