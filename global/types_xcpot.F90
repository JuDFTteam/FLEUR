MODULE m_types_xcpot
  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=4),PARAMETER:: xc_names(19)=[&
       'l91 ','x-a ','wign','mjw ','hl  ','bh  ','vwn ','pz  ', &  
       'pw91','pbe ','rpbe','Rpbe','wc  ','PBEs', & 
       'hse ','vhse','lhse','exx ','hf  '] 
  
  LOGICAL,PARAMETER:: priv_gga(19)=[&
       .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,&
       .TRUE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.]

  LOGICAL,PARAMETER:: priv_hybrid(19)=[&
       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
       .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
       .TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.]

  REAL, PARAMETER       ::  amix_pbe0 = 0.25
  REAL, PARAMETER       ::  amix_hse  = 0.25
  REAL, PARAMETER       ::  amix_hf   = 1.00
       
  TYPE t_xcpot
#ifdef CPP_MPI     
     INTEGER             :: icorr=0 !not private to allow bcasting it around
#else
     INTEGER,PRIVATE     :: icorr=0
#endif

     !in the pbe case (exchpbe.F) lots of test are made
     !in addition some constants are set
     !to speed up this code precalculate things in init
     LOGICAL             :: is_rpbe !Rpbe
     LOGICAL             :: is_wc
     LOGICAL             :: is_hse !hse,lhse,vhse
     REAL                :: uk,um
     
     LOGICAL,ALLOCATABLE :: lda_atom(:)
     REAL                :: gmaxxc
     INTEGER             :: krla !relativistic corrections
     
   CONTAINS
     PROCEDURE        :: is_gga=>xcpot_is_gga
     PROCEDURE        :: get_name=>xcpot_get_name
     PROCEDURE        :: init=>xcpot_init
     PROCEDURE        :: is_hybrid=>xcpot_is_hybrid 
     PROCEDURE        :: is_name=>xcpot_is_name
     PROCEDURE        :: get_exchange_weight=>xcpot_get_exchange_weight
  END TYPE t_xcpot
  PUBLIC t_xcpot
CONTAINS
  CHARACTER(len=4) FUNCTION xcpot_get_name(xcpot)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN)    :: xcpot
    IF (xcpot%icorr==0) CALL judft_error("xc-potential not initialized",calledby="types_xcpot.F90")
    xcpot_get_name=xc_names(xcpot%icorr)
  END FUNCTION xcpot_get_name

  SUBROUTINE xcpot_init(xcpot,namex,relcor)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(INOUT)    :: xcpot
    CHARACTER(len=*),INTENT(IN)  :: namex
    LOGICAL,INTENT(IN)           :: relcor
    INTEGER:: n
    !Determine icorr from name
    xcpot%icorr=0
    DO n=1,SIZE(xc_names)
       IF (TRIM(ADJUSTL(namex))==TRIM(xc_names(n))) THEN
          xcpot%icorr=n
       ENDIF
    ENDDO
    if (xcpot%icorr==0) CALL judft_error("Unkown xc-potential:"//namex,calledby="types_xcpot.F90")
    xcpot%krla=MERGE(1,0,relcor)

    
    !Code from exchpbe to speed up determination of constants
    IF (xcpot%is_name("rpbe")) THEN
       xcpot%uk=1.2450
    ELSE
       xcpot%uk=0.8040
    ENDIF
    IF (xcpot%is_name("PBEs")) THEN     ! pbe_sol
       xcpot%um=0.123456790123456d0
    ELSE
       xcpot%um=0.2195149727645171e0
    ENDIF
    xcpot%is_hse=xcpot%is_name("hse").OR.xcpot%is_name("lhse").OR.xcpot%is_name("vhse")
    xcpot%is_rpbe=xcpot%is_name("Rpbe") !Rpbe
    xcpot%is_wc=xcpot%is_name("wc")
      
      
      
  END SUBROUTINE xcpot_init
  
  LOGICAL FUNCTION xcpot_is_name(xcpot,name)
    CLASS(t_xcpot),INTENT(IN):: xcpot
    CHARACTER(len=*),INTENT(IN)  :: name
    xcpot_is_name=(trim(xc_names(xcpot%icorr))==trim((name)))
  END FUNCTION xcpot_is_name

  LOGICAL FUNCTION xcpot_is_gga(xcpot,icorr)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    INTEGER,OPTIONAL,INTENT(IN):: icorr
    IF (PRESENT(icorr)) THEN
       xcpot_is_gga=priv_gga(icorr)
    ELSE
       xcpot_is_gga=priv_gga(xcpot%icorr)
    ENDIF
  END FUNCTION xcpot_is_gga

  LOGICAL FUNCTION xcpot_is_hybrid(xcpot,icorr)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    INTEGER,OPTIONAL,INTENT(IN):: icorr
    IF (PRESENT(icorr)) THEN
       xcpot_is_hybrid=priv_hybrid(icorr)
    ELSE
       xcpot_is_hybrid=priv_hybrid(xcpot%icorr)
    ENDIF
  END FUNCTION xcpot_is_hybrid

  FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    
    REAL:: a_ex

    a_ex=-1
    IF (xcpot%is_name("pbe0")) a_ex=amix_pbe0
    IF (xcpot%is_name("hf")) a_ex=amix_hf
    IF (xcpot%is_name("hse")) a_ex=amix_hse
    IF (xcpot%is_name("vhse")) a_ex=amix_hse

    IF (a_ex==-1) CALL judft_error('xc functional can not be identified')
  END FUNCTION xcpot_get_exchange_weight
END MODULE m_types_xcpot
