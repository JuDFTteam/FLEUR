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
       
  INTEGER,PARAMETER:: ILLEGAL_XCPOT=0
  TYPE t_xcpot
     INTEGER,PRIVATE     :: icorr
     INTEGER,ALLOCATABLE :: icorr_mt(:)
     REAL                :: gmaxxc
     INTEGER             :: krla !relativistic corrections
     
   CONTAINS
     PROCEDURE        :: is_gga=>xcpot_is_gga
     PROCEDURE,NOPASS :: from_name=>xcpot_from_name
     PROCEDURE        :: init=>xcpot_init
     PROCEDURE        :: is_hybrid=>xcpot_is_hybrid
     PROCEDURE        :: is_name=>xcpot_is_name
     PROCEDURE        :: get_exchange_weight=>xcpot_get_exchange_weight
  END TYPE t_xcpot
  PUBLIC t_xcpot,ILLEGAL_XCPOT
CONTAINS
  PURE INTEGER FUNCTION xcpot_from_name(name)
    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)  :: name
    INTEGER :: n
    xcpot_from_name=ILLEGAL_XCPOT
    DO n=1,SIZE(xc_names)
       IF (TRIM(ADJUSTL(name))==TRIM(xc_names(n))) xcpot_from_name=n
    ENDDO
  END FUNCTION xcpot_from_name

  SUBROUTINE xcpot_init(xcpot,namex,relcor)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(INOUT)    :: xcpot
    CHARACTER(len=*),INTENT(IN)  :: namex
    LOGICAL,INTENT(IN)           :: relcor

    xcpot%icorr=xcpot%from_name(namex)
    xcpot%krla=MERGE(1,0,relcor)
  END SUBROUTINE xcpot_init
  
  PURE LOGICAL FUNCTION xcpot_is_name(xcpot,name)
    CLASS(t_xcpot),INTENT(IN):: xcpot
    CHARACTER(len=*),INTENT(IN)  :: name
    xcpot_is_name=(xcpot%icorr/=ILLEGAL_XCPOT.AND.xcpot%icorr==xcpot_from_name(name))
  END FUNCTION xcpot_is_name

  ELEMENTAL LOGICAL FUNCTION xcpot_is_gga(xcpot,icorr)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    INTEGER,OPTIONAL,INTENT(IN):: icorr
    IF (PRESENT(icorr)) THEN
       xcpot_is_gga=priv_gga(icorr)
    ELSE
       xcpot_is_gga=priv_gga(xcpot%icorr)
    ENDIF
  END FUNCTION xcpot_is_gga

  ELEMENTAL LOGICAL FUNCTION xcpot_is_hybrid(xcpot,icorr)
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
