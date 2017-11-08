MODULE m_bmt
contains
  SUBROUTINE bmt(&
       & stars,input,noco,atoms,sphhar,vacuum,&
       & cell,sym,oneD)

    USE m_constants
    USE m_types
    USE m_juDFT
    USE m_cdn_io
    USE m_wrtdop
    IMPLICIT NONE
    !     ..
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_input),INTENT(INOUT) :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_oneD),INTENT(IN)     :: oneD

    TYPE(t_potden) :: den

    INTEGER k,i,ivac
    INTEGER type,typmag, archiveType
    REAL fermiEnergyTemp
    LOGICAL l_qfix
    CHARACTER(len=8) filename

    typmag= atoms%ntype 
    ! only muffin-tins with type <= typmag remain magnetic  


    IF (input%jspins/=2) THEN
       CALL juDFT_error("Stop in bmt:  jspins/=2",calledby="bmt")
    ENDIF

    !atoms%jmtd = maxval(atoms%jri(:))
    !sphhar%nlhd = maxval(sphhar%nlh(:))

    CALL den%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,.FALSE.,POTDEN_TYPE_DEN)
    IF(noco%l_noco) THEN
       archiveType = CDN_ARCHIVE_TYPE_NOCO_const
    ELSE
       archiveType = CDN_ARCHIVE_TYPE_CDN1_const
    END IF

    CALL readDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,&
                     CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den)

    IF ( typmag < atoms%ntype ) THEN 
       DO type= typmag+1,atoms%ntype 
          DO k= 0,sphhar%nlhd
             DO i= 1,atoms%jmtd
                den%mt(i,k,type,1)= (den%mt(i,k,type,1) + den%mt(i,k,type,2))/2. 
                den%mt(i,k,type,2)= den%mt(i,k,type,1) 
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    DO k= 1,stars%ng3
       den%pw(k,1)= (den%pw(k,1) + den%pw(k,2))/2.0
       den%pw(k,2)= den%pw(k,1)
    ENDDO
    IF (input%film) THEN
       DO ivac= 1,vacuum%nvac
          DO i= 1,vacuum%nmz
             den%vacz(i,ivac,1)= (den%vacz(i,ivac,1) + den%vacz(i,ivac,2))/2.0
             den%vacz(i,ivac,2)= den%vacz(i,ivac,1)
          ENDDO
          DO k= 2,stars%ng2
             DO i= 1,vacuum%nmzxy
                den%vacxy(i,k-1,ivac,1)= (den%vacxy(i,k-1,ivac,1) + den%vacxy(i,k-1,ivac,2))/2.0
                den%vacxy(i,k-1,ivac,2)= den%vacxy(i,k-1,ivac,1)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    filename= 'cdnbmtXX'
    IF ( typmag < atoms%ntype ) THEN
       filename(7:7)= ACHAR(IACHAR('0')+ MOD(typmag,100)/10 )
       filename(8:8)= ACHAR(IACHAR('0')+ MOD(typmag,10) )
       OPEN(98,file=filename(1:8),form='unformatted',status='replace')
    ELSE
       OPEN(98,file=filename(1:6),form='unformatted',status='replace')
    ENDIF
    CALL wrtdop(&
         & stars,vacuum,atoms,sphhar,input,sym,&
         & 98,&
         & den%iter,den%mt,den%pw,den%vacz,den%vacxy)
    CLOSE(98) 

  END SUBROUTINE bmt
END MODULE m_bmt
