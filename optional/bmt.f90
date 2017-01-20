MODULE m_bmt
contains
  SUBROUTINE bmt(&
       & stars,input,noco,atoms,sphhar,vacuum,&
       & cell,sym,oneD)
    !
    use m_types
    use m_juDFT
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
    INTEGER k,i,ivac  ,it 
    INTEGER type,typmag, archiveType
    CHARACTER(len=8) filename 
    COMPLEX, ALLOCATABLE :: fpw(:,:),fzxy(:,:,:,:)
    REAL,    ALLOCATABLE :: fz(:,:,:),fr(:,:,:,:)
    COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)

    !     ..
    !     ..


    typmag= atoms%ntype 
    ! only muffin-tins with type <= typmag remain magnetic  


    IF (input%jspins/=2) THEN
       CALL juDFT_error("Stop in bmt:  jspins/=2",calledby="bmt")
    ENDIF

    !atoms%jmtd = maxval(atoms%jri(:))
    !sphhar%nlhd = maxval(sphhar%nlh(:))

    ALLOCATE(fpw(stars%ng3,input%jspins),fzxy(vacuum%nmzxy,stars%ng2-1,2,input%jspins))
    ALLOCATE(fr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),fz(vacuum%nmz,2,input%jspins))
    ALLOCATE(cdom(stars%n3d),cdomvz(vacuum%nmzd,2),cdomvxy(vacuum%nmzxyd,oneD%odi%n2d-1,2))

    archiveType = CDN_ARCHIVE_TYPE_CDN1_const
    IF (noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_NOCO_const

    CALL readDensity(stars,vacuum,atoms,sphhar,input,sym,oneD,archiveType,&
                     CDN_INPUT_DEN_const,0,it,fr,fpw,fz,fzxy,cdom,cdomvz,cdomvxy)

    IF ( typmag < atoms%ntype ) THEN 
       DO type= typmag+1,atoms%ntype 
          DO k= 0,sphhar%nlhd
             DO i= 1,atoms%jmtd
                fr(i,k,type,1)= ( fr(i,k,type,1) + fr(i,k,type,2) )/2. 
                fr(i,k,type,2)= fr(i,k,type,1) 
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    DO k= 1,stars%ng3
       fpw(k,1)= ( fpw(k,1) + fpw(k,2) )/2.
       fpw(k,2)= fpw(k,1)
    ENDDO
    IF (input%film) THEN
       DO ivac= 1,vacuum%nvac
          DO i= 1,vacuum%nmz
             fz(i,ivac,1)= ( fz(i,ivac,1) + fz(i,ivac,2) )/2.  
             fz(i,ivac,2)= fz(i,ivac,1)
          ENDDO
          DO k= 2,stars%ng2
             DO i= 1,vacuum%nmzxy
                fzxy(i,k-1,ivac,1)= &
                     &         ( fzxy(i,k-1,ivac,1) + fzxy(i,k-1,ivac,2) )/2. 
                fzxy(i,k-1,ivac,2)= fzxy(i,k-1,ivac,1)
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
         & it,fr,fpw,fz,fzxy)
    CLOSE(98) 

    DEALLOCATE(cdom,cdomvz,cdomvxy)
    DEALLOCATE(fpw,fzxy,fr,fz)

  END SUBROUTINE bmt
END MODULE m_bmt
