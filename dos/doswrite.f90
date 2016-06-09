MODULE m_doswrite
  USE m_juDFT
  !
  !-- now write cdninf for all kpts if on T3E
  !-- now read data from tmp_dos and write to vacdos&dosinp .. dw
  !
CONTAINS
  SUBROUTINE doswrite(&
       &                   eig_id,DIMENSION,kpts,atoms,vacuum,&
       &                   input,banddos,&
       &                   sliceplot,noco,sym,&
       &                   cell,&
       &                   l_mcd,ncored,ncore,e_mcd,&
       &                   efermi,nsld,oneD)
    USE m_eig66_io,ONLY:read_dos,read_eig
    USE m_evaldos
    USE m_cdninf
    USE m_types
    IMPLICIT NONE
  
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_banddos),INTENT(IN)   :: banddos
    TYPE(t_sliceplot),INTENT(IN) :: sliceplot
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_atoms),INTENT(IN)     :: atoms

    !     .. Scalar Arguments ..
    INTEGER,PARAMETER :: n2max=13 
    INTEGER, INTENT (IN) :: nsld,eig_id 
    INTEGER, INTENT (IN) :: ncored
    REAL,    INTENT (IN) :: efermi
    LOGICAL, INTENT (IN) :: l_mcd
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN)  :: ncore(atoms%ntypd)
    REAL, INTENT(IN)      :: e_mcd(atoms%ntypd,input%jspins,ncored)
    !-odim
    !+odim

    !    locals
    INTEGER :: jsym(DIMENSION%neigd),ksym(DIMENSION%neigd)
    REAL    :: wk,bkpt(3)
    REAL   :: eig(DIMENSION%neigd)
    REAL   :: qal(0:3,atoms%ntypd,DIMENSION%neigd,DIMENSION%jspd)
    REAL   :: qis(DIMENSION%neigd,kpts%nkptd,DIMENSION%jspd)
    REAL   :: qvac(DIMENSION%neigd,2,kpts%nkptd,DIMENSION%jspd)
    REAL   :: qvlay(DIMENSION%neigd,vacuum%layerd,2)
    COMPLEX :: qstars(vacuum%nstars,DIMENSION%neigd,vacuum%layerd,2)
    INTEGER :: ne,ikpt,kspin,j,i,n
    COMPLEX, ALLOCATABLE :: ac(:,:),bc(:,:)


    !     check if there is anything todo here
    IF (.NOT.(banddos%dos.OR.input%cdinf.OR.banddos%vacdos.OR.(vacuum%nstm.EQ.3))) RETURN
    !     check if settings in inp-file make any sense
    IF (banddos%vacdos.AND..NOT.banddos%dos) THEN
       WRITE(6,*) "STOP DOS: only set banddos%vacdos = .true. if banddos%dos=.true."
       CALL juDFT_error("DOS",calledby ="doswrite")
    ENDIF
    IF (banddos%vacdos.AND.(.NOT.vacuum%starcoeff.AND.(vacuum%nstars.NE.1)))THEN
       WRITE(6,*) "STOP DOS: if stars = f set vacuum%nstars=1"
       CALL juDFT_error("DOS",calledby ="doswrite")
    ENDIF

    IF (banddos%dos.AND.(banddos%ndir.GE.0)) THEN
       !---  >    open files for bandstucture+ old style vacdos
       OPEN (85,file='dosinp')
       IF (banddos%vacdos) THEN
          OPEN (86,file='vacDOS')
       ENDIF
    ENDIF

    IF ((banddos%dos.AND.(banddos%ndir.GE.0)).OR.input%cdinf) THEN
       !
       !      write bandstructure or cdn-info to output-file
       !         
       DO kspin = 1,input%jspins
          IF (banddos%dos.AND.(banddos%ndir.GE.0)) THEN
             !---  >       write header information to vacdos & dosinp
             !     
             IF (input%film) THEN
                WRITE (85,FMT=8080) vacuum%nvac,kpts%nkpt
             ELSE
                WRITE (85,FMT=8080) input%jspins,kpts%nkpt
             ENDIF
8080         FORMAT (12i6)
             WRITE (85,FMT=8080) atoms%ntype, (atoms%neq(n),n=1,atoms%ntype)
             IF (banddos%vacdos) THEN
                WRITE (86,FMT=8080) vacuum%nvac,kpts%nkpt
                WRITE (86,FMT=8080) vacuum%layers
                WRITE (86,'(20(i3,1x))') (vacuum%izlay(i,1),i=1,vacuum%layers)
             ENDIF
          ENDIF

          DO ikpt=1,kpts%nkpt
             call read_eig(eig_id,ikpt,kspin,&
                                 bk=bkpt,wk=wk,neig=ne,eig=eig)
             call read_dos(eig_id,ikpt,kspin,&
                  &              qal(:,:,:,kspin),qvac(:,:,ikpt,kspin),&
                  &              qis(:,ikpt,kspin),&
                  &              qvlay(:,:,:),qstars,ksym,jsym)

             CALL cdninf(&
                  &              input,sym,noco,kspin,atoms,&
                  &              vacuum,sliceplot,banddos,ikpt,bkpt,&
                  &              wk,cell,kpts,&
                  &              ne,eig,qal(0:,:,:,kspin),qis,qvac,&
                  &              qvlay(:,:,:),&
                  &              qstars,ksym,jsym)
          ENDDO

       ENDDO                  ! end spin loop (kspin = 1,input%jspins)

    ENDIF

    IF (banddos%dos.AND.(banddos%ndir.GE.0)) THEN
       CLOSE(85)
       RETURN
       !     ok, all done in the bandstructure/cdninf case
    ENDIF
    !
    !     write DOS/VACDOS
    !
    !     
    IF (banddos%dos.AND.(banddos%ndir.LT.0)) THEN
       CALL evaldos(&
            &                eig_id,input,banddos,vacuum,kpts,atoms,sym,noco,oneD,cell,&
            &                DIMENSION,efermi,&
            &                l_mcd,ncored,ncore,e_mcd,nsld)
    ENDIF
    !
    !      Now write to vacwave if nstm=3 
    !     all data
    !     has been written to tmp_vacwave and must be written now
    !     by PE=0 only!
    !
    IF (vacuum%nstm.EQ.3) THEN
       call juDFT_error("nstm=3 not implemented in doswrite")
       OPEN (89,file='tmp_vacwave',status='old',access='direct')!, recl=reclength_vw)
       ALLOCATE ( ac(n2max,DIMENSION%neigd),bc(n2max,DIMENSION%neigd) )
       DO ikpt = 1,kpts%nkpt
          WRITE(*,*) 'Read rec',ikpt,'from vacwave'
          READ(89,rec=ikpt) wk,ne,bkpt(1),bkpt(2),&
               &                 eig,ac,bc
          WRITE (87,'(i3,1x,f12.6)') ikpt,wk
          i=0
          DO n = 1, ne
             IF (ABS(eig(n)-vacuum%tworkf).LE.banddos%e2_dos) i=i+1
          END DO
          WRITE (87,FMT=990) bkpt(1), bkpt(2), i, n2max
          DO n = 1, ne
             IF (ABS(eig(n)-vacuum%tworkf).LE.banddos%e2_dos) THEN
                WRITE (87,FMT=1000) eig(n)
                DO j=1,n2max
                   WRITE (87,FMT=1010) ac(j,n),bc(j,n)
                END DO
             END IF
          END DO
990       FORMAT(2(f8.4,1x),i3,1x,i3)
1000      FORMAT(e10.4)
1010      FORMAT(2(2e20.8,1x))
       END DO
       DEALLOCATE ( ac,bc )
       !
       CLOSE(89)

    ENDIF
    RETURN
  END SUBROUTINE doswrite
END MODULE m_doswrite
