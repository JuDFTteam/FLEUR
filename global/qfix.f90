MODULE m_qfix
  USE m_juDFT
  !Calculate total charge
  !Depending on variable input%qfix, the following will be done to fix the charge
  !Input qfix can be 1 or 2
  ! qfix=0 (no qfix in inp.xml) means we usually do not run the code
  !                       if force_fix is .true. we run the code and assume qfix=2
  !                       in the call to qfix we will always run it
  ! qfix=1 (qfix=f in inp.xml) means we fix only in INT (only done in firstcall)
  ! qfix=2 (qfix=t in inp.xml) means we fix total charge
  ! qfix file no longer supported!
  ! If l_par=.TRUE., MPI parallelization in the cdntot will be used.
  ! Be carefull not to set it to .TRUE. if you are calling only from one MPI
  ! rank.

CONTAINS
  SUBROUTINE qfix(fmpi,stars,atoms,sym,vacuum,sphhar,input,cell ,&
                  den,l_noco,l_printData,l_par,force_fix,fix,fix_pw_only)

    USE m_types
    USE m_constants
    USE m_cdntot
    USE m_xmlOutput

    IMPLICIT NONE

    !     .. Scalar Arguments ..
    TYPE(t_mpi),INTENT(IN)       :: fmpi
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_input),INTENT(IN)     :: input

    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_potden),INTENT(INOUT) :: den
    LOGICAL,INTENT(IN)           :: l_noco,l_printData,force_fix,l_par
    REAL,    INTENT (OUT)        :: fix
    LOGICAL,INTENT(IN),OPTIONAL  :: fix_pw_only
    !     .. Local Scalars ..
    LOGICAL :: l_qfixfile,fixtotal
    LOGICAL :: l_firstcall=.true.
    REAL    :: qtot,qis,zc
    INTEGER :: jm,lh,n,na
    !     ..
    fixtotal=.true. !this is the default
    IF (PRESENT(fix_pw_only)) fixtotal=.NOT.fix_pw_only
    fix=1.0
    IF (l_firstcall) THEN
       l_firstcall=.false.
    ELSE
       IF (MOD(input%qfix,2)==0.AND..NOT.force_fix) RETURN
    ENDIF
    ! qfix==0 means no qfix was given in inp.xml.
    ! In this case do nothing except when forced to fix!

    CALL cdntot(stars,atoms,sym,vacuum,input,cell ,den,l_printData,qtot,qis,fmpi,l_par)

    IF (fmpi%irank.EQ.0) THEN
       !The total nucleii charge
       zc=SUM(atoms%neq(:)*atoms%zatom(:))
       !zc = zc + 2*input%sigma !TODO : reactivate fields

       IF (fixtotal) THEN
          !-roa
          fix = zc/qtot
          na = 1
          DO n = 1,atoms%ntype
             lh = sphhar%nlh(sym%ntypsy(na))
             jm = atoms%jri(n)
             den%mt(:jm,0:lh,n,:) = fix*den%mt(:jm,0:lh,n,:)
             na = na + atoms%neq(n)
          ENDDO
          den%pw(:stars%ng3,:) = fix*den%pw(:stars%ng3,:)
          IF (input%film) THEN
             den%vacz(:vacuum%nmz,:vacuum%nvac,:) = fix*den%vacz(:vacuum%nmz,:vacuum%nvac,:)
             den%vacxy(:vacuum%nmzxy,:stars%ng2-1,:vacuum%nvac,:) = fix*&
                den%vacxy(:vacuum%nmzxy,:stars%ng2-1,:vacuum%nvac,:)
          END IF
          WRITE (oUnit,FMT=8000) zc,fix
       ELSE
          fix = (zc - qtot) / qis + 1.
          den%pw(:stars%ng3,:) = fix*den%pw(:stars%ng3,:)
          WRITE (oUnit,FMT=8001) zc,fix
       ENDIF

       ! TODO: This looks spooky.
       ! a) All noco quantities are already included in the fix above.
       ! b) Below, MT is missing.
       !IF (l_noco) THEN
          !fix also the off-diagonal part of the density matrix
          !den%pw(:stars%ng3,3) = fix*den%pw(:stars%ng3,3)
          !IF (input%film.AND.fixtotal) THEN
             !den%vacz(:,:,3:4) = fix*den%vacz(:,:,3:4)
             !den%vacxy(:,:,:,3) = fix*den%vacxy(:,:,:,3)
          !END IF
       !END IF

       IF (ABS(fix-1.0)<1.E-6) RETURN !no second calculation of cdntot as nothing was fixed

       IF(l_printData) CALL openXMLElementNoAttributes('fixedCharges')
       CALL cdntot(stars,atoms,sym,vacuum,input,cell ,den,l_printData,qtot,qis,fmpi,.FALSE.)
       IF(l_printData) CALL closeXMLElement('fixedCharges')

       IF (fix>1.1) CALL juDFT_WARN("You lost too much charge")
       IF (fix<.9) CALL juDFT_WARN("You gained too much charge")

8000   FORMAT (/,10x,'zc= ',f12.6,5x,'qfix=  ',f10.6)
8001   FORMAT (/,' > broy only qis: ','zc= ',f12.6,5x,'qfix=  ',f10.6)
       !-roa
    ENDIF

  END SUBROUTINE qfix
END MODULE m_qfix
