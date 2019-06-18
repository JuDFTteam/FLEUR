!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rinpXML
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!
   !!! The routine r_inpXML reads in the inp.xml file
   !!!
   !!!                               GM'16
   !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
   SUBROUTINE r_inpXML(&
      atoms,vacuum,input,stars,sliceplot,banddos,forcetheo,field,&
      cell,sym,xcpot,noco,oneD,hybrid,kpts,enpara,coreSpecInput,wann,&
      noel,namex,relcor,a1,a2,a3,dtild)

      USE iso_c_binding
      USE m_juDFT
      USE m_types
      USE m_types_forcetheo_extended
      USE m_symdata , ONLY : nammap, ord2, l_c2
      !USE m_rwsymfile
      !USE m_xmlIntWrapFort
      USE m_inv3
      !USE m_spg2set
      !USE m_closure, ONLY : check_close
      !USE m_symproperties
      USE m_calculator
      USE m_constants
      !USE m_inpeig
      USE m_sort
      USE m_types_xcpot_inbuild
#ifdef CPP_LIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE

      TYPE(t_input),INTENT(INOUT)   :: input!ok
      TYPE(t_sym),INTENT(INOUT)     :: sym !ok
      TYPE(t_stars),INTENT(INOUT)   :: stars !ok
      TYPE(t_atoms),INTENT(INOUT)   :: atoms
      TYPE(t_vacuum),INTENT(INOUT)   :: vacuum !ok
      TYPE(t_kpts),INTENT(INOUT)     :: kpts!ok
      TYPE(t_oneD),INTENT(INOUT)     :: oneD !ok
      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid !ok
      TYPE(t_cell),INTENT(INOUT)     :: cell!ok
      TYPE(t_banddos),INTENT(INOUT)  :: banddos!ok
      TYPE(t_sliceplot),INTENT(INOUT):: sliceplot !ok
      CLASS(t_xcpot),INTENT(INOUT),ALLOCATABLE :: xcpot!ok
      TYPE(t_noco),INTENT(INOUT)     :: noco!ok
      TYPE(t_enpara)   ,INTENT(OUT)  :: enpara !ok
      TYPE(t_field), INTENT(INOUT)   :: field !ok
      CLASS(t_forcetheo),ALLOCATABLE,INTENT(OUT):: forcetheo!ok
      TYPE(t_coreSpecInput),INTENT(OUT) :: coreSpecInput !ok
      TYPE(t_wann)   ,INTENT(INOUT)  :: wann !ok

      !+odim
      INTEGER MM,vM,m_cyl
      LOGICAL invs1,zrfs1
      INTEGER chi,rot
      LOGICAL d1,band
      NAMELIST /odim/ d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
      !-odim
      ! ..
      ! ..  Local Variables
      REAL     :: scpos  ,zc
      INTEGER ieq,i,k,na,n,ii,vxc_id_c,vxc_id_x,exc_id_c,exc_id_x
      REAL s3,ah,a,hs2,rest,thetaj
      LOGICAL l_hyb,l_sym,ldum
      INTEGER :: ierr
      ! ..
      !...  Local Arrays
      !   CHARACTER :: helpchar(atoms%ntype)
      CHARACTER(len=  4) :: chntype
      CHARACTER(len= 41) :: chform
      CHARACTER(len=100) :: line

      CHARACTER(len=20) :: tempNumberString
      CHARACTER(len=150) :: format
      CHARACTER(len=20) :: mixingScheme
      CHARACTER(len=10) :: loType
      LOGICAL :: kptGamma, l_relcor,ldummy
      INTEGER :: iAtomType
      CHARACTER(len=100) :: xPosString, yPosString, zPosString
      CHARACTER(len=20),ALLOCATABLE :: speciesName(:)
      
      !   REAL :: tempTaual(3,atoms%nat)
      TYPE(t_econfig)  :: econf

      INTEGER            :: iType, iLO, iSpecies, lNumCount, nNumCount, iLLO, jsp, j, l, absSum, numTokens
      INTEGER            :: numberNodes, nodeSum, numSpecies, n2spg, n1, n2, ikpt, iqpt
      INTEGER            :: atomicNumber,  gridPoints, lmax, lnonsphr, lmaxAPW
      INTEGER            :: latticeDef, symmetryDef, nop48, firstAtomOfType, errorStatus
      INTEGER            :: loEDeriv, ntp1, ios, ntst, jrc, minNeigd
      INTEGER            :: nv, nv2, kq1, kq2, kq3, nprncTemp, kappaTemp, tempInt
      INTEGER            :: ldau_l(4), numVac, numU
      INTEGER            :: speciesEParams(0:3)
      INTEGER            :: mrotTemp(3,3,48)
      REAL               :: tauTemp(3,48)
      REAL               :: bk(3)
      LOGICAL            :: flipSpin, l_eV, invSym, l_qfix, relaxX, relaxY, relaxZ
      LOGICAL            :: coreConfigPresent, l_enpara, l_orbcomp, tempBool, l_nocoinp
      REAL               :: magMom, radius, logIncrement, qsc(3), latticeScale, dr
      REAL               :: aTemp, zp, rmtmax, sumWeight, ldau_u(4), ldau_j(4), tempReal
      REAL               :: ldau_phi(4),ldau_theta(4)
      REAL               :: weightScale, eParamUp, eParamDown
      LOGICAL            :: l_amf(4)
      INTEGER            :: lcutm,lcutwf,hybSelect(4)
      REAL               :: evac0Temp(2,2)

      CHARACTER(LEN=200,KIND=c_char) :: schemaFilename, docFilename
      CHARACTER(LEN=255) :: valueString, lString, nString, token
      CHARACTER(LEN=255) :: xPathA, xPathB, xPathC, xPathD, xPathE
      CHARACTER(LEN=11)  :: latticeType
      CHARACTER(LEN=50)  :: versionString
      CHARACTER(LEN=150) :: kPointsPrefix

      INTEGER            :: altKPointSetIndex,  altKPointSetIndices(2)
      LOGICAL            :: ldaSpecies
      REAL               :: socscaleSpecies,b_field_mtspecies,vcaspecies

      INTEGER, ALLOCATABLE :: lNumbers(:), nNumbers(:), speciesLLO(:)
      INTEGER, ALLOCATABLE :: loOrderList(:)
      INTEGER, ALLOCATABLE :: speciesNLO(:)
      INTEGER, ALLOCATABLE :: multtab(:,:), invOps(:), optype(:)
      INTEGER, ALLOCATABLE :: lmx1(:), nq1(:), nlhtp1(:)
      INTEGER, ALLOCATABLE :: speciesLOEDeriv(:)
      REAL,    ALLOCATABLE :: speciesLOeParams(:), speciesLLOReal(:)
      LOGICAL, ALLOCATABLE :: wannAtomList(:)

      ! Variables for MT radius testing:

      REAL                 :: dtild1,kmax1,dvac1
      LOGICAL              :: l_test
      INTEGER, ALLOCATABLE :: jri1(:), lmax1(:)
      REAL, ALLOCATABLE    :: rmt1(:), dx1(:)


    

     
      IF (ABS(input%fixed_moment)>1E-8.AND.(input%jspins==1.OR.noco%l_noco)) CALL judft_error("Fixed moment only in collinear calculations with two spins")

 
    
   
      

     
    
      ! Construction of amat requires additional information about the lattice
      ! and is done later (scroll down)!

         IF(input%film) THEN
         evac0Temp = eVac0Default_const
            xPathB = TRIM(ADJUSTL(xPathA))//'/vacuumEnergyParameters'
            numberNodes = xmlGetNumberOfNodes(xPathB)
            IF(numberNodes.GE.1) THEN
               DO i = 1, numberNodes
                  xPathC = ''
                  WRITE(xPathC,'(a,i0,a)') TRIM(ADJUSTL(xPathB))//'[',i,']'
                  numVac = evaluateFirstIntOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@vacuum'))
                  eParamUp = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@spinUp'))
                  eParamDown = eParamUp
                  IF (xmlGetNumberOfNodes(TRIM(ADJUSTL(xPathC))//'/@spinDown').GE.1) THEN
                     eParamDown = evaluateFirstOnly(xmlGetAttributeValue(TRIM(ADJUSTL(xPathC))//'/@spinDown'))
                  END IF
                  evac0Temp(numVac,1) = eParamUp
                  IF(input%jspins.GT.1) evac0Temp(numVac,2) = eParamDown
                  IF(i == 1) THEN
                     evac0Temp(3-numVac,1) = eParamUp
                     IF(input%jspins.GT.1) evac0Temp(3-numVac,2) = eParamDown
                  END IF
               END DO
            END IF
      

      IF (.NOT.input%film) vacuum%dvac = a3(3)
      vacuum%dvac = latticeScale*vacuum%dvac
      dtild = latticeScale*dtild


  !some settings for film calculations
     
      IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1
      IF (oneD%odd%d1) vacuum%nvac = 1
    
      IF (oneD%odd%d1) vacuum%delz = 20.0 / vacuum%nmz
    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Start of force-theorem section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      xPathA = '/fleurInput/forceTheorem'
      numberNodes = xmlGetNumberOfNodes(xPathA)
      IF (numberNodes.EQ.1) THEN
         !Magnetic anisotropy...
         xPathA = '/fleurInput/forceTheorem/MAE'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF
         !spin-spiral dispersion
         xPathA = '/fleurInput/forceTheorem/spinSpiralDispersion'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF
         !dmi
         
         !jij
         xPathA = '/fleurInput/forceTheorem/Jij'
         numberNodes = xmlGetNumberOfNodes(xPathA)
         IF (numberNodes.EQ.1) THEN
         ENDIF

      ELSE
         ALLOCATE(t_forcetheo::forcetheo) !default no forcetheorem type
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! End of force-theorem section
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
         

         

      !CALL xmlFreeResources()
      
      !WRITE(*,*) 'Reading of inp.xml file finished'

      DEALLOCATE(speciesNLO)

   END SUBROUTINE r_inpXML

   SUBROUTINE setXCParameters(atoms,namex,relcor,jspins,vxc_id_x,vxc_id_c,exc_id_x,exc_id_c,xcpot)
      USE m_juDFT
      USE m_types
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_libxc

      IMPLICIT NONE
      TYPE(t_atoms),INTENT(IN)          :: atoms
      CHARACTER(LEN=*),     INTENT(IN)  :: namex
      LOGICAL,              INTENT(IN)  :: relcor
      INTEGER,              INTENT(IN)  :: jspins,vxc_id_c,vxc_id_x,exc_id_x,exc_id_c
      CLASS(t_xcpot),INTENT(OUT),ALLOCATABLE      :: xcpot

      IF (namex(1:5)=='LibXC') THEN
         ALLOCATE(t_xcpot_libxc::xcpot)
      ELSE
         ALLOCATE(t_xcpot_inbuild::xcpot)
      ENDIF

      SELECT TYPE(xcpot)
      TYPE IS(t_xcpot_inbuild)
         CALL xcpot%init(namex(1:4),relcor,atoms%ntype)
      TYPE IS(t_xcpot_libxc)
         CALL xcpot%init(jspins,vxc_id_x,vxc_id_c,exc_id_x,exc_id_c)

      END SELECT

      CALL set_xcpot_usage(xcpot)
   END SUBROUTINE setXCParameters

   SUBROUTINE set_xcpot_usage(xcpot)
      use m_judft_usage
      USE m_types
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_libxc
      implicit none
      class(t_xcpot), intent(in)    :: xcpot

      ! give some information about XC functional to usage.json
      ! 1 -> LDA
      ! 2 -> GGA
      ! 3 -> MetaGGA
      ! 4 -> Hybrid functional
      if(xcpot%vxc_is_lda()) then
         call add_usage_data("XC-treatment", 1)
         return
      endif

      if(xcpot%exc_is_MetaGGA()) then
         call add_usage_data("XC-treatment", 3)
         return
      endif

      if(xcpot%vxc_is_GGA()) then
         call add_usage_data("XC-treatment", 2)
         return
      endif

      if(xcpot%is_hybrid()) then
         call add_usage_data("XC-treatment", 4)
         return
      endif

   END SUBROUTINE set_xcpot_usage

   SUBROUTINE getIntegerSequenceFromString(string, sequence, count)
      use m_juDFT_stop
      IMPLICIT NONE

      CHARACTER(*),         INTENT(IN)  :: string
      INTEGER, ALLOCATABLE, INTENT(OUT) :: sequence(:)
      INTEGER,              INTENT(OUT) :: count

      INTEGER :: i, length, start, lastNumber, currentNumber, index
      LOGICAL singleNumber, comma, dash

      ! 3 cases: 1. a single number
      !          2. number - number
      !          3. comma separated numbers

      length = LEN(string)
      count = 0
      start = 1
      singleNumber = .TRUE.
      comma = .FALSE.
      dash = .FALSE.
      lastNumber = 0
      count = 0

      ! 1. Determine number count

      DO i = 1, length
         SELECT CASE (string(i:i))
         CASE ('0':'9')
         CASE (',')
            IF ((start.EQ.i).OR.(dash)) THEN
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            comma = .TRUE.
            READ(string(start:i-1),*) lastNumber
            count = count + 1
            start = i+1
         CASE ('-')
            IF ((start.EQ.i).OR.(dash).OR.(comma)) THEN
               CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
            END IF
            singleNumber = .FALSE.
            dash = .TRUE.
            READ(string(start:i-1),*) lastNumber
            start = i+1
         CASE DEFAULT
            CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
         END SELECT
      END DO
      IF(start.GT.length) THEN
         CALL judft_error('String has wrong syntax (in getIntegerSequenceFromString)')
      END IF
      READ(string(start:length),*) currentNumber
      IF (dash) THEN
         count = currentNumber - lastNumber + 1
      ELSE
         count = count + 1
      END IF

      IF (ALLOCATED(sequence)) THEN
         DEALLOCATE(sequence)
      END IF
      ALLOCATE(sequence(count))

      ! 2. Read in numbers iff comma separation ...and store numbers in any case

      IF (singleNumber) THEN
         sequence(1) = currentNumber
      ELSE IF (dash) THEN
         DO i = 1, count
            sequence(i) = lastNumber + i - 1
         END DO
      ELSE
         index = 1
         start = 1
         DO i = 1, length
            SELECT CASE (string(i:i))
            CASE (',')
               comma = .TRUE.
               READ(string(start:i-1),*) lastNumber
               start = i+1
               SEQUENCE(index) = lastNumber
               index = index + 1
            END SELECT
         END DO
         sequence(index) = currentNumber
      END IF

   END SUBROUTINE getIntegerSequenceFromString

  
END MODULE m_rinpXML
