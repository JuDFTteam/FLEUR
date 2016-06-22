MODULE m_winpXML

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!   XML input file generator
!!!
!!!   This subroutine is supposed to write out a file inp.xml
!!!   containing all required input data.
!!!                                         GM'16
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
SUBROUTINE w_inpXML(&
&                   atoms,obsolete,vacuum,input,stars,sliceplot,banddos,&
&                   cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,div,l_gamma,&
&                   noel,namex,relcor,a1,a2,a3,scale,dtild_opt,name_opt,&
&                   xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
&                   atomTypeSpecies,speciesRepAtomType,l_outFile,numSpecies,&
&                   enpara)

   USE m_types
   USE m_juDFT_init
   USE m_constants
   USE m_xmlOutput

   IMPLICIT NONE

! arguments

   TYPE(t_input),INTENT(IN)   :: input
   TYPE(t_sym),INTENT(IN)     :: sym
   TYPE(t_stars),INTENT(IN)   :: stars 
   TYPE(t_atoms),INTENT(IN)   :: atoms
   TYPE(t_vacuum),INTENT(IN)   :: vacuum
   TYPE(t_obsolete),INTENT(IN) :: obsolete
   TYPE(t_kpts),INTENT(IN)     :: kpts
   TYPE(t_oneD),INTENT(IN)     :: oneD
   TYPE(t_hybrid),INTENT(IN)   :: hybrid
   TYPE(t_Jij),INTENT(IN)      :: Jij
   TYPE(t_cell),INTENT(IN)     :: cell
   TYPE(t_banddos),INTENT(IN)  :: banddos
   TYPE(t_sliceplot),INTENT(IN):: sliceplot
   TYPE(t_xcpot),INTENT(IN)    :: xcpot
   TYPE(t_noco),INTENT(IN)     :: noco
   TYPE(t_enpara),INTENT(IN)   :: enpara
   INTEGER, INTENT (IN)        :: numSpecies
   INTEGER, INTENT (IN)        :: div(3)
   INTEGER, INTENT (IN)        :: atomTypeSpecies(atoms%ntype)
   INTEGER, INTENT (IN)        :: speciesRepAtomType(numSpecies)
   LOGICAL, INTENT (IN)        :: l_gamma, l_outFile
   REAL,    INTENT (IN)        :: a1(3),a2(3),a3(3),scale
   REAL, INTENT (IN)     :: xmlCoreOccs(2,29,atoms%ntype)
   INTEGER, INTENT (IN)  :: xmlElectronStates(29,atoms%ntype)
   LOGICAL, INTENT (IN)  :: xmlPrintCoreStates(29,atoms%ntype)
   CHARACTER(len=3),INTENT(IN) :: noel(atoms%ntypd)
   CHARACTER(len=4),INTENT(IN) :: namex
   CHARACTER(len=12),INTENT(IN):: relcor
   REAL,INTENT(IN),OPTIONAL    :: dtild_opt
   CHARACTER(len=8),INTENT(IN),OPTIONAL:: name_opt(10)


   INTEGER :: iSpecies, fileNum
   CHARACTER(len=8) :: name(10)

!+lda+u
   REAL    u,j
   INTEGER l
   LOGICAL l_amf
   CHARACTER(len=3) ch_test
   NAMELIST /ldaU/ l,u,j,l_amf
!-lda+u
!+odim
   INTEGER MM,vM,m_cyl
   LOGICAL invs1,zrfs1
   INTEGER chi,rot
   LOGICAL d1,band
   NAMELIST /odim/ d1,MM,vM,m_cyl,chi,rot,invs1,zrfs1
!-odim
! ..
! ..  Local Variables
   REAL     ::dtild ,scpos, zc, sumWeight
   INTEGER  ::nw,idsprs, n1, n2
   INTEGER ieq,i,k,na,n,ilo
   REAL s3,ah,a,hs2,rest
   LOGICAL l_hyb,l_sym,ldum
   INTEGER :: ierr
! ..
!...  Local Arrays
   CHARACTER :: helpchar(atoms%ntypd)
   CHARACTER(len=  4) :: chntype
   CHARACTER(len= 41) :: chform
   CHARACTER(len=100) :: line

!     added for HF and hybrid functionals
   REAL                  ::  aMix,omega
   INTEGER               :: idum
   CHARACTER (len=1)     ::  check

   CHARACTER(len=20) :: tempNumberString, speciesName
   CHARACTER(len=150) :: format
   CHARACTER(len=20) :: mixingScheme
   CHARACTER(len=10) :: loType
   CHARACTER(len=10) :: bzIntMode
   CHARACTER(len=200) :: symFilename
   LOGICAL :: kptGamma, l_relcor, l_explicit
   INTEGER :: iAtomType, startCoreStates, endCoreStates
   CHARACTER(len=100) :: xPosString, yPosString, zPosString
   CHARACTER(len=200) :: coreStatesString, valenceStatesString
   REAL :: tempTaual(3,atoms%nat)
   REAL :: a1Temp(3),a2Temp(3),a3Temp(3)
   REAL :: amatTemp(3,3), bmatTemp(3,3)
   CHARACTER(len=7) :: coreStateList(29) !'(1s1/2)'
   CHARACTER(len=4) :: nobleGasConfigList(6) !'[He]'

   DATA coreStateList / '(1s1/2)','(2s1/2)','(2p1/2)','(2p3/2)','(3s1/2)',&
&                       '(3p1/2)','(3p3/2)','(3d3/2)','(3d5/2)','(4s1/2)',&
&                       '(4p1/2)','(4p3/2)','(5s1/2)','(4d3/2)','(4d5/2)',&
&                       '(5p1/2)','(5p3/2)','(6s1/2)','(4f5/2)','(4f7/2)',&
&                       '(5d3/2)','(5d5/2)','(6p1/2)','(6p3/2)','(7s1/2)',&
&                       '(5f5/2)','(5f7/2)','(6d3/2)','(6d5/2)' /

   DATA nobleGasConfigList / '[He]','[Ne]','[Ar]','[Kr]','[Xe]','[Rn]' /

   IF (PRESENT(dtild_opt)) dtild=dtild_opt
   IF (PRESENT(name_opt)) name=name_opt

   symFilename = 'sym.out'
   kptGamma = l_gamma
   band = .false.
   nw=1
   IF (TRIM(ADJUSTL(namex)).EQ.'hf'.OR.TRIM(ADJUSTL(namex)).EQ.'exx'.OR.&
       TRIM(ADJUSTL(namex)).EQ.'hse'.OR.TRIM(ADJUSTL(namex)).EQ.'vhse') l_hyb = .true.
   l_relcor=.true.
   IF(relcor.EQ.'relativi') THEN
      l_relcor=.true.
   ELSE 
      l_relcor=.false.
   END IF

   DO i = 1, 3
      a1Temp(i) = a1(i)
      a2Temp(i) = a2(i)
      a3Temp(i) = a3(i)
   END DO

   fileNum = -1
   IF(l_outFile) THEN
      fileNum = getXMLOutputUnitNumber()
   ELSE
      fileNum = 5
      OPEN (fileNum,file='inp.xml',form='formatted',status='unknown')
      REWIND (fileNum)

      WRITE (fileNum,'(a)') '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
      WRITE (fileNum,'(a)') '<fleurInput fleurInputVersion="0.27">'
   END IF

   IF(PRESENT(name_opt)) THEN
      WRITE (fileNum,'(a)') '   <comment>'
      WRITE (fileNum,'(a6,10a8)') '      ',name
      WRITE (fileNum,'(a)') '   </comment>'
   END IF

   WRITE (fileNum,'(a)') '   <calculationSetup>'

!      <cutoffs Kmax="3.60000" Gmax="11.000000" GmaxXC="9.200000" numbands="0"/>
   110 FORMAT('      <cutoffs Kmax="',f0.8,'" Gmax="',f0.8,'" GmaxXC="',f0.8,'" numbands="',i0,'"/>')
   WRITE (fileNum,110) input%rkmax,stars%gmax,xcpot%gmaxxc,input%gw_neigd

!      <scfLoop itmax="9" maxIterBroyd="99" imix="Anderson" alpha="0.05" spinf="2.00"/>
   120 FORMAT('      <scfLoop itmax="',i0,'" maxIterBroyd="',i0,'" imix="',a,'" alpha="',f0.8,'" spinf="',f0.8,'"/>')
   SELECT CASE (input%imix)
      CASE (1) 
         mixingScheme='straight'
      CASE (3) 
         mixingScheme='Broyden1'
      CASE (5) 
         mixingScheme='Broyden2'
      CASE (7) 
         mixingScheme='Anderson'
      CASE DEFAULT 
         mixingScheme='errorUnknownMixing'
   END SELECT
   WRITE (fileNum,120) input%itmax,input%maxiter,TRIM(mixingScheme),input%alpha,input%spinf

!      <coreElectrons ctail="T" frcor="F" kcrel="0"/>
   130 FORMAT('      <coreElectrons ctail="',l1,'" frcor="',l1,'" kcrel="',i0,'"/>')
   WRITE (fileNum,130) input%ctail,input%frcor,input%kcrel

!      <magnetism jspins="1" l_noco="F" l_J="F" swsp="F" lflip="F"/>
   140 FORMAT('      <magnetism jspins="',i0,'" l_noco="',l1,'" l_J="',l1,'" swsp="',l1,'" lflip="',l1,'"/>')
   WRITE (fileNum,140) input%jspins,noco%l_noco,jij%l_J,input%swsp,input%lflip

!      <soc theta="0.00000" phi="0.00000" l_soc="F" spav="F" off="F" soc66="F"/>
   150 FORMAT('      <soc theta="',f0.8,'" phi="',f0.8,'" l_soc="',l1,'" spav="',l1,'" off="',l1,'" soc66="',l1,'"/>')
   WRITE (fileNum,150) noco%theta,noco%phi,noco%l_soc,noco%soc_opt(atoms%ntype+2),noco%soc_opt(atoms%ntype+1),obsolete%eig66(2)

   IF (noco%l_noco) THEN
      160 FORMAT('      <nocoParams l_ss="',l1,'" l_mperp="',l1,'" l_constr="',l1,'" l_disp="',l1,'" sso_opt="',a3,'" mix_b="',f0.8,'" thetaJ="',f0.8,'" nsh="',i0,'"/>')
      STOP 'Output of Noco input not yet implemented!'
   END IF

   IF (oneD%odd%d1) THEN
      170 FORMAT('      <oneDParams d1="',l1,'" MM="',i0,'" vM="',i0,'" m_cyl="',i0,'" chi="',i0,'" rot="',i0,'" invs1="',l1,'" zrfs1="',l1,'"/>')
      WRITE (fileNum,170) oneD%odd%d1,oneD%odd%M,oneD%odd%mb,oneD%odd%m_cyl,oneD%odd%chi,oneD%odd%rot,oneD%odd%invs,oneD%odd%zrfs
   END IF

!      <expertModes gw="0" pot8="F" eig66="F" lpr="0" isec1="99" secvar="F" />
   180 FORMAT('      <expertModes gw="',i0,'" pot8="',l1,'" eig66="',l1,'" lpr="',i0,'" isec1="',i0,'" secvar="',l1,'"/>')
   WRITE (fileNum,180) input%gw,obsolete%pot8,obsolete%eig66(1),obsolete%lpr,input%isec1,input%secvar

!      <geometryOptimization l_f="F" xa="2.00000" thetad="330.00000" epsdisp="0.00001" epsforce="0.00001"/>
   190 FORMAT('      <geometryOptimization l_f="',l1,'" xa="',f0.8,'" thetad="',f0.8,'" epsdisp="',f0.8,'" epsforce="',f0.8,'"/>')
   WRITE (fileNum,190) input%l_f,input%xa,input%thetad,input%epsdisp,input%epsforce

   IF(input%gauss.AND.input%tria) THEN
      STOP 'Error: bz integration modes gauss AND tria selected!'
   END IF

   bzIntMode = 'hist'
   IF(input%gauss) THEN
      bzIntMode = 'gauss'
   ELSE IF(input%tria) THEN
      bzIntMode = 'tria'
   END IF
!      <bzIntegration valenceElectrons="8.00000" mode="hist" fermiSmearingEnergy="0.00100">
   200 FORMAT('      <bzIntegration valenceElectrons="',f0.8,'" mode="',a,'" fermiSmearingEnergy="',f0.8,'">')
   WRITE (fileNum,200) input%zelec,TRIM(ADJUSTL(bzIntMode)),input%tkb

   l_explicit = juDFT_was_argument("-explicit").OR.l_outFile
   IF(l_explicit) THEN
      sumWeight = 0.0
      DO i = 1, kpts%nkpt
         sumWeight = sumWeight + kpts%weight(i)
      END DO
      205 FORMAT('         <kPointList posScale="',f0.8,'" weightScale="',f0.8,'" count="',i0,'">')
      WRITE (fileNum,205) kpts%posScale, sumWeight, kpts%nkpt
      DO i = 1, kpts%nkpt
         206 FORMAT('            <kPoint weight="',f12.6,'">',f12.6,' ',f12.6,' ',f12.6,'</kPoint>')
         WRITE (fileNum,206) kpts%weight(i), kpts%bk(1,i), kpts%bk(2,i), kpts%bk(3,i)
      END DO
      WRITE (fileNum,'(a)')('         </kPointList>')
   ELSE IF( (div(1) == 0).OR.(div(2) == 0) ) THEN
!            <kPointCount count="100" gamma="F"/>
      208 FORMAT('         <kPointCount count="',i0,'" gamma="',l1,'"/>')
      WRITE (fileNum,208) kpts%nkpt,kptGamma
   ELSE
!            <kPointMesh nx="10" ny="10" nz="10" gamma="F"/>
      210 FORMAT('         <kPointMesh nx="',i0,'" ny="',i0,'" nz="',i0,'" gamma="',l1,'"/>')
      WRITE (fileNum,210) div(1),div(2),div(3),kptGamma
   END IF
   WRITE (fileNum,'(a)') '      </bzIntegration>'

!      <energyParameterLimits ellow="-2.00000" elup="2.00000"/>
   220 FORMAT('      <energyParameterLimits ellow="',f0.8,'" elup="',f0.8,'"/>')
   WRITE (fileNum,220) input%ellow,input%elup

   WRITE (fileNum,'(a)') '   </calculationSetup>'
   WRITE (fileNum,'(a)') '   <cell>'

   IF(l_explicit) THEN
      WRITE(fileNum,'(a)') '      <symmetryOperations>'
      DO i = 1, sym%nop
      WRITE(fileNum,'(a)') '         <symOp>'
      224 FORMAT('            <row-1>',i0,' ',i0,' ',i0,' ',f0.15,'</row-1>')
      WRITE(fileNum,224) sym%mrot(1,1,i), sym%mrot(1,2,i), sym%mrot(1,3,i), sym%tau(1,i)
      225 FORMAT('            <row-2>',i0,' ',i0,' ',i0,' ',f0.15,'</row-2>')
      WRITE(fileNum,225) sym%mrot(2,1,i), sym%mrot(2,2,i), sym%mrot(2,3,i), sym%tau(2,i)
      226 FORMAT('            <row-3>',i0,' ',i0,' ',i0,' ',f0.15,'</row-3>')
      WRITE(fileNum,226) sym%mrot(3,1,i), sym%mrot(3,2,i), sym%mrot(3,3,i), sym%tau(3,i)
      WRITE(fileNum,'(a)') '         </symOp>'
      END DO
      WRITE(fileNum,'(a)') '      </symmetryOperations>'
   ELSE IF(TRIM(ADJUSTL(sym%namgrp)).EQ.'any') THEN
      228 FORMAT('      <symmetryFile filename="',a,'"/>')
      WRITE(fileNum,228) TRIM(ADJUSTL(symFilename))
   ELSE
!      <symmetry spgrp="any" invs="T" zrfs="F"/>
      230 FORMAT('      <symmetry spgrp="',a,'" invs="',l1,'" zrfs="',l1,'"/>')
      WRITE (fileNum,230) TRIM(ADJUSTL(sym%namgrp)),sym%invs,sym%zrfs

   END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note: Different options for the cell definition!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF (cell%latnam.EQ.'c-b') THEN
      a1Temp(1) = sqrt(2.)* a1Temp(1)
   END IF
   IF (cell%latnam.EQ.'hex') THEN
      s3 = sqrt(3.)
      a1Temp(1) = 2*a1Temp(1)/sqrt(3.)
   END IF
   IF (cell%latnam.EQ.'hx3') THEN
      a1Temp(1) = 2*a1Temp(1)
   END IF

   IF (input%film) THEN
!      <xsd:attribute name="dVac" type="xsd:double" use="required"/>
!      <xsd:attribute name="dTilda" type="xsd:double" use="required"/>
!      <filmLattice ...>
      241 FORMAT('      <filmLattice scale="',f0.8,'" latnam="',a,'" dVac="',f0.8,'" dTilda="',f0.8,'">')
      WRITE(fileNum,241) scale, TRIM(ADJUSTL(cell%latnam)), vacuum%dvac, dtild
      IF (cell%latnam.EQ.'any') THEN
         WRITE (fileNum,'(a)') '         <bravaisMatrix>'
         255 FORMAT('            <row-1>',f0.12,' ',f0.12,' ',f0.12,'</row-1>')
         WRITE (fileNum,255) a1Temp(1),a1Temp(2),a1Temp(3)
         265 FORMAT('            <row-2>',f0.12,' ',f0.12,' ',f0.12,'</row-2>')
         WRITE (fileNum,265) a2Temp(1),a2Temp(2),a2Temp(3)
         275 FORMAT('            <row-3>',f0.12,' ',f0.12,' ',f0.12,'</row-3>')
         WRITE (fileNum,275) a3Temp(1),a3Temp(2),a3Temp(3)
         WRITE (fileNum,'(a)') '         </bravaisMatrix>'
      ELSE
         IF ((cell%latnam.EQ.'squ').OR.(cell%latnam.EQ.'hex').OR.&
     &       (cell%latnam.EQ.'c-b').OR.(cell%latnam.EQ.'hx3').OR.&
     &       (cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
            256 FORMAT('         <a1>',f0.12,'</a1>')
            WRITE (fileNum,256) a1Temp(1)
         END IF
         IF ((cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
            266 FORMAT('         <a2>',f0.12,'</a2>')
            WRITE (fileNum,266) a2Temp(2)
         END IF

         IF (cell%latnam.EQ.'obl') THEN
            257 FORMAT('         <row-1>',f0.12,' ',f0.12,'</row-1>')
            WRITE (fileNum,257) a1Temp(1), a1Temp(2)
            267 FORMAT('         <row-2>',f0.12,' ',f0.12,'</row-2>')
            WRITE (fileNum,267) a2Temp(1), a2Temp(2)
         END IF
      END IF

      WRITE (fileNum,'(a)') '      </filmLattice>'
   ELSE

      242 FORMAT('      <bulkLattice scale="',f0.12,'" latnam="',a,'">')
      WRITE (fileNum,242) scale, TRIM(ADJUSTL(cell%latnam))

      IF (cell%latnam.EQ.'any') THEN

!         <bravaisMatrix scale="1.0000000">
         WRITE (fileNum,'(a)') '         <bravaisMatrix>'

!            <row-1>0.00000 5.13000 5.13000</row-1>
         250 FORMAT('            <row-1>',f0.12,' ',f0.12,' ',f0.12,'</row-1>')
         WRITE (fileNum,250) a1Temp(1),a1Temp(2),a1Temp(3)
!            <row-2>5.13000 0.00000 5.13000</row-2>
         260 FORMAT('            <row-2>',f0.12,' ',f0.12,' ',f0.12,'</row-2>')
         WRITE (fileNum,260) a2Temp(1),a2Temp(2),a2Temp(3)
!            <row-3>5.13000 5.13000 0.00000</row-3>
         270 FORMAT('            <row-3>',f0.12,' ',f0.12,' ',f0.12,'</row-3>')
         WRITE (fileNum,270) a3Temp(1),a3Temp(2),a3Temp(3)

         WRITE (fileNum,'(a)') '         </bravaisMatrix>'
      END IF

      IF ((cell%latnam.EQ.'squ').OR.(cell%latnam.EQ.'hex').OR.&
     &    (cell%latnam.EQ.'c-b').OR.(cell%latnam.EQ.'hx3').OR.&
     &    (cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
         252 FORMAT('         <a1>',f0.12,'</a1>')
         WRITE (fileNum,252) a1Temp(1)

         IF ((cell%latnam.EQ.'c-r').OR.(cell%latnam.EQ.'p-r')) THEN
            262 FORMAT('         <a2>',f0.12,'</a2>')
            WRITE (fileNum,262) a2Temp(2)
         END IF

         272 FORMAT('         <c>',f0.12,'</c>')
         WRITE (fileNum,272) dtild
      END IF

      IF (cell%latnam.EQ.'obl') THEN
         254 FORMAT('         <row-1>',f0.12,' ',f0.12,'</row-1>')
         WRITE (fileNum,254) a1Temp(1), a1Temp(2)

         264 FORMAT('         <row-2>',f0.12,' ',f0.12,'</row-2>')
         WRITE (fileNum,264) a2Temp(1), a2Temp(2)

         274 FORMAT('         <c>',f0.12,'</c>')
         WRITE (fileNum,274) dtild
      END IF

      WRITE (fileNum,'(a)') '      </bulkLattice>'
   END IF
   WRITE (fileNum,'(a)') '   </cell>'

!   <xcFunctional name="pbe" relativisticCorrections="F">
   280 FORMAT('   <xcFunctional name="',a,'" relativisticCorrections="',l1,'"/>')
   WRITE (fileNum,280) TRIM(namex), l_relcor

!      <xcParams igrd="1" lwb="F" ndvgrd="6" idsprs="0" chng="-0.100e-11"/>

!   290 FORMAT('      <xcParams igrd="',i0,'" lwb="',l1,'" ndvgrd="',i0,'" idsprs="',i0,'" chng="',e,'"/>')
!   WRITE (fileNum,290) xcpot%igrd,obsolete%lwb,obsolete%ndvgrd,0,obsolete%chng
!   WRITE (fileNum,'(a)') '   </xcFunctional>'

   WRITE (fileNum,'(a)') '   <atomSpecies>'
   DO iSpecies=1, numSpecies
      iAtomType = speciesRepAtomType(iSpecies)
      IF(iAtomType.EQ.-1) THEN
         EXIT
      END IF
!      <species name="Si-1" element="Si" atomicNumber="14" coreStates="4" magMom="0.0" flipSpin="F">
      300 FORMAT('      <species name="',a,'" element="',a,'" atomicNumber="',i0,'" coreStates="',i0,'" magMom="',f0.8,'" flipSpin="',l1,'">')
      tempNumberString = ''
      WRITE(tempNumberString,'(i0)') iSpecies
      speciesName = TRIM(ADJUSTL(noel(iAtomType))) // '-' // TRIM(ADJUSTL(tempNumberString))
      WRITE (fileNum,300) TRIM(ADJUSTL(speciesName)),TRIM(ADJUSTL(noel(iAtomType))),atoms%nz(iAtomType),atoms%ncst(iAtomType),atoms%bmu(iAtomType),atoms%nflip(iAtomType)

!         <mtSphere radius="2.160000" gridPoints="521" logIncrement="0.022000"/>
      310 FORMAT('         <mtSphere radius="',f0.8,'" gridPoints="',i0,'" logIncrement="',f0.8,'"/>')
      WRITE (fileNum,310) atoms%rmt(iAtomType),atoms%jri(iAtomType),atoms%dx(iAtomType)

!         <atomicCutoffs lmax="8" lnonsphr="6"/>
      320 FORMAT('         <atomicCutoffs lmax="',i0,'" lnonsphr="',i0,'"/>')
      WRITE (fileNum,320) atoms%lmax(iAtomType),atoms%lnonsph(iAtomType)

      IF (ALL((enpara%el0(0:3,iAtomType,1)-INT(enpara%el0(0:3,iAtomType,1))).LE.0.00000001)) THEN
!         <energyParameters s="3" p="3" d="3" f="4"/>
         321 FORMAT('         <energyParameters s="',i0,'" p="',i0,'" d="',i0,'" f="',i0,'"/>')
         WRITE (fileNum,321) INT(enpara%el0(0,iAtomType,1)),INT(enpara%el0(1,iAtomType,1)),&
                             INT(enpara%el0(2,iAtomType,1)),INT(enpara%el0(3,iAtomType,1))
      END IF

      IF(ANY(xmlElectronStates(:,iAtomType).NE.noState_const)) THEN
         endCoreStates = 1
         startCoreStates = 1
         coreStatesString = ''
         valenceStatesString = ''
         DO i = 1, 29
            IF (xmlElectronStates(i,iAtomType).EQ.coreState_const) endCoreStates = i
         END DO
         IF ((endCoreStates.GE.24).AND.&
&            (ALL(xmlPrintCoreStates(1:24,iAtomType).EQV..FALSE.)).AND.&
&            (ALL(xmlElectronStates(1:24,iAtomType).EQ.coreState_const)) ) THEN
            coreStatesString = nobleGasConfigList(6)
            startCoreStates = 25
         ELSE IF ((endCoreStates.GE.17).AND.&
&                 (ALL(xmlPrintCoreStates(1:17,iAtomType).EQV..FALSE.)).AND.&
&                 (ALL(xmlElectronStates(1:17,iAtomType).EQ.coreState_const))) THEN
            coreStatesString = nobleGasConfigList(5)
            startCoreStates = 18
         ELSE IF ((endCoreStates.GE.12).AND.&
&                 (ALL(xmlPrintCoreStates(1:12,iAtomType).EQV..FALSE.)).AND.&
&                 (ALL(xmlElectronStates(1:12,iAtomType).EQ.coreState_const))) THEN
            coreStatesString = nobleGasConfigList(4)
            startCoreStates = 13
         ELSE IF ((endCoreStates.GE.7).AND.&
&                 (ALL(xmlPrintCoreStates(1:7,iAtomType).EQV..FALSE.)).AND.&
&                 (ALL(xmlElectronStates(1:7,iAtomType).EQ.coreState_const))) THEN
            coreStatesString = nobleGasConfigList(3)
            startCoreStates = 8
         ELSE IF ((endCoreStates.GE.4).AND.&
&                 (ALL(xmlPrintCoreStates(1:4,iAtomType).EQV..FALSE.)).AND.&
&                 (ALL(xmlElectronStates(1:4,iAtomType).EQ.coreState_const))) THEN
            coreStatesString = nobleGasConfigList(2)
            startCoreStates = 5
         ELSE IF ((endCoreStates.GE.1).AND.&
&                 (ALL(xmlPrintCoreStates(1:1,iAtomType).EQV..FALSE.)).AND.&
&                 (ALL(xmlElectronStates(1:1,iAtomType).EQ.coreState_const))) THEN
            coreStatesString = nobleGasConfigList(1)
            startCoreStates = 2
         END IF
         DO i = startCoreStates, endCoreStates
            IF(xmlElectronStates(i,iAtomType).EQ.coreState_const) THEN
               coreStatesString = TRIM(ADJUSTL(coreStatesString)) // ' ' // coreStateList(i)
            END IF
         END DO
         DO i = 1, 29
            IF(xmlElectronStates(i,iAtomType).EQ.valenceState_const) THEN
               valenceStatesString = TRIM(ADJUSTL(valenceStatesString)) // ' ' // coreStateList(i)
            END IF
         END DO
         WRITE (fileNum,'(a)') '         <electronConfig>'
!         <coreConfig>[He] (2s1/2) (2p1/2) (2p3/2)</coreConfig>
         322 FORMAT('            <coreConfig>',a,'</coreConfig>')
         WRITE(fileNum,322) TRIM(ADJUSTL(coreStatesString))
         323 FORMAT('            <valenceConfig>',a,'</valenceConfig>')
         WRITE(fileNum,323) TRIM(ADJUSTL(valenceStatesString))
         DO i = startCoreStates, 29
            IF ((xmlElectronStates(i,iAtomType).NE.noState_const).AND.(xmlPrintCoreStates(i,iAtomType))) THEN
!         <coreStateOccupation state="(2s1/2)" spinUp="1.0" spinDown="1.0"/>
               325 FORMAT('            <stateOccupation state="',a,'" spinUp="',f0.8,'" spinDown="',f0.8,'"/>')
               WRITE(fileNum,325) coreStateList(i), xmlCoreOccs(1,i,iAtomType), xmlCoreOccs(2,i,iAtomType)
            END IF
         END DO
         WRITE (fileNum,'(a)') '         </electronConfig>'
      END IF

      DO ilo = 1, atoms%nlo(iAtomType)
!         <lo type="HELO" l="0" n="4"/>
         l = atoms%llo(ilo,iAtomType)
         n = INT(enpara%ello0(ilo,iAtomType,1))
         loType = 'SCLO'
         IF(n.LT.0) THEN
            loType = 'HELO'
         END IF
         n = ABS(n)
         324 FORMAT('         <lo type="',a,'" l="',i0,'" n="',i0,'" eDeriv="',i0,'"/>')
         WRITE (fileNum,324) TRIM(ADJUSTL(loType)), l, n, atoms%ulo_der(ilo,iAtomType)
      END DO

      WRITE (fileNum,'(a)') '      </species>'
   END DO
   WRITE (fileNum,'(a)') '   </atomSpecies>'
   WRITE (fileNum,'(a)') '   <atomGroups>'
   na = 0
   DO iAtomType=1, atoms%ntype
      iSpecies = atomTypeSpecies(iAtomType)
!      <atomGroup species="Si-1">
      330 FORMAT('      <atomGroup species="',a,'">')
      tempNumberString = ''
      WRITE(tempNumberString,'(i0)') iSpecies
      speciesName = TRIM(ADJUSTL(noel(iAtomType))) // '-' // TRIM(ADJUSTL(tempNumberString))
      WRITE (fileNum,330) TRIM(ADJUSTL(speciesName))

      DO ieq=1,atoms%neq(iAtomType)
         na = na + 1
         tempTaual(1,na) = atoms%taual(1,na)
         tempTaual(2,na) = atoms%taual(2,na)
         tempTaual(3,na) = atoms%taual(3,na)
         DO i = 2,9
            rest = ABS(i*tempTaual(1,na) - NINT(i*tempTaual(1,na)) ) + ABS(i*tempTaual(2,na) - NINT(i*tempTaual(2,na)))
            IF (.not.input%film) THEN
               rest = rest + ABS(i*tempTaual(3,na) - NINT(i*tempTaual(3,na)) )
            END IF
            IF (rest.LT.(i*0.000001)) EXIT
         END DO
         scpos = 1.0
         IF (i.LT.10) scpos = real(i)  ! common factor found (x,y)
         DO i = 1,2
            tempTaual(i,na) = tempTaual(i,na)*scpos
         ENDDO
         IF (.not.input%film) tempTaual(3,na) = tempTaual(3,na)*scpos
         IF (input%film) THEN
            tempTaual(3,na) = dtild*tempTaual(3,na)/scale
         END IF
!+odim in 1D case all the coordinates are given in cartesian YM
         IF (oneD%odd%d1) THEN
            tempTaual(1,na) = tempTaual(1,na)*a1(1)
            tempTaual(2,na) = tempTaual(2,na)*a2(2)
         END IF
!-odim
         IF (oneD%odd%d1) THEN
            STOP '1D position output not implemented!'
         ELSE IF (input%film) THEN
!         <filmPos> x/myConstant  y/myConstant  1/myConstant</filmPos>
            340 FORMAT('         <filmPos>',a,' ',a,' ',a,'</filmPos>')
            xPosString = ''
            yPosString = ''
            zPosString = ''
            IF((scpos.NE.1.0).AND.((tempTaual(1,na).NE.0.0).OR.(tempTaual(2,na).NE.0.0).OR.(tempTaual(3,na).NE.0.0))) THEN
               WRITE(xPosString,'(f0.12,a1,f0.12)') tempTaual(1,na), '/', scpos
               WRITE(yPosString,'(f0.12,a1,f0.12)') tempTaual(2,na), '/', scpos
            ELSE
               WRITE(xPosString,'(f0.12)') tempTaual(1,na)
               WRITE(yPosString,'(f0.12)') tempTaual(2,na)
            END IF
            WRITE(zPosString,'(f0.12)') tempTaual(3,na)
            WRITE (fileNum,340) TRIM(ADJUSTL(xPosString)),TRIM(ADJUSTL(yPosString)),TRIM(ADJUSTL(zPosString))
         ELSE
!         <relPos> x/myConstant  y/myConstant  z/myConstant</relPos>
            350 FORMAT('         <relPos>',a,' ',a,' ',a,'</relPos>')
            xPosString = ''
            yPosString = ''
            zPosString = ''
            IF((scpos.NE.1.0).AND.((tempTaual(1,na).NE.0.0).OR.(tempTaual(2,na).NE.0.0).OR.(tempTaual(3,na).NE.0.0))) THEN
               WRITE(xPosString,'(f0.12,a1,f0.12)') tempTaual(1,na), '/', scpos
               WRITE(yPosString,'(f0.12,a1,f0.12)') tempTaual(2,na), '/', scpos
               WRITE(zPosString,'(f0.12,a1,f0.12)') tempTaual(3,na), '/', scpos
            ELSE
               WRITE(xPosString,'(f0.12)') tempTaual(1,na)
               WRITE(yPosString,'(f0.12)') tempTaual(2,na)
               WRITE(zPosString,'(f0.12)') tempTaual(3,na)
            END IF
            WRITE (fileNum,350) TRIM(ADJUSTL(xPosString)),TRIM(ADJUSTL(yPosString)),TRIM(ADJUSTL(zPosString))
         END IF
      END DO
!         <force calculate="F" relaxX="T" relaxY="T" relaxZ="T"/>
      360 FORMAT('         <force calculate="',l1,'" relaxXYZ="',3l1,'"/>')
      WRITE (fileNum,360) atoms%l_geo(iAtomType),atoms%relax(1,iAtomType),atoms%relax(2,iAtomType),atoms%relax(3,iAtomType)

      WRITE (fileNum,'(a)') '      </atomGroup>'
   END DO
   WRITE (fileNum,'(a)') '   </atomGroups>'

   368 FORMAT('   <output dos="',l1,'" band="',l1,'" vacdos="',l1,'" slice="',l1,'">')
   WRITE (fileNum,368) banddos%dos,band,banddos%vacdos,sliceplot%slice

!      <checks vchk="F" cdinf="F" disp="F"/>
   370 FORMAT('      <checks vchk="',l1,'" cdinf="',l1,'" disp="',l1,'"/>')
   WRITE (fileNum,370) input%vchk,input%cdinf,obsolete%disp

!      <densityOfStates ndir="0" minEnergy="-0.50000" maxEnergy="0.50000" sigma="0.01500"/>  
   380 FORMAT('      <densityOfStates ndir="',i0,'" minEnergy="',f0.8,'" maxEnergy="',f0.8,'" sigma="',f0.8,'"/>')
   WRITE (fileNum,380) banddos%ndir,banddos%e2_dos,banddos%e1_dos,banddos%sig_dos

!      <vacuumDOS layers="0" integ="F" star="F" nstars="0" locx1="0.00" locy1="0.00" locx2="0.00" locy2="0.00" nstm="0" tworkf="0.000000"/>
   390 FORMAT('      <vacuumDOS layers="',i0,'" integ="',l1,'" star="',l1,'" nstars="',i0,'" locx1="',f0.8,'" locy1="',f0.8,'" locx2="',f0.8,'" locy2="',f0.8,'" nstm="',i0,'" tworkf="',f0.8,'"/>')
   WRITE (fileNum,390) vacuum%layers,input%integ,vacuum%starcoeff,vacuum%nstars,vacuum%locx(1),vacuum%locy(1),vacuum%locx(2),vacuum%locy(2),vacuum%nstm,vacuum%tworkf

!      <plotting iplot="F" score="F" plplot="F"/>
   400 FORMAT('      <plotting iplot="',l1,'" score="',l1,'" plplot="',l1,'"/>')
   WRITE (fileNum,400) sliceplot%iplot,input%score,sliceplot%plpot

!      <chargeDensitySlicing numkpt="0" minEigenval="0.000000" maxEigenval="0.000000" nnne="0" pallst="F"/>
   410 FORMAT('      <chargeDensitySlicing numkpt="',i0,'" minEigenval="',f0.8,'" maxEigenval="',f0.8,'" nnne="',i0,'" pallst="',l1,'"/>')
   WRITE (fileNum,410) sliceplot%kk,sliceplot%e1s,sliceplot%e2s,sliceplot%nnne,input%pallst

!      <specialOutput form66="F" eonly="F" bmt="F"/>
   420 FORMAT('      <specialOutput form66="',l1,'" eonly="',l1,'" bmt="',l1,'"/>')
   WRITE (fileNum,420) obsolete%form66,input%eonly,input%l_bmt

   WRITE (fileNum,'(a)') '   </output>'
   IF(.NOT.l_outFile) THEN
      WRITE (fileNum,'(a)') '</fleurInput>'
      CLOSE (fileNum)
   END IF

END SUBROUTINE w_inpXML
END MODULE m_winpXML
