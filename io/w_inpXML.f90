!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

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
   SUBROUTINE w_inpXML( &
      atoms, vacuum, input, stars, sliceplot, forcetheo, banddos, &
      cell, sym, xcpot, noco, oneD, hybrid, kpts, enpara, &
      l_explicitIn, l_includeIn, filename)

      use m_types_input
      use m_types_sym
      use m_types_stars
      use m_types_atoms
      use m_types_vacuum
      use m_types_kpts
      use m_types_oneD
      use m_types_hybrid
      use m_types_cell
      use m_types_banddos
      use m_types_sliceplot
      USE m_types_xcpot
      USE m_types_xcpot_inbuild_nofunction
      USE m_types_noco
      use m_types_enparaxml
      USE m_types_forcetheo

      USE m_juDFT
      USE m_constants
      USE m_xmlOutput

      IMPLICIT NONE

! arguments

      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_sym), INTENT(IN)     :: sym
      TYPE(t_stars), INTENT(IN)   :: stars
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_vacuum), INTENT(IN)   :: vacuum
      TYPE(t_kpts), INTENT(IN)     :: kpts
      TYPE(t_oneD), INTENT(IN)     :: oneD
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_cell), INTENT(IN)     :: cell
      TYPE(t_banddos), INTENT(IN)  :: banddos
      TYPE(t_sliceplot), INTENT(IN):: sliceplot
      CLASS(t_xcpot), INTENT(IN)   :: xcpot
      TYPE(t_noco), INTENT(IN)     :: noco
      CLASS(t_enparaxml), INTENT(IN)   :: enpara
      CLASS(t_forcetheo), INTENT(IN):: forcetheo !nothing is done here so far....
      LOGICAL, INTENT(IN)        :: l_explicitIn, l_includeIn(4)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

      INTEGER          :: iSpecies, fileNum
      CHARACTER(len=8) :: name(10)
      INTEGER          :: numSpecies
      INTEGER          :: speciesRepAtomType(atoms%ntype)
      CHARACTER(len=20):: speciesNames(atoms%ntype)
      LOGICAL          :: known_species

!+lda+u
      REAL u, j
      INTEGER l, i_u
      INTEGER uIndices(2, atoms%ntype)
      LOGICAL l_amf
      CHARACTER(len=3) ch_test
      NAMELIST /ldaU/ l, u, j, l_amf
!-lda+u
!+odim
      INTEGER MM, vM, m_cyl
      LOGICAL invs1, zrfs1
      INTEGER chi, rot
      LOGICAL d1, band
      NAMELIST /odim/ d1, MM, vM, m_cyl, chi, rot, invs1, zrfs1
!-odim
! ..
! ..  Local Variables
      REAL     :: zc, sumWeight, occ(2)
      INTEGER  ::nw, idsprs, n1, n2
      INTEGER ieq, i, k, na, n, ilo
      REAL s3, ah, a, hs2, rest
      LOGICAL l_hyb, ldum
      INTEGER :: ierr
! ..
!...  Local Arrays
      CHARACTER :: helpchar(atoms%ntype)
      CHARACTER(len=4) :: chntype
      CHARACTER(len=41) :: chform
      CHARACTER(len=100) :: line

!     added for HF and hybrid functionals
      REAL                  ::  aMix, omega
      INTEGER               :: idum
      CHARACTER(len=1)     ::  check

      CHARACTER(len=20) :: speciesName
      CHARACTER(len=150) :: format
      CHARACTER(len=20) :: mixingScheme
      CHARACTER(len=10) :: loType
      CHARACTER(len=10) :: bzIntMode
      LOGICAL ::   l_explicit, l_nocoOpt, l_include(4)
      INTEGER :: iAtomType, startCoreStates, endCoreStates
      CHARACTER(len=100) :: posString(3)
      CHARACTER(len=7) :: str
      REAL :: tempTaual(3, atoms%nat), scpos(3)
      REAL :: amatTemp(3, 3), bmatTemp(3, 3)

      l_include = l_includeIn .or. .not. present(filename)
      l_explicit = l_explicitIn .OR. .not. present(filename)
      l_nocoOpt = noco%l_noco .OR. juDFT_was_argument("-noco")

      band = .false.
      nw = 1

      IF (PRESENT(filename)) THEN
         filenum = 98
         OPEN (fileNum, file=TRIM(ADJUSTL(filename)), form='formatted', status='replace')
         WRITE (fileNum, '(a)') '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
         WRITE (fileNum, '(a)') '<fleurInput fleurInputVersion="0.30">'
      ELSE
         fileNum = getXMLOutputUnitNumber()
         CALL openXMLElementNoAttributes('inputData')
      END IF

      WRITE (fileNum, '(a)') '   <comment>'
      WRITE (fileNum, '(a6,10a8)') '      ', input%comment
      WRITE (fileNum, '(a)') '   </comment>'

      WRITE (fileNum, '(a)') '   <calculationSetup>'

!      <cutoffs Kmax="3.60000" Gmax="11.000000" GmaxXC="9.200000" numbands="0"/>
110   FORMAT('      <cutoffs Kmax="', f0.8, '" Gmax="', f0.8, '" GmaxXC="', f0.8, '" numbands="', i0, '"/>')
      WRITE (fileNum, 110) input%rkmax, input%gmax, xcpot%gmaxxc, input%gw_neigd

!      <scfLoop itmax="9" maxIterBroyd="99" imix="Anderson" alpha="0.05" precondParam="0.0" spinf="2.00"/>
120   FORMAT('      <scfLoop itmax="', i0, '" minDistance="', f0.8, '" maxIterBroyd="', i0, '" imix="', a, '" alpha="', f0.8, '" precondParam="', f3.1, '" spinf="', f0.8, '"/>')
      SELECT CASE (input%imix)
      CASE (1)
         mixingScheme = 'straight'
      CASE (3)
         mixingScheme = 'Broyden1'
      CASE (5)
         mixingScheme = 'Broyden2'
      CASE (7)
         mixingScheme = 'Anderson'
      CASE DEFAULT
         mixingScheme = 'errorUnknownMixing'
      END SELECT
      WRITE (fileNum, 120) input%itmax, input%minDistance, input%maxiter, TRIM(mixingScheme), input%alpha, input%preconditioning_param, input%spinf

!      <coreElectrons ctail="T" frcor="F" kcrel="0"/>
130   FORMAT('      <coreElectrons ctail="', l1, '" frcor="', l1, '" kcrel="', i0, '" coretail_lmax="', i0, '"/>')
      WRITE (fileNum, 130) input%ctail, input%frcor, input%kcrel, input%coretail_lmax

!      <magnetism jspins="1" l_noco="F" l_J="F" swsp="F" lflip="F"/>
140   FORMAT('      <magnetism jspins="', i0, '" l_noco="', l1, '" swsp="', l1, '" lflip="', l1, '"/>')
      WRITE (fileNum, 140) input%jspins, noco%l_noco, input%swsp, input%lflip

      !      <soc theta="0.00000" phi="0.00000" l_soc="F" spav="F" off="F" soc66="F"/>
150   FORMAT('      <soc theta="', f0.8, '" phi="', f0.8, '" l_soc="', l1, '" spav="', l1, '"/>')
      WRITE (fileNum, 150) noco%theta, noco%phi, noco%l_soc, noco%l_spav

      IF (l_explicit .OR. hybrid%l_hybrid) THEN
155      FORMAT('      <prodBasis gcutm="', f0.8, '" tolerance="', f0.8, '" ewaldlambda="', i0, '" lexp="', i0, '" bands="', i0, '"/>')
         WRITE (fileNum, 155) hybrid%g_cutoff, hybrid%linear_dep_tol, hybrid%ewaldlambda, hybrid%lexp, hybrid%bands1
      END IF

      IF (l_nocoOpt .OR. l_explicit) THEN
160      FORMAT('      <nocoParams l_ss="', l1, '" l_mperp="', l1, '" l_constr="', l1, &
                '" mix_b="', f0.8, '">')
         WRITE (fileNum, 160) noco%l_ss, noco%l_mperp, noco%l_constr, noco%mix_b
162      FORMAT('         <qss>', f0.10, ' ', f0.10, ' ', f0.10, '</qss>')
         WRITE (fileNum, 162) noco%qss(1), noco%qss(2), noco%qss(3)
         WRITE (fileNum, '(a)') '      </nocoParams>'
      END IF

      IF (oneD%odd%d1) THEN
170      FORMAT('      <oneDParams d1="', l1, '" MM="', i0, '" vM="', i0, '" m_cyl="', i0, '" chi="', i0, '" rot="', i0, '" invs1="', l1, '" zrfs1="', l1, '"/>')
         WRITE (fileNum, 170) oneD%odd%d1, oneD%odd%M, oneD%odd%mb, oneD%odd%m_cyl, oneD%odd%chi, oneD%odd%rot, oneD%odd%invs, oneD%odd%zrfs
      END IF

!      <expertModes gw="0"  eig66="F" lpr="0" secvar="F" />
180   FORMAT('      <expertModes gw="', i0, '" secvar="', l1, '"/>')
      WRITE (fileNum, 180) input%gw, input%secvar

!      <geometryOptimization l_f="F" xa="2.00000" thetad="330.00000" epsdisp="0.00001" epsforce="0.00001"/>
190   FORMAT('      <geometryOptimization l_f="', l1, '" forcealpha="', f0.8, '" forcemix="', a, '" epsdisp="', f0.8, '" epsforce="', f0.8, '"/>')
      SELECT CASE (input%forcemix)
      CASE (0)
         mixingScheme = 'Straight'
      CASE (1)
         mixingScheme = 'CG'
      CASE (2)
         mixingScheme = 'BFGS'
      CASE DEFAULT
         mixingScheme = 'errorUnknownMixing'
      END SELECT
      WRITE (fileNum, 190) input%l_f, input%forcealpha, TRIM(mixingScheme), input%epsdisp, input%epsforce

      IF (input%gauss .AND. input%tria) THEN
         call juDFT_error('Error: bz integration modes gauss AND tria selected!')
      END IF

      bzIntMode = 'hist'
      IF (input%gauss) THEN
         bzIntMode = 'gauss'
      ELSE IF (input%tria) THEN
         bzIntMode = 'tria'
      END IF

!      <ldaU l_linMix="F" mixParam="0.05" spinf="1.0" />
195   FORMAT('      <ldaU l_linMix="', l1, '" mixParam="', f0.6, '" spinf="', f0.6, '"/>')
      WRITE (fileNum, 195) input%ldauLinMix, input%ldauMixParam, input%ldauSpinf

!      <bzIntegration valenceElectrons="8.00000" mode="hist" fermiSmearingEnergy="0.00100">
200   FORMAT('      <bzIntegration valenceElectrons="', f0.8, '" mode="', a, '" fermiSmearingEnergy="', f0.8, '">')
      WRITE (fileNum, 200) input%zelec, TRIM(ADJUSTL(bzIntMode)), input%tkb

      if (l_include(1)) THEN
         call kpts%print_xml(fileNum)
      else
         WRITE (fileNum, '(a)') '         <!-- k-points included here -->'
         WRITE (fileNum, '(a)') '         <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="kpts.xml"> </xi:include>'
      end if
      WRITE (fileNum, '(a)') '      </bzIntegration>'

!      <energyParameterLimits ellow="-2.00000" elup="2.00000"/>
220   FORMAT('      <energyParameterLimits ellow="', f0.8, '" elup="', f0.8, '"/>')
      WRITE (fileNum, 220) input%ellow, input%elup
      if (l_include(2)) THEN
         call sym%print_xml(fileNum)
      else
         WRITE (fileNum, '(a)') '      <!-- symmetry operations included here -->'
         WRITE (fileNum, '(a)') '      <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="sym.xml"> </xi:include>'
      end if
      WRITE (fileNum, '(a)') '   </calculationSetup>'
      WRITE (fileNum, '(a)') '   <cell>'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note: Different options for the cell definition!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (input%film) THEN
!      <xsd:attribute name="dVac" type="xsd:double" use="required"/>
!      <xsd:attribute name="dTilda" type="xsd:double" use="required"/>
!      <filmLattice ...>
241      FORMAT('      <filmLattice scale="', f0.8, '" dVac="', f0.8, '" dTilda="', f0.8, '">')
         WRITE (fileNum, 241) 1.0, vacuum%dvac, cell%amat(3, 3)
      ELSE
242      FORMAT('      <bulkLattice scale="', f0.10, '">')
         WRITE (fileNum, 242) 1.0
      ENDIF
      !         <bravaisMatrix>
      WRITE (fileNum, '(a)') '         <bravaisMatrix>'
      !            <row-1>0.00000 5.13000 5.13000</row-1>
250   FORMAT('            <row-1>', f0.10, ' ', f0.10, ' ', f0.10, '</row-1>')
      WRITE (fileNum, 250) cell%amat(:, 1)
      !            <row-2>5.13000 0.00000 5.13000</row-2>
260   FORMAT('            <row-2>', f0.10, ' ', f0.10, ' ', f0.10, '</row-2>')
      WRITE (fileNum, 260) cell%amat(:, 2)
      !            <row-3>5.13000 5.13000 0.00000</row-3>
270   FORMAT('            <row-3>', f0.10, ' ', f0.10, ' ', f0.10, '</row-3>')
      WRITE (fileNum, 270) cell%amat(:, 3)
      WRITE (fileNum, '(a)') '         </bravaisMatrix>'

      IF (input%film) THEN
268      FORMAT('         <vacuumEnergyParameters vacuum="', i0, '" spinUp="', f0.8, '" spinDown="', f0.8, '"/>')
         DO i = 1, vacuum%nvac
            WRITE (fileNum, 268) i, enpara%evac0(i, 1), enpara%evac0(i, input%jspins)
         END DO

         WRITE (fileNum, '(a)') '      </filmLattice>'
      ELSE
         WRITE (fileNum, '(a)') '      </bulkLattice>'
      END IF
      WRITE (fileNum, '(a)') '   </cell>'

      SELECT TYPE (xcpot)
      CLASS IS (t_xcpot_inbuild_nf)
         !   <xcFunctional name="pbe" relativisticCorrections="F">
280      FORMAT('   <xcFunctional name="', a, '" relativisticCorrections="', l1, '"/>')
         WRITE (fileNum, 280) trim(xcpot%get_name()), xcpot%relativistic_correction()
      END SELECT

      uIndices = -1
      DO i_u = 1, atoms%n_u
         IF (uIndices(1, atoms%lda_u(i_u)%atomType) .EQ. -1) uIndices(1, atoms%lda_u(i_u)%atomType) = i_u
         uIndices(2, atoms%lda_u(i_u)%atomType) = i_u
      END DO

      !Build list of species
      speciesNames = ''
      numSpecies = 0
      DO n = 1, atoms%ntype
         known_species = ANY(atoms%speciesname(n) == speciesNames(:numSpecies))
         if (.not. known_species) THEN
            numSpecies = numSpecies + 1
            speciesNames(numSpecies) = atoms%speciesname(n)
            speciesRepAtomType(numSpecies) = n
         end if
      enddo

      if (.not. l_include(3)) then
         open (99, file='species.xml')
         WRITE (fileNum, '(a)') '      <!-- species included here -->'
         WRITE (fileNum, '(a)') '      <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="species.xml"> </xi:include>'
         fileNum = 99
      endif

      WRITE (fileNum, '(a)') '   <atomSpecies>'
      DO iSpecies = 1, numSpecies
         iAtomType = speciesRepAtomType(iSpecies)
         IF (iAtomType .EQ. -1) THEN
            EXIT
         END IF
!      <species name="Si-1" element="Si" atomicNumber="14" coreStates="4" magMom="0.0" flipSpin="F">
300      FORMAT('      <species name="', a, '" element="', a, '" atomicNumber="', i0, '" flipSpinPhi="', f0.8, '" flipSpinTheta="', f0.8, '" flipSpinScale="', l1, '">')
         speciesName = TRIM(ADJUSTL(atoms%speciesName(iSpecies)))
         WRITE (fileNum, 300) TRIM(ADJUSTL(speciesName)), TRIM(ADJUSTL(namat_const(atoms%nz(iAtomType)))), atoms%nz(iAtomType), atoms%flipSpinPhi(iAtomType), atoms%flipSpinTheta(iAtomType), atoms%flipSpinScale(iAtomType)

!         <mtSphere radius="2.160000" gridPoints="521" logIncrement="0.022000"/>
310      FORMAT('         <mtSphere radius="', f0.8, '" gridPoints="', i0, '" logIncrement="', f0.8, '"/>')
         WRITE (fileNum, 310) atoms%rmt(iAtomType), atoms%jri(iAtomType), atoms%dx(iAtomType)

!         <atomicCutoffs lmax="8" lnonsphr="6"/>
320      FORMAT('         <atomicCutoffs lmax="', i0, '" lnonsphr="', i0, '"/>')
         WRITE (fileNum, 320) atoms%lmax(iAtomType), atoms%lnonsph(iAtomType)

         WRITE (fileNum, '(a)') '         <electronConfig>'
!         <coreConfig>[He] (2s1/2) (2p1/2) (2p3/2)</coreConfig>
322      FORMAT('            <coreConfig>', a, '</coreConfig>')
         WRITE (fileNum, 322) TRIM(ADJUSTL(atoms%econf(iAtomType)%coreconfig))
323      FORMAT('            <valenceConfig>', a, '</valenceConfig>')
         IF (len_TRIM(atoms%econf(iAtomType)%valenceconfig) > 1) THEN
            WRITE (fileNum, 323) TRIM(ADJUSTL(atoms%econf(iAtomType)%valenceconfig))
         END IF
         DO i = 1, MERGE(atoms%econf(iAtomType)%num_states, atoms%econf(iAtomType)%num_core_states, (len_TRIM(atoms%econf(iAtomType)%valenceconfig) > 1))
            occ = atoms%econf(iAtomType)%occupation(i, :)
            IF (ABS(occ(1) - occ(2)) > 1E-5 .OR. ABS(occ(1) - ABS(atoms%econf(iAtomType)%kappa(i))) > 1E-5) THEN
               !State not fully occupied
325            FORMAT('            <stateOccupation state="', a, '" spinUp="', f0.8, '" spinDown="', f0.8, '"/>')
               str = atoms%econf(iAtomType)%get_state_string(i)
               WRITE (fileNum, 325) str, occ(1), occ(2)
            END IF
         END DO
         WRITE (fileNum, '(a)') '         </electronConfig>'

         !IF (ALL(enpara%qn_el(0:3,iAtomType,1).ne.0)) THEN
!!         <energyParameters s="3" p="3" d="3" f="4"/>
321      FORMAT('         <energyParameters s="', i0, '" p="', i0, '" d="', i0, '" f="', i0, '"/>')
         WRITE (fileNum, 321) enpara%qn_el(0:3, iAtomType, 1)
         !END IF
         IF (l_explicit .OR. hybrid%l_hybrid) THEN
315         FORMAT('         <prodBasis lcutm="', i0, '" lcutwf="', i0, '" select="', a, '"/>')
            line = ''
            WRITE (line, '(i0,1x,i0,1x,i0,1x,i0)') hybrid%select1(1:4, iAtomType)
            WRITE (fileNum, 315) hybrid%lcutm1(iAtomType), hybrid%lcutwf(iAtomType), TRIM(ADJUSTL(line))
         END IF

         IF (uIndices(1, iAtomType) .NE. -1) THEN
!         <ldaU l="2" U="5.5" J="0.9" l_amf="F"/>
            DO i_u = uIndices(1, iAtomType), uIndices(2, iAtomType)
326            FORMAT('         <ldaU l="', i0, '" U="', f0.5, '" J="', f0.5, '" l_amf="', l1, '"/>')
               WRITE (fileNum, 326) atoms%lda_u(i_u)%l, atoms%lda_u(i_u)%u, atoms%lda_u(i_u)%j, atoms%lda_u(i_u)%l_amf
            END DO
         END IF

         DO ilo = 1, atoms%nlo(iAtomType)
!         <lo type="HELO" l="0" n="4"/>
            l = atoms%llo(ilo, iAtomType)
            n = enpara%qn_ello(ilo, iAtomType, 1)
            loType = 'SCLO'
            IF (n .LT. 0) THEN
               loType = 'HELO'
            END IF
            n = ABS(n)
324         FORMAT('         <lo type="', a, '" l="', i0, '" n="', i0, '" eDeriv="', i0, '"/>')
            WRITE (fileNum, 324) TRIM(ADJUSTL(loType)), l, n, atoms%ulo_der(ilo, iAtomType)
         END DO

         WRITE (fileNum, '(a)') '      </species>'
      END DO
      WRITE (fileNum, '(a)') '   </atomSpecies>'

      if (.not. l_include(3)) then
         close (99)
         fileNum = 98
      endif

      if (.not. l_include(4)) then
         open (99, file='atoms.xml')
         WRITE (fileNum, '(a)') '      <!-- atoms group included here -->'
         WRITE (fileNum, '(a)') '      <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="atoms.xml"> </xi:include>'
         fileNum = 98
      endif

      WRITE (fileNum, '(a)') '   <atomGroups>'
      na = 0
      DO iAtomType = 1, atoms%ntype
!      <atomGroup species="Si-1">
330      FORMAT('      <atomGroup species="', a, '">')
         speciesName = TRIM(ADJUSTL(atoms%speciesName(iAtomType)))
         WRITE (fileNum, 330) TRIM(ADJUSTL(speciesName))

         DO ieq = 1, atoms%neq(iAtomType)
            na = na + 1
            tempTaual(1, na) = atoms%taual(1, na)
            tempTaual(2, na) = atoms%taual(2, na)
            tempTaual(3, na) = atoms%taual(3, na)
            scpos = 1.0
            DO i = 2, 40
               rest = ABS(i*tempTaual(1, na) - NINT(i*tempTaual(1, na)))
               IF ((scpos(1) .EQ. 1.0) .AND. (rest .LT. (i*0.000001))) scpos(1) = real(i)
               rest = ABS(i*tempTaual(2, na) - NINT(i*tempTaual(2, na)))
               IF ((scpos(2) .EQ. 1.0) .AND. (rest .LT. (i*0.000001))) scpos(2) = real(i)
               IF (.not. input%film) THEN
                  rest = ABS(i*tempTaual(3, na) - NINT(i*tempTaual(3, na)))
                  IF ((scpos(3) .EQ. 1.0) .AND. (rest .LT. (i*0.000001))) scpos(3) = real(i)
               END IF
            END DO
            DO i = 1, 2
               tempTaual(i, na) = tempTaual(i, na)*scpos(i)
            END DO
            IF (.not. input%film) tempTaual(3, na) = tempTaual(3, na)*scpos(3)
            IF (input%film) THEN
               tempTaual(3, na) = cell%amat(3, 3)*tempTaual(3, na)
            END IF
!+odim in 1D case all the coordinates are given in cartesian YM
            IF (oneD%odd%d1) THEN
               !tempTaual(1,na) = tempTaual(1,na)*a1(1)
               !tempTaual(2,na) = tempTaual(2,na)*a2(2)
            END IF
!-odim
            IF (oneD%odd%d1) THEN
               call judft_error('1D position output not implemented!')
            ELSE IF (input%film) THEN
!         <filmPos> x/myConstant  y/myConstant  1/myConstant</filmPos>
340            FORMAT('         <filmPos label="', a20, '">', a, ' ', a, ' ', a, '</filmPos>')
               posString(:) = ''
               DO i = 1, 2
                  IF ((scpos(i) .NE. 1.0) .AND. (tempTaual(i, na) .NE. 0.0)) THEN
                     WRITE (posString(i), '(f0.3,a1,f0.3)') tempTaual(i, na), '/', scpos(i)
                  ELSE
                     WRITE (posString(i), '(f0.10)') tempTaual(i, na)
                  END IF
               END DO
               WRITE (posString(3), '(f0.10)') tempTaual(3, na)
               WRITE (fileNum, 340) TRIM(ADJUSTL(atoms%label(na))), &
                  TRIM(ADJUSTL(posString(1))), TRIM(ADJUSTL(posString(2))), TRIM(ADJUSTL(posString(3)))
            ELSE
!         <relPos> x/myConstant  y/myConstant  z/myConstant</relPos>
350            FORMAT('         <relPos label="', a20, '">', a, ' ', a, ' ', a, '</relPos>')
               posString(:) = ''
               DO i = 1, 3
                  IF ((scpos(i) .NE. 1.0) .AND. (tempTaual(i, na) .NE. 0.0)) THEN
                     WRITE (posString(i), '(f0.3,a1,f0.3)') tempTaual(i, na), '/', scpos(i)
                  ELSE
                     WRITE (posString(i), '(f0.10)') tempTaual(i, na)
                  END IF
               END DO
               WRITE (fileNum, 350) TRIM(ADJUSTL(atoms%label(na))), &
                  TRIM(ADJUSTL(posString(1))), TRIM(ADJUSTL(posString(2))), TRIM(ADJUSTL(posString(3)))
            END IF
         END DO
!         <force calculate="F" relaxX="T" relaxY="T" relaxZ="T"/>
360      FORMAT('         <force calculate="', l1, '" relaxXYZ="', 3l1, '"/>')
         WRITE (fileNum, 360) atoms%l_geo(iAtomType), atoms%relax(1, iAtomType), atoms%relax(2, iAtomType), atoms%relax(3, iAtomType)

         IF (l_nocoOpt .OR. l_explicit) THEN
362         FORMAT('         <nocoParams l_relax="', l1, '" alpha="', f0.8, '" beta="', &
                   f0.8, '" b_cons_x="', f0.8, '" b_cons_y="', f0.8, '"/>')
            WRITE (fileNum, 362) noco%l_relax(iAtomType), noco%alphInit(iAtomType), &
               noco%beta(iAtomType), noco%b_con(1, iAtomType), noco%b_con(2, iAtomType)
         END IF

         WRITE (fileNum, '(a)') '      </atomGroup>'
      END DO
      WRITE (fileNum, '(a)') '   </atomGroups>'
      if (.not. l_include(4)) then
         close (99)
         fileNum = 98
      endif

368   FORMAT('   <output dos="', l1, '" band="', l1, '" vacdos="', l1, '" slice="', l1, '" mcd="', l1, '">')
      WRITE (fileNum, 368) banddos%dos, band, banddos%vacdos, sliceplot%slice, banddos%l_mcd

!      <checks vchk="F" cdinf="F" disp="F"/>
370   FORMAT('      <checks vchk="', l1, '" cdinf="', l1, '"/>')
      WRITE (fileNum, 370) input%vchk, input%cdinf

!      <densityOfStates ndir="0" minEnergy="-0.50000" maxEnergy="0.50000" sigma="0.01500"/>
380   FORMAT('      <densityOfStates ndir="', i0, '" minEnergy="', f0.8, '" maxEnergy="', f0.8, '" sigma="', f0.8, '"/>')
      WRITE (fileNum, 380) banddos%ndir, banddos%e2_dos, banddos%e1_dos, banddos%sig_dos

!      <vacuumDOS layers="0" integ="F" star="F" nstars="0" locx1="0.00" locy1="0.00" locx2="0.00" locy2="0.00" nstm="0" tworkf="0.000000"/>
390   FORMAT('      <vacuumDOS layers="', i0, '" integ="', l1, '" star="', l1, '" nstars="', i0, '" locx1="', f0.5, '" locy1="', f0.5, '" locx2="', f0.5, '" locy2="', f0.5, '" nstm="', i0, '" tworkf="', f0.5, '"/>')
      WRITE (fileNum, 390) vacuum%layers, input%integ, vacuum%starcoeff, vacuum%nstars, vacuum%locx(1), vacuum%locy(1), vacuum%locx(2), vacuum%locy(2), vacuum%nstm, vacuum%tworkf

!      <unfoldingBand unfoldBand="F" supercellX="1" supercellY="1" supercellZ="1"/>
395   FORMAT('      <unfoldingBand unfoldBand="', l1, '" supercellX="', i0, '" supercellY="', i0, '" supercellZ="', i0, '"/>')
      WRITE (fileNum, 395) banddos%unfoldband, banddos%s_cell_x, banddos%s_cell_y, banddos%s_cell_z

!      <plotting iplot="0" />
400   FORMAT('      <plotting iplot="', i0, '"/>')
      WRITE (fileNum, 400) sliceplot%iplot

!      <chargeDensitySlicing numkpt="0" minEigenval="0.000000" maxEigenval="0.000000" nnne="0" pallst="F"/>
410   FORMAT('      <chargeDensitySlicing numkpt="', i0, '" minEigenval="', f0.8, '" maxEigenval="', f0.8, '" nnne="', i0, '" pallst="', l1, '"/>')
      WRITE (fileNum, 410) sliceplot%kk, sliceplot%e1s, sliceplot%e2s, sliceplot%nnne, input%pallst

!      <specialOutput form66="F" eonly="F" bmt="F"/>
420   FORMAT('      <specialOutput eonly="', l1, '" bmt="', l1, '"/>')
      WRITE (fileNum, 420) input%eonly, input%l_bmt

!      <magneticCircularDichroism energyLo="-10.0" energyUp="0.0"/>
430   FORMAT('      <magneticCircularDichroism energyLo="', f0.8, '" energyUp="', f0.8, '"/>')
      WRITE (fileNum, 430) banddos%e_mcd_lo, banddos%e_mcd_up

      WRITE (fileNum, '(a)') '   </output>'
      IF (present(filename)) THEN
         WRITE (fileNum, '(a)') '  <!-- We include the file relax.inp here to enable relaxations (see documentation) -->'
         WRITE (fileNum, '(a)') '  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="relax.xml"> <xi:fallback/> </xi:include>'
         WRITE (fileNum, '(a)') '</fleurInput>'
         CLOSE (fileNum)
      ELSE
         CALL closeXMLElement('inputData')
      END IF

   END SUBROUTINE w_inpXML
END MODULE m_winpXML
