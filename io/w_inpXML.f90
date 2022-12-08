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
      atoms, vacuum, input, stars, sliceplot, forcetheo, banddos, juPhon, &
      cell, sym, xcpot, noco,   mpinp, hybinp, kptsArray, kptsSelection, enpara, &
      gfinp, hub1inp, l_explicitIn, l_includeIn, filename)

      use m_types_input
      use m_types_sym
      use m_types_stars
      use m_types_atoms
      use m_types_vacuum
      use m_types_kpts
       
      use m_types_mpinp
      use m_types_hybinp
      use m_types_gfinp
      use m_types_hub1inp
      use m_types_cell
      use m_types_banddos
      use m_types_sliceplot
      USE m_types_xcpot
      USE m_types_xcpot_inbuild_nofunction
      USE m_types_noco
      use m_types_enparaxml
      USE m_types_forcetheo
      USE m_types_juPhon

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
      TYPE(t_kpts), INTENT(IN)     :: kptsArray(:)
       

      TYPE(t_mpinp), INTENT(IN)    :: mpinp
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_cell), INTENT(IN)     :: cell
      TYPE(t_banddos), INTENT(IN)  :: banddos
      TYPE(t_juPhon), INTENT(IN)   :: juPhon
      TYPE(t_sliceplot), INTENT(IN):: sliceplot
      CLASS(t_xcpot), INTENT(IN)   :: xcpot
      TYPE(t_noco), INTENT(IN)     :: noco
      TYPE(t_gfinp), INTENT(IN)    :: gfinp
      TYPE(t_hub1inp), INTENT(IN)  :: hub1inp
      CLASS(t_enparaxml), INTENT(IN)   :: enpara
      CLASS(t_forcetheo), INTENT(IN):: forcetheo !nothing is done here so far....
      CHARACTER(LEN=40)          :: kptsSelection(3) ! 1: default selection, 2: alternative for band structures, alternative for GW
      LOGICAL, INTENT(IN)        :: l_explicitIn, l_includeIn(4)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

      INTEGER          :: iSpecies, fileNum
      CHARACTER(len=8) :: name(10)
      INTEGER          :: numSpecies
      INTEGER          :: speciesRepAtomType(atoms%ntype)
      CHARACTER(len=20):: speciesNames(atoms%ntype)
      CHARACTER(LEN=50):: tempStringA, tempStringB
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
      INTEGER ieq, i, k, na, n, ilo,iContour, iKpts
      REAL s3, ah, a, hs2, rest
      LOGICAL l_hyb, ldum
      INTEGER :: ierr
! ..
!...  Local Arrays
      CHARACTER :: helpchar(atoms%ntype)
      CHARACTER(len=4) :: chntype
      CHARACTER(len=41) :: chform
      CHARACTER(len=100) :: line

!     added for HF and hybinp functionals
      REAL                  ::  aMix, omega
      INTEGER               :: idum
      CHARACTER(len=1)     ::  check

      CHARACTER(len=20) :: speciesName
      CHARACTER(len=150) :: format
      CHARACTER(len=20) :: mixingScheme
      CHARACTER(len=10) :: loType
      CHARACTER(len=10) :: bzIntMode
      LOGICAL ::   l_explicit, l_nocoOpt, l_gfOpt, l_include(4)
      INTEGER :: iAtomType, startCoreStates, endCoreStates
      CHARACTER(len=100) :: posString(3)
      CHARACTER(len=7) :: str
      REAL :: tempTaual(3, atoms%nat), scpos(3)
      REAL :: amatTemp(3, 3), bmatTemp(3, 3)

      l_include = l_includeIn .or. .not. present(filename)
      l_explicit = l_explicitIn .OR. .not. present(filename)
      l_nocoOpt = noco%l_noco .OR. juDFT_was_argument("-noco")
      l_gfOpt = gfinp%n>0 .OR. juDFT_was_argument("-greensf")

      band = .false.
      nw = 1

      IF (PRESENT(filename)) THEN
         filenum = 98
         OPEN (fileNum, file=TRIM(ADJUSTL(filename)), form='formatted', status='replace')
         WRITE (fileNum, '(a)') '<?xml version="1.0" encoding="UTF-8" standalone="no"?>'
         WRITE (fileNum, '(a)') '<fleurInput fleurInputVersion="0.36">'
      ELSE
         fileNum = getXMLOutputUnitNumber()
         CALL openXMLElementNoAttributes('inputData')
      END IF

      WRITE (fileNum, '(a)') '   <comment>'
      WRITE (fileNum, '(a6,10a)') '      ', input%comment
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

!      <coreElectrons ctail="T" frcor="F" kcrel="0" coretail_lmax="0" l_core_confpot="T"/>
130   FORMAT('      <coreElectrons ctail="', l1, '" frcor="', l1, '" kcrel="', i0, '" coretail_lmax="', i0, '"/>')
      WRITE (fileNum, 130) input%ctail, input%frcor, input%kcrel, input%coretail_lmax

      SELECT TYPE (xcpot)
      CLASS IS (t_xcpot_inbuild_nf)
         !   <xcFunctional name="pbe" relativisticCorrections="F">
135      FORMAT('      <xcFunctional name="', a, '" relativisticCorrections="', l1, '"/>')
         WRITE (fileNum, 135) trim(xcpot%get_name()), xcpot%relativistic_correction()
      END SELECT

!      <magnetism jspins="1" l_noco="F" l_J="F" swsp="F" lflip="F"/>
140   FORMAT('      <magnetism jspins="', i0, '" l_noco="', l1, '" l_ss="', l1, '">')
141   FORMAT('      <magnetism jspins="', i0, '"/>')
      IF(l_explicit.OR.l_nocoOpt) THEN
         WRITE (fileNum, 140) input%jspins, noco%l_noco, noco%l_ss
162      FORMAT('         <qss>', f0.10, ' ', f0.10, ' ', f0.10, '</qss>')
         WRITE (fileNum, 162) noco%qss_inp
164      FORMAT('         <mtNocoParams l_mperp="', l1, '" l_mtNocoPot="', l1,'" l_relaxSQA="', l1,'" mag_mixing_scheme="', i1, '" mix_RelaxWeightOffD="',f0.8,'" l_constrained="', l1,'" mix_constr="', f0.8,'"/>')
         WRITE (fileNum, 164) noco%l_mperp,any(noco%l_unrestrictMT), any(noco%l_alignMT), noco%mag_mixing_scheme, minval(noco%mix_RelaxWeightOffD), any(noco%l_constrained), noco%mix_b
166      FORMAT('         <sourceFreeMag l_sourceFree="', l1, '" l_scaleMag="', l1, '" mag_scale="', f0.8,'"/>')
         WRITE (fileNum, 166) noco%l_sourceFree, noco%l_scaleMag, noco%mag_scale
         WRITE (fileNum, '(a)') '      </magnetism>'
      ELSE
         WRITE (fileNum, 141) input%jspins
      END IF

      !      <soc theta="0.00000" phi="0.00000" l_soc="F" spav="F" off="F" soc66="F"/>
150   FORMAT('      <soc l_soc="', l1, '" theta="', f0.8, '" phi="', f0.8, '" spav="', l1, '"/>')
      WRITE (fileNum, 150) noco%l_soc, noco%theta_inp, noco%phi_inp, noco%l_spav

      IF (l_explicit .OR. hybinp%l_hybrid) THEN
155      FORMAT('      <prodBasis gcutm="', f0.8, '" tolerance="', f0.8, '" ewaldlambda="', i0, '" lexp="', i0, '" bands="', i0, '" fftcut="', f0.8, '"/>')
         WRITE (fileNum, 155) mpinp%g_cutoff, mpinp%linear_dep_tol, hybinp%ewaldlambda, hybinp%lexp, hybinp%bands1, hybinp%fftcut
      END IF


!      <expertModes spex="0"  eig66="F" lpr="0" secvar="F" />
180   FORMAT('      <expertModes spex="', i0, '" secvar="', l1, '"/>')
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

      SELECT CASE (input%bz_integration)
         CASE (BZINT_METHOD_HIST)
            bzIntMode = 'hist'
         CASE (BZINT_METHOD_GAUSS)
            bzIntMode = 'gauss'
         CASE (BZINT_METHOD_TRIA)
            bzIntMode = 'tria'
         CASE (BZINT_METHOD_TETRA)
            bzIntMode = 'tetra'
         CASE DEFAULT
            CALL judft_error("Invalid brillouin zone integration mode",calledby="w_inpXML")
      END SELECT

!      <ldaU l_linMix="F" mixParam="0.05" spinf="1.0" />
195   FORMAT('      <ldaU l_linMix="', l1, '" mixParam="', f0.6, '" spinf="', f0.6, '"/>')
      WRITE (fileNum, 195) input%ldauLinMix, input%ldauMixParam, input%ldauSpinf

      IF(atoms%n_hia>0 .OR. l_explicit) THEN
196      FORMAT('      <ldaHIA itmaxHubbard1="', i0, '" minoccDistance="', f0.6, '" minmatDistance="', f0.6, '" beta="', f0.1, '" dftspinpol="', l1, '"/>')
         WRITE (fileNum, 196) hub1inp%itmax, hub1inp%minoccDistance, hub1inp%minmatDistance, hub1inp%beta, hub1inp%l_dftspinpol
      ENDIF

      IF(l_gfOpt) THEN
205      FORMAT('      <greensFunction l_mperp="', l1'">')
         WRITE(fileNum, 205) gfinp%l_mperp
206      FORMAT('         <realAxis ne="', i0, '" ellow="', f0.8, '" elup="', f0.8, '"/>')
         WRITE(fileNum, 206) gfinp%ne, gfinp%ellow, gfinp%elup
         IF(gfinp%numberContours>0) THEN
            DO iContour = 1, gfinp%numberContours
               SELECT CASE(gfinp%contour(iContour)%shape)
               CASE(CONTOUR_RECTANGLE_CONST)
207               FORMAT('         <contourRectangle n1="', i0, '" n2="', i0, '" n3="', i0, '" nmatsub="', i0,&
                         '" sigma="', f0.8, '" eb="', f0.8, '" label="', a,'"/>')
                  WRITE(fileNum, 207) gfinp%contour(iContour)%n1, gfinp%contour(iContour)%n2, gfinp%contour(iContour)%n3,&
                                      gfinp%contour(iContour)%nmatsub, gfinp%contour(iContour)%sigma, gfinp%contour(iContour)%eb,&
                                      gfinp%contour(iContour)%label
               CASE(CONTOUR_SEMICIRCLE_CONST)
208               FORMAT('         <contourSemicircle n="', i0, '" eb="', f0.8, '" et="', f0.8, '" alpha="', f0.8, '" label="', a,'"/>')
                  WRITE(fileNum, 208) gfinp%contour(iContour)%ncirc, gfinp%contour(iContour)%eb, gfinp%contour(iContour)%et,&
                                      gfinp%contour(iContour)%alpha,gfinp%contour(iContour)%label
               CASE(CONTOUR_DOS_CONST)
209               FORMAT('         <contourDOS n="', i0, '" sigma="', f0.8, '" eb="', f0.8, '" et="', f0.8, &
                         '" analytical_cont="', l1, '" l_fermi="', l1, '" label="', a,'"/>')
                  WRITE(fileNum, 209) gfinp%contour(iContour)%nDOS, gfinp%contour(iContour)%sigmaDOS, gfinp%contour(iContour)%eb,&
                                      gfinp%contour(iContour)%et, gfinp%contour(iContour)%l_anacont, gfinp%contour(iContour)%l_dosfermi,&
                                      gfinp%contour(iContour)%label
               CASE DEFAULT
                  CALL judft_error("Unknown green's function contour mode", calledby="w_inpXML")
               END SELECT
            ENDDO
         ELSE
            !Write out a default contour (Semicircle)
            WRITE(fileNum, 208) 128, -1.0, 0.0,1.0,"default"
         ENDIF
         WRITE(fileNum, '(a)') '      </greensFunction>'
      ENDIF

! 

      WRITE (fileNum, '(a)') '   </calculationSetup>'
      WRITE (fileNum, '(a)') '   <cell>'

!      <bzIntegration valenceElectrons="8.00000" mode="hist" fermiSmearingEnergy="0.00100">
200   FORMAT('      <bzIntegration valenceElectrons="', f0.8, '" mode="', a, '" fermiSmearingEnergy="', f0.8, '">')
      WRITE (fileNum, 200) input%zelec, TRIM(ADJUSTL(bzIntMode)), input%tkb

210   FORMAT('         <kPointListSelection listName="', a, '"/>')
      WRITE (filenum, 210) TRIM(ADJUSTL(kptsSelection(1)))

!211   FORMAT('         <altKPointList listName="', a, '" purpose="', a, '"/>')
!      IF(kptsSelection(2).NE.'') THEN
!         WRITE (filenum, 211) TRIM(ADJUSTL(kptsSelection(2))), 'bands'
!      END IF
!      IF(kptsSelection(3).NE.'') THEN
!         WRITE (filenum, 211) TRIM(ADJUSTL(kptsSelection(3))), 'GW'
!      END IF

      if (l_include(1)) THEN
         WRITE (fileNum, '(a)') "         <kPointLists>"
         DO iKpts = 1, SIZE(kptsArray)
            CALL kptsArray(iKpts)%print_XML(fileNum)
         END DO
         WRITE (fileNum, '(a)') "         </kPointLists>"
      else
         WRITE (fileNum, '(a)') '         <!-- k-points included here -->'
         WRITE (fileNum, '(a)') '         <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="kpts.xml"> </xi:include>'
      end if
      WRITE (fileNum, '(a)') '      </bzIntegration>'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note: Different options for the cell definition!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (l_include(2)) THEN
         call sym%print_xml(fileNum)
      else
         WRITE (fileNum, '(a)') '      <!-- symmetry operations included here -->'
         WRITE (fileNum, '(a)') '      <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="sym.xml"> </xi:include>'
      end if
      IF (input%film) THEN
!      <xsd:attribute name="dVac" type="xsd:double" use="required"/>
!      <xsd:attribute name="dTilda" type="xsd:double" use="required"/>
!      <filmLattice ...>
241      FORMAT('      <filmLattice scale="', f0.8, '" dVac="', f0.8, '" dTilda="', f0.8, '">')
         WRITE (fileNum, 241) 1.0, vacuum%dvac, cell%amat(3, 3)
         !       <bravaisMatrixFilm>
         WRITE (fileNum, '(a)') '         <bravaisMatrixFilm>'
      !            <row-1>0.00000 5.13000 </row-1>
251      FORMAT('            <row-1>', f22.16, ' ', f22.16, '</row-1>')
         WRITE (fileNum, 251) cell%amat(:2, 1)
      !            <row-2>5.13000 0.00000 5.13000</row-2>
252      FORMAT('            <row-2>', f22.16, ' ', f22.16, '</row-2>')
         WRITE (fileNum, 252) cell%amat(:2, 2)
         WRITE (fileNum, '(a)') '         </bravaisMatrixFilm>'
      ELSE
242      FORMAT('      <bulkLattice scale="', f0.10, '">')
WRITE (fileNum, 242) 1.0
!         <bravaisMatrix>
   WRITE (fileNum, '(a)') '         <bravaisMatrix>'
!            <row-1>0.00000 5.13000 5.13000</row-1>
250      FORMAT('            <row-1>', f22.16, ' ', f22.16, ' ', f22.16, '</row-1>')
   WRITE (fileNum, 250) cell%amat(:, 1)
!            <row-2>5.13000 0.00000 5.13000</row-2>
260      FORMAT('            <row-2>', f22.16, ' ', f22.16, ' ', f22.16, '</row-2>')
   WRITE (fileNum, 260) cell%amat(:, 2)
!            <row-3>5.13000 5.13000 0.00000</row-3>
270      FORMAT('            <row-3>', f22.16, ' ', f22.16, ' ', f22.16, '</row-3>')
   WRITE (fileNum, 270) cell%amat(:, 3)
   WRITE (fileNum, '(a)') '         </bravaisMatrix>'
   ENDIF

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

      uIndices = -1
      DO i_u = 1, atoms%n_u
         IF (uIndices(1, atoms%lda_u(i_u)%atomType) .EQ. -1) uIndices(1, atoms%lda_u(i_u)%atomType) = i_u
         uIndices(2, atoms%lda_u(i_u)%atomType) = i_u
      END DO

      !Build list of species
      speciesNames = ''
      numSpecies = 0
      DO n = 1, atoms%ntype
         known_species = ANY(trim(atoms%speciesname(n)) == speciesNames(:numSpecies))
         if (.not. known_species) THEN
            numSpecies = numSpecies + 1
            speciesNames(numSpecies) = trim(atoms%speciesname(n))
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
300      FORMAT('      <species name="', a, '" element="', a, '" atomicNumber="', i0, '">')
         speciesName = TRIM(ADJUSTL(speciesNames(iSpecies)))
         WRITE (fileNum, 300) TRIM(ADJUSTL(speciesName)), TRIM(ADJUSTL(namat_const(atoms%nz(iAtomType)))), atoms%nz(iAtomType)

!         <mtSphere radius="2.160000" gridPoints="521" logIncrement="0.022000"/>
310      FORMAT('         <mtSphere radius="', f0.8, '" gridPoints="', i0, '" logIncrement="', f0.8, '"/>')
         WRITE (fileNum, 310) atoms%rmt(iAtomType), atoms%jri(iAtomType), atoms%dx(iAtomType)

!         <atomicCutoffs lmax="8" lnonsphr="6"/>
320      FORMAT('         <atomicCutoffs lmax="', i0, '" lnonsphr="', i0, '"/>')
         WRITE (fileNum, 320) atoms%lmax(iAtomType), atoms%lnonsph(iAtomType)

         WRITE (fileNum, '(a)') '         <electronConfig flipSpins="F">'
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
         IF (l_explicit .OR. hybinp%l_hybrid) THEN
315         FORMAT('         <prodBasis lcutm="', i0, '" lcutwf="', i0, '" select="', a, '"/>')
            line = ''
            WRITE (line, '(i0,1x,i0,1x,i0,1x,i0)') hybinp%select1(1:4, iAtomType)
            WRITE (fileNum, 315) hybinp%lcutm1(iAtomType), hybinp%lcutwf(iAtomType), TRIM(ADJUSTL(line))
         END IF

         IF (l_explicit) THEN
328         FORMAT('         <modInitDen flipSpinPhi="', f0.8, '" flipSpinTheta="', f0.8, '" flipSpinScale="', l1, '"/>')
            WRITE (fileNum, 328) atoms%flipSpinPhi(iAtomType), atoms%flipSpinTheta(iAtomType), atoms%flipSpinScale(iAtomType)
         END IF

         IF (uIndices(1, iAtomType) .NE. -1) THEN
!         <ldaU l="2" U="5.5" J="0.9" l_amf="F"/>
            DO i_u = uIndices(1, iAtomType), uIndices(2, iAtomType)
326            FORMAT('         <ldaU l="', i0, '" U="', f0.5, '" J="', f0.5, '" l_amf="', l1, '"/>')
               WRITE (fileNum, 326) atoms%lda_u(i_u)%l, atoms%lda_u(i_u)%u, atoms%lda_u(i_u)%j, atoms%lda_u(i_u)%l_amf
            END DO
         END IF

         IF(l_gfOpt) THEN
            WRITE (fileNum,316) .TRUE., "default",0,"calc"
            WRITE (fileNum,318) .FALSE.,.FALSE.,.FALSE.,.FALSE.
            WRITE (fileNum, '(a)') '         </greensfCalculation>'
316         FORMAT('         <greensfCalculation l_sphavg="', l1, '" label="', a, '" nshells="', i0, '" kkintgrCutoff="', a, '">')
318         FORMAT('            <diagElements s="', l1, '" p="', l1, '" d="', l1, '" f="', l1, '"/>')
         ENDIF

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
            IF (input%film) THEN
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
362         FORMAT('         <nocoParams  alpha="', f0.8, '" beta="', &
                   f0.8,  '"/>')
            WRITE (fileNum, 362)  noco%alph_inp(iAtomType), &
               noco%beta_inp(iAtomType)
         END IF

         WRITE (fileNum, '(a)') '      </atomGroup>'
      END DO
      WRITE (fileNum, '(a)') '   </atomGroups>'
      if (.not. l_include(4)) then
         close (99)
         fileNum = 98
      endif

368   FORMAT('   <output dos="', l1, '" band="', l1,  '" slice="', l1, '">')
      WRITE (fileNum, 368) banddos%dos, band, sliceplot%slice

!      <checks vchk="F" cdinf="F" disp="F"/>
370   FORMAT('      <checks vchk="', l1, '" cdinf="', l1, '"/>')
      WRITE (fileNum, 370) input%vchk, input%cdinf

!      <densityOfStates ndir="0" minEnergy="-0.50000" maxEnergy="0.50000" sigma="0.01500"/>
      WRITE(tempStringA,'(f0.8,a)') banddos%e2_dos, '*Htr'
      WRITE(tempStringB,'(f0.8,a)') banddos%e1_dos, '*Htr'
380   FORMAT('      <bandDOS minEnergy="', a, '" maxEnergy="', a, '" sigma="', f0.8, '" storeEVData="', l1, '"/>')
      WRITE (fileNum, 380)  TRIM(ADJUSTL(tempStringA)), TRIM(ADJUSTL(tempStringB)), banddos%sig_dos, banddos%l_storeEVData

!      <vacuumDOS layers="0" integ="F" star="F" nstars="0" locx1="0.00" locy1="0.00" locx2="0.00" locy2="0.00" nstm="0" tworkf="0.000000"/>
390   FORMAT('      <vacuumDOS vacdos="', l1, '" integ="', l1, '" star="', l1, '" nstars="', i0, '" locx1="', f0.5, '" locy1="', f0.5, '" locx2="', f0.5, '" locy2="', f0.5, '" nstm="', i0, '" tworkf="', f0.5, '"/>')
      WRITE (fileNum, 390) banddos%vacdos, input%integ, banddos%starcoeff, banddos%nstars, banddos%locx(1), banddos%locy(1), banddos%locx(2), banddos%locy(2), 0, 0.0

!      <unfoldingBand unfoldBand="F" supercellX="1" supercellY="1" supercellZ="1"/>
395   FORMAT('      <unfoldingBand unfoldBand="', l1, '" supercellX="', i0, '" supercellY="', i0, '" supercellZ="', i0, '"/>')
      WRITE (fileNum, 395) banddos%unfoldband, banddos%s_cell_x, banddos%s_cell_y, banddos%s_cell_z

!!      <juPhon l_potout="F" l_eigout="F"/>
!396   FORMAT('      <juPhon l_potout="', l1, '" l_eigout="', l1, '"/>')
!      WRITE (fileNum, 396) juPhon%l_potout, juPhon%l_eigout

!      <plotting iplot="0" />
      IF(SIZE(sliceplot%plot)>0) THEN
400      FORMAT('      <plotting iplot="', i0, '" polar="', l1, '">')
         WRITE (fileNum, 400) sliceplot%iplot, sliceplot%polar
401      FORMAT('         <plot TwoD="', l1, '" vec1="', 3f5.1,  '" vec2="', 3f5.1, '" vec3="', 3f5.1, '" zero="', 3f5.1, '" file="', a, '"/>')
         WRITE (fileNum, 401) sliceplot%plot(1)%twodim, sliceplot%plot(1)%vec1(:), sliceplot%plot(1)%vec2(:), sliceplot%plot(1)%vec3(:),&
                              sliceplot%plot(1)%zero(:), TRIM(ADJUSTL(sliceplot%plot(1)%filename))
         WRITE (fileNum, '(a)') '      </plotting>'
      ELSE
402      FORMAT('      <plotting iplot="', i0, '" polar="', l1, '"/>')
         WRITE (fileNum, 402) sliceplot%iplot, sliceplot%polar
      ENDIF
      
!      <chargeDensitySlicing numkpt="0" minEigenval="0.000000" maxEigenval="0.000000" nnne="0" pallst="F"/>
410   FORMAT('      <chargeDensitySlicing numkpt="', i0, '" minEigenval="', f0.8, '" maxEigenval="', f0.8, '" nnne="', i0, '" pallst="', l1, '"/>')
      WRITE (fileNum, 410) sliceplot%kk, sliceplot%e1s, sliceplot%e2s, sliceplot%nnne, input%pallst

!      <specialOutput form66="F" eonly="F" bmt="F"/>
420   FORMAT('      <specialOutput eonly="', l1, '" bmt="', l1, '"/>')
      WRITE (fileNum, 420) input%eonly, input%l_bmt

!      <magneticCircularDichroism energyLo="-10.0" energyUp="0.0"/>
430   FORMAT('      <magneticCircularDichroism mcd="',l1,'" energyLo="', f0.8, '" energyUp="', f0.8, '"/>')
      WRITE (fileNum, 430) banddos%l_mcd,banddos%e_mcd_lo, banddos%e_mcd_up

      WRITE (fileNum, '(a)') '   </output>'
      IF (present(filename)) THEN
         WRITE (fileNum, '(a)') '  <!-- We include the file relax.xml here to enable relaxations (see documentation) -->'
         WRITE (fileNum, '(a)') '  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="relax.xml"> <xi:fallback/> </xi:include>'
         WRITE (fileNum, '(a)') '</fleurInput>'
         CLOSE (fileNum)
      ELSE
         CALL closeXMLElement('inputData')
      END IF

   END SUBROUTINE w_inpXML
END MODULE m_winpXML
