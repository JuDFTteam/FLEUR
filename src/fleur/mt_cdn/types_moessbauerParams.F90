!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_moessbauerParams
   implicit none

   PRIVATE

   PUBLIC :: t_moessbauerParams

   TYPE t_moessbauerParams

      ! Per-category on/off switches, from the optional <moessbauerParams>
      ! tag in the <output> section of inp.xml (all on by default). Since
      ! the orbital contribution to the valence hyperfine field is only
      ! meaningful with spin-orbit coupling, l_hyperfine additionally
      ! requires noco%l_soc -- with SOC off, enabling valenceHyperfine in
      ! inp.xml is equivalent to leaving it off. It is also restricted to
      ! kcrel=1, so that valence and core hyperfine fields are only ever
      ! reported together, never with just one of the two available.
      LOGICAL :: l_hyperfine
      LOGICAL :: l_calcEFG
      LOGICAL :: l_calcIsomerShift
      LOGICAL :: l_calcCoreHyperfine

      ! Valence hyperfine field: contact/dipolar/orbital contributions.
      REAL, ALLOCATABLE :: hypFineContribs(:,:,:,:) ! (-1:3,ntype,jspins,3)

      ! Electric field gradient tensor.
      REAL, ALLOCATABLE    :: efgTensor(:,:,:) ! (3,3,ntype)
      LOGICAL, ALLOCATABLE :: l_efgValid(:)

      ! Valence+core isomer shift (average charge density over the nucleus).
      REAL, ALLOCATABLE    :: isomerShiftNucRad(:)
      REAL, ALLOCATABLE    :: isomerShiftAverageDen(:)
      LOGICAL, ALLOCATABLE :: l_isomerShiftValid(:)

      ! Relativistic core hyperfine field + core isomer shift (per core shell).
      INTEGER, ALLOCATABLE :: nCoreShells(:)
      INTEGER, ALLOCATABLE :: coreShellN(:,:), coreShellL(:,:) ! (15,ntype)
      REAL,    ALLOCATABLE :: coreHFF(:,:), coreIsomerShift(:,:) ! (15,ntype)
      LOGICAL              :: l_coreHFFAvailable
      LOGICAL              :: l_coreIsomerShiftAvailable

      CONTAINS

      PROCEDURE, PASS :: init => mPInit
      PROCEDURE, PASS :: calcEFG => mPCalcEFG
      PROCEDURE, PASS :: calcIS => mPCalcIS
      PROCEDURE, PASS :: setCoreHFF => mPSetCoreHFF
      PROCEDURE, PASS :: accumCoreIS => mPAccumCoreIS
      PROCEDURE, PASS :: printAll => mPPrintAll

   END TYPE t_moessbauerParams

   CONTAINS

   SUBROUTINE mPInit(this, input, noco, atoms)

      USE m_types

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_noco),              INTENT(IN)    :: noco
      TYPE(t_atoms),             INTENT(IN)    :: atoms

      this%l_calcEFG = input%l_moessbauerEFG
      this%l_calcIsomerShift = input%l_moessbauerIsomerShift
      this%l_calcCoreHyperfine = input%l_moessbauerCoreHyperfine
      ! Valence hyperfine field needs jspins==2 (spin density) and, for a
      ! meaningful orbital term, noco%l_soc; restricted to kcrel=1 so it is
      ! only ever calculated/reported alongside the core hyperfine field.
      this%l_hyperfine = input%l_moessbauerValenceHyperfine .AND. noco%l_soc &
                         .AND. (input%jspins.EQ.2) .AND. (input%kcrel.EQ.1)

      IF (.NOT.ALLOCATED(this%hypFineContribs)) THEN
         ALLOCATE(this%hypFineContribs(-1:3,atoms%ntype,input%jspins,3))
         ALLOCATE(this%efgTensor(3,3,atoms%ntype))
         ALLOCATE(this%l_efgValid(atoms%ntype))
         ALLOCATE(this%isomerShiftNucRad(atoms%ntype))
         ALLOCATE(this%isomerShiftAverageDen(atoms%ntype))
         ALLOCATE(this%l_isomerShiftValid(atoms%ntype))
         ALLOCATE(this%nCoreShells(atoms%ntype))
         ALLOCATE(this%coreShellN(15,atoms%ntype))
         ALLOCATE(this%coreShellL(15,atoms%ntype))
         ALLOCATE(this%coreHFF(15,atoms%ntype))
         ALLOCATE(this%coreIsomerShift(15,atoms%ntype))
      END IF

      this%hypFineContribs = 0.0
      this%efgTensor = 0.0
      this%l_efgValid = .FALSE.
      this%isomerShiftNucRad = 0.0
      this%isomerShiftAverageDen = 0.0
      this%l_isomerShiftValid = .FALSE.
      this%nCoreShells = 0
      this%coreShellN = 0
      this%coreShellL = 0
      this%coreHFF = 0.0
      this%coreIsomerShift = 0.0
      this%l_coreHFFAvailable = .FALSE.
      this%l_coreIsomerShiftAvailable = .FALSE.

   END SUBROUTINE mPInit

   SUBROUTINE mPCalcEFG(this, atoms, sym, sphhar, fmpi, vCoul)

      USE m_types
      USE m_efg

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_sym),               INTENT(IN)    :: sym
      TYPE(t_sphhar),            INTENT(IN)    :: sphhar
      TYPE(t_mpi),               INTENT(IN)    :: fmpi
      TYPE(t_potden),            INTENT(IN)    :: vCoul

      IF (.NOT.this%l_calcEFG) RETURN

      CALL calc_efg(atoms, sym, sphhar, fmpi, vCoul, this%efgTensor, this%l_efgValid)

   END SUBROUTINE mPCalcEFG

   SUBROUTINE mPCalcIS(this, input, atoms, fmpi, den)

      USE m_constants
      USE m_types
      USE m_juDFT
      USE m_intgr, ONLY : intgr2

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      TYPE(t_input),             INTENT(IN)    :: input
      TYPE(t_atoms),             INTENT(IN)    :: atoms
      TYPE(t_mpi),               INTENT(IN)    :: fmpi
      TYPE(t_potden),            INTENT(IN)    :: den

      INTEGER, PARAMETER :: nFitPoints = 4

      INTEGER :: iType, i, iRad, iStart, iEnd
      REAL    :: nucRad
      REAL    :: indefInteg(atoms%jmtd), sphrDen(atoms%jmtd)
      REAL    :: avgDenAtMesh(nFitPoints)
      LOGICAL :: isSmaller
      CHARACTER(LEN=100) :: warnMsg

      IF (fmpi%irank /= 0) RETURN
      IF (.NOT.this%l_calcIsomerShift) RETURN

      DO iType = 1, atoms%ntype
         nucRad = r0_const*(atomicMasses_const(atoms%nz(iType))**(1.0/3.0))
         this%isomerShiftNucRad(iType) = nucRad

         IF (atoms%rmsh(1,iType) >= nucRad) THEN
            this%l_isomerShiftValid(iType) = .FALSE.
            WRITE(warnMsg,'(a,i0,a)') 'Nuclear radius smaller than smallest radial mesh point for atom type ', &
                                       iType, ' -- skipping its isomer shift'
            CALL juDFT_warn(TRIM(warnMsg), calledby='types_moessbauerParams')
            CYCLE
         END IF

         ! Locate the mesh window around r=nucRad: since the mesh is
         ! logarithmic, nucRad is typically many mesh points away from the
         ! innermost one (r(1)), so the fit below must use a local window
         ! straddling nucRad, not simply the mesh's first nFitPoints points.
         iRad = 0
         isSmaller = .TRUE.
         DO WHILE (isSmaller)
            iRad = iRad + 1
            IF (atoms%rmsh(iRad,iType).GE.nucRad) isSmaller = .FALSE.
         END DO
         iStart = MAX(1, iRad-2)
         iEnd = iStart + nFitPoints - 1
         IF (iEnd > atoms%jri(iType)) THEN
            iEnd = atoms%jri(iType)
            iStart = iEnd - nFitPoints + 1
         END IF

         IF (input%jspins.EQ.1) THEN
            sphrDen(:) = den%mt(:,0,iType,1)
         ELSE
            sphrDen(:) = den%mt(:,0,iType,1) + den%mt(:,0,iType,2)
         END IF

         ! Average density enclosed within the nuclear sphere: quadratic
         ! least-squares fit (mirrors extrapolate_r0 in efg.F90) of the
         ! enclosed-charge/volume ratio against r, using the nFitPoints mesh
         ! points straddling nucRad, evaluated at r=nucRad. A multi-point fit
         ! (as opposed to a 2-point bracket-and-interpolate) is less
         ! sensitive to noise at any single mesh point.
         indefInteg = 0.0
         CALL intgr2(sphrDen(:),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),indefInteg)

         DO i = 1, nFitPoints
            avgDenAtMesh(i) = indefInteg(iStart+i-1)*sfp_const / ((4.0/3.0)*pi_const*(atoms%rmsh(iStart+i-1,iType)**3.0))
         END DO
         this%isomerShiftAverageDen(iType) = quadFitAtR(avgDenAtMesh, atoms%rmsh(iStart:iEnd,iType), nucRad)

         this%l_isomerShiftValid(iType) = .TRUE.
      END DO

   END SUBROUTINE mPCalcIS

   SUBROUTINE mPSetCoreHFF(this, iType, nshell, nqntab, lqntab, bhff, isomerShift)

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      INTEGER,                   INTENT(IN)    :: iType, nshell
      INTEGER,                   INTENT(IN)    :: nqntab(15), lqntab(15)
      REAL,                      INTENT(IN)    :: bhff(15), isomerShift(15)

      IF (.NOT.this%l_calcCoreHyperfine) RETURN

      this%nCoreShells(iType) = nshell
      this%coreShellN(1:nshell,iType) = nqntab(1:nshell)
      this%coreShellL(1:nshell,iType) = lqntab(1:nshell)
      this%coreHFF(1:nshell,iType) = bhff(1:nshell)
      this%coreIsomerShift(1:nshell,iType) = isomerShift(1:nshell)
      this%l_coreHFFAvailable = .TRUE.
      this%l_coreIsomerShiftAvailable = .TRUE.

   END SUBROUTINE mPSetCoreHFF

   SUBROUTINE mPAccumCoreIS(this, iType, nshell, nqntab, lqntab, isomerShift, l_valid)

      !-------------------------------------------------------------------------
      ! Accumulates one spin channel's contribution to the core isomer shift
      ! (per shell) for the kcrel=0 core solver (src/fleur/core/cored.F90),
      ! which calls this once per (atom type, spin) combination.
      !-------------------------------------------------------------------------

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      INTEGER,                   INTENT(IN)    :: iType, nshell
      INTEGER,                   INTENT(IN)    :: nqntab(15), lqntab(15)
      REAL,                      INTENT(IN)    :: isomerShift(15)
      LOGICAL,                   INTENT(IN)    :: l_valid

      IF (.NOT.this%l_calcCoreHyperfine) RETURN

      IF ((.NOT.l_valid) .OR. (nshell == 0)) THEN
         this%nCoreShells(iType) = 0
         RETURN
      END IF

      this%nCoreShells(iType) = nshell
      this%coreShellN(1:nshell,iType) = nqntab(1:nshell)
      this%coreShellL(1:nshell,iType) = lqntab(1:nshell)
      this%coreIsomerShift(1:nshell,iType) = this%coreIsomerShift(1:nshell,iType) + isomerShift(1:nshell)
      this%l_coreIsomerShiftAvailable = .TRUE.

   END SUBROUTINE mPAccumCoreIS

   SUBROUTINE mPPrintAll(this, fmpi, atoms)

      USE m_constants
      USE m_types
      USE m_xmlOutput

      CLASS(t_moessbauerParams), INTENT(INOUT) :: this
      TYPE(t_mpi),               INTENT(IN)    :: fmpi
      TYPE(t_atoms),             INTENT(IN)    :: atoms

      CHARACTER, PARAMETER :: txtl(0:3) = (/'s','p','d','f'/)

      INTEGER :: iType, iShell, l
      REAL    :: a0, e0, cautog, bohrMagInCGS
      REAL    :: hyperfineResults(-1:3), hyperfineResultsTotal(-1:3)
      REAL    :: eig(3), eta, btot, istot
      REAL    :: Vij(3,3)
      REAL    :: work(9)
      INTEGER :: info
      EXTERNAL :: dsyev

      ! Stashed away here as they are computed below (in the same units/scale
      ! as the values written to the "out" file) so the XML section at the
      ! end of this routine can report exactly what was printed there,
      ! without redoing any of the unit conversions a second time.
      REAL    :: contactKG(-1:3,atoms%ntype), dipolarKG(-1:3,atoms%ntype)
      REAL    :: orbitalKG(-1:3,atoms%ntype), allTermsKG(-1:3,atoms%ntype)
      REAL    :: eigStash(3,atoms%ntype), etaStash(atoms%ntype)

      CHARACTER(LEN=24) :: efgMatChar(3,3)
      INTEGER           :: efgLengths(5,2)
      INTEGER           :: i, j

      IF (fmpi%irank /= 0) RETURN

      contactKG = 0.0; dipolarKG = 0.0; orbitalKG = 0.0; allTermsKG = 0.0
      eigStash = 0.0; etaStash = 0.0

      IF (this%l_hyperfine) THEN
         a0 = bohr_to_angstrom_const * 1.0e-8
         e0 = 1.6021892e-19 * 2.997930e+09
         cautog = e0 / (a0*a0)
         bohrMagInCGS = 1.0/(2.0*c_light(1.0))
         WRITE(oUnit,*) ''
         WRITE(oUnit,*) ' Hyperfine field valence contributions in kG '
         WRITE(oUnit,*) ' ========================================================== '
         WRITE(oUnit,*) ' atom type                          contribution'
         WRITE(oUnit,*) '                    total            s            p            d            f'
         DO iType = 1, atoms%ntype
            this%hypFineContribs(:,iType,1,1) = this%hypFineContribs(:,iType,1,1) - this%hypFineContribs(:,iType,2,1)
            hyperfineResults(:) = this%hypFineContribs(:,iType,1,1) * cautog * 0.001 * sfp_const * bohrMagInCGS * 8.0 * pi_const / 3.0
            WRITE(oUnit,'(i7,5x,5f13.5,5x,a)') iType, hyperfineResults(-1:3), 'contact term'
            hyperfineResultsTotal(:) = hyperfineResults(:)
            contactKG(:,iType) = hyperfineResults(:)
            this%hypFineContribs(:,iType,1,2) = this%hypFineContribs(:,iType,1,2) - this%hypFineContribs(:,iType,2,2)
            hyperfineResults(:) = this%hypFineContribs(:,iType,1,2) * cautog * 0.001 * bohrMagInCGS
            WRITE(oUnit,'(i7,5x,5f13.5,5x,a)') iType, hyperfineResults(-1:3), 'dipolar term'
            hyperfineResultsTotal(:) = hyperfineResultsTotal(:) + hyperfineResults(:)
            dipolarKG(:,iType) = hyperfineResults(:)
            this%hypFineContribs(:,iType,1,3) = this%hypFineContribs(:,iType,1,3) + this%hypFineContribs(:,iType,2,3)
            hyperfineResults(:) = this%hypFineContribs(:,iType,1,3) * cautog * 0.001 / c_light(1.0)
            WRITE(oUnit,'(i7,5x,5f13.5,5x,a)') iType, hyperfineResults(-1:3), 'orbital term'
            hyperfineResultsTotal(:) = hyperfineResultsTotal(:) + hyperfineResults(:)
            orbitalKG(:,iType) = hyperfineResults(:)
            WRITE(oUnit,'(i7,5x,5f13.5,5x,a)') iType, hyperfineResultsTotal(-1:3), 'all terms'
            allTermsKG(:,iType) = hyperfineResultsTotal(:)
         END DO
         WRITE(oUnit,*) ' ========================================================== '
      END IF

      IF (this%l_coreHFFAvailable) THEN
         WRITE(oUnit,*) ''
         WRITE(oUnit,*) ' Core hyperfine field / isomer shift per shell '
         WRITE(oUnit,*) ' ========================================================== '
         DO iType = 1, atoms%ntype
            IF (this%nCoreShells(iType) == 0) CYCLE
            WRITE(oUnit,'(a,i0)') ' atom type ', iType
            btot = 0.0
            istot = 0.0
            DO iShell = 1, this%nCoreShells(iType)
               l = this%coreShellL(iShell,iType)
               WRITE(oUnit,'(I4,A1,23X,F13.5,3X,F13.5)') this%coreShellN(iShell,iType), txtl(l), &
                  this%coreHFF(iShell,iType), this%coreIsomerShift(iShell,iType)
               btot = btot + this%coreHFF(iShell,iType)
               istot = istot + this%coreIsomerShift(iShell,iType)
            END DO
            WRITE(oUnit,'('' HFFTOT, ISTOT:'',13X,F13.5,3X,F13.5)') btot, istot
         END DO
         WRITE(oUnit,*) ' ========================================================== '
      ELSE IF (this%l_coreIsomerShiftAvailable) THEN
         ! kcrel=0: no coupled-branch Dirac solver, so no physically
         ! meaningful core hyperfine field is available (see
         ! mPCalcCoreISNonRel) -- only the isomer shift per shell.
         WRITE(oUnit,*) ''
         WRITE(oUnit,*) ' Core isomer shift per shell '
         WRITE(oUnit,*) ' ========================================== '
         DO iType = 1, atoms%ntype
            IF (this%nCoreShells(iType) == 0) CYCLE
            WRITE(oUnit,'(a,i0)') ' atom type ', iType
            istot = 0.0
            DO iShell = 1, this%nCoreShells(iType)
               l = this%coreShellL(iShell,iType)
               WRITE(oUnit,'(I4,A1,23X,F13.5)') this%coreShellN(iShell,iType), txtl(l), &
                  this%coreIsomerShift(iShell,iType)
               istot = istot + this%coreIsomerShift(iShell,iType)
            END DO
            WRITE(oUnit,'('' ISTOT:'',27X,F13.5)') istot
         END DO
         WRITE(oUnit,*) ' ========================================== '
      END IF

      IF (this%l_calcIsomerShift) THEN
         WRITE(oUnit,*) ''
         WRITE(oUnit,*) ' Isomer shift output (average charge density over nucleus) '
         WRITE(oUnit,*) ' ============================================================== '
         WRITE(oUnit,*) ' atom type     radius of nucleus (Bohr)   smallest radial mesh point (Bohr)   average charge density (e/Bohr^3) '
         DO iType = 1, atoms%ntype
            IF (.NOT.this%l_isomerShiftValid(iType)) CYCLE
            WRITE(oUnit,'(i7,10x,f15.7,15x,f15.7,20x,f13.5)') iType, this%isomerShiftNucRad(iType), &
               atoms%rmsh(1,iType), this%isomerShiftAverageDen(iType)
         END DO
         WRITE(oUnit,*) ' ============================================================== '
      END IF

      IF (this%l_calcEFG) THEN
         WRITE(oUnit,*) ''
         WRITE(oUnit,*) ' Electric field gradient (EFG) tensor at the nuclei (Hartree/bohr^2) '
         WRITE(oUnit,*) ' ================================================================== '
         DO iType = 1, atoms%ntype
            IF (.NOT.this%l_efgValid(iType)) CYCLE

            Vij = this%efgTensor(:,:,iType)
            CALL dsyev('N','U',3,Vij,3,eig,work,9,info)
            CALL sort_by_abs_descending(eig)

            eta = 0.0
            IF (ABS(eig(3)) > 1.0e-12) eta = (eig(1)-eig(2)) / eig(3)
            eigStash(:,iType) = eig
            etaStash(iType) = eta

            WRITE(oUnit,'(a,i0)') ' atom type ', iType
            WRITE(oUnit,'(a,3f13.5)') '   V_xx, V_yy, V_zz : ', this%efgTensor(1,1,iType), this%efgTensor(2,2,iType), this%efgTensor(3,3,iType)
            WRITE(oUnit,'(a,3f13.5)') '   V_xy, V_xz, V_yz : ', this%efgTensor(1,2,iType), this%efgTensor(1,3,iType), this%efgTensor(2,3,iType)
            WRITE(oUnit,'(a,3f13.5)') '   principal values V_XX, V_YY, V_ZZ : ', eig(1), eig(2), eig(3)
            WRITE(oUnit,'(a,f13.5)')  '   asymmetry parameter               : ', eta
         END DO
         WRITE(oUnit,*) ' ================================================================== '
      END IF

      CALL openXMLElement('moessbauerParams',(/'unitsHFF','unitsIS ','unitsEFG'/),&
                           (/'kG        ','e/bohr^3  ','Htr/bohr^2'/))

      IF (this%l_hyperfine) THEN
         DO iType = 1, atoms%ntype
            CALL openXMLElement('valenceHyperfineField',(/'atomType'/),&
                                 (/i2s(iType)/))
            CALL writeXMLElement('contact',(/'total','s    ','p    ','d    ','f    '/),&
                 (/r2s(contactKG(-1,iType)),r2s(contactKG(0,iType)),r2s(contactKG(1,iType)),&
                   r2s(contactKG(2,iType)),r2s(contactKG(3,iType))/))
            CALL writeXMLElement('dipolar',(/'total','s    ','p    ','d    ','f    '/),&
                 (/r2s(dipolarKG(-1,iType)),r2s(dipolarKG(0,iType)),r2s(dipolarKG(1,iType)),&
                   r2s(dipolarKG(2,iType)),r2s(dipolarKG(3,iType))/))
            CALL writeXMLElement('orbital',(/'total','s    ','p    ','d    ','f    '/),&
                 (/r2s(orbitalKG(-1,iType)),r2s(orbitalKG(0,iType)),r2s(orbitalKG(1,iType)),&
                   r2s(orbitalKG(2,iType)),r2s(orbitalKG(3,iType))/))
            CALL writeXMLElement('all',(/'total','s    ','p    ','d    ','f    '/),&
                 (/r2s(allTermsKG(-1,iType)),r2s(allTermsKG(0,iType)),r2s(allTermsKG(1,iType)),&
                   r2s(allTermsKG(2,iType)),r2s(allTermsKG(3,iType))/))
            CALL closeXMLElement('valenceHyperfineField')
         END DO
      END IF

      IF (this%l_coreHFFAvailable) THEN
         DO iType = 1, atoms%ntype
            IF (this%nCoreShells(iType) == 0) CYCLE
            btot = SUM(this%coreHFF(1:this%nCoreShells(iType),iType))
            istot = SUM(this%coreIsomerShift(1:this%nCoreShells(iType),iType))
            CALL openXMLElement('coreHyperfineField',(/'atomType','hffTotal','isTotal '/),&
                                 (/i2s(iType),r2s(btot),r2s(istot)/))
            DO iShell = 1, this%nCoreShells(iType)
               CALL writeXMLElement('shell',(/'n          ','l          ','hff        ','isomerShift'/),&
                    (/i2s(this%coreShellN(iShell,iType)),i2s(this%coreShellL(iShell,iType)),&
                      r2s(this%coreHFF(iShell,iType)),r2s(this%coreIsomerShift(iShell,iType))/))
            END DO
            CALL closeXMLElement('coreHyperfineField')
         END DO
      ELSE IF (this%l_coreIsomerShiftAvailable) THEN
         DO iType = 1, atoms%ntype
            IF (this%nCoreShells(iType) == 0) CYCLE
            istot = SUM(this%coreIsomerShift(1:this%nCoreShells(iType),iType))
            CALL openXMLElement('coreIsomerShift',(/'atomType','isTotal '/),&
                                 (/i2s(iType),r2s(istot)/))
            DO iShell = 1, this%nCoreShells(iType)
               CALL writeXMLElement('shell',(/'n          ','l          ','isomerShift'/),&
                    (/i2s(this%coreShellN(iShell,iType)),i2s(this%coreShellL(iShell,iType)),&
                      r2s(this%coreIsomerShift(iShell,iType))/))
            END DO
            CALL closeXMLElement('coreIsomerShift')
         END DO
      END IF

      IF (this%l_calcIsomerShift) THEN
         DO iType = 1, atoms%ntype
            IF (.NOT.this%l_isomerShiftValid(iType)) CYCLE
            CALL writeXMLElement('isomerShift',&
                 (/'atomType      ','nucleusRadius ',&
                   'averageDensity'/),&
                 (/i2s(iType),r2sHP(this%isomerShiftNucRad(iType)),&
                   r2s(this%isomerShiftAverageDen(iType))/))
         END DO
      END IF

      IF (this%l_calcEFG) THEN
         DO iType = 1, atoms%ntype
            IF (.NOT.this%l_efgValid(iType)) CYCLE
            DO j = 1, 3
               DO i = 1, 3
                  efgMatChar(i,j) = r2s(this%efgTensor(i,j,iType))
               END DO
            END DO
            efgLengths = 0
            CALL writeXMLElementMatrixForm('electricFieldGradient',&
                 (/'atomType ','VXX      ',&
                   'VYY      ','VZZ      ','asymParam'/),&
                 (/i2s(iType),r2s(eigStash(1,iType)),&
                   r2s(eigStash(2,iType)),r2s(eigStash(3,iType)),r2s(etaStash(iType))/),&
                 efgLengths, efgMatChar)
         END DO
      END IF

      CALL closeXMLElement('moessbauerParams')

   END SUBROUTINE mPPrintAll

   SUBROUTINE sort_by_abs_descending(eig)
      ! Sorts 3 eigenvalues in place so that |eig(3)| >= |eig(2)| >= |eig(1)|.
      REAL, INTENT(INOUT) :: eig(3)
      REAL :: tmp
      INTEGER :: i, j

      DO i = 1, 2
         DO j = 1, 3-i
            IF (ABS(eig(j)) > ABS(eig(j+1))) THEN
               tmp = eig(j); eig(j) = eig(j+1); eig(j+1) = tmp
            END IF
         END DO
      END DO

   END SUBROUTINE sort_by_abs_descending

   REAL FUNCTION quadFitAtR(vals, rmsh, r0)
      ! Quadratic (y = a + b*x + c*x^2) least-squares fit of vals(i) against
      ! rmsh(i), evaluated at r=r0. Used to get a value at (or extrapolated
      ! to) a small radius r0 from a few points of a radial mesh, more
      ! robustly than a 2-point bracket-and-interpolate.
      REAL, INTENT(IN) :: vals(:), rmsh(:)
      REAL, INTENT(IN) :: r0

      REAL :: mat(3,3), rhs(3), coef(3)

      mat(1,:) = (/ REAL(SIZE(vals)), SUM(rmsh),         SUM(rmsh**2)        /)
      mat(2,:) = (/ SUM(rmsh),        SUM(rmsh**2),      SUM(rmsh**3)        /)
      mat(3,:) = (/ SUM(rmsh**2),     SUM(rmsh**3),      SUM(rmsh**4)        /)
      rhs      = (/ SUM(vals),        SUM(rmsh*vals),    SUM(rmsh**2*vals)   /)

      CALL solve3x3(mat, rhs, coef)

      quadFitAtR = coef(1) + coef(2)*r0 + coef(3)*r0**2

   END FUNCTION quadFitAtR

   SUBROUTINE solve3x3(a, b, x)
      ! Solves the 3x3 linear system a*x=b via Gaussian elimination with
      ! partial pivoting.
      REAL, INTENT(IN)  :: a(3,3), b(3)
      REAL, INTENT(OUT) :: x(3)

      REAL    :: m(3,4), factor, temp(4)
      INTEGER :: i, k, p

      m(:,1:3) = a
      m(:,4)   = b

      DO k = 1, 2
         p = k
         DO i = k+1, 3
            IF (ABS(m(i,k)) > ABS(m(p,k))) p = i
         END DO
         IF (p /= k) THEN
            temp = m(k,:); m(k,:) = m(p,:); m(p,:) = temp
         END IF
         DO i = k+1, 3
            factor = m(i,k) / m(k,k)
            m(i,:) = m(i,:) - factor*m(k,:)
         END DO
      END DO

      x(3) = m(3,4) / m(3,3)
      x(2) = (m(2,4) - m(2,3)*x(3)) / m(2,2)
      x(1) = (m(1,4) - m(1,2)*x(2) - m(1,3)*x(3)) / m(1,1)

   END SUBROUTINE solve3x3

   CHARACTER(LEN=24) FUNCTION r2s(x)
      ! Formats a real value as an xml attribute string: up to 6 digits
      ! (plus sign) before the decimal separator, 5 digits after it.
      REAL, INTENT(IN) :: x
      r2s = ''
      WRITE(r2s,'(f13.5)') x
   END FUNCTION r2s

   CHARACTER(LEN=24) FUNCTION r2sHP(x)
      ! Higher-precision variant of r2s (7 digits after the decimal
      ! separator instead of 5).
      REAL, INTENT(IN) :: x
      r2sHP = ''
      WRITE(r2sHP,'(f15.7)') x
   END FUNCTION r2sHP

   CHARACTER(LEN=24) FUNCTION i2s(i)
      ! Formats an integer value as an xml attribute string.
      INTEGER, INTENT(IN) :: i
      i2s = ''
      WRITE(i2s,'(i0)') i
   END FUNCTION i2s

END MODULE m_types_moessbauerParams
