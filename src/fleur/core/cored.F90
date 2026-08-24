MODULE m_cored
CONTAINS
   SUBROUTINE cored(input, jspin, iType, atoms, rho,  sphhar, l_CoreDenPresent, vr, qint, rhc, tec, seig, l_useOtherCoreSolver, EnergyDen, &
                     moessbauerParams)
      !     *******************************************************
      !     *****   set up the core densities for compounds.  *****
      !     *****                      d.d.koelling           *****
      !     *******************************************************
      USE m_juDFT
      USE m_intgr, ONLY : intgr3,intgr0,intgr1,intgr2
      USE m_constants
      !USE m_setcor
      USE m_differ
      USE m_types
      USE m_types_moessbauerParams
      USE m_cdn_io
      USE m_xmlOutput
      IMPLICIT NONE

      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_sphhar),INTENT(IN)      :: sphhar
      TYPE(t_atoms),INTENT(IN)       :: atoms
      !
      !     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: jspin, iType
      LOGICAL, INTENT (IN) :: l_CoreDenPresent
      REAL,    INTENT (OUT) :: seig
      !     ..
      !     .. Array Arguments ..
      REAL, INTENT(IN)              :: vr(atoms%jmtd,atoms%ntype)
      REAL, INTENT(INOUT)           :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      REAL, INTENT(INOUT)           :: rhc(atoms%msh,atoms%ntype,input%jspins)
      REAL, INTENT(INOUT)           :: qint(atoms%ntype,input%jspins)
      REAL, INTENT(INOUT)           :: tec(atoms%ntype,input%jspins)
      LOGICAL, INTENT (INOUT), OPTIONAL :: l_useOtherCoreSolver
      REAL, INTENT(INOUT), OPTIONAL :: EnergyDen(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins)
      TYPE(t_moessbauerParams), OPTIONAL, INTENT(INOUT) :: moessbauerParams

      !     ..
      !     .. Local Scalars ..
      REAL eig,fj,fl,fn,qOutside,rad,rhos,rhs,sea,sume,t2
      REAL d,dxx,rn,rnot,z,t1,rr,r,lambd,c,bmu,weight, aux_weight
      INTEGER i,j,korb,ncmsh,nm,nm1,nst ,l,ierr
      INTEGER ilshell, nShQ, lShQ, n_prev, l_prev
      INTEGER iRad
      LOGICAL isSmaller, l_isoShiftOK
      REAL nucRad, alpha, smallIsoVal, largeIsoVal, stateIsoShift
      CHARACTER(LEN=150) :: warnMsg
      !     ..
      !     .. Local Arrays ..

      REAL rhcs(atoms%msh),rhoc(atoms%msh),rhoss(atoms%msh),vrd(atoms%msh),f(0:3)
      REAL rhcs_aux(atoms%msh), rhoss_aux(atoms%msh) !> quantities for energy density calculations
      REAL occ(maxval(atoms%econf%num_states)),a(atoms%msh),b(atoms%msh),ain(atoms%msh),ahelp(atoms%msh)
      REAL occ_h(maxval(atoms%econf%num_states),2)
      INTEGER kappa(maxval(atoms%econf%num_states)),nprnc(maxval(atoms%econf%num_states))
      CHARACTER(LEN=20) :: attributes(6)
      REAL stateEnergies(29)
      LOGICAL l_useThisCoreSolver
      INTEGER shellN(15), shellL(15)
      REAL shellIsomerShift(15)
      REAL indefInteg(atoms%msh)
      !     ..

      l_useThisCoreSolver = .TRUE.
      IF(PRESENT(l_useOtherCoreSolver)) THEN
         l_useThisCoreSolver = l_useOtherCoreSolver
      END IF

      c = c_light(1.0)
      seig = 0.0
      !
      IF(.NOT.l_useThisCoreSolver) THEN
         IF (PRESENT(moessbauerParams)) moessbauerParams%nCoreShells(iType) = 0
         RETURN
      END IF

      IF (input%frcor.and.l_CoreDenPresent) THEN
         rnot = atoms%rmsh(1,iType) ; dxx = atoms%dx(iType)
         ncmsh = NINT( LOG( (atoms%rmt(iType)+10.0)/rnot ) / dxx + 1 )
         ncmsh = MIN( ncmsh, atoms%msh )
         !     --->    update spherical charge density
         DO i = 1,atoms%jri(iType)
            rhoc(i) = rhc(i,iType,jspin)
            rho(i,0,iType,jspin) = rho(i,0,iType,jspin) + rhoc(i)/sfp_const
         ENDDO
         !     ---> for total energy calculations, determine the sum of the
         !     ---> eigenvalues by requiring that the core kinetic energy
         !     ---> remains constant.
         DO i = 1,atoms%jri(iType)
            rhoc(i) = rhoc(i)*vr(i,iType)/atoms%rmsh(i,iType)
         ENDDO
         nm = atoms%jri(iType)
         CALL intgr3(rhoc,atoms%rmsh(1,iType),atoms%dx(iType),nm,rhos)
         sea = tec(iType,jspin) + rhos
         WRITE (oUnit,FMT=8030) iType,jspin,tec(iType,jspin),sea
         seig = seig + atoms%neq(iType)*sea
      END IF
      
      IF (input%frcor.and.l_CoreDenPresent) THEN
         IF (PRESENT(moessbauerParams)) moessbauerParams%nCoreShells(iType) = 0
         RETURN
      END IF

      !     ---> set up densities

      sume = 0.
      z = atoms%zatom(iType)
      !         rn = rmt(iType)
      dxx = atoms%dx(iType)

      ! Locate the mesh bracket around the nuclear radius, for the
      ! shell-wise isomer shift below (same 2-point linear interpolation as
      ! the kcrel=1 core solver in core.f90).
      nucRad = r0_const*(atomicMasses_const(atoms%nz(iType))**(1.0/3.0))
      l_isoShiftOK = (atoms%rmsh(1,iType) < nucRad)
      IF (.NOT.l_isoShiftOK) THEN
         IF (jspin==1) THEN
            WRITE(warnMsg,'(a,i0,a)') 'Nuclear radius smaller than smallest radial mesh point for atom type ', &
                                       iType, ' -- skipping its core isomer shift'
            WRITE(*,'(a)') TRIM(warnMsg)
         END IF
      ELSE
         iRad = 0
         isSmaller = .TRUE.
         DO WHILE (isSmaller)
            iRad = iRad + 1
            IF (atoms%rmsh(iRad,iType).GE.nucRad) isSmaller = .FALSE.
         END DO
         alpha = (nucRad - atoms%rmsh(iRad-1,iType)) / (atoms%rmsh(iRad,iType) - atoms%rmsh(iRad-1,iType))
      END IF
      bmu = 0.0
      !CALL setcor(iType,input%jspins,atoms,input,bmu,nst,kappa,nprnc,occ_h)
      CALL atoms%econf(iType)%get_core(nst,nprnc,kappa,occ_h)


      occ(1:nst) = occ_h(1:nst,jspin)
         
      rnot = atoms%rmsh(1,iType)
      d = EXP(atoms%dx(iType))
      ncmsh = NINT( LOG( (atoms%rmt(iType)+10.0)/rnot ) / dxx + 1 )
      ncmsh = MIN( ncmsh, atoms%msh )
      rn = rnot* (d** (ncmsh-1))
      WRITE (oUnit,FMT=8000) z,rnot,dxx,atoms%jri(iType)
      DO j = 1,atoms%jri(iType)
         rhoss(j)     = 0.0
         if(present(EnergyDen)) rhoss_aux(j) = 0.0
         vrd(j) = vr(j,iType)
      ENDDO
      !
      IF (input%l_core_confpot) THEN
         !--->    linear extension of the potential with slope t1 / a.u.
         t1=0.125
         t1 = MAX( (vrd(atoms%jri(iType)) - vrd(atoms%jri(iType)-1)*d)*&
                  d / (atoms%rmt(iType)**2 * (d-1) ) , t1)
         t2=vrd(atoms%jri(iType))/atoms%rmt(iType)-atoms%rmt(iType)*t1
         rr = atoms%rmt(iType)
      ELSE
         t2 = vrd(atoms%jri(iType)) / ( atoms%jri(iType) - ncmsh )
      ENDIF
      IF (atoms%jri(iType) < ncmsh) THEN
         DO i = atoms%jri(iType) + 1,ncmsh
            rhoss(i) = 0.
            if(present(EnergyDen)) rhoss_aux(i) = 0.0
            IF (input%l_core_confpot) THEN
               rr = d*rr
               vrd(i) = rr*( t2 + rr*t1 )
               !               vrd(i) = 2*vrd(jri(iType)) - rr*( t2 + rr*t1 )
            ELSE
               vrd(i) = vrd(atoms%jri(iType)) + t2* (i-atoms%jri(iType))
            ENDIF
            !
         ENDDO
      END IF

      nst = atoms%econf(iType)%num_core_states        ! for lda+U

      IF (input%gw==1 .OR. input%gw==3) WRITE(15) nst,atoms%rmsh(1:atoms%jri(iType),iType)

      stateEnergies = 0.0
      ilshell = 0
      n_prev = -999
      l_prev = -999
      shellN = 0
      shellL = 0
      shellIsomerShift = 0.0
      DO korb = 1,nst
         IF (occ(korb) /= 0.0) THEN
            fn = nprnc(korb)
            fj = iabs(kappa(korb)) - .5e0

!            weight = 2*fj + 1.e0
!            IF (bmu > 99.) weight = occ(korb)
            weight = 2*occ(korb)

            fl = fj + (.5e0)*isign(1,kappa(korb))

            eig        = -2* (z/ (fn+fl))**2

            CALL differ(fn,fl,fj,c,z,dxx,rnot,rn,d,ncmsh,vrd, eig, a,b,ierr)
            stateEnergies(korb) = eig
            WRITE (oUnit,FMT=8010) fn,fl,fj,eig,weight

            IF (ierr/=0)  CALL juDFT_error("error in core-level routine" ,calledby ="cored")
            IF (input%gw==1 .OR. input%gw==3) WRITE (15) NINT(fl),weight,eig,&
               a(1:atoms%jri(iType)),b(1:atoms%jri(iType))

            sume = sume + weight*eig/input%jspins
            DO j = 1,ncmsh
               rhcs(j)  = weight* (a(j)**2+b(j)**2)
               rhoss(j) = rhoss(j) + rhcs(j)
            ENDDO

            nShQ = nprnc(korb)
            lShQ = NINT(fl)
            IF (nShQ /= n_prev .OR. lShQ /= l_prev) THEN
               ilshell = ilshell + 1
               shellN(ilshell) = nShQ
               shellL(ilshell) = lShQ
               n_prev = nShQ
               l_prev = lShQ
            END IF

            ! Shell-wise isomer shift: average density enclosed within the
            ! nuclear sphere, from this state's own density rhcs, using the
            ! same 2-point linear bracket interpolation at r=nucRad as the
            ! kcrel=1 core solver (core.f90). Same /jspins normalization as
            ! used above for rhoc=rhoss/input%jspins. Skipped entirely if the
            ! nuclear radius is smaller than the innermost mesh point (see
            ! l_isoShiftOK above).
            IF (l_isoShiftOK) THEN
               CALL intgr2(rhcs(1:atoms%jri(iType)), atoms%rmsh(1:atoms%jri(iType),iType), atoms%dx(iType), &
                           atoms%jri(iType), indefInteg(1:atoms%jri(iType)))
               smallIsoVal = indefInteg(iRad-1) / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad-1,iType)**3.0))
               largeIsoVal = indefInteg(iRad) / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad,iType)**3.0))
               stateIsoShift = (1.0-alpha)*smallIsoVal + alpha*largeIsoVal
               shellIsomerShift(ilshell) = shellIsomerShift(ilshell) + stateIsoShift/input%jspins
            END IF

            IF(present(EnergyDen)) THEN
               !rhoss_aux = rhoss
               DO j = 1,ncmsh
                  ! for energy density we want to multiply the weights
                  ! with the eigenenergies
                  rhoss_aux(j) = rhoss_aux(j) + (rhcs(j) * eig)
               ENDDO
            ENDIF
         ENDIF
      ENDDO

      IF (PRESENT(moessbauerParams)) THEN
         CALL moessbauerParams%accumCoreIS(iType, ilshell, shellN, shellL, shellIsomerShift, l_isoShiftOK)
      END IF

      !     ---->update spherical charge density rho with the core density.
      !     ---->for spin-polarized (jspins=2), take only half the density
      nm = atoms%jri(iType)
      DO  j = 1,nm
         rhoc(j) = rhoss(j)/input%jspins
         rho(j,0,iType,jspin) = rho(j,0,iType,jspin) + rhoc(j)/sfp_const
      ENDDO

      IF(present(EnergyDen)) then
         DO j = 1,nm
            EnergyDen(j,0,iType,jspin) = EnergyDen(j,0,iType,jspin) &
                                        + rhoss_aux(j) /(input%jspins * sfp_const)
         ENDDO
      ENDIF

      rhc(1:ncmsh,iType,jspin) = rhoss(1:ncmsh) / input%jspins
      rhc(ncmsh+1:atoms%msh,iType,jspin) = 0.0

      seig = seig + atoms%neq(iType)*sume
      DO i = 1,nm
         rhoc(i) = rhoc(i)*vr(i,iType)/atoms%rmsh(i,iType)
      ENDDO
      CALL intgr3(rhoc,atoms%rmsh(1,iType),atoms%dx(iType),nm,rhs)
      tec(iType,jspin) = sume - rhs
      WRITE (oUnit,FMT=8030) iType,jspin,tec(iType,jspin),sume

      !     ---> simpson integration
      rad = atoms%rmt(iType)
      ! qOutside is the charge outside a single MT sphere of the considered atom type
      qOutside = rad*rhoss(nm)/2.
      DO nm1 = nm + 1,ncmsh - 1,2
         rad = d*rad
         qOutside = qOutside + 2*rad*rhoss(nm1)
         rad = d*rad
         qOutside = qOutside + rad*rhoss(nm1+1)
      ENDDO
      qOutside = 2*qOutside*dxx/3
      !+sb
      WRITE (oUnit,FMT=8020) qOutside/input%jspins
      !-sb
      qint(iType,jspin) = qOutside*atoms%neq(iType)
      attributes = ''
      WRITE(attributes(1),'(i0)') iType
      WRITE(attributes(2),'(i0)') NINT(z)
      WRITE(attributes(3),'(i0)') jspin
      WRITE(attributes(4),'(f18.10)') tec(iType,jspin)
      WRITE(attributes(5),'(f18.10)') sume
      WRITE(attributes(6),'(f9.6)') qOutside/input%jspins
      CALL openXMLElementForm('coreStates',(/'atomType     ','atomicNumber ','spin         ','kinEnergy    ',&
                                             'eigValSum    ','lostElectrons'/),&
                              attributes,RESHAPE((/8,12,4,9,9,13,6,3,1,18,18,9/),(/6,2/)))
      DO korb = 1, atoms%econf(iType)%num_core_states
         fj = iabs(kappa(korb)) - .5e0
!         weight = 2*fj + 1.e0
!         IF (bmu > 99.) weight = occ(korb)
         weight = occ(korb)
         fl = fj + (.5e0)*isign(1,kappa(korb))
         attributes = ''
         WRITE(attributes(1),'(i0)') nprnc(korb)
         WRITE(attributes(2),'(i0)') NINT(fl)
         WRITE(attributes(3),'(f4.1)') fj
         WRITE(attributes(4),'(f20.10)') stateEnergies(korb)
         WRITE(attributes(5),'(f15.10)') weight
         CALL writeXMLElementForm('state',(/'n     ','l     ','j     ','energy','weight'/),&
                                  attributes(1:5),RESHAPE((/1,1,1,6,6,2,2,4,20,15/),(/5,2/)))
      END DO
      CALL closeXMLElement('coreStates')

      RETURN

8000  FORMAT (/,/,10x,'z=',f4.0,5x,'r(1)=',e14.6,5x,'dx=',f9.6,5x,&
          &       'm.t.index=',i4,/,15x,'n',4x,'l',5x,'j',4x,'energy',7x,&
          &       'weight')
8010  FORMAT (12x,2f5.0,f6.1,f10.4,f12.4)
8020  FORMAT (f20.8,'  electrons lost from core.')
8030  FORMAT (10x,'atom type',i5,'  (spin',i2,') ',/,10x,&
          &       'kinetic energy=',e20.12,5x,'sum of the eigenvalues=',&
          &       e20.12)
   END SUBROUTINE cored
END MODULE m_cored
