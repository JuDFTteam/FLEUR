!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_totale
CONTAINS
  SUBROUTINE totale(fmpi,atoms,sphhar,stars,vacuum, &
       sym,input,noco,cell,oneD, xcpot,hybdat,vTot,vCoul,it,den,results)
    !
    !     ***************************************************
    !     subroutine calculates the total energy
    !     ***************************************************
    !     single particle energies
    !     SEIGC  sum of the eigenvalues of the core states
    !            calculated in cdngen.f
    !     SEIGSCV  sum of the eigenvalues of the semicore and valence states
    !              calculated in fermie.f
    !     TS         : entropy contribution to the free energy
    !     SEIGC,SEIGSCV, TS are calculated in fermie.f
    !     ***************************************************
    !     TE_VCOUL  :   charge density-coulomb potential integral
    !     TE_VEFF:   charge density-effective potential integral
    !     TE_EXC :   charge density-ex-corr.energy density integral
    !                 exchange-correlation energy
    !     TE_VCOUL,TE_VEFF,TE_EXC are calculated in vgen.f
    !     VMD :   Madelung term
    !     ***************************************************
    !     TOTE    :   total energy due to all electrons
    !     TOTE = SEIGC + SEIGSCV + TE_VCOUL/2 -TE_VEFF + TE_EXC + VMD
    !
    !     if HF calculation/hybinp-functional calculation :
    !     TOTE = SEIGC + SEIGSCV + TE_VCOUL/2 -TE_VEFF + TE_EXC_loc + VMD - 1/2 E_FOCK
    !
    !     E_FOCK: sum of diagonal elements of fock matrix
    !
    !     ***************************************************
    !     FREE ENRGY: F = TOTE - TS
    !     total electron energy extrapolated for T->0
    !     E0 = TOTE - TS/2
    !     ***************************************************
    !
    USE m_intgr    , ONLY : intgr3
    USE m_constants
    USE m_force_a4
    USE m_force_a3
    USE m_force_a4_add ! Klueppelberg (force level 1)
    USE m_force_sf ! Klueppelberg (force level 3)
    USE m_forcew
    USE m_cdn_io
    USE m_types
    USE m_xmlOutput
    use m_judft
    USE m_vdWfleur_grimme
    
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)          :: fmpi
    TYPE(t_results),INTENT(INOUT)   :: results
    CLASS(t_xcpot),INTENT(IN)       :: xcpot
    TYPE(t_oneD),INTENT(IN)         :: oneD
    TYPE(t_hybdat),INTENT(IN)       :: hybdat
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_vacuum),INTENT(IN)       :: vacuum
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_stars),INTENT(IN)        :: stars
    TYPE(t_cell),INTENT(IN)         :: cell
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_atoms),INTENT(IN)        :: atoms

    TYPE(t_potden),INTENT(IN)       :: vTot,vCoul
    TYPE(t_potden),INTENT(IN)       :: den
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT (IN) :: it

    ! Local type instances

    !     .. Local Scalars ..
    REAL rhs,totz, eigSum, fermiEnergyTemp
    INTEGER n,j,nt,i, archiveType,jsp
    LOGICAL l_qfix

    !     .. Local Arrays ..
    REAL vmd(atoms%ntype),zintn_r(atoms%ntype)
    REAL dpj(atoms%jmtd),mt(atoms%jmtd,atoms%ntype)
    CHARACTER(LEN=20) :: attributes(3)

    !CALL den%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,.FALSE.,POTDEN_TYPE_DEN)
    IF (fmpi%irank==0) THEN
       WRITE (oUnit,FMT=8000)
8000   FORMAT (/,/,/,5x,'t o t a l  e n e r g y')
       !
       !      ---> sum of eigenvalues (core, semicore and valence states)
       !
       eigSum = results%seigscv + results%seigc
       results%tote = eigSum
       WRITE (oUnit,FMT=8010) results%tote
8010   FORMAT (/,10x,'sum of eigenvalues =',t40,f20.10)
       !
       !      ---> add contribution of coulomb potential
       !
       results%tote = results%tote + 0.5e0*results%te_vcoul
       WRITE (oUnit,FMT=8020) results%te_vcoul
8020   FORMAT (/,10x,'density-coulomb potential integral =',t40,f20.10)
       !
       !      ---> add contribution of effective potential
       !
       results%tote = results%tote - results%te_veff
       WRITE (oUnit,FMT=8030) results%te_veff
8030   FORMAT (/,10x,'density-effective potential integral =',t40,f20.10)
       !
       !      ---> add contribution of exchange-correlation energy
       !
       results%tote = results%tote + results%te_exc
       WRITE (oUnit,FMT=8040) results%te_exc
8040   FORMAT (/,10x,'charge density-ex.-corr.energy density integral=', t40,f20.10)
       !
       !      ---> Fock exchange contribution
       !
       IF (xcpot%is_hybrid()) THEN
          !IF (xcpot%is_name("exx")) THEN
          !   results%tote = results%tote + 0.5e0*results%te_hfex%valence
          !ELSE
          results%tote = results%tote - 0.5e0*results%te_hfex%valence + 0.5e0*results%te_hfex%core
          !END IF
          write (oUnit,*)  'Fock-exchange energy (valence)= ' // float2str(0.5e0*results%te_hfex%valence)
          write (oUnit,*)  'Fock-exchange energy (core)=    ' // float2str(0.5e0*results%te_hfex%core)
       ENDIF


       !     ----> VM terms
       !     ---> reload the density
       !
       !archiveType = CDN_ARCHIVE_TYPE_CDN1_const
       !IF (noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_CDN_const

       !CALL readDensity(stars,noco,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,&
       !                 CDN_INPUT_DEN_const,0,fermiEnergyTemp,l_qfix,den)


       ! CLASSICAL HELLMAN-FEYNMAN FORCE
       CALL force_a3(atoms,sym,sphhar, input, den%mt,vCoul%mt, results%force)

       IF (input%l_f) THEN
          ! core contribution to force: needs TOTAL POTENTIAL and core charge

          IF (input%ctail.AND.(input%f_level.GE.1)) THEN
             ! Add core correction to forces from tails of core states
             ! Klueppelberg, Sep'12 (force level 1)
             CALL force_a4_add(atoms,input,results)
          END IF

          CALL force_a4(atoms,sym,sphhar,input, vTot%mt, results%force)

       END IF

       !-for
       !     ---> add spin-up and spin-down charge density for lh=0
       !
       mt=0.0
       DO  n = 1,atoms%ntype
          DO  i = 1,atoms%jri(n)
             mt(i,n) = den%mt(i,0,n,1) + den%mt(i,0,n,input%jspins)
          ENDDO
       ENDDO
       IF (input%jspins.EQ.1) mt=mt/2 !we just added the same value twice

       !
       ! ----> coulomb interaction between electrons and nuclei of different m.t.s
       !
       DO  n = 1,atoms%ntype
          DO  j = 1,atoms%jri(n)
             dpj(j) = mt(j,n)/atoms%rmsh(j,n)
          ENDDO
          CALL intgr3(dpj,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),rhs)
          !
          zintn_r(n) = atoms%neq(n)*atoms%zatom(n)*sfp_const*rhs/2.
          results%tote = results%tote - zintn_r(n)
          !
          WRITE (oUnit,FMT=8045) zintn_r(n)
          CALL intgr3(mt(1,n),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),totz)
          vmd(n) = atoms%rmt(n)*vCoul%mt(atoms%jri(n),0,n,1)/sfp_const + atoms%zatom(n) - totz*sfp_const
          vmd(n) = -atoms%neq(n)*atoms%zatom(n)*vmd(n)/ (2.*atoms%rmt(n))
          WRITE (oUnit,FMT=8050) n,vmd(n)
          results%tote = results%tote + vmd(n)
       ENDDO
       IF (atoms%n_u+atoms%n_hia.GT.0) THEN
          WRITE (oUnit,FMT=8090) results%e_ldau
          results%tote = results%tote - results%e_ldau             ! gu test
       ENDIF
       ! print 'HF' before total energy to make it grepable
       IF ( .NOT. hybdat%l_calhf ) THEN
          WRITE (oUnit,FMT=8060) results%tote
       ELSE
          WRITE (oUnit,FMT=8061) results%tote
       END IF
       !
       !     ---> calculate the free energy and the ground state energy,
       !          extrapolated for T->0
       !
       ! print 'HF' before all energies to make them grepable
       IF ( .NOT. hybdat%l_calhf ) THEN
          WRITE (oUnit,FMT=8065) results%ts
          WRITE (oUnit,FMT=8070) results%tote-results%ts
          WRITE (oUnit,FMT=8080) results%tote-0.5e0*results%ts
       ELSE
          WRITE (oUnit,FMT=8066) results%ts
          WRITE (oUnit,FMT=8071) results%tote-results%ts
          WRITE (oUnit,FMT=8081) results%tote-0.5e0*results%ts
       END IF

       ! vdW D3 Grimme contribution
       IF (btest(input%vdW,0)) THEN
         if (.not.allocated(results%force_vdw)) ALLOCATE(results%force_vdW(3,atoms%ntype))
         call vdW_fleur_grimme(input,atoms,sym,cell,results%e_vdW,results%force_vdw)
      ENDIF
      !ADD vdW contribution to energy (is zero if no vdW was evaluated)
      results%tote=results%tote+results%e_vdW
      WRITE(attributes(1),'(f20.10)') results%tote
       WRITE(attributes(2),'(a)') 'Htr'
       WRITE(attributes(3),'(a)') 'HF'
       IF (hybdat%l_calhf) THEN
          CALL openXMLElementForm('totalEnergy',(/'value  ','units  ','comment'/),attributes,reshape((/40,20/),(/1,2/)))
       ELSE
          CALL openXMLElementForm('totalEnergy',(/'value','units'/),attributes(1:2),reshape((/40,20/),(/1,2/)))
       END IF
       CALL openXMLElementFormPoly('sumOfEigenvalues',(/'value'/),(/eigSum/),reshape((/32,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('coreElectrons',(/'value'/),(/results%seigc/),reshape((/32,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('valenceElectrons',(/'value'/),(/results%seigscv/),reshape((/29,20/),(/1,2/)))
       CALL closeXMLElement('sumOfEigenvalues')
       CALL writeXMLElementFormPoly('densityCoulombPotentialIntegral',(/'value'/),(/results%te_vcoul/),reshape((/17,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('densityEffectivePotentialIntegral',(/'value'/),(/results%te_veff/),reshape((/15,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('chargeDenXCDenIntegral',(/'value'/),(/results%te_exc/),reshape((/26,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('FockExchangeEnergyValence',(/'value'/),(/0.5e0*results%te_hfex%valence/),reshape((/23,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('FockExchangeEnergyCore',(/'value'/),(/0.5e0*results%te_hfex%core/),reshape((/26,20/),(/1,2/)))
       if (btest(input%vdw,0).or.btest(input%vdW,1)) call writeXMLElementFormPoly('vdWEnergy',(/'value'/),(/results%e_vdW/),reshape((/17,20/),(/1,2/)))
       DO  n = 1,atoms%ntype
          CALL openXMLElementPoly('atomTypeDependentContributions',(/'atomType'/),(/n/))
          CALL writeXMLElementFormPoly('electronNucleiInteractionDifferentMTs',(/'value'/),(/zintn_r(n)/),reshape((/8,20/),(/1,2/)))
          CALL writeXMLElementFormPoly('MadelungTerm',(/'value'/),(/vmd(n)/),reshape((/33,20/),(/1,2/)))
          CALL closeXMLElement('atomTypeDependentContributions')
       END DO
       IF (atoms%n_u+atoms%n_hia.GT.0) THEN
          CALL writeXMLElementFormPoly('dftUCorrection',(/'value'/),(/results%e_ldau/),reshape((/34,20/),(/1,2/)))
       END IF
       CALL writeXMLElementFormPoly('tkbTimesEntropy',(/'value'/),(/results%ts/),reshape((/33,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('freeEnergy',(/'value'/),(/results%tote-results%ts/),reshape((/38,20/),(/1,2/)))
       CALL writeXMLElementFormPoly('extrapolationTo0K',(/'value'/),(/results%tote-0.5e0*results%ts/),reshape((/31,20/),(/1,2/)))
       CALL closeXMLElement('totalEnergy')
8060   FORMAT (/,/,' ---->    total energy=',t40,f20.10,' htr')
8061   FORMAT (/,/,' ----> HF total energy=',t40,f20.10,' htr')
8050   FORMAT (/,10x,'Madelung term for atom type:',i3,t40,f20.10)
8045   FORMAT (/,10x,'el.-nucl. inter. diff. m.t.',t40,f20.10)
8065   FORMAT (/,/,' ---->    (input%tkb*entropy) TS=',t40,f20.10,' htr')
8066   FORMAT (/,/,' ----> HF (input%tkb*entropy) TS=',t40,f20.10,' htr')
8070   FORMAT (/,/,' ---->    free energy=',t40,f20.10,' htr')
8071   FORMAT (/,/,' ----> HF free energy=',t40,f20.10,' htr')
8080   FORMAT (/,/,'      extrapolation for T->0',&
            /,' ---->    total electron energy=',t40,f20.10,' htr')
8081   FORMAT (/,/,'      extrapolation for T->0',&
            /,' ----> HF total electron energy=',t40,f20.10,' htr')
8090   FORMAT (/,/,' ---->    correction for lda+U =',t40,f20.10,' htr')
    ENDIF

    ! Klueppelberg (force level 3)
    IF (input%l_f.AND.(input%f_level.GE.3)) THEN 
       DO jsp=1,input%jspins
          CALL exit_sf(jsp,atoms,results%force)
       END DO
    END IF
    CALL force_w(fmpi,input,atoms,sym,results,cell,oneD,vacuum)

  END SUBROUTINE totale
END MODULE m_totale
