MODULE m_totale
CONTAINS
  SUBROUTINE totale(atoms,sphhar,stars,vacuum, &
       sym,input,noco,cell,oneD, xcpot,hybrid, it,results)
    !
    !     ***************************************************
    !     subroutine calculates the total energy of the slab
    !                                  c.l.fu
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
    !     if HF calculation/hybrid-functional calculation :
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
    USE m_constants, ONLY : sfp_const
    USE m_force_a4
    USE m_force_a3
    USE m_forcew
    USE m_loddop
    USE icorrkeys
    USE m_types
    IMPLICIT NONE

    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_xcpot),INTENT(IN)        :: xcpot
    TYPE(t_oneD),INTENT(IN)         :: oneD
    TYPE(t_hybrid),INTENT(IN)       :: hybrid
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_vacuum),INTENT(IN)       :: vacuum
    TYPE(t_noco),INTENT(IN)         :: noco
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_stars),INTENT(IN)        :: stars
    TYPE(t_cell),INTENT(IN)         :: cell
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_atoms),INTENT(IN)        :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT (IN) :: it      
    !     ..
    !     .. Local Scalars ..
    REAL rhs,totz,vmd,zintn_r
    INTEGER n,j,nt,iter,i

    !     .. Local Arrays ..
    REAL dpj(atoms%jmtd)
    !.....density
    REAL,    ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
    COMPLEX, ALLOCATABLE :: qpw(:,:),rhtxy(:,:,:,:)
    !.....potential
    REAL,    ALLOCATABLE :: vr(:,:,:,:),vz(:,:,:)
    COMPLEX, ALLOCATABLE :: vpw(:,:),vxy(:,:,:,:)
    !     ..
    ALLOCATE ( rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins),rht(vacuum%nmzd,2,input%jspins),&
         qpw(stars%n3d,input%jspins),rhtxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,input%jspins),&
         vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins),vz(vacuum%nmzd,2,input%jspins),&
         vpw(stars%n3d,input%jspins),vxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,input%jspins) )
    !
    WRITE (6,FMT=8000)
    WRITE (16,FMT=8000)
8000 FORMAT (/,/,/,5x,'t o t a l  e n e r g y')
    !
    !      ---> sum of eigenvalues (core, semicore and valence states)
    ! 
    results%tote = results%seigscv + results%seigc
    WRITE (6,FMT=8010) results%tote
    WRITE (16,FMT=8010) results%tote
8010 FORMAT (/,10x,'sum of eigenvalues =',t40,f20.10)
    !
    !      ---> add contribution of coulomb potential
    !
    results%tote = results%tote + 0.5e0*results%te_vcoul
    WRITE (6,FMT=8020) results%te_vcoul
    WRITE (16,FMT=8020) results%te_vcoul
8020 FORMAT (/,10x,'density-coulomb potential integral =',t40,f20.10)
    !
    !      ---> add contribution of effective potential
    !
    results%tote = results%tote - results%te_veff
    WRITE (6,FMT=8030) results%te_veff
    WRITE (16,FMT=8030) results%te_veff
8030 FORMAT (/,10x,'density-effective potential integral =',t40,f20.10)
    !
    !      ---> add contribution of exchange-correlation energy
    !
    results%tote = results%tote + results%te_exc
    WRITE (6,FMT=8040) results%te_exc
    WRITE (16,FMT=8040) results%te_exc
8040 FORMAT (/,10x,'charge density-ex.-corr.energy density integral=', t40,f20.10)
    !
    !      ---> Fock exchange contribution 
    !
    IF( xcpot%icorr .EQ. icorr_hf .OR. xcpot%icorr .EQ. icorr_pbe0 &
             .OR. xcpot%icorr .EQ. icorr_hse .OR. xcpot%icorr .EQ. icorr_vhse ) THEN
       results%tote = results%tote - 0.5e0*results%te_hfex%valence + 0.5e0*results%te_hfex%core
    ELSE IF ( xcpot%icorr .EQ. icorr_exx ) THEN
       results%tote = results%tote + 0.5e0*results%te_hfex%valence
    END IF
    WRITE (6,FMT=8100)  0.5e0*results%te_hfex%valence
    WRITE (16,FMT=8100) 0.5e0*results%te_hfex%valence
    WRITE (6,FMT=8101)  0.5e0*results%te_hfex%core
    WRITE (16,FMT=8101) 0.5e0*results%te_hfex%core
8100 FORMAT (/,10x,'Fock-exchange energy (valence)=',t40,f20.10)
8101 FORMAT (10x,'Fock-exchange energy (core)=',t40,f20.10)


    !     ----> VM terms
    !     ---> reload the density
    !
    IF (noco%l_noco) THEN
       nt = 70
       OPEN (nt,file='cdn',form='unformatted',status='old')
    ELSE
       nt = 71
       OPEN (nt,file='cdn1',form='unformatted',status='old')
    ENDIF
    CALL loddop(stars,vacuum,atoms,sphhar, input,sym,&
                     nt, iter,rho,qpw,rht,rhtxy)
    CLOSE (nt)
    !+for
    !     ---> reload the COULOMB potential
    !
    OPEN (11,file='potcoul',form='unformatted',status='old')
    REWIND 11
    CALL loddop(stars,vacuum,atoms,sphhar, input,sym,&
                     11, iter,vr,vpw,vz,vxy)
    CLOSE (11)
    !
    !     CLASSICAL HELLMAN-FEYNMAN FORCE
    !
    CALL force_a3(atoms,sphhar, input, rho,vr, results%force)
    !
    IF (input%l_f) THEN
       !
       !       core contribution to force: needs TOTAL POTENTIAL and core charge
       OPEN (8,file='pottot',form='unformatted',status='old')
       REWIND 8
       CALL loddop(stars,vacuum,atoms,sphhar, input,sym,&
                          8, iter,vr,vpw,vz,vxy)
       CLOSE (8)
       !
       CALL force_a4(atoms,sphhar,input, vr, results%force)
       !
    ENDIF
    !

    !-for
    !     ---> add spin-up and spin-down charge density for lh=0
    !
    IF (input%jspins.EQ.2) THEN
       DO  n = 1,atoms%ntype
          DO  i = 1,atoms%jri(n)
             rho(i,0,n,1) = rho(i,0,n,1) + rho(i,0,n,input%jspins)
          ENDDO
       ENDDO
    END IF
    !
    ! ----> coulomb interaction between electrons and nuclei of different m.t.s
    !
    DO  n = 1,atoms%ntype
       DO  j = 1,atoms%jri(n)
          dpj(j) = rho(j,0,n,1)/atoms%rmsh(j,n)
       ENDDO
       CALL intgr3(dpj,atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),rhs)
       !
       results%tote = results%tote - atoms%neq(n)*atoms%zatom(n)*sfp_const*rhs/2.
       !
       zintn_r = atoms%neq(n)*atoms%zatom(n)*sfp_const*rhs/2.
       WRITE (6,FMT=8045) zintn_r
       WRITE (16,FMT=8045) zintn_r
       CALL intgr3(rho(1,0,n,1),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),totz)
       vmd = atoms%rmt(n)*atoms%vr0(n)/sfp_const + atoms%zatom(n) - totz*sfp_const
       vmd = -atoms%neq(n)*atoms%zatom(n)*vmd/ (2.*atoms%rmt(n))
       WRITE (6,FMT=8050) n,vmd
       WRITE (16,FMT=8050) n,vmd
       results%tote = results%tote + vmd
    ENDDO
    IF (atoms%n_u.GT.0) THEN
       WRITE ( 6,FMT=8090) results%e_ldau
       WRITE (16,FMT=8090) results%e_ldau
       results%tote = results%tote - results%e_ldau             ! gu test
    ENDIF
    ! print 'HF' before total energy to make it grepable
    IF ( .NOT. hybrid%l_calhf ) THEN
       WRITE ( 6,FMT=8060) results%tote
       WRITE (16,FMT=8060) results%tote
    ELSE
       WRITE ( 6,FMT=8061) results%tote
       WRITE (16,FMT=8061) results%tote
    END IF

    CALL force_w(input,atoms,sym,results,cell,oneD)
    !
    !     ---> calculate the free energy and the ground state energy,
    !          extrapolated for T->0
    !
    ! print 'HF' before all energies to make them grepable
    IF ( .NOT. hybrid%l_calhf ) THEN
       WRITE ( 6,FMT=8065) results%ts
       WRITE (16,FMT=8065) results%ts
       WRITE ( 6,FMT=8070) results%tote-results%ts
       WRITE (16,FMT=8070) results%tote-results%ts
       WRITE ( 6,FMT=8080) results%tote-0.5e0*results%ts
       WRITE (16,FMT=8080) results%tote-0.5e0*results%ts
    ELSE
       WRITE ( 6,FMT=8066) results%ts
       WRITE (16,FMT=8066) results%ts
       WRITE ( 6,FMT=8071) results%tote-results%ts
       WRITE (16,FMT=8071) results%tote-results%ts
       WRITE ( 6,FMT=8081) results%tote-0.5e0*results%ts
       WRITE (16,FMT=8081) results%tote-0.5e0*results%ts
    END IF
8060 FORMAT (/,/,' ---->    input%total energy=',t40,f20.10,' htr')
8061 FORMAT (/,/,' ----> HF input%total energy=',t40,f20.10,' htr')
8050 FORMAT (/,10x,'Madelung term for atom type:',i3,t40,f20.10)
8045 FORMAT (/,10x,'el.-nucl. inter. diff. m.t.',t40,f20.10)
8065 FORMAT (/,/,' ---->    (input%tkb*entropy) TS=',t40,f20.10,' htr')
8066 FORMAT (/,/,' ----> HF (input%tkb*entropy) TS=',t40,f20.10,' htr')
8070 FORMAT (/,/,' ---->    free energy=',t40,f20.10,' htr')
8071 FORMAT (/,/,' ----> HF free energy=',t40,f20.10,' htr')
8080 FORMAT (/,/,'      extrapolation for T->0',&
               /,' ---->    input%total electron energy=',t40,f20.10,' htr')
8081 FORMAT (/,/,'      extrapolation for T->0',&
               /,' ----> HF input%total electron energy=',t40,f20.10,' htr')
8090 FORMAT (/,/,' ---->    correction for lda+U =',t40,f20.10,' htr')

    DEALLOCATE ( rho,rht,qpw,rhtxy,vr,vz,vpw,vxy )

  END SUBROUTINE totale
END MODULE m_totale
