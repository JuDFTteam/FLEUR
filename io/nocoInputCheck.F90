!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_nocoInputCheck

   CONTAINS

   SUBROUTINE nocoInputCheck(atoms,input,sym,vacuum,noco)

      USE m_juDFT
      USE m_constants
      USE m_types_atoms
      USE m_types_input
      USE m_types_sym
      USE m_types_vacuum
      USE m_types_noco
      USE m_sssym

      IMPLICIT NONE

      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_noco),   INTENT(IN)    :: noco

      INTEGER itype
      LOGICAL l_relax_any,error(sym%nop)

!---> make sure second variation is switched off
      IF (input%secvar) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the second variation!!'
         CALL juDFT_error("Second variation cannot be used!!!" ,calledby="nocoInputCheck")
      END IF

!---> make sure histogram method is used
      IF (input%bz_integration==BZINT_METHOD_GAUSS) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the Gaussian smearing for '
         WRITE (oUnit,*) 'the Brillouin zone integration!!'
         WRITE (oUnit,*) 'Please use the histogram method.'
         CALL juDFT_error("Only histogram Brillouin zone integration can be used!!!",calledby ="nocoInputCheck")
      END IF

!---> make sure force is switched off
      IF (input%l_f) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'does not support force calculations.'
         CALL juDFT_error("force calculations not supported!!!",calledby="nocoInputCheck")
      END IF


!---> make sure starcoeff is switched off
!      IF (starcoeff) THEN
!         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
!         WRITE (oUnit,*) 'does not support starcoefficients output.'
!     CALL juDFT_error("starcoefficients output (for STM) cannot be !!!"
!     generated
!      ENDIF

!---> make sure coretails are switched off
      IF (input%ctail) THEN
         WRITE (oUnit,*) 'This non-collinear version of the flapw program'
         WRITE (oUnit,*) 'cannot be used with the coretail option!! '
         CALL juDFT_error("Coretail option cannot be used!!!",hint="Set /calculationSetup/coreElectrons/@ctail to F in the FLEUR input file.",calledby="nocoInputCheck")
      END IF

!---> make sure that moments are not relaxed and constrained
      l_relax_any = any(noco%l_alignMt.and.noco%l_constrained)
      IF (l_relax_any) THEN
         WRITE (oUnit,*)'The relaxation of the moment is switched on for at'
         WRITE (oUnit,*)'least one atom. At the same time the constrained'
         WRITE (oUnit,*)'moment option has been switched on!!!'
        CALL juDFT_error("You can not constrain and relax a magnetic moment simultaniously")
      ENDIF
!---> make sure that perp. component of mag. is calculated if needed
      IF ( (l_relax_any .or. any(noco%l_constrained)) .and. (.not. noco%l_mperp) ) THEN
         WRITE (oUnit,*)'The relaxation of the moment is switched on for at'
         WRITE (oUnit,*)'least one atom or the constrained moment option is'
         WRITE (oUnit,*)'switched on. In either case, you need to set'
         WRITE (oUnit,*)'l_mperp=T !!'
         CALL juDFT_error("Stop: Set l_mperp = T to relax or constrain the moments!!",calledby ="nocoInputCheck")
      ENDIF
!---> make sure l_constr is switched off in the case of spin spirals
      IF (any(noco%l_constrained) .and. noco%l_ss) THEN
         WRITE (oUnit,*)'The constraint moment option is not implemeted'
         WRITE (oUnit,*)'for spin spirals.'
         CALL juDFT_error("Stop: constraint not implemented for spin spirals!!",calledby ="nocoInputCheck")
      ENDIF

      IF((any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau)) &
         .AND.atoms%n_hia+atoms%n_u>0.AND.sym%nop.NE.1) THEN
         CALL juDFT_warn("LDA+U and FullyFullyNoco with symmetries is not correctly implemented at the moment",calledby="nocoInputCheck")
      ENDIF

    IF(any(noco%l_unrestrictMT).AND.sym%nop.NE.1) THEN
       CALL juDFT_warn("FullyFullyNoco with symmetries might not deliver the desired results for you. This probably would require the implementation of magnetic space groups. ",calledby="nocoInputCheck")
    END IF

    IF(any(noco%l_unrestrictMT).AND.atoms%n_hia+atoms%n_u>0.AND.(.NOT.any(noco%l_alignMT))) THEN
         CALL juDFT_warn("LDA+U and FullyFullyNoco should only be used together with the l_RelaxSQA=T setting to achieve reasonable results.",calledby="nocoInputCheck")
      ENDIF

          !Warning on strange choice of switches
    IF (any(noco%l_unrestrictMT).AND..NOT.noco%l_mperp) THEN
    	CALL juDFT_error("l_mperp='F' and l_mtNocoPot='T' makes no sense.",calledby='nocoInputCheck')
    END IF

    if (noco%l_ss) then
       CALL ss_sym(sym%nop,sym%mrot,noco%qss_inp,error)
       IF (ANY(error)) CALL judft_warn("Symmetry incompatible with Spin-Spiral")
       IF (ANY(noco%qss_inp(:).NE.0.0)) THEN
          IF(ALL(noco%beta_inp(:).EQ.0.0)) CALL juDFT_warn("No spin-spiral cone has a finite opening angle. Is this wanted?", hint='This is the beta angle.')
       END IF
    endif

    IF (any(noco%l_spinoffd_ldau).AND..NOT.noco%l_mperp) THEN
      CALL juDFT_error("l_spinoffd='T' for ldaU and l_mperp='F' makes no sense.",calledby='nocoInputCheck')
    END IF

    IF (any(noco%l_spinoffd_ldau).AND..NOT.noco%l_noco) THEN
      CALL juDFT_error("l_spinoffd='T' for ldaU and l_noco='F' makes no sense.",calledby='nocoInputCheck')
    END IF

    IF (noco%l_ss .and. noco%l_soc) THEN
      CALL juDFT_error("You use l_soc='T' and l_ss='T'.",hint="In a spin-spiral calculation SOC cannot be used. These are incompatible features. Please see the documentation for details.",calledby='nocoInputCheck')
    END IF

   END SUBROUTINE nocoInputCheck

END MODULE m_nocoInputCheck
