!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_nocoInputCheck

   CONTAINS

   SUBROUTINE nocoInputCheck(atoms,input,vacuum,jij,noco)

      USE m_juDFT
      USE m_types

      IMPLICIT NONE

      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_vacuum), INTENT(IN)    :: vacuum
      TYPE(t_Jij),    INTENT(IN)    :: Jij
      TYPE(t_noco),   INTENT(IN)    :: noco

      INTEGER itype
      LOGICAL l_relax_any



!---> make sure Wu-diagonalization is switched off
      IF (input%isec1 .LE. input%itmax) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Wu-diagonalization!!'
         WRITE (6,*) 'itmax = ',input%itmax,'isec1 = ',input%isec1
         CALL juDFT_error("Wu-diagonalization cannot be used!!!",calledby="nocoInputCheck")
      END IF

!---> make sure second variation is switched off
      IF (input%secvar) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the second variation!!'
         CALL juDFT_error("Second variation cannot be used!!!" ,calledby="nocoInputCheck")
      END IF

!---> make sure histogram method is used
      IF (input%gauss) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Gaussian smearing for '
         WRITE (6,*) 'the Brillouin zone integration!!'
         WRITE (6,*) 'Please use the histogram method.'
         CALL juDFT_error("Only histogram Brillouin zone integration can be used!!!",calledby ="nocoInputCheck")
      END IF

!---> make sure force is switched off
      IF (input%l_f) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support force calculations.'
         CALL juDFT_error("force calculations not supported!!!",calledby="nocoInputCheck")
      END IF

!---> make sure nstm equals zero
      IF (vacuum%nstm.NE.0) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support STM calculations(nstm .NE. 0).'
         CALL juDFT_error("nstm /= 0 not supported!",calledby ="nocoInputCheck")
      END IF

!---> make sure starcoeff is switched off
!      IF (starcoeff) THEN
!         WRITE (6,*) 'This non-collinear version of the flapw program'
!         WRITE (6,*) 'does not support starcoefficients output.'
!     CALL juDFT_error("starcoefficients output (for STM) cannot be !!!"
!     generated
!      ENDIF

!---> make sure coretails are switched off
      IF (input%ctail) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the coretail option!! '
         CALL juDFT_error("Coretail option cannot be used!!!",calledby="nocoInputCheck")
      END IF

!---> make sure score is false
      IF (input%score) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the score option!! '
         CALL juDFT_error("score must be false!!!",calledby ="nocoInputCheck")
      END IF

!---> make sure that moments are not relaxed and constrained
      l_relax_any = .FALSE.
      DO itype = 1,atoms%ntype
         l_relax_any = l_relax_any.OR.noco%l_relax(itype)
      END DO
      IF (l_relax_any.AND.noco%l_constr) THEN
         WRITE (6,*)'The relaxation of the moment is switched on for at'
         WRITE (6,*)'least one atom. At the same time the constrained'
         WRITE (6,*)'moment option has been switched on!!!'
!          CALL juDFT_error("relaxation of moments and constraint are sw
      ENDIF
!---> make sure that perp. component of mag. is calculated if needed
      IF ( (l_relax_any .or. noco%l_constr) .and. (.not. noco%l_mperp) ) THEN
         WRITE (6,*)'The relaxation of the moment is switched on for at'
         WRITE (6,*)'least one atom or the constrained moment option is'
         WRITE (6,*)'switched on. In either case, you need to set'
         WRITE (6,*)'l_mperp=T !!'
         CALL juDFT_error("Stop: Set l_mperp = T to relax or constrain the moments!!",calledby ="nocoInputCheck")
      ENDIF
!---> make sure l_constr is switched off in the case of spin spirals
      IF (noco%l_constr .and. noco%l_ss) THEN
         WRITE (6,*)'The constraint moment option is not implemeted'
         WRITE (6,*)'for spin spirals.'
         CALL juDFT_error("Stop: constraint not implemented for spin spirals!!",calledby ="nocoInputCheck")
      ENDIF

   END SUBROUTINE nocoInputCheck

END MODULE m_nocoInputCheck
