!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_inpnoco
      use m_juDFT
!**********************************************************************
!     This subroutine reads the  Euler angles of the  magnetic field 
!     directions and other noncollinear parameters from the file nocoinp
!
!                                                 Philipp Kurz 98/01/28
!******** ABBREVIATIONS ***********************************************
!     alpha,beta:Euler angles of the local magnetic field direction of
!                each atom(-type). 
!**********************************************************************
      CONTAINS
      SUBROUTINE inpnoco(atoms,input,vacuum,jij,noco)

      USE m_constants, ONLY : tpi_const
      USE m_rwnoco
      USE m_types
      IMPLICIT NONE
      TYPE(t_atoms),INTENT(INOUT) ::atoms
      TYPE(t_input),INTENT(INOUT) ::input
      TYPE(t_vacuum),INTENT(IN)   ::vacuum
      TYPE(t_Jij),INTENT(INOUT)   ::Jij
      TYPE(t_noco),INTENT(INOUT)  ::noco

!     ..
!     .. Local Scalars ..
      INTEGER itype,iatom
      LOGICAL l_relax_any

!---> make sure that the complex program has been compiled
#ifdef CPP_INVERSION
         WRITE (6,*) 'Non-collinear calculations can only be done with'
         WRITE (6,*) 'the complex program. Please compile without the'
         WRITE (6,*) 'option "CPP_INVERSION".'
         CALL juDFT_error&
     &        ("for l_noco = T: recompile without CPP_INVERSION!"&
     &        ,calledby ="inpnoco")
#endif

!---> make sure Wu-diagonalization is switched off
      IF (input%isec1 .LE. input%itmax) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Wu-diagonalization!!'
         WRITE (6,*) 'itmax = ',input%itmax,'isec1 = ',input%isec1
         CALL juDFT_error("Wu-diagonalization cannot be used!!!",&
     &        calledby="inpnoco")
      ENDIF

!---> make sure second variation is switched off
      IF (input%secvar) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the second variation!!'
         CALL juDFT_error("Second variation cannot be used!!!" ,calledby&
     &        ="inpnoco")
      ENDIF

!---> make sure histogram method is used
      IF (input%gauss) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the Gaussian smearing for '
         WRITE (6,*) 'the Brillouin zone integration!!'
         WRITE (6,*) 'Please use the histogram method.'
         CALL juDFT_error&
     &       ("Only histogram Brillouin zone integration can be used!!!"&
     &    ,calledby ="inpnoco")
      ENDIF

!---> make sure force is switched off
      IF (input%l_f) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support force calculations.'
         CALL juDFT_error("force calculations not supported!!!",calledby&
     &        ="inpnoco")
      ENDIF

!---> make sure nstm equals zero
      IF (vacuum%nstm.NE.0) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'does not support STM calculations(nstm .NE. 0).'
         CALL juDFT_error("nstm /= 0 not supported!",calledby ="inpnoco"&
     &        )
      ENDIF

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
         CALL juDFT_error("Coretail option cannot be used!!!",calledby&
     &        ="inpnoco")
      ENDIF

!---> make sure score is false
      IF (input%score) THEN
         WRITE (6,*) 'This non-collinear version of the flapw program'
         WRITE (6,*) 'cannot be used with the score option!! '
         CALL juDFT_error("score must be false!!!",calledby ="inpnoco")
      ENDIF

      OPEN (24,file='nocoinp',form='formatted',status='old')

      WRITE (6,*)'This is a non-collinear calculation. The magnetic'
      WRITE (6,*)'moments of the atoms have a fixed direction.'
      WRITE (6,*)'The Euler-angles alpha and beta of this direction'
      WRITE (6,*)'are equal to the polar angles of the magnetic'
      WRITE (6,*)'moment vector phi and theta respectively.'
      WRITE (6,*)

      CALL rw_noco_read(&
     &             atoms,jij,noco,input)

!---> make sure that moments are not relaxed and constrained
      l_relax_any = .false.
      DO itype = 1,atoms%ntype
         l_relax_any = l_relax_any.OR.noco%l_relax(itype)
      ENDDO
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
        CALL juDFT_error&
     &      ("Stop: Set l_mperp = T to relax or constrain the moments!!"&
     &,calledby ="inpnoco")
      ENDIF
!---> make sure l_constr is switched off in the case of spin spirals
      IF (noco%l_constr .and. noco%l_ss) THEN
        WRITE (6,*)'The constraint moment option is not implemeted'
        WRITE (6,*)'for spin spirals.'
        CALL juDFT_error&
     &       ("Stop: constraint not implemented for spin spirals!!"&
     &       ,calledby ="inpnoco")
      ENDIF
         
      IF (.not.jij%l_j.and.noco%l_ss) THEN
!
!--->    the angle beta is relative to the spiral in a spin-spiral
!--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
!--->    that means that the moments are "in line" with the spin-spiral
!--->    (beta = qss * taual). note: this means that only atoms within
!--->    a plane perpendicular to qss can be equivalent!
         iatom = 1
         DO itype = 1,atoms%ntype
            noco%phi = tpi_const*dot_product(noco%qss,atoms%taual(:,iatom))
            noco%alph(itype) = noco%alph(itype) + noco%phi
            iatom = iatom + atoms%neq(itype)
         ENDDO
      ENDIF

      WRITE (6,*)
      CLOSE (24)

      END SUBROUTINE inpnoco
      END MODULE m_inpnoco
