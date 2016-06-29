!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_stmix
  !
  !      straight mixing,     r.pentcheva, iff, 1996
  !
  !     sm   : input charge density of iteration m
  !     sm1  : input charge density of iteration m+1
  !     fsm  : output minus input charge densityof iteration m
  !
CONTAINS
  SUBROUTINE stmix(&
       &                 atoms,input,noco,&
       &                 nmap,nmaph,fsm,sm)

    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nmaph,nmap  
    !     ..
    !     .. Array Arguments ..
    REAL fsm(:),sm(:)
    !     ..
    !     .. Local Scalars ..
    INTEGER imap
    REAL,PARAMETER:: tol_6=1.0e-6
    !     ..
    !
    WRITE (16,FMT='(a)') 'STRAIGHT MIXING'
    IF (input%jspins.EQ.1) WRITE (16,FMT='(a,2f10.5)')&
         &    'charge density mixing parameter:',input%alpha
    IF (input%jspins.EQ.2) WRITE (16,FMT='(a,2f10.5)')&
         &    'spin density mixing parameter:',input%alpha*input%spinf
    IF ( ABS(input%spinf-1.0e0).LE.tol_6 .OR. input%jspins.EQ.1 ) THEN
       !     --> perform simple mixing 
       !
       !        sm1 = sm + alpha * F(sm)

       sm(:nmap) = sm(:nmap) + input%alpha*fsm(:nmap)
       RETURN
    ELSE
       !     -->perform simple mixing with the mixing parameters
       !        for charge and spin
       !
       !       sm1+/_ = (sm+/_) + alpha* F(sm)
       !                +/-0.5alpha(spinf-1)( F(sm+) + F(sm-) )

       DO imap = 1,nmaph
          sm(imap) = sm(imap) + input%alpha*fsm(imap) &
               &            + input%alpha/2.0*(input%spinf-1.0)*(fsm(imap) - fsm(imap+nmaph))
       ENDDO

       DO imap = nmaph+1,2*nmaph
          sm(imap) = sm(imap) + input%alpha*fsm(imap) &
               &            + input%alpha/2.0*(input%spinf-1.0)*(fsm(imap) - fsm(imap-nmaph))
       ENDDO
       IF (noco%l_noco) THEN
          DO imap = 2*nmaph+1, nmap - 98*input%jspins*atoms%n_u
             sm(imap) = sm(imap) + input%alpha*input%spinf*fsm(imap) 
          ENDDO
       ENDIF
       IF ( atoms%n_u > 0 )  THEN
          DO imap = nmap - 98*input%jspins*atoms%n_u + 1, nmap - 98*atoms%n_u 
             sm(imap) = sm(imap) + input%alpha*fsm(imap) &
                  &            + input%alpha/2.0*(input%spinf-1.0)*(fsm(imap) - fsm(imap+98*atoms%n_u))
          ENDDO
          DO imap = nmap - 98*atoms%n_u + 1, nmap
             sm(imap) = sm(imap) + input%alpha*fsm(imap) &
                  &            + input%alpha/2.0*(input%spinf-1.0)*(fsm(imap) - fsm(imap-98*atoms%n_u))
          ENDDO
       ENDIF
    END IF

  END SUBROUTINE stmix
END MODULE m_stmix
