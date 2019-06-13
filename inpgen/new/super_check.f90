MODULE m_supercheck
!*********************************************************************
!     checks whether this is a supercell and determines the
!     translations that take the crystal in to itself
!*********************************************************************
CONTAINS
  SUBROUTINE super_check(nat,pos,ity,ntypm, ns,trs)
    
    IMPLICIT NONE
    
    !==> Arguments
    INTEGER, INTENT (IN)  :: nat,ntypm,ity(nat)
    REAL,    INTENT (IN)  :: pos(3,nat)
    INTEGER, INTENT (OUT) :: ns       ! size of supercell
    REAL,    INTENT (OUT) :: trs(3,*) ! translations
    !==> Locals
    INTEGER i,j,n,ntysmin
    INTEGER ntys(ntypm)
    REAL    tr(3)
    LOGICAL l_f
    
    REAL, PARAMETER :: eps=1.e-7
    
    
    ns = 1
    trs(1:3,ns) = (/ 0.000 , 0.000 , 0.000 /)
    
    !---> check if possible to be a supercell
    ntys(1:ntypm) = 0
    DO n = 1, nat
       ntys( ity(n) ) = ntys( ity(n) ) + 1
    ENDDO
    ntysmin = MINVAL( ntys )
    
    IF ( ntysmin == 1 ) RETURN  ! not a supercell if one atom different
    
    DO i = ntysmin,2,-1                 ! check ratios of atoms
       IF ( MOD(ntysmin,i).NE.0 ) CYCLE ! only factors of ntysmin allowed
       l_f = .TRUE.
       DO n = 1, ntypm
          IF ( MOD( ntys(n), i ) .NE. 0 ) THEN
             l_f = .FALSE.
             EXIT
          ENDIF
       ENDDO
       IF (l_f) THEN
          ns = i    ! possible value
          EXIT
       ENDIF
    ENDDO
    IF (ns == 1) RETURN
    
    !---> based on number of atoms (and type), possibly a supercell;
    !---> now need to check by doing translations
    
    ns = 1
    !---> get possible shifts
    DO j = 1, nat-1
       shift_i:  DO i = j+1, nat  ! -ve shifts will come through naturally
          
          tr(:) = pos(:,i)-pos(:,j) - ANINT( pos(:,i)-pos(:,j) - eps )
          
          !--->       check if already done
          DO n=1,ns
             IF ( ALL( ABS( tr(:)-trs(:,n) ) < eps ) ) CYCLE shift_i
          ENDDO
          
          IF ( l_shiftm(tr,pos,nat) ) THEN
             ns = ns + 1
             trs(:,ns) = tr(:)
          ENDIF
          
       ENDDO shift_i
    ENDDO
    
    IF ( ns > 1 ) THEN
       WRITE(6,'(/," The system appears to be a supercell"," containing",i4," primitive cells:")') ns
       DO n = 1, ns
          WRITE (6,'(i8,3f12.6)') n,trs(1:3,n)
       ENDDO
    ENDIF
    
  CONTAINS ! internal function
    
    LOGICAL FUNCTION l_shiftm(tr,pos,nat)
      !********************************************************************
      !     determines whether the vector tr is a translation of the
      !     crystal (non-primitive for supercell)
      !********************************************************************
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nat
      REAL,    INTENT(IN) :: tr(3)
      REAL,    INTENT(IN) :: pos(3,nat)
      
      REAL    rp(3)
      INTEGER i,j,in
      
      l_shiftm = .FALSE.
      
      DO i = 1, nat
         !--->    rotated and shifted atom, reduced to (-1/2,1/2]
         rp(:) = pos(:,i) + tr(:) - ANINT( pos(:,i) + tr(:) - eps )
         !--->    find which atom, if any, this matches
         in = 0
         DO j = 1, nat
            IF ( ity(i).NE.ity(j) ) CYCLE
            !            if( all( abs(pos(:,j)-rp(:) ) < eps ) ) then
            ! causes problem with intel compiler ifc Version 5.0.1 (gs2001-11-07)
            IF ( ABS(pos(1,j)-rp(1) ) < eps  .AND. ABS(pos(2,j)-rp(2) ) < eps  .AND.&
                 ABS(pos(3,j)-rp(3) ) < eps  ) THEN
               in = j
               EXIT
            ENDIF
         ENDDO
         IF (in == 0 ) RETURN
      ENDDO
      
      l_shiftm = .TRUE.  !  only if everything matches
      
    END FUNCTION l_shiftm
    
  END SUBROUTINE super_check
END MODULE m_supercheck
