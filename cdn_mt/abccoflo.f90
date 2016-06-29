!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abccoflo
  USE m_juDFT
  !*********************************************************************
  ! Calculates the (upper case) A, B and C coefficients for the local
  ! orbitals.
  ! Philipp Kurz 99/04
  !*********************************************************************
CONTAINS
  SUBROUTINE abccoflo(atoms, con1,rph,cph,ylm,ntyp,na,k,nv, l_lo1,alo1,blo1,&
       clo1, nkvec, enough,alo,blo,clo,kvec)
    !
    !*************** ABBREVIATIONS ***************************************
    ! kvec    : stores the number of the G-vectors, that have been used to
    !           construct the local orbitals
    ! nkvec   : stores the number of G-vectors that have been found and
    !           accepted during the construction of the local orbitals.
    ! enough  : enough is set to .true. when enough G-vectors have been
    !           accepted.
    ! linindq : if the norm of that part of a local orbital (contructed 
    !           with a trial G-vector) that is orthogonal to the previous
    !           ones is larger than linindq, then this G-vector is 
    !           accepted.
    !*********************************************************************
    !
    USE m_constants
    USE m_types
    IMPLICIT NONE

    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     .. 
    !     .. Scalar Arguments ..
    REAL,    INTENT (IN) :: con1,cph ,rph
    INTEGER, INTENT (IN) :: k,na,ntyp,nv
    LOGICAL, INTENT (IN) :: l_lo1
    LOGICAL, INTENT (OUT):: enough
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN)::  kvec(2* (2*atoms%llod+1),atoms%nlod )
    REAL,    INTENT (IN) :: alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    COMPLEX, INTENT (IN) :: ylm( (atoms%lmaxd+1)**2 )
    COMPLEX, INTENT (OUT):: alo(-atoms%llod:atoms%llod,2* (2*atoms%llod+1),atoms%nlod)  
    COMPLEX, INTENT (OUT):: blo(-atoms%llod:atoms%llod,2* (2*atoms%llod+1),atoms%nlod)  
    COMPLEX, INTENT (OUT):: clo(-atoms%llod:atoms%llod,2* (2*atoms%llod+1),atoms%nlod)  
    INTEGER,INTENT (INOUT):: nkvec(atoms%nlod)
    !     ..
    !     .. Local Scalars ..
    COMPLEX term1
    REAL,PARAMETER:: linindq=1.e-4
    INTEGER l,lo ,mind,ll1,lm,m
    LOGICAL linind
    !     ..
    !
    !---> the whole program is in hartree units, therefore 1/wronskian is
    !---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
    !---> and c coefficients, is included in the t-matrices. thus, it does
    !---> not show up in the formula above.
    !
    !-abccoflo1
    IF ( l_lo1) THEN
       DO lo = 1,atoms%nlo(ntyp)
          IF ( (nkvec(lo).EQ.0).AND.(atoms%llo(lo,ntyp).EQ.0) ) THEN
             enough = .FALSE.
             nkvec(lo) = 1
             m = 0
             clo(m,nkvec(lo),lo) = con1* ((atoms%rmt(ntyp)**2)/2) / SQRT(fpi_const)
             alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*alo1(lo)
             blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*blo1(lo)
             clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*clo1(lo)
             IF (kvec(nkvec(lo),lo)/=k)  CALL juDFT_error("abccoflo:1"&
                  &           ,calledby ="abccoflo")

          ENDIF
       ENDDO
    ELSE
       enough = .TRUE.
       term1 = con1* ((atoms%rmt(ntyp)**2)/2)*CMPLX(rph,cph)
       DO lo = 1,atoms%nlo(ntyp)
          IF (atoms%invsat(na).EQ.0) THEN
             IF ((nkvec(lo)).LT. (2*atoms%llo(lo,ntyp)+1)) THEN
                enough = .FALSE.
                nkvec(lo) = nkvec(lo) + 1
                l = atoms%llo(lo,ntyp)
                ll1 = l*(l+1) + 1
                DO m = -l,l
                   lm = ll1 + m
                   clo(m,nkvec(lo),lo) = term1*ylm(lm)
                END DO
                IF ( kvec(nkvec(lo),lo) == k ) THEN
                   DO m = -l,l
                      alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*alo1(lo)
                      blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*blo1(lo)
                      clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*clo1(lo)
                   END DO
                   !                  WRITE(6,9000) nkvec(lo),k,lo,na,
                   !     +                          (clo(m,nkvec(lo),lo),m=-l,l)
                   ! 9000             format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
                ELSE
                   nkvec(lo) = nkvec(lo) - 1
                ENDIF
             ENDIF
          ELSE
             IF ((atoms%invsat(na).EQ.1) .OR. (atoms%invsat(na).EQ.2)) THEN
                !           only invsat=1 is needed invsat=2 for testing
                IF ((nkvec(lo)).LT. (2* (2*atoms%llo(lo,ntyp)+1))) THEN
                   enough = .FALSE.
                   nkvec(lo) = nkvec(lo) + 1
                   l = atoms%llo(lo,ntyp)
                   ll1 = l*(l+1) + 1
                   DO m = -l,l
                      lm = ll1 + m
                      clo(m,nkvec(lo),lo) = term1*ylm(lm)
                   END DO
                   IF ( kvec(nkvec(lo),lo) == k ) THEN
                      DO m = -l,l
                         !                            if(l.eq.1) then
                         !               WRITE(*,*)'k=',k,' clotmp=',clo(m,nkvec(lo),lo)
                         !               WRITE(*,*)'clo1=',clo1(lo),' term1=',term1
                         !                            endif
                         alo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*alo1(lo)
                         blo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*blo1(lo)
                         clo(m,nkvec(lo),lo) = clo(m,nkvec(lo),lo)*clo1(lo)
                         !                        kvec(nkvec(lo),lo) = k
                      END DO
                   ELSE
                      nkvec(lo) = nkvec(lo) - 1
                   END IF
                END IF
             END IF
          END IF
       END DO
       IF ((k.EQ.nv) .AND. (.NOT.enough)) THEN
          WRITE (6,FMT=*)&
               &     'abccoflo did not find enough linearly independent'
          WRITE (6,FMT=*)&
               &     'clo coefficient-vectors. the linear independence'
          WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
          WRITE (6,FMT=*) 'this value might be to large.'
          CALL juDFT_error&
               &        ("abccoflo: did not find enough lin. ind. clo-vectors"&
               &        ,calledby ="abccoflo")
       END IF
    ENDIF  ! abccoflo1

  END SUBROUTINE abccoflo
END MODULE m_abccoflo
