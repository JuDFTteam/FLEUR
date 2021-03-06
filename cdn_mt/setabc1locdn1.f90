!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_setabc1locdn1
  !***********************************************************************
  ! calculates the (lower case) a, b and c coefficients for the local
  ! orbitals. The radial function of the local orbital is a linear 
  ! combination of the apw radial function and its derivative and the
  ! extra radial funtion (a*u + b*udot + c*ulo). This function is zero
  ! and has zero derivative at the muffin tin boundary.
  ! In addition the the total number of basisfuntions (apw + lo) nbasf and
  ! the number of the first basisfunction of each local orbital nbasf0 is
  ! determined.
  ! Philipp Kurz 99/04
  !***********************************************************************
CONTAINS
  SUBROUTINE setabc1locdn1(jsp,atoms,lapw,sym,usdus,&
        enough,nkvec,kvec,nbasf0, alo1,blo1,clo1)
    !
    !*************** ABBREVIATIONS *****************************************
    ! nbasf   : total number of basisfunctions (apw + lo)
    ! nbasf0  : number of the first basisfunction of each local orbital
    ! nkvec   : stores the number of G-vectors that have been found and
    !           accepted during the construction of the local orbitals.
    !***********************************************************************
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_usdus),INTENT(IN)   :: usdus
    TYPE(t_lapw),INTENT(IN)    :: lapw
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    INTEGER,INTENT(IN)         :: jsp
    !     .. Array Arguments ..
    INTEGER, INTENT (OUT) :: nbasf0(atoms%nlod,atoms%nat),nkvec(atoms%nlod,atoms%nat)
    INTEGER, INTENT (OUT) :: kvec(2*(2*atoms%llod+1),atoms%nlod,atoms%nat  )
    REAL,    INTENT (OUT) :: alo1(atoms%nlod,atoms%ntype),blo1(atoms%nlod,atoms%ntype)
    REAL,    INTENT (OUT) :: clo1(atoms%nlod,atoms%ntype)
    LOGICAL, INTENT (OUT) :: enough(atoms%nat)
    !     ..
    !     .. Local Scalars ..
    REAL ka,kb,ws
    INTEGER i,l,lo ,natom,nbasf,nn,ntyp,lm,m
    LOGICAL apw_at
    !     ..
    enough(:) = .true.
    DO ntyp = 1,atoms%ntype
       !     ..
       ! look, whether 'ntyp' is a APW atom; then set apw_at=.true.
       !
       apw_at = .false.
       DO lo = 1,atoms%nlo(ntyp)
          IF (atoms%l_dulo(lo,ntyp)) apw_at = .true.
       ENDDO

       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          IF (apw_at) THEN
             IF (atoms%l_dulo(lo,ntyp)) THEN
                ! udot lo
                ka = sqrt( 1+(usdus%us(l,ntyp,jsp)/usdus%uds(l,ntyp,jsp))**2 * usdus%ddn(l,ntyp,jsp))
                alo1(lo,ntyp)=1.00 / ka
                blo1(lo,ntyp)=-usdus%us(l,ntyp,jsp)/ (usdus%uds(l,ntyp,jsp) * ka)
                clo1(lo,ntyp)=0.00
             ELSE
                ! u2 lo
                alo1(lo,ntyp)=1.00
                blo1(lo,ntyp)=0.00
                clo1(lo,ntyp)=-usdus%us(l,ntyp,jsp)/usdus%ulos(lo,ntyp,jsp)
             ENDIF
          ELSE
             ws = usdus%uds(l,ntyp,jsp)*usdus%dus(l,ntyp,jsp) - usdus%us(l,ntyp,jsp)*usdus%duds(l,ntyp,jsp)
             ka = 1.0/ws* (usdus%duds(l,ntyp,jsp)*usdus%ulos(lo,ntyp,jsp)- usdus%uds(l,ntyp,jsp)*usdus%dulos(lo,ntyp,jsp))
             kb = 1.0/ws* (usdus%us(l,ntyp,jsp)*usdus%dulos(lo,ntyp,jsp)- usdus%dus(l,ntyp,jsp)*usdus%ulos(lo,ntyp,jsp))
             clo1(lo,ntyp) = 1.0/sqrt(ka**2+ (kb**2)*usdus%ddn(l,ntyp,jsp)+1.0+&
                  2.0*ka*usdus%uulon(lo,ntyp,jsp)+2.0*kb*usdus%dulon(lo,ntyp,jsp))
             alo1(lo,ntyp) = ka*clo1(lo,ntyp)
             blo1(lo,ntyp) = kb*clo1(lo,ntyp)
          ENDIF
       END DO
    END DO
    !---> set up enough, nbasf0 and initialize nkvec
    natom = 0
    nbasf = lapw%nv(jsp)
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          DO lo = 1,atoms%nlo(ntyp)
             enough(natom) = .false.
             nkvec(lo,natom) = 0
             l = atoms%llo(lo,ntyp)
             IF (sym%invsat(natom).EQ.0) THEN
                nbasf0(lo,natom) = nbasf
                nbasf = nbasf + 2*l + 1
             END IF
             IF (sym%invsat(natom).EQ.1) THEN
                nbasf0(lo,natom) = nbasf
                nbasf0(lo,sym%invsatnr(natom)) = nbasf
                nbasf = nbasf + 2* (2*l+1)
             END IF
          END DO
       END DO
    END DO


    !      write (*,*) 'in setabc1locdn: nmat = ',nmat,' nbasf = ',nbasf
    !      write (*,*) 'array nbasf0 :'
    !      do natom = 1,natd
    !         write (*,fmt='(15i4)') (nbasf0(lo,natom),lo=1,nlod)
    !      enddo
    !      write (*,*)
    IF ((lapw%nmat).NE.nbasf) THEN
       write (*,*) 'in setabc1locdn: lapw%nmat = ',lapw%nmat,' nbasf = ',nbasf
       STOP 'setabc1locdn: number of bas.-fcn.'
    ENDIF
    !
    !--> sort the k-vectors used for the LO's according to atom & lo:
    !
    natom = 0
    lm = 0
    DO ntyp = 1, atoms%ntype
       DO nn = 1, atoms%neq(ntyp)
          natom = natom + 1
          IF ((sym%invsat(natom).EQ.0) .OR. (sym%invsat(natom).EQ.1)) THEN
             DO lo = 1,atoms%nlo(ntyp)
                m = ( sym%invsat(natom) +1 ) * ( 2 * atoms%llo(lo,ntyp) + 1 )
                DO l = 1, m
                   lm = lm + 1
                   kvec(l,lo,natom) =  lapw%kvec(l,lo,natom)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE setabc1locdn1
END MODULE m_setabc1locdn1
