!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_slomat
  !***********************************************************************
  ! updates the overlap matrix with the contributions from the local
  ! orbitals.
  !                                                p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE slomat(&
       input,atoms,mpi,lapw,cell,noco,ntyp,na,&
       isp,ud, alo1,blo1,clo1,fj,gj,&
       iintsp,jintsp,chi,smat)
    !***********************************************************************
    ! locol stores the number of columns already processed; on parallel
    !       computers this decides, whether the LO-contribution is
    !       done on this node                                          gb00
    !
    ! function legpol() at end of module
    !***********************************************************************
    USE m_constants,ONLY: fpi_const
    USE m_types
    USE m_apws
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: mpi
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)      :: na,ntyp 
    INTEGER, INTENT (IN)      :: iintsp,jintsp
    COMPLEX, INTENT (IN)      :: chi
    INTEGER, INTENT(IN)       :: isp
    !     ..
    !     .. Array Arguments ..
    REAL,   INTENT (IN)       :: alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    REAL,    INTENT (IN) :: fj(:,0:,:),gj(:,0:,:)
    TYPE(t_usdus),INTENT(IN)  :: ud
    TYPE(t_lapwmat),INTENT(INOUT) :: smat
    
    !     ..
    !     .. Local Scalars ..
    REAL con,dotp,fact1,fact2,fact3,fl2p1
    INTEGER invsfct,k ,l,lo,lop,lp,nkvec,nkvecp,kp,i
    INTEGER locol,lorow
    !     ..
    !     ..

    COMPLEX,   ALLOCATABLE  :: cph(:,:)
    ALLOCATE(cph(MAXVAL(lapw%nv),2))
    DO i=MIN(iintsp,jintsp),MAX(iintsp,jintsp)
       CALL lapw_phase_factors(lapw,i,atoms%taual(:,na),noco%qss,cph(:,i))
    ENDDO

    IF ((atoms%invsat(na).EQ.0) .OR. (atoms%invsat(na).EQ.1)) THEN
       !--->    if this atom is the first of two atoms related by inversion,
       !--->    the contributions to the overlap matrix of both atoms are added
       !--->    at once. where it is made use of the fact, that the sum of
       !--->    these contributions is twice the real part of the contribution
       !--->    of each atom. note, that in this case there are twice as many
       !--->    (2*(2*l+1)) k-vectors (compare abccoflo and comments there).
       IF (atoms%invsat(na).EQ.0) invsfct = 1
       IF (atoms%invsat(na).EQ.1) invsfct = 2
   
       con = fpi_const/SQRT(cell%omtil)* ((atoms%rmt(ntyp))**2)/2.0

       DO lo = 1,atoms%nlo(ntyp) !loop over all LOs for this atom
          l = atoms%llo(lo,ntyp)
          fl2p1 = (2*l+1)/fpi_const
          fact1 = (con**2)* fl2p1 * (&
               alo1(lo)* (  alo1(lo) + &
               2*clo1(lo) * ud%uulon(lo,ntyp,isp) ) +&
               blo1(lo)* (  blo1(lo) * ud%ddn(l, ntyp,isp) +&
               2*clo1(lo) * ud%dulon(lo,ntyp,isp) ) +&
               clo1(lo)*    clo1(lo) )
          DO nkvec = 1,invsfct* (2*l+1) !Each LO can have several functions
             !+t3e
             locol = lapw%nv(iintsp)+lapw%index_lo(lo,ntyp)+nkvec !this is the column of the matrix
             IF (MOD(locol-1,mpi%n_size).EQ.mpi%n_rank) THEN
                locol=(locol-1)/mpi%n_size+1 !this is the column in local storage
                !-t3e
                k = lapw%kvec(nkvec,lo,ntyp)
                !--->          calculate the overlap matrix elements with the regular
                !--->          flapw basis-functions
                DO kp = 1,lapw%nv(jintsp)
                   fact2 = con * fl2p1 * (&
                        fj(kp,l,jintsp)* ( alo1(lo) + &
                        clo1(lo)*ud%uulon(lo,ntyp,isp))+&
                        gj(kp,l,jintsp)* ( blo1(lo) * ud%ddn(l,ntyp,isp)+&
                        clo1(lo)*ud%dulon(lo,ntyp,isp)))
                   dotp = dot_PRODUCT(lapw%gk(:,k,iintsp),lapw%gk(:,kp,jintsp))
                   IF (smat%l_real) THEN
                      smat%data_r(kp,locol) = smat%data_r(kp,locol) + chi*invsfct*fact2 * legpol(l,dotp) *&
                           cph(k,iintsp)*CONJG(cph(kp,jintsp))
                   ELSE
                      smat%data_c(kp,locol) = smat%data_c(kp,locol) + chi*invsfct*fact2 * legpol(l,dotp) *&
                           cph(k,iintsp)*CONJG(cph(kp,jintsp))
                   ENDIF
                END DO
                !--->          calculate the overlap matrix elements with other local
                !--->          orbitals at the same atom, if they have the same l
                DO lop = 1, (lo-1)
                   lp = atoms%llo(lop,ntyp)
                   IF (l.EQ.lp) THEN
                      fact3 = con**2 * fl2p1 * (&
                           alo1(lop)*(alo1(lo) + &
                           clo1(lo)*ud%uulon(lo,ntyp,isp))+&
                           blo1(lop)*(blo1(lo)*ud%ddn(l,ntyp,isp) +&
                           clo1(lo)*ud%dulon(lo,ntyp,isp))+&
                           clo1(lop)*(alo1(lo)*ud%uulon(lop,ntyp,isp)+&
                           blo1(lo)*ud%dulon(lop,ntyp,isp)+&
                           clo1(lo)*ud%uloulopn(lop,lo,ntyp,isp)))
                      DO nkvecp = 1,invsfct* (2*lp+1)
                         kp = lapw%kvec(nkvecp,lop,ntyp)
                         lorow=lapw%nv(jintsp)+lapw%index_lo(lop,ntyp)+nkvecp
                         dotp = dot_PRODUCT(lapw%gk(:,k,iintsp),lapw%gk(:,kp,jintsp))
                         IF (smat%l_real) THEN
                            smat%data_r(lorow,locol) =smat%data_r(lorow,locol)+chi*invsfct*fact3*legpol(l,dotp)* &
                                 cph(k,iintsp)*conjg(cph(kp,jintsp))
                         ELSE
                            smat%data_c(lorow,locol) =smat%data_c(lorow,locol)+chi*invsfct*fact3*legpol(l,dotp)*&
                                 cph(k,iintsp)*CONJG(cph(kp,jintsp)) 
                         ENDIF
                      END DO
                   ELSE
                   END IF
                END DO
                !--->          calculate the overlap matrix elements of one local
                !--->          orbital with itself
                DO nkvecp = 1,nkvec
                   kp = lapw%kvec(nkvecp,lo,ntyp)
                   lorow=lapw%nv(jintsp)+lapw%index_lo(lo,ntyp)+nkvecp
                   dotp = dot_PRODUCT(lapw%gk(:,k,iintsp),lapw%gk(:,kp,jintsp))
                   IF (smat%l_real) THEN
                      smat%data_r(lorow,locol) = smat%data_r(lorow,locol) + chi*invsfct*fact1*legpol(l,dotp) *&
                           cph(k,iintsp)*CONJG(cph(kp,jintsp))
                   ELSE
                      smat%data_c(lorow,locol) = smat%data_c(lorow,locol) + chi*invsfct*fact1*legpol(l,dotp)*&
                           cph(k,iintsp)*CONJG(cph(kp,jintsp))
                   ENDIF
                END DO
             ENDIF ! mod(locol-1,n_size) = nrank 
             !-t3e
          END DO
       END DO
    END IF
  END SUBROUTINE slomat
  !===========================================================================
  PURE REAL FUNCTION legpol(l,arg)
    !
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    REAL,INTENT(IN)   :: arg
    INTEGER,INTENT(IN):: l
    !     ..
    !     .. Local Scalars ..
    INTEGER lp
    !     ..
    !     .. Local Arrays ..
    REAL plegend(0:l)
    !     ..
    plegend(0) = 1.0
    IF (l.GE.1) THEN
       plegend(1) = arg
       DO lp = 1,l - 1
          plegend(lp+1) = (lp+lp+1)*arg*plegend(lp)/ (lp+1) -lp*plegend(lp-1)/ (lp+1)
       END DO
    END IF
    legpol = plegend(l)
  END FUNCTION legpol
  !===========================================================================
END MODULE m_slomat
