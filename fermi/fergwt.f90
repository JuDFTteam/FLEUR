MODULE  m_fergwt
  USE m_juDFT
  !****************************************************************
  !     determines the fermi energy and weights for the k-space
  !     integration using gaussing-smearing method.
  !                                               c.l.fu
  !*****************************************************************
CONTAINS
  SUBROUTINE fergwt(kpts,input,mpi, ne,eig, ef,w_iks,seigv)

    USE m_constants
    USE m_types
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_kpts),INTENT(IN)      :: kpts
    REAL,        INTENT(INOUT)   :: ef,seigv
    REAL,INTENT(INOUT)           :: w_iks(:,:,:)
    !     ..
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: ne(:,:)    !(kpts%nkpt,dimension%jspd)
    REAL,    INTENT (IN) :: eig(:,:,:) !dimension%neigd,kpts%nkpt,dimension%jspd)
    !     ..
    !     .. Local Scalars ..
    REAL chmom,de,ef0,ef1,elow,en,eps,eup,fac,fact1,s,s0,s1,s2,&
         workf,wt,wtk,zcdiff,zero,seigv1
    INTEGER i,ifl,it,jspin,k,nbnd
    !     ..
    !     .. External Functions ..
    REAL  erf
    !     ..
    !     .. Data statements ..
    DATA zero/0.e0/,eps/1.e-5/,eup/3.0e0/,elow/-3.0e0/
    !     ..
    fact1 = input%delgau/SQRT(pi_const)
    !     ---> determines ef
    ifl = 0
    conv_loop:DO WHILE (.TRUE.)
       DO  it = 1,50
          s = 0.
          DO  jspin = 1,input%jspins
             DO  k = 1,kpts%nkpt
                wtk = kpts%wtkpt(k)
                nbnd = ne(k,jspin)
                DO  i = 1,nbnd
                   en = eig(i,k,jspin)
                   de = (en-ef)/input%delgau
                   wt = 2.0
                   IF (de.GT.eup) wt = 0.0
                   IF (de.GE.elow .AND. de.LE.eup) THEN
                      IF (de.LT.zero) THEN
                         wt = 1. + ERF(-de)
                      ELSE
                         wt = 1. - ERF(de)
                      END IF
                   END IF
                   s = s + wt*wtk
                   w_iks(i,k,jspin) = wt/2.
                ENDDO
             ENDDO
          ENDDO
          s = s/REAL(input%jspins)
          zcdiff = input%zelec - s
          IF (ABS(zcdiff).LT.eps) EXIT conv_loop
          IF (ifl.EQ.0) THEN
             ifl = 1
             ef0 = ef
             ef = ef + 0.003
             s0 = s
          ELSE
             fac = (s0-s)/ (input%zelec-s)
             IF (ABS(fac).LT.1.0e-1) THEN
                ef0 = ef
                s0 = s
                IF (zcdiff.GE.zero) THEN
                   ef = ef + 0.003
                ELSE
                   ef = ef - 0.003
                END IF
             ELSE
                ef1 = ef
                ef = ef + (ef0-ef)/fac
                ef0 = ef1
                s0 = s
             END IF
          END IF
       ENDDO
       eps = 1.25*eps
       IF ( mpi%irank == 0 ) WRITE (6,FMT=8000) eps
8000   FORMAT (10x,'warning: eps has been increased to',e12.5)
    ENDDO conv_loop
    workf = -hartree_to_ev_const*ef
    IF ( mpi%irank == 0 ) THEN
       WRITE (16,FMT=8010) ef,workf,s
       WRITE (6,FMT=8010) ef,workf,s
    END IF
8010 FORMAT (/,10x,'fermi energy=',f10.5,' har',3x,'work function=',&
                f10.5,' ev',/,10x,'number of valence electrons=',f10.5)
    IF (ABS(zcdiff).GT.5.0e-4) THEN
       CALL juDFT_error('Fermi-level determination did not converge'&
            ,hint ="change temperature or set input = F" ,calledby ="fergwt")
    ENDIF
    DO  jspin = 1,input%jspins
       IF ( mpi%irank == 0 ) WRITE (6,FMT=8020) jspin
8020   FORMAT (/,/,5x,'band-weighting factor for spin=',i5)
       DO  k = 1,kpts%nkpt
          nbnd = ne(k,jspin)
          IF ( mpi%irank == 0 ) WRITE (6,FMT=8030) k
8030      FORMAT (/,5x,'k-point=',i5,/)
          w_iks(:,k,jspin) = kpts%wtkpt(k)*w_iks(:,k,jspin)
          IF ( mpi%irank == 0) WRITE (6,FMT=8040) (w_iks(i,k,jspin),i=1,nbnd)
8040      FORMAT (5x,16f6.3)
       ENDDO
    ENDDO
    s1 = 0.
    s2 = 0.
    DO  jspin = 1,input%jspins
       s = 0.
       DO  k = 1,kpts%nkpt
          DO  i = 1,ne(k,jspin)
             s = s + w_iks(i,k,jspin)
             seigv = seigv + w_iks(i,k,jspin)*eig(i,k,jspin)
             en = eig(i,k,jspin)
             de = (en-ef)/input%delgau
             !     ---> correction term
             IF (ABS(de).LT.3.) THEN
                de = de*de
                s2 = s2 + EXP(-de)*kpts%wtkpt(k)
             END IF
          ENDDO
       ENDDO
       s1 = s1 + s
    ENDDO
    seigv = (2/input%jspins)*seigv
    seigv1 = (1/input%jspins)*fact1*s2
    chmom = s1 - input%jspins*s
    IF ( mpi%irank == 0 ) THEN
       WRITE (6,FMT=8050) seigv - seigv1,s1,chmom
       WRITE (16,FMT=8050) seigv - seigv1,s1,chmom
    END IF
8050 FORMAT (/,10x,'sum of eigenvalues-correction=',f12.5,/,10x,&
          'sum of weight                =',f12.5,/,10x,&
          'total moment                 =',f12.5,/)

  END SUBROUTINE fergwt
END MODULE m_fergwt

