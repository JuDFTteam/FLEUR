MODULE m_convn
  use m_juDFT
CONTAINS
  SUBROUTINE convn(atoms,stars)
    !
    !     ***********************************************************
    !     determines the optimum values for the convergence parameter
    !     for each atom type using the criterion discussed in
    !     m. weinert, j. math. phys. 22, 2433 (1981).  each sphere
    !     and l component may have different values.  (psqpw changed
    !     to allow this option).
    !          m. weinert july 1982
    !     ***********************************************************
    USE m_types
    USE m_constants
    IMPLICIT NONE
    !     ..
    TYPE(t_atoms),INTENT(INOUT)  :: atoms
    TYPE(t_stars),INTENT(IN)     :: stars
    !     .. Local Scalars ..
    REAL sck,z0
    INTEGER i,l,n,n1,nc
    !     ..
    !     .. Local Arrays ..
    REAL z(17)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC min0
    !     ..
    !     .. Data statements ..
    DATA z/6.9e0,8.1e0,9.3e0,10.5e0,11.6e0,12.7e0,13.9e0,15.0e0,&
         &     16.1e0,17.2e0,18.3e0,19.4e0,20.5e0,21.6e0,22.7e0,23.7e0,&
         &     24.8e0/,z0/5.7e0/
    !     ..
    !--->    read in values of ncv (if ncv(1).le.0, calculate best values)
    !      read(5,1000) (ncv(n),n=1,ntype)
    !      if(ncv(1).le.0) go to 2
    !      n1=ncv(1)
    !      do 1 n=2,ntype
    !    1 if(ncv(n).le.0) ncv(n)=n1
    !      go to 5
    !--->    calculate values
    !    2 continue
    !
    DO 20 n = 1,atoms%ntype
       sck = stars%gmax*atoms%rmt(n)
       IF (sck.LT.z0) GO TO 60
       DO 10 i = 1,17
          IF (sck.GT.z(i)) GO TO 10
          atoms%ncv(n) = i
          GO TO 20
10        CONTINUE
          n1 = 0.9e0* (sck-z(17))
          atoms%ncv(n) = 18 + n1
20        CONTINUE
          !--->    output and make sure ncv(n).le.ncvd
30        CONTINUE
          WRITE (6,FMT=8010)
          WRITE (16,FMT=8010)
          DO 40 n = 1,atoms%ntype
             nc = atoms%ncv(n)
             l = nc - 1
             WRITE (6,FMT=8020) n,nc,l
             WRITE (16,FMT=8020) n,nc,l
40           CONTINUE
             l = ncvd_dim - 1
             WRITE (6,FMT=8030) ncvd_dim,l
             WRITE (16,FMT=8030) ncvd_dim,l
             DO 50 n = 1,atoms%ntype
                atoms%ncv(n) = min0(atoms%ncv(n),ncvd_dim)
50              CONTINUE
                RETURN
60              WRITE (6,FMT=8040) n,sck
                WRITE (16,FMT=8040) n,sck
                CALL juDFT_error("ncv",calledby="convn")
8000            FORMAT (10i5)
8010            FORMAT (/,/,10x,'convergence parameters for the pseudocharge',&
                     &       ' density expansion',/,10x,'atom',5x,'parameter',5x,&
                     &       'max. l to include',/)
8020            FORMAT (10x,i3,9x,i3,13x,i3)
8030            FORMAT (10x,'max values allowed: ncvd=',i3,', l=',i3,/)
8040            FORMAT (/,/,10x,'atom type',i3,' has rkmax=',f6.4,/,10x,&
                     &       '$$$ stop ncv error')
              END SUBROUTINE convn
            END MODULE m_convn
