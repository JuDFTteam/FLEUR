      MODULE m_chkmt
      use m_juDFT
      private
      public chkmt
      INTEGER,PARAMETER:: MAX_CHECK=100
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
      CONTAINS
      SUBROUTINE chkmt(&
     &                 atoms,input,vacuum,cell,oneD,&
     &                 l_gga,noel,l_test,&
     &                 kmax,dtild,dvac1,lmax1,jri1,rmt1,dx1)

      USE m_types
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_input),INTENT(IN) :: input
      TYPE(t_vacuum),INTENT(IN):: vacuum
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_oneD),INTENT(IN)  :: oneD
      CHARACTER*3, INTENT (IN) :: noel(atoms%ntype)
      LOGICAL, INTENT (IN)     :: l_gga,l_test
      REAL,    INTENT (OUT)    :: kmax,dtild,dvac1
!     ..
!     .. Array Arguments ..
      INTEGER, INTENT (OUT)    :: lmax1(atoms%ntype),jri1(atoms%ntype)
      REAL,    INTENT (OUT)    :: rmt1(atoms%ntype),dx1(atoms%ntype)
!     ..
!     .. Local Scalars ..
      INTEGER na,n,nna,nn,n1,nn1,k1,k2,k3
      INTEGER i,j,jri11,lmax11
      REAL    sum,sss,xmin,dx11,rkm,fac,sum_r,fac_1,fac_2
      LOGICAL error
      REAL    distvec(3),vec(3,26)
!     ..
!     .. Local Arrays ..
      INTEGER minni(2)
      REAL    dist(atoms%ntype,atoms%ntype),dist1(atoms%ntype,atoms%ntype),t_rmt(0:103)

      IF (l_test.and.atoms%ntype>MAX_CHECK) THEN
          write(6,*) "No test for MT radii performed"
          write(6,*) "Change MAX_CHECK in chkmt if needed"
          return
      ENDIF
!
! typical muffin-tin radii
!
      i=0
      DO k1=-1,1
         DO k2=-1,1
            DO k3=-1,1
               if (k1==0.and.k2==0.and.k3==0) cycle
               i=i+1
               vec(:,i)=matmul(cell%amat,(/k1,k2,k3/))
            ENDDO
         ENDDO
      ENDDO
      t_rmt(0:103) = 2.3 ! default value
      t_rmt(1) = 1.0 ; t_rmt(5:9) = 1.3 ; t_rmt(16:17) = 1.8

      error=.false.
      dist(:,:) = 9.99e19
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(atoms,dist,error,l_test,input,oneD,vacuum,vec)&
      !$OMP PRIVATE(n,n1,na,nn,nn1,nna,sss,distvec,k1)
      DO n=1,atoms%ntype
        DO n1=1,atoms%neq(n)
          na=n1+sum(atoms%neq(:n-1))
!
! check distance to other atoms:
!
          DO nn=n,atoms%ntype
            DO nn1=1,atoms%neq(nn)
              nna=nn1+sum(atoms%neq(:nn-1))
              sss=1.E19
              distvec=atoms%pos(:,na)-atoms%pos(:,nna)
              if (nna.ne.na) sss=distvec(1)**2+distvec(2)**2+distvec(3)**2
              DO k1=1,26
                 sss=min(sss,(distvec(1)+vec(1,k1))**2+(distvec(2)+vec(2,k1))**2+(distvec(3)+vec(3,k1))**2)
              ENDDO
              dist(n,nn) = min( dist(n,nn),sqrt(sss) ) 
              dist(nn,n) = dist(n,nn) 
              IF ( sss.LE.(atoms%rmt(nn)+atoms%rmt(n))**2 ) THEN
                 error=.true.
                 IF (l_test)&
     &           WRITE(6,240) nn,nn1,(atoms%pos(i,nna),i=1,3),atoms%rmt(nn),&
     &                        n ,n1 ,(atoms%pos(i,na),i=1,3),atoms%rmt(n )
              ENDIF

            ENDDO ! nn1
          ENDDO   ! nn
!
! distance to vacuum
!
          IF (input%film) THEN
             IF (oneD%odd%d1) THEN
                IF ((sqrt(atoms%pos(1,na)**2+atoms%pos(2,na)**2)+&
     &               atoms%rmt(n)).GT.vacuum%dvac/2.) THEN
                   error=.true.
                   WRITE(6,241) n ,n1
                   WRITE(6,*) sqrt(atoms%pos(1,na)**2+atoms%pos(2,na)**2),&
     &                  atoms%rmt(n),vacuum%dvac/2.
                END IF
             ELSE
                IF ( ( (atoms%pos(3,na)+atoms%rmt(n) ).GT. vacuum%dvac/2.).OR.&
     &               ( (atoms%pos(3,na)-atoms%rmt(n) ).LT.-vacuum%dvac/2.) ) THEN
                   error=.true.
                   WRITE(6,241) n ,n1
                   WRITE(6,*) atoms%pos(3,na),atoms%rmt(n),vacuum%dvac/2.
                ENDIF
             ENDIF
          END IF
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
       
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
!     DO n = 1, ntype
!       WRITE (*,'(12f12.6)') dist(:,n)
!     ENDDO
      dist1 = dist
      rmt1(:) = 999.
      WRITE (6,*) '----------------------------------------------------'
      WRITE (6,*) 'Suggested values for input: '
      WRITE (6,*) 

      IF (input%film) THEN
        fac = 0.95
      ELSE
        fac = 0.975
      ENDIF

      minni = minloc(dist)       ! minni(1) and minni(2) are the indices of the closest atoms
      xmin  = minval(dist)       ! xmin is their distance

      DO WHILE ( xmin < 999.0 )

        sum_r = 1.0 / ( t_rmt(atoms%nz(minni(1))) + t_rmt(atoms%nz(minni(2))) )
        fac_1 = t_rmt(atoms%nz(minni(1))) * sum_r
        fac_2 = t_rmt(atoms%nz(minni(2))) * sum_r

        IF (rmt1(minni(1)) > 990.) THEN         ! if not set, determine MTR 
          IF (rmt1(minni(2)) > 990.) THEN       ! both not set, choose in between
            rmt1(minni(1)) = fac * xmin * fac_1 ! / 2
            rmt1(minni(2)) = fac * xmin * fac_2 ! / 2
          ELSE
            rmt1(minni(1)) = fac * ( xmin - rmt1(minni(2)) )
            IF (2*rmt1(minni(1)).GT.dist(minni(1),minni(1))) THEN
              rmt1(minni(1)) = fac * dist(minni(1),minni(1)) / 2
            ENDIF
          ENDIF
        ELSEIF (rmt1(minni(2)) > 990.) THEN
          rmt1(minni(2)) = fac * ( xmin - rmt1(minni(1)) )
          IF (2*rmt1(minni(2)).GT.dist(minni(2),minni(2))) THEN
            rmt1(minni(2)) = fac * dist(minni(2),minni(2)) / 2
          ENDIF
        ENDIF

        dist(minni(1),minni(1)) = 999.0
        dist(minni(2),minni(1)) = 999.0
        dist(minni(1),minni(2)) = 999.0
        dist(minni(2),minni(2)) = 999.0

        DO j = 1, 2
          DO n = 1, atoms%ntype
            IF (atoms%nz(n) == atoms%nz(minni(j))) THEN
              IF (rmt1(n) > 990.) THEN 
                rmt1(n) = rmt1(minni(j))
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        minni = minloc(dist)
        xmin  = minval(dist)

      ENDDO

      dvac1 = 0.0
      rkm = 0.0
      na = 0
      WRITE (6,230)
      DO n= 1,atoms%ntype
!
!--> determine M.T. radii 
!
        DO j= 1,atoms%ntype
           dist(j,n) = rmt1(n)+rmt1(j)
           IF ( dist1(j,n)-dist(j,n) < 0.0 ) THEN
             WRITE(*,*) j,n,dist1(j,n)-dist(j,n)
             rmt1(n) = fac * dist1(j,n) / 2.
             rmt1(j) = fac * dist1(j,n) / 2.
           ENDIF
        ENDDO
        IF (input%film) THEN
          DO nn = 1, atoms%neq(n)
            na = na + 1
            IF (oneD%odd%d1) THEN
               dvac1 = max( dvac1, sqrt(atoms%pos(1,na)**2+atoms%pos(2,na)**2)&
     &              +rmt1(n) )
            ELSE
               dvac1 = max( dvac1, abs(atoms%pos(3,na))+rmt1(n) )
            END IF
          ENDDO
        ENDIF 
!
!--> calculate jri1, dx1 and lmax1
!
        IF (rmt1(n).LT.1.8) THEN
          lmax11 = 6
        ELSEIF (rmt1(n).LT.2.4) THEN
          lmax11 = 8
        ELSEIF (rmt1(n).LT.2.8) THEN
          lmax11 = 10
        ELSE
          WRITE (6,'("Atom Nr.",i3,"( ",a3,") has a M.T. radius of",&
     &                                     f8.4)') n,noel(n),rmt1(n)
          WRITE (6,'("that was truncated to 2.8")')
          rmt1(n) = 2.8
          lmax11 = 10
        ENDIF
        IF (l_gga) THEN
          jri11 = nint( 330*rmt1(n) ) 
        ELSE
          jri11 = nint( 220*rmt1(n) ) 
        ENDIF
        jri11 = nint( jri11*0.5 ) * 2 + 1
        dx11 =  log(3200*atoms%nz(n)*rmt1(n))/(jri11-1)
        rkm = max( rkm , lmax11/rmt1(n) )
           
        WRITE (6,9070) noel(n),atoms%nz(n),lmax11,jri11,rmt1(n),dx11
 9070   FORMAT (a3,i3,2i5,2f10.6)
        dx1(n) = dx11 ; lmax1(n) = lmax11 ; jri1(n) = jri11
        
      ENDDO ! loop over atom types
  230 FORMAT ('Atom Z  lmax jri    rmt         dx')

      IF (input%film) THEN
        dvac1 = 2* (dvac1+0.3)
        dtild = dvac1 + 1.5 * maxval( rmt1(:) )
        WRITE (6,'("vacuum distance dvac =",f10.5)') dvac1
        WRITE (6,'("extra vac.dist. dtild=",f10.5)') dtild
      ENDIF
      WRITE (6,'("k_max =",f8.5)') rkm
      WRITE (6,'("G_max =",f8.5)') 3*rkm
      kmax = rkm

      WRITE (6,*) '----------------------------------------------------'
 
      IF ( error.AND.l_test )  CALL juDFT_error&
     &     ("Error checking M.T. radii",calledby ="chkmt")
  240 FORMAT('   error in muffin tin radii  , pair  ',2(/,2i5,4f10.5))
  241 FORMAT('   error: atom ',i3,' # ',i3,'reaches out into vaccuum')

      END SUBROUTINE chkmt
      END MODULE m_chkmt
