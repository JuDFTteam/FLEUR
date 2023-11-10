!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_spg2set
      use m_juDFT
!-------------------------------------------------------------------------+
! determine the rotation martrices (mrot) and non-symorphic translations  !
! for a given 2D symmetry (defined by a number n2spg) plus eventually     !
! inversion (invs) or z-reflection (zrfs) symmetry.                       !
! Plane groups 1-17 are according to Int. Tables of Crystallography       !
! 18-20 : pm, pg and cm with mirror plane y -> -y  (as used in old spgset)!
! 21-25 : p3, p3m1, p31m, p6 and p6mm  with sharp angle (lattice type hex)!
!                                                                   gb`02 !
!-------------------------------------------------------------------------+
      CONTAINS
      SUBROUTINE spg2set(
     >                   nop,zrfs,invs,namgrp,latnam,
     <                   mrot,tau,nop2,symor)

      USE m_constants
      USE m_symdata, ONLY : gen2,tau2,spg2,gnt2,namgr2,nammap,ord2
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nop
      LOGICAL, INTENT (IN) :: invs,zrfs
      CHARACTER(len=4),INTENT (IN) :: namgrp
      CHARACTER(len=3),INTENT (IN) :: latnam

      INTEGER, INTENT (OUT) :: nop2
      INTEGER, INTENT (OUT) :: mrot(3,3,nop)
      REAL,    INTENT (OUT) :: tau(3,nop)
      LOGICAL, INTENT (OUT) :: symor


      INTEGER n2spg, ngen, igen, n1, n2, n, i, j, d, t
      INTEGER mt(3,3), multab(nop)
      LOGICAL l_new
      CHARACTER(len=3) type

!mod: INTEGER spg2(3,25)             ! generators for 2d space groups
!mod: INTEGER gnt2(3,25)             ! translations for 2d space groups
!mod: INTEGER gen2(2,2,9)            ! rotation matrices for the generators
!mod: REAL    tau2(2,3)              ! translations for the generators
!mod: CHARACTER(len=4) :: namgr2(25) ! names of 2d space groups
!
! Determine number of 2d space group
!
      n2spg=0
      DO i = 1, 20
        IF (namgrp.EQ.nammap(i)) n2spg = i
      ENDDO
      IF ((latnam.EQ.'hex').AND.(n2spg.GT.12).
     +                      AND.(n2spg.LT.18)) THEN
        n2spg = n2spg + 8
      ENDIF
      IF (n2spg == 0)  THEN
        print *,latnam,namgrp 
        CALL juDFT_error("2D-symmetry group not found!"
     +     ,calledby ="spg2set")
      ENDIF
!
! Determine number of generators ngen
!
      ngen = 1
      DO igen = 1, 3
        IF (spg2(igen,n2spg).NE.0) ngen = ngen + 1
      ENDDO
      nop2 = ord2(n2spg)
!
! make 3d rotation matrices for the generators
!
      symor = .true.
      mrot(:,:,1) = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))  ! Identity is the first generator
      tau(:,1) = 0.0
      DO igen = 2, ngen
        mrot(:,:,igen) = reshape((/0,0,0,0,0,0,0,0,1/),(/3,3/))
        mrot(1:2,1:2,igen) = gen2(1:2,1:2,spg2(igen-1,n2spg))
        tau(:,igen) =  0.0
        IF (gnt2(igen-1,n2spg).NE.0) THEN
          symor = .false.
          tau(1:2,igen) = tau2(1:2,gnt2(igen-1,n2spg))
        ENDIF
      ENDDO
!
! now close the group
!
      mrot(:,:,ngen+1:nop) = 0 ; tau(:,ngen+1:nop) = 0.0
  10  CONTINUE
      rm1 : DO n1 = 1,ngen
         rm2 : DO  n2 = 1,ngen

            !CALL matmul1(mrot(1,1,n1),mrot(1,1,n2),mt)
            mt=matmul(mrot(:,:,n1),mrot(:,:,n2))
            rm3 : DO n = 1,nop
               DO i = 1,3
                  DO j = 1,3
                     IF (mt(i,j).NE.mrot(i,j,n)) CYCLE rm3
                  ENDDO
               ENDDO
               CYCLE rm2
            ENDDO rm3

! -->       new element found
            ngen = ngen + 1
            IF (ngen.gt.nop) THEN
             WRITE(oUnit,'(a7,i4,a7,i4)') 'ngen = ',ngen,' nop = ',nop
                CALL juDFT_error("ngen > nop",calledby="spg2set")
            ENDIF

            CALL matmul4(mrot(1,1,n1),tau(1,n1),
     >                   mrot(1,1,n2),tau(1,n2),
     <                   mrot(1,1,ngen),tau(1,ngen))
            GOTO 10

         ENDDO rm2
      ENDDO rm1
!
! add inversion or z-reflection symmetry
!
      IF ( (nop2.EQ.ngen).AND.(zrfs.OR.invs) ) THEN
        l_new = .true.
        IF (invs) THEN
          ngen = ngen + 1
          mrot(:,:,ngen) = - mrot(:,:,1)                     ! I
          tau(:,ngen) =  0.0
          IF (spg2(1,n2spg).EQ.2) l_new = .false.            ! if c_2 & I are generators
        ENDIF                                                ! m_z is no new generator
        IF (zrfs.AND.l_new) THEN
          ngen = ngen + 1
          mrot(:,:,ngen) =  reshape((/1,0,0,0,1,0,0,0,-1/),(/3,3/))  !m_z
          tau(:,ngen) =  0.0
        ENDIF
        GOTO 10                                              ! now close the 3D group
      ENDIF
!
! Output the symmetry elements
!
      IF (nop.NE.ngen) THEN
          WRITE (oUnit,*) 'nop =',nop,' =/= ngen = ',ngen
           CALL juDFT_error("nop =/= ngen",calledby="spg2set")
      ENDIF
      WRITE (oUnit,FMT=8010) namgr2(n2spg),invs,zrfs,nop,symor
 8010 FORMAT (/,/,' space group: ',a4,' invs=',l1,' zrfs=',l1,/,
     +  ' number of operations=',i3,/,' symmorphic=',l1,/,/)
      WRITE (oUnit,'("Number of 2D operations=",i3)') nop2
      DO n = 1,nop
!
! Determine the kind of symmetry operation we have here
!
         d = mrot(1,1,n)*mrot(2,2,n)*mrot(3,3,n) +
     +       mrot(1,2,n)*mrot(2,3,n)*mrot(3,1,n) +
     +       mrot(2,1,n)*mrot(3,2,n)*mrot(1,3,n) -
     +       mrot(1,3,n)*mrot(2,2,n)*mrot(3,1,n) -
     +       mrot(2,3,n)*mrot(3,2,n)*mrot(1,1,n) -
     +       mrot(2,1,n)*mrot(1,2,n)*mrot(3,3,n)
         t =  mrot(1,1,n) + mrot(2,2,n) + mrot(3,3,n)

         IF (d.EQ.-1) THEN
           type = 'm  '
           IF (t.EQ.-3) type = 'I  '
         ELSEIF (d.EQ.1) THEN
           IF (t.EQ.-1) type = 'c_2'
           IF (t.EQ. 0) type = 'c_3'
           IF (t.EQ. 1) type = 'c_4'
           IF (t.EQ. 2) type = 'c_6'
           IF (t.EQ. 3) type = 'E  '
         ELSE
            CALL juDFT_error("determinant=/=+/- 1",calledby ="spg2set")
         ENDIF

         WRITE (oUnit,FMT=8020) n, type
 8020    FORMAT (/,1x,i3,' : ',a3)
         DO i = 1,3
            WRITE (oUnit,FMT=8030) (mrot(i,j,n),j=1,3),tau(i,n)
         ENDDO
 8030    FORMAT (5x,3i3,3x,f4.1)
      ENDDO
c
c     check closure
c
      WRITE (oUnit,FMT=8040)
 8040 FORMAT (/,/,' multiplication table',/,/)

      op1 : DO n1 = 1,nop
         op2 : DO  n2 = 1,nop

            !CALL matmul1(mrot(1,1,n1),mrot(1,1,n2),mt)
            mt=matmul(mrot(:,:,n1),mrot(:,:,n2))
            op3 : DO n = 1,nop
               DO i = 1,3
                  DO j = 1,3
                     IF (mt(i,j).NE.mrot(i,j,n)) CYCLE op3
                  ENDDO
               ENDDO
               multab(n2) = n
               CYCLE op2
            ENDDO op3

            WRITE (oUnit,FMT=8050) n1,n2
 8050       FORMAT (' error - n1,n2=',2i3)
             CALL juDFT_error("mult",calledby="spg2set")
         ENDDO op2
         WRITE (oUnit,FMT=8060) (multab(n),n=1,nop)
 8060    FORMAT (1x,48i2)
      ENDDO op1
      CONTAINS
      SUBROUTINE matmul4 (ma,ta,mb,tb,mc,tc)

         IMPLICIT NONE
         INTEGER, INTENT (IN) :: ma(3,3),mb(3,3)
         REAL,    INTENT (IN) :: ta(3),tb(3)
         INTEGER, INTENT (OUT):: mc(3,3)
         REAL,    INTENT (OUT):: tc(3)
   
         INTEGER :: i,j,k
         REAL x,xa(4,4),xb(4,4),xc(4,4)
   
         xa(:,:) = 0.0 ; xa(4,4) = 1.0 ; xb(:,:) = xa(:,:)
         xa(1:3,1:3) = real(ma(1:3,1:3)) ; xa(1:3,4) = ta(:) 
         xb(1:3,1:3) = real(mb(1:3,1:3)) ; xb(1:3,4) = tb(:)
         DO i=1,4
            DO k=1,4
            x=0.e0
               DO j=1,4
                  x = x + xa(i,j)*xb(j,k)
               ENDDO
            xc(i,k) = x
            END DO
         END DO
         mc(1:3,1:3) = nint(xc(1:3,1:3)) ; tc(:) = xc(1:3,4) 
         DO i = 1,3
           IF (tc(i).GT.1.0) THEN
             tc(i) = tc(i) - int(tc(i))
           ELSEIF (tc(i).LT.0.0) THEN
             tc(i) = tc(i) + int(tc(i)) + 1
           ENDIF
         ENDDO
   
         END SUBROUTINE matmul4
      END SUBROUTINE spg2set
      END MODULE m_spg2set
