      MODULE m_writestruct
!********************************************************************
!     write structure in povray format
!     usage:  povray +i struct.pov
!     maybe:   +H750 +W1000 +L/usr/local/lib/povray3/include 
!********************************************************************
      CONTAINS
      SUBROUTINE write_struct(
     >                        ntype,nat,neq,
     >                        rmt,pos,natmap,amat)

     
      IMPLICIT NONE

!===> Arguments

      INTEGER, INTENT (IN) :: ntype,nat
      INTEGER, INTENT (IN) :: natmap(nat),neq(ntype)
      REAL,    INTENT (IN) :: rmt(ntype),pos(3,nat),amat(3,3)

!===> Local Variables

      INTEGER       :: i,n,na,nn
      REAL          :: posc(3),col(3,12)

      DATA col/ 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0,
     +          0.0,1.0,1.0, 1.0,0.0,1.0, 1.0,1.0,0.0,
     +          0.0,0.5,1.0, 0.5,0.0,1.0, 0.5,1.0,0.0,
     +          0.0,1.0,0.5, 1.0,0.0,0.5, 1.0,0.5,0.0/

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      OPEN (45,file='struct.pov',form='formatted',status='unknown')
!
! --> header information
!
      WRITE (45,*) '#include "colors.inc"'
      WRITE (45,*) '#include "shapes.inc"'
      WRITE (45,*) 'global_settings { max_trace_level 20 '
      WRITE (45,*) '                  assumed_gamma 2.2 }'
      WRITE (45,*) 'light_source { <60,-10,-10>'
      WRITE (45,*) '              color rgb <2.5,2.5,2.5> }'
      WRITE (45,*) 'camera { location <90,10,10>'
      WRITE (45,*) '         look_at <0.0,0.0,0.0> angle 20 }'
      WRITE (45,*) 'background {color White}'
!
! --> colors etc.
!
      WRITE (45,*) ' #declare R1 =  pigment{ color Black } '
      WRITE (45,*) ' #declare Rd =  0.05 ;'
      DO nn = 1,ntype
        i = mod(nn-1,12) + 1
        IF ( nn < 10 ) THEN
          WRITE (45,1005) nn,col(:,i)
        ELSE
          WRITE (45,1006) nn,col(:,i)
        ENDIF
      ENDDO
 1005 FORMAT ('#declare Acol',i1,'= color rgb <',3f4.1,'>;')     
 1006 FORMAT ('#declare Acol',i2,'= color rgb <',3f4.1,'>;')     

!---> output the atomic definitions 

      DO nn = 1,ntype
        IF ( nn < 10 ) THEN
          WRITE (45,1010) nn,nn
        ELSE
          WRITE (45,1011) nn,nn
        ENDIF
      ENDDO
 1010 FORMAT ('#declare Atom',i1,' = pigment { Acol',i1,' }')
 1011 FORMAT ('#declare Atom',i2,' = pigment { Acol',i2,' }')
      
      WRITE (45,1015)
 1015 FORMAT (/,'#declare Ascale = 0.80 ;')

      DO nn = 1,ntype
        IF ( nn < 10 ) THEN
          WRITE (45,1030) nn,rmt(nn)
        ELSE
          WRITE (45,1031) nn,rmt(nn)
        ENDIF
      ENDDO
 1030 FORMAT ('#declare Asize',i1,' = ',f6.3,'*Ascale ;')
 1031 FORMAT ('#declare Asize',i2,' = ',f6.3,'*Ascale ;')


!---> output the atomic positions

      na = 0
      DO nn = 1, ntype
         DO n = 1, neq(nn)
           !CALL cotra0(pos(:,natmap(na+n)),posc,amat)
           posc=matmul(amat,pos(:,natmap(na+n)))
!           DO i = 1, 2
!             IF (posc(i).LT.0) posc(i) = posc(i) + amat(i,i)
!           ENDDO
           IF ( nn < 10 ) THEN
            WRITE (45,1020)
     &           posc(:),nn,nn,nn,natmap(na+n)
          ELSE
            WRITE (45,1021)
     &           posc(:),nn,nn,nn,natmap(na+n)
          ENDIF
         ENDDO
         na = na + neq(nn)
      ENDDO
 1020 FORMAT('sphere { <',3(f8.3,','),'>, Asize',i1,' texture { Atom',
     &       i1,' } } // ',2i4)
 1021 FORMAT('sphere { <',3(f8.3,','),'>, Asize',i2,' texture { Atom',
     &       i2,' } } // ',2i4)


!---> output the primitive cell

      DO i = 1,3
         WRITE (45,2000) 0.0, 0.0, 0.0, amat(:,i)
         WRITE (45,2000) -amat(:,i)/2, amat(:,i)/2
      ENDDO
 2000 FORMAT('cylinder { ',2('<',3(f8.3,','),'>,'),
     &                   ' Rd texture { R1 } }')

      CLOSE (45)

      END SUBROUTINE write_struct
      END MODULE m_writestruct
