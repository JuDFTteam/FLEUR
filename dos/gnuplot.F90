    !----------------------------------------------------------------------
! once the file "bands.1" and "bands.2" are created, activate with:
! gnuplot < band.gnu > band.ps
!----------------------------------------------------------------------
      SUBROUTINE write_gnu(
     >                     nosyp,d,ssy,name,jspins)
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: nosyp,jspins
      REAL,    INTENT (IN) :: d(nosyp)
      CHARACTER(len=1), INTENT (IN) :: ssy(nosyp)
      CHARACTER(len=8), INTENT (IN) :: name(10)

      INTEGER n,aoff,adel
      aoff = iachar('a')-1
      adel = iachar('a')-iachar('A')
      write(*,*) aoff,adel

      OPEN (27,file='band.gnu',status='unknown')
      WRITE (27,900)
      WRITE (27,901)
      WRITE (27,902)
      WRITE (27,903)
      WRITE (27,904) name(:)
      DO n = 1, nosyp
        WRITE (27,905) d(n),d(n)
      ENDDO
      WRITE (27,906) d(1),d(nosyp)
!
! nomal labels
!
      IF (iachar(ssy(1)) < aoff ) THEN
        WRITE (27,907) ssy(1),d(1)
      ELSE
       WRITE (27,907) " ",d(1)
      ENDIF
      DO n = 2, nosyp-1
        IF (iachar(ssy(n)) < aoff ) THEN
          WRITE (27,908) ssy(n),d(n)
        ELSE
          WRITE (27,908) " ",d(n)
        ENDIF
      ENDDO
      IF (iachar(ssy(nosyp)) < aoff ) THEN
        WRITE (27,909) ssy(nosyp),d(nosyp)
      ELSE
        WRITE (27,909) " ",d(nosyp)
      ENDIF
!
! greek labels
!
      DO n = 1, nosyp
        IF (iachar(ssy(n)) > aoff ) THEN
          WRITE (27,914) achar(iachar(ssy(n))-adel),d(n)
        ENDIF
      ENDDO
!
! now write the rest
      WRITE (27,910)
      WRITE (27,911) d(nosyp)+0.00001
      IF (jspins == 2) WRITE (27,912)
      WRITE (27,913)
      CLOSE (27)

 900  FORMAT ('set terminal postscript enhanced "Times-Roman" 20')
 901  FORMAT ('set xlabel ""')
 902  FORMAT ('set ylabel "E - E_F (eV)"')
 903  FORMAT ('set nokey')
 904  FORMAT ('set title "',10a8,'"')
 905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
 906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
#ifdef CPP_AIX
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \\')
 908  FORMAT ('           "',a1,'"',f9.5,', \\')
#else
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \')
 908  FORMAT ('           "',a1,'"',f9.5,', \')
#endif
 909  FORMAT ('           "',a1,'"',f9.5,'  )')
 910  FORMAT ('set ytics -8,2,4')
#ifdef CPP_AIX
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \\')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \\')
#else
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \')
#endif
 913  FORMAT ('"bands.1" using 1:($2+0.00)  w p pt  7 ps 0.5')
 914  FORMAT ('set label "',a1,'" at ',f9.5,
     +        ', -9.65 center font "Symbol,20"')

      END SUBROUTINE write_gnu
