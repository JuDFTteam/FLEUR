module m_gnuplot_BS
CONTAINS

!----------------------------------------------------------------------
! once the file "bands.1" and "bands.2" are created, activate with:
! gnuplot < band.gnu > band.ps
!----------------------------------------------------------------------
      SUBROUTINE gnuplot_BS(kpts,cell,jspins)
      use m_types_cell
      use m_types_kpts
      IMPLICIT NONE
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      INTEGER, INTENT (IN)     :: jspins
      CHARACTER(len=1) :: ssy
      real   :: d(kpts%numSpecialPoints),vkr(3),vkr_prev(3)

      INTEGER n,aoff,adel,k
      aoff = iachar('a')-1
      adel = iachar('a')-iachar('A')
      !write(*,*) aoff,adel

      !Generate distances
      d(1)=0.0
      vkr_prev=matmul(kpts%specialPoints(:,1),cell%bmat)
      DO k=2,kpts%numSpecialPoints
        vkr=matmul(kpts%specialPoints(:,k),cell%bmat)
        d(k)=d(k-1)+sqrt(dot_product(vkr-vkr_prev,vkr-vkr_prev))
        vkr_prev=vkr
      ENDDO


      OPEN (27,file='band.gnu',status='unknown')
      WRITE (27,900)
      WRITE (27,901)
      WRITE (27,902)
      WRITE (27,903)
      WRITE (27,904) kpts%kptsName
      DO n = 1, kpts%numSpecialPoints
        WRITE (27,905) d(n),d(n)
      ENDDO
      WRITE (27,906) d(1),d(size(d))
!
! nomal labels
!
      ssy=kpts%specialPointNames(1)
      IF (iachar(ssy) < aoff ) THEN
        WRITE (27,907) ssy,d(1)
      ELSE
       WRITE (27,907) " ",d(1)
     ENDIF
      DO n = 2, kpts%numSpecialPoints-1
        ssy=kpts%specialPointNames(n)
        IF (iachar(ssy) < aoff ) THEN
          WRITE (27,908) ssy,d(n)
        ELSE
          WRITE (27,908) " ",d(n)
        ENDIF
      ENDDO
      ssy=kpts%specialPointNames(n)
      IF (iachar(ssy) < aoff ) THEN
        WRITE (27,909) ssy,d(n)
      ELSE
        WRITE (27,909) " ",d(n)
      ENDIF

      DO n=1,kpts%numSpecialPoints
        ssy=kpts%specialPointNames(n)
        IF (iachar(ssy) > aoff ) THEN
          WRITE (27,914) achar(iachar(ssy)-adel),d(n)
        ENDIF
      ENDDO

!
! now write the rest
      WRITE (27,910)
      WRITE (27,911) d(size(d))+0.00001
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
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \ ')
 908  FORMAT ('           "',a1,'"',f9.5,', \ ')
#endif
 909  FORMAT ('           "',a1,'"',f9.5,'  )')
 910  FORMAT ('set ytics -8,2,4')
#ifdef CPP_AIX
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \\')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \\')
#else
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \ ')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \ ')
#endif
 913  FORMAT ('"bands.1" using 1:($2+0.00)  w p pt  7 ps 0.5')
 914  FORMAT ('set label "',a1,'" at ',f9.5,', -9.65 center font "Symbol,20"')

    !call system("gnuplot <band.gnu >band.ps")

  END SUBROUTINE gnuplot_BS
end MODULE
