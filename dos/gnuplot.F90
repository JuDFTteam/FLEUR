module m_gnuplot_BS
CONTAINS

!----------------------------------------------------------------------
! once the file "bands.1" and "bands.2" are created, activate with:
! gnuplot < band.gnu > band.ps
!----------------------------------------------------------------------
      SUBROUTINE gnuplot_BS(kpts,title,cell,jspins)
      use m_types_cell
      use m_types_kpts
      IMPLICIT NONE
      type(t_kpts),intent(in)  :: kpts
      type(t_cell),intent(in)  :: cell
      INTEGER, INTENT (IN)     :: jspins
      CHARACTER(LEN=*), INTENT(IN) :: title
      CHARACTER(len=1) :: ssy
      real   :: vkr(3),vkr_prev(3)

      INTEGER n,aoff,adel,k
      INTEGER numSpecialPoints, addFront, addBack

      REAL,              ALLOCATABLE :: d(:)
      REAL,              ALLOCATABLE :: specialPoints(:,:)
      CHARACTER(LEN=50), ALLOCATABLE :: specialPointNames(:)

      aoff = iachar('a')-1
      adel = iachar('a')-iachar('A')
      !write(*,*) aoff,adel

      numSpecialPoints = kpts%numSpecialPoints
      addFront = 0
      addBack = 0
      
      IF (kpts%numSpecialPoints.EQ.0) THEN
         numSpecialPoints = 2
         addFront = 1
         addBack = 1
      ELSE
         IF (kpts%specialPointIndices(1).NE.1) THEN
            addFront = 1
            numSpecialPoints = numSpecialPoints + 1
         END IF
         IF (kpts%specialPointIndices(kpts%numSpecialPoints).NE.kpts%nkpt) THEN
            addBack = 1
            numSpecialPoints = numSpecialPoints + 1
         END IF
      END IF

      ALLOCATE (specialPoints(3,numSpecialPoints))
      ALLOCATE (specialPointNames(numSpecialPoints))
      ALLOCATE (d(numSpecialPoints))
      IF(addFront.EQ.1) THEN
         specialPoints(:,1) = kpts%bk(:,1)
         specialPointNames(1) = '-'
      END IF
      IF(ALLOCATED(kpts%specialPoints)) THEN
         specialPoints(:,1+addFront:kpts%numSpecialPoints+addFront) = kpts%specialPoints(:,1:kpts%numSpecialPoints)
         specialPointNames(1+addFront:kpts%numSpecialPoints+addFront) = kpts%specialPointNames(1:kpts%numSpecialPoints)
      END IF
      IF(addBack.EQ.1) THEN
         specialPoints(:,kpts%numSpecialPoints+addFront+1) = kpts%bk(:,kpts%nkpt)
         specialPointNames(kpts%numSpecialPoints+addFront+1) = '-'
      END IF



      !Generate distances
      d(1)=0.0
      vkr_prev=matmul(specialPoints(:,1),cell%bmat)
      DO k=2,numSpecialPoints
        vkr=matmul(specialPoints(:,k),cell%bmat)
        d(k)=d(k-1)+sqrt(dot_product(vkr-vkr_prev,vkr-vkr_prev))
        vkr_prev=vkr
      ENDDO


      OPEN (27,file='band.gnu',status='unknown')
      WRITE (27,900)
      WRITE (27,901)
      WRITE (27,902)
      WRITE (27,903)
      WRITE (27,'(3a)') 'set title "', TRIM(ADJUSTL(title)), '"'
      DO n = 1, numSpecialPoints
        WRITE (27,905) d(n),d(n)
      ENDDO
      WRITE (27,906) d(1),d(size(d))
!
! nomal labels
!
      ssy=specialPointNames(1)
      IF (iachar(ssy) < aoff ) THEN
        WRITE (27,907) ssy,d(1)
      ELSE
       WRITE (27,907) " ",d(1)
     ENDIF
      DO n = 2, numSpecialPoints-1
        ssy=specialPointNames(n)
        IF (iachar(ssy) < aoff ) THEN
          WRITE (27,908) ssy,d(n)
        ELSE
          WRITE (27,908) " ",d(n)
        ENDIF
      ENDDO
      ssy=specialPointNames(n)
      IF (iachar(ssy) < aoff ) THEN
        WRITE (27,909) ssy,d(n)
      ELSE
        WRITE (27,909) " ",d(n)
      ENDIF

      DO n=1,numSpecialPoints
        ssy=specialPointNames(n)
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

      DEALLOCATE (d,specialPoints,specialPointNames)

 900  FORMAT ('set terminal postscript enhanced color "Times-Roman" 20')
 901  FORMAT ('set xlabel ""')
 902  FORMAT ('set ylabel "E - E_F (eV)"')
 903  FORMAT ('set nokey')
 905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
 906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
#if (defined(CPP_AIX) || defined(__PGI))
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \\')
 908  FORMAT ('           "',a1,'"',f9.5,', \\')
#else
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', \')
 908  FORMAT ('           "',a1,'"',f9.5,', \')
#endif
 909  FORMAT ('           "',a1,'"',f9.5,'  )')
 910  FORMAT ('set ytics -8,2,4')
#if (defined(CPP_AIX) || defined(__PGI))
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \\')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \\')
#else
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] \')
 912  FORMAT ('"bands.2" using 1:($2+0.00)  w p pt 12 ps 0.5, \')
#endif
 913  FORMAT ('"bands.1" using 1:($2+0.00)  w p pt  7 ps 0.5')
 914  FORMAT ('set label "',a1,'" at ',f9.5,', -9.65 center font "Symbol,20"')

    !call system("gnuplot <band.gnu >band.ps")

  END SUBROUTINE gnuplot_BS

  SUBROUTINE write_gnu_sc(banddos,kpts,title,cell,jspins)
    use m_types_cell
    use m_types_kpts
    use m_types_banddos
    USE m_inv3
    USE m_constants
    IMPLICIT NONE
    type(t_kpts),intent(in)  :: kpts
    type(t_cell),intent(in)  :: cell
    INTEGER, INTENT (IN)     :: jspins
    type(t_banddos),INTENT(IN)   :: banddos
    CHARACTER(LEN=*), INTENT(IN) :: title
    CHARACTER(len=1) :: ssy
    type(t_cell) :: p_cell
    real   :: vkr(3),vkr_prev(3)

    INTEGER n,aoff,adel,k

      INTEGER numSpecialPoints, addFront, addBack

      REAL,              ALLOCATABLE :: d(:)
      REAL,              ALLOCATABLE :: specialPoints(:,:)
      CHARACTER(LEN=50), ALLOCATABLE :: specialPointNames(:)


      numSpecialPoints = kpts%numSpecialPoints
      addFront = 0
      addBack = 0
      
      IF (kpts%numSpecialPoints.EQ.0) THEN
         numSpecialPoints = 2
         addFront = 1
         addBack = 1
      ELSE
         IF (kpts%specialPointIndices(1).NE.1) THEN
            addFront = 1
            numSpecialPoints = numSpecialPoints + 1
         END IF
         IF (kpts%specialPointIndices(kpts%numSpecialPoints).NE.kpts%nkpt) THEN
            addBack = 1
            numSpecialPoints = numSpecialPoints + 1
         END IF
      END IF

      ALLOCATE (specialPoints(3,numSpecialPoints))
      ALLOCATE (specialPointNames(numSpecialPoints))
      ALLOCATE (d(numSpecialPoints))
      IF(addFront.EQ.1) THEN
         specialPoints(:,1) = kpts%bk(:,1)
         specialPointNames(1) = '-'
      END IF
      IF(ALLOCATED(kpts%specialPoints)) THEN
         specialPoints(:,1+addFront:kpts%numSpecialPoints+addFront) = kpts%specialPoints(:,1:kpts%numSpecialPoints)
         specialPointNames(1+addFront:kpts%numSpecialPoints+addFront) = kpts%specialPointNames(1:kpts%numSpecialPoints)
      END IF
      IF(addBack.EQ.1) THEN
         specialPoints(:,kpts%numSpecialPoints+addFront+1) = kpts%bk(:,kpts%nkpt)
         specialPointNames(kpts%numSpecialPoints+addFront+1) = '-'
      END IF

    aoff = iachar('a')-1
    adel = iachar('a')-iachar('A')
    !write(*,*) aoff,adel
    
    p_cell=cell
    DO k =1,3
      p_cell%amat(1,k)=cell%amat(1,k)/banddos%s_cell_x
      p_cell%amat(2,k)=cell%amat(2,k)/banddos%s_cell_y
      p_cell%amat(3,k)=cell%amat(3,k)/banddos%s_cell_z
    END DO
      CALL inv3(p_cell%amat,p_cell%bmat,p_cell%omtil)
      p_cell%bmat=p_cell%bmat*tpi_const


    !Generate distances
    d(1)=0.0
    vkr_prev=matmul(specialPoints(:,1),p_cell%bmat)
    DO k=2,numSpecialPoints
      vkr=matmul(specialPoints(:,k),p_cell%bmat)
      d(k)=d(k-1)+sqrt(dot_product(vkr-vkr_prev,vkr-vkr_prev))
      vkr_prev=vkr
    ENDDO


    OPEN (27,file='band_sc.gnu',status='unknown')
    WRITE (27,900)
    WRITE (27,901)
    WRITE (27,902)
    WRITE (27,903)
    WRITE (27,'(3a)') 'set title "', TRIM(ADJUSTL(title)), '"'
    DO n = 1, numSpecialPoints
      WRITE (27,905) d(n),d(n)
    ENDDO
    WRITE (27,906) d(1),d(size(d))
!
! nomal labels
!
    ssy=specialPointNames(1)
    IF (iachar(ssy) < aoff ) THEN
      WRITE (27,907) ssy,d(1)
    ELSE
     WRITE (27,907) " ",d(1)
   ENDIF
    DO n = 2, numSpecialPoints-1
      ssy=specialPointNames(n)
      IF (iachar(ssy) < aoff ) THEN
        WRITE (27,908) ssy,d(n)
      ELSE
        WRITE (27,908) " ",d(n)
      ENDIF
    ENDDO
    ssy=specialPointNames(n)
    IF (iachar(ssy) < aoff ) THEN
      WRITE (27,909) ssy,d(n)
    ELSE
      WRITE (27,909) " ",d(n)
    ENDIF

    DO n=1,numSpecialPoints
      ssy=specialPointNames(n)
      IF (iachar(ssy) > aoff ) THEN
        WRITE (27,914) achar(iachar(ssy)-adel),d(n)
      ENDIF
    ENDDO

!
! now write the rest
    WRITE (27,910)
    WRITE (27,*) 'set palette model RGB'
    WRITE (27,*) 'set palette defined (-2 "black", -1 "white" ,0 "white",',achar(92)
    WRITE (27,*) '0.67 "light-blue",1 "blue")'
    WRITE (27,*) 'set cbrange [-2:1]'
    WRITE (27,*) 'unset colorbox'
    WRITE (27,*) 'size1(x)=0.9*x**(0.4)'
    WRITE (27,*) 'color1(x)=0.3+x/2.4'
    WRITE (27,*) 'size2(x)=0.35*(1-x**(0.01))'
    WRITE (27,*) 'color2(x)=1.15*(x-1)'
    WRITE (27,*) 'e_f=0.000000 #fermi energy is already corrected when using hdf5'
    WRITE (27,911) d(size(d))+0.00001
    IF (jspins == 2) THEN
       WRITE (27,912)
       WRITE (27,916)
    END IF
    WRITE (27,913)
    WRITE (27,915)
    CLOSE (27)

      DEALLOCATE (d,specialPoints,specialPointNames)

900  FORMAT ('set terminal postscript enhanced color "Times-Roman" 20')
901  FORMAT ('set xlabel ""')
902  FORMAT ('set ylabel "E - E_F (eV)"')
903  FORMAT ('set nokey')
905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
#if (defined(CPP_AIX) || defined(__PGI))
907  FORMAT ('set xtics ("',a1,'"',f9.5,', \\')
908  FORMAT ('           "',a1,'"',f9.5,', \\')
#else
907  FORMAT ('set xtics ("',a1,'"',f9.5,', \')
908  FORMAT ('           "',a1,'"',f9.5,', \')
#endif
909  FORMAT ('           "',a1,'"',f9.5,'  )')
910  FORMAT ('set ytics -8,2,4')
#if (defined(CPP_AIX) || defined(__PGI))
911  FORMAT ('plot [0:',f9.5,'] [-9:5] \\')
912  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, \\')
916  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette,\\')
#else
911  FORMAT ('plot [0:',f9.5,'] [-9:5] \')
912  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, \')
916  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette,\')
#endif
#if (defined(CPP_AIX) || defined(__PGI))
913  FORMAT ('"bands_sc.1" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, \\')
#else
913  FORMAT ('"bands_sc.1" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, \')
#endif
915  FORMAT ('"bands_sc.1" using 1:($2-e_f):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette')
914  FORMAT ('set label "',a1,'" at ',f9.5,', -9.65 center font "Symbol,20"')

  !call system("gnuplot <band_sc.gnu >band_sc.ps")
  END SUBROUTINE write_gnu_sc
end MODULE
