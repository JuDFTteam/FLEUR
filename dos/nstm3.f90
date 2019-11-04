MODULE m_nstm3
  USE m_juDFT
  !***********************************************************************
  !     included writing to vacwave!
  !     set up mapping array to general G_parallel(j)=(gvac1(j),gvac2(j))
  !             for vacuum density in order to write out information
  !             on electronic structure for calculation of tunneling current    
  !                            change by shz, Jan.99
  !
  !***********************************************************************
CONTAINS
  SUBROUTINE nstm3(sym,atoms,vacuum,stars,lapw,ikpt,input,jspin,kpts,&
                   cell,evac,vz,gvac1d,gvac2d)

    USE m_sort
    USE m_types_setup

    USE m_types_lapw
    USE m_types_kpts
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_vacuum),INTENT(IN)   :: vacuum
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_kpts),INTENT(IN)     :: kpts
    TYPE(t_atoms),INTENT(IN)    :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: ikpt   
    INTEGER, INTENT (IN) :: jspin      
    !     ..
    !     .. Array  Arguments ..
    REAL,    INTENT (IN) :: evac(2)
    REAL,    INTENT (IN) :: vz(:,:)!(vacuum%nmzd,2)
    INTEGER, INTENT (OUT) :: gvac1d(:),gvac2d(:) !(lapw%dim_nv2d())
    !     ..
    !     .. Local Scalars
    INTEGER n2,k,j,i,ivac
    REAL    dz0
    !     ..
    !     .. Local Arrays ..
    INTEGER gvac1(SIZE(gvac1d)),gvac2(SIZE(gvac1d)),gindex(SIZE(gvac1d))
    REAL gvacl(SIZE(gvac1d)),gvac(2)
    !     ..
    !
    IF (ikpt.EQ.1) THEN
       n2 = 0
       k_loop: DO  k = 1,lapw%nv(jspin)
          DO j = 1,n2
             IF (lapw%k1(k,jspin).EQ.gvac1(j).AND.lapw%k2(k,jspin).EQ.gvac2(j)) THEN
                CYCLE k_loop
             END IF
          ENDDO
          n2 = n2 + 1
          gvac1(n2) = lapw%k1(k,jspin)
          gvac2(n2) = lapw%k2(k,jspin)
          DO i=1,2
             gvac(i)=lapw%k1(k,jspin)*cell%bmat(1,i)+lapw%k2(k,jspin)*cell%bmat(2,i)
          END DO
          gvacl(n2) = SQRT(REAL(gvac(1)**2+gvac(2)**2))
       ENDDO k_loop
       CALL sort(gindex(:n2),gvacl)
       DO j = 1,n2
          !  gvac1d, gvac2d are now ordered by increasing length
          gvac1d(j)=gvac1(gindex(j))
          gvac2d(j)=gvac2(gindex(j))
       END DO
       ! 
       IF (jspin.EQ.1) THEN
          WRITE (87,'(f10.6,1x,i1,1x,f10.6)') vacuum%tworkf,input%jspins,cell%area
          WRITE (87,'(2(f10.6,1x))') cell%amat(1,1), cell%amat(2,1)
          WRITE (87,'(2(f10.6,1x))') cell%amat(1,2), cell%amat(2,2)
          WRITE (87,'(2(f10.6,1x))') cell%bmat(1,1), cell%bmat(2,1)
          WRITE (87,'(2(f10.6,1x))') cell%bmat(1,2), cell%bmat(2,2)
          WRITE (87,'(i2)') sym%nop2
          DO j = 1, sym%nop2
             WRITE (87,'(i2,1x,i2)') sym%mrot(1,1,j), sym%mrot(1,2,j)
             WRITE (87,'(i2,1x,i2)') sym%mrot(2,1,j), sym%mrot(2,2,j)
          END DO
          WRITE (87,'(i3)') n2
          DO j = 1,n2
             WRITE (87,'(3(i3,1x),f10.6)') j, gvac1(gindex(j)), &
                  &              gvac2(gindex(j)),gvacl(gindex(j))
          END DO
          !
          !     Write info on 2D-starfunctions

          WRITE (87,'(i2,1x,i2,1x,i2)') stars%mx1,stars%mx2, stars%ng2
          DO i=1, stars%ng2
             WRITE (87,'(i2)') stars%nstr2(i)
          END DO
          DO i=-stars%mx1, stars%mx1
             DO j=-stars%mx2,stars%mx2
                WRITE (87,'(i2,1x,e12.4)') stars%ig2(stars%ig(i,j,0)),stars%rgphs(i,j,0)
             END DO
          END DO
       END IF
       WRITE (87,'(i1,1x,i1)') jspin, vacuum%nvac
       WRITE (87,'(2(e16.8,1x))') (evac(i), i=1,vacuum%nvac)
       WRITE (87,'(2(e16.8,1x))') (vz(vacuum%nmz,i), i=1,vacuum%nvac)
       dz0=0.0
       DO i=1, atoms%nat
          IF (ABS(atoms%taual(3,i)).GT.dz0) dz0=ABS(atoms%taual(3,i))
       END DO
       dz0=cell%z1-dz0*cell%amat(3,3)
       WRITE (87,'(i3,1x,f6.4,1x,f12.6)') vacuum%nmz,vacuum%delz,dz0   
       DO ivac=1,vacuum%nvac
          DO i=1, vacuum%nmz
             WRITE (87,'(e16.8)') vz(i,ivac)
          END DO
       END DO
       WRITE (87,'(i4)') kpts%nkpt
    END IF

    !  only write here if not on T3E


    WRITE (87,'(i3,1x,f12.6)') ikpt,kpts%wtkpt(ikpt)

  END SUBROUTINE nstm3
END MODULE m_nstm3
