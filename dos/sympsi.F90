!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_sympsi

  ! Calculates the irreducible represetantions of the wave functions.
  ! if k-point is in Brillouin zone boundary results are correct only for
  ! non-symmorphic groups (factor groups would be needed for that...). 
  ! jsym contains the number of irreducible rep., corresponding character
  ! tables are given in the file syminfo.
  !
  ! Double groups work only with non-collinear calculations, for normal spin-orbit 
  ! calculations both spin up and down components would be needed...

  ! Jussi Enkovaara, Juelich 2004

CONTAINS
  SUBROUTINE sympsi(bkpt,nv,kx,ky,kz,sym,DIMENSION,ne,cell,eig,noco, ksym,jsym,z_r,z_c,l_real)

    USE m_grp_k
    USE m_inv3
    USE m_types
    IMPLICIT NONE

    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    !
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nv,ne
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: kx(:),ky(:),kz(:)!(nvd) 
    REAL,    INTENT (IN) :: bkpt(3),eig(DIMENSION%neigd)
    REAL,    INTENT (IN) :: z_r(DIMENSION%nbasfcn,DIMENSION%neigd)
    COMPLEX, INTENT (IN) :: z_c(DIMENSION%nbasfcn,DIMENSION%neigd)
    INTEGER, INTENT (OUT):: jsym(DIMENSION%neigd),ksym(DIMENSION%neigd)
    LOGICAL,INTENT(IN)   :: l_real
    !     ..
    !     .. Local Scalars ..
    REAL degthre
    INTEGER i,k,n,c
    INTEGER nclass,nirr,n1,n2 ,ndeg
    LOGICAL soc, char_written
    !     ..
    !     .. Local Arrays ..
    INTEGER mrot_k(3,3,2*sym%nop)
    INTEGER :: mtmpinv(3,3),d
    INTEGER :: gmap(DIMENSION%nvd,sym%nop)
    REAL ::    kv(3),kvtest(3)
    INTEGER :: deg(ne)

    REAL :: norm(ne)
    LOGICAL :: symdone(ne)

    COMPLEX, ALLOCATABLE :: csum(:,:,:),chars(:,:)
    COMPLEX, SAVE, ALLOCATABLE :: char_table(:,:)
    CHARACTER(LEN=7) :: grpname
    CHARACTER(LEN=5) :: irrname(2*sym%nop)
    COMPLEX          :: c_table(2*sym%nop,2*sym%nop)
    COMPLEX, ALLOCATABLE :: su(:,:,:)
    !
    REAL,PARAMETER:: small=1.0e-4

    soc=noco%l_soc.AND.noco%l_noco
    jsym=0
    ksym=0
    IF (noco%l_soc.AND.(.NOT.noco%l_noco)) RETURN

    IF (soc) THEN
       ALLOCATE(su(2,2,2*sym%nop))
       CALL grp_k(sym,mrot_k,cell,bkpt,nclass,nirr,c_table, grpname,irrname,su)
    ELSE
       CALL grp_k(sym,mrot_k,cell,bkpt,nclass,nirr,c_table, grpname,irrname)
    ENDIF
    ALLOCATE(csum(ne,ne,nclass))
    ALLOCATE(chars(ne,nclass))
    chars=0.0
    !>

    IF (ALLOCATED(char_table)) THEN
       IF (SIZE(char_table,2).NE.nclass) THEN
          DEALLOCATE(char_table)
          ALLOCATE(char_table(nirr,nclass))
          char_written=.FALSE.
       ENDIF
    ELSE
       ALLOCATE(char_table(nirr,nclass))
       char_written=.FALSE.
    ENDIF
    char_table(:,:) = c_table(1:nirr,1:nclass)

    !<--map the (k+g)-vectors related by inv(rot)
    gmap=0
    DO c=1,nclass
       CALL inv3(mrot_k(:,:,c),mtmpinv,d)
       kloop: DO k=1,nv
          kv(1)=kx(k)
          kv(2)=ky(k)
          kv(3)=kz(k)
          kv=kv+bkpt
          kvtest=MATMUL(kv,mtmpinv)
          !         kvtest=MATMUL(kv,mrot_k(:,:,c))
          DO i = 1,nv
             kv(1)=kx(i)
             kv(2)=ky(i)
             kv(3)=kz(i)
             kv=kv+bkpt
             IF (ABS(kvtest(1)-kv(1)).LT.small.AND.&
                  ABS(kvtest(2)-kv(2)).LT.small.AND. ABS(kvtest(3)-kv(3)).LT.small) THEN
                gmap(k,c)=i
                CYCLE kloop
             ENDIF
          ENDDO
          WRITE(6,*) 'Problem in symcheck, cannot find rotated kv for', k,kx(k),ky(k),kz(k)
          RETURN
       ENDDO kloop
    ENDDO

    !norms
    DO i=1,ne
       norm(i)=0.0
       IF (soc) THEN
          DO k=1,nv*2
             norm(i)=norm(i)+ABS(z_c(k,i))**2
          ENDDO
       ELSE
          IF (l_real) THEN
             DO k=1,nv
                norm(i)=norm(i)+ABS(z_r(k,i))**2
             ENDDO
          ELSE
             DO k=1,nv
                norm(i)=norm(i)+ABS(z_c(k,i))**2
             ENDDO
          ENDIF
       ENDIF
       norm(i)=SQRT(norm(i))
    ENDDO


    !<-- Calculate the characters
    symdone=.FALSE.
    stateloop: DO i=1,ne
       IF (symdone(i)) CYCLE stateloop
       ndeg=0
       deg=0
       degthre=0.0001
       DO n=1,ne
          IF (ABS(eig(i)-eig(n)).LT.degthre) THEN
             ndeg=ndeg+1
             deg(ndeg)=n
          ENDIF
       ENDDO

       csum=0.0
       DO c=1,nclass
          DO n1=1,ndeg
             DO n2=1,ndeg
                IF (l_real) THEN
                   DO k=1,nv
                      csum(n1,n2,c)=csum(n1,n2,c)+z_r(k,deg(n1))*&
                           z_r(gmap(k,c),deg(n2))/(norm(deg(n1))*norm(deg(n2)))
                   END DO
                ELSE
                   IF (soc) THEN  
                      DO k=1,nv

                         csum(n1,n2,c)=csum(n1,n2,c)+(CONJG(z_c(k,deg(n1)))*&
                              (su(1,1,c)*z_c(gmap(k,c),deg(n2))+ su(1,2,c)*z_c(gmap(k,c)+nv,deg(n2)))+&
                              CONJG(z_c(k+nv,deg(n1)))* (su(2,1,c)*z_c(gmap(k,c),deg(n2))+&
                              su(2,2,c)*z_c(gmap(k,c)+nv,deg(n2))))/ (norm(deg(n1))*norm(deg(n2)))
                      END DO
                   ELSE
                      DO k=1,nv
                         csum(n1,n2,c)=csum(n1,n2,c)+CONJG(z_c(k,deg(n1)))*&
                              z_c(gmap(k,c),deg(n2))/(norm(deg(n1))*norm(deg(n2)))
                      END DO
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       ! We might have taken degenerate states which are not degenerate due to symmetry
       ! so look for irreducible reps
       DO n1=1,ndeg
          chars(deg(n1),:)=0.0
          DO n2=1,ndeg
             IF (ANY(ABS(csum(n1,n2,:)).GT.0.01)) THEN
                chars(deg(n1),:)=chars(deg(n1),:)+csum(n2,n2,:)
             ENDIF
          ENDDO
          symdone(deg(n1))=.TRUE.
       ENDDO


       ! determine the irreducible presentation
       irrloop: DO n1=1,ndeg 
          !        write(*,'(2i3,6(2f6.3,2x))') n1,i,chars(deg(n1),1:nclass)
          DO c=1,nirr
             IF (ALL(ABS(chars(deg(n1),1:nclass)-&
                  &             char_table(c,1:nclass)).LT.0.001)) THEN
                jsym(deg(n1))=c
                CYCLE irrloop
             ELSE IF (ALL(ABS(char_table(c,1:nclass)).LT.0.001)) THEN
                char_table(c,:)=chars(deg(n1),:)
                jsym(deg(n1))=c
                CYCLE irrloop
             ENDIF
          ENDDO
       ENDDO irrloop

    ENDDO stateloop
    !>

    IF (.NOT.char_written) THEN
       WRITE(444,124) bkpt
       WRITE(444,*) 'Group is ' ,grpname
       DO c=1,nirr
          IF (l_real)THEN
             IF (ANY(ABS(char_table).GT.0.001)) THEN
                WRITE(444,123) c,irrname(c),(char_table(c,n),n=1,nclass)
             ELSE
                WRITE(444,123) c,irrname(c),(REAL(char_table(c,n)),n=1,nclass)
             ENDIF
          ELSE
             IF (ANY(AIMAG(char_table).GT.0.001)) THEN
                WRITE(444,123) c,irrname(c),(char_table(c,n),n=1,nclass)
             ELSE
                WRITE(444,123) c,irrname(c),(REAL(char_table(c,n)),n=1,nclass)
             ENDIF
          ENDIF
       ENDDO
       char_written=.TRUE.
    ENDIF
123 FORMAT(i3,1x,a5,1x,20f7.3)
124 FORMAT('Character table for k: ',3f8.4)

    DEALLOCATE(csum)
    DEALLOCATE(chars)

  END SUBROUTINE sympsi

END  MODULE m_sympsi
