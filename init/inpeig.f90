      MODULE m_inpeig
      CONTAINS
      SUBROUTINE inpeig(&
     &                  atoms,cell,input,l_is_oneD,kpts,enpara,kptsFilename,latnam)
!*********************************************************************
!     inputs the necessary quantities for the eigenvalue part (energy
!     parameters, k-points, wavefunction cutoffs, etc.).
!                  m. weinert   jan. 1987
!     modification dec. 1990:
!     dummyline before reading l-dependent energies to make reading of
!     input easier (e.g. insert name of atom etc.)
!     modification dec. 93:
!     for step-forward diagonalization a la wu in case of more
!     than 1 window we read now
!     number of occupied states for EACH window
!*********************************************************************
      USE m_gkptwgt
      USE m_constants
      USE m_types_atoms
      USE m_types_cell
      USE m_types_input
      USE m_types_kpts
      USE m_types_enpara
      USE m_juDFT
     
      IMPLICIT NONE
!     ..
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_input),INTENT(IN)     :: input
      LOGICAL,INTENT(IN)           :: l_is_oneD
      TYPE(t_kpts),INTENT(INOUT)   :: kpts
      TYPE(t_enpara),OPTIONAL,INTENT(INOUT) :: enpara
      
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: kptsFilename
      CHARACTER(len=*),INTENT(IN)    :: latnam

!     ..
!     .. Local Scalars ..
      REAL      :: wt,scale
      INTEGER   :: i,j,nk,jsp,n
      LOGICAL   :: xyu,l_enpara,l_clf,l_k
      CHARACTER(LEN=255) :: fname
!     ..
!
     
!---> input energy parameters for each atom.
!---> the energy parameters for l.ge.3 have the same value
!---> read from file 40='enpara'  shz Jan.96
!

      IF(PRESENT(enpara)) THEN
         IF (.NOT.input%l_inpXML) THEN
            !read enpara file if present!
            CALL enpara%init_enpara(atoms,input%jspins,input%film)
         END IF
      END IF
!
!---> read k-points from file 41='kpts'
!

      IF(PRESENT(kptsFilename)) THEN
         fname = TRIM(ADJUSTL(kptsFilename))
      ELSE
         fname = 'kpts'
      END IF

      INQUIRE(file=TRIM(ADJUSTL(fname)),exist=l_k)
      if (.not.l_k) return 

      OPEN (41,file=TRIM(ADJUSTL(fname)),form='formatted',status='old')
!
!---> k-mesh: given in units of the reciprocal lattice basis vectors
!---> scale is a factor to make input easier (default=1.0). k-pt
!---> weights can be relative weights since they are renormalized.
!---> input: for bulk - k1,k2,k3,wtkpt
!--->        for film - k1,k2,wtkpt
!--->           weights are calculated for films, if wtkpt=0
!     for film calculation k1,k2 may also be read in xy - units : xyu=T
!     1 = boundery of BZ on kx/ky axis
!                                                  shz Feb.96
         READ (41,FMT=8110,ERR=911,END=911) kpts%nkpt,scale,xyu
         GOTO 912
  911    CONTINUE
         xyu = .false.
  912    CONTINUE
         
         IF (kpts%nkpt.GT.kpts%nkpt)  THEN
           CALL juDFT_error('nkptd too small',calledby='inpeig')
         ENDIF
 8100    FORMAT (i5,f20.10)
 8110    FORMAT (i5,f20.10,3x,l1)
         IF (scale.EQ.0.0) scale = 1.0
         DO nk = 1,kpts%nkpt
            READ (41,FMT=8040) (kpts%bk(i,nk),i=1,3),kpts%wtkpt(nk)
 8040       FORMAT (4f10.5)
            IF (input%film .AND. .NOT.l_is_oneD) THEN
               kpts%wtkpt(nk) = kpts%bk(3,nk)
               kpts%bk(3,nk) = 0.0
               IF (xyu) THEN
!           transform to cartesian coordinates
                  IF (latnam.EQ.'hex') THEN
                     kpts%bk(1,nk) = kpts%bk(1,nk)*tpi_const/cell%amat(2,2)
                     kpts%bk(2,nk) = kpts%bk(2,nk)*pi_const/cell%amat(1,1)
                  ELSE
                     kpts%bk(1,nk) = kpts%bk(1,nk)*pi_const/cell%amat(1,1)
                     kpts%bk(2,nk) = kpts%bk(2,nk)*pi_const/cell%amat(2,2)
                  END IF
!           transform to internal coordinates
                  kpts%bk(1:2,nk)=matmul(kpts%bk(1:2,nk),cell%amat(1:2,1:2))/tpi_const
               END IF
            ELSEIF (.NOT.input%film) THEN
              IF (xyu) THEN
                call juDFT_warn("The xyu feature is not tested",calledby="inpeig")
                kpts%bk(:,nk) = kpts%bk(:,nk)!*2.0/sc !TODO what is this scaling?
                kpts%bk(:,nk) = matmul( cell%amat, kpts%bk(:,nk) )
              ENDIF
            END IF
            DO  i = 1,3
               kpts%bk(i,nk) = kpts%bk(i,nk)/scale
            ENDDO
!-odim
            IF (l_is_oneD) THEN
!--> trapezoidal
               IF (kpts%bk(3,nk).EQ.0. .OR. kpts%bk(3,nk).EQ.0.5) THEN
                  kpts%wtkpt(nk) = 1.
               ELSE
                  kpts%wtkpt(nk) = 2.
               END IF
            END IF
!-odim
         ENDDO
         wt = sum(kpts%wtkpt(:kpts%nkpt))

         IF (wt.EQ.0.0) THEN
            IF (input%film) THEN
!
!---> generate k-point weights for 2d BZ: squ, rec, cen, hex
!     determine new wt
!
               CALL gkptwgt(&
     &                      kpts,cell,latnam)
               wt=sum(kpts%wtkpt)
             ELSE
               CALL juDFT_error("wtkpts",calledby ="inpeig",hint&
     &              ="The sum of weights in the kpts file is zero")
            END IF
         END IF
         IF (l_is_oneD)  kpts%bk(1:2,:) = 0.0    
         kpts%wtkpt(:) = kpts%wtkpt(:)/wt
     
         WRITE (6,FMT=8120)  kpts%nkpt
         DO  nk = 1,kpts%nkpt
            WRITE (6,FMT=8040)  kpts%bk(:,nk),kpts%wtkpt(nk)
         ENDDO
 8120    FORMAT (1x,/,' number of k-points for this window =',i5,/,t12,&
     &          'coordinates',t34,'weights')
      CLOSE (41)

      END SUBROUTINE inpeig
      END MODULE m_inpeig
