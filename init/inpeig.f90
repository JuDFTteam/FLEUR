      MODULE m_inpeig
      CONTAINS
      SUBROUTINE inpeig(&
     &                  atoms,cell,input,l_is_oneD,kpts,enpara)
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
      USE m_enpara,    ONLY : r_enpara,default_enpara
      USE m_types
      use m_juDFT

      IMPLICIT NONE
!     ..
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_input),INTENT(IN)     :: input
      LOGICAL,INTENT(IN)           :: l_is_oneD
      TYPE(t_kpts),INTENT(INOUT)   :: kpts
      TYPE(t_enpara),INTENT(INOUT) :: enpara

!     ..
!     .. Local Scalars ..
      REAL      :: wt,scale
      INTEGER   :: i,j,nk,jsp,n
      LOGICAL   :: xyu,l_enpara
!     ..
!
     
!---> input energy parameters for each atom.
!---> the energy parameters for l.ge.3 have the same value
!---> read from file 40='enpara'  shz Jan.96
!
      l_enpara = .FALSE.
      INQUIRE (file ='enpara',exist= l_enpara)
      IF (l_enpara) THEN
         OPEN (40,file ='enpara',form='formatted',status='old')
           DO jsp = 1,input%jspins
            CALL r_enpara(&
     &                    atoms,input,jsp,enpara)
           ENDDO !dimension%jspd
         CLOSE (40)
      ELSE IF (.NOT.input%l_inpXML) THEN
         WRITE(6,*) "No enpara file found, using default values"
         enpara%el0(:,:,1)=0.0
         enpara%el0(0,:,1)=-999999.0
         DO n = 1, atoms%ntype
            enpara%skiplo(n,:) = 0
            DO i = 1, atoms%nlo(n)
               enpara%skiplo(n,:) = enpara%skiplo(n,1) + (2*atoms%llo(i,n)+1)
            END DO
         END DO
         CALL default_enpara(1,atoms,enpara)
         IF (input%jspins>1) THEN
           enpara%el0(:,:,2)=enpara%el0(:,:,1)
           enpara%ello0(:,:,2)=enpara%ello0(:,:,1)
         ENDIF
         IF (input%film) enpara%evac0 = eVac0Default_const
      END IF
!
!---> read k-points from file 41='kpts'
!
      OPEN (41,file='kpts',form='formatted',status='old')
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
         
         IF (kpts%nkpt.GT.kpts%nkptd)  THEN
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
                  IF (cell%latnam.EQ.'hex') THEN
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
     &                      kpts,cell)
               wt=sum(kpts%weight)
             ELSE
               CALL juDFT_error("wtkpts",calledby ="inpeig",hint&
     &              ="The sum of weights in the kpts file is zero")
            END IF
         END IF
         IF (l_is_oneD)  kpts%bk(1:2,:) = 0.0    
         kpts%wtkpt(:) = kpts%wtkpt(:)/wt
     
         WRITE (6,FMT=8120)  kpts%nkpt
         WRITE (16,FMT=8120) kpts%nkpt
         DO  nk = 1,kpts%nkpt
            WRITE (6,FMT=8040)  kpts%bk(:,nk),kpts%wtkpt(nk)
            WRITE (16,FMT=8040) kpts%bk(:,nk),kpts%wtkpt(nk)
         ENDDO
 8120    FORMAT (1x,/,' number of k-points for this window =',i5,/,t12,&
     &          'coordinates',t34,'weights')
      CLOSE (41)

      END SUBROUTINE inpeig
      END MODULE m_inpeig
