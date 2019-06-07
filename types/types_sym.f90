 !--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sym
  IMPLICIT NONE
  !symmetry information
   TYPE t_sym
      !No of sym ops
      INTEGER ::nop
      !Rot-matrices (3,3,nop)
      INTEGER, ALLOCATABLE::mrot(:, :, :)
      !translation vectors (3,nop)
      REAL, ALLOCATABLE::tau(:, :)
      !Symophic group
      LOGICAL ::symor
      !2D-inv-sym
      LOGICAL ::invs2
      !Inversion-sym
      LOGICAL ::invs
      !Z-refls. sym
      LOGICAL ::zrfs
      !inverse operation (nop)
      INTEGER, ALLOCATABLE::invtab(:)
      !multiplication table
      INTEGER, ALLOCATABLE :: multab(:,:) !(nop,nop)
      !No of 2D-sym ops
      INTEGER ::nop2
      !Wigner matrix for lda+u
      COMPLEX, ALLOCATABLE:: d_wgn(:, :, :, :)
      !
      ! Atom sepecific stuff
      !
      INTEGER, ALLOCATABLE :: invsatnr(:) !(atoms%nat)
      INTEGER, ALLOCATABLE :: invarop(:,:)!(atoms%nat,nop)
      INTEGER, ALLOCATABLE :: invarind(:) !(atoms%nat)

      !
      ! Hybrid specific stuff TODO
      !
      INTEGER ::nsymt
      INTEGER :: nsym
      !
      ! Description of initalization
      !
      INTEGER :: symSpecType
      !Name of lattice type
      CHARACTER*3   :: latnam
      !Name of sym
      CHARACTER*4   :: namgrp


    CONTAINS
      PROCEDURE :: init
   END TYPE t_sym
 CONTAINS
   SUBROUTINE init(sym,amat,film)
     !Generates missing symmetry info.
     !tau,mrot and nop have to be specified already
     USE m_closure
     USE m_dwigner
     CLASS(t_sym),INTENT(INOUT):: sym
     REAL,INTENT(in)           :: amat(3,3)
     LOGICAL,INTENT(in)        :: film


     INTEGER :: invsop, zrfsop, invs2op, magicinv,n,i,j,nn
     INTEGER :: optype(sym%nop),indtwo(sym%nop),usedop(sym%nop)
     INTEGER :: mrotaux(3,3,sym%nop)
     REAL    :: tauaux(3,sym%nop)
     LOGICAL :: zorth           ! true, if z-axis is othorgonal
     REAL,PARAMETER  :: eps12   = 1.0e-12

     IF (sym%nop==0) CALL judft_error("BUG in calling sym%init. mrot,tau&nop have to be set before")

     IF (ALLOCATED(sym%invtab)) DEALLOCATE(sym%intab)
     IF (ALLOCATED(sym%multab)) DEALLOCATE(sym%multab)
     ALLOCATE ( sym%invtab(sym%nop),sym%multab(sym%nop,sym%nop) )
     CALL check_close(sym%nop,sym%mrot,sym%tau, sym%multab,sym%invtab,optype)

     !---> determine properties of symmmetry operations, 
     ! Code previously in symproperties
     sym%symor=.NOT.(ANY(ABS(sym%tau(:,:sym%nop))>eps12))

     sym%invs    = .FALSE.
     sym%zrfs    = .FALSE.
     sym%invs2   = .FALSE.
     invsop  = 0
     zrfsop  = 0
     invs2op = 0
     
     DO n = 1, sym%nop
        IF ( optype(n) == -1 )THEN
            !--->  check if we have inversion as a symmetry operation
            IF ( ALL ( ABS( tau(:,n) ) < eps12 ) ) THEN
               invsop = n
               sym%invs = .TRUE.
            ENDIF
         ENDIF
         IF ( optype(n) == -2 )THEN
            !---> check for z-reflection
            IF ( mrot(3,3,n) == -1 .AND. ALL(ABS(tau(:,n))<eps12) ) THEN
               zrfsop = n
               sym%zrfs = .TRUE.
            ENDIF
         ENDIF
         IF ( optype(n) == 2 )THEN
            !---> check for 2d inversion
            IF ( mrot(3,3,n) == 1 .AND. ALL(ABS(tau(:,n))<eps12) ) THEN
               invs2op = n
               sym%invs2 = .TRUE.
            ENDIF
         ENDIF     
      ENDDO !nops

      !if z-axis is not orthogonal we will not use z-reflect and 2d-invs
      IF ( amat(3,1)==0.00 .AND. amat(3,2)==0.00 .AND.amat(1,3)==0.00 .AND. amat(2,3)==0.00 ) THEN
         zorth= .TRUE.
      ELSE       
         zorth= .FALSE.
         ! reset the following...
         sym%zrfs    = .FALSE.
         sym%invs2   = .FALSE.
      ENDIF
      IF (film.AND.((.NOT.zorth) )) &
           CALL juDFT_error("film = t and z-axis not orthogonal",calledby ="types_sym")

      IF (film) THEN
      !---> now we have to sort the ops to find the two-dimensional ops
      !---> and their 3-dim inverted or z-reflected counterparts
 
      mrotaux(:,:,1:sym%nop) = sym%mrot(:,:,1:sym%nop)
      tauaux(:,1:sym%nop) = sym%tau(:,1:sym%nop)

      DO i=1,sym%nop
         indtwo(i)= i
      ENDDO

      !Find number of pure 2D-operations
      sym%nop2=0
      DO i = 1, sym%nop
         IF ( sym%mrot(3,3,i) == 1 ) THEN
            sym%nop2 = sym%nop2 + 1
            indtwo(sym%nop2)= i
         ENDIF
      ENDDO


      magicinv = 0
      IF (sym%zrfs) magicinv = zrfsop
      IF (sym%invs) magicinv = invsop
      usedop = 1

      IF ( magicinv > 0 ) THEN
        DO i = 1, sym%nop2
          j = indtwo(i)
          sym%mrot(:,:,i) = mrotaux(:,:,j)
          sym%tau(:,i)    = tauaux(:,j)
          usedop(j) = usedop(j) - 1
          j = sym%multab(magicinv,indtwo(i))
          sym%mrot(:,:,i+nop2) =  mrotaux(:,:,j)
          sym%tau(:,i+nop2) = tauaux(:,j)
          usedop(j) = usedop(j) - 1
        ENDDO
        IF ( ANY( usedop(1:nops) < 0 ) )  CALL juDFT_error("Fatal Error! #01",calledby="types_sym")
 
        IF ( 2*sym%nop2.ne.sym%nop ) THEN
          n = 0
          DO i = 1, sym%nop
            IF ( usedop(i) == 1 ) THEN
              n = n + 1
              sym%mrot(:,:,2*sym%nop2+n) =  mrotaux(:,:,i)
              sym%tau(:,2*sym%nop2+n) = tauaux(:,i)
            ENDIF
          ENDDO
          IF ( n+2*sym%nop2 /= sym%nop )  CALL juDFT_error("Fatal Error! #02",calledby="types_sym")
        ENDIF

      ENDIF

!---> check for nonsymmorphic translations in z-direction in
!---> the film (oldfleur=t) case
      n = 1
      DO WHILE (n <= sym%nop)
         IF (ABS(sym%tau(3,n)) > 0.000001) THEN
            WRITE(6,'(/," Full space group has",i3," operations.",/)') nops
            WRITE(6,'(i3,"th operation violate the 2d symmetry in fleur and has been removed.",/)') n
             DO nn = n+1, sym%nop
                sym%mrot(:,:,nn-1) = sym%mrot(:,:,nn)
                sym%tau(:,nn-1) =sym%tau(:,nn)
             ENDDO
             sym%nop = sym%nop - 1
          ELSE
             n = n + 1
          ENDIF
       ENDDO
    ELSE !film
       sym%nop2 = 0
    ENDIF

    !Generated wigner symbols for LDA+U
    IF (ALLOCATED(sym%d_wgn)) DEALLOCATE(sym%d_wgn)
    ALLOCATE(sym%d_wgn(-3:3,-3:3,3,sym%nop))
    CALL d_wigner(sym%nop,sym%mrot,cell%bmat,3,sym%d_wgn)
    
    !---> redo to ensure proper mult. table and mapping functions
    CALL check_close(sym%nops,sym%mrot,sym%tau, sym%multab,sym%invtab,optype)
  END SUBROUTINE init
END MODULE m_types_sym
