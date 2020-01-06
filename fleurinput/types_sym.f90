 !--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_sym
  USE m_juDFT
  USE m_types_fleurinput_base
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sym
  !> Type to contain the symmetry information
  TYPE,EXTENDS(t_fleurinput_base):: t_sym
     !No of sym ops
     INTEGER :: nop = -1
     !Rot-matrices (3,3,nop)
     INTEGER, ALLOCATABLE::mrot(:, :, :)
     !translation vectors (3,nop)
     REAL, ALLOCATABLE::tau(:, :)
     !Symophic group
     LOGICAL ::symor=.false.
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
        !no of op that maps atom into
      INTEGER, ALLOCATABLE::ngopr(:)
      !symetry of atoms(nat)
      INTEGER, ALLOCATABLE::ntypsy(:)
      !atom mapped to by inversion (nat
      INTEGER, ALLOCATABLE ::invsat(:)

     !
     ! hybinp specific stuff TODO
     !
     INTEGER ::nsymt
     INTEGER :: nsym
   CONTAINS
     PROCEDURE :: init
     PROCEDURE :: print_xml
     PROCEDURE :: closure
     PROCEDURE :: read_xml=>read_xml_sym
     PROCEDURE :: mpi_bc => mpi_bc_sym
     PROCEDURE,PRIVATE :: check_close
  END TYPE t_sym
CONTAINS

    subroutine mpi_bc_sym(this,mpi_comm,irank)
      use m_mpi_bc_tool
      class(t_sym),INTENT(INOUT)::this
      integer,INTENT(IN):: mpi_comm
      INTEGER,INTENT(IN),OPTIONAL::irank
      INTEGER ::rank
      if (present(irank)) THEN
         rank=irank
      else
         rank=0
      end if

      call mpi_bc(this%nop,rank,mpi_comm)
      call mpi_bc(this%mrot,rank,mpi_comm)
      call mpi_bc(this%tau,rank,mpi_comm)
      call mpi_bc(this%symor,rank,mpi_comm)
      call mpi_bc(this%invs2,rank,mpi_comm)
      call mpi_bc(this%invs,rank,mpi_comm)
      call mpi_bc(this%zrfs,rank,mpi_comm)
      call mpi_bc(this%invtab,rank,mpi_comm)
      call mpi_bc(this%multab,rank,mpi_comm)
      call mpi_bc(this%nop2,rank,mpi_comm)
      call mpi_bc(this%d_wgn,rank,mpi_comm)
      call mpi_bc(this%invsatnr,rank,mpi_comm)
      call mpi_bc(this%invarop,rank,mpi_comm)
      call mpi_bc(this%invarind,rank,mpi_comm)
      call mpi_bc(this%nsymt,rank,mpi_comm)
      call mpi_bc(this%nsym,rank,mpi_comm)
      call mpi_bc(this%invsat,rank,mpi_comm)
      CALL mpi_bc(this%ngopr,rank,mpi_comm)
     call mpi_bc(this%ntypsy,rank,mpi_comm)
    end subroutine mpi_bc_sym

  SUBROUTINE read_xml_sym(this,xml)
    USE m_types_xml
    USE m_calculator
    CLASS(t_sym),INTENT(inout):: this
    TYPE(t_xml),INTENT(IN)  :: xml

    INTEGER:: number_sets,n
    CHARACTER(len=200)::str,path,path2



    this%nop = xml%GetNumberOfNodes('/fleurInput/calculationSetup/symmetryOperations/symOp')
    this%nop2=this%nop !might be changed later in film case
    ALLOCATE(this%mrot(3,3,this%nop))
    ALLOCATE(this%tau(3,this%nop))

    IF (this%nop<1) CALL judft_error("No symmetries in inp.xml")

    DO n=1,this%nop
       WRITE(path,"(a,i0,a)") '/fleurInput/calculationSetup/symmetryOperations/symOp[',n,']'
       str=xml%GetAttributeValue(TRIM(path)//'/row-1')
       READ(str,*) this%mrot(1,:,n),this%tau(1,n)
       str=xml%GetAttributeValue(TRIM(path)//'/row-2')
       READ(str,*) this%mrot(2,:,n),this%tau(2,n)
       str=xml%GetAttributeValue(TRIM(path)//'/row-3')
       READ(str,*) this%mrot(3,:,n),this%tau(3,n)
    ENDDO
  END SUBROUTINE read_xml_sym


  SUBROUTINE print_xml(sym,fh,filename)
    CLASS(t_sym),INTENT(IN)   :: sym
    INTEGER,INTENT(in)        ::fh
    CHARACTER(len=*),INTENT(in),OPTIONAL::filename

    INTEGER::i

    IF (PRESENT(filename)) OPEN(fh,file=filename,status='replace',action='write')

    WRITE(fh,'(a)') '      <symmetryOperations>'
    DO i = 1, sym%nop
       WRITE(fh,'(a)') '         <symOp>'
224    FORMAT('            <row-1>',i0,' ',i0,' ',i0,' ',f0.10,'</row-1>')
       WRITE(fh,224) sym%mrot(1,1,i), sym%mrot(1,2,i), sym%mrot(1,3,i), sym%tau(1,i)
225    FORMAT('            <row-2>',i0,' ',i0,' ',i0,' ',f0.10,'</row-2>')
       WRITE(fh,225) sym%mrot(2,1,i), sym%mrot(2,2,i), sym%mrot(2,3,i), sym%tau(2,i)
226    FORMAT('            <row-3>',i0,' ',i0,' ',i0,' ',f0.10,'</row-3>')
       WRITE(fh,226) sym%mrot(3,1,i), sym%mrot(3,2,i), sym%mrot(3,3,i), sym%tau(3,i)
       WRITE(fh,'(a)') '         </symOp>'
    END DO
    WRITE(fh,'(a)') '      </symmetryOperations>'
    IF (PRESENT(filename)) CLOSE(fh)
  END SUBROUTINE print_xml



  SUBROUTINE init(sym,cell,film)
    !Generates missing symmetry info.
    !tau,mrot and nop have to be specified alread
    USE m_types_cell
    use m_types_input
    CLASS(t_sym),INTENT(INOUT):: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    LOGICAL,INTENT(IN)        :: film



    INTEGER :: invsop, zrfsop, invs2op, magicinv,n,i,j,nn
    INTEGER :: optype(sym%nop),indtwo(sym%nop),usedop(sym%nop)
    INTEGER :: mrotaux(3,3,sym%nop)
    REAL    :: tauaux(3,sym%nop)
    LOGICAL :: zorth           ! true, if z-axis is othorgonal
    REAL,PARAMETER  :: eps12   = 1.0e-12

    IF (sym%nop==0) CALL judft_error("BUG in calling sym%init. mrot,tau&nop have to be set before")

    IF (ALLOCATED(sym%invtab)) DEALLOCATE(sym%invtab)
    IF (ALLOCATED(sym%multab)) DEALLOCATE(sym%multab)
    ALLOCATE ( sym%invtab(sym%nop),sym%multab(sym%nop,sym%nop) )
    CALL sym%check_close(optype)

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
          IF ( ALL ( ABS( sym%tau(:,n) ) < eps12 ) ) THEN
             invsop = n
             sym%invs = .TRUE.
          ENDIF
       ENDIF
       IF ( optype(n) == -2 )THEN
          !---> check for z-reflection
          IF ( sym%mrot(3,3,n) == -1 .AND. ALL(ABS(sym%tau(:,n))<eps12) ) THEN
             zrfsop = n
             sym%zrfs = .TRUE.
          ENDIF
       ENDIF
       IF ( optype(n) == 2 )THEN
          !---> check for 2d inversion
          IF ( sym%mrot(3,3,n) == 1 .AND. ALL(ABS(sym%tau(:,n))<eps12) ) THEN
             invs2op = n
             sym%invs2 = .TRUE.
          ENDIF
       ENDIF
    ENDDO !nops

    !if z-axis is not orthogonal we will not use z-reflect and 2d-invs
    IF ( cell%amat(3,1)==0.00 .AND. cell%amat(3,2)==0.00 .AND.cell%amat(1,3)==0.00 .AND.cell%amat(2,3)==0.00 ) THEN
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
             sym%mrot(:,:,i+sym%nop2) =  mrotaux(:,:,j)
             sym%tau(:,i+sym%nop2) = tauaux(:,j)
             usedop(j) = usedop(j) - 1
          ENDDO
          IF ( ANY( usedop(1:sym%nop) < 0 ) )  CALL juDFT_error("Fatal Error! #01",calledby="types_sym")

          IF ( 2*sym%nop2.NE.sym%nop ) THEN
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
             WRITE(6,'(/," Full space group has",i3," operations.",/)') sym%nop
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
    ENDIF


    !---> redo to ensure proper mult. table and mapping functions
    CALL sym%check_close(optype)


  END SUBROUTINE init

  FUNCTION closure(sym)RESULT(lclose)
    CLASS(t_sym),INTENT(IN):: sym
    LOGICAL                :: lclose

   !INTEGER, INTENT (IN)  :: mops           ! number of operations of the bravais lattice
   !INTEGER, INTENT (IN)  :: sym%nop           ! number of operations in space group
   !INTEGER, INTENT (IN)  :: sym%mrot(3,3,mops) ! refer to the operations of the
   !REAL,    INTENT (IN)  :: sym%tau(3,mops)    ! bravais lattice
    !LOGICAL, INTENT (OUT) :: lclose

   REAL    ttau(3),eps7
   INTEGER i,j,k,mp(3,3),map(sym%nop)

   eps7 = 1.0e-7

   ! loop over all operations
   DO j = 1, sym%nop

      map(1:sym%nop) = 0

      ! multiply {R_j|t_j}{R_i|t_i}
      DO i = 1, sym%nop
         mp = matmul( sym%mrot(:,:,j) , sym%mrot(:,:,i) )
         ttau = sym%tau(:,j) + matmul( sym%mrot(:,:,j) , sym%tau(:,i) )
         ttau = ttau - anint( ttau - eps7 )

         ! determine which operation this is
         DO k=1,sym%nop
            IF ( all( mp(:,:) == sym%mrot(:,:,k) ) .AND. all( abs( ttau(:)-sym%tau(:,k) ) < eps7 ) ) THEN
               IF ( map(i) .eq. 0 ) THEN
                  map(i) = k
               ELSE
                  write(6,*)'ERROR Closure: Multiplying ', j,' with ',k, ' and with ',map(i)
                  write(6,*) 'yields the same matrix'
                  lclose = .false.
                  RETURN
               END IF
            END IF
         END DO

         IF (map(i).eq.0) THEN
            write(6,*)'ERROR Closure:',i,' times',j,' leaves group'
            lclose = .false.
            RETURN
         END IF
      END DO
   END DO

   lclose = .true.

 END FUNCTION closure

SUBROUTINE check_close(sym,optype)
  CLASS(t_sym),INTENT(inout)::sym
  INTEGER, INTENT (OUT) :: optype(sym%nop)

   REAL    ttau(3)
   INTEGER i,j,n,k,mp(3,3),mdet,mtr

   REAL,    PARAMETER :: eps=1.0e-7
   INTEGER, PARAMETER :: cops(-1:3)=(/ 2, 3, 4, 6, 1 /)

   sym%invtab(1:sym%nop) = 0

   sym%multab = 0

   ! loop over all operations
   DO j = 1, sym%nop

      ! multiply {R_j|t_j}{R_i|t_i}
      DO i = 1, sym%nop
         mp = matmul( sym%mrot(:,:,j) , sym%mrot(:,:,i) )
         ttau = sym%tau(:,j) + matmul( sym%mrot(:,:,j) , sym%tau(:,i) )
         ttau = ttau - anint( ttau - eps )

         ! determine which operation this is
         DO k=1,sym%nop
            IF ( all( mp(:,:) == sym%mrot(:,:,k) ) .and. all( abs( ttau(:)-sym%tau(:,k) ) < eps ) ) THEN
               IF ( sym%multab(j,i) .eq. 0 ) THEN
                  sym%multab(j,i) = k
                  IF (k .eq. 1) sym%invtab(j)=i
               ELSE
                  WRITE(6,'(" Symmetry error: multiple ops")')
                  CALL juDFT_error("check_close: Multiple ops",calledby ="closure")
               END IF
            END IF
         END DO

         IF (sym%multab(j,i).eq.0) THEN
            WRITE (6,'(" Group not closed")')
            WRITE (6,'("  j , i =",2i4)') j,i
            CALL juDFT_error("check_close: Not closed",calledby="closure")
         END IF
      END DO
   END DO

   ! determine the type of each operation
   DO n = 1, sym%nop
      mtr = sym%mrot(1,1,n) + sym%mrot(2,2,n) + sym%mrot(3,3,n)
      mdet = sym%mrot(1,1,n)*(sym%mrot(2,2,n)*sym%mrot(3,3,n)-sym%mrot(3,2,n)*sym%mrot(2,3,n)) +&
             sym%mrot(1,2,n)*(sym%mrot(3,1,n)*sym%mrot(2,3,n)-sym%mrot(2,1,n)*sym%mrot(3,3,n)) +&
             sym%mrot(1,3,n)*(sym%mrot(2,1,n)*sym%mrot(3,2,n)-sym%mrot(3,1,n)*sym%mrot(2,2,n))

      optype(n) = mdet*cops(mdet*mtr)

   END DO

END SUBROUTINE check_close
END MODULE m_types_sym
