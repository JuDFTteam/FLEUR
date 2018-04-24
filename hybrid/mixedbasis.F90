!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the mixed basis set used to evaluate the  !
! exchange term in HF/hybrid functional calculations or EXX           !
! calculations. In the latter case a second mixed basis set is setup  !
! for the OEP integral equation.                                      !
! In all cases the mixed basis consists of IR plane waves             !
!                                                                     !
! IR:                                                                 !
!    M_{\vec{k},\vec{G}} = 1/\sqrt{V} \exp{i(\vec{k}+\vec{G})}        !
!                                                                     !
! which are zero in the MT spheres and radial functions times         !
! spherical harmonics in the MT spheres                               !
!                                                                     !
! MT:                                                                 !
!     a                a                                              !
!    M              = U   Y                                           !
!     PL               PL  LM                                         !
!                                                                     !
!            where     a    a  a                                      !
!                     U  = u  u                                       !
!                      PL   pl  p'l'                                  !
!                                                                     !
!               and    L \in {|l-l'|,...,l+l'}                        !
!                                                                     !
!               and    P counts the different combinations of         !
!                      pl,p'l' which contribute to L                  !
!                                                                     !
!                                               M.Betzinger (09/07)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE m_mixedbasis
  use m_judft
CONTAINS

  SUBROUTINE mixedbasis( atoms,kpts,DIMENSION,input, &
       cell,sym,xcpot,hybrid, eig_id, mpi,v,l_restart)

    USE m_radfun,   ONLY : radfun
    USE m_radflo,   ONLY : radflo
    USE m_loddop,   ONLY : loddop
    USE m_util,     ONLY : intgrf_init,intgrf,rorderpf
    USE m_read_core
    USE m_wrapper
    USE m_eig66_io
    USE m_types
    IMPLICIT NONE

    TYPE(t_xcpot_inbuild),INTENT(IN)     :: xcpot
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_potden),INTENT(IN)    :: v

    ! - scalars -
    INTEGER,INTENT(IN)              :: eig_id
    LOGICAL,INTENT(INOUT)           :: l_restart



    ! - local scalars -
    INTEGER                         ::  ilo
    INTEGER                         ::  ikpt
    INTEGER                         ::  ispin,itype,l1,l2,l,n,igpt,n1,n2,nn,i,j,ic ,ng      
    INTEGER                         ::  jsp
    INTEGER                         ::  nodem,noded
    INTEGER                         ::  m,nk,ok
    INTEGER                         ::  x,y,z
    INTEGER                         ::  maxindxc,lmaxcd
    INTEGER                         ::  divconq ! use Divide & Conquer algorithm for array sorting (>0: yes, =0: no)
    REAL                            ::  gcutm
    REAL                            ::  wronk
    REAL                            ::  rdum,rdum1,rdum2

    LOGICAL                         ::  ldum,ldum1
    LOGICAL                         ::  exx


    ! - local arrays -
    INTEGER                         ::  g(3)
    INTEGER                         ::  lmaxc(atoms%ntype)
    INTEGER,ALLOCATABLE             ::  nindxc(:,:)
    INTEGER,ALLOCATABLE             ::  ihelp(:) 
    INTEGER,ALLOCATABLE             ::  ptr(:)           ! pointer for array sorting
    INTEGER,ALLOCATABLE             ::  unsrt_pgptm(:,:) ! unsorted pointers to g vectors

    REAL                            ::  kvec(3)
    REAL                            ::  flo(atoms%jmtd,2,atoms%nlod)
    TYPE(t_usdus)                   ::  usdus
    REAL                            ::  uuilon(atoms%nlod,atoms%ntype), duilon(atoms%nlod,atoms%ntype)
    REAL                            ::  ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
    REAL                            ::  potatom(atoms%jmtd,atoms%ntype)
    REAL                            ::  el(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd)
    REAL                            ::  ello(atoms%nlod,atoms%ntype,DIMENSION%jspd)
    REAL                            ::  bashlp(atoms%jmtd)

    REAL   ,ALLOCATABLE             ::  f(:,:,:),df(:,:,:)
    REAL   ,ALLOCATABLE             ::  olap(:,:),work(:),eig(:), eigv(:,:)      
    REAL   ,ALLOCATABLE             ::  bas1(:,:,:,:,:), bas2(:,:,:,:,:)
    REAL   ,ALLOCATABLE             ::  basmhlp(:,:,:,:)
    REAL   ,ALLOCATABLE             ::  gridf(:,:),vr0(:,:,:)
    REAL   ,ALLOCATABLE             ::  core1(:,:,:,:,:), core2(:,:,:,:,:)
    REAL   ,ALLOCATABLE             ::  eig_c(:,:,:,:)
    REAL   ,ALLOCATABLE             ::  length_kg(:,:) ! length of the vectors k + G

 
    LOGICAL,ALLOCATABLE             ::  selecmat(:,:,:,:)
    LOGICAL,ALLOCATABLE             ::  seleco(:,:),selecu(:,:)

    CHARACTER, PARAMETER            :: lchar(0:38) =&
         &      (/'s','p','d','f','g','h','i','j','k','l','m','n','o',&
         &        'x','x','x','x','x','x','x','x','x','x','x','x','x',&
         &        'x','x','x','x','x','x','x','x','x','x','x','x','x' /)
    
    CHARACTER(len=2)                ::  nchar
    CHARACTER(len=2)                ::  noel(atoms%ntype)
    CHARACTER(len=10)               ::  fname(atoms%ntype)
    ! writing to a file
    INTEGER, PARAMETER              ::  iounit = 125
    CHARACTER(10), PARAMETER        ::  ioname = 'mixbas'
    LOGICAL                         ::  l_found

    IF ( mpi%irank == 0 ) WRITE(6,'(//A,I2,A)') '### subroutine: mixedbasis ###'


    exx = xcpot%is_name("exx")
    if (exx) call judft_error("EXX is not implemented in this version",calledby='mixedbasis.F90')
    
    ! Deallocate arrays which might have been allocated in a previous run of this subroutine
    IF(ALLOCATED(hybrid%ngptm))    DEALLOCATE(hybrid%ngptm)
    IF(ALLOCATED(hybrid%ngptm1))   DEALLOCATE(hybrid%ngptm1)
    IF(ALLOCATED(hybrid%nindxm1))  DEALLOCATE(hybrid%nindxm1)
    IF(ALLOCATED(hybrid%pgptm))    DEALLOCATE(hybrid%pgptm)
    IF(ALLOCATED(hybrid%pgptm1))   DEALLOCATE(hybrid%pgptm1)
    IF(ALLOCATED(hybrid%gptm))     DEALLOCATE(hybrid%gptm)
    IF(ALLOCATED(hybrid%basm1) )   DEALLOCATE(hybrid%basm1)
   
    call usdus%init(atoms,dimension%jspd)

    ! If restart is specified read file if it already exists
    ! create it otherwise
    IF ( l_restart ) THEN

       ! Test if file exists
       INQUIRE(FILE=ioname,EXIST=l_found)

       IF ( l_found .and..false.) THEN !reading not working yet

          ! Open file
          OPEN(UNIT=iounit,FILE=ioname,FORM='unformatted',STATUS='old')

          ! Read array sizes
          !READ(iounit) kpts%nkptf,hybrid%gptmd
          READ(iounit) hybrid%maxgptm,hybrid%maxindx
          ! Allocate necessary array size and read arrays
          ALLOCATE ( hybrid%ngptm(kpts%nkptf),hybrid%gptm(3,hybrid%gptmd),hybrid%pgptm(hybrid%maxgptm,kpts%nkptf) )
          READ(iounit) hybrid%ngptm,hybrid%gptm,hybrid%pgptm,hybrid%nindx

          ! Read array sizes
          READ(iounit) hybrid%maxgptm1,hybrid%maxlcutm1,hybrid%maxindxm1,hybrid%maxindxp1
          ! Allocate necessary array size and read arrays
          ALLOCATE ( hybrid%ngptm1(kpts%nkptf),hybrid%pgptm1(hybrid%maxgptm1,kpts%nkptf),&
               &               hybrid%nindxm1(0:hybrid%maxlcutm1,atoms%ntype) )
          ALLOCATE ( hybrid%basm1(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
          READ(iounit) hybrid%ngptm1,hybrid%pgptm1,hybrid%nindxm1
          READ(iounit) hybrid%basm1

       
          CLOSE(iounit)

          RETURN
       END IF
    END IF

    hybrid%maxindx = 0
    DO itype = 1,atoms%ntype
       DO l = 0,atoms%lmax(itype)
          hybrid%maxindx = MAX(hybrid%maxindx,2+COUNT( atoms%llo(:atoms%nlo(itype),itype) .EQ. l))
       END DO
    END DO
    !       maxindx   = maxval( nlo   ) + 2



    ! initialize gridf for radial integration
    CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,gridf)

    !
    ! read in energy parameters from file eig
    ! to avoid meaningless energy parameters which occur in the case
    ! that the energy parameter is set to the atomic prinicipal 
    ! quantum number
    ! (el0 and ello0 are just the values in the enpara file)
    ! 


    DO jsp=1,DIMENSION%jspd
       CALL judft_error("TODO,mixedbasis")
!       CALL read_eig(eig_id,1,jsp,el=el(:,:,jsp),ello=ello(:,:,jsp))
    ENDDO

    ALLOCATE ( vr0(atoms%jmtd,atoms%ntype,DIMENSION%jspd) )

    vr0(:,:,:) = v%mt(:,0,:,:) 
   
    ! calculate radial basisfunctions u and u' with
    ! the spherical part of the potential vr0 and store them in
    ! bas1 = large component ,bas2 = small component

    ALLOCATE( f(atoms%jmtd,2,0:atoms%lmaxd), df(atoms%jmtd,2,0:atoms%lmaxd) )
    ALLOCATE( bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd) )
    ALLOCATE( bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd) )

    DO itype=1,atoms%ntype
       ng = atoms%jri(itype)
       DO ispin=1,DIMENSION%jspd
          DO l=0,atoms%lmax(itype)
             CALL radfun(l,itype,ispin,el(l,itype,ispin),vr0(:,itype,ispin),atoms,f(:,:,l),df(:,:,l),&
                  &                  usdus,nodem,noded,wronk)
          END DO
          bas1(1:ng,1,0:atoms%lmaxd,itype,ispin) =  f(1:ng,1,0:atoms%lmaxd)
          bas2(1:ng,1,0:atoms%lmaxd,itype,ispin) =  f(1:ng,2,0:atoms%lmaxd)
          bas1(1:ng,2,0:atoms%lmaxd,itype,ispin) = df(1:ng,1,0:atoms%lmaxd)
          bas2(1:ng,2,0:atoms%lmaxd,itype,ispin) = df(1:ng,2,0:atoms%lmaxd)

          hybrid%nindx(:,itype) = 2
          ! generate radial functions for local orbitals
          IF (atoms%nlo(itype).GE.1) THEN
             CALL radflo( atoms,itype,ispin,&
                  &                   ello(1,1,ispin),vr0(:,itype,ispin),&
                  &                   f,df,mpi,&
                  &                   usdus,&
                  &                   uuilon,duilon,ulouilopn,flo)

             DO ilo=1,atoms%nlo(itype)
                hybrid%nindx(atoms%llo(ilo,itype),itype) =&
                     &          hybrid%nindx(atoms%llo(ilo,itype),itype) + 1
                bas1(1:ng,hybrid%nindx(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),&
                     &                                    itype,ispin) = flo(1:ng,1,ilo)
                bas2(1:ng,hybrid%nindx(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),&
                     &                                    itype,ispin) = flo(1:ng,2,ilo)
             END DO

          END IF

       END DO
    END DO

    DEALLOCATE(f,df)

    ! the radial functions are normalized
    DO ispin=1,DIMENSION%jspd
       DO itype=1,atoms%ntype
          DO l=0,atoms%lmax(itype)
             DO i=1,hybrid%nindx(l,itype)
                rdum = intgrf(bas1(:,i,l,itype,ispin)**2&
                     &                     +bas2(:,i,l,itype,ispin)**2,&
                     &                      atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

                bas1(:atoms%jri(itype),i,l,itype,ispin)&
                     &         = bas1(:atoms%jri(itype),i,l,itype,ispin)/SQRT(rdum)
                bas2(:atoms%jri(itype),i,l,itype,ispin)&
                     &         = bas2(:atoms%jri(itype),i,l,itype,ispin)/SQRT(rdum)
             END DO
          END DO
       END DO
    END DO

  
    !
    ! - - - - - - SETUP OF THE MIXED BASIS IN THE IR - - - - - - -
    !


    !
    ! construct G-vectors with cutoff smaller than gcutm
    !
    gcutm = hybrid%gcutm1
    ALLOCATE ( hybrid%ngptm(kpts%nkptf) )

    hybrid%ngptm = 0
    i     = 0
    n     =-1

    rdum1 = MAXVAL( (/ (SQRT(SUM(MATMUL(kpts%bkf(:,ikpt),cell%bmat)**2)),ikpt=1,kpts%nkptf ) /) )

    ! a first run for the determination of the dimensions of the fields
    ! gptm,pgptm

    DO
       n    = n + 1
       ldum = .FALSE.
       DO x = -n,n
          n1 = n-ABS(x)
          DO y = -n1,n1
             n2 = n1-ABS(y)
             DO z = -n2,n2,MAX(2*n2,1)
                g     = (/x,y,z/) 
                rdum  = SQRT(SUM(MATMUL(g,cell%bmat)**2))-rdum1
                IF(rdum.GT. gcutm) CYCLE
                ldum1 = .FALSE.
                DO ikpt = 1,kpts%nkptf
                   kvec = kpts%bkf(:,ikpt)
                   rdum = SUM(MATMUL(kvec+g,cell%bmat)**2)

                   IF(rdum.LE.gcutm**2) THEN
                      IF(.NOT.ldum1) THEN
                         i                     = i + 1
                         ldum1                 = .TRUE.
                      END IF

                      hybrid%ngptm(ikpt) = hybrid%ngptm(ikpt) + 1
                      ldum        = .TRUE.
                   END IF
                END DO
             END DO
          END DO
       END DO
       IF(.NOT.ldum) EXIT
    END DO

    hybrid%gptmd   = i
    hybrid%maxgptm = MAXVAL(hybrid%ngptm)

    ALLOCATE ( hybrid%gptm(3,hybrid%gptmd) )
    ALLOCATE ( hybrid%pgptm(hybrid%maxgptm,kpts%nkptf) )

    hybrid%gptm  = 0
    hybrid%pgptm = 0
    hybrid%ngptm = 0

    i     = 0
    n     =-1

    !     
    ! Allocate and initialize arrays needed for G vector ordering
    !
    ALLOCATE ( unsrt_pgptm(hybrid%maxgptm,kpts%nkptf) )
    ALLOCATE ( length_kG(hybrid%maxgptm,kpts%nkptf) )
    length_kG   = 0
    unsrt_pgptm = 0

    DO
       n    = n + 1
       ldum = .FALSE.
       DO x = -n,n
          n1 = n-ABS(x)
          DO y = -n1,n1
             n2 = n1-ABS(y)
             DO z = -n2,n2,MAX(2*n2,1)
                g     = (/x,y,z/) 
                rdum  = SQRT(SUM(MATMUL(g,cell%bmat)**2))-rdum1
                IF(rdum.GT. gcutm) CYCLE
                ldum1 = .FALSE.
                DO ikpt = 1,kpts%nkptf
                   kvec = kpts%bkf(:,ikpt)
                   rdum = SUM(MATMUL(kvec+g,cell%bmat)**2)

                   IF(rdum.LE.(gcutm)**2) THEN
                      IF(.NOT.ldum1) THEN
                         i          = i + 1
                         hybrid%gptm(:,i)  = g
                         ldum1      = .TRUE.
                      END IF

                      hybrid%ngptm(ikpt)              = hybrid%ngptm(ikpt) + 1
                      ldum                     = .TRUE.

                      ! Position of the vector is saved as pointer
                      unsrt_pgptm(hybrid%ngptm(ikpt),ikpt) = i
                      ! Save length of vector k + G for array sorting
                      length_kG(hybrid%ngptm(ikpt),ikpt)   = rdum
                   END IF
                END DO
             END DO
          END DO
       END DO
       IF(.NOT.ldum) EXIT
    END DO

    !
    ! Sort pointers in array, so that shortest |k+G| comes first
    !
    DO ikpt = 1,kpts%nkptf
       ALLOCATE( ptr(hybrid%ngptm(ikpt)) )
       ! Divide and conquer algorithm for arrays > 1000 entries
       divconq = MAX( 0, INT( 1.443*LOG( 0.001*hybrid%ngptm(ikpt) ) ) )
       ! create pointers which correspond to a sorted array
       CALL rorderpf(ptr, length_kG(1:hybrid%ngptm(ikpt),ikpt),hybrid%ngptm(ikpt), divconq )
       ! rearrange old pointers
       DO igpt = 1,hybrid%ngptm(ikpt)
          hybrid%pgptm(igpt,ikpt) = unsrt_pgptm(ptr(igpt),ikpt)
       END DO
       DEALLOCATE( ptr )
    END DO
    DEALLOCATE( unsrt_pgptm )
    DEALLOCATE( length_kG )

    !
    ! construct IR mixed basis set for the representation of the non local exchange elements
    ! with cutoff gcutm1
    !

    ! first run to determine dimension of pgptm1
    ALLOCATE( hybrid%ngptm1(kpts%nkptf) )
    hybrid%ngptm1 = 0
    DO igpt = 1,hybrid%gptmd
       g = hybrid%gptm(:,igpt)
       DO ikpt = 1,kpts%nkptf
          kvec = kpts%bkf(:,ikpt)
          rdum = SUM(MATMUL(kvec+g,cell%bmat)**2)
          IF( rdum .LE. hybrid%gcutm1**2 ) THEN
             hybrid%ngptm1(ikpt) = hybrid%ngptm1(ikpt) + 1
          END IF
       END DO
    END DO
    hybrid%maxgptm1 = MAXVAL( hybrid%ngptm1 )

    ALLOCATE( hybrid%pgptm1(hybrid%maxgptm1,kpts%nkptf) )
    hybrid%pgptm1 = 0
    hybrid%ngptm1 = 0
    !     
    ! Allocate and initialize arrays needed for G vector ordering
    !
    ALLOCATE ( unsrt_pgptm(hybrid%maxgptm1,kpts%nkptf) )
    ALLOCATE ( length_kG(hybrid%maxgptm1,kpts%nkptf) )
    length_kG   = 0
    unsrt_pgptm = 0
    DO igpt = 1,hybrid%gptmd
       g = hybrid%gptm(:,igpt)
       DO ikpt = 1,kpts%nkptf
          kvec = kpts%bkf(:,ikpt)
          rdum = SUM(MATMUL(kvec+g,cell%bmat)**2)
          IF( rdum .LE. hybrid%gcutm1**2 ) THEN
             hybrid%ngptm1(ikpt)                   = hybrid%ngptm1(ikpt) + 1
             unsrt_pgptm(hybrid%ngptm1(ikpt),ikpt) = igpt
             length_kG(hybrid%ngptm1(ikpt),ikpt)   = rdum
          END IF
       END DO
    END DO

    !
    ! Sort pointers in array, so that shortest |k+G| comes first
    !
    DO ikpt = 1,kpts%nkptf
       ALLOCATE( ptr(hybrid%ngptm1(ikpt)) )
       ! Divide and conquer algorithm for arrays > 1000 entries
       divconq = MAX( 0, INT( 1.443*LOG( 0.001*hybrid%ngptm1(ikpt) ) ) )
       ! create pointers which correspond to a sorted array
       CALL rorderpf(ptr, length_kG(1:hybrid%ngptm1(ikpt),ikpt),hybrid%ngptm1(ikpt), divconq )
       ! rearrange old pointers
       DO igpt = 1,hybrid%ngptm1(ikpt)
          hybrid%pgptm1(igpt,ikpt) = unsrt_pgptm(ptr(igpt),ikpt)
       END DO
       DEALLOCATE( ptr )
    END DO
    DEALLOCATE( unsrt_pgptm )
    DEALLOCATE( length_kG )

    
    IF ( mpi%irank == 0 ) THEN
       WRITE(6,'(/A)')       'Mixed basis'
       WRITE(6,'(A,I5)')     'Number of unique G-vectors: ',hybrid%gptmd
       WRITE(6,*)
       WRITE(6,'(3x,A)') 'IR Plane-wave basis with cutoff of gcutm (hybrid%gcutm1/2*input%rkmax):'
       WRITE(6,'(5x,A,I5)')  'Maximal number of G-vectors:',hybrid%maxgptm
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(3x,A)') 'IR Plane-wave basis for non-local exchange potential:'
       WRITE(6,'(5x,A,I5)')  'Maximal number of G-vectors:',hybrid%maxgptm1
       WRITE(6,*)
    END IF

    !
    ! - - - - - - - - Set up MT product basis for the non-local exchange potential  - - - - - - - - - -
    !

    IF ( mpi%irank == 0 ) THEN
       WRITE(6,'(A)') 'MT product basis for non-local exchange potential:'
       WRITE(6,'(A)') 'Reduction due to overlap (quality of orthonormality, should be < 1.0E-06)'
    END IF

    hybrid%maxlcutm1 = MAXVAL( hybrid%lcutm1 )


    ALLOCATE ( hybrid%nindxm1(0:hybrid%maxlcutm1,atoms%ntype) )
    ALLOCATE ( seleco(hybrid%maxindx,0:atoms%lmaxd) , selecu(hybrid%maxindx,0:atoms%lmaxd) )
    ALLOCATE ( selecmat(hybrid%maxindx,0:atoms%lmaxd,hybrid%maxindx,0:atoms%lmaxd) )
    hybrid%nindxm1 = 0    !!! 01/12/10 jij%M.b.

    !
    ! determine maximal indices of (radial) mixed-basis functions (->nindxm1) 
    ! (will be reduced later-on due to overlap)
    !

    hybrid%maxindxp1  = 0
    DO itype = 1,atoms%ntype
       seleco = .FALSE.
       selecu = .FALSE.
       seleco(1,0:hybrid%select1(1,itype)) = .TRUE.
       selecu(1,0:hybrid%select1(3,itype)) = .TRUE.
       seleco(2,0:hybrid%select1(2,itype)) = .TRUE.
       selecu(2,0:hybrid%select1(4,itype)) = .TRUE.

       ! include local orbitals 
       IF(hybrid%maxindx.GE.3) THEN
          seleco(3:,:) = .TRUE. 
          selecu(3:,:) = .TRUE.
       END IF

       DO l=0,hybrid%lcutm1(itype) 
          n = 0
          M = 0
        
          !
          ! valence * valence
          !

          ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
          selecmat = RESHAPE ( (/ ((((seleco(n1,l1).AND.selecu(n2,l2),&
               n1=1,hybrid%maxindx),l1=0,atoms%lmaxd), n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,&
               (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

          DO l1=0,atoms%lmax(itype)
             DO l2=0,atoms%lmax(itype)
                IF( l.GE.ABS(l1-l2) .AND. l.LE.l1+l2) THEN
                   DO n1=1,hybrid%nindx(l1,itype)
                      DO n2=1,hybrid%nindx(l2,itype)
                         M = M + 1
                         IF(selecmat(n1,l1,n2,l2)) THEN
                            n = n + 1
                            selecmat(n2,l2,n1,l1) = .FALSE. ! prevent double counting of products (a*b = b*a)
                         END IF
                      END DO
                   END DO
                END IF
             END DO
          END DO
          IF(n.EQ.0 .AND. mpi%irank==0) &
               WRITE(6,'(A)') 'mixedbasis: Warning!  No basis-function product of ' // lchar(l) //&
               '-angular momentum defined.'
          hybrid%maxindxp1        = MAX(hybrid%maxindxp1,M)
          hybrid%nindxm1(l,itype) = n*DIMENSION%jspd
       END DO
    END DO
    hybrid%maxindxm1 = MAXVAL(hybrid%nindxm1)

    ALLOCATE ( hybrid%basm1(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
    hybrid%basm1 = 0

    ! Define product bases and reduce them according to overlap

    DO itype=1,atoms%ntype
       seleco = .FALSE.
       selecu = .FALSE.
       seleco(1,0:hybrid%select1(1,itype)) = .TRUE.
       selecu(1,0:hybrid%select1(3,itype)) = .TRUE.
       seleco(2,0:hybrid%select1(2,itype)) = .TRUE.
       selecu(2,0:hybrid%select1(4,itype)) = .TRUE. 
       ! include lo's
       IF(hybrid%maxindx.GE.3) THEN
          seleco(3:,:) = .TRUE. 
          selecu(3:,:) = .TRUE.
       END IF
       IF(atoms%ntype.GT.1 .AND. mpi%irank==0)&
            &    WRITE(6,'(6X,A,I3)') 'Atom type',itype
       ng = atoms%jri(itype)
       DO l=0,hybrid%lcutm1(itype)
          n = hybrid%nindxm1(l,itype)
          ! allow for zero product-basis functions for
          ! current l-quantum number
          IF(n.EQ.0) THEN
             IF ( mpi%irank==0 ) WRITE(6,'(6X,A,'':   0 ->   0'')') lchar(l)
             CYCLE 
          END IF

          ! set up the overlap matrix
          ALLOCATE (olap(n,n),eigv(n,n),work(3*n),eig(n),ihelp(n))
          ihelp    = 1 ! initialize to avoid a segfault
          i        = 0

       
          ! valence*valence

          ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
          selecmat = RESHAPE ( (/ ((((seleco(n1,l1).AND.selecu(n2,l2),&
               n1=1,hybrid%maxindx),l1=0,atoms%lmaxd), n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,  &
               (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

          DO l1=0,atoms%lmax(itype)
             DO l2=0,atoms%lmax(itype)
                IF( l .LT. ABS(l1-l2) .OR. l .GT. l1+l2 ) CYCLE

                DO n1=1,hybrid%nindx(l1,itype)
                   DO n2=1,hybrid%nindx(l2,itype)

                      IF(selecmat(n1,l1,n2,l2)) THEN
                         DO ispin=1,DIMENSION%jspd
                            i = i + 1
                            IF(i.GT.n) call judft_error('got too many product functions',hint='This is a BUG, please report',calledby='mixedbasis')

                            hybrid%basm1(:ng,i,l,itype)  &
                                 = (  bas1(:ng,n1,l1,itype,ispin)&
                                 * bas1(:ng,n2,l2,itype,ispin) &
                                 + bas2(:ng,n1,l1,itype,ispin)&
                                 * bas2(:ng,n2,l2,itype,ispin) )/atoms%rmsh(:ng,itype)

                            !normalize basm1
                            rdum=SQRT(intgrf(hybrid%basm1(:,i,l,itype)**2,&
                                 atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )

                            hybrid%basm1(:ng,i,l,itype) = hybrid%basm1(:ng,i,l,itype)/rdum

                         END DO !ispin
                         ! prevent double counting of products (a*b = b*a)
                         selecmat(n2,l2,n1,l1) = .FALSE. 
                      END IF

                   END DO !n2
                END DO !n1


             END DO !l2
          END DO  !l1

          IF(i .NE. n) call judft_error('counting error for product functions',hint='This is a BUG, please report',calledby='mixedbasis')

          ! In order to get ride of the linear dependencies in the 
          ! radial functions basm1 belonging to fixed l and itype
          ! the overlap matrix is diagonalized and those eigenvectors
          ! with a eigenvalue greater then hybrid%tolerance1 are retained


          ! Calculate overlap
          olap = 0
          DO n2=1,n
             DO n1=1,n2
                olap(n1,n2) = intgrf(hybrid%basm1(:,n1,l,itype) *hybrid%basm1(:,n2,l,itype),&
                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
                olap(n2,n1) = olap(n1,n2)
             END DO
          END DO

          ! Diagonalize
          CALL diagonalize(eigv,eig,olap)

          ! Get rid of linear dependencies (eigenvalue <= hybrid%tolerance1)
          nn = 0
          DO i=1,hybrid%nindxm1(l,itype)
             IF(eig(i).GT.hybrid%tolerance1) THEN 
                nn        = nn + 1
                ihelp(nn) = i
             END IF
          END DO
          hybrid%nindxm1(l,itype) = nn
          eig              = eig(ihelp)
          eigv(:,:)        = eigv(:,ihelp)

          DO i=1,ng
             hybrid%basm1(i,1:nn,l,itype) = MATMUL(hybrid%basm1(i,1:n,l,itype), eigv(:,1:nn))/SQRT(eig(:nn))
          END DO

          ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
          IF(l.EQ.0) THEN

             ! Check if basm1 must be reallocated
             IF(nn+1.GT.SIZE(hybrid%basm1,2)) THEN
                ALLOCATE ( basmhlp(atoms%jmtd,nn+1,0:hybrid%maxlcutm1,atoms%ntype) )
                basmhlp(:,1:nn,:,:) = hybrid%basm1
                DEALLOCATE (hybrid%basm1)
                ALLOCATE ( hybrid%basm1(atoms%jmtd,nn+1,0:hybrid%maxlcutm1,atoms%ntype) )
                hybrid%basm1(:,1:nn,:,:) = basmhlp(:,1:nn,:,:) 
                DEALLOCATE ( basmhlp )
             END IF

             hybrid%basm1(:ng,nn+1,0,itype) = atoms%rmsh(:ng,itype) / SQRT(atoms%rmsh(ng,itype)**3/3)
             DO i = nn,1,-1
                DO j = i+1,nn+1
                   hybrid%basm1(:ng,i,0,itype) = hybrid%basm1(:ng,i,0,itype)&
                        - intgrf(hybrid%basm1(:ng,i,0,itype)*hybrid%basm1(:ng,j,0,itype),&
                        atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
                        * hybrid%basm1(:ng,j,0,itype)

                END DO
                hybrid%basm1(:ng,i,0,itype) = hybrid%basm1(:ng,i,0,itype)&
                     / SQRT(intgrf(hybrid%basm1(:ng,i,0,itype)**2,&
                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf))
             END DO
             nn              = nn + 1
             DEALLOCATE ( olap )
             ALLOCATE   ( olap(nn,nn) )
             hybrid%nindxm1(l,itype) = nn
          END IF

          ! Check orthonormality of product basis
          rdum = 0
          DO i=1,nn
             DO j=1,i
                rdum1 = intgrf(hybrid%basm1(:,i,l,itype)*hybrid%basm1(:,j,l,itype),&
                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

                IF(i.EQ.j) THEN
                   rdum1 = ABS(1-rdum1)
                   rdum  = rdum + rdum1**2
                ELSE
                   rdum1 = ABS(rdum1)
                   rdum  = rdum + 2*rdum1**2
                END IF

                IF(rdum1.GT.1d-6) THEN
                   IF ( mpi%irank == 0 ) THEN
                      WRITE(6,'(A)') 'mixedbasis: Bad orthonormality of ' &
                           //lchar(l)//'-product basis. Increase tolerance.'
                      WRITE(6,'(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of',&
                           rdum1,' occurred for (',i,',',j,')-overlap.'
                   END IF
                   CALL judft_error("Bad orthonormality of product basis",hint='Increase tolerance',calledby='mixedbasis')
                END IF

             END DO
          END DO

          IF ( mpi%irank == 0 ) THEN
             WRITE(6,'(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')') lchar(l)//':',n,nn,SQRT(rdum)/nn
          END IF

          DEALLOCATE(olap,eigv,work,eig,ihelp)

       END DO !l
       IF ( mpi%irank == 0 ) WRITE(6,'(6X,A,I7)') 'Total:', SUM(hybrid%nindxm1(0:hybrid%lcutm1(itype),itype))
    END DO ! itype

    hybrid%maxindxm1 = MAXVAL(hybrid%nindxm1)

    ALLOCATE( basmhlp(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
    basmhlp(1:atoms%jmtd,1:hybrid%maxindxm1,0:hybrid%maxlcutm1,1:atoms%ntype)&
         = hybrid%basm1(1:atoms%jmtd,1:hybrid%maxindxm1,0:hybrid%maxlcutm1,1:atoms%ntype)
    DEALLOCATE(hybrid%basm1)
    ALLOCATE( hybrid%basm1(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
    hybrid%basm1=basmhlp

    DEALLOCATE(basmhlp,seleco,selecu,selecmat)


    !
    ! now we build linear combinations of the radial functions
    ! such that they possess no moment except one radial function in each l-channel
    !
    IF ( mpi%irank == 0 ) THEN
       WRITE(6,'(/,A,/,A)')'Build linear combinations of radial '//&
            'functions in each l-channel,',&
            'such that they possess no multipolmoment'//&
            ' except the last function:'

       WRITE(6,'(/,17x,A)') 'moment  (quality of orthonormality)'
    END IF
    DO itype = 1,atoms%ntype
       ng = atoms%jri(itype)

       IF(atoms%ntype.GT.1 .AND. mpi%irank==0) WRITE(6,'(6X,A,I3)') 'Atom type',itype

       DO l = 0,hybrid%lcutm1(itype)
          ! determine radial function with the largest moment
          ! this function is used to build the linear combinations
          rdum = 0
          DO i = 1,hybrid%nindxm1(l,itype)
             rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
                  atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
             IF( ABS(rdum1) .GT.rdum ) THEN
                n = i
                rdum = rdum1
             END IF
          END DO

          ! rearrange order of radial functions such that the last function possesses the largest moment
          j = 0
          bashlp(:ng) = hybrid%basm1(:ng,n,l,itype)
          DO i = 1,hybrid%nindxm1(l,itype)
             IF( i .EQ. n ) CYCLE
             j = j + 1
             hybrid%basm1(:ng,j,l,itype) = hybrid%basm1(:ng,i,l,itype)
          END DO
          hybrid%basm1(:ng,hybrid%nindxm1(l,itype),l,itype) = bashlp(:ng)

       END DO


       DO l = 0,hybrid%lcutm1(itype)
          IF ( mpi%irank == 0 ) WRITE(6,'(6X,A)') lchar(l)//':'

          IF(hybrid%nindxm1(l,itype).EQ.0) THEN
             IF( mpi%irank == 0 ) WRITE(6,'(6X,A,'':   0 ->    '')') lchar(l)
             CYCLE
          END IF

          n = hybrid%nindxm1(l,itype)
          DO i = 1,hybrid%nindxm1(l,itype)
             IF(i.EQ.n) CYCLE
             ! calculate moment of radial function i
             rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
                  atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

             rdum  = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,n,l,itype),&
                  atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

             bashlp(:ng) = hybrid%basm1(:ng,n,l,itype)

             IF( SQRT(rdum**2 + rdum1**2) .LE. 1E-06 .AND. mpi%irank == 0)&
                  WRITE(6,*) 'Warning: Norm is smaller thann 1E-06!'

             ! change function n such that n is orthogonal to i
             ! since the functions basm1 have been orthogonal on input
             ! the linear combination does not destroy the orthogonality to the residual functions
             hybrid%basm1(:ng,n,l,itype) = rdum /SQRT(rdum**2+rdum1**2) * bashlp(:ng)&
                  + rdum1/SQRT(rdum**2+rdum1**2) * hybrid%basm1(:ng,i,l,itype)

             ! combine basis function i and n so that they possess no momemt
             hybrid%basm1(:ng,i,l,itype) = rdum1/SQRT(rdum**2+rdum1**2) * bashlp(:ng)&
                  - rdum /SQRT(rdum**2+rdum1**2) * hybrid%basm1(:ng,i,l,itype)

             rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
                  atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

             IF( rdum1 .GT. 1E-10 ) call judft_error('moment of radial function does not vanish',calledby='mixedbasis')

             IF ( mpi%irank == 0 ) WRITE(6,'(6x,I4,'' ->  '',ES8.1)') i,rdum1
          END DO

          ! test orthogonality
          rdum = 0
          DO i = 1,hybrid%nindxm1(l,itype)
             DO j = 1,hybrid%nindxm1(l,itype) 
                rdum1 = intgrf(hybrid%basm1(:ng,i,l,itype)*hybrid%basm1(:ng,j,l,itype),&
                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
                IF( i .NE. j ) THEN
                   rdum = rdum + rdum1
                ELSE
                   rdum = rdum + ABS(1 - rdum1)
                END IF
             END DO
          END DO
          IF ( mpi%irank == 0 )&
               WRITE(6,'(6x,I4,'' ->'',f10.5,''   ('',ES8.1,'' )'')')  n,&
               intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,n,l,itype),atoms%jri,&
               atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf),rdum
       END DO
    END DO



    DO itype = 1,atoms%ntype
       IF ( ANY(hybrid%nindxm1(0:hybrid%lcutm1(itype),itype) == 0) ) call judft_error('any hybrid%nindxm1 eq 0',calledby='mixedbasis')
    END DO

    !count basis functions
    hybrid%nbasp   = 0
    DO itype=1,atoms%ntype
       DO i=1,atoms%neq(itype)
          DO l=0,hybrid%lcutm1(itype)
             DO M=-l,l
                DO j=1,hybrid%nindxm1(l,itype)
                   hybrid%nbasp = hybrid%nbasp + 1
                END DO
             END DO
          END DO
       END DO
    END DO
    hybrid%maxbasm1  = hybrid%nbasp  + hybrid%maxgptm
    DO nk = 1,kpts%nkptf
       hybrid%nbasm(nk) = hybrid%nbasp + hybrid%ngptm(nk)
    END DO

     hybrid%maxlmindx = MAXVAL((/ ( SUM( (/ (hybrid%nindx(l,itype)*(2*l+1), l=0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )

    !
    !
    ! Write all information to a file to enable restarting
    !
    IF ( l_restart .AND. mpi%irank == 0 ) THEN

       OPEN(UNIT=iounit,FILE=ioname,FORM='unformatted', STATUS='replace')
 
       WRITE(iounit) kpts%nkptf,hybrid%gptmd
       WRITE(iounit) hybrid%maxgptm,hybrid%maxindx
       WRITE(iounit) hybrid%ngptm,hybrid%gptm,hybrid%pgptm,hybrid%nindx

       WRITE(iounit) hybrid%maxgptm1,hybrid%maxlcutm1,hybrid%maxindxm1,hybrid%maxindxp1
       WRITE(iounit) hybrid%ngptm1,hybrid%pgptm1,hybrid%nindxm1
       WRITE(iounit) hybrid%basm1

   
       CLOSE(iounit)

       l_restart = .FALSE.

    END IF

  END SUBROUTINE mixedbasis


  SUBROUTINE read_atompotential(atoms,fname, potential)

    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms

    ! - scalars -
    CHARACTER(10),INTENT(IN) ::   fname(atoms%ntype)

    ! - arrays -
    REAL        ,INTENT(OUT)::   potential(atoms%jmtd,atoms%ntype)
    ! - local scalars -
    INTEGER                 ::   i,j,itype,ic,m     
    REAL                    ::   rdum,r ,n
    ! - local arrays -
    REAL    ,ALLOCATABLE    ::   v_atom (:),r_atom(:)

    potential =0
    DO itype = 1,atoms%ntype
       OPEN(331,file=fname(itype),form='formatted')

       ic = 0
       DO 
          READ(331,*) rdum
          ic = ic + 1
          IF( rdum .GT. (atoms%rmt(itype)+0.2) ) EXIT
       END DO

       ALLOCATE( v_atom(ic),r_atom(ic) )
       REWIND(331)

       DO i = 1,ic
          READ(331,*) r_atom(i),v_atom(i)
       END DO

       ! transfer potential from grid used in the atom program
       ! to those used here
       DO i = 1,atoms%jri(itype)
          r = atoms%rmsh(i,itype)
          DO j = 1,ic-1
             IF( r_atom(j) .GT. r  ) THEN
                M = (v_atom(j)-v_atom(j+1))/(r_atom(j)-r_atom(j+1))
                n = v_atom(j) - M*r_atom(j)
                potential(i,itype) =n + M*r
                !               WRITE(333,'(4f25.10)') n + m*r_atom(j),v_atom(j),n + m*r_atom(j+1),v_atom(j+1)
                EXIT
             END IF
          END DO
       END DO

       DEALLOCATE( v_atom,r_atom )
       CLOSE(331) 
    END DO

  END SUBROUTINE read_atompotential

END MODULE m_mixedbasis

