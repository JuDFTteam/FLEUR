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
 
      CONTAINS

      SUBROUTINE mixedbasis( atoms,kpts,&
     &                       sphhar,&
     &                       dimension,input,&
     &                       enpara,&
     &                       cell,vacuum,sym,&
     &                       oneD,stars,&
     &                       xcpot,&
     &                       hybrid,&
     &                       eig_id,&
     &                       mpi,l_restart)



      USE m_radfun,   ONLY : radfun
      USE m_radflo,   ONLY : radflo
      USE m_loddop,   ONLY : loddop
      USE m_util,     ONLY : intgrf_init,intgrf,rorderpf
      USE m_gen_bz
      USE m_read_core
      USE m_wrapper
      USE m_icorrkeys
      USE m_eig66_io
      USE m_types
      IMPLICIT NONE

      TYPE(t_xcpot),INTENT(IN)   :: xcpot

      TYPE(t_mpi),INTENT(IN)   :: mpi

      TYPE(t_dimension),INTENT(IN)   :: dimension

      TYPE(t_oneD),INTENT(IN)   :: oneD

      TYPE(t_hybrid),INTENT(INOUT)   :: hybrid

      TYPE(t_enpara),INTENT(IN)   :: enpara

      TYPE(t_input),INTENT(IN)   :: input

      TYPE(t_vacuum),INTENT(IN)   :: vacuum

      TYPE(t_sym),INTENT(IN)   :: sym

      TYPE(t_stars),INTENT(IN)   :: stars

      TYPE(t_cell),INTENT(IN)   :: cell

      TYPE(t_kpts),INTENT(INOUT)   :: kpts

      TYPE(t_sphhar),INTENT(IN)   :: sphhar

      TYPE(t_atoms),INTENT(IN)   :: atoms


      ! - scalars -
      INTEGER,INTENT(IN)              :: eig_id
      
      LOGICAL,INTENT(INOUT)           ::  l_restart
      
      

      ! - local scalars -
      INTEGER                         ::  idum,iter,ilo
      INTEGER                         ::  ikpt
      INTEGER                         ::  ispin,itype,l1,l2,l,n,igpt,&
     &                                    n1,n2,nn,i,j,ic ,ng      
      INTEGER                         ::  jsp
      INTEGER                         ::  nodem,noded
      INTEGER                         ::  nrec,m
      INTEGER                         ::  x,y,z
      INTEGER                         ::  maxindxc,lmaxcd
      INTEGER                         ::  divconq ! use Divide & Conquer algorithm for array sorting (>0: yes, =0: no)

      REAL                            ::  gcutm
      REAL                            :: wronk
      REAL                            ::  rdum,rdum1,rdum2
    
      COMPLEX                         ::  cdum

      LOGICAL                         ::  ldum,ldum1
      LOGICAL                         ::  exx
  

      ! - local arrays -
      INTEGER                         ::  g(3)
      INTEGER                         ::  lmaxc(atoms%ntype)
      INTEGER,ALLOCATABLE             ::  nindxc(:,:)
      INTEGER,ALLOCATABLE             ::  ihelp(:) 
      INTEGER,ALLOCATABLE             ::  ptr(:)           ! pointer for array sorting
      INTEGER,ALLOCATABLE             ::  unsrt_pgptm(:,:) ! unsorted pointers to g vectors

      REAL                            ::  kvec(3),bkpt(3)
      REAL                            ::  flo(atoms%jmtd,2,atoms%nlod)
      TYPE(t_usdus)                   ::  usdus
      REAL                            ::  uuilon(atoms%nlod,atoms%ntype),&
     &                                    duilon(atoms%nlod,atoms%ntype)
      REAL                            ::  ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
      REAL                            ::  potatom(atoms%jmtd,atoms%ntype)
      REAL                            ::  el(0:atoms%lmaxd,atoms%ntype,dimension%jspd)
      REAL                            ::  ello(atoms%nlod,atoms%ntype,dimension%jspd)
      REAL                            ::  bashlp(atoms%jmtd)

      REAL   ,ALLOCATABLE             ::  f(:,:,:),df(:,:,:),vr0(:,:,:)
      REAL   ,ALLOCATABLE             ::  vr(:,:,:,:),vz(:,:,:)
      REAL   ,ALLOCATABLE             ::  olap(:,:),work(:),eig(:),&
     &                                    eigv(:,:)      
      REAL   ,ALLOCATABLE             ::  bas1(:,:,:,:,:),&
     &                                    bas2(:,:,:,:,:)
      REAL   ,ALLOCATABLE             ::  basmhlp(:,:,:,:)
      REAL   ,ALLOCATABLE             ::  gridf(:,:)
      REAL   ,ALLOCATABLE             ::  core1(:,:,:,:,:),&
     &                                    core2(:,:,:,:,:)
      REAL   ,ALLOCATABLE             ::  eig_c(:,:,:,:)
      REAL   ,ALLOCATABLE             ::  length_kg(:,:) ! length of the vectors k + G

      COMPLEX,ALLOCATABLE             ::  vpw(:,:),vzxy(:,:,:,:) 

      LOGICAL,ALLOCATABLE             ::  selecmat(:,:,:,:)
      LOGICAL,ALLOCATABLE             ::  seleco(:,:),selecu(:,:)

      CHARACTER, PARAMETER            :: lchar(0:38) =&
     &      (/'s','p','d','f','g','h','i','j','k','l','m','n','o',&
     &        'x','x','x','x','x','x','x','x','x','x','x','x','x',&
     &        'x','x','x','x','x','x','x','x','x','x','x','x','x' /)

      CHARACTER(len=2)                ::  nchar
      CHARACTER(len=2)                ::  noel(atoms%ntype)
      CHARACTER(len=10)               ::  fname(atoms%ntype)
      CHARACTER(8)                    ::  dop,iop,name(10)
      ! writing to a file
      INTEGER, PARAMETER              ::  iounit = 125
      CHARACTER(10), PARAMETER        ::  ioname = 'mixbas'
      LOGICAL                         ::  l_found

      IF ( mpi%irank == 0 )&
     &   WRITE(6,'(//A,I2,A)') '### subroutine: mixedbasis ###'

    
      exx = ( xcpot%icorr .EQ. icorr_exx )

      ! Deallocate arrays which might have been allocated in a previous run of this subroutine
      IF(allocated(hybrid%ngptm))    DEALLOCATE(hybrid%ngptm)
      IF(allocated(hybrid%ngptm1))   DEALLOCATE(hybrid%ngptm1)
      IF(allocated(hybrid%ngptm2))   DEALLOCATE(hybrid%ngptm2)
      IF(allocated(hybrid%nindxm1))  DEALLOCATE(hybrid%nindxm1)
      IF(allocated(hybrid%nindxm2))  DEALLOCATE(hybrid%nindxm2)
      IF(allocated(hybrid%pgptm))    DEALLOCATE(hybrid%pgptm)
      IF(allocated(hybrid%pgptm1))   DEALLOCATE(hybrid%pgptm1)
      IF(allocated(hybrid%pgptm2))   DEALLOCATE(hybrid%pgptm2)
      IF(allocated(vr0))      DEALLOCATE(vr0)
      IF(allocated(hybrid%gptm))     DEALLOCATE(hybrid%gptm)
      IF(allocated(hybrid%basm1) )   DEALLOCATE(hybrid%basm1)
      IF(allocated(hybrid%basm2) )   DEALLOCATE(hybrid%basm2)

      ! If restart is specified read file if it already exists
      ! create it otherwise
      IF ( l_restart ) THEN

        ! Test if file exists
        INQUIRE(FILE=ioname,EXIST=l_found)

        IF ( l_found ) THEN

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

          IF ( exx ) THEN
            ! Read array sizes
            READ(iounit) hybrid%maxgptm2,hybrid%maxlcutm2,hybrid%maxindxm2,hybrid%maxindxp2
            ! Allocate necessary array size and read arrays
            ALLOCATE ( hybrid%ngptm2(kpts%nkptf),hybrid%pgptm2(hybrid%maxgptm2,kpts%nkptf),&
     &                 hybrid%nindxm2(0:hybrid%maxlcutm2,atoms%ntype) )
            ALLOCATE ( hybrid%basm2(atoms%jmtd,hybrid%maxindxm2,0:hybrid%maxlcutm2,atoms%ntype) )
            READ(iounit) hybrid%ngptm2,hybrid%pgptm2,hybrid%nindxm2
            READ(iounit) hybrid%basm2
          END IF

          CLOSE(iounit)

          RETURN
        END IF
      END IF

      hybrid%maxindx = 0
      DO itype = 1,atoms%ntype
        DO l = 0,atoms%lmax(itype)
          hybrid%maxindx = max(hybrid%maxindx,2+count( atoms%llo(:atoms%nlo(itype),itype) .eq. l))
         END DO
      END DO
!       maxindx   = maxval( nlo   ) + 2

      

      ! initialize gridf for radial integration
      CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,gridf)

      ! generate whole BZ form k-points given in kpts
      CALL gen_bz(kpts,sym)

      !
      ! read in energy parameters from file eig
      ! to avoid meaningless energy parameters which occur in the case
      ! that the energy parameter is set to the atomic prinicipal 
      ! quantum number
      ! (el0 and ello0 are just the values in the enpara file)
      ! 


      DO jsp=1,dimension%jspd
         STOP "TODO"
 !       CALL read_eig(eig_id,nkpti,jsp,el=el(:,:,jsp),&
  !   &   ello=ello(:,:,jsp),kpts=bkpt)
      ENDDO

      ! load potential from file pottot
      ALLOCATE ( vpw(stars%ng3,dimension%jspd),vzxy(vacuum%nmzxyd,oneD%odi%n2d-1,2,dimension%jspd) )
      ALLOCATE ( vz(vacuum%nmzd,2,4), vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,dimension%jspd) )
      ALLOCATE ( vr0(atoms%jmtd,atoms%ntype,dimension%jspd) )

      OPEN (8,file='pottot',form='unformatted',status='old')
      STOP "TODO"
      !      CALL loddop(&
!     &            dimension,stars,oneD%n2d,vacuum,atoms,sphhar,&
!     &            input%nq2,sym,sym,&
!     &            8,&
!     &            iop,dop,iter,vr,vpw,vz,vzxy,name)
      CLOSE(8)
      vr0(:,:,:) = vr(:,0,:,:) 
      DEALLOCATE (vpw,vzxy,vz,vr)

      ! calculate radial basisfunctions u and u' with
      ! the spherical part of the potential vr0 and store them in
      ! bas1 = large component ,bas2 = small component

      ALLOCATE( f(atoms%jmtd,2,0:atoms%lmaxd), df(atoms%jmtd,2,0:atoms%lmaxd) )
      ALLOCATE( bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype,dimension%jspd) )
      ALLOCATE( bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype,dimension%jspd) )
      
      DO itype=1,atoms%ntype
        ng = atoms%jri(itype)
        DO ispin=1,dimension%jspd
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
      DO ispin=1,dimension%jspd
        DO itype=1,atoms%ntype
          DO l=0,atoms%lmax(itype)
            DO i=1,hybrid%nindx(l,itype)
              rdum = intgrf(bas1(:,i,l,itype,ispin)**2&
     &                     +bas2(:,i,l,itype,ispin)**2,&
     &                      atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

              bas1(:atoms%jri(itype),i,l,itype,ispin)&
     &         = bas1(:atoms%jri(itype),i,l,itype,ispin)/sqrt(rdum)
              bas2(:atoms%jri(itype),i,l,itype,ispin)&
     &         = bas2(:atoms%jri(itype),i,l,itype,ispin)/sqrt(rdum)
            END DO
          END DO
        END DO
      END DO

      IF( exx ) THEN
        ! generate core wave functions (-> core1/2(jmtd,nindxc,0:lmaxc,ntype,jspd) )   
      
        CALL core_init( dimension,input,atoms,&
     &                  lmaxcd,maxindxc)
      
   
        ALLOCATE( nindxc(0:lmaxcd,atoms%ntype) )
        ALLOCATE( core1(atoms%jmtd,maxindxc,0:lmaxcd,atoms%ntype,dimension%jspd) )      
        ALLOCATE( core2(atoms%jmtd,maxindxc,0:lmaxcd,atoms%ntype,dimension%jspd) )
        ALLOCATE( eig_c(maxindxc,0:lmaxcd,atoms%ntype,dimension%jspd)      )
      
        DO ispin = 1,input%jspins
          CALL corewf( atoms,ispin,input,dimension,vr0,&
     &                 lmaxcd,maxindxc,mpi,&
     &                 lmaxc,nindxc,core1(:,:,:,:,ispin),&
     &                 core2(:,:,:,:,ispin),eig_c(:,:,:,ispin))

        END DO
      END IF ! exx

      !
      ! - - - - - - SETUP OF THE MIXED BASIS IN THE IR - - - - - - -
      !


      !
      ! construct G-vectors with cutoff smaller than gcutm
      !
      IF( exx ) THEN
        gcutm = 2*input%rkmax
      ELSE
        gcutm = hybrid%gcutm1
      END IF
      ALLOCATE ( hybrid%ngptm(kpts%nkptf) )
     
      hybrid%ngptm = 0
      i     = 0
      n     =-1
      
      rdum1 = maxval( (/ (sqrt(sum(matmul(kpts%bk(:,ikpt),cell%bmat)**2)),&
     &                   ikpt=1,kpts%nkptf ) /) )
      
      ! a first run for the determination of the dimensions of the fields
      ! gptm,pgptm
      
      DO
        n    = n + 1
        ldum = .false.
        DO x = -n,n
          n1 = n-abs(x)
          DO y = -n1,n1
            n2 = n1-abs(y)
            DO z = -n2,n2,max(2*n2,1)
              g     = (/x,y,z/) 
              rdum  = sqrt(sum(matmul(g,cell%bmat)**2))-rdum1
              IF(rdum.gt. gcutm) CYCLE
              ldum1 = .false.
              DO ikpt = 1,kpts%nkptf
                kvec = kpts%bk(:,ikpt)
                rdum = sum(matmul(kvec+g,cell%bmat)**2)

                IF(rdum.le.gcutm**2) THEN
                  IF(.not.ldum1) THEN
                    i                     = i + 1
                    ldum1                 = .true.
                  END IF
                   
                  hybrid%ngptm(ikpt) = hybrid%ngptm(ikpt) + 1
                  ldum        = .true.
                END IF
              END DO
            END DO
          END DO
        END DO
        IF(.not.ldum) EXIT
      END DO
 
      hybrid%gptmd   = i
      hybrid%maxgptm = maxval(hybrid%ngptm)

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
        ldum = .false.
        DO x = -n,n
          n1 = n-abs(x)
          DO y = -n1,n1
            n2 = n1-abs(y)
            DO z = -n2,n2,max(2*n2,1)
              g     = (/x,y,z/) 
              rdum  = sqrt(sum(matmul(g,cell%bmat)**2))-rdum1
              IF(rdum.gt. gcutm) CYCLE
              ldum1 = .false.
              DO ikpt = 1,kpts%nkptf
                kvec = kpts%bk(:,ikpt)
                rdum = sum(matmul(kvec+g,cell%bmat)**2)

                IF(rdum.le.(gcutm)**2) THEN
                  IF(.not.ldum1) THEN
                    i          = i + 1
                    hybrid%gptm(:,i)  = g
                    ldum1      = .true.
                  END IF
                   
                  hybrid%ngptm(ikpt)              = hybrid%ngptm(ikpt) + 1
                  ldum                     = .true.

                  ! Position of the vector is saved as pointer
                  unsrt_pgptm(hybrid%ngptm(ikpt),ikpt) = i
                  ! Save length of vector k + G for array sorting
                  length_kG(hybrid%ngptm(ikpt),ikpt)   = rdum
                END IF
              END DO
            END DO
          END DO
        END DO
        IF(.not.ldum) EXIT
      END DO

      !
      ! Sort pointers in array, so that shortest |k+G| comes first
      !
      do ikpt = 1,kpts%nkptf
         ALLOCATE( ptr(hybrid%ngptm(ikpt)) )
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = max( 0, int( 1.443*log( 0.001*hybrid%ngptm(ikpt) ) ) )
         ! create pointers which correspond to a sorted array
         call rorderpf(&
     &        ptr,&
     &        length_kG(1:hybrid%ngptm(ikpt),ikpt),hybrid%ngptm(ikpt), divconq )
         ! rearrange old pointers
         do igpt = 1,hybrid%ngptm(ikpt)
            hybrid%pgptm(igpt,ikpt) = unsrt_pgptm(ptr(igpt),ikpt)
         end do
         DEALLOCATE( ptr )
      end do
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
          kvec = kpts%bk(:,ikpt)
          rdum = sum(matmul(kvec+g,cell%bmat)**2)
          IF( rdum .le. hybrid%gcutm1**2 ) THEN
            hybrid%ngptm1(ikpt) = hybrid%ngptm1(ikpt) + 1
          END IF          
        END DO
      END DO
      hybrid%maxgptm1 = maxval( hybrid%ngptm1 )
      
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
          kvec = kpts%bk(:,ikpt)
          rdum = sum(matmul(kvec+g,cell%bmat)**2)
          IF( rdum .le. hybrid%gcutm1**2 ) THEN
            hybrid%ngptm1(ikpt)                   = hybrid%ngptm1(ikpt) + 1
            unsrt_pgptm(hybrid%ngptm1(ikpt),ikpt) = igpt
            length_kG(hybrid%ngptm1(ikpt),ikpt)   = rdum
          END IF          
        END DO
      END DO
      
      !
      ! Sort pointers in array, so that shortest |k+G| comes first
      !
      do ikpt = 1,kpts%nkptf
         ALLOCATE( ptr(hybrid%ngptm1(ikpt)) )
         ! Divide and conquer algorithm for arrays > 1000 entries
         divconq = max( 0, int( 1.443*log( 0.001*hybrid%ngptm1(ikpt) ) ) )
         ! create pointers which correspond to a sorted array
         call rorderpf(&
     &        ptr,&
     &        length_kG(1:hybrid%ngptm1(ikpt),ikpt),hybrid%ngptm1(ikpt), divconq )
         ! rearrange old pointers
         do igpt = 1,hybrid%ngptm1(ikpt)
            hybrid%pgptm1(igpt,ikpt) = unsrt_pgptm(ptr(igpt),ikpt)
         end do
         DEALLOCATE( ptr )
      end do
      DEALLOCATE( unsrt_pgptm )
      DEALLOCATE( length_kG )

      IF( exx ) THEN
        !
        ! construct IR mixed basis set for the representation of the response function
        ! with cutoff gcutm2 
        !


        ALLOCATE( hybrid%ngptm2(kpts%nkptf) )
        hybrid%ngptm2 = 0

        DO igpt = 1,hybrid%gptmd
          g = hybrid%gptm(:,igpt)
          DO ikpt = 1,kpts%nkptf
            kvec = kpts%bk(:,ikpt)
            rdum = sum(matmul(kvec+g,cell%bmat)**2)
            IF( rdum .le. hybrid%gcutm2**2 ) THEN
              hybrid%ngptm2(ikpt) = hybrid%ngptm2(ikpt) + 1
            END IF          
          END DO
        END DO
        hybrid%maxgptm2 = maxval( hybrid%ngptm2 )

        ALLOCATE( hybrid%pgptm2(hybrid%maxgptm2,kpts%nkptf) )
        hybrid%pgptm2 = 0
        hybrid%ngptm2 = 0
        !     
        ! Allocate and initialize arrays needed for G vector ordering
        !
        ALLOCATE ( unsrt_pgptm(hybrid%maxgptm2,kpts%nkptf) )
        ALLOCATE ( length_kG(hybrid%maxgptm2,kpts%nkptf) )
        length_kG   = 0
        unsrt_pgptm = 0
        DO igpt = 1,hybrid%gptmd
          g = hybrid%gptm(:,igpt)
          DO ikpt = 1,kpts%nkptf
            kvec = kpts%bk(:,ikpt)
            rdum = sum(matmul(kvec+g,cell%bmat)**2)
            IF( rdum .le. hybrid%gcutm2**2 ) THEN
              hybrid%ngptm2(ikpt)                   = hybrid%ngptm2(ikpt) + 1
              unsrt_pgptm(hybrid%ngptm2(ikpt),ikpt) = igpt
              length_kG(hybrid%ngptm2(ikpt),ikpt)   = rdum
            END IF          
          END DO
        END DO

        !
        ! Sort pointers in array, so that shortest |k+G| comes first
        !
        do ikpt = 1,kpts%nkptf
           ALLOCATE( ptr(hybrid%ngptm2(ikpt)) )
           ! Divide and conquer algorithm for arrays > 1000 entries
           divconq = max( 0, int( 1.443*log( 0.001*hybrid%ngptm2(ikpt) ) ) )
           ! create pointers which correspond to a sorted array
        call rorderpf(&
     &        ptr,&
     &        length_kG(1:hybrid%ngptm2(ikpt),ikpt),hybrid%ngptm2(ikpt), divconq )
           ! rearrange old pointers
           do igpt = 1,hybrid%ngptm2(ikpt)
              hybrid%pgptm2(igpt,ikpt) = unsrt_pgptm(ptr(igpt),ikpt)
           end do
           DEALLOCATE( ptr )
        end do
        DEALLOCATE( unsrt_pgptm )
        DEALLOCATE( length_kG )

      END IF

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(/A)')       'Mixed basis'
        WRITE(6,'(A,I5)')     'Number of unique G-vectors: ',hybrid%gptmd
        WRITE(6,*)
        WRITE(6,'(3x,A)')&
     &     'IR Plane-wave basis with cutoff of gcutm (hybrid%gcutm1/2*input%rkmax):'
        WRITE(6,'(5x,A,I5)')  'Maximal number of G-vectors:',hybrid%maxgptm
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,'(3x,A)')&
     &     'IR Plane-wave basis for non-local exchange potential:'
        WRITE(6,'(5x,A,I5)')  'Maximal number of G-vectors:',hybrid%maxgptm1
        WRITE(6,*)
        IF( exx ) THEN
          WRITE(6,'(3x,A)')&
     &     'IR Plane-wave basis for response function:'
          WRITE(6,'(5x,A,I5)')  'Maximal number of G-vectors:',hybrid%maxgptm2
          WRITE(6,*)
        END IF
      END IF

      !
      ! - - - - - - - - Set up MT product basis for the non-local exchange potential  - - - - - - - - - -
      !

      IF ( mpi%irank == 0 ) THEN
        WRITE(6,'(A)')&
     &       'MT product basis for non-local exchange potential:'
        WRITE(6,'(A)') 'Reduction due to overlap (quality of &&
     &orthonormality, should be < 1.0E-06)'
      END IF

      hybrid%maxlcutm1 = maxval( hybrid%lcutm1 )


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
        seleco = .false.
        selecu = .false.
        seleco(1,0:hybrid%select1(1,itype)) = .true.
        selecu(1,0:hybrid%select1(3,itype)) = .true.
        seleco(2,0:hybrid%select1(2,itype)) = .true.
        selecu(2,0:hybrid%select1(4,itype)) = .true.

        ! include local orbitals 
        IF(hybrid%maxindx.ge.3) THEN
          seleco(3:,:) = .true. 
          selecu(3:,:) = .true.
        END IF

        DO l=0,hybrid%lcutm1(itype) 
          n = 0
          M = 0
          IF( exx ) THEN
            !
            ! core * valence products
            !
            DO l1 = 0,lmaxc(itype)
              DO l2 = 0,atoms%lmax(itype)
                IF( l .lt. abs(l1-l2) .or. l .gt. l1+l2 ) CYCLE
                DO n1 = 1,nindxc(l1,itype)             
                  DO n2 = 1,hybrid%nindx(l2,itype)
                    M = M + 1
                    IF( .not. seleco(n2,l2) ) CYCLE
                    n = n + 1
                  END DO
                END DO
              END DO
            END DO
            
          END IF

          !
          ! valence * valence
          !

          ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
          selecmat = reshape ( (/ ((((&
     &               seleco(n1,l1).and.selecu(n2,l2),&
     &               n1=1,hybrid%maxindx),l1=0,atoms%lmaxd),           &
     &               n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,&
     &               (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

          DO l1=0,atoms%lmax(itype)
            DO l2=0,atoms%lmax(itype)
              IF( l.ge.abs(l1-l2) .and. l.le.l1+l2) THEN
                DO n1=1,hybrid%nindx(l1,itype)
                  DO n2=1,hybrid%nindx(l2,itype)
                    M = M + 1
                    IF(selecmat(n1,l1,n2,l2)) THEN
                      n = n + 1
                      selecmat(n2,l2,n1,l1) = .false. ! prevent double counting of products (a*b = b*a)
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO 
          IF(n.eq.0 .AND. mpi%irank==0) &
     &      WRITE(6,'(A)') 'mixedbasis: Warning!  No basis-function '//&
     &                     'product of ' // lchar(l) //&
     &                     '-angular momentum defined.'
          hybrid%maxindxp1        = max(hybrid%maxindxp1,M)
          hybrid%nindxm1(l,itype) = n*dimension%jspd
        END DO
      END DO
      hybrid%maxindxm1 = maxval(hybrid%nindxm1)

      ALLOCATE ( hybrid%basm1(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
      hybrid%basm1 = 0

      ! Define product bases and reduce them according to overlap

      DO itype=1,atoms%ntype
        seleco = .false.
        selecu = .false.
        seleco(1,0:hybrid%select1(1,itype)) = .true.
        selecu(1,0:hybrid%select1(3,itype)) = .true.
        seleco(2,0:hybrid%select1(2,itype)) = .true.
        selecu(2,0:hybrid%select1(4,itype)) = .true. 
        ! include lo's
        IF(hybrid%maxindx.ge.3) THEN
          seleco(3:,:) = .true. 
          selecu(3:,:) = .true.
        END IF
        IF(atoms%ntype.gt.1 .AND. mpi%irank==0)&
     &    WRITE(6,'(6X,A,I3)') 'Atom type',itype
        ng = atoms%jri(itype)
        DO l=0,hybrid%lcutm1(itype)
          n = hybrid%nindxm1(l,itype)
          ! allow for zero product-basis functions for
          ! current l-quantum number
          IF(n.eq.0) THEN
            IF ( mpi%irank==0 ) WRITE(6,'(6X,A,'':   0 ->   0'')') lchar(l)
            CYCLE 
          END IF

          ! set up the overlap matrix
          ALLOCATE (olap(n,n),eigv(n,n),work(3*n),eig(n),ihelp(n))
          ihelp    = 1 ! initialize to avoid a segfault
          i        = 0

          IF( exx ) THEN
          ! core*valence
            DO l1 = 0,lmaxc(itype)
              DO l2 = 0,atoms%lmax(itype)
                IF( l .lt. abs(l1-l2) .or. l .gt. l1+l2 ) CYCLE
                DO n1 = 1,nindxc(l1,itype)
                  DO n2 = 1,hybrid%nindx(l2,itype)
                    IF( .not. seleco(n2,l2) ) CYCLE
                    DO ispin = 1,dimension%jspd  
                      i = i + 1
                      IF(i.gt.n)&
     &STOP 'mixedbasis:got too many product functions (bug).'
                      hybrid%basm1(:ng,i,l,itype)  &
     &              = (  core1(:ng,n1,l1,itype,ispin)&
     &                 * bas1(:ng,n2,l2,itype,ispin) &
     &                 + core2(:ng,n1,l1,itype,ispin)&
     &                 * bas2(:ng,n2,l2,itype,ispin)  )/atoms%rmsh(:ng,itype)

                      !normalize basm1
                      rdum=sqrt(intgrf( hybrid%basm1(:,i,l,itype)**2,&
     &                          atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )

                      hybrid%basm1(:ng,i,l,itype)=hybrid%basm1(:ng,i,l,itype)/rdum

                    END DO  !ispin
                  END DO  !n2
                END DO  !n1
              END DO  !l2
            END DO  !l1
          END IF !exx

          ! valence*valence

          ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
          selecmat = reshape ( (/ ((((&
     &               seleco(n1,l1).and.selecu(n2,l2),&
     &               n1=1,hybrid%maxindx),l1=0,atoms%lmaxd),      &
     &               n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,  &
     &               (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

          DO l1=0,atoms%lmax(itype)
            DO l2=0,atoms%lmax(itype)
              IF( l .lt. abs(l1-l2) .or. l .gt. l1+l2 ) CYCLE

              DO n1=1,hybrid%nindx(l1,itype)
                DO n2=1,hybrid%nindx(l2,itype)

                  IF(selecmat(n1,l1,n2,l2)) THEN
                    DO ispin=1,dimension%jspd
                      i = i + 1
                      IF(i.gt.n)&
     & STOP 'mixedbasis:got too many product functions (bug).'

                      hybrid%basm1(:ng,i,l,itype)  &
     &              = (  bas1(:ng,n1,l1,itype,ispin)&
     &                  * bas1(:ng,n2,l2,itype,ispin) &
     &                 + bas2(:ng,n1,l1,itype,ispin)&
     &                  * bas2(:ng,n2,l2,itype,ispin) )/atoms%rmsh(:ng,itype)

                      !normalize basm1
                      rdum=sqrt(intgrf(hybrid%basm1(:,i,l,itype)**2,&
     &                             atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )

                      hybrid%basm1(:ng,i,l,itype) = hybrid%basm1(:ng,i,l,itype)/rdum

                    END DO !ispin
                    ! prevent double counting of products (a*b = b*a)
                    selecmat(n2,l2,n1,l1) = .false. 
                  END IF

                END DO !n2
              END DO !n1

        
            END DO !l2
          END DO  !l1

          IF(i .ne. n)&
     &    STOP 'mixedbasis: counting error for product functions (bug).'

          ! In order to get ride of the linear dependencies in the 
          ! radial functions basm1 belonging to fixed l and itype
          ! the overlap matrix is diagonalized and those eigenvectors
          ! with a eigenvalue greater then hybrid%tolerance1 are retained


          ! Calculate overlap
          olap = 0
          DO n2=1,n
            DO n1=1,n2
              olap(n1,n2) = intgrf(hybrid%basm1(:,n1,l,itype)&
     &                               *hybrid%basm1(:,n2,l,itype),&
     &                             atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
              olap(n2,n1) = olap(n1,n2)
            END DO
          END DO

          ! Diagonalize
          CALL diagonalize(eigv,eig,olap)

          ! Get rid of linear dependencies (eigenvalue <= hybrid%tolerance1)
          nn = 0
          DO i=1,hybrid%nindxm1(l,itype)
            IF(eig(i).gt.hybrid%tolerance1) THEN 
              nn        = nn + 1
              ihelp(nn) = i
            END IF
          END DO
          hybrid%nindxm1(l,itype) = nn
          eig              = eig(ihelp)
          eigv(:,:)        = eigv(:,ihelp)

          DO i=1,ng
            hybrid%basm1(i,1:nn,l,itype) = matmul(hybrid%basm1(i,1:n,l,itype),&
     &                                     eigv(:,1:nn))/sqrt(eig(:nn))
          END DO

          ! Add constant function to l=0 basis and then do a Gram-Schmidt orthonormalization
          IF(l.eq.0) THEN
            
            ! Check if basm1 must be reallocated
            IF(nn+1.gt.size(hybrid%basm1,2)) THEN
              ALLOCATE ( basmhlp(atoms%jmtd,nn+1,0:hybrid%maxlcutm1,atoms%ntype) )
              basmhlp(:,1:nn,:,:) = hybrid%basm1
              DEALLOCATE (hybrid%basm1)
              ALLOCATE ( hybrid%basm1(atoms%jmtd,nn+1,0:hybrid%maxlcutm1,atoms%ntype) )
              hybrid%basm1(:,1:nn,:,:) = basmhlp(:,1:nn,:,:) 
              DEALLOCATE ( basmhlp )
            END IF
            
            hybrid%basm1(:ng,nn+1,0,itype) = atoms%rmsh(:ng,itype)&
     &                              / sqrt(atoms%rmsh(ng,itype)**3/3)
            DO i = nn,1,-1
              DO j = i+1,nn+1
                hybrid%basm1(:ng,i,0,itype) = hybrid%basm1(:ng,i,0,itype)&
     &               - intgrf(hybrid%basm1(:ng,i,0,itype)*hybrid%basm1(:ng,j,0,itype),&
     &                        atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
     &                 * hybrid%basm1(:ng,j,0,itype)

              END DO
              hybrid%basm1(:ng,i,0,itype) = hybrid%basm1(:ng,i,0,itype)&
     &                            / sqrt(intgrf(hybrid%basm1(:ng,i,0,itype)**2,&
     &                              atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf))
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
     &                       atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

              IF(i.eq.j) THEN
                rdum1 = abs(1-rdum1)
                rdum  = rdum + rdum1**2
              ELSE
                rdum1 = abs(rdum1)
                rdum  = rdum + 2*rdum1**2
              END IF

              IF(rdum1.gt.1d-6) THEN
                  IF ( mpi%irank == 0 ) THEN
                    WRITE(6,'(A)') 'mixedbasis: Bad orthonormality of ' &
     &                //lchar(l)//'-product basis. Increase tolerance.'
                    WRITE(6,'(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of',&
     &                    rdum1,' occurred for (',i,',',j,')-overlap.'
                  END IF
                  STOP 'mixedbasis: Bad orthonormality of product basis.&
     &                  Increase tolerance.'
              END IF

            END DO
          END DO

          IF ( mpi%irank == 0 ) THEN
            WRITE(6,'(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')')&
     &                           lchar(l)//':',n,nn,sqrt(rdum)/nn
          END IF

          DEALLOCATE(olap,eigv,work,eig,ihelp)

        END DO !l
        IF ( mpi%irank == 0 ) WRITE(6,'(6X,A,I7)') 'Total:',&
     &                               sum(hybrid%nindxm1(0:hybrid%lcutm1(itype),itype))
      END DO ! itype

      hybrid%maxindxm1 = maxval(hybrid%nindxm1)

      ALLOCATE( basmhlp(atoms%jmtd,hybrid%maxindxm1,0:hybrid%maxlcutm1,atoms%ntype) )
      basmhlp(1:atoms%jmtd,1:hybrid%maxindxm1,0:hybrid%maxlcutm1,1:atoms%ntype)&
     &                   = hybrid%basm1(1:atoms%jmtd,1:hybrid%maxindxm1,0:hybrid%maxlcutm1,1:atoms%ntype)
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
     &                      'functions in each l-channel,',&
     &                      'such that they possess no multipolmoment'//&
     &                      ' except the last function:'

        WRITE(6,'(/,17x,A)') 'moment  (quality of orthonormality)'
      END IF
      DO itype = 1,atoms%ntype
        ng = atoms%jri(itype)

        IF(atoms%ntype.gt.1 .AND. mpi%irank==0)&
     &    WRITE(6,'(6X,A,I3)') 'Atom type',itype

        DO l = 0,hybrid%lcutm1(itype)
          ! determine radial function with the largest moment
          ! this function is used to build the linear combinations
          rdum = 0
          DO i = 1,hybrid%nindxm1(l,itype)
            rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
     &                       atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
            IF( abs(rdum1) .gt.rdum ) THEN
              n = i
              rdum = rdum1
            END IF
          END DO
          
          ! rearrange order of radial functions such that the last function possesses the largest moment
          j = 0
          bashlp(:ng) = hybrid%basm1(:ng,n,l,itype)
          DO i = 1,hybrid%nindxm1(l,itype)
            IF( i .eq. n ) CYCLE
            j = j + 1
            hybrid%basm1(:ng,j,l,itype) = hybrid%basm1(:ng,i,l,itype)
          END DO
          hybrid%basm1(:ng,hybrid%nindxm1(l,itype),l,itype) = bashlp(:ng)
          
        END DO
        
        
        DO l = 0,hybrid%lcutm1(itype)
          IF ( mpi%irank == 0 ) WRITE(6,'(6X,A)') lchar(l)//':'

          IF(hybrid%nindxm1(l,itype).eq.0) THEN
            IF( mpi%irank == 0 ) WRITE(6,'(6X,A,'':   0 ->    '')') lchar(l)
            CYCLE
          END IF
          
          n = hybrid%nindxm1(l,itype)
          DO i = 1,hybrid%nindxm1(l,itype)
            IF(i.eq.n) CYCLE
            ! calculate moment of radial function i
            rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
     &                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

            rdum  = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,n,l,itype),&
     &                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

            bashlp(:ng) = hybrid%basm1(:ng,n,l,itype)

            IF( sqrt(rdum**2 + rdum1**2) .le. 1E-06 .AND. mpi%irank == 0)&
     &                WRITE(6,*) 'Warning: Norm is smaller thann 1E-06!'

            ! change function n such that n is orthogonal to i
            ! since the functions basm1 have been orthogonal on input
            ! the linear combination does not destroy the orthogonality to the residual functions
            hybrid%basm1(:ng,n,l,itype) &
     &    = rdum /sqrt(rdum**2+rdum1**2) * bashlp(:ng)&
     &    + rdum1/sqrt(rdum**2+rdum1**2) * hybrid%basm1(:ng,i,l,itype)

            ! combine basis function i and n so that they possess no momemt
            hybrid%basm1(:ng,i,l,itype) &
     &    = rdum1/sqrt(rdum**2+rdum1**2) * bashlp(:ng)&
     &    - rdum /sqrt(rdum**2+rdum1**2) * hybrid%basm1(:ng,i,l,itype)

            rdum1 = intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,i,l,itype),&
     &                     atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

            IF( rdum1 .gt. 1E-10 )&
     &      STOP 'mixedbasis: moment of radial function does not vanish'

            IF ( mpi%irank == 0 ) WRITE(6,'(6x,I4,'' ->  '',ES8.1)') i,rdum1
          END DO

          ! test orthogonality
          rdum = 0
          DO i = 1,hybrid%nindxm1(l,itype)
            DO j = 1,hybrid%nindxm1(l,itype) 
              rdum1 = intgrf(hybrid%basm1(:ng,i,l,itype)*hybrid%basm1(:ng,j,l,itype),&
     &                       atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
              IF( i .ne. j ) THEN
                rdum = rdum + rdum1
              ELSE
                rdum = rdum + abs(1 - rdum1)
              END IF    
            END DO
          END DO
          IF ( mpi%irank == 0 )&
     &        WRITE(6,'(6x,I4,'' ->'',f10.5,''   ('',ES8.1,'' )'')')  n,&
     &           intgrf(atoms%rmsh(:ng,itype)**(l+1)*hybrid%basm1(:ng,n,l,itype),atoms%jri,&
     &                  atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf),rdum
        END DO
      END DO

      IF( exx ) THEN
        !
        ! - - - - - - Set up MT product basis for the representation of the response function  - - - - - -
        !

        IF ( mpi%irank == 0 ) THEN
          WRITE(6,*)
          WRITE(6,'(A)')       'MT product basis for response function:'
          WRITE(6,'(A)')       'Reduction due to overlap (quality of'//&
     &                         ' orthonormality, should be < 1.0E-06)'
        END IF

        hybrid%maxlcutm2  = maxval(hybrid%lcutm2)

        !
        ! read in atomic exchange potential
        !

        !
        ! to get the name of each element type read in inp file
        !
        OPEN (123,file='inp',form='formatted',status='old')
        DO
          READ(123,*) nchar
          IF(nchar .eq. '**') EXIT
        END DO

        DO itype = 1,atoms%ntype
          READ(123,*) noel(itype)
          IF( noel(itype)(2:2) .eq. " " ) THEN
            WRITE(fname(itype),*) 'vx-',noel(itype)(1:1),'.tab'
          ELSE
            WRITE(fname(itype),*) 'vx-',noel(itype),'.tab'
          END IF
          DO
            READ(123,*) nchar
            IF(nchar .eq. '**') EXIT
          END DO
        END DO
        CLOSE(123)

        CALL read_atompotential(atoms,&
     &                          (/(fname(i),i=1,atoms%ntype)/),potatom)

        ALLOCATE ( hybrid%nindxm2(0:hybrid%maxlcutm2,atoms%ntype) )
        ALLOCATE ( seleco(hybrid%maxindx,0:atoms%lmaxd) , selecu(hybrid%maxindx,0:atoms%lmaxd) )
        ALLOCATE ( selecmat(hybrid%maxindx,0:atoms%lmaxd,hybrid%maxindx,0:atoms%lmaxd) )


        ! determine maximal indices of (radial) mixed-basis functions (->nindxm) 
        ! (will be reduced later-on due to overlap)

        hybrid%nindxm2    = 0
        hybrid%maxindxp2  = 0
        DO itype=1,atoms%ntype
          seleco = .false.
          selecu = .false.
          seleco(1,0:hybrid%select2(1,itype)) = .true. !0!1
          selecu(1,0:hybrid%select2(3,itype)) = .true. !0!1
          seleco(2,0:hybrid%select2(2,itype)) = .true. !1!1
          selecu(2,0:hybrid%select2(4,itype)) = .true. !0!1
          DO l=0,hybrid%lcutm2(itype) 
            n = 0
            M = 0
!             !
!             ! core * valence products
!             !
!             DO l1 = 0,lmaxc(itype)
!               DO l2 = 0,lmax(itype)
!                 IF( l .lt. abs(l1-l2) .or. l .gt. l1+l2 ) CYCLE
!                 DO n1 = 1,nindxc(l1,itype)             
!                   DO n2 = 1,nindx(l2,itype)
!                     m = m + 1
!                     IF( .not. seleco(n2,l2) ) CYCLE
!                     n = n + 1
!                   END DO
!                 END DO
!               END DO
!             END DO

            ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
            selecmat = reshape ( (/ ((((&
     &                 seleco(n1,l1).and.selecu(n2,l2),&
     &                 n1=1,hybrid%maxindx),l1=0,atoms%lmaxd),           &
     &                 n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,&
     &                 (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

            DO l1=0,atoms%lmax(itype)
              DO l2=0,atoms%lmax(itype)
                IF(l.ge.abs(l1-l2).and.l.le.l1+l2) THEN
                  DO n1=1,hybrid%nindx(l1,itype)
                    DO n2=1,hybrid%nindx(l2,itype)
                      M = M + 1
                      IF(selecmat(n1,l1,n2,l2)) THEN
                        n = n + 1
                        selecmat(n2,l2,n1,l1) = .false. ! prevent double counting of products (a*b = b*a)
                      END IF
                    END DO !n2
                  END DO !n1
                END IF
              END DO !l2
            END DO !l1
            IF(n.eq.0 .AND. mpi%irank==0) &
     &         WRITE(6,'(A)') 'mixedbasis: Warning! No basis-function'//&
     &                        ' product of ' // lchar(l) //&
     &                        '-angular momentum defined.'
            hybrid%maxindxp2        = max(hybrid%maxindxp2,M)
            hybrid%nindxm2(l,itype) = n*dimension%jspd
          END DO
        END DO

        hybrid%maxindxm2 = maxval(hybrid%nindxm2)


        ALLOCATE (hybrid%basm2(atoms%jmtd,hybrid%maxindxm2,0:hybrid%maxlcutm2,atoms%ntype) )
        hybrid%basm2 = 0

        ! Define product bases and reduce them according to overlap
       
        DO itype=1,atoms%ntype
          seleco = .false.
          selecu = .false.
          seleco(1,0:hybrid%select2(1,itype)) = .true. !0!1
          selecu(1,0:hybrid%select2(3,itype)) = .true. !0!1
          seleco(2,0:hybrid%select2(2,itype)) = .true. !1!1
          selecu(2,0:hybrid%select2(4,itype)) = .true. !0!1
          IF(atoms%ntype.gt.1 .AND. mpi%irank==0)&
     &      WRITE(6,'(6X,A,I3)') 'Atom type',itype
          ng = atoms%jri(itype)
          DO l=0,hybrid%lcutm2(itype)
            n = hybrid%nindxm2(l,itype)
            ! allow for zero product-basis functions for
            ! current l-quantum number
            IF(n.eq.0) THEN
              IF ( mpi%irank == 0 ) WRITE(6,'(6X,A,'':   0 ->   0'')')&
     &                                                          lchar(l)
              CYCLE 
            END IF

            ! set up the overlap matrix
            ALLOCATE (olap(n,n),eigv(n,n),work(3*n),eig(n),ihelp(n))
            ihelp    = 1 ! initialize to avoid a segfault
            i        = 0

!             ! core*valence
!             DO l1 = 0,lmaxc(itype)
!               DO l2 = 0,lmax(itype)
!                 IF( l .lt. abs(l1-l2) .or. l .gt. l1+l2 ) CYCLE
!                 DO n1 = 1,nindxc(l1,itype)
!                   DO n2 = 1,nindx(l2,itype)
!                     IF( .not. seleco(n2,l2) ) CYCLE
!                     DO ispin = 1,jspd  
!                       i = i + 1
!                       IF(i.gt.n) THEN
!                          STOP 'mixedbasis:got too many product functions (bug).'
!                       END IF
!                       basm2(:ng,i,l,itype)  
!      &              = ( core1(:ng,n1,l1,itype,ispin) * bas1(:ng,n2,l2,itype,ispin) 
!      &                 +core2(:ng,n1,l1,itype,ispin) * bas2(:ng,n2,l2,itype,ispin) )/rmsh(:ng,itype)
!                     
!                       !normalize basm1
!                       rdum=sqrt(intgrf( basm1(:,i,l,itype)**2,
!      &                                  jri,jmtd,rmsh,dx,ntype,itype,gridf) ) 
! 
!                       basm2(:ng,i,l,itype)=basm2(:ng,i,l,itype)/rdum
! 
!                     END DO  !ispin
!                   END DO  !n2
!                 END DO  !n1
!               END DO  !l2
!             END DO  !l1



            ! Condense seleco and seleco into selecmat (each product corresponds to a matrix element)
            selecmat = reshape ( (/ ((((&
     &                 seleco(n1,l1).and.selecu(n2,l2),&
     &                 n1=1,hybrid%maxindx),l1=0,atoms%lmaxd),      &
     &                 n2=1,hybrid%maxindx),l2=0,atoms%lmaxd) /) ,  &
     &                 (/hybrid%maxindx,atoms%lmaxd+1,hybrid%maxindx,atoms%lmaxd+1/) )

            DO l1=0,atoms%lmax(itype)
              DO n1 = 1,hybrid%nindx(l1,itype)
                DO l2=0,atoms%lmax(itype)
                  DO n2 = 1,hybrid%nindx(l2,itype)
                
                    IF(l.ge.abs(l1-l2).and.l.le.l1+l2) THEN
                      IF(selecmat(n1,l1,n2,l2)) THEN
                        DO ispin=1,dimension%jspd
                          i = i + 1
                          IF(i.gt.n) STOP 'mixedbasis:&&
     &                    got too many product functions (bug).'

                          hybrid%basm2(:ng,i,l,itype) &
     &                  = ( bas1(:ng,n1,l1,itype,ispin)&
     &                      * bas1(:ng,n2,l2,itype,ispin)&
     &                  +   bas2(:ng,n1,l1,itype,ispin)&
     &                      * bas2(:ng,n2,l2,itype,ispin) )&
     &                    / atoms%rmsh(:ng,itype)

                          !normalize basm2
                          rdum=sqrt(intgrf(hybrid%basm2(:,i,l,itype)**2,&
     &                             atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )

                          hybrid%basm2(:ng,i,l,itype)=hybrid%basm2(:ng,i,l,itype)/rdum

                        END DO !ispin
                        ! prevent double counting of products (a*b = b*a)
                        selecmat(n2,l2,n1,l1) = .false. 
                      END IF
                    END IF
                  
                  END DO  !n2
                END DO  !l2
              END DO !n1
            END DO  !l1

            IF(i.ne.n) STOP 'mixedbasis: counting error for product&&
     &                       functions (bug).'


            !Add constant function and atomic EXX potential to l=0 basis
          
            IF(l.eq.0) THEN
            
              ! Check if basm2 must be reallocated
              IF(n+2.gt.size(hybrid%basm2,2)) THEN
                ALLOCATE ( basmhlp(atoms%jmtd,n,0:hybrid%maxlcutm2,atoms%ntype) )
                basmhlp(:,1:n,:,:) = hybrid%basm2
                DEALLOCATE (hybrid%basm2)
                ALLOCATE ( hybrid%basm2(atoms%jmtd,n+2,0:hybrid%maxlcutm2,atoms%ntype) )
                hybrid%basm2(:,1:n,:,:) = basmhlp(:,1:n,:,:) 
                DEALLOCATE ( basmhlp )
              END IF

              ! add normalized atomic exx potential
              hybrid%basm2(:ng,n+2,0,itype)= potatom(:ng,itype)*atoms%rmsh(:ng,itype)

              rdum=intgrf(hybrid%basm2(:ng,n+2,0,itype)*hybrid%basm2(:ng,n+2,0,itype),&
     &                    atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

              hybrid%basm2(:ng,n+2,0,itype) = hybrid%basm2(:ng,n+2,0,itype)/sqrt(rdum)
            
            
              ! add normalized constant function
              hybrid%basm2(:ng,n+1,0,itype) = atoms%rmsh(:ng,itype)&
     &                               / sqrt(atoms%rmsh(ng,itype)**3/3)
            
              ! orthogonalize the constant function with respect the atomic potential
              hybrid%basm2(:ng,n+1,0,itype) = hybrid%basm2(:ng,n+1,0,itype)&
     &          - intgrf(hybrid%basm2(:ng,n+1,0,itype)* hybrid%basm2(:ng,n+2,0,itype),&
     &                   atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)&
     &            * hybrid%basm2(:ng,n+2,0,itype)

              
              hybrid%basm2(:ng,n+1,0,itype) = hybrid%basm2(:ng,n+1,0,itype)&
     &              / sqrt( intgrf(hybrid%basm2(:ng,n+1,0,itype)**2,&
     &                             atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )
            
              n               = n + 2
              DEALLOCATE( olap )
              ALLOCATE  ( olap(n,n) )
              hybrid%nindxm2(l,itype)= n
            END IF

      
            ! Calculate overlap
            olap = 0
            DO n2=1,n
              DO n1=1,n2
                olap(n1,n2) = &
     &                   intgrf(hybrid%basm2(:,n1,l,itype)*hybrid%basm2(:,n2,l,itype),&
     &                          atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
                olap(n2,n1) = olap(n1,n2)
              END DO
            END DO

            IF( l .eq. 0 ) THEN
         

              ! orthogonalize all function ( 1:n-2 ) with respect to the constant and atomic function
              DO i = n-2,1,-1
            
                DO j = n-1,n
                  hybrid%basm2(:ng,i,0,itype) = hybrid%basm2(:ng,i,0,itype)&
     &                                 - olap(j,i)*hybrid%basm2(:ng,j,0,itype)
                END DO
                hybrid%basm2(:ng,i,0,itype) = hybrid%basm2(:ng,i,0,itype) &
     &              / sqrt( intgrf(hybrid%basm2(:ng,i,0,itype)**2,&
     &                             atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf) )
              END DO
  
              ! Calculate overlap
              olap = 0
              DO n2=1,n
                DO n1=1,n2
                  olap(n1,n2) = &
     &               intgrf(hybrid%basm2(:,n1,l,itype)*hybrid%basm2(:,n2,l,itype),&
     &                      atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
                  olap(n2,n1) = olap(n1,n2)
                END DO
              END DO

              n = n - 2

            END IF

            ! In order to get ride of the linear dependencies in the
            ! radial functions basm2 belonging to fixed l and itype
            ! the overlap matrix is diagonalized and those eigenvectors
            ! with a eigenvalue greater then tolerance are retained


            ! Diagonalize

            CALL diagonalize(eigv,eig,olap(:n,:n))

            IF( l .eq. 0 ) THEN
              OPEN (313,file='eig_basis',form='formatted')
              WRITE(313,'(f30.10)') eig
              CLOSE(313)
            END IF

            ! Get rid of linear dependencies (eigenvalue <= tolerance)

            IF( .false. ) THEN
              IF( l .eq. 0 ) rdum2 = hybrid%tolerance2*1E07
              ic = int(rdum2/10.**(8-l))
              rdum2 = rdum2 - ic*10.**(8-l)
              WRITE(*,*) ic
              nn = 0
              DO i= n,1,-1
                nn        = nn + 1
                ihelp(nn) = i
                WRITE(*,*) nn,i
                IF( nn .eq. ic ) EXIT
              END DO

            ELSE
              nn = 0
              DO i=1,n !hybrid%nindxm2(l,itype)
                IF(eig(i).gt.hybrid%tolerance2) THEN 
                  nn        = nn + 1
                  ihelp(nn) = i
                END IF
              END DO
            END IF
            
            hybrid%nindxm2(l,itype) = nn
            eig              = eig(ihelp)
            eigv(:,:)        = eigv(:,ihelp)

            DO i=1,ng
              hybrid%basm2(i,1:nn,l,itype) = matmul(hybrid%basm2(i,1:n,l,itype),&
     &                                    eigv(:,1:nn)) / sqrt(eig(:nn))
            END DO

            IF( l .eq. 0 ) THEN
              hybrid%nindxm2(l,itype)         = nn + 2
              hybrid%basm2(:,nn+1,l,itype) = hybrid%basm2(:,n+1,l,itype)
              hybrid%basm2(:,nn+2,l,itype) = hybrid%basm2(:,n+2,l,itype)
              nn = nn + 2
            END IF

 
            ! Check orthonormality of product basis
            rdum = 0
            DO i=1,nn
              DO j=1,i
                rdum1 = intgrf(hybrid%basm2(:,i,l,itype)*hybrid%basm2(:,j,l,itype),&
     &                         atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)

                IF(i.eq.j) THEN
                  rdum1 = abs(1-rdum1)
                  rdum  = rdum + rdum1**2
                ELSE
                  rdum1 = abs(rdum1)
                  rdum  = rdum + 2*rdum1**2
                ENDIF

                IF(rdum1.gt.1d-6) THEN
                  IF ( mpi%irank == 0 ) THEN
                    WRITE(6,'(A)') 'mixedbasis: Bad orthonormality of ' &
     &               //lchar(l)//'-product basis. Increase tolerance.'
                    WRITE(6,'(12X,A,F9.6,A,2(I3.3,A))') 'Deviation of',&
     &                  rdum1,' occurred for (',i,',',j,')-overlap.'
                  END IF
                  STOP 'mixedbasis: Bad orthonormality of product basis.&
     &                Increase tolerance.'
                END IF

              END DO
            END DO

            IF ( mpi%irank == 0 )&
     &             WRITE(6,'(6X,A,I4,'' ->'',I4,''   ('',ES8.1,'' )'')')&
     &                                  lchar(l)//':',n,nn,sqrt(rdum)/nn

            DEALLOCATE(olap,eigv,work,eig,ihelp)

          END DO !l
          IF ( mpi%irank == 0 ) WRITE(6,'(6X,A,I7)')&
     &                      'Total:',sum(hybrid%nindxm2(0:hybrid%lcutm2(itype),itype))
        END DO ! itype

        hybrid%maxindxm2 = maxval(hybrid%nindxm2)

        ALLOCATE( basmhlp(atoms%jmtd,hybrid%maxindxm2,0:hybrid%maxlcutm2,atoms%ntype) )
        basmhlp(1:atoms%jmtd,1:hybrid%maxindxm2,0:hybrid%maxlcutm2,1:atoms%ntype)&
     &    = hybrid%basm2(1:atoms%jmtd,1:hybrid%maxindxm2,0:hybrid%maxlcutm2,1:atoms%ntype)
        DEALLOCATE(hybrid%basm2)
        ALLOCATE( hybrid%basm2(atoms%jmtd,hybrid%maxindxm2,0:hybrid%maxlcutm2,atoms%ntype) )
        hybrid%basm2=basmhlp
        DEALLOCATE(basmhlp)
      END IF

      
      DO itype = 1,atoms%ntype
        IF ( ANY(hybrid%nindxm1(0:hybrid%lcutm1(itype),itype) == 0) )&
     &    STOP 'mixedbasis: any hybrid%nindxm1 eq 0'
      END DO

!       ic = 0
!       DO itype = 1,ntype
!         DO l = 0,lcutm2(itype)
!           DO i = 1,nindxm2(l,itype)
!             ic = ic + 1
!             DO j = 1,jmtd
!               WRITE(200+ic,'(2f15.10)') rmsh(j,itype),basm2(j,i,l,itype,1)
!             END DO
!           END DO
!         END DO
!       END DO
!       
!       STOP

      !
      ! Write all information to a file to enable restarting
      !
      IF ( l_restart .AND. mpi%irank == 0 ) THEN

        OPEN(UNIT=iounit,FILE=ioname,FORM='unformatted',&
     &       STATUS='replace')

        WRITE(iounit) kpts%nkptf,hybrid%gptmd
        WRITE(iounit) hybrid%maxgptm,hybrid%maxindx
        WRITE(iounit) hybrid%ngptm,hybrid%gptm,hybrid%pgptm,hybrid%nindx

        WRITE(iounit) hybrid%maxgptm1,hybrid%maxlcutm1,hybrid%maxindxm1,hybrid%maxindxp1
        WRITE(iounit) hybrid%ngptm1,hybrid%pgptm1,hybrid%nindxm1
        WRITE(iounit) hybrid%basm1

        IF ( exx ) THEN
          WRITE(iounit) hybrid%maxgptm2,hybrid%maxlcutm2,hybrid%maxindxm2,hybrid%maxindxp2
          WRITE(iounit) hybrid%ngptm2,hybrid%pgptm2,hybrid%nindxm2
          WRITE(iounit) hybrid%basm2
        END IF

        CLOSE(iounit)

        l_restart = .false.

      END IF

      END SUBROUTINE mixedbasis


      SUBROUTINE read_atompotential(atoms,fname,&
     &                              potential)

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
          IF( rdum .gt. (atoms%rmt(itype)+0.2) ) EXIT
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
            IF( r_atom(j) .gt. r  ) THEN
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

      END SUBROUTINE

      END MODULE m_mixedbasis

