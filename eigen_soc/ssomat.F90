MODULE m_ssomat
  USE m_judft
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ssomat(seigvso,theta,phi,eig_id,DIMENSION,atoms,kpts,sym,&
       cell,noco, input,mpi, oneD,enpara,v,results )
    
    USE m_eig66_io
    USE m_spnorb 
    USE m_abcof 
    USE m_types_mat
    USE m_types_setup
    USE m_types_mpi
    USE m_types_enpara
    USE m_types_potden
    USE m_types_misc
    USE m_types_kpts
    USE m_types_tlmplm
    USE m_types_usdus
    USE m_types_lapw
    USE m_constants
    IMPLICIT NONE

    TYPE(t_mpi),INTENT(IN)         :: mpi
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_oneD),INTENT(IN)        :: oneD
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_enpara),INTENT(IN)      :: enpara
    TYPE(t_potden),INTENT(IN)      :: v
    TYPE(t_results),INTENT(IN)     :: results
    INTEGER,INTENT(IN)             :: eig_id
    REAL,INTENT(in)                :: theta(:),phi(:) ! more than a single angle at once...
    REAL,INTENT(OUT)               :: seigvso(:)
    !     ..
    !     .. Locals ..
#ifdef CPP_MPI
    INTEGER:: ierr
    include 'mpif.h'
#endif
    INTEGER :: neigf=1  !not full-matrix
    INTEGER :: ilo,js,jsloc,nk,n,l ,lm,band,nr,ne,nat,m
    INTEGER :: na 
    REAL    :: r1,r2  
    COMPLEX :: c1,c2  

    COMPLEX, ALLOCATABLE :: matel(:,:,:) 
    REAL,    ALLOCATABLE :: eig_shift(:,:,:) 

    COMPLEX, ALLOCATABLE :: acof(:,:,:,:,:), bcof(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: ccof(:,:,:,:,:,:)
    COMPLEX,ALLOCATABLE  :: soangl(:,:,:,:,:,:,:)
  
    TYPE(t_rsoc) :: rsoc
    TYPE(t_mat)  :: zmat
    TYPE(t_usdus):: usdus
    TYPE(t_lapw) :: lapw

    IF (ANY(atoms%neq/=1)) CALL judft_error('(spin spiral + soc) does not work'//&
         ' properly for more than one atom per type!',calledby="ssomat")


    
    ! needed directly for calculating matrix elements  
    seigvso=0.0
    ALLOCATE(eig_shift(DIMENSION%neigd,kpts%nkpt,SIZE(theta)));eig_shift=0.0
    ALLOCATE( acof(dimension%neigd,0:dimension%lmd,atoms%nat,2,2),&
         bcof(dimension%neigd,0:dimension%lmd,atoms%nat,2,2) )
    ALLOCATE( ccof(-atoms%llod:atoms%llod,dimension%neigd,atoms%nlod,atoms%nat,2,2) )

    ALLOCATE( matel(neigf,DIMENSION%neigd,0:atoms%ntype) )

  

    CALL usdus%init(atoms,2)


    !Calculate radial and angular matrix elements of SOC
    !many directions of SOC at once...
    CALL spnorb(atoms,noco,input,mpi, enpara, v%mt, usdus, rsoc,.FALSE.)

    ALLOCATE(soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,&
         atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,SIZE(theta)))
    DO nr=1,SIZE(theta)
       CALL spnorb_angles(atoms,mpi,theta(nr),phi(nr),soangl(:,:,:,:,:,:,nr))
    ENDDO
    
    DO nk=mpi%irank+1,kpts%nkpt,mpi%isize
       CALL lapw%init(input,noco, kpts,atoms,sym,nk,cell,.false., mpi)
       zMat%matsize1=lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot
       zmat%matsize2=DIMENSION%neigd
       zmat%l_real=.FALSE.
       IF (ALLOCATED(zmat%data_c)) DEALLOCATE(zmat%data_c)
       ALLOCATE(zmat%data_c(zMat%matsize1,zmat%matsize2))
       CALL read_eig(eig_id,nk,1,neig=ne,eig=eig_shift(:,nk,1),zmat=zmat)
       DO jsloc= 1,2 
          eig_shift(:,nk,1)=0.0 !not needed
          CALL abcof(input,atoms,sym, cell,lapw,ne,usdus,noco,jsloc,oneD, &
               acof(:,:,:,jsloc,1),bcof(:,:,:,jsloc,1),ccof(:,:,:,:,jsloc,1),zMat)
       ENDDO

       ! rotate abcof into global spin coordinate frame 
       nat= 0 
       DO n= 1,atoms%ntype  
          DO na= 1,atoms%neq(n)
             nat= nat+1 
             r1= noco%alph(n) 
             r2= noco%beta(n) 
             DO lm= 0,DIMENSION%lmd
                DO band= 1,DIMENSION%neigd
                   c1= acof(band,lm,nat,1,1)
                   c2= acof(band,lm,nat,2,1)
                   acof(band,lm,nat,1,1)= CMPLX(COS(r1/2.),-SIN(r1/2.)) *CMPLX( COS(r2/2.),0.) *c1
                   acof(band,lm,nat,2,1)= CMPLX(COS(r1/2.),-SIN(r1/2.)) *CMPLX(-SIN(r2/2.),0.) *c2
                   acof(band,lm,nat,1,2)= CMPLX(COS(r1/2.),+SIN(r1/2.)) *CMPLX(+SIN(r2/2.),0.) *c1
                   acof(band,lm,nat,2,2)= CMPLX(COS(r1/2.),+SIN(r1/2.)) *CMPLX( COS(r2/2.),0.) *c2
                   c1= bcof(band,lm,nat,1,1)
                   c2= bcof(band,lm,nat,2,1)
                   bcof(band,lm,nat,1,1)= CMPLX(COS(r1/2.),-SIN(r1/2.)) *CMPLX( COS(r2/2.),0.) *c1
                   bcof(band,lm,nat,2,1)= CMPLX(COS(r1/2.),-SIN(r1/2.)) *CMPLX(-SIN(r2/2.),0.) *c2
                   bcof(band,lm,nat,1,2)= CMPLX(COS(r1/2.),+SIN(r1/2.)) *CMPLX(+SIN(r2/2.),0.) *c1
                   bcof(band,lm,nat,2,2)= CMPLX(COS(r1/2.),+SIN(r1/2.)) *CMPLX( COS(r2/2.),0.) *c2
                ENDDO ! band
             ENDDO   ! lm
             DO ilo = 1,atoms%nlo(n)
                l = atoms%llo(ilo,n)
                DO band= 1,DIMENSION%neigd
                   DO m = -l, l
                      c1= ccof(m,band,ilo,nat,1,1)
                      c2= ccof(m,band,ilo,nat,2,1)
                      ccof(m,band,ilo,nat,1,1)= CMPLX(COS(r1/2.),-SIN(r1/2.))*CMPLX( COS(r2/2.),0.)*c1
                      ccof(m,band,ilo,nat,2,1)= CMPLX(COS(r1/2.),-SIN(r1/2.))*CMPLX(-SIN(r2/2.),0.)*c2
                      ccof(m,band,ilo,nat,1,2)= CMPLX(COS(r1/2.),+SIN(r1/2.))*CMPLX(+SIN(r2/2.),0.)*c1
                      ccof(m,band,ilo,nat,2,2)= CMPLX(COS(r1/2.),+SIN(r1/2.))*CMPLX( COS(r2/2.),0.)*c2
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO nr=1,size(theta) !loop over angles
          ! matrix elements within k
          CALL ssomatel(neigf,dimension,atoms, noco, &
               soangl(:,:,:,:,:,:,nr),rsoc%rsopp(:,:,:,:),rsoc%rsoppd(:,:,:,:),&
               rsoc%rsopdp(:,:,:,:),rsoc%rsopdpd(:,:,:,:),rsoc%rsoplop(:,:,:,:), &
               rsoc%rsoplopd(:,:,:,:),rsoc%rsopdplo(:,:,:,:),rsoc%rsopplo(:,:,:,:),&
               rsoc%rsoploplop(:,:,:,:,:),&
               .TRUE.,&
               acof,bcof, ccof,&
               acof,bcof, ccof,&
               matel )
          eig_shift(:,nk,nr)=matel(1,:,0)
       ENDDO
    ENDDO

    !Collect data from distributed k-loop
#ifdef CPP_MPI
    IF (mpi%irank==0) THEN
       CALL MPI_REDUCE(MPI_IN_PLACE,eig_shift,SIZE(eig_shift),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi%mpi_comm,ierr)
    ELSE
       CALL MPI_REDUCE(eig_shift,eig_shift,SIZE(eig_shift),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi%mpi_comm,ierr)
    ENDIF
#endif

    IF (mpi%irank==0) THEN
       !Sum all shift using weights
       DO nr=1,SIZE(theta)
          DO nk=1,kpts%nkpt
             seigvso(nr)=seigvso(nr)+dot_PRODUCT(results%w_iks(:,nk,1),eig_shift(:,nk,nr))
          ENDDO
       ENDDO
       seigvso= results%seigv+seigvso
    ENDIF
  END SUBROUTINE ssomat

  ! ==================================================================== ! 

  SUBROUTINE ssomatel(neigf,dimension,atoms, noco,&
       soangl,rsopp,rsoppd,rsopdp,rsopdpd,rsoplop,&
       rsoplopd,rsopdplo,rsopplo,rsoploplop,&
       diag, &
       acof1,bcof1,ccof1,acof2,bcof2,ccof2,&
       matel )
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_atoms),INTENT(IN)       :: atoms

    LOGICAL, INTENT(IN)  :: diag 
    INTEGER, INTENT(IN)  :: neigf   
    REAL,    INTENT(IN)  :: &
         rsopp(:,:,:,:), rsoppd(:,:,:,:),&
         rsopdp(:,:,:,:), rsopdpd(:,:,:,:),  &
         rsoplop(:,:,:,:),rsoplopd(:,:,:,:),&
         rsopdplo(:,:,:,:),rsopplo(:,:,:,:),&
         rsoploplop(:,:,:,:,:)
    COMPLEX, INTENT(IN)  :: &
         soangl(:,-atoms%lmaxd:,:,:,-atoms%lmaxd:,:),  &
         acof1(:,0:,:,:,:), &
         bcof1(:,0:,:,:,:),&
         ccof1(-atoms%llod:,:,:,:,:,:),&
         acof2(:,0:,:,:,:), &
         bcof2(:,0:,:,:,:),&
         ccof2(-atoms%llod:,:,:,:,:,:)
 
    Complex, INTENT(OUT) :: matel(neigf,dimension%neigd,0:atoms%ntype)

    INTEGER :: band1,band2,bandf, n ,na, l,m1,m2,lm1,lm2,&
         jsloc1,jsloc2, js1,js2,jsnumber,ilo,ilop,nat
    COMPLEX, ALLOCATABLE :: sa(:,:),sb(:,:),sc(:,:,:),ral(:,:,:)
    COMPLEX, ALLOCATABLE :: ra(:,:),rb(:,:),rc(:,:,:),rbl(:,:,:)

    ! with the following nesting of loops the calculation of the 
    ! matrix-elements is of order
    ! natall*lmd*neigd*(lmd+neigd) ; note that  lmd+neigd << lmd*neigd

    matel(:,:,:)= CMPLX(0.,0.) 
    ALLOCATE ( sa(2,0:dimension%lmd),sb(2,0:dimension%lmd),ra(2,0:dimension%lmd),rb(2,0:dimension%lmd) )
    ALLOCATE ( sc(2,-atoms%llod:atoms%llod,atoms%nlod),rc(2,-atoms%llod:atoms%llod,atoms%nlod) )
    ALLOCATE ( ral(2,-atoms%llod:atoms%llod,atoms%nlod),rbl(2,-atoms%llod:atoms%llod,atoms%nlod) )

    ! within one k-point loop over global spin 
    IF (diag) THEN 
       jsnumber= 2
    ELSE 
       jsnumber= 1
    ENDIF
    DO js2= 1,jsnumber 
       IF (diag) THEN
          js1= js2
       ELSE
          js1= 2
       ENDIF

       ! loop over MT 
       na= 0 
       DO n= 1,atoms%ntype  
          DO nat= 1,atoms%neq(n) 
             na= na+1 

             DO band2= 1,dimension%neigd ! loop over eigenstates 2

                DO l= 1,atoms%lmax(n) ! loop over l
                   DO m1= -l,l   ! loop over m1
                      lm1= l*(l+1) + m1

                      DO jsloc2= 1,2
                         sa(jsloc2,lm1) = CMPLX(0.,0.)
                         sb(jsloc2,lm1) = CMPLX(0.,0.)
                         DO m2= -l,l
                            lm2= l*(l+1) + m2

                            sa(jsloc2,lm1)= sa(jsloc2,lm1) + &
                                 CONJG(acof2(band2,lm2,na,jsloc2,js2))&
                                 * soangl(l,m2,js2,l,m1,js1)
                            sb(jsloc2,lm1)= sb(jsloc2,lm1) + &
                                 CONJG(bcof2(band2,lm2,na,jsloc2,js2))&
                                 * soangl(l,m2,js2,l,m1,js1)

                         ENDDO ! m2  
                      ENDDO   ! jsloc2

                   ENDDO ! m1
                ENDDO   ! l

                DO ilo = 1, atoms%nlo(n) ! LO-part
                   l = atoms%llo(ilo,n)
                   DO m1 = -l, l
                      DO jsloc2= 1,2
                         sc(jsloc2,m1,ilo) = CMPLX(0.,0.)
                         DO m2= -l, l
                            sc(jsloc2,m1,ilo) = sc(jsloc2,m1,ilo) +&
                                 CONJG(ccof2(m2,band2,ilo,na,jsloc2,js2))&
                                 * soangl(l,m2,js2,l,m1,js1)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO ! ilo

                DO l= 1,atoms%lmax(n) ! loop over l
                   DO m1= -l,l   ! loop over m1
                      lm1= l*(l+1) + m1

                      DO jsloc1= 1,2
                         ra(jsloc1,lm1)= CMPLX(0.,0.) 
                         rb(jsloc1,lm1)= CMPLX(0.,0.) 
                         DO jsloc2= 1,2
                            ra(jsloc1,lm1)= ra(jsloc1,lm1) +  &
                                 sa(jsloc2,lm1) * rsopp(n,l,jsloc1,jsloc2) &
                                 + sb(jsloc2,lm1) * rsoppd(n,l,jsloc1,jsloc2) 
                            rb(jsloc1,lm1)= rb(jsloc1,lm1) +&
                                 sa(jsloc2,lm1) * rsopdp(n,l,jsloc1,jsloc2)&
                                 + sb(jsloc2,lm1) * rsopdpd(n,l,jsloc1,jsloc2)
                         ENDDO ! jsloc2
                      ENDDO   ! jsloc1 

                   ENDDO ! m1
                ENDDO   ! l

                DO ilo = 1, atoms%nlo(n) ! LO-part
                   l = atoms%llo(ilo,n)
                   DO m1 = -l, l
                      lm1= l*(l+1) + m1
                      DO jsloc1= 1,2
                         ral(jsloc1,m1,ilo) = CMPLX(0.,0.)
                         rbl(jsloc1,m1,ilo) = CMPLX(0.,0.)
                         rc(jsloc1,m1,ilo)  = CMPLX(0.,0.)
                         DO jsloc2= 1,2
                            ral(jsloc1,m1,ilo) = ral(jsloc1,m1,ilo) +&
                                 sc(jsloc2,m1,ilo) * rsopplo(n,ilo,jsloc1,jsloc2)
                            rbl(jsloc1,m1,ilo) = rbl(jsloc1,m1,ilo) +&
                                 sc(jsloc2,m1,ilo) * rsopdplo(n,ilo,jsloc1,jsloc2)
                            rc(jsloc1,m1,ilo) = rc(jsloc1,m1,ilo) +&
                                 sa(jsloc2,lm1) * rsoplop(n,ilo,jsloc1,jsloc2)&
                                 + sb(jsloc2,lm1) * rsoplopd(n,ilo,jsloc1,jsloc2)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO ! ilo

                DO l= 1,atoms%lmax(n) ! loop over l
                   DO m1= -l,l   ! loop over m1
                      lm1= l*(l+1) + m1

                      DO jsloc1= 1,2
                         DO bandf= 1,neigf
                            IF (neigf==dimension%neigd) THEN
                               band1= bandf
                            ELSE
                               band1= band2
                            ENDIF
                            matel(bandf,band2,n)= matel(bandf,band2,n) +&
                                 acof1(band1,lm1,na,jsloc1,js1)*ra(jsloc1,lm1)   &
                                 + bcof1(band1,lm1,na,jsloc1,js1)*rb(jsloc1,lm1)   
                         ENDDO ! band1
                      ENDDO   ! jsloc1 

                   ENDDO ! m1,lm1 
                ENDDO   ! l

                DO ilo = 1, atoms%nlo(n) ! LO-part
                   l = atoms%llo(ilo,n)
                   DO m1 = -l, l
                      lm1= l*(l+1) + m1

                      DO jsloc1= 1,2
                         DO bandf= 1,neigf
                            IF (neigf==dimension%neigd) THEN
                               band1= bandf
                            ELSE
                               band1= band2
                            ENDIF
                            matel(bandf,band2,n)= matel(bandf,band2,n) +&
                                 ccof1(m1,band1,ilo,na,jsloc1,js1)*rc(jsloc1,m1,ilo)&
                                 + acof1(band1,lm1,na,jsloc1,js1)*ral(jsloc1,m1,ilo)&
                                 + bcof1(band1,lm1,na,jsloc1,js1)*rbl(jsloc1,m1,ilo)
                         ENDDO ! band1
                      ENDDO   ! jsloc1 

                      DO ilop = 1,atoms%nlo(n)
                         IF (atoms%llo(ilop,n).EQ.l) THEN
                            DO jsloc1= 1,2
                               DO bandf= 1,neigf
                                  IF (neigf==dimension%neigd) THEN
                                     band1= bandf
                                  ELSE
                                     band1= band2
                                  ENDIF
                                  DO jsloc2= 1,2 
                                     matel(bandf,band2,n)= matel(bandf,band2,n) +&
                                          ccof1(m1,band1,ilo,na,jsloc1,js1)*&
                                          rsoploplop(n,ilo,ilop,jsloc1,jsloc2)*&
                                          sc(jsloc2,m1,ilop)
                                  ENDDO   ! jsloc2
                               ENDDO     ! band1
                            ENDDO   ! jsloc1
                         ENDIF
                      ENDDO ! ilop

                   ENDDO   ! m1 
                ENDDO     ! ilo

             ENDDO     ! band2
          ENDDO       ! nat,na 
       ENDDO         ! n        
    ENDDO           ! js2,js1

    DO n= 1,atoms%ntype 
          DO band2= 1,dimension%neigd 
             DO bandf= 1,neigf 
                matel(bandf,band2,0)= matel(bandf,band2,0) + matel(bandf,band2,n) 
             ENDDO
          ENDDO
    ENDDO

    IF (diag) THEN 
       DO n= 1,atoms%ntype
          DO band2= 1,dimension%neigd
             IF (neigf==dimension%neigd) THEN
                bandf= band2 
             ELSE
                bandf= 1 
             ENDIF
             IF (ABS(AIMAG(matel(bandf,band2,n)))>1.e-10) THEN
                PRINT *,bandf,band2,n,AIMAG(matel(bandf,band2,n))
                CALL judft_error('Stop in ssomatel:  diagonal matrix element not real')
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    DEALLOCATE ( sa,sb,ra,rb )

  END SUBROUTINE ssomatel
END MODULE m_ssomat
