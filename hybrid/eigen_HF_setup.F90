
module m_eigen_hf_setup
  contains
    subroutine eigen_hf_setup(hybrid,input,sym,kpts,dimension,atoms,mpi,noco,cell,oneD,results,jsp,eig_id_hf,&
         hybdat,irank2,it,vr0)
    use m_types
    use m_eig66_io
    use m_util
    use m_apws
    use m_checkolap
    use m_read_core
    use m_gen_wavf
    implicit none
    TYPE(t_hybrid),INTENT(INOUT)   :: hybrid
    TYPE(t_kpts),INTENT(IN)     :: kpts
    TYPE(t_dimension),INTENT(IN):: dimension
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_results),INTENT(IN)  :: results
    INTEGER,INTENT(IN)          :: irank2(:),it
 
    INTEGER,INTENT(IN)          :: jsp,eig_id_hf
    REAL, INTENT(IN)            :: vr0(:,:,:)
    TYPE(t_hybdat),INTENT(INOUT) ::hybdat

     
    INTEGER:: ok,nk,nrec1,i,j,bands,ll,l1,l2,ng,itype,n,l,n1,n2,nn
    

    TYPE(t_zmat),allocatable :: zmat(:)
      REAL,    ALLOCATABLE    ::  eig_irr(:,:)
      REAL,    ALLOCATABLE    ::  basprod(:)
      REAL                    ::  el_eig(0:atoms%lmaxd,atoms%ntype), ello_eig(atoms%nlod,atoms%ntype),bk(3)
      INTEGER                 ::  degenerat(dimension%neigd2+1,kpts%nkpt)
      INTEGER                 :: matind(dimension%nbasfcn,2),nred
      TYPE(t_lapw):: lapw

      LOGICAL:: skip_kpt(kpts%nkpt)
      REAL   :: g(3)
      skip_kpt=.false.
      
    IF( hybrid%l_calhf ) THEN
!
       !             preparations for HF and hybrid functional calculation
!
       CALL timestart("gen_bz and gen_wavf")
       STOP "TODO allocate"
       !ALLOCATE ( zmat%z_c(dimension%nbasfcn,dimension%neigd2,kpts%nkpt),stat=ok )
       !zmat(:)%z_c = 0
       IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation z_c'
       ALLOCATE ( eig_irr(dimension%neigd2,kpts%nkpt)      ,stat=ok )
       IF( ok .ne. 0 ) STOP'eigen_hf: failure allocation eig_irr'
       ALLOCATE ( hybdat%kveclo_eig(atoms%nlotot,kpts%nkpt)  ,stat=ok )
       IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%kveclo_eig'
       eig_irr = 0 ; hybdat%kveclo_eig = 0

              ! Reading the eig file
              DO nk = 1,kpts%nkpt
#               ifdef CPP_MPI
                  ! jump to next k-point if this process is not present in communicator
                  IF ( skip_kpt(nk) ) CYCLE
#               endif
                  nrec1 = kpts%nkpt*(jsp-1) + nk
                 CALL read_eig(eig_id_hf,nk,jsp,el=el_eig,ello=ello_eig, neig=hybdat%ne_eig(nk),eig=eig_irr(:,nk), kveclo=hybdat%kveclo_eig(:,nk)) !TODO introduce zmat!!,z=z_irr(:,:,nk))
                

              END DO

              !
              !determine degenerate states at each k-point
              !
              ! degenerat(i) =1  band i  is not degenerat ,
              ! degenerat(i) =j  band i  has j-1 degenart states ( i, i+1, ..., i+j)
              ! degenerat(i) =0  band i  is  degenerat, but is not the lowest band
              !                  of the group of degenerate states
              IF ( mpi%irank == 0 ) THEN
                WRITE(6,*)
                WRITE(6,'(A)') "   k-point      |   number of occupied&
      &bands  |   maximal number of bands"
              END IF
              degenerat = 1
              hybdat%nobd      = 0
              DO nk=1 ,kpts%nkpt
#               ifdef CPP_MPI
                  ! jump to next k-point if this k-point is not treated at this process
                  IF ( skip_kpt(nk) ) CYCLE
#               endif
                DO i=1,hybdat%ne_eig(nk)
                  DO j=i+1,hybdat%ne_eig(nk)
                    IF( abs(eig_irr(i,nk)-eig_irr(j,nk)) < 1E-07) THEN !0.015
                      degenerat(i,nk) = degenerat(i,nk) + 1
                    END IF
                  END DO
                END DO

                DO i=1,hybdat%ne_eig(nk)
                  IF( degenerat(i,nk) .ne. 1 .or. degenerat(i,nk) .ne. 0 ) &
                    degenerat(i+1:i+degenerat(i,nk)-1,nk) = 0
                END DO


                ! set the size of the exchange matrix in the space of the wavefunctions

                hybdat%nbands(nk)=bands
                IF(hybdat%nbands(nk).gt.hybdat%ne_eig(nk)) THEN
                  IF ( mpi%irank == 0 ) THEN
                    WRITE(*,*) ' maximum for hybdat%nbands is', hybdat%ne_eig(nk)
                    WRITE(*,*) ' increase energy window to obtain enough eigenvalues'
                    WRITE(*,*) ' set hybdat%nbands equal to hybdat%ne_eig'
                  END IF
                  hybdat%nbands(nk)=hybdat%ne_eig(nk)
                END IF

                DO i = hybdat%nbands(nk)-1,1,-1
                  IF( (degenerat(i,nk) .ge. 1) .and. (degenerat(i,nk)+i-1 .ne. hybdat%nbands(nk) ) ) THEN
                    hybdat%nbands(nk) = i + degenerat(i,nk) - 1
                    EXIT
                  END IF
                END DO

                DO i = 1,hybdat%ne_eig(nk)
                  IF( results%w_iks(i,nk,jsp) .gt. 0d0 ) hybdat%nobd(nk) = hybdat%nobd(nk) + 1
                 
                END DO


              END DO

#             ifdef CPP_MPI
                ! send results for occupied bands to all processes
                sndreqd = 0 ; rcvreqd = 0
                DO nk = 1,kpts%nkpt
                  IF ( skip_kpt(nk) ) THEN
                    rcvreqd = rcvreqd + 1
                    CALL MPI_IRECV(hybdat%nobd(nk),1,MPI_INTEGER4, MPI_ANY_SOURCE,TAG_SNDRCV_HYBDAT%NOBD+nk, mpi,rcvreq(rcvreqd),ierr(1))
                  ELSE
                    i = mod( mpi%irank + isize2(nk), mpi%isize )
                    DO WHILE ( i < mpi%irank-irank2(nk) .OR. i >= mpi%irank-irank2(nk)+isize2(nk) )
                      sndreqd = sndreqd + 1
                      CALL MPI_ISSEND(hybdat%nobd(nk),1,MPI_INTEGER4,i, TAG_SNDRCV_HYBDAT%NOBD+nk,mpi, sndreq(sndreqd),ierr(1) )
                      i = mod( i + isize2(nk), mpi%isize )
                    END DO
                  END IF
                END DO
                CALL MPI_WAITALL( rcvreqd, rcvreq, MPI_STATUSES_IGNORE, ierr(1) )
                ! Necessary to avoid compiler optimization
                ! Compiler does not know that hybdat%nobd is modified in mpi_waitall
                CALL MPI_GET_ADDRESS( hybdat%nobd, addr, ierr(1) )
                rcvreqd = 0

#             endif

              ! spread hybdat%nobd from IBZ to whole BZ
              DO nk = 1,kpts%nkptf
                i       = kpts%bkp(nk)
                hybdat%nobd(nk)= hybdat%nobd(i)
              END DO

              !
              ! generate eigenvectors z and MT coefficients from the previous iteration
              ! at all k-points
              !
              CALL gen_wavf(&
     &                 kpts%nkpt,kpts,it,sym,&
     &                 atoms,el_eig,&
     &                 ello_eig,cell,dimension,&
     &                 hybrid,vr0,&
     &                 hybdat%kveclo_eig,&
     &                 noco,oneD,mpi,irank2,&
     &                 hybdat%nbands,input,jsp,&
     &                 zmat,&
     &                 hybdat%bas1,hybdat%bas2,hybdat%bas1_MT,hybdat%drbas1_MT)

              ! generate core wave functions (-> core1/2(jmtd,hybdat%nindxc,0:lmaxc,ntype) )
              CALL corewf(atoms,jsp,input,dimension,vr0,&
     &                    hybdat%lmaxcd,hybdat%maxindxc,mpi,&
     &                    hybdat%lmaxc,hybdat%nindxc,hybdat%core1,hybdat%core2,hybdat%eig_c)

#             ifdef CPP_MPI
                ! wait until all files are written in gen_wavf
                CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#             endif

              !
              ! check olap between core-basis/core-valence/basis-basis
              !
              CALL checkolap(atoms,hybdat%lmaxc,hybdat%lmaxcd,hybdat%nindxc,hybdat%maxindxc,&
     &                      hybdat%core1,hybdat%core2,hybrid,hybdat%bas1,&
     &                      hybdat%bas2,kpts%nkpt,kpts,&
     &                      dimension,mpi,irank2,skip_kpt,&
     &                      input,sym,&
     &                      noco,cell,lapw,jsp,&
     &                       maxval(hybdat%nbands),hybdat%nbands)

              !
              ! set up pointer pntgpt
              !

              ! setup dimension of pntgpt
              hybdat%pntgptd = 0
              DO nk = 1,kpts%nkptf
                CALL apws(dimension,input,noco, kpts,nk,cell,sym%zrfs,&
     &                    1,jsp, bk,lapw,matind,nred)
                hybdat%pntgptd(1) = maxval( (/ ( abs(lapw%k1(i,jsp)),i=1,lapw%nv(jsp)), hybdat%pntgptd(1) /) )
                hybdat%pntgptd(2) = maxval( (/ ( abs(lapw%k2(i,jsp)),i=1,lapw%nv(jsp)), hybdat%pntgptd(2) /) )
                hybdat%pntgptd(3) = maxval( (/ ( abs(lapw%k3(i,jsp)),i=1,lapw%nv(jsp)), hybdat%pntgptd(3) /) )
              END DO

              ALLOCATE( hybdat%pntgpt(-hybdat%pntgptd(1):hybdat%pntgptd(1), -hybdat%pntgptd(2):hybdat%pntgptd(2),&
     &                         -hybdat%pntgptd(3):hybdat%pntgptd(3),kpts%nkptf),stat=ok )
              IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation pntgpt'
              hybdat%pntgpt = 0
              DO nk = 1,kpts%nkptf
                CALL apws( dimension,input,noco, kpts,nk,cell,sym%zrfs,&
     &                     1,jsp, bk,lapw,matind,nred)
                DO i = 1,lapw%nv(jsp)
                  g = (/ lapw%k1(i,jsp),lapw%k2(i,jsp),lapw%k3(i,jsp) /)
                  hybdat%pntgpt(g(1),g(2),g(3),nk) = i
                END DO
              END DO

              ALLOCATE ( basprod(atoms%jmtd),stat=ok )
              IF( ok .ne. 0 )STOP 'eigen_hf: failure allocation basprod'
              ALLOCATE ( hybdat%prodm(hybrid%maxindxm1,hybrid%maxindxp1,0:hybrid%maxlcutm1,atoms%ntype), stat= ok )
              IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%prodm'
              ALLOCATE ( hybdat%prod(hybrid%maxindxp1,0:hybrid%maxlcutm1,atoms%ntype),stat= ok )
              IF( ok .ne. 0 ) STOP 'eigen_hf: failure allocation hybdat%prod'
              basprod = 0 ; hybdat%prodm = 0 ; hybdat%prod%l1 = 0 ; hybdat%prod%l2 = 0
              hybdat%prod%n1 = 0 ; hybdat%prod%n2 = 0

              hybrid%nindxp1 = 0
              DO itype = 1,atoms%ntype
                ng = atoms%jri(itype)
                DO l2 = 0,min(atoms%lmax(itype),hybrid%lcutwf(itype))
                  ll = l2
                  DO l1 = 0,ll
                    IF(abs(l1-l2).le.hybrid%lcutm1(itype)) THEN
                      DO n2 = 1,hybrid%nindx(l2,itype)
                        nn = hybrid%nindx(l1,itype)
                        IF(l1.eq.l2) nn = n2
                        DO n1 = 1,nn
                          ! Calculate all basis-function hybdat%products to obtain
                          ! the overlaps with the hybdat%product-basis functions (hybdat%prodm)
                          basprod(:ng) = ( hybdat%bas1(:ng,n1,l1,itype)*hybdat%bas1(:ng,n2,l2,itype) +hybdat%bas2(:ng,n1,l1,itype)*hybdat%bas2(:ng,n2,l2,itype)) / atoms%rmsh(:ng,itype)
                          DO l = abs(l1-l2),min(hybrid%lcutm1(itype),l1+l2)
                            IF(mod(l1+l2+l,2).eq.0) THEN
                              hybrid%nindxp1(l,itype)    = hybrid%nindxp1(l,itype) + 1
                              n                  = hybrid%nindxp1(l,itype)
                              hybdat%prod(n,l,itype)%l1 = l1
                              hybdat%prod(n,l,itype)%l2 = l2
                              hybdat%prod(n,l,itype)%n1 = n1
                              hybdat%prod(n,l,itype)%n2 = n2
                              DO i = 1,hybrid%nindxm1(l,itype)
                                hybdat%prodm(i,n,l,itype) = intgrf(basprod(:ng)*hybrid%basm1(:ng,i,l,itype), atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)
                              END DO
                            END IF
                          END DO

                        END DO
                      END DO
                    END IF
                  END DO
                END DO
              END DO
              DEALLOCATE(basprod)
              CALL timestop("gen_bz and gen_wavf")
             
         
            ELSE IF (hybrid%l_hybrid ) THEN ! hybrid%l_calhf is false

              ! Reading the eig file
              ALLOCATE ( eig_irr(dimension%neigd2,1) )
              ALLOCATE ( hybdat%kveclo_eig(atoms%nlotot,1) )
              !DO nk = n_start,kpts%nkpt,n_stride
              DO nk = 1,kpts%nkpt,1
                nrec1 = kpts%nkpt*(jsp-1) + nk
                CALL read_eig(eig_id_hf,nk,jsp,el=el_eig, ello=ello_eig,neig=hybdat%ne_eig(nk))
                hybdat%nobd(nk) = count( results%w_iks(:hybdat%ne_eig(nk),nk,jsp) > 0.0 )
              END DO
              DEALLOCATE ( eig_irr , hybdat%kveclo_eig )

              hybrid%maxlmindx = maxval((/ ( sum( (/ (hybrid%nindx(l,itype)*(2*l+1), l=0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
              hybdat%nbands    = min( bands, dimension%neigd )

           ENDIF ! hybrid%l_calhf
         end subroutine eigen_hf_setup
       end module m_eigen_hf_setup
