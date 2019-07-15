      MODULE m_checkolap

      CONTAINS

      SUBROUTINE checkolap(atoms,hybdat,&
     &                  hybrid,&
     &                  nkpti,kpts,&
     &                  dimension,mpi,skip_kpt,&
     &                  input,sym,noco,&
     &                  cell,lapw,jsp)
      USE m_util     , ONLY: intgrf,intgrf_init,chr,sphbessel,harmonicsr
      USE m_constants
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE

      TYPE(t_hybdat),INTENT(IN)   :: hybdat

      TYPE(t_mpi),INTENT(IN)         :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)      :: hybrid
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_noco),INTENT(IN)        :: noco
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_kpts),INTENT(IN)        :: kpts
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_lapw),INTENT(INOUT)     :: lapw

      ! - scalars -
      INTEGER, INTENT(IN)     :: jsp
      INTEGER, INTENT(IN)     ::  nkpti 


      ! - arrays -
      LOGICAL, INTENT(IN)     ::  skip_kpt(nkpti)

      ! - local scalars -
      INTEGER                 ::  i,itype,iatom,ikpt,ineq,igpt,iband
      INTEGER                 ::  irecl_cmt
      INTEGER                 ::  j,m
      INTEGER                 ::  l
      INTEGER                 :: lm,lm1
      INTEGER                 ::  n,nred, nbasfcn

      REAL                    ::  rdum,rdum1
      REAL                    ::  qnorm

      COMPLEX                 ::  cexp,cdum
      COMPLEX , PARAMETER     ::  img = (0.0,1.0)

      ! -local arrays -
      INTEGER                 ::  iarr(2),gpt(3)
      INTEGER , ALLOCATABLE   ::  olapcv_loc(:,:,:,:,:)

      REAL                    ::  sphbes(0:atoms%lmaxd)
      REAL                    ::  q(3)
      REAL                    ::  integrand(atoms%jmtd)
      REAL                    ::  bkpt(3) 
      REAL                    ::  rarr(maxval(hybrid%nbands))
      REAL                    ::  rtaual(3)
      REAL    , ALLOCATABLE   ::  olapcb(:)
      REAL    , ALLOCATABLE   :: olapcv_avg(:,:,:,:),olapcv_max(:,:,:,:)
      TYPE(t_mat),ALLOCATABLE :: z(:)

      COMPLEX                 ::  cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat,nkpti)
      COMPLEX                 ::  y((atoms%lmaxd+1)**2)
      COMPLEX , ALLOCATABLE   ::  olapcv(:,:),olapww(:,:)
      COMPLEX , ALLOCATABLE   ::  carr1(:,:),carr2(:,:),carr3(:,:)

      CHARACTER, PARAMETER    ::  lchar(0:38) =&
     &          (/'s','p','d','f','g','h','i','j','k','l','m','n','o',&
     &            'x','x','x','x','x','x','x','x','x','x','x','x','x',&
     &            'x','x','x','x','x','x','x','x','x','x','x','x','x' /)
      LOGICAL                 ::  l_mism = .true.

      ALLOCATE(z(nkpti))
      DO ikpt=1,nkpti
         CALL lapw%init(input,noco,kpts,atoms,sym,ikpt,cell,sym%zrfs)
         nbasfcn = MERGE(lapw%nv(1)+lapw%nv(2)+2*atoms%nlotot,lapw%nv(1)+atoms%nlotot,noco%l_noco)
         call z(ikpt)%alloc(sym%invs,nbasfcn,dimension%neigd)
      ENDDO
      
      IF ( mpi%irank == 0 ) WRITE(6,'(//A)') '### checkolap ###'

      cmt = 0
    
      ! initialize gridf -> was done in eigen_HF_init
      !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

      ! read in cmt
      DO ikpt = 1,nkpti
#ifdef CPP_MPI
        IF ( skip_kpt(ikpt) ) CYCLE
#endif
        call read_cmt(cmt(:,:,:,ikpt),ikpt)
     END DO
     
      IF ( mpi%irank == 0 ) WRITE(6,'(/A)') ' Overlap <core|core>'
      DO itype=1,atoms%ntype
        IF(atoms%ntype.gt.1 .AND. mpi%irank==0) &
     &     WRITE(6,'(A,I3)') ' Atom type',itype
        DO l=0,hybdat%lmaxc(itype)
          DO i=1,hybdat%nindxc(l,itype)
            IF ( mpi%irank == 0 )&
     &        WRITE(6,'(1x,I1,A,2X)',advance='no') i+l,lchar(l)
            DO j=1,i
              integrand =  hybdat%core1(:,i,l,itype)*hybdat%core1(:,j,l,itype)&
     &                  +  hybdat%core2(:,i,l,itype)*hybdat%core2(:,j,l,itype)
              IF ( mpi%irank == 0 ) WRITE(6,'(F10.6)', advance='no')&
     &           intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)
            END DO
            IF ( mpi%irank == 0 ) WRITE(6,*)
          END DO 
        END DO
      END DO

      IF ( mpi%irank == 0 ) WRITE(6,'(/A)') ' Overlap <core|basis>'
      ALLOCATE( olapcb(hybrid%maxindx),olapcv(maxval(hybrid%nbands),nkpti),&
     &          olapcv_avg(  -hybdat%lmaxcd:hybdat%lmaxcd,hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype),&
     &          olapcv_max(  -hybdat%lmaxcd:hybdat%lmaxcd,hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype),&
     &          olapcv_loc(2,-hybdat%lmaxcd:hybdat%lmaxcd,hybdat%maxindxc,0:hybdat%lmaxcd,atoms%ntype) )

      DO itype=1,atoms%ntype
        IF(atoms%ntype.gt.1 .AND. mpi%irank==0) &
     &     WRITE(6,'(A,I3)') ' Atom type',itype
        DO l=0,hybdat%lmaxc(itype)
          IF(l.gt.atoms%lmax(itype)) EXIT ! very improbable case
          IF ( mpi%irank == 0 ) &
     &        WRITE(6,"(9X,'u(',A,')',4X,'udot(',A,')',:,3X,'ulo(',A,"//&
     &                "') ...')") (lchar(l),i=1,min(3,hybrid%nindx(l,itype)))
          DO i=1,hybdat%nindxc(l,itype)
            IF ( mpi%irank == 0 )&
     &        WRITE(6,'(1x,I1,A,2X)',advance='no') i+l,lchar(l)
            DO j=1,hybrid%nindx(l,itype)

              integrand = hybdat%core1(:,i,l,itype)*hybdat%bas1(:,j,l,itype)&
     &                  + hybdat%core2(:,i,l,itype)*hybdat%bas2(:,j,l,itype)

              olapcb(j) = &
     &              intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)

              IF ( mpi%irank == 0 )&
     &          WRITE(6,'(F10.6)',advance='no') olapcb(j)
            END DO

            lm    = sum ( (/ (hybrid%nindx(j,itype)*(2*j+1),j=0,l-1) /) )
            iatom = sum(atoms%neq(1:itype-1))+1 ! take first of group of equivalent atoms
            DO m=-l,l
              olapcv = 0
              DO j=1,hybrid%nindx(l,itype)
                lm = lm + 1
                olapcv(:,:) = olapcv(:,:) + &
     &                        olapcb(j)*cmt(:maxval(hybrid%nbands),lm,iatom,:nkpti)
              END DO
              rdum                      = sum    ( abs(olapcv(:,:))**2 )
              rdum1                     = maxval ( abs(olapcv(:,:))    )
              iarr                      = maxloc ( abs(olapcv(:,:))    )
              olapcv_avg(  m,i,l,itype) = &
     &                sqrt( rdum / nkpti / sum(hybrid%nbands(:nkpti)) * nkpti )
              olapcv_max(  m,i,l,itype) = rdum1
              olapcv_loc(:,m,i,l,itype) = iarr
            END DO
            IF ( mpi%irank == 0 ) WRITE(6,*)

          END DO
        END DO
      END DO

      IF( mpi%irank == 0 ) THEN
        WRITE(6,'(/A)') ' Average overlap <core|val>'
        DO itype=1,atoms%ntype
          IF(atoms%ntype.gt.1) write(6,'(A,I3)') ' Atom type',itype
          DO l=0,hybdat%lmaxc(itype)
            DO i=1,hybdat%nindxc(l,itype)
              WRITE(6,'(1x,I1,A,2X)',advance='no') i+l,lchar(l)
              WRITE(6,'('//chr(2*l+1)//'F10.6)') &
     &                                        olapcv_avg(-l:l,i,l,itype)
            END DO
          END DO
        END DO

        WRITE(6,'(/A)') ' Maximum overlap <core|val> at (band/kpoint)'
        DO itype=1,atoms%ntype
          IF(atoms%ntype.gt.1) write(6,'(A,I3)') ' Atom type',itype
          DO l=0,hybdat%lmaxc(itype)
            DO i=1,hybdat%nindxc(l,itype)
              WRITE(6,'(1x,I1,A,2X)',advance='no') i+l,lchar(l)
              WRITE(6,'('//chr(2*l+1)//&
     &                 '(F10.6,'' ('',I3.3,''/'',I4.3,'')''))')&
     &                          (olapcv_max(  m,i,l,itype),&
     &                           olapcv_loc(:,m,i,l,itype),m=-l,l)
            END DO
          END DO
        END DO
      END IF ! mpi%irank == 0

      DEALLOCATE (olapcb,olapcv,olapcv_avg,olapcv_max,olapcv_loc)

      IF ( mpi%irank == 0 ) WRITE(6,'(/A)') ' Overlap <basis|basis>'

      DO itype=1,atoms%ntype
        IF(atoms%ntype.gt.1 .AND. mpi%irank==0) &
     &     WRITE(6,'(A,I3)') ' Atom type',itype
        DO l=0,atoms%lmax(itype)
          DO i=1,hybrid%nindx(l,itype)
            IF ( mpi%irank == 0 ) THEN
              SELECT CASE(i)
               CASE(1)
                 WRITE(6,'(1x,''   u('',A,'')'')',advance='no') lchar(l)
               CASE(2)
                 WRITE(6,'(1x,''udot('',A,'')'')',advance='no') lchar(l)
               CASE DEFAULT
                 WRITE(6,'(1x,'' ulo('',A,'')'')',advance='no') lchar(l)
              END SELECT
            END IF
            DO j=1,i
              integrand = hybdat%bas1(:,i,l,itype)*hybdat%bas1(:,j,l,itype)&
     &                  + hybdat%bas2(:,i,l,itype)*hybdat%bas2(:,j,l,itype)

              IF ( mpi%irank == 0 ) WRITE(6,'(F10.6)',advance='no')&
     &              intgrf(integrand,atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,hybdat%gridf)
            END DO
            IF ( mpi%irank == 0 ) WRITE(6,*)
          END DO
        END DO
      END DO

      IF( .not. l_mism ) RETURN

      IF ( mpi%irank == 0 ) WRITE(6,'(/A)') &
     &          'Mismatch of wave functions at the MT-sphere boundaries'
      ALLOCATE (carr1(maxval(hybrid%nbands),(atoms%lmaxd+1)**2))
      ALLOCATE (carr2(maxval(hybrid%nbands),(atoms%lmaxd+1)**2))
      ALLOCATE (carr3(maxval(hybrid%nbands),(atoms%lmaxd+1)**2))
      DO ikpt = 1,nkpti
         call read_z(z(ikpt),ikpt)
      END DO

      iatom = 0
      DO itype = 1,atoms%ntype
        DO ineq = 1,atoms%neq(itype)
          iatom = iatom + 1
          IF ( mpi%irank == 0 ) THEN
            if(atoms%nat.gt.1) WRITE(6,'(2X,A,I3)') 'Atom',iatom
            WRITE(6,'(2X,A)') 'k-point    average      (   maximum    )'
          END IF

          DO ikpt = 1,nkpti
#ifdef CPP_MPI
            IF ( skip_kpt(ikpt) ) CYCLE
#endif
            carr1 = 0; carr2 = 0; carr3 = 0

            ! calculate k1,k2,k3
            CALL lapw%init(input,noco,kpts,atoms,sym,ikpt,cell,sym%zrfs)

            ! PW part
            DO igpt = 1,lapw%nv(jsp)
              gpt(1) = lapw%k1(igpt,jsp)
              gpt(2) = lapw%k2(igpt,jsp)
              gpt(3) = lapw%k3(igpt,jsp)

              cexp = exp(img * 2*pi_const * &
     &                   dot_product(kpts%bkf(:,ikpt)+gpt,atoms%taual(:,iatom)) )
              q    = matmul(kpts%bkf(:,ikpt)+gpt,cell%bmat)

              qnorm = sqrt(sum(q**2))
              call sphbessel(sphbes,atoms%rmt(itype)*qnorm,atoms%lmax(itype))
              call harmonicsr(y,q,atoms%lmax(itype))
              y  = conjg(y)
              lm = 0
              DO l = 0,atoms%lmax(itype)
                cdum = 4*pi_const*img**l/sqrt(cell%omtil) * sphbes(l) * cexp
                DO m = -l,l
                  lm = lm + 1
                  DO iband = 1,hybrid%nbands(ikpt)
                     if (z(1)%l_real) THEN
                        carr2(iband,lm) = carr2(iband,lm) + cdum * z(ikpt)%data_r(igpt,iband) * y(lm)
                     Else
                        carr2(iband,lm) = carr2(iband,lm) + cdum * z(ikpt)%data_c(igpt,iband) * y(lm)
                     END if
                  end DO
                END DO
              END DO
            END DO

            ! MT
            lm  = 0
            lm1 = 0
            DO l = 0,atoms%lmax(itype)
              DO m = -l,l
                lm = lm + 1
                DO n = 1,hybrid%nindx(l,itype)
                  lm1  = lm1 + 1
                  rdum = hybdat%bas1(atoms%jri(itype),n,l,itype) / atoms%rmt(itype)
                  DO iband = 1,hybrid%nbands(ikpt)
                    carr3(iband,lm) = carr3(iband,lm) + cmt(iband,lm1,iatom,ikpt) * rdum
                  END DO
                END DO
              END DO
            END DO
            carr1 = carr2 - carr3

            rarr = 0
            lm   = 0
            DO l = 0,atoms%lmax(itype)
              DO m = -l,l
                lm   = lm + 1
                rarr = rarr + abs(carr1(:,lm))**2
              END DO
            END DO
            rarr = sqrt ( rarr / (4*pi_const) )
 !             WRITE(outtext,'(I6,4X,F14.12,''  ('',F14.12,'')'')') &
 !    &              ikpt,sum(rarr(:1)**2/nbands(ikpt)),maxval(rarr(:1))
 !             CALL writeout(outtext,mpi%irank)
!             IF( iatom .eq. 6 ) THEN
!               cdum = exp(2*pi*img*dot_product(bkf(:,ikpt),(/0.0,0.0,1.0/) ))
!               lm = 0 
!               DO l = 0,lmax(itype)
!                 DO m = -l,l
!                   lm = lm + 1
!                   DO iband = 1,nbands(ikpt)
!                     WRITE(700+ikpt,'(3i4,6f15.10)') iband,l,m,carr2(iband,lm),carr3(iband,lm),
!      &                                              carr2(iband,lm)/(carr3(iband,lm))
!                   END DO
!                 END DO
!               END DO
!             END IF

          END DO
        END DO
      END DO


      END SUBROUTINE checkolap

      END MODULE m_checkolap
