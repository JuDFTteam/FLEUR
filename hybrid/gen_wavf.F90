!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This module generates the cmt coefficients and eigenvectors z   !
!     at all kpoints nkpt from the irreducible kpoints nkpti          !
!     and writes them out in cmt and z, respectively.                 !
!                                                 M.Betzinger(09/07)  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
MODULE m_gen_wavf

CONTAINS

   SUBROUTINE gen_wavf (nkpti,kpts,it,sym,atoms,el_eig,ello_eig,cell,dimension,hybrid,vr0,&
                        hybdat,noco,oneD,mpi,irank2,input,jsp,zmat)

      ! nkpti      ::     number of irreducible k-points
      ! nkpt       ::     number of all k-points 

      USE m_radfun
      USE m_radflo
      USE m_abcof
      USE m_trafo     ,ONLY: waveftrafo_genwavf
      USE m_util      ,ONLY: modulo1
      USE m_olap
      USE m_types
      USE m_hyb_abcrot
      USE m_io_hybrid

      IMPLICIT NONE

      TYPE(t_hybdat),    INTENT(INOUT) :: hybdat
      TYPE(t_mpi),       INTENT(IN)    :: mpi
      TYPE(t_dimension), INTENT(IN)    :: dimension
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_hybrid),    INTENT(IN)    :: hybrid
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_kpts),      INTENT(IN)    :: kpts
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_mat),       INTENT(IN)    :: zmat(:) !for all kpoints 

      INTEGER,           INTENT(IN)    :: nkpti, it
      INTEGER,           INTENT(IN)    :: jsp

      INTEGER,           INTENT(IN)    :: irank2(nkpti)
      REAL,              INTENT(IN)    :: vr0(:,:,:)!(jmtd,ntype,jspd)
      REAL,              INTENT(IN)    :: el_eig(0:atoms%lmaxd,atoms%ntype)
      REAL,              INTENT(IN)    :: ello_eig(atoms%nlod,atoms%ntype)

      ! local scalars
      INTEGER                 :: ilo,idum,m,irecl_cmt,irecl_z
      COMPLEX                 :: cdum
      TYPE(t_mat)             :: zhlp
      INTEGER                 :: ikpt0,ikpt,itype,iop,ispin,ieq,indx,iatom
      INTEGER                 :: i,j,l ,ll,lm,ng,ok
      COMPLEX                 :: img=(0d0,1d0)

      INTEGER                 :: nodem,noded
      REAL                    :: wronk

      INTEGER                 :: lower, upper
      LOGICAL                 :: found

      ! local arrays
      INTEGER                 :: rrot(3,3,sym%nsym)
      INTEGER                 :: map_lo(atoms%nlod)
      INTEGER                 :: iarr(0:atoms%lmaxd,atoms%ntype)
      COMPLEX,ALLOCATABLE     :: acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)
      
      COMPLEX,ALLOCATABLE     :: cmt(:,:,:),cmthlp(:,:,:)

      REAL                    :: vr(atoms%jmtd,atoms%ntype,input%jspins)
      REAL,ALLOCATABLE        :: f(:,:,:),df(:,:,:)

      REAL                    :: flo(atoms%jmtd,2,atoms%nlod)
      REAL                    :: uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
      REAL                    :: ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
   
      REAL                    :: bkpt(3)

!     local arrays for abcof1
!      COMPLEX                 ::  a(nvd,0:lmd,natd,nkpti),b(nvd,0:lmd,natd,nkpti)


      TYPE(t_lapw)  :: lapw(kpts%nkptf)
      TYPE(t_usdus) :: usdus

      CALL usdus%init(atoms,input%jspins)
      CALL zhlp%alloc(zmat(1)%l_real,zmat(1)%matsize1,zmat(1)%matsize2)

      
      ! setup rotations in reciprocal space
      DO iop = 1, sym%nsym
         IF(iop.LE.sym%nop) THEN
            rrot(:,:,iop) = transpose(sym%mrot(:,:,sym%invtab(iop)))
         ELSE
            rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
         END IF
      END DO

      ! generate G-vectors, which fulfill |k+G|<rkmax
      ! for all k-points
      DO ikpt = 1, kpts%nkptf
         CALL lapw(ikpt)%init(input,noco,kpts,atoms,sym,ikpt,cell,sym%zrfs)
      END DO

      ! set spherical component of the potential from the previous iteration vr
      vr = vr0


!       ALLOCATE ( z_out(nbasfcn,neigd,nkpti),stat=ok )
!       IF ( ok .ne. 0) STOP 'gen_wavf: failure allocation z'
!       z_out = 0
!       z_out(:,:,:nkpti) = z_in


      ! calculate radial basis functions belonging to the
      ! potential vr stored in bas1 and bas2
      ! bas1 denotes the large component
      ! bas2    "     "  small component

      ALLOCATE(f(atoms%jmtd,2,0:atoms%lmaxd),df(atoms%jmtd,2,0:atoms%lmaxd))
      f = 0
      df = 0
      iarr = 2
      DO itype = 1, atoms%ntype
         IF (mpi%irank == 0) WRITE (6,FMT=8000) itype
         ng = atoms%jri(itype)
         DO l=0,atoms%lmax(itype)
            CALL radfun(l,itype,1,el_eig(l,itype),vr(:,itype,jsp),atoms,f(:,:,l),df(:,:,l),usdus,nodem,noded,wronk)
            IF (mpi%irank == 0 ) WRITE (6,FMT=8010) l,el_eig(l,itype),usdus%us(l,itype,1),usdus%dus(l,itype,1),nodem,&
                                                    usdus%uds(l,itype,1),usdus%duds(l,itype,1),noded,usdus%ddn(l,itype,1),wronk

            hybdat%bas1(1:ng,1,l,itype) =  f(1:ng,1,l)
            hybdat%bas2(1:ng,1,l,itype) =  f(1:ng,2,l)
            hybdat%bas1(1:ng,2,l,itype) = df(1:ng,1,l)
            hybdat%bas2(1:ng,2,l,itype) = df(1:ng,2,l)

            hybdat%bas1_MT(1,l,itype)   = usdus%us(l,itype,1)
            hybdat%drbas1_MT(1,l,itype) = usdus%dus(l,itype,1)
            hybdat%bas1_MT(2,l,itype)   = usdus%uds(l,itype,1)
            hybdat%drbas1_MT(2,l,itype) = usdus%duds(l,itype,1)
         END DO

         IF (atoms%nlo(itype).GE.1) THEN
            CALL radflo(atoms,itype,jsp,ello_eig,vr(:,itype,jsp),f,df,mpi,usdus,uuilon,duilon,ulouilopn,flo)

            DO ilo=1,atoms%nlo(itype)
               iarr(atoms%llo(ilo,itype),itype) = iarr(atoms%llo(ilo,itype),itype) + 1
               hybdat%bas1(1:ng,iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype) = flo(1:ng,1,ilo)
               hybdat%bas2(1:ng,iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype) = flo(1:ng,2,ilo)
               hybdat%bas1_MT(iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype) = usdus%ulos(ilo,itype,1)
               hybdat%drbas1_MT(iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype) = usdus%dulos(ilo,itype,1)
            END DO
         END IF
      END DO
      DEALLOCATE (f,df)

#if CPP_DEBUG
      ! consistency check
      IF(.not.all(iarr.eq.hybrid%nindx)) STOP 'gen_wavf: counting error'
#endif

      8000 FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',/,t32,'radial function',t79,&
                   'energy derivative',/,t3,'l',t8,'energy',t26,'value',t39,'derivative',t53,&
                   'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,'norm',t119,'wronskian')
      8010 FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

      ! determine boundaries for parallel calculations
      lower = 1
      upper = nkpti
      found = .false.
#ifdef CPP_MPI
      DO ikpt = 1, nkpti
         IF (irank2(ikpt) == 0 .AND. .NOT.found) THEN
            lower = ikpt
            found = .true.
         ELSE IF (irank2(ikpt) /= 0 .AND. found) THEN
            upper = ikpt-1
            EXIT
         END IF
      END DO
#else
      found = .true.
#endif
      IF (.NOT.found) THEN
         upper = 0
      END IF

      ! calculate wavefunction expansion in the the MT region
      ! (acof,bcof,ccof) and APW-basis coefficients
      ! (a,b,bascofold_lo) at irred. kpoints

      ALLOCATE(acof(dimension%neigd,0:dimension%lmd,atoms%nat),stat=ok)
      IF (ok.NE.0) STOP 'gen_wavf: failure allocation acof'
      ALLOCATE(bcof(dimension%neigd,0:dimension%lmd,atoms%nat),stat=ok)
      IF (ok.NE.0) STOP 'gen_wavf: failure allocation bcof'
      ALLOCATE(ccof(-atoms%llod:atoms%llod,dimension%neigd,atoms%nlod,atoms%nat),stat=ok)
      IF (ok.NE.0) STOP 'gen_wavf: failure allocation ccof'
      ALLOCATE (cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat),stat=ok)
      IF (ok.NE.0) STOP 'gen_wavf: Failure allocation cmt'
      ALLOCATE (cmthlp(dimension%neigd,hybrid%maxlmindx,atoms%nat),stat=ok)
      IF (ok.NE.0) STOP 'gen_wavf: failure allocation cmthlp'

      DO ikpt0 = lower, upper

         acof = 0; bcof = 0; ccof = 0

         ! abcof calculates the wavefunction coefficients
         ! stored in acof,bcof,ccof
         lapw(ikpt0)%nmat = lapw(ikpt0)%nv(jsp) + atoms%nlotot
         CALL abcof(input,atoms,sym,cell,lapw(ikpt0),hybrid%nbands(ikpt0),usdus,noco,jsp,&!hybdat%kveclo_eig(:,ikpt0),&
                   oneD,acof(: hybrid%nbands(ikpt0),:,:),bcof(: hybrid%nbands(ikpt0),:,:),&
                   ccof(:,: hybrid%nbands(ikpt0),:,:),zmat(ikpt0))

! call was ...
          ! gpt(1,:,:,ikpt0),gpt(2,:,:,ikpt0),&
          ! gpt(3,:,:,ikpt0),ngpt(:,ikpt0),&!k1hlp,k2hlp,k3hlp,nvhlp,&
          !    ngpt(jsp,ikpt0)+nbands(ikpt0),z(:,:,ikpt0),&!nvhlp(jsp)+ &
          !   &usdus,&
          !    noco,&
          !    jsp,kveclo_eig(:ikpt0),oneD,oneD,&
          !    acof(:nbands(ikpt0),:,:),&
          !    bcof(:nbands(ikpt0),:,:),ccof(:,:nbands(ikpt0),:,:) )

         ! MT wavefunction coefficients are calculated in a local coordinate system rotate them in the global one

         CALL hyb_abcrot(hybrid,atoms,hybrid%nbands(ikpt0),sym,cell,oneD,acof(: hybrid%nbands(ikpt0),:,:),&
                         bcof(: hybrid%nbands(ikpt0),:,:),ccof(:,: hybrid%nbands(ikpt0),:,:))

         ! decorate acof, bcof, ccof with coefficient i**l and store them
         ! in the field cmt(neigd,nkpt,maxlmindx,nat), i.e.
         ! where maxlmindx subsumes l,m and nindx

         cmt = 0
         iatom = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx = 0
               DO l = 0, atoms%lmax(itype)
                  ll = l*(l+1)
                  cdum = img**l

                  ! determine number of local orbitals with quantum number l
                  ! map returns the number of the local orbital of quantum
                  ! number l in the list of all local orbitals of the atom type
                  idum   = 0
                  map_lo = 0
                  IF (hybrid%nindx(l,itype).GT.2) THEN
                     DO j = 1, atoms%nlo(itype)
                        IF (atoms%llo(j,itype).EQ.l) THEN
                           idum = idum + 1
                           map_lo(idum) = j
                        END IF
                     END DO
                  END IF

                  DO M = -l, l
                     lm = ll + M
                     DO i = 1, hybrid%nindx(l,itype)
                        indx = indx + 1
                        IF(i.EQ.1) THEN
                           cmt(:,indx,iatom) = cdum * acof(:,lm,iatom)
                        ELSE IF (i.EQ.2) THEN
                           cmt(:,indx,iatom) = cdum * bcof(:,lm,iatom)
                        ELSE
                           idum = i - 2
                           cmt(:,indx,iatom) = cdum * ccof(M,:,map_lo(idum),iatom)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! write cmt at irreducible k-points in direct-access file cmt
         CALL write_cmt(cmt,ikpt0)
         CALL zhlp%alloc(zmat(ikpt0)%l_real,zmat(ikpt0)%matsize1,zmat(ikpt0)%matsize2)
        
         IF (zhlp%l_real) THEN
            zhlp%data_r = zmat(ikpt0)%data_r
         ELSE
            zhlp%data_c = zmat(ikpt0)%data_c
         END IF
         CALL write_z(zhlp,ikpt0)

         ! generate wavefunctions coefficients at all k-points from
         ! irreducible k-points

         DO ikpt = 1, kpts%nkptf
            IF ((kpts%bkp(ikpt).EQ.ikpt0).AND.(ikpt0.NE.ikpt)) THEN
               iop = kpts%bksym(ikpt)
               CALL waveftrafo_genwavf(cmthlp,zhlp%data_r,zhlp%data_c,cmt(:,:,:),zmat(1)%l_real,zmat(ikpt0)%data_r(:,:),&
                                       zmat(ikpt0)%data_c(:,:),ikpt0,iop,atoms,hybrid,kpts,sym,jsp,dimension,&
                                       hybrid%nbands(ikpt0),cell,lapw(ikpt0),lapw(ikpt),.true.)

               CALL write_cmt(cmthlp,ikpt)
               CALL write_z(zhlp,ikpt)
            END IF
         END DO  !ikpt
      END DO !ikpt0

      DEALLOCATE(acof,bcof,ccof)
      DEALLOCATE(cmt,cmthlp)

   END SUBROUTINE gen_wavf

END MODULE m_gen_wavf
