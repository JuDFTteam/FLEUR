MODULE m_stden

USE m_juDFT
!     ************************************************************
!     generate flapw starting density by superposition of
!     atomic densities. the non-spherical terms inside
!     the spheres are obtained by a least squares fit
!     and the interstitial and vacuum warping terms are calculated
!     by a fast fourier transform.
!     e. wimmer   nov. 1984       c.l.fu diagonized 1987
!     *************************************************************

CONTAINS

SUBROUTINE stden(mpi,sphhar,stars,atoms,sym,vacuum,&
                 input,cell,xcpot,noco,oneD)

   USE m_constants
   USE m_qsf
   USE m_checkdopall
   USE m_cdnovlp
   USE m_cdn_io
   USE m_qfix
   USE m_atom2
   USE m_types
   USE m_types_xcpot_inbuild
   USE m_juDFT_init

   IMPLICIT NONE

   TYPE(t_mpi),INTENT(IN)      :: mpi
   TYPE(t_atoms),INTENT(IN)    :: atoms
   
   TYPE(t_sphhar),INTENT(IN)   :: sphhar
   TYPE(t_sym),INTENT(IN)      :: sym
   TYPE(t_stars),INTENT(IN)    :: stars
   TYPE(t_noco),INTENT(IN)     :: noco
   TYPE(t_oneD),INTENT(IN)     :: oneD
   TYPE(t_input),INTENT(IN)    :: input
   TYPE(t_vacuum),INTENT(IN)   :: vacuum
   TYPE(t_cell),INTENT(IN)     :: cell
   CLASS(t_xcpot),INTENT(IN)   :: xcpot

   ! Local type instances
   TYPE(t_potden)   :: den
   TYPE(t_enpara)   :: enpara
   ! Local Scalars
   REAL d,del,fix,h,r,rnot,z,bm,qdel,va
   REAL denz1(1,1),vacxpot(1,1),vacpot(1,1)
   INTEGER i,ivac,iza,j,jr,k,n,n1,ispin
   INTEGER nw,ilo,natot,nat

   ! Local Arrays
   REAL,    ALLOCATABLE :: vbar(:,:)
   REAL,    ALLOCATABLE :: rat(:,:),eig(:,:,:),sigm(:)
   REAL,    ALLOCATABLE :: rh(:,:,:),rh1(:,:,:),rhoss(:,:)
   REAL     :: vacpar(2)
   INTEGER lnum(29,atoms%ntype),nst(atoms%ntype)
   INTEGER jrc(atoms%ntype)
   LOGICAL l_found(0:3),llo_found(atoms%nlod),l_st
   REAL,ALLOCATABLE   :: occ(:,:)
   ! Data statements
   DATA del/1.e-6/
   PARAMETER (l_st=.true.)

   IF (input%jspins > input%jspins) CALL juDFT_error("input%jspins > input%jspins", calledby = "stden")

   CALL den%init(stars,atoms,sphhar,vacuum,noco,input%jspins,POTDEN_TYPE_DEN)

   ALLOCATE ( rat(atoms%msh,atoms%ntype),eig(29,input%jspins,atoms%ntype) )
   ALLOCATE ( rh(atoms%msh,atoms%ntype,input%jspins),rh1(atoms%msh,atoms%ntype,input%jspins) )
   ALLOCATE ( vbar(2,atoms%ntype),sigm(vacuum%nmzd) )
   ALLOCATE ( rhoss(atoms%msh,input%jspins) )

   rh = 0.0
   rhoss = 0.0

   IF (mpi%irank == 0) THEN
      ! if sigma is not 0.0, then divide this charge among all atoms
      !TODO: reactivate efields
      !IF ( ABS(input%sigma).LT. 1.e-6) THEN
      IF (1==1) THEN
         qdel = 0.0
      ELSE
         natot = 0
         DO n = 1, atoms%ntype
            IF (atoms%zatom(n).GE.1.0) natot = natot + atoms%neq(n)
         END DO
         !qdel = 2.*input%sigma/natot
      END IF

      WRITE (6,FMT=8000)
      8000 FORMAT (/,/,/,' superposition of atomic densities',/,/,' original atomic densities:',/)
      DO n = 1,atoms%ntype
         r = atoms%rmsh(1,n)
         d = EXP(atoms%dx(n))
         jrc(n) = 0
         DO WHILE (r < atoms%rmt(n) + 20.0)
            IF (jrc(n) > atoms%msh) CALL juDFT_error("increase msh in fl7para!",calledby ="stden")
            jrc(n) = jrc(n) + 1
            rat(jrc(n),n) = r
            r = r*d
         END DO
      END DO

      ! Generate the atomic charge densities
      DO n = 1,atoms%ntype
         z = atoms%zatom(n)
         r = atoms%rmt(n)
         h = atoms%dx(n)
         jr = atoms%jri(n)
         IF (input%jspins.EQ.2) THEN
            bm = atoms%bmu(n)
         ELSE
            bm = 0.
         END IF
         occ=atoms%econf(n)%Occupation(:,:)
         ! check whether this atom has been done already
         DO n1 = 1, n - 1
            IF (ABS(z-atoms%zatom(n1)).GT.del) CYCLE
            IF (ABS(r-atoms%rmt(n1)).GT.del) CYCLE
            IF (ABS(h-atoms%dx(n1)).GT.del) CYCLE
            IF (ABS(bm-atoms%bmu(n1)).GT.del) CYCLE
            IF (ANY(ABS(occ(:,:)-atoms%econf(n1)%Occupation(:,:))>del)) CYCLE
            IF (jr.NE.atoms%jri(n1)) CYCLE
            DO ispin = 1, input%jspins
               DO i = 1,jrc(n) ! atoms%msh
                  rh(i,n,ispin) = rh(i,n1,ispin)
               END DO
            END DO
            nst(n) = nst(n1)
            vbar(1,n) = vbar(1,n1)
            vbar(input%jspins,n) = vbar(input%jspins,n1)
            DO i = 1, nst(n1)
               lnum(i,n)  = lnum(i,n1)
               eig(i,1,n) = eig(i,1,n1)
               eig(i,input%jspins,n) = eig(i,input%jspins,n1)
            END DO
            GO TO 70
         END DO
         !--->    new atom
         rnot = atoms%rmsh(1,n)
         IF (z.LT.1.0) THEN
            va = max(z,1.e-8)/(input%jspins*sfp_const*atoms%volmts(n))
            DO ispin = 1, input%jspins
               DO i = 1,jrc(n) ! atoms%msh
                  rh(i,n,ispin) = va/rat(i,n)**2
               END DO
            END DO
         ELSE
            CALL atom2(atoms,xcpot,input,n,jrc(n),rnot,qdel,&
                       rhoss,nst(n),lnum(1,n),eig(1,1,n),vbar(1,n),.true.)
            DO ispin = 1, input%jspins
               DO i = 1, jrc(n) ! atoms%msh
                  rh(i,n,ispin) = rhoss(i,ispin)
               END DO
            END DO
         END IF
         ! list atomic density
         iza = atoms%zatom(n) + 0.0001
         WRITE (6,FMT=8030) namat_const(iza)
         8030 FORMAT (/,/,' atom: ',a2,/)
         8040 FORMAT (4 (3x,i5,f8.5,f12.6))
         70 CONTINUE
      END DO

      !roa+
      ! use cdnovlp to generate total density out of atom densities
      DO ispin = 1, input%jspins
         nat = 1
         DO n = 1,atoms%ntype
            DO i = 1, jrc(n)
               rh1(i,n,ispin) = rh(i,n,ispin)*fpi_const*rat(i,n)**2
            END DO
            rh1(jrc(n):atoms%msh,n,ispin) = 0.0
            ! prepare spherical mt charge
            DO i = 1,atoms%jri(n)
               den%mt(i,0,n,ispin) = rh(i,n,ispin)*sfp_const*atoms%rmsh(i,n)**2
            END DO
            ! reset nonspherical mt charge
            DO k = 1,sphhar%nlh(sym%ntypsy(nat))
               DO j = 1,atoms%jri(n)
                  den%mt(j,k,n,ispin) = 0.e0
               END DO
            END DO
            nat = nat + atoms%neq(n)
         END DO
      END DO ! ispin
   END IF ! mpi%irank == 0

   DO ispin = 1, input%jspins
      CALL cdnovlp(mpi,sphhar,stars,atoms,sym,vacuum,&
                   cell,input,oneD,l_st,ispin,rh1(:,:,ispin),&
                   den%pw,den%vacxy,den%mt,den%vacz)
      !roa-
   END DO

   IF (mpi%irank == 0) THEN

      ! Check the normalization of total density
      CALL qfix(mpi,stars,atoms,sym,vacuum,sphhar,input,cell,oneD,den,.FALSE.,.FALSE.,.true.,fix)
      z=SUM(atoms%neq(:)*atoms%zatom(:))
      IF (ABS(fix*z-z)>0.5) THEN
         CALL judft_warn("Starting density not charge neutral",hint= &
                         "Your electronic configuration might be broken",calledby="stden.f90")
      END IF

      ! Write superposed density onto density file
      den%iter = 0
      CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,CDN_ARCHIVE_TYPE_CDN1_const,CDN_INPUT_DEN_const,&
                        1,-1.0,0.0,.TRUE.,den)

      ! Check continuity
      IF (input%vchk) THEN
         DO ispin = 1, input%jspins
            WRITE (6,'(a8,i2)') 'spin No.',ispin
            CALL checkDOPAll(input,sphhar,stars,atoms,sym,vacuum,oneD,&
                           cell,den,ispin)
         END DO ! ispin = 1, input%jspins
      END IF ! input%vchk

      ! set up parameters for enpara-file
      IF ((juDFT_was_argument("-genEnpara"))) THEN
         CALL enpara%init_enpara(atoms,input%jspins,input%film)

         enpara%lchange = .TRUE.
         enpara%llochg = .TRUE.

         DO ispin = 1, input%jspins
            ! vacpar is taken as highest occupied level from atomic eigenvalues
            ! vacpar (+0.3)  serves as fermi level for film or bulk
            vacpar(1) = -999.9
            DO n = 1,atoms%ntype
               vacpar(1) = MAX( vacpar(1),eig(nst(n),ispin,n) )
            END DO
            IF (.NOT.input%film) vacpar(1) = vacpar(1) + 0.4
            vacpar(2) = vacpar(1)
            DO n = 1,atoms%ntype
               enpara%skiplo(n,ispin) = 0
               DO i = 0, 3
                  l_found(i) = .FALSE.
                  enpara%el0(i,n,ispin) = vacpar(1)
               END DO
               DO i = 1, atoms%nlo(n)
                  llo_found(i) = .FALSE.
                  enpara%ello0(i,n,ispin) = vacpar(1) - 1.5
               END DO

               ! take energy-parameters from atomic calculation
               DO i = nst(n), 1, -1
                  IF (.NOT.input%film) eig(i,ispin,n) = eig(i,ispin,n) + 0.4
                  IF (.NOT.l_found(lnum(i,n)).AND.(lnum(i,n).LE.3)) THEN
                     enpara%el0(lnum(i,n),n,ispin) = eig(i,ispin,n)
                     IF (enpara%el0(lnum(i,n),n,ispin).LT.input%ellow) THEN
                        enpara%el0(lnum(i,n),n,ispin) = vacpar(1)
                        l_found(lnum(i,n))  = .TRUE.
                     END IF
                     IF (enpara%el0(lnum(i,n),n,ispin).LT.input%elup) THEN
                        l_found(lnum(i,n))  = .TRUE.
                     END IF
                  ELSE
                     IF (l_found(lnum(i,n)).AND.(atoms%nlo(n).GT.0)) THEN
                        DO ilo = 1, atoms%nlo(n)
                           IF (atoms%llo(ilo,n).EQ.lnum(i,n)) THEN
                              IF (.NOT.llo_found(ilo)) THEN
                                 enpara%ello0(ilo,n,ispin) = eig(i,ispin,n)
                                 IF ((enpara%ello0(ilo,n,ispin).GT.input%elup).OR.&
                                     (enpara%ello0(ilo,n,ispin).LT.input%ellow)) THEN
                                    enpara%ello0(ilo,n,ispin)= vacpar(1)
                                 ELSE IF (atoms%l_dulo(ilo,n)) THEN
                                    enpara%ello0(ilo,n,ispin)= enpara%el0(atoms%llo(ilo,n),n,ispin)
                                 ELSE
                                    enpara%skiplo(n,ispin) = enpara%skiplo(n,ispin) + 2*atoms%llo(ilo,n)+1
                                 END IF
                                 llo_found(ilo) = .TRUE.
                                 IF ((enpara%el0(atoms%llo(ilo,n),n,ispin)-enpara%ello0(ilo,n,ispin).LT.0.5)&
                                     .AND.(atoms%llo(ilo,n).GE.0)) THEN
                                    enpara%el0(atoms%llo(ilo,n),n,ispin) = vacpar(1)
                                    IF (atoms%l_dulo(ilo,n)) enpara%ello0(ilo,n,ispin)= enpara%el0(atoms%llo(ilo,n),n,ispin)
                                 END IF
                              END IF
                           END IF
                        END DO ! ilo = 1, atoms%nlo(n)
                     END IF
                  END IF ! .NOT.l_found(lnum(i,n)).AND.(lnum(i,n).LE.3)
               END DO ! i = nst(n), 1, -1
            END DO ! atom types

            IF (input%film) THEN
               ! get guess for vacuum parameters
               ! YM : in 1D case should be modified, but also not that bad like this
               !
               ! generate coulomb potential by integrating inward to z1

               DO ivac = 1, vacuum%nvac
                  DO i=1,vacuum%nmz
                     sigm(i) = (i-1)*vacuum%delz*den%vacz(i,ivac,ispin)
                  END DO
                  CALL qsf(vacuum%delz,sigm,vacpar(ivac),vacuum%nmz,0)
                  denz1 = den%vacz(1,ivac,ispin)          ! get estimate for potential at vacuum boundary
                  CALL xcpot%get_vxc(1,denz1,vacpot,vacxpot)
                  ! seems to be the best choice for 1D not to substract vacpar
                  IF (.NOT.oneD%odi%d1) THEN
                     vacpot = vacpot - fpi_const*vacpar(ivac)
                  END IF
                  vacpar(ivac) = vacpot(1,1)
               END DO
               IF (vacuum%nvac.EQ.1) vacpar(2) = vacpar(1)
            END IF

            enpara%enmix = 1.0


            enpara%evac0(:,ispin)=vacpar(:SIZE(enpara%evac0,1))

         END DO ! ispin
         CALL enpara%WRITE(atoms,input%jspins,input%film)
      END IF
   END IF ! mpi%irank == 0
   DEALLOCATE ( rat,eig )
   DEALLOCATE ( rh,rh1)
   DEALLOCATE ( vbar,sigm )
   DEALLOCATE ( rhoss )
   deallocate(den%vacz)
   deallocate(den%vacxy)

 END SUBROUTINE stden

END MODULE m_stden
