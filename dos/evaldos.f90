      MODULE m_evaldos
      CONTAINS
      SUBROUTINE evaldos(eig_id,input,banddos,vacuum,kpts,atoms,sym,noco,oneD,cell,&
           dimension,efermiarg, l_mcd,ncored,ncore,e_mcd,nsld)
!----------------------------------------------------------------------
!
!     vk: k-vectors
!     wk weight of k-point (not used)
!     nevk no of eigenvalues
!     ev eigenvalue
!     qal partial charges
!             partial charge interstitial: qal(lmax*ntype+1...
!             partial charge vacuum      : qal(lmax*ntype+2...
!     qlay,qstar read in vacuum charges
!     qval partial charges in the vacuum
!             qval(m*nstars,neigd,nkptd):charge in m'th layer
!             qval(m*nstars+nstars,... ):star-resolved charge of that layer
!             qval(layers*nstars+....  ):same for second vacuum
!     ntb=max(nevk)
!
!----------------------------------------------------------------------
      USE m_eig66_io,ONLY:read_dos,read_eig
      USE m_triang
      USE m_maketetra
      USE m_tetrados
      USE m_dosbin
      USE m_ptdos
      USE m_smooth
      USE m_types
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: eig_id
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_banddos),INTENT(IN)     :: banddos
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_noco),INTENT(IN)        :: noco
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_kpts),INTENT(IN)        :: kpts
      TYPE(t_atoms),INTENT(IN)       :: atoms

      INTEGER, INTENT(IN) :: ncored 
      INTEGER, INTENT(IN) :: nsld
      REAL,    INTENT(IN) :: efermiarg
      LOGICAL, INTENT(IN) :: l_mcd 

      INTEGER, INTENT(IN) :: ncore(:)!(ntype)
      REAL,    INTENT(IN) :: e_mcd(atoms%ntypd,input%jspins,ncored) 
!-odim
!+odim 
!    locals
      INTEGER, PARAMETER ::  lmax= 4, ned = 1301
      REAL,    PARAMETER :: factor = 27.2
      INTEGER  i,s,v,index,jspin,k,l,l1,l2,ln,n,nl,ntb,ntria,ntetra
      INTEGER  icore,qdim,n_orb
      REAL     as,de,efermi,emax,emin,qmt,sigma,totdos
      REAL     e_up,e_lo,e_test1,e_test2,fac,sumwei,dk
      LOGICAL  l_tria,l_orbcomp

      INTEGER  itria(3,2*kpts%nkpt),nevk(kpts%nkpt),itetra(4,6*kpts%nkpt)
      INTEGER, ALLOCATABLE :: ksym(:),jsym(:)
      REAL     vk(3,kpts%nkptd),wt(kpts%nkptd),voltet(6*kpts%nkpt),kx(kpts%nkpt),vkr(3,kpts%nkpt)
      REAL     ev(dimension%neigd,kpts%nkptd),e(ned),gpart(ned,atoms%ntype),atr(2*kpts%nkpt)
      REAL     e_grid(ned+1),spect(ned,3*atoms%ntypd),ferwe(dimension%neigd,kpts%nkptd)
      REAL,    ALLOCATABLE :: qal(:,:,:),qval(:,:,:),qlay(:,:,:),g(:,:),qal_tmp(:,:,:),qis(:),qvlay(:,:,:)
      REAL,    ALLOCATABLE :: mcd(:,:,:),orbcomp(:,:,:),qmtp(:,:)
      REAL,    ALLOCATABLE :: qintsl(:,:),qmtsl(:,:),qvac(:,:)
      COMPLEX, ALLOCATABLE :: qstars(:,:,:,:)
      CHARACTER(len=2) :: spin12(2),ch_mcd(3)
      CHARACTER(len=8) :: chntype*2,chform*19
      DATA spin12/'.1' , '.2'/
      DATA ch_mcd/'.+' , '.-' , '.0'/
!
      IF  (banddos%ndir.NE.-3) THEN
        qdim = lmax*atoms%ntype+3
      ELSE
        qdim = 2*nsld 
        n_orb = 0
        l_orbcomp = .false.
        INQUIRE(file='orbcomp',exist=l_orbcomp)
        IF (l_orbcomp) THEN
          OPEN (1,file='orbcomp',form='formatted')
          READ (1,*) n_orb
          WRITE (*,*) 'DOS: orbcomp',n_orb
          CLOSE (1)
          qdim = 23
        ENDIF
      ENDIF
      ALLOCATE( qal(qdim,dimension%neigd,kpts%nkptd),&
     &          qval(vacuum%nstars*vacuum%layers*vacuum%nvac,dimension%neigd,kpts%nkptd),&
     &          qlay(dimension%neigd,vacuum%layerd,2),qstars(vacuum%nstars,dimension%neigd,vacuum%layerd,2))
      IF (l_mcd) ALLOCATE( mcd(3*atoms%ntypd*ncored,dimension%neigd,kpts%nkptd) )
!
! scale energies
      sigma = banddos%sig_dos*factor
      emin =min(banddos%e1_dos*factor,banddos%e2_dos*factor)
      emax =max(banddos%e1_dos*factor,banddos%e2_dos*factor)
      efermi = efermiarg*factor
 
      WRITE (6,'(a)') 'DOS-Output is generated!'

      IF ( NINT((emax - emin)/sigma) > ned ) THEN
        WRITE(6,*) 'sig_dos too small for DOS smoothing:'   
        WRITE(6,*) 'Reduce energy window or enlarge banddos%sig_dos!'
        WRITE(6,*) 'For now: setting sigma to zero !'
        sigma = 0.0
      ENDIF

      WRITE (6,*) 'sigma=   ' , sigma
      WRITE (6,*) 'emax=   ' , emax
      WRITE (6,*) 'emin=   ' , emin
      WRITE (6,*) 'ef_inp=   ' , efermi
!
!     create energy grid
      emax = emax - efermi
      emin = emin - efermi
      de = (emax-emin)/(ned-1)
      DO i=1,ned
         e(i) = emin + (i-1)*de
      ENDDO
 
      IF ( l_mcd ) THEN ! create an energy grid for mcd-spectra
        e_lo =  9.9d+9 
        e_up = -9.9d+9     
        DO jspin = 1,input%jspins
          DO n = 1,atoms%ntype
            DO icore = 1 , ncore(n)
              e_lo = min(e_mcd(n,jspin,icore),e_lo)
              e_up = max(e_mcd(n,jspin,icore),e_up)
            ENDDO
          ENDDO
        ENDDO
        e_lo = e_lo*factor - efermi - emax 
        e_up = e_up*factor - efermi
        de = (e_up-e_lo)/(ned-1)
        DO i=1,ned
          e_grid(i) = e_lo + (i-1)*de
          spect(i,:) = 0.0
        ENDDO
        e_grid(ned+1) = e_lo + ned*de
      ENDIF

      DO jspin = 1,input%jspins
         ntb = 0
         DO k = 1,kpts%nkpt
!
!     initialize arrays
!
            DO n = 1,dimension%neigd
               DO i = 1,qdim
                  qal(i,n,k) = 0.
               ENDDO
               DO i = 1,vacuum%nstars*vacuum%layers*vacuum%nvac
                  qval(i,n,k) = 0.
               ENDDO
            ENDDO
!
!     read data from file!
!
            ALLOCATE( ksym(dimension%neigd),jsym(dimension%neigd) )
            ALLOCATE( qal_tmp(1:lmax,atoms%ntype,dimension%neigd))
            ALLOCATE( orbcomp(dimension%neigd,23,atoms%natd),qintsl(nsld,dimension%neigd))
            ALLOCATE( qmtsl(nsld,dimension%neigd),qmtp(dimension%neigd,atoms%natd),qvac(dimension%neigd,2))
            ALLOCATE( qis(dimension%neigd),qvlay(dimension%neigd,vacuum%layerd,2))
            CALL read_eig(eig_id,k,jspin,bk=vk(:,k),wk=wt(k),neig=nevk(k),eig=ev(:,k))
            CALL read_dos(eig_id,k,jspin,qal_tmp,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
            qal(1:lmax*atoms%ntype,:,k)=reshape(qal_tmp,(/lmax*atoms%ntype,size(qal_tmp,3)/))
            qal(lmax*atoms%ntype+2,:,k)=qvac(:,1)         ! vacuum 1
            qal(lmax*atoms%ntype+3,:,k)=qvac(:,2)         ! vacuum 2
            qal(lmax*atoms%ntype+1,:,k)=qis              ! interstitial
            DEALLOCATE( ksym,jsym )
            IF (l_orbcomp)THEN
               IF (n_orb == 0) THEN
                  qal(1:nsld,:,k)        = qintsl(:,:)
                  qal(nsld+1:2*nsld,:,k) = qmtsl(:,:)
               ELSE 
                  DO i = 1, 23
                     DO l = 1, nevk(k)
                        qal(i,l,k) = orbcomp(l,i,n_orb)*qmtp(l,n_orb)/10000.
                     ENDDO
                     DO l = nevk(k)+1, dimension%neigd
                        qal(i,l,k) = 0.0
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            DEALLOCATE( orbcomp,qintsl,qmtsl,qmtp,qvac,qis,qal_tmp,qvlay)
            ntb = max(ntb,nevk(k))
!
!     set vacuum partial charge zero, if bulk calculation
!     otherwise, write vacuum charge in correct arrays
!
            IF ((.NOT.input%film).AND.(banddos%ndir.NE.-3)) THEN
               DO n = 1,dimension%neigd
                  qal(lmax*atoms%ntype+2,n,k) = 0.0
                  qal(lmax*atoms%ntype+3,n,k) = 0.0
               ENDDO
            ELSEIF ( banddos%vacdos .and. input%film ) THEN
               DO i = 1,nevk(k)
                  DO v = 1,vacuum%nvac
                     DO l = 1,vacuum%layers
                        index = (l-1)*vacuum%nstars + (v-1)*(vacuum%nstars*vacuum%layers) + 1
                        qval(index,i,k) = qlay(i,l,v)
                        DO s = 1,vacuum%nstars - 1
                           qval(index+s,i,k) = real(qstars(s,i,l,v))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
!
!     calculate interstitial dos if not noco
!     in the noco case, qis has been calculated in pwden and is read in from tmp_dos
!
            IF ((.NOT.noco%l_noco).AND.(banddos%ndir.NE.-3)) THEN
               DO i = 1 , dimension%neigd
                  qal(lmax*atoms%ntype+1,i,k) = 1.
                  DO nl = 1 , atoms%ntype
                     l1 = lmax*(nl-1) + 1
                     l2 = lmax*nl
                     qmt=0.0
                     DO l = l1 , l2
                        qmt = qmt + qal(l,i,k)*atoms%neq(nl)
                     ENDDO
                     qal(lmax*atoms%ntype+1,i,k) = qal(lmax*atoms%ntype+1,i,k) - qmt
                  ENDDO
                  qal(lmax*atoms%ntype+1,i,k) = qal(lmax*atoms%ntype+1,i,k)&
                       -qal(lmax*atoms%ntype+2,i,k)*(3-vacuum%nvac) -qal(lmax*atoms%ntype+3,i,k)*(vacuum%nvac-1) 
               ENDDO
            ENDIF
!
!---- >     convert eigenvalues to ev and shift them by efermi
!
            DO i = 1 , nevk(k)
               ev(i,k) = ev(i,k)*factor - efermi
            ENDDO
            DO i = nevk(k) + 1, dimension%neigd
               ev(i,k) = 9.9e+99
            ENDDO
!
!
         ENDDO                                                 ! end of k-point loop
!
!     calculate the triangles!
!
         IF ( jspin.EQ.1 ) THEN
           l_tria=.true.
           IF (input%film .AND. .NOT.oneD%odi%d1) THEN
             CALL triang(vk,kpts%nkpt,itria,ntria,atr,as,l_tria)
             IF (sym%invs) THEN
               IF (abs(sym%nop2*as-0.5).GT.0.000001) l_tria=.false.
             ELSE
               IF (abs(sym%nop2*as-1.0).GT.0.000001) l_tria=.false.
             ENDIF
             write(*,*) as,sym%nop2,l_tria
!             l_tria=.true.
           ELSE
             IF (input%l_inpXML) GOTO 66
             OPEN (41,file='kpts',FORM='formatted',STATUS='old')
             DO i = 1, kpts%nkpt+1
                READ (41,*,END=66,ERR=66)
             ENDDO
             READ (41,'(i5)',END=66,ERR=66) ntetra
             READ (41,'(4(4i6,4x))') ((itetra(i,k),i=1,4),k=1,ntetra)
             READ (41,'(4f20.13)') (voltet(k),k=1,ntetra)
             CLOSE(41)
             voltet(1:ntetra) = voltet(1:ntetra) / ntetra
             l_tria=.true.
             GOTO 67
 66          CONTINUE                       ! no tetrahedron-information of file
             CALL triang(vk,kpts%nkpt,itria,ntria,atr,as,l_tria)
             l_tria=.true.
! YM: tetrahedrons is not the way in 1D
             IF (oneD%odi%d1) as = 0.0         
             IF (sym%invs) THEN
               IF (abs(sym%nop2*as-1.0).GT.0.000001) l_tria=.false.
             ELSE
               IF (abs(sym%nop2*as-0.5).GT.0.000001) l_tria=.false.
             ENDIF

             IF (l_tria) THEN
               CALL make_tetra(kpts%nkpt,vk,ntria,itria,atr,&
                    ntetra,itetra,voltet)
             ELSE
               WRITE (6,*) 'no tetrahedron method with these k-points!'
               WRITE (6,*) sym%nop2,as
             ENDIF
 67          CONTINUE                       ! tetrahedron-information read or created
           ENDIF
         ENDIF
!
        IF ( .not.l_mcd ) THEN
         ALLOCATE (g(ned,qdim))
        ELSE
         ALLOCATE (g(ned,3*atoms%ntypd*ncored))
        ENDIF
!
         IF ( l_tria.and.(.not.l_mcd).and.(banddos%ndir.NE.-3) ) THEN
!
!     DOS calculation: use triangular method!!
!
            IF ( input%film ) THEN
!             CALL ptdos(
!    >                  emin,emax,jspins,ned,qdim,neigd,
!    >                  ntria,as,atr,2*nkpt,itria,nkpt,ev,qal,e,
!    <                  g)
              CALL ptdos(emin,emax,input%jspins,ned,qdim,ntb,ntria,as,&
                        atr,2*kpts%nkpt,itria,kpts%nkpt,ev(1:ntb,1:kpts%nkpt),&
                        qal(:,1:ntb,1:kpts%nkpt),e, g)
            ELSE
              write(*,*) efermi
              CALL tetra_dos(lmax,atoms%ntype,dimension%neigd,ned,ntetra,kpts%nkpt,&
                            itetra,efermi,voltet,e,nevk, ev,qal, g)
              IF (input%jspins.EQ.1) g(:,:) = 2 * g(:,:)
            ENDIF
         ELSE
!
!     DOS calculation: use histogram method
!
            IF ( .not.l_mcd ) THEN
            CALL dos_bin(input%jspins,qdim,ned,emin,emax,dimension%neigd,kpts%nkpt,&
                 nevk,wt,ev,qal, g)
            ELSE
            CALL dos_bin(input%jspins,3*atoms%ntypd*ncored,ned,emin,emax,ntb,kpts%nkpt,&
                 nevk(1:kpts%nkpt),wt(1:kpts%nkpt),ev(1:ntb,1:kpts%nkpt), mcd(1:3*atoms%ntypd*ncored,1:ntb,1:kpts%nkpt), g)
            ENDIF
         ENDIF
!
!---- >     smoothening
!
         IF ( .not.l_mcd ) THEN
            IF ( sigma.GT.0.0 ) THEN
              DO ln = 1 , qdim
                CALL smooth(e,g(1,ln),sigma,ned)
              ENDDO
            ENDIF
 
!*** sum up for all atoms
 
         IF (banddos%ndir.NE.-3) THEN
            DO l = 1 , atoms%ntype
               l1 = lmax*(l-1) + 1
               l2 = lmax*l
               DO i = 1 , ned
                  gpart(i,l) = 0.0
                  DO nl = l1 , l2
                     gpart(i,l) = gpart(i,l) + g(i,nl)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF (n_orb == 0) THEN
            DO l = 1 , nsld
               nl = nsld+l
               DO i = 1 , ned
                  gpart(i,l) = g(i,l) + g(i,nl)
               ENDDO
            ENDDO
         ENDIF
    
!**** write out DOS
         OPEN (18,FILE='DOS'//spin12(jspin))

         DO i = 1 , ned
           totdos = 0.0
           IF (banddos%ndir.NE.-3) THEN
             DO nl = 1 , atoms%ntype
                totdos = totdos + gpart(i,nl)*atoms%neq(nl)
             ENDDO
             totdos = totdos + g(i,lmax*atoms%ntype+1) + g(i,lmax*atoms%ntype+2) *&
                  (3 - vacuum%nvac) + g(i,lmax*atoms%ntype+3)*(vacuum%nvac - 1)
             IF (atoms%ntype < 20) THEN
                WRITE (18,99001)  e(i),totdos,g(i,lmax*atoms%ntype+1), &
                     g(i,lmax*atoms%ntype+2),g(i,lmax*atoms%ntype+3),&
                     (gpart(i,l),l=1,atoms%ntype), (g(i,l),l=1,atoms%ntype*lmax)
             ELSE
             WRITE (18,99001)  e(i),totdos,g(i,lmax*atoms%ntype+1), &
                  g(i,lmax*atoms%ntype+2),g(i,lmax*atoms%ntype+3), (gpart(i,l),l=1,atoms%ntype)
          ENDIF
       ELSEIF (n_orb == 0) THEN
          DO nl = 1 , nsld
             totdos = totdos + gpart(i,nl)
          ENDDO
          WRITE (18,99001)  e(i),totdos,(gpart(i,nl),nl=1,nsld), (g(i,l),l=1,2*nsld)
           ELSE
             DO nl = 1 , 23
                totdos = totdos + g(i,nl)
             ENDDO
             WRITE (18,99001)  e(i),totdos,(g(i,l),l=1,23)
           ENDIF
         ENDDO
         CLOSE (18)

         ELSE
           write(*,'(4f15.8)') ((e_mcd(n,jspin,i),n=1,atoms%ntype),i=1,ncored)
           write(*,*)
           write(*,'(4f15.8)') (g(800,n),n=1,3*atoms%ntype*ncored)
           write(*,*)
           write(*,'(4f15.8)') (mcd(n,10,8),n=1,3*atoms%ntype*ncored)
           DO n = 1,atoms%ntype
             DO l = 1 , ned
               DO icore = 1 , ncore(n)
                 DO i = 1 , ned-1
                   IF (e(i).GT.0) THEN     ! take unoccupied part only
                   e_test1 = -e(i) - efermi +e_mcd(n,jspin,icore)*factor
                   e_test2 = -e(i+1)-efermi +e_mcd(n,jspin,icore)*factor
                   IF ((e_test2.LE.e_grid(l)).AND. (e_test1.GT.e_grid(l))) THEN
                     fac = (e_grid(l)-e_test1)/(e_test2-e_test1)
                     DO k = 3*(n-1)+1,3*(n-1)+3
                       spect(l,k) = spect(l,k)+ g(i,3*atoms%ntypd*(icore-1)+k)&
                           *(1.-fac) + fac * g(i+1,3*atoms%ntypd*(icore-1)+k)
                     ENDDO
                   ENDIF
                   ENDIF
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
           CLOSE (18) 
         ENDIF
         DEALLOCATE (g)
!         
!------------------------------------------------------------------------------
!     now calculate the VACOS
!------------------------------------------------------------------------------
            
         IF ( banddos%vacdos .and. input%film ) THEN
            ALLOCATE(g(ned,vacuum%nstars*vacuum%layers*vacuum%nvac))
!            CALL ptdos(
!     >                 emin,emax,jspins,ned,nstars*nvac*layers,neigd,
!     >                 ntria,as,atr,2*nkpt,itria,nkptd,ev,qval,e,
!     <                 g)
            CALL ptdos(emin,emax,input%jspins,ned,vacuum%nstars*vacuum%nvac*vacuum%layers,ntb,ntria&
                ,as,atr,2*kpts%nkpt,itria,kpts%nkptd,ev(1:ntb,1:kpts%nkptd), qval(:,1:ntb,1:kpts%nkptd),e,g)
            
!---- >     smoothening
            IF ( sigma.GT.0.0 ) THEN
               DO ln = 1 , vacuum%nstars*vacuum%nvac*vacuum%layers
                  CALL smooth(e,g(1,ln),sigma,ned)
               ENDDO
            ENDIF
            
!     write VACDOS
            
            OPEN (18,FILE='VACDOS'//spin12(jspin))
!            WRITE (18,'(i2,25(2x,i3))') Layers , (Zlay(l),l=1,Layers)
            DO i = 1 , ned
             WRITE (18,99001) e(i) , (g(i,l),l=1,vacuum%Layers*vacuum%Nstars*vacuum%Nvac)
            ENDDO
            CLOSE (18)
            DEALLOCATE(g)
         ENDIF
!
!------------------------------------------------------------------------------
!     for bandstructures
!------------------------------------------------------------------------------

         IF (banddos%ndir == -4) THEN
            OPEN (18,FILE='bands'//spin12(jspin))
            ntb = minval(nevk(:))    
            kx(1) = 0.0
            vkr(:,1)=matmul(cell%bmat,vk(:,1))
            DO k = 2, kpts%nkpt
              
            vkr(:,k)=matmul(cell%bmat,vk(:,k))
              dk = ( vkr(1,k) - vkr(1,k-1) )**2 + ( vkr(2,k) - vkr(2,k-1) )**2 +&
                   ( vkr(3,k) - vkr(3,k-1) )**2
              kx(k) = kx(k-1) + sqrt(dk)
            ENDDO
            DO i = 1, ntb
              DO k = 1, kpts%nkpt
                write(18,'(2f15.9)') kx(k),ev(i,k)
              ENDDO
            ENDDO
            CLOSE (18)
         ENDIF

      ENDDO
!         
!------------------------------------------------------------------------------
!     for MCD calculations ...
!------------------------------------------------------------------------------

      IF (l_mcd) THEN
        WRITE (chntype,'(i2)') atoms%ntype+1
        chform = '('//chntype//'f15.8)'
        IF ( sigma.GT.0.0 ) THEN
           IF ( l_mcd ) THEN
             DO ln = 1 , 3*atoms%ntypd
               CALL smooth(e_grid,spect(1,ln),sigma,ned)
             ENDDO
           ENDIF
        ENDIF
        DO l = 1,3
          OPEN (18,FILE='MCD_SPEC'//ch_mcd(l))
          DO i = 1 , ned
          WRITE (18,FMT=chform) e_grid(i),(spect(i,3*(n-1)+l),n=1,atoms%ntype)
          ENDDO
          CLOSE (18)
        ENDDO
      ENDIF

      DEALLOCATE(qal,qval,qlay,qstars)
      IF (l_mcd) DEALLOCATE( mcd )
99001 FORMAT (f10.5,110(1x,e10.3))

      END SUBROUTINE evaldos
      END MODULE m_evaldos 
