      MODULE m_flipcdn
!     *******************************************************
!     this subroutine reads the charge density and flips the 
!     magnetic moment within the m.t.sphere for each atom 
!     according to the variable nflip. This variable is read in
!     the main program
!             nflip = -1 : flip spin in sphere
!             nflip = -2 : scale spin by bmu(n)
!             nflip = any: no spin flip
!                            r.pentcheva,kfa,Feb'96
!     *******************************************************
      CONTAINS
        SUBROUTINE flipcdn(&
             &                   atoms,input,vacuum,sphhar,&
             &                   stars,sym,cell,&
             &                   l_noco)
          USE m_loddop
          USE m_wrtdop
          USE m_types
          IMPLICIT NONE
          TYPE(t_stars),INTENT(IN)  :: stars
          TYPE(t_vacuum),INTENT(IN) :: vacuum
          TYPE(t_atoms),INTENT(IN)  :: atoms
          TYPE(t_sphhar),INTENT(IN) :: sphhar
          TYPE(t_input),INTENT(IN)  :: input
          TYPE(t_sym),INTENT(IN)    :: sym
          TYPE(t_cell),INTENT(IN)   :: cell
          LOGICAL,INTENT(IN)        :: l_noco

          !     .. Local Scalars ..
          REAL    rhodummy,rhodumms
          INTEGER i,iter,n,nt,j,lh,na ,mp,ispin,n_ldau,urec,itype,m
            
          CHARACTER(len=8) iop,dop
          LOGICAL n_exist
          !     ..
          !     .. Local Arrays ..
          COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:),qpw(:,:),rhtxy(:,:,:,:)
          REAL   , ALLOCATABLE :: rho(:,:,:,:),rht(:,:,:)
          COMPLEX, ALLOCATABLE :: cdom(:),cdomvz(:,:),cdomvxy(:,:,:)
          CHARACTER(len=80), ALLOCATABLE :: clines(:)
          CHARACTER(len=8) name(10)
          !     ..
          !atoms%jmtd = MAXVAL(atoms%jri(:))
          !sphhar%nlhd = MAXVAL(sphhar%nlh(:))
          ALLOCATE ( rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,input%jspins),qpw(stars%ng3,input%jspins) )
          ALLOCATE ( rhtxy(vacuum%nmzxy,stars%ng2-1,2,input%jspins),rht(vacuum%nmz,2,input%jspins) )
          IF (l_noco) THEN
             ALLOCATE( cdom(stars%ng3) )
             ALLOCATE( cdomvz(vacuum%nmz,2),cdomvxy(vacuum%nmzxy,stars%ng2-1,2) )
          ENDIF

          nt = 71                   !gs, see sub mix
          IF (l_noco) THEN
             OPEN (nt,file='rhomat_inp',form='unformatted',status='old')
          ELSE
             OPEN (nt,file='cdn1',form='unformatted',status='old')
          ENDIF
          !     ---> read the charge density 
          CALL loddop(&
               &            stars,vacuum,atoms,sphhar,input,sym,&
               &            nt,&
               &            iter,rho,qpw,rht,rhtxy)
          IF (l_noco) THEN
             READ (nt) (cdom(j),j=1,stars%ng3)
             IF (input%film) THEN
                READ (nt) ((cdomvz(n,i),n=1,vacuum%nmz),i=1,vacuum%nvac)
                READ (nt) (((cdomvxy(n,j-1,i),n=1,vacuum%nmzxy),j=2,stars%ng2),i=1,vacuum%nvac) 
             ENDIF
          ENDIF
          !     ---> flip cdn for each atom with nflip=-1
          !
          na = 1
          DO n = 1, atoms%ntype
             IF (atoms%nflip(n).EQ.-1) THEN
                !     ---> spherical and non-spherical m.t. charge density
                DO lh = 0,sphhar%nlh(atoms%ntypsy(na))
                   DO j = 1,atoms%jri(n)
                      rhodummy = rho(j,lh,n,1)
                      rho(j,lh,n,1) = rho(j,lh,n,input%jspins)
                      rho(j,lh,n,input%jspins) = rhodummy
                   ENDDO
                ENDDO
             ELSEIF (atoms%nflip(n).EQ.-2) THEN
                DO lh = 0,sphhar%nlh(atoms%ntypsy(na))
                   DO j = 1,atoms%jri(n)
                      rhodummy = rho(j,lh,n,1) + rho(j,lh,n,input%jspins)
                      rhodumms = rho(j,lh,n,1) - rho(j,lh,n,input%jspins)
                      rho(j,lh,n,1) = 0.5 * ( rhodummy + atoms%bmu(n) * rhodumms )
                      rho(j,lh,n,input%jspins) = 0.5*(rhodummy - atoms%bmu(n)*rhodumms )
                   ENDDO
                ENDDO
             END IF
             na = na + atoms%neq(n)
          ENDDO
          !     ----> write the spin-polarized density on unit nt
          REWIND nt
          CALL wrtdop(&
               &            stars,vacuum,atoms,sphhar,input,sym,&
               &            nt,&
               &            iter,rho,qpw,rht,rhtxy)
          IF (l_noco) THEN
             WRITE (nt) (cdom(j),j=1,stars%ng3)
             IF (input%film) THEN
                WRITE (nt) ((cdomvz(n,i),n=1,vacuum%nmz),i=1,vacuum%nvac)
                WRITE (nt) (((cdomvxy(n,j-1,i),n=1,vacuum%nmzxy),j=2,stars%ng2),i=1,vacuum%nvac)
             ENDIF
          ENDIF
          CLOSE (nt)
          !
          ! for lda+U: flip n-matrix 
          !
          IF (atoms%n_u.GT.0) THEN
             INQUIRE (file='n_mmp_mat',exist=n_exist)
             IF (n_exist) THEN
                OPEN (69,file='n_mmp_mat',status='old',form='formatted')
                ALLOCATE (  n_mmp(-3:3,-3:3,atoms%n_u,2) )

                READ (69,9000) n_mmp
                !   flip    ...
                n_ldau = 0
                DO n = 1,atoms%ntype
                   IF (atoms%lda_u(n)%l.GE.0) THEN
                      n_ldau = n_ldau + 1
                      IF (atoms%nflip(n).EQ.-1) THEN
                         DO m = -3,3
                            DO mp = -3,3
                               rhodummy = n_mmp(m,mp,n_ldau,1)
                               n_mmp(m,mp,n_ldau,1) = n_mmp(m,mp,n_ldau,input%jspins)
                               n_mmp(m,mp,n_ldau,input%jspins) = rhodummy
                            ENDDO
                         ENDDO
                      ELSEIF (atoms%nflip(n).EQ.-2) THEN
                         DO m = -3,3
                            DO mp = -3,3
                               rhodummy = n_mmp(m,mp,n_ldau,1) + &
                                    &                           n_mmp(m,mp,n_ldau,input%jspins)
                               rhodumms = n_mmp(m,mp,n_ldau,1) - &
                                    &                           n_mmp(m,mp,n_ldau,input%jspins)
                               n_mmp(m,mp,n_ldau,1) = 0.5 * ( rhodummy + &
                                    &                                      atoms%bmu(n) * rhodumms )
                               n_mmp(m,mp,n_ldau,input%jspins) = 0.5*( rhodummy - &
                                    &                                         atoms%bmu(n) * rhodumms )
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO
                !   flip    ...
                REWIND (69)
                WRITE (69,9000) n_mmp
9000            FORMAT(7f20.13)
                !
                DEALLOCATE ( n_mmp )
             ENDIF
          ENDIF
          !-lda+U
          !
          !--->   read enpara and  flip lines
          !
          INQUIRE(file='enpara',exist=n_exist)
          IF (n_exist) THEN
             OPEN(40,file ='enpara',status='old',form='formatted')

             n = 2
             DO itype = 1 , atoms%ntype
                n         = n + 1
                IF (atoms%nlo(itype)>0) n = n + 2
             ENDDO
             IF (input%film) n = n + 1
             ALLOCATE (clines(2*n))
             DO i = 1,2*n
                READ (40,'(a)') clines(i)
             ENDDO

             REWIND 40
             i = 0 
             DO ispin = 1,input%jspins
                i         = i + 2
                WRITE (40,'(a)') TRIM(clines(i-1))
                WRITE (40,'(a)') TRIM(clines(i))
                DO itype = 1 , atoms%ntype
                   i          = i + 1
                   m               = i
                   IF (atoms%nflip(itype)==-1) m = MOD(i+n,2*n)
                   IF (m==0) m        = 2*n
                   WRITE (40,'(a)') TRIM(clines(m))
                   IF (atoms%nlo(itype)>0) THEN
                      WRITE (40,'(a)') TRIM(clines(m+1))
                      WRITE (40,'(a)') TRIM(clines(m+2))
                      i = i + 2
                   ENDIF
                ENDDO
                IF (input%film) THEN
                   i = i + 1
                   WRITE (40,'(a)') TRIM(clines(i))
                ENDIF
             ENDDO

             DEALLOCATE (clines,rho,qpw,rhtxy,rht)
             IF (l_noco) THEN
                DEALLOCATE (cdom,cdomvz,cdomvxy)
             ENDIF
             CLOSE(40)
          ENDIF
        END SUBROUTINE flipcdn
      END MODULE m_flipcdn
