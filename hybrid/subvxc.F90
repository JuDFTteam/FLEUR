MODULE m_subvxc
CONTAINS
  SUBROUTINE subvxc(lapw,bk, DIMENSION,input,jsp,vr0, atoms,usdus, hybrid, el,ello,sym,&
       nlot_d,kveclo, cell, sphhar, stars,xcpot,mpi,oneD,hamovlp,vx)


    USE m_intgr,     ONLY : intgr3
    USE m_constants
    USE m_gaunt,     ONLY : gaunt1
    USE m_wrapper
    USE m_loddop
    USE m_radflo
    USE m_radfun
    USE m_abcof3
    USE m_icorrkeys
    USE m_hybridmix
    USE m_types
    IMPLICIT NONE
    TYPE(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_lapw),INTENT(IN)      :: lapw
    TYPE(t_usdus),INTENT(INOUT)  :: usdus
    TYPE(t_potden),INTENT(IN)    :: vx
    !     .. Scalar Arguments ..

    INTEGER, INTENT (IN) :: jsp 
    INTEGER, INTENT (IN) :: nlot_d




    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: kveclo(nlot_d)

    REAL,    INTENT (IN) :: vr0(atoms%jmtd,atoms%ntype,DIMENSION%jspd)               ! just for radial functions
    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (IN) :: ello(atoms%nlod,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (IN) :: bk(3)
    TYPE(t_hamovlp),INTENT(INOUT)::hamovlp

    !     .. Local Scalars ..
    INTEGER               ::  ic,indx,m,ig1,ig2
    INTEGER               ::  nlharm,nnbas,typsym,lm
    INTEGER               ::  noded,nodeu
    INTEGER               ::  nbasf0
    INTEGER               ::  i,j,l,ll,l1,l2 ,m1,m2  ,j1,j2
    INTEGER               ::  ok,p1,p2,lh,mh,pp1,pp2
    INTEGER               ::  igrid,itype,ilharm,istar
    INTEGER               ::  ineq,iatom,ilo,ilop,ieq,icentry
    INTEGER               ::  ikvecat,ikvecprevat,invsfct,ikvec,ikvecp
    INTEGER               ::  lp,mp,pp
    REAL                  ::  a_ex
    REAL                  ::  wronk
    COMPLEX               ::  rc,rr

    !     .. Local Arrays ..
    INTEGER               ::  gg(3)
    INTEGER               ::  pointer_lo(atoms%nlod,atoms%ntype)

    REAL                  ::  integ(0:sphhar%nlhd,hybrid%maxindx,0:atoms%lmaxd,hybrid%maxindx,0:atoms%lmaxd)
    REAL                  ::  grid(atoms%jmtd)
    REAL                  ::  vr(atoms%jmtd,0:sphhar%nlhd)
    REAL                  ::  f(atoms%jmtd,2,0:atoms%lmaxd),g(atoms%jmtd,2,0:atoms%lmaxd)
    REAL                  ::  flo(atoms%jmtd,2,atoms%nlod)
    REAL                  ::  uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
    REAL                  ::  ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)


    REAL                  ::  bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
         bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)

    COMPLEX               ::  vpw(stars%ng3)
    COMPLEX               ::  vxc(hamovlp%matsize)
    COMPLEX               ::  vrmat(hybrid%maxlmindx,hybrid%maxlmindx)
    COMPLEX               ::  carr(hybrid%maxlmindx,DIMENSION%nvd),carr1(DIMENSION%nvd,DIMENSION%nvd)
    COMPLEX ,ALLOCATABLE  ::  ahlp(:,:,:),bhlp(:,:,:)
    COMPLEX, ALLOCATABLE  ::  bascof(:,:,:)
    COMPLEX               ::  bascof_lo(3,-atoms%llod:atoms%llod,4*atoms%llod+2,atoms%nlod, atoms%nat)



    CALL timestart("subvxc")
    vxc=0


    !  calculate radial functions
    hybrid%nindx      = 2
    DO itype = 1,atoms%ntype


       !
       !--->    generate the radial basis-functions for each l
       !

       WRITE(6,'(a,i3,a)') new_LINE('n')//new_LINE('n')//' wavefunction parameters for atom type',itype,':'
       WRITE(6,'(31x,a,32x,a)') 'radial function','energy derivative'
       WRITE(6,'(a)') '  l    energy            value        '//&
            'derivative    nodes          value        derivative    nodes       norm        wronskian'
       DO l = 0,atoms%lmax(itype)
          CALL radfun(l,itype,jsp,el(l,itype,jsp),vr0(1,itype,jsp),atoms,&
               f(1,1,l),g(1,1,l),usdus,&
               nodeu,noded,wronk)
          WRITE (6,FMT=8010) l,el(l,itype,jsp),usdus%us(l,itype,jsp),&
               usdus%dus(l,itype,jsp),nodeu,usdus%uds(l,itype,jsp),usdus%duds(l,itype,jsp),noded,&
               usdus%ddn(l,itype,jsp),wronk
       END DO
       !  8000    FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',
       !      +          /,t32,'radial function',t79,'energy derivative',/,t3,
       !      +          'l',t8,'energy',t26,'value',t39,'derivative',t53,
       !      +          'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,
       !      +          'norm',t119,'wronskian')
8010   FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

       bas1(:,1,:,itype)=f(:,1,:)
       bas1(:,2,:,itype)=g(:,1,:)
       bas2(:,1,:,itype)=f(:,2,:)
       bas2(:,2,:,itype)=g(:,2,:)

       !
       !--->   generate the extra radial basis-functions for the local orbitals,
       !--->   if there are any.
       !        
       IF (atoms%nlo(itype).GE.1) THEN

          CALL radflo(&
               atoms,itype,jsp,ello(1,1,jsp),vr0(1,itype,jsp),&
               f,g,mpi,&
               usdus,&
               uuilon,duilon,ulouilopn,flo,.TRUE.)

          DO i=1,atoms%nlo(itype)
             hybrid%nindx(atoms%llo(i,itype),itype) = hybrid%nindx(atoms%llo(i,itype),itype) + 1
             pointer_lo(i,itype)       = hybrid%nindx(atoms%llo(i,itype),itype)
             bas1(:,hybrid%nindx(atoms%llo(i,itype),itype),atoms%llo(i,itype),itype)=&
                  flo(:,1,i)
             bas2(:,hybrid%nindx(atoms%llo(i,itype),itype),atoms%llo(i,itype),itype)=&
                  flo(:,2,i)
          END DO
       END IF
    END DO


    ! compute APW coefficients

    !  calculate bascof
    ALLOCATE( ahlp(DIMENSION%nvd,0:DIMENSION%lmd,atoms%nat),bhlp(DIMENSION%nvd,0:DIMENSION%lmd,atoms%nat),stat=ok)
    IF( ok .NE. 0 ) STOP 'subvxc: error in allocation of ahlp/bhlp'

    CALL abcof3( input,atoms,sym,jsp,cell, bk,lapw,&
         usdus, kveclo,oneD,ahlp,bhlp,bascof_lo)

    ALLOCATE( bascof(DIMENSION%nvd,2*(DIMENSION%lmd+1),atoms%nat), stat=ok )
    IF( ok .NE. 0 ) STOP 'subvxc: error in allocation of bascof'
    bascof = 0
    ic     = 0

    DO itype=1,atoms%ntype
       DO ieq=1,atoms%neq(itype)
          ic   = ic + 1
          indx = 0
          DO l=0,atoms%lmax(itype)
             ll = l*(l+1)
             DO M=-l,l
                lm=ll+M 
                DO i=1,2
                   indx = indx + 1
                   IF( i .EQ. 1) THEN
                      bascof(:,indx,ic) = ahlp(:,lm,ic)
                   ELSE IF( i .EQ. 2 ) THEN
                      bascof(:,indx,ic) = bhlp(:,lm,ic)
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE( ahlp,bhlp )

    ! Loop over atom types
    iatom = 0
    DO itype = 1,atoms%ntype

       typsym = atoms%ntypsy( SUM(atoms%neq(:itype-1))+1 )
       nlharm = sphhar%nlh(typsym)

       ! Calculate vxc = vtot - vcoul
       DO l=0,nlharm
          DO i=1,atoms%jri(itype)
             IF(l.EQ.0) THEN
                !               vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) -  vrcou(i,0,itype,jsp)   
                vr(i,0)=  vx%mt(i,0,itype,jsp)*sfp_const/atoms%rmsh(i,itype)  !
             ELSE                                              ! vxc = vtot - vcoul
                !               vr(i,l)= vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                vr(i,l)=  vx%mt(i,l,itype,jsp)      
             END IF
          END DO
       END DO


       ! Calculate MT contribution to vxc matrix elements
       ! Precompute auxiliary radial integrals
       DO ilharm = 0,nlharm
          i = 0
          DO l1 = 0,atoms%lmax(itype)
             DO p1 = 1,2
                i = i + 1
                j = 0
                DO l2 = 0,atoms%lmax(itype)     
                   DO p2 = 1,2
                      j = j + 1
                      IF( j .LE. i) THEN
                         DO igrid = 1,atoms%jri(itype)
                            grid(igrid)=vr(igrid,ilharm)*(bas1(igrid,p1,l1,itype)*bas1(igrid,p2,l2,itype)+ bas2(igrid,p1,l1,itype)*bas2(igrid,p2,l2,itype) )
                         END DO

                         CALL intgr3(grid,atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),integ(ilharm,p1,l1,p2,l2) ) ! numerical integration

                         integ(ilharm,p2,l2,p1,l1)=integ(ilharm,p1,l1,p2,l2)
                      END IF
                   END DO
                END DO

             END DO
          END DO
       END DO

       ! Calculate muffin tin contribution to vxc matrix
       vrmat=0

       j1=0
       DO l1 = 0,atoms%lmax(itype) ! loop: left basis function
          DO m1 = -l1,l1
             DO p1 = 1,2
                j1 = j1+1
                j2 = 0
                DO l2 = 0,atoms%lmax(itype) ! loop: right basis function
                   DO m2 = -l2,l2
                      DO p2 = 1,2
                         j2 = j2+1
                         rr = 0
                         DO ilharm = 0,nlharm ! loop: lattice harmonics of vxc
                            l = sphhar%llh(ilharm,typsym)
                            DO i = 1,sphhar%nmem(ilharm,typsym)
                               M  = sphhar%mlh(i,ilharm,typsym)
                               rc = sphhar%clnu(i,ilharm,typsym)* gaunt1(l1,l,l2,m1,M,m2,atoms%lmaxd)
                               rr = rr+integ(ilharm,p1,l1,p2,l2)*rc
                            END DO
                         END DO

                         rc           = CMPLX(0,1)**(l2-l1) ! adjusts to a/b/ccof-scaling
                         vrmat(j1,j2) = rr*rc

                      END DO
                   END DO
                END DO

             END DO
          END DO
       END DO
       nnbas = j1

       !        ! Project on bascof
       DO ineq = 1,atoms%neq(itype)
          iatom = iatom+1

          carr (:nnbas,:lapw%nv(jsp)) = CONJG(MATMUL(vrmat(:nnbas,:nnbas), TRANSPOSE(bascof(:lapw%nv(jsp),:nnbas,iatom)) ))

          carr1(:lapw%nv(jsp),:lapw%nv(jsp)) = MATMUL(bascof(:lapw%nv(jsp),:nnbas,iatom),carr(:nnbas,:lapw%nv(jsp)) )
          ic    = 0
          DO j = 1,lapw%nv(jsp)
             !            carr(:nnbas) =  matmul(vrmat(:nnbas,:nnbas),
             !     +                             bascof(j,:nnbas,iatom) )
             DO i = 1,j
                ic = ic + 1
                vxc(ic) = vxc(ic) + carr1(i,j)
                !             vxc(ic) = vxc(ic) + conjg(dotprod ( bascof(i,:nnbas,iatom),
                !     +                                           carr(:nnbas) ))
             END DO
          END DO
       END DO

    END DO ! End loop over atom types

    ! ---------------------------------------------------------------
    ! Calculate plane wave contribution
    DO i=1,stars%ng3
       vpw(i)= vx%pw(i,jsp) 
       !         vpw(i)=vpwtot(i)-vpwcou(i,jsp)      
    END DO

    ! Calculate vxc-matrix,  left basis function (ig1)
    !                        right basis function (ig2)
    ic = 0
    DO ig1=1,lapw%nv(jsp)
       DO ig2=1,ig1
          ic = ic + 1
          gg(1)=lapw%k1(ig1,jsp)-lapw%k1(ig2,jsp)
          gg(2)=lapw%k2(ig1,jsp)-lapw%k2(ig2,jsp)
          gg(3)=lapw%k3(ig1,jsp)-lapw%k3(ig2,jsp)
          istar=stars%ig(gg(1),gg(2),gg(3))
          IF(istar.NE.0) THEN
             vxc(ic)= vxc(ic) + stars%rgphs(gg(1),gg(2),gg(3))*vpw(istar)
          ELSE
             IF ( mpi%irank == 0 ) WRITE(6,'(A,/6I5)') 'Warning: Gi-Gj not in any star:',&
                  lapw%k1(ig1,jsp),lapw%k2(ig1,jsp),lapw%k3(ig1,jsp),&
                  lapw%k1(ig2,jsp),lapw%k2(ig2,jsp),lapw%k3(ig2,jsp)
          ENDIF
       ENDDO
    ENDDO

    !    
    ! -------------------------------------------------------------------
    ! Calculate local orbital contribution

    IF( ANY( atoms%nlo .NE. 0) ) THEN 

       nbasf0      = lapw%nv(jsp)*(lapw%nv(jsp)+1)/2    ! number of pure APW contributions
       icentry     = nbasf0                   ! icentry counts the entry in the matrix vxc
       iatom       = 0
       ikvecat     = 0
       ikvecprevat = 0

       DO itype=1,atoms%ntype

          typsym = atoms%ntypsy(SUM(atoms%neq(:itype-1))+1)
          nlharm = sphhar%nlh(typsym)

          ! Calculate vxc = vtot - vcoul
          DO l=0,nlharm
             DO i=1,atoms%jri(itype)
                IF(l.EQ.0) THEN
                   !                 vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) -  vrcou(i,0,itype,jsp)
                   vr(i,0)=  vx%mt(i,0,itype,jsp)*sfp_const/atoms%rmsh(i,itype)  !
                ELSE                                              ! vxc = vtot - vcoul
                   vr(i,l)=  vx%mt(i,l,itype,jsp)                    !
                   !                 vr(i,l)=  vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                END IF
             END DO
          END DO

          ! Precompute auxiliary radial integrals
          DO ilharm=0,nlharm
             i = 0
             DO l1=0,atoms%lmax(itype)
                DO p1=1,hybrid%nindx(l1,itype)
                   i = i + 1
                   j = 0
                   DO l2=0,atoms%lmax(itype)     
                      DO p2=1,hybrid%nindx(l2,itype)
                         j = j + 1
                         IF( j .LE. i) THEN
                            DO igrid=1,atoms%jri(itype)
                               grid(igrid)=vr(igrid,ilharm)*(bas1(igrid,p1,l1,itype)*bas1(igrid,p2,l2,itype)+ bas2(igrid,p1,l1,itype)*bas2(igrid,p2,l2,itype))
                            END DO

                            CALL intgr3(grid,atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),integ(ilharm,p1,l1,p2,l2) ) ! numerical integration

                            integ(ilharm,p2,l2,p1,l1) = integ(ilharm,p1,l1,p2,l2)
                         END IF
                      END DO
                   END DO

                END DO
             END DO
          END DO


          DO ieq = 1,atoms%neq(itype)
             iatom = iatom + 1
             IF( (atoms%invsat(iatom).EQ.0) .OR. (atoms%invsat(iatom) .EQ. 1) ) THEN

                IF( atoms%invsat(iatom) .EQ. 0 ) invsfct = 1
                IF( atoms%invsat(iatom) .EQ. 1 ) invsfct = 2


                DO ilo = 1,atoms%nlo(itype)
                   l1 = atoms%llo(ilo,itype)
                   DO ikvec = 1,invsfct*(2*l1+1)

                      DO m1 = -l1,l1
                         DO p1 = 1,3
                            IF( p1 .EQ. 3) THEN
                               pp1 = pointer_lo(ilo,itype)
                            ELSE
                               pp1 = p1
                            END IF

                            IF( hybrid%nindx(l1,itype) .LE. 2) STOP 'subvxc: error hybrid%nindx'

                            lm = 0

                            !loop over APW
                            DO l2 = 0,atoms%lmax(itype)
                               DO m2 = -l2,l2
                                  DO p2 = 1,2
                                     lm = lm + 1

                                     rr = 0
                                     DO ilharm = 0,nlharm
                                        lh = sphhar%llh(ilharm,typsym)
                                        DO i = 1,sphhar%nmem(ilharm,typsym)
                                           mh = sphhar%mlh(i,ilharm,typsym)
                                           rc = sphhar%clnu(i,ilharm,typsym)* gaunt1(l1,lh,l2,m1,mh,m2,atoms%lmaxd)
                                           rr = rr+integ(ilharm,p2,l2,pp1,l1)*rc
                                        END DO
                                     END DO

                                     rc = CMPLX(0d0,1d0)**(l2-l1) ! adjusts to a/b/ccof-scaling

                                     ! ic counts the entry in vxc
                                     ic = icentry
                                     DO i=1,lapw%nv(jsp)
                                        ic = ic + 1
                                        IF (hamovlp%l_real) THEN
                                           vxc(ic) = vxc(ic) + invsfct * REAL(rr*rc*bascof(i,lm,iatom) * CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom)))
                                        ELSE
                                           vxc(ic) = vxc(ic) + rr*rc*bascof(i,lm,iatom) *CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom))
                                        ENDIF
                                     END DO

                                  END DO  !p2
                               END DO  ! m2
                            END DO ! l2 ->  loop over APW


                            ! calcualte matrix-elements with local orbitals at the same atom
                            IF( ic .NE. icentry + lapw%nv(jsp) ) STOP 'subvxc: error counting ic'

                            ic = ic + ikvecprevat

                            DO ilop = 1,(ilo-1)
                               lp = atoms%llo(ilop,itype)

                               DO ikvecp = 1,invsfct*(2*lp+1)

                                  ic = ic + 1

                                  DO mp = -lp,lp
                                     DO pp = 1,3
                                        IF ( pp .EQ. 3) THEN
                                           pp2 = pointer_lo(ilop,itype)
                                        ELSE
                                           pp2 = pp
                                        END IF

                                        rr = 0
                                        DO ilharm = 0,nlharm
                                           lh = sphhar%llh(ilharm,typsym)
                                           DO i = 1,sphhar%nmem(ilharm,typsym)
                                              mh = sphhar%mlh(i,ilharm,typsym)
                                              rc = sphhar%clnu(i,ilharm,typsym)* gaunt1(l1,lh,lp,m1,mh,mp,atoms%lmaxd)
                                              rr = rr+integ(ilharm,pp2,lp,pp1,l1)*rc
                                           END DO
                                        END DO

                                        rc = CMPLX(0d0,1d0)**(lp-l1) ! adjusts to a/b/ccof-scaling

                                        IF (hamovlp%l_real) THEN

                                           vxc(ic) = vxc(ic) + invsfct * REAL( rr*rc*bascof_lo(pp,mp, ikvecp,ilop,iatom) * CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom)) )
                                        ELSE
                                           vxc(ic) = vxc(ic) + rr*rc*bascof_lo(pp,mp,ikvecp, ilop,iatom) *CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom))
                                        ENDIF
                                     END DO ! pp
                                  END DO ! mp

                               END DO !ikvecp
                            END DO ! ilop

                            ! calculate matrix-elements of one local orbital with itself

                            DO ikvecp = 1,ikvec
                               ic = ic + 1

                               lp   = l1
                               ilop = ilo
                               DO mp = -lp,lp
                                  DO pp = 1,3
                                     IF ( pp .EQ. 3) THEN
                                        pp2 = pointer_lo(ilop,itype)
                                     ELSE
                                        pp2 = pp
                                     END IF

                                     rr = 0
                                     DO ilharm = 0,nlharm
                                        lh = sphhar%llh(ilharm,typsym)
                                        DO i = 1,sphhar%nmem(ilharm,typsym)
                                           mh = sphhar%mlh(i,ilharm,typsym)
                                           rc = sphhar%clnu(i,ilharm,typsym)* gaunt1(l1,lh,lp,m1,mh,mp,atoms%lmaxd)
                                           rr = rr+integ(ilharm,pp2,lp,pp1,l1)*rc
                                        END DO
                                     END DO

                                     rc = CMPLX(0d0,1d0)**(lp-l1) ! adjusts to a/b/ccof-scaling

                                     IF (hamovlp%l_real) THEN
                                        vxc(ic) = vxc(ic) + invsfct*REAL( rr*rc* bascof_lo(pp,mp,ikvecp,ilop,iatom) * CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom)) )
                                     ELSE
                                        vxc(ic) = vxc(ic) + rr*rc*bascof_lo(pp,mp,ikvecp,ilop, iatom) * CONJG(bascof_lo(p1,m1,ikvec,ilo, iatom))
                                     ENDIF
                                  END DO ! pp
                               END DO ! mp

                            END DO ! ikvecp


                         END DO  ! p1
                      END DO  ! m1

                      icentry = ic       
                   END DO !ikvec
                   ikvecat = ikvecat + invsfct*(2*l1+1)
                END DO  ! ilo
                ikvecprevat = ikvecprevat + ikvecat
                ikvecat     = 0
             END IF  ! atoms%invsat(iatom)

          END DO ! ieq
       END DO !itype

    END IF ! if any atoms%llo

    !initialize weighting factor
    IF( xcpot%icorr .EQ. icorr_hf) THEN
       a_ex = amix_hf
    ELSE  IF( xcpot%icorr .EQ. icorr_pbe0 ) THEN
       a_ex = amix_pbe0
    ELSE IF ( xcpot%icorr .EQ. icorr_hse ) THEN
       a_ex = aMix_HSE
    ELSE IF ( xcpot%icorr .EQ. icorr_vhse ) THEN
       a_ex = aMix_VHSE()
    ELSE
       STOP 'subvxc: error icorr'
    END IF

    IF (hamovlp%l_real) THEN
       DO i=1,hamovlp%matsize
          hamovlp%a_r(i) = hamovlp%a_r(i) - a_ex*REAL(vxc(i))
       ENDDO
    ELSE
       DO i=1,hamovlp%matsize
          hamovlp%a_c(i) = hamovlp%a_c(i) - a_ex*vxc(i)
       ENDDO
    ENDIF

    CALL timestop("subvxc")

    DEALLOCATE( bascof )

  END SUBROUTINE subvxc

END MODULE m_subvxc
