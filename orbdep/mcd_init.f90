MODULE m_mcdinit
CONTAINS
  SUBROUTINE mcd_init(atoms,input,DIMENSION,vr,g,f,mcd,itype,jspin)

    !-----------------------------------------------------------------------
    !
    ! For a given atom-type 'itype' look, whether a core state is in the
    ! energy interval [emcd_lo,emcd_up] and, if found, calculate the 
    ! MCD-matrix elements 'm_mcd'.
    !          
    !-----------------------------------------------------------------------

    USE m_nabla
    USE m_dr2fdr
    USE m_constants, ONLY : c_light
    !USE m_setcor
    USE m_differ
    USE m_types
    IMPLICIT NONE

    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_mcd),INTENT(INOUT)    :: mcd

    INTEGER, PARAMETER :: l_max = 3

    ! Arguments ...

    INTEGER, INTENT (IN)  :: itype
    INTEGER, INTENT (IN)  :: jspin
    REAL,    INTENT (IN)  :: vr(atoms%jmtd,atoms%ntype,input%jspins)
    REAL,    INTENT (IN)  :: f(atoms%jmtd,2,0:atoms%lmaxd,jspin:jspin)
    REAL,    INTENT (IN)  :: g(atoms%jmtd,2,0:atoms%lmaxd,jspin:jspin)

    ! Locals ...

    INTEGER kap,mue,iri,l,ispin,i,icore,korb,nst,n_core,ierr
    REAL  c,t2,e,fj,fl,fn ,d,ms,rn 
    INTEGER kappa(maxval(atoms%econf%num_states)),nprnc(maxval(atoms%econf%num_states)),l_core(maxval(atoms%econf%num_states))
    REAL vrd(atoms%msh),occ(maxval(atoms%econf%num_states),2),a(atoms%msh),b(atoms%msh),j_core(maxval(atoms%econf%num_states)),e_mcd1(maxval(atoms%econf%num_states))
    REAL gv1(atoms%jmtd)
    REAL, ALLOCATABLE :: gc(:,:,:),fc(:,:,:)
    REAL, ALLOCATABLE :: gv(:,:,:,:),fv(:,:,:,:),dgv(:,:,:,:)

    !-----------------------------------------------------------------------

    c = c_light(1.0)
    ALLOCATE ( gc(atoms%jri(itype),atoms%econf(itype)%num_core_states,input%jspins) )
    ALLOCATE ( fc(atoms%jri(itype),atoms%econf(itype)%num_core_states,input%jspins) )

    ! core setup

    mcd%ncore(itype) = 0
    CALL atoms%econf(itype)%get_core(nst,nprnc,kappa,occ)

    DO ispin = jspin, jspin

       ! extend core potential

       DO iri = 1, atoms%jri(itype)
          vrd(iri) = vr(iri,itype,ispin)
       ENDDO
       t2 = vrd(atoms%jri(itype)) / (atoms%jri(itype) - atoms%msh)
       DO iri = atoms%jri(itype) + 1, atoms%msh
          vrd(iri) =  vrd(atoms%jri(itype))  + t2* ( iri-atoms%jri(itype) )
       ENDDO

       ! calculate core

       n_core = 0
       DO korb = 1, atoms%econf(itype)%num_core_states
          IF (occ(korb,1).GT.0) THEN
             fn = nprnc(korb)
             fj = iabs(kappa(korb)) - .5e0
             fl = fj + (.5e0)*isign(1,kappa(korb))
             e = -2* (atoms%zatom(itype)/ (fn+fl))**2
             d = EXP(atoms%dx(itype))
             rn = atoms%rmsh(1,itype)*( d**(atoms%msh-1) )
             CALL differ(fn,fl,fj,c,atoms%zatom(itype),atoms%dx(itype),atoms%rmsh(1,itype),&
                  rn,d,atoms%msh,vrd, e, a,b,ierr)
             IF (ierr/=0)  CALL juDFT_error("error in core-levels", calledby="mcd_init")
             IF ( (e.LE.mcd%emcd_up).AND.(e.GE.mcd%emcd_lo) ) THEN
                WRITE(*,*) 'good    ev = ',e
                n_core = n_core + 1
                j_core(n_core) = fj
                l_core(n_core) = NINT( fl )
                e_mcd1(n_core) = e
                DO iri = 1, atoms%jri(itype)
                   gc(iri,n_core,ispin) = a(iri)
                   fc(iri,n_core,ispin) = b(iri)
                ENDDO
             ENDIF
          ENDIF
       ENDDO

    ENDDO

    !-----------------------------------------------------------------------

    IF (n_core.GT.0) THEN

       ALLOCATE ( gv(atoms%jri(itype),0:l_max,input%jspins,2) )
       ALLOCATE (dgv(atoms%jri(itype),0:l_max,input%jspins,2) )
       ALLOCATE ( fv(atoms%jri(itype),0:l_max,input%jspins,2) )
       DO i = 1, 2
          DO iri = 3*(itype-1)+1 , 3*(itype-1)+3
             DO l = 1, (l_max+1)**2
                DO icore = 1, maxval(atoms%econf%num_states)
                   mcd%m_mcd(icore,l,iri,i) = CMPLX(0.0,0.0)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
       ! bring LAPW wavefunctions in a proper form:
       !
       DO ispin = jspin, jspin
          ms = ispin - 1.5
          DO l = 0, l_max
             DO iri = 1, atoms%jri(itype)
                gv(iri,l,ispin,1) = f(iri,1,l,ispin)   ! large component of u
                fv(iri,l,ispin,1) = f(iri,2,l,ispin)   ! small              .
                gv(iri,l,ispin,2) = g(iri,1,l,ispin)   ! large component of u
                fv(iri,l,ispin,2) = g(iri,2,l,ispin)   ! small
             ENDDO
             gv1(:) = atoms%rmsh(:,itype) * gv(:,l,ispin,1)
             CALL dr2fdr(&                                          ! deriative of u (large)&
                  gv1,atoms%rmsh(1,itype),atoms%jri(itype), dgv(1,l,ispin,1) )
             gv1(:) = atoms%rmsh(:,itype) * gv(:,l,ispin,2)              !              .
             CALL dr2fdr(&                                          ! deriative of u (large)&
                  gv1,atoms%rmsh(1,itype),atoms%jri(itype), dgv(1,l,ispin,2) )
          ENDDO
          !
          !
          !
          DO icore = 1, n_core

             DO i = 1, 2
                !              write(*,*) j_core(icore),l_core(icore),l_max,ms
                CALL nabla(itype,icore,atoms%jri(itype),atoms%dx(itype),maxval(atoms%econf%num_states),atoms%ntype,&
                     j_core(icore),l_core(icore),l_max,ms,atoms%rmsh(:,itype),gc(:,icore,ispin),&
                     gv(:,0:,ispin,i),dgv(:,0:,ispin,i), mcd%m_mcd(:,:,:,i) )
             ENDDO

             DO i = 1, 2*icore*l_core(icore)
                mcd%ncore(itype) = mcd%ncore(itype) + 1
                IF (mcd%ncore(itype)>maxval(atoms%econf%num_states))  CALL juDFT_error("maxval(atoms%econf%num_states) too small" ,calledby ="mcd_init")
                mcd%e_mcd(itype,ispin,mcd%ncore(itype)) = e_mcd1(icore)
             ENDDO
          ENDDO
       ENDDO

       DEALLOCATE (gv,fv,dgv)
    ENDIF
    DEALLOCATE (gc,fc)


    !      DO i = 1, 2
    !       DO iri = 3*(itype-1)+1 , 3*(itype-1)+3
    !         write (*,*) iri
    !         DO icore = 1, mcd%ncore(itype)
    !           write (*,'(10f10.5)') (mcd%m_mcd(icore,l,iri,i),l=1,9)
    !         ENDDO
    !       ENDDO
    !      ENDDO
  END SUBROUTINE mcd_init
END MODULE m_mcdinit
