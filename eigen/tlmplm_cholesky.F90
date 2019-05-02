MODULE m_tlmplm_cholesky
  use m_judft
  IMPLICIT NONE
  !*********************************************************************
  !     sets up the local Hamiltonian, i.e. the Hamiltonian in the
  !     l',m',l,m,u- basis which is independent from k!
  !     shifts this local Hamiltonian to make it positive definite
  !     and does a cholesky decomposition
  !*********************************************************************
CONTAINS
  SUBROUTINE tlmplm_cholesky(sphhar,atoms,noco,enpara,&
       jspin,jsp,mpi,v,input,td,ud,hub1)
    USE m_tlmplm
    USE m_types
    USE m_gaunt, ONLY: gaunt2
    USE m_constants, ONLY: lmaxU_const
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_hub1ham),INTENT(IN)  :: hub1
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin,jsp !physical spin&spin index for data
    INTEGER, INTENT (IN) :: iterHIA
    !     ..
    TYPE(t_potden),INTENT(IN)    :: v
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud

    !     ..
    !     .. Local Scalars ..
    REAL temp
    INTEGER i,l,lm,lmin,lmin0,lmp,lmplm,lp,info,in
    INTEGER lpl ,mp,n,m,s,i_u
    LOGICAL OK
    LOGICAL l_hia(0:lmaxU_const)
    !     ..
    !     .. Local Arrays ..


    REAL,PARAMETER:: e_shift_min=0.5
    REAL,PARAMETER:: e_shift_max=65.0


    !     ..e_shift
    td%e_shift(:,jsp)=0.0
    IF (jsp<3) td%e_shift(:,jsp)=e_shift_min


    td%tdulo(:,:,:,jsp) = CMPLX(0.0,0.0)
    td%tuulo(:,:,:,jsp) = CMPLX(0.0,0.0)
    td%tuloulo(:,:,:,jsp) = CMPLX(0.0,0.0)



    td%h_off=0.0
    !$    call gaunt2(atoms%lmaxd)
    !$OMP PARALLEL DO DEFAULT(NONE)&
    !$OMP PRIVATE(temp,i,l,lm,lmin,lmin0,lmp)&
    !$OMP PRIVATE(lmplm,lp,m,mp,n)&
    !$OMP PRIVATE(OK,s,in,info)&
    !$OMP SHARED(atoms,jspin,jsp,sphhar,enpara,td,ud,v,mpi,input,hub1)
    DO  n = 1,atoms%ntype
       l_hia = .FALSE.
       DO i = 1, hub1%n_hia
          IF(n.EQ.hub1%lda_u(i)%atomType) THEN
             l_hia(hub1%lda_u(i)%l) = .TRUE..AND.(iterHIA.GT.0)
          ENDIF
       ENDDO
       CALL tlmplm(n,sphhar,atoms,enpara,jspin,jsp,mpi,v,input,td,ud,l_hia)
       OK=.FALSE.
       cholesky_loop:DO WHILE(.NOT.OK)
          td%h_loc(:,:,n,jsp)=0.0
          OK=.TRUE.
          !
          !--->    generate the wavefunctions for each l
          !
          s=atoms%lnonsph(n)*(atoms%lnonsph(n)+2)+1
          !Setup local hamiltonian
          DO lmp=0,atoms%lnonsph(n)*(atoms%lnonsph(n)+2)
             lp=FLOOR(SQRT(1.0*lmp))
             mp=lmp-lp*(lp+1)
             IF (lp>atoms%lmax(n).OR.ABS(mp)>lp) STOP "BUG"
             !--->             loop over l,m
             DO l = 0,atoms%lnonsph(n)
                DO m = -l,l
                   lm = l* (l+1) + m
                   in = td%ind(lmp,lm,n,jsp)
                   IF (in/=-9999) THEN
                      IF (in>=0) THEN
                         td%h_loc(lm,lmp,n,jsp)    = CONJG(td%tuu(in,n,jsp))
                         td%h_loc(lm+s,lmp,n,jsp)  = CONJG(td%tud(in,n,jsp))
                         td%h_loc(lm,lmp+s,n,jsp)  = CONJG(td%tdu(in,n,jsp))
                         td%h_loc(lm+s,lmp+s,n,jsp)= CONJG(td%tdd(in,n,jsp))
                      ELSE
                         td%h_loc(lm,lmp,n,jsp)    = td%tuu(-in,n,jsp)
                         td%h_loc(lm+s,lmp,n,jsp)  = td%tdu(-in,n,jsp)
                         td%h_loc(lm,lmp+s,n,jsp)  = td%tud(-in,n,jsp)
                         td%h_loc(lm+s,lmp+s,n,jsp)= td%tdd(-in,n,jsp)
                      END IF
                   END IF
                END DO
             END DO
          ENDDO
          !Include contribution from LDA+U
          DO i_u=1,atoms%n_u
             IF (n.NE.atoms%lda_u(i_u)%atomtype) CYCLE
             !Found a "U" for this atom type
             l=atoms%lda_u(i_u)%l
             lp=atoms%lda_u(i_u)%l
             DO m = -l,l
                lm = l* (l+1) + m
                DO mp = -lp,lp
                   lmp = lp* (lp+1) + mp
                   td%h_loc(lm,lmp,n,jsp)     =td%h_loc(lm,lmp,n,jsp) + v%mmpMat(m,mp,i_u,jsp)
                   td%h_loc(lm+s,lmp+s,n,jsp) =td%h_loc(lm+s,lmp+s,n,jsp)+ v%mmpMat(m,mp,i_u,jsp)*ud%ddn(lp,n,jsp)
                ENDDO
             ENDDO
          END DO

         !Contribution from LDA+HIA
         DO i_u=1,hub1%n_hia
            IF (n.NE.hub1%lda_u(i_u)%atomtype) CYCLE
            !Found a "U" for this atom type
            l=hub1%lda_u(i_u)%l
            lp=hub1%lda_u(i_u)%l
            DO m = -l,l
               lm = l* (l+1) + m
               DO mp = -lp,lp
                  lmp = lp* (lp+1) + mp
                  td%h_loc(lm,lmp,n,jsp)     =td%h_loc(lm,lmp,n,jsp) + v%mmpMat(m,mp,atoms%n_u+i_u,jsp)
                  td%h_loc(lm+s,lmp+s,n,jsp) =td%h_loc(lm+s,lmp+s,n,jsp)+ v%mmpMat(m,mp,atoms%n_u+i_u,jsp)*ud%ddn(lp,n,jsp)
               ENDDO
            ENDDO
         END DO

          !Now add diagonal contribution to matrices
          IF (jsp<3) THEN
             DO l = 0,atoms%lmax(n)
                DO  m = -l,l
                   lm = l* (l+1) + m
                   lmplm = (lm* (lm+3))/2
                   td%tuu(lmplm,n,jsp)=td%tuu(lmplm,n,jsp) + enpara%el0(l,n,jsp)
                   td%tdd(lmplm,n,jsp)=td%tdd(lmplm,n,jsp) + enpara%el0(l,n,jsp)*ud%ddn(l,n,jsp)
                   td%tud(lmplm,n,jsp)=td%tud(lmplm,n,jsp) + 0.5
                   td%tdu(lmplm,n,jsp)=td%tdu(lmplm,n,jsp) + 0.5
                ENDDO
             ENDDO
             !Create Cholesky decomposition of local hamiltonian
             
             !--->    Add diagonal terms to make matrix positive definite
             DO lp = 0,atoms%lnonsph(n)
                DO mp = -lp,lp
                   lmp = lp* (lp+1) + mp
                   td%h_loc(lmp,lmp,n,jsp)=td%e_shift(n,jsp)+td%h_loc(lmp,lmp,n,jsp)
                   td%h_loc(lmp+s,lmp+s,n,jsp)=td%e_shift(n,jsp)*ud%ddn(lp,n,jsp)+td%h_loc(lmp+s,lmp+s,n,jsp)
                END DO
             END DO
             IF (lmp+1.NE.s) CALL judft_error("BUG in tlmpln_cholesky")
             !Perform cholesky decomposition
             info=0
             CALL zpotrf("L",2*s,td%h_loc(:,:,n,jsp),SIZE(td%h_loc,1),info)

             !Upper part to zero
             DO l=0,2*s-1
                DO lp=0,l-1
                   td%h_loc(lp,l,n,jsp)=0.0
                ENDDO
             ENDDO

             IF (info.NE.0) THEN
                td%e_shift(n,jsp)=td%e_shift(n,jsp)*2.0
                PRINT *,"Potential shift to small, increasing the value to:",td%e_shift(n,jsp)
                IF (td%e_shift(n,jsp)>e_shift_max) THEN
                   CALL judft_error("Potential shift at maximum")
                ENDIF
                OK=.FALSE.
             ENDIF
          ENDIF
       ENDDO cholesky_loop
    ENDDO
    !$OMP END PARALLEL DO
    IF (noco%l_constr) CALL tlmplm_constrained(atoms,v,enpara,input,ud,noco,td)



  END SUBROUTINE tlmplm_cholesky


   


  SUBROUTINE tlmplm_constrained(atoms,v,enpara,input,ud,noco,td)
    USE m_radovlp
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_enpara),INTENT(IN)   :: enpara
    TYPE(t_potden),INTENT(IN)   :: v
    TYPE(t_tlmplm),INTENT(INOUT):: td
    TYPE(t_usdus),INTENT(INOUT) :: ud
    TYPE(t_noco),INTENT(IN)     :: noco
    
    REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)
    COMPLEX :: c
    INTEGER :: n,l,s  
      
    
    ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntype),udn21(0:atoms%lmaxd,atoms%ntype),&
         dun21(0:atoms%lmaxd,atoms%ntype),ddn21(0:atoms%lmaxd,atoms%ntype) )
    CALL rad_ovlp(atoms,ud,input,v%mt,enpara%el0, uun21,udn21,dun21,ddn21)
    
    DO  n = 1,atoms%ntype
       !If we do  a constraint calculation, we have to calculate the
       !local spin off-diagonal contributions
       s=atoms%lnonsph(n)+1
       !first ispin=2,jspin=1 case
       DO l=0,atoms%lnonsph(n)
          c=(-0.5)*CMPLX(noco%b_con(1,n),noco%b_con(2,n))
          td%h_off(l  ,l  ,n,1)     =td%h_off(l  ,l  ,n,1) + uun21(l,n)*c
          td%h_loc(l  ,l+s,n,1)     =td%h_off(l  ,l+s,n,1) + udn21(l,n)*c
          td%h_loc(l+s,l  ,n,1)     =td%h_off(l+s,l  ,n,1) + dun21(l,n)*c
          td%h_loc(l+s,l+s,n,1)     =td%h_off(l+s,l+s,n,1) + ddn21(l,n)*c
       ENDDO
       
       
       !then ispin=2,jspin=1 case
       DO l=0,atoms%lnonsph(n)
          c=(-0.5)*CMPLX(noco%b_con(1,n),-noco%b_con(2,n))
          td%h_off(l  ,l  ,n,2)     =td%h_off(l  ,l  ,n,2) + uun21(l,n)*c
          td%h_loc(l  ,l+s,n,2)     =td%h_off(l  ,l+s,n,2) + udn21(l,n)*c
          td%h_loc(l+s,l  ,n,2)     =td%h_off(l+s,l  ,n,2) + dun21(l,n)*c
          td%h_loc(l+s,l+s,n,2)     =td%h_off(l+s,l+s,n,2) + ddn21(l,n)*c
       ENDDO
    END DO
  END SUBROUTINE tlmplm_constrained
  
  END MODULE m_tlmplm_cholesky
  
