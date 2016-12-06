MODULE m_forcea21
CONTAINS
  SUBROUTINE force_a21(&
       input,atoms,DIMENSION,nobd,sym,oneD,cell,&
       we,jsp,epar,ne,eig,usdus,&
       acof,bcof,ccof,aveccof,bveccof,cveccof, results,f_a21,f_b4)

    ! ************************************************************
    ! Pulay 2nd and 3rd (A17+A20) term force contribution a la Rici
    ! combined
    ! NOTE: we do NOT include anymore  the i**l factors
    ! in the alm,blm coming from to_pulay. Therefore, we can
    ! use matrixelements from file 28,38 DIRECTLY
    ! note: present version only yields forces for
    ! highest energy window (=valence states)
    ! if also semicore forces are wanted the tmas and tmat files
    ! have to be saved, indexed and properly used here in force_a21
    ! 22/june/97: probably we found symmetrization error replacing
    ! now S^-1 by S (IS instead of isinv)
    ! ************************************************************
    !
    ! Force contribution B4 added following
    ! Madsen, Blaha, Schwarz, Sjostedt, Nordstrom
    ! GMadsen FZJ 20/3-01
    !
    USE m_forcea21lo
    USE m_forcea21U
    USE m_tlmplm_store
    USE m_types
    USE m_constants
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_usdus),INTENT(IN)     :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nobd
    INTEGER, INTENT (IN) :: ne,jsp
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(nobd),epar(0:atoms%lmaxd,atoms%ntypd)
    REAL,    INTENT (IN) :: eig(DIMENSION%neigd)  
    COMPLEX, INTENT (INOUT) :: f_a21(3,atoms%ntypd),f_b4(3,atoms%ntypd)
    COMPLEX, INTENT (IN) :: acof(nobd,0:atoms%lmaxd*(atoms%lmaxd+2) ,atoms%natd)
    COMPLEX, INTENT (IN) :: bcof(nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd )
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
    COMPLEX, INTENT (IN) :: aveccof(3,nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd )
    COMPLEX, INTENT (IN) :: bveccof(3,nobd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%natd )
    COMPLEX, INTENT (IN) :: cveccof(3,-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%natd)
    !     ..
    !     .. Local Scalars ..
    INTEGER, PARAMETER :: lmaxb=3
    COMPLEX dtd,dtu,utd,utu
    INTEGER lo, mlotot, mlolotot, mlot_d, mlolot_d
    INTEGER i,ie,im,in,l1,l2,ll1,ll2,lm1,lm2,m1,m2,n,natom,m
    INTEGER natrun,is,isinv,j,irinv,it
    REAL   ,PARAMETER:: zero=0.0
    COMPLEX,PARAMETER:: czero=CMPLX(0.,0.)
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: v_mmp(:,:)
    REAL,    ALLOCATABLE :: a21(:,:),b4(:,:)
    COMPLEX forc_a21(3),forc_b4(3)
    REAL starsum(3),starsum2(3),gvint(3),gvint2(3)
    REAL vec(3),vec2(3),vecsum(3),vecsum2(3)

    TYPE(t_tlmplm)::tlmplm
    !     ..
    !     ..
    !dimension%lmplmd = (dimension%lmd* (dimension%lmd+3))/2
    mlotot = 0 ; mlolotot = 0
    DO n = 1, atoms%ntype
       mlotot = mlotot + atoms%nlo(n)
       mlolotot = mlolotot + atoms%nlo(n)*(atoms%nlo(n)+1)/2
    ENDDO
    mlot_d = MAX(mlotot,1)
    mlolot_d = MAX(mlolotot,1)
    ALLOCATE ( tlmplm%tdd(0:DIMENSION%lmplmd,atoms%ntype,1),tlmplm%tuu(0:DIMENSION%lmplmd,atoms%ntype,1),&
         tlmplm%tdu(0:DIMENSION%lmplmd,atoms%ntype,1),tlmplm%tud(0:DIMENSION%lmplmd,atoms%ntype,1),&
         tlmplm%tuulo(0:DIMENSION%lmd,-atoms%llod:atoms%llod,mlot_d,1),&
         tlmplm%tdulo(0:DIMENSION%lmd,-atoms%llod:atoms%llod,mlot_d,1),&
         tlmplm%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,mlolot_d,1),&
         v_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb),&
         a21(3,atoms%natd),b4(3,atoms%natd),tlmplm%ind(0:DIMENSION%lmd,0:DIMENSION%lmd,atoms%ntype,1) )
    !
    natom = 1
    DO  n = 1,atoms%ntype
       IF (atoms%l_geo(n)) THEN
          forc_a21(:) = czero
          forc_b4(:) = czero


          CALL read_tlmplm(n,jsp,atoms%nlo,atoms%lda_u%l.GE.0,&
               tlmplm%tuu(:,n,1),tlmplm%tud(:,n,1),tlmplm%tdu(:,n,1),tlmplm%tdd(:,n,1),&
               tlmplm%ind(:,:,n,1),tlmplm%tuulo(:,:,:,1),tlmplm%tuloulo(:,:,:,1),tlmplm%tdulo(:,:,:,1),v_mmp)

          DO natrun = natom,natom + atoms%neq(n) - 1
             a21(:,natrun) = zero
             b4(:,natrun) = zero

          END DO
          !
          DO ie = 1,ne
             !
             !
             DO l1 = 0,atoms%lmax(n)
                ll1 = l1* (l1+1)
                DO m1 = -l1,l1
                   lm1 = ll1 + m1
                   DO l2 = 0,atoms%lmax(n)
                      !
                      ll2 = l2* (l2+1)
                      DO m2 = -l2,l2
                         lm2 = ll2 + m2
                         DO natrun = natom,natom + atoms%neq(n) - 1
                            in = tlmplm%ind(lm1,lm2,n,1)
                            IF (in.NE.-9999) THEN
                               IF (in.GE.0) THEN
                                  !
                                  ! ATTENTION: the matrix elements tuu,tdu,tud,tdd
                                  ! as calculated in tlmplm are the COMPLEX CONJUGATE
                                  ! of the non-spherical matrix elements because in the
                                  ! matrix building routine hssphn (or similar routines)
                                  ! the COMPLEX CONJUGATE of alm,blm is calculated (to
                                  ! save complex operations presumably)
                                  ! Her, A20 is formulated in the usual way therefore
                                  ! we have to take the COMPLEX CONJUGATE versions
                                  ! of tuu,tdu,tud,tdd as compared to hssphn!
                                  !
                                  utu = tlmplm%tuu(in,n,1)
                                  dtu = tlmplm%tdu(in,n,1)
                                  utd = tlmplm%tud(in,n,1)
                                  dtd = tlmplm%tdd(in,n,1)
                               ELSE
                                  im = -in
                                  utu = CONJG(tlmplm%tuu(im,n,1))
                                  dtd = CONJG(tlmplm%tdd(im,n,1))
                                  utd = CONJG(tlmplm%tdu(im,n,1))
                                  dtu = CONJG(tlmplm%tud(im,n,1))
                               END IF
                               DO i = 1,3
                                  a21(i,natrun) = a21(i,natrun) + 2.0*&
                                       AIMAG( CONJG(acof(ie,lm1,natrun)) *utu*aveccof(i,ie,lm2,natrun)&
                                       +CONJG(acof(ie,lm1,natrun)) *utd*bveccof(i,ie,lm2,natrun)&
                                       +CONJG(bcof(ie,lm1,natrun)) *dtu*aveccof(i,ie,lm2,natrun)&
                                       +CONJG(bcof(ie,lm1,natrun)) *dtd*bveccof(i,ie,lm2,natrun))*we(ie)/atoms%neq(n)
                                  !   END i loop
                               END DO
                            END IF
                            !   END natrun
                         END DO
                         !
                         !   END m2 loop
                      END DO
                      !   END l2 loop
                   END DO
                   !+gu 20.11.97
                   utu = epar(l1,n)-eig(ie)
                   utd = 0.5
                   dtu = 0.5
                   dtd = utu*usdus%ddn(l1,n,jsp)
                   DO i = 1,3
                      DO natrun = natom,natom + atoms%neq(n) - 1
                         a21(i,natrun) = a21(i,natrun) + 2.0*&
                              AIMAG(CONJG(acof(ie,lm1,natrun)) *utu*aveccof(i,ie,lm1,natrun)&
                              +CONJG(acof(ie,lm1,natrun)) *utd*bveccof(i,ie,lm1,natrun)&
                              +CONJG(bcof(ie,lm1,natrun)) *dtu*aveccof(i,ie,lm1,natrun)&
                              +CONJG(bcof(ie,lm1,natrun)) *dtd*bveccof(i,ie,lm1,natrun)&
                              )*we(ie) /atoms%neq(n)
                      END DO
                      !
                      !-gu
                      ! END  i loop
                   END DO
                   !   END m1 loop
                END DO
                !   END l1 loop
             END DO
             !   END ie loop
          END DO
          !
          !--->    add the local orbital and U contribution to a21
          !
          CALL force_a21_lo(nobd,atoms,jsp,n,we,eig,ne,&
               acof,bcof,ccof,aveccof,bveccof,&
               cveccof, tlmplm,usdus, a21)

          CALL force_a21_U(nobd,atoms,lmaxb,n,jsp,we,ne,&
               usdus,v_mmp,acof,bcof,ccof,&
               aveccof,bveccof,cveccof, a21)
          IF (input%l_useapw) THEN
             ! -> B4 force
             DO ie = 1,ne
                DO l1 = 0,atoms%lmax(n)
                   ll1 = l1* (l1+1)
                   DO m1 = -l1,l1
                      lm1 = ll1 + m1
                      DO i = 1,3
                         DO natrun = natom,natom + atoms%neq(n) - 1
                            b4(i,natrun) = b4(i,natrun) + 0.5 *&
                                 we(ie)/atoms%neq(n)*atoms%rmt(n)**2*AIMAG(&
                                 CONJG(acof(ie,lm1,natrun)*usdus%us(l1,n,jsp)&
                                 +bcof(ie,lm1,natrun)*usdus%uds(l1,n,jsp))*&
                                 (aveccof(i,ie,lm1,natrun)*usdus%dus(l1,n,jsp)&
                                 +bveccof(i,ie,lm1,natrun)*usdus%duds(l1,n,jsp) )&
                                 -CONJG(aveccof(i,ie,lm1,natrun)*usdus%us(l1,n,jsp)&
                                 +bveccof(i,ie,lm1,natrun)*usdus%uds(l1,n,jsp) )*&
                                 (acof(ie,lm1,natrun)*usdus%dus(l1,n,jsp)&
                                 +bcof(ie,lm1,natrun)*usdus%duds(l1,n,jsp)) )
                         END DO
                      END DO
                   END DO
                END DO
                DO lo = 1,atoms%nlo(n)
                   l1 = atoms%llo(lo,n)
                   DO m = -l1,l1
                      lm1 = l1* (l1+1) + m
                      DO i=1,3
                         DO natrun = natom,natom + atoms%neq(n) - 1
                            b4(i,natrun) = b4(i,natrun) + 0.5 *&
                                 we(ie)/atoms%neq(n)*atoms%rmt(n)**2*AIMAG(&
                                 CONJG( acof(ie,lm1,natrun)* usdus%us(l1,n,jsp)&
                                 + bcof(ie,lm1,natrun)* usdus%uds(l1,n,jsp) ) *&
                                 cveccof(i,m,ie,lo,natrun)*usdus%dulos(lo,n,jsp)&
                                 + CONJG(ccof(m,ie,lo,natrun)*usdus%ulos(lo,n,jsp)) *&
                                 ( aveccof(i,ie,lm1,natrun)* usdus%dus(l1,n,jsp)&
                                 + bveccof(i,ie,lm1,natrun)* usdus%duds(l1,n,jsp)&
                                 + cveccof(i,m,ie,lo,natrun)*usdus%dulos(lo,n,jsp) )  &
                                 - (CONJG( aveccof(i,ie,lm1,natrun) *usdus%us(l1,n,jsp)&
                                 + bveccof(i,ie,lm1,natrun) *usdus%uds(l1,n,jsp) ) *&
                                 ccof(m,ie,lo,natrun)  *usdus%dulos(lo,n,jsp)&
                                 + CONJG(cveccof(i,m,ie,lo,natrun)*usdus%ulos(lo,n,jsp)) *&
                                 ( acof(ie,lm1,natrun)*usdus%dus(l1,n,jsp)&
                                 + bcof(ie,lm1,natrun)*usdus%duds(l1,n,jsp)&
                                 + ccof(m,ie,lo,natrun)*usdus%dulos(lo,n,jsp) ) ) )  
                         END DO
                      ENDDO
                   ENDDO
                ENDDO
             END DO
          ENDIF
          !
          DO natrun = natom,natom + atoms%neq(n) - 1
             !
             !  to complete summation over stars of k now sum
             !  over all operations which leave (k+G)*R(natrun)*taual(natrun)
             !  invariant. We sum over ALL these operations and not only
             !  the ones needed for the actual star of k. Should be
             !  ok if we divide properly by the number of operations
             !  First, we find operation S where RS=T. T -like R- leaves
             !  the above scalar product invariant (if S=1 then R=T).
             !  R is the operation which generates position of equivalent atom
             !  out of position of representative
             !  S=R^(-1) T
             !  number of ops which leave (k+G)*op*taual invariant: invarind
             !  index of inverse operation of R: irinv
             !  index of operation T: invarop
             !  now, we calculate index of operation S: is
             !
             !  note, that vector in expression A17,A20 + A21 is a
             !  reciprocal lattice vector! other transformation rules
             !
             !  transform recip vector g-g' into internal coordinates

             vec(:) = a21(:,natrun)
             vec2(:) = b4(:,natrun)

             gvint=MATMUL(cell%bmat,vec)/tpi_const
             gvint2=MATMUL(cell%bmat,vec2)/tpi_const
             vecsum(:) = zero
             vecsum2(:) = zero

             !-gb2002
             !            irinv = invtab(ngopr(natrun))
             !            DO it = 1,invarind(natrun)
             !               is = multab(irinv,invarop(natrun,it))
             !c  note, actually we need the inverse of S but -in principle
             !c  because {S} is a group and we sum over all S- S should also
             !c  work; to be lucid we take the inverse:
             !                isinv = invtab(is)
             !!               isinv = is
             ! Rotation is alreadt done in to_pulay, here we work only in the
             ! coordinate system of the representative atom (natom):
             !!        
             DO it = 1,sym%invarind(natom)
                is =sym%invarop(natom,it)
                isinv = sym%invtab(is)
                IF (oneD%odi%d1) isinv = oneD%ods%ngopr(natom)
                !-gb 2002
                !  now we have the wanted index of operation with which we have
                !  to rotate gv. Note gv is given in cart. coordinates but
                !  mrot acts on internal ones
                DO i = 1,3
                   vec(i) = zero
                   vec2(i) = zero
                   DO j = 1,3
                      IF (.NOT.oneD%odi%d1) THEN
                         vec(i) = vec(i) + sym%mrot(i,j,isinv)*gvint(j)
                         vec2(i) = vec2(i) + sym%mrot(i,j,isinv)*gvint2(j)
                      ELSE
                         vec(i) = vec(i) + oneD%ods%mrot(i,j,isinv)*gvint(j)
                         vec2(i) = vec2(i) + oneD%ods%mrot(i,j,isinv)*gvint2(j)
                      END IF
                   END DO
                END DO
                DO i = 1,3
                   vecsum(i) = vecsum(i) + vec(i)
                   vecsum2(i) = vecsum2(i) + vec2(i)
                END DO
                !   end operator loop
             END DO
             !
             !   transform from internal to cart. coordinates
             starsum=MATMUL(cell%amat,vecsum)
             starsum2=MATMUL(cell%amat,vecsum2)
             DO i = 1,3
                forc_a21(i) = forc_a21(i) + starsum(i)/sym%invarind(natrun)
                forc_b4(i) = forc_b4(i) + starsum2(i)/sym%invarind(natrun)
             END DO
             !
             !  natrun loop end
          END DO
          !
          !     sum to existing forces
          !
          !  NOTE: force() IS REAL AND THEREFORE TAKES ONLY THE
          !  REAL PART OF forc_a21(). IN GENERAL, FORCE MUST BE
          !  REAL AFTER k-STAR SUMMATION. NOW, WE PUT THE PROPER
          !  OPERATIONS INTO REAL SPACE. PROBLEM: WHAT HAPPENS
          !  IF IN REAL SPACE THERE IS NO INVERSION ANY MORE?
          !  BUT WE HAVE INVERSION IN k-SPACE DUE TO TIME REVERSAL
          !  SYMMETRY, E(k)=E(-k)
          !  WE ARGUE THAT k-SPACE INVERSION IS AUTOMATICALLY TAKEN
          !  INTO ACCOUNT IF FORCE = (1/2)(forc_a21+conjg(forc_a21))
          !  BECAUSE TIME REVERSAL SYMMETRY MEANS THAT conjg(PSI)
          !  IS ALSO A SOLUTION OF SCHR. EQU. IF PSI IS ONE.
          DO i = 1,3
             results%force(i,n,jsp) = results%force(i,n,jsp) + REAL(forc_a21(i) + forc_b4(i))
             f_a21(i,n)     = f_a21(i,n)     + forc_a21(i)
             f_b4(i,n)      = f_b4(i,n)      + forc_b4(i)
          END DO
          !
          !     write result moved to force_a8
          !
          !         write(*,*) a21(:,n) 
       ENDIF                                            !  IF (atoms%l_geo(n)) ...
       natom = natom + atoms%neq(n)
    ENDDO
    !
    DEALLOCATE (tlmplm%tdd,tlmplm%tuu,tlmplm%tdu,tlmplm%tud,tlmplm%tuulo,tlmplm%tdulo,tlmplm%tuloulo,v_mmp,tlmplm%ind,a21,b4)

  END SUBROUTINE force_a21
END MODULE m_forcea21
