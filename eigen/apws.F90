MODULE m_apws
  use m_juDFT
  !*********************************************************************
  !     determines the lapw list such that |k+G|<rkmax.
  !     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.
  !        m. weinert  1986
  !     unit 29 removed gb 2004
  !*********************************************************************
  !     modified for explicit use of z-reflection symmetry in seclr4.f
  !        g. bihlmayer '96
  !     subroutine boxdim added to treat non-orthogonal lattice vectors
  !        s.bluegel, IFF, 18.Nov.97
  !*********************************************************************
CONTAINS
  SUBROUTINE apws(dimension,input,noco,kpts,&
       nk,cell,l_zref,n_size,jspin,bkpt,lapw,nred)

    USE m_types
    USE m_sort
    USE m_boxdim
    IMPLICIT NONE

    TYPE(t_dimension),INTENT(IN)   :: dimension
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_lapw),INTENT(INOUT)     :: lapw
    !     .. 
    !     .. Scalar Arguments ..
    INTEGER, INTENT  (IN) :: nk,n_size,jspin
    INTEGER, INTENT (OUT) :: nred
    LOGICAL, INTENT (IN)  :: l_zref
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (OUT) :: bkpt(3)
    !     ..
    !     .. Local Scalars ..
    REAL arltv1,arltv2,arltv3,r2,rk2,rkm,t,r2q,gla,eps
    INTEGER i,itt,j,j1,j2,j3,k,l ,mk1,mk2,mk3,n,ispin,jsp_start,jsp_end,gmi,m
    LOGICAL :: done
    !     ..
    !     .. Local Arrays ..
    REAL s(3),sq(3),rkq(dimension%nvd),gsk3(dimension%nvd)
    INTEGER k1rev(dimension%nvd),k2rev(dimension%nvd),k3rev(dimension%nvd),index3(dimension%nvd)
#ifdef CPP_MPI
    INTEGER              :: n_inner,n_bound
    REAL,    ALLOCATABLE :: rk_help(:)
    INTEGER, ALLOCATABLE :: k_help(:,:) ,pos(:)
#endif
    IF (.not.allocated(lapw%k1)) THEN
       ALLOCATE ( lapw%k1(DIMENSION%nvd,DIMENSION%jspd),lapw%k2(DIMENSION%nvd,DIMENSION%jspd),&
            lapw%k3(DIMENSION%nvd,DIMENSION%jspd),lapw%rk(DIMENSION%nvd,DIMENSION%jspd) )
       ALLOCATE(lapw%gvec(3,DIMENSION%nvd,DIMENSION%jspd))
       ALLOCATE(lapw%vk(3,DIMENSION%nvd,DIMENSION%jspd))
       ALLOCATE(lapw%gk(3,DIMENSION%nvd,DIMENSION%jspd))
    ENDIF
    lapw%rk = 0 ; lapw%k1 = 0 ; lapw%k2 = 0 ; lapw%k3 = 0 ;lapw%nv=0
    !     ..
    !     ..
    !---> in a spin-spiral calculation different basis sets are used for
    !---> the two spin directions, because the cutoff radius is defined
    !---> by |G + k +/- qss/2| < rkmax.
    IF (nk>kpts%nkpt) THEN
       bkpt(:)=kpts%bkf(:,nk)
    ELSE
       bkpt(:) = kpts%bk(:,nk)
    ENDIF
    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !     arltv(i) length of reciprical lattice vector along direction (i)
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int( input%rkmax/arltv1 ) + 4
    mk2 = int( input%rkmax/arltv2 ) + 4
    mk3 = int( input%rkmax/arltv3 ) + 4

    IF (noco%l_ss) THEN
       jsp_start = 1
       jsp_end   = input%jspins
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
    ENDIF
    rkm = input%rkmax
    !---> if too many basis functions, reduce rkmax
    DO ispin = jsp_start,jsp_end
       done=.false.
       rkm_reduction_loop:DO while(.not.done)
          rk2 = rkm*rkm
          !--->    obtain vectors
          n = 0
          DO  j1 = -mk1,mk1
             s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*noco%qss(1)
             sq(1) = bkpt(1) + j1
             DO  j2 = -mk2,mk2
                s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*noco%qss(2)
                sq(2) = bkpt(2) + j2 
                DO  j3 = -mk3,mk3
                   s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*noco%qss(3)
                   sq(3) = bkpt(3) + j3
                   r2 = dot_product(s,matmul(s,cell%bbmat))
                   r2q = dot_product(sq,matmul(sq,cell%bbmat))
                   IF (r2.LE.rk2) THEN
                      n = n + 1
                      IF (n.GT.dimension%nvd) THEN
                         rkm = rkm - 0.1
                         WRITE (*,FMT=8000) (bkpt(i),i=1,3),rkm
8000                     FORMAT (' $$$ k=(',3f10.6,'): rkm truncated to',f12.6)
                         cycle rkm_reduction_loop
                      endif
                      lapw%k1(n,ispin) = j1
                      lapw%k2(n,ispin) = j2
                      lapw%k3(n,ispin) = j3
                      lapw%rk(n,ispin) = sqrt(r2)
                      rkq(n) = sqrt(r2q)
                   END IF
                enddo
             enddo
          enddo
          done=.true.
       enddo rkm_reduction_loop
       lapw%nv(ispin) = n
       !
       !--->    sort by shell-metzner
       !
       ! (for spin-spirals & LO's we have to sort according to the k+G's (rkq), not
       !  the k+G+q's (rk). Otherwise we might couple an LO to k+G1+q and k+G2-q !)
       !                                                                       gb01
       m = n
80     m = m/2
       IF (m.LE.0) GO TO 130
       k = n - m
       j = 1
90     i = j
100    l = i + m
       IF (rkq(i).GT.rkq(l)) GO TO 120
110    j = j + 1
       IF (j.GT.k) GO TO 80
       GO TO 90
120    t = rkq(i)
       rkq(i) = rkq(l)
       rkq(l) = t
       t = lapw%rk(i,ispin)
       lapw%rk(i,ispin) = lapw%rk(l,ispin)
       lapw%rk(l,ispin) = t
       itt = lapw%k1(i,ispin)
       lapw%k1(i,ispin) = lapw%k1(l,ispin)
       lapw%k1(l,ispin) = itt
       itt = lapw%k2(i,ispin)
       lapw%k2(i,ispin) = lapw%k2(l,ispin)
       lapw%k2(l,ispin) = itt
       itt = lapw%k3(i,ispin)
       lapw%k3(i,ispin) = lapw%k3(l,ispin)
       lapw%k3(l,ispin) = itt
       i = i - m
       IF (i.LT.1) GO TO 110
       GO TO 100
130    CONTINUE
       !+gu
       !--->    determine pairs of K-vectors, where K_z = K'_-z to use 
       !--->    z-reflection
       IF (l_zref) THEN
          n=0
          DO i=1,lapw%nv(ispin)
             DO j=1,i
                IF (((lapw%k1(i,ispin).EQ.lapw%k1(j,ispin)).AND.&
                     (lapw%k2(i,ispin).EQ.lapw%k2(j,ispin))).AND.&
                     (lapw%k3(i,ispin).EQ.-lapw%k3(j,ispin))) THEN
                   n=n+1 
                   lapw%matind(n,1)=i
                   lapw%matind(n,2)=j
                ENDIF
             ENDDO
          ENDDO
          nred=n

#ifdef CPP_MPI
          IF (n_size.GT.1) THEN
             !
             !--->     order K's in sequence K_1,...K_n | K_0,... | K_-1,....K_-n
             !
             ALLOCATE (pos(lapw%nv(ispin)))
             n_inner = lapw%nv(ispin) - nred
             IF (mod(nred,n_size).EQ.0) THEN
                n_bound = nred
             ELSE
                n_bound = (1+int( nred/n_size ))*n_size
             ENDIF
             IF (lapw%nv(ispin) - nred + n_bound.GT.dimension%nvd) THEN
                WRITE ( 6,*) 'increase dimension%nvd by:', lapw%nv(ispin)-nred+n_bound-dimension%nvd
                WRITE (16,*) 'increase dimension%nvd by:', lapw%nv(ispin)-nred+n_bound-dimension%nvd
                CALL juDFT_error("z-ref & ev || : dimension too small!" ,calledby ="apws")
             ENDIF

             i = 1
             j = 1
             DO n = 1, nred 
                IF (lapw%matind(n,1).EQ.lapw%matind(n,2)) THEN
                   pos(lapw%matind(n,1)) = n_inner + i
                   i = i + 1
                ELSE
                   pos(lapw%matind(n,1)) = j
                   pos(lapw%matind(n,2)) = j + n_bound
                   j = j + 1
                ENDIF
             ENDDO
             !--->          resort the rk,k1,k2,k3 and lapw%matind arrays:
             ALLOCATE (rk_help(lapw%nv(ispin)),k_help(3,lapw%nv(ispin)))
             DO n = 1, lapw%nv(ispin)
                rk_help(n)  = lapw%rk(n,ispin)
                k_help(1,n) = lapw%k1(n,ispin)
                k_help(2,n) = lapw%k2(n,ispin)
                k_help(3,n) = lapw%k3(n,ispin)
             ENDDO
             DO n = lapw%nv(ispin), 1, -1
                lapw%rk(pos(n),ispin) = rk_help(n)
                lapw%k1(pos(n),ispin) = k_help(1,n)
                lapw%k2(pos(n),ispin) = k_help(2,n)
                lapw%k3(pos(n),ispin) = k_help(3,n)
             ENDDO
             DO n = nred + 1, n_bound
                lapw%rk(n,ispin) = lapw%rk(lapw%nv(ispin),ispin)
                lapw%k1(n,ispin) = lapw%k1(lapw%nv(ispin),ispin)
                lapw%k2(n,ispin) = lapw%k2(lapw%nv(ispin),ispin)
                lapw%k3(n,ispin) = lapw%k3(lapw%nv(ispin),ispin)
             ENDDO
             DEALLOCATE (rk_help,k_help)
             DEALLOCATE (pos)
             lapw%nv(ispin) = lapw%nv(ispin) - nred + n_bound
          ENDIF
#endif
       ENDIF

       IF (noco%l_ss) THEN  ! sort additionally like in strgn1... gb
          i = 1
          gla = 0.
          gsk3(1) = 0.0
          eps=1.e-10
          DO  k = 1,lapw%nv(ispin)
             IF (rkq(k)-gla.GE.eps) i=i+1
             gla = rkq(k)
             gmi = (mk1+lapw%k1(k,ispin)) + (mk2+lapw%k2(k,ispin))*(2*mk1+1) +&
                  (mk3+lapw%k3(k,ispin))*(2*mk1+1)*(2*mk2+1)
             gsk3(k) = i * (9.+(2*mk1+1)*(2*mk2+1)*(2*mk3+1)) + gmi
          ENDDO
          CALL sort(lapw%nv(ispin),gsk3,index3)
          DO  k = 1,lapw%nv(ispin)
             k1rev(k) = lapw%k1(index3(k),ispin)
             k2rev(k) = lapw%k2(index3(k),ispin)
             k3rev(k) = lapw%k3(index3(k),ispin)
             gsk3(k) =  lapw%rk(index3(k),ispin)
          ENDDO
          DO  k = 1,lapw%nv(ispin)
             lapw%k1(k,ispin) = k1rev(k)
             lapw%k2(k,ispin) = k2rev(k)
             lapw%k3(k,ispin) = k3rev(k)
             lapw%rk(k,ispin) = gsk3(k)
          ENDDO
       ENDIF
       !-gu
    ENDDO

    IF ((.NOT. noco%l_ss) .AND. (input%jspins.EQ.2) ) THEN
       lapw%nv(input%jspins-(jspin-1)) = lapw%nv(jspin)
       DO i = 1,lapw%nv(jspin)
          lapw%rk(i,input%jspins-(jspin-1)) = lapw%rk(i,jspin)
          lapw%k1(i,input%jspins-(jspin-1)) = lapw%k1(i,jspin)
          lapw%k2(i,input%jspins-(jspin-1)) = lapw%k2(i,jspin)
          lapw%k3(i,input%jspins-(jspin-1)) = lapw%k3(i,jspin)
       ENDDO
       !      ELSE
       !         DO i = 1, min(nv(1),nv(jspd))
       !            WRITE(*,'(3i10,3x,3i10,3f12.6)') k1(i,1),k2(i,1),k3(i,1),
       !     +      k1(i,jspd),k2(i,jspd),k3(i,jspd),rk(i,1),rk(i,jspd),rkq(i)
       !         ENDDO
    ENDIF
    lapw%gvec(1,:,:)=lapw%k1
    lapw%gvec(2,:,:)=lapw%k2
    lapw%gvec(3,:,:)=lapw%k3

    DO ispin=1,SIZE(lapw%vk,3)
       DO k=1,lapw%nv(ispin)
          lapw%vk(:,k,ispin)=bkpt+lapw%gvec(:,k,ispin)+(ispin-1.5)*noco%qss
          lapw%gk(:,k,ispin)=MATMUL(TRANSPOSE(cell%bmat),lapw%vk(:,k,ispin))/MAX (lapw%rk(k,ispin),1.0e-30)
       ENDDO
    END DO

    lapw%num_local_cols=lapw%nv
    
  END SUBROUTINE apws
END MODULE m_apws
