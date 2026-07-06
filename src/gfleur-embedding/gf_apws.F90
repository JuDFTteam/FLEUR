!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_apws
      use m_juDFT
      USE m_constants, ONLY: oUnit
      IMPLICIT NONE

      PRIVATE
      PUBLIC gf_apws, gf_apws_DIM
      CONTAINS
      SUBROUTINE gf_apws(jspins,jspin,bk,lpe0,atoms,input,sym,noco,      &
     &     nococonv,embinp,cell,lapw,lapw_gf,layer)

!*********************************************************************
!     determines the lapw list such that |k+G|<rkmax.
!     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.
!        m. weinert  1986
!*********************************************************************
!     changed subroutine to generate g-vectors out of a zylinder not
!     the usual sphere.
!
!     added additional index kp of parallel g-vectors to handle
!     matrixes expanded in terms of the 2D Basis as used in Embedding
!     New variables:
!     kp :index (now in the modern t_lapw)
!     nv2:no of 2D Basis fcts (in t_lapw_gf)
!            Daniel Wortmann, Tokyo, 2001
!
!     In the port to current FLEUR this routine fills a modern t_lapw
!     (gvec/rk/gk/vk/nv/nmat + the standard LO bookkeeping via
!     lapw%init_lo) so that the standard kernels (hsmt, hsmt_ab,
!     fjgj%calculate) consume the cylinder basis unchanged; the
!     GF-specific 2D bookkeeping lives in t_lapw_gf.
!*******************************************************************
      USE m_gf_types
      USE m_boxdim
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT  (IN)      :: jspins,jspin
      LOGICAL , INTENT  (IN)     :: lpe0
      TYPE(t_lapw),INTENT(INOUT) :: lapw
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_embinp),INTENT(IN)  :: embinp
      TYPE(t_cell),INTENT(IN)    :: cell
      TYPE(t_noco),INTENT(IN)    :: noco
      TYPE(t_nococonv),INTENT(IN):: nococonv
      REAL   ,INTENT(IN)         :: bk(3)
      INTEGER,INTENT(IN)         :: layer

!     ..
!     .. Local Scalars ..
      REAL arltv1,arltv2,arltv3,s3,s2,rk2,rkm,t
      INTEGER i,itt,j,j1,j2,j3,ki,l,m,mk1,mk2,mk3,n,                    &
     &        ispin,jsp_start,jsp_end,n2,k
      LOGICAL l_sorted
!     ..
!     .. Local Arrays ..
      REAL s(3)
!     ..
      !<-- allocate the basis arrays if needed
      IF (.NOT.ALLOCATED(lapw%gvec)) THEN
         ALLOCATE(lapw%gvec(3,lapw_gf%nvd,jspins));lapw%gvec=0
         ALLOCATE(lapw%rk(lapw_gf%nvd,jspins));lapw%rk=0.
         ALLOCATE(lapw%gk(3,lapw_gf%nvd,jspins));lapw%gk=0.
         ALLOCATE(lapw%vk(3,lapw_gf%nvd,jspins));lapw%vk=0.
         ALLOCATE(lapw%kp(lapw_gf%nvd,jspins));lapw%kp=0
         ALLOCATE(lapw_gf%rkp(lapw_gf%nv2d,jspins));lapw_gf%rkp=0.
         ALLOCATE(lapw_gf%k1p(lapw_gf%nv2d,jspins));lapw_gf%k1p=0
         ALLOCATE(lapw_gf%k2p(lapw_gf%nv2d,jspins));lapw_gf%k2p=0
      ENDIF
      !>
      !initialize
      lapw%nv=0
      lapw%nlotot=atoms%nlotot
      lapw%bkpt=bk
      lapw%qphon=0.0
      lapw_gf%nv2=0


      !<-- Calculate the Box

!---> Determine rkmax box of size mk1, mk2, mk3,
!     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
!     arltv(i) length of reciprical lattice vector along direction (i)
!
      !cylinder setup is used if the projected Green fct's are
      !calculated


      CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

      mk1 = int( lapw_gf%rkmax/arltv1 ) + 4
      mk2 = int( lapw_gf%rkmax/arltv2 ) + 4

      IF (lapw_gf%l_cylinder) THEN
         mk3 = embinp%napw(layer)
      ELSE
         mk3=int( lapw_gf%rkmax/arltv3 ) + 4
      ENDIF

      !>
      ! No noco-setup in this subroutine !!
      jsp_start = jspin
      jsp_end   = jspin
      rkm = lapw_gf%rkmax

      DO ispin = jsp_start,jsp_end
         rk2 = rkm*rkm
!
!--->    2D obtain vectors
!     First calculate only 2D Vectors, sort them and then add z-komponen
!
         n2= 0
         s(3)=0
         DO  j1 = -mk1,mk1
            s(1) = bk(1) + j1 + (2*ispin - 3)/2.0*nococonv%qss(1)
            DO  j2 = -mk2,mk2
               s(2) = bk(2) + j2 + (2*ispin - 3)/2.0*nococonv%qss(2)
               s2=dot_product(s,matmul(s,cell%bbmat))
               IF (s2<=rk2) THEN
                  n2=n2+1
                  lapw_gf%rkp(n2,ispin)=sqrt(s2)
                  lapw_gf%k1p(n2,ispin)=j1
                  lapw_gf%k2p(n2,ispin)=j2
               ENDIF
            ENDDO
         ENDDO

         lapw_gf%nv2(ispin)=n2
         IF (n2>lapw_gf%nv2d)CALL juDFT_error('gf_apws: wrong 2D-dimensions')
         IF (lapw_gf%l_cylinder)  lapw%nv(ispin) = n2*(2*mk3+1)
!
! Sort the 2D-PW
!
         l_sorted=.FALSE.
         DO WHILE(.NOT.l_sorted)
            l_sorted=.TRUE.
            DO i=1,n2-1
               !reorderer such that |G(i)|<|G(i+1)| and
               !if |G(i)|=|G(i+1)| and G_x(i)<G_x(i+1) and
               !if |G(i)|=|G(i+1)| and G_x(i)=G_x(i+1) G_y(i)<G_y(i+1)
               IF (                                                     &
     &              (lapw_gf%rkp(i,ispin)-lapw_gf%rkp(i+1,ispin)        &
     &              )>EPSILON(1.0).OR.((ABS(lapw_gf%rkp(i,ispin)        &
     &              -lapw_gf%rkp(i+1,ispin))<EPSILON(1.0)).AND.         &
     &               ((lapw_gf%k1p(i,ispin)>lapw_gf%k1p(i+1,ispin)).OR. &
     &               (lapw_gf%k1p(i,ispin) == lapw_gf%k1p(i+1,ispin).AND.&
     &               (lapw_gf%k2p(i,ispin)>lapw_gf%k2p(i+1,ispin))))))  &
     &              THEN
                  t = lapw_gf%rkp(i,ispin)
                  lapw_gf%rkp(i,ispin) = lapw_gf%rkp(i+1,ispin)
                  lapw_gf%rkp(i+1,ispin) = t
                  itt = lapw_gf%k1p(i,ispin)
                  lapw_gf%k1p(i,ispin)=lapw_gf%k1p(i+1,ispin)
                  lapw_gf%k1p(i+1,ispin)=itt
                  itt =lapw_gf%k2p(i,ispin)
                  lapw_gf%k2p(i,ispin) = lapw_gf%k2p(i+1,ispin)
                  lapw_gf%k2p(i+1,ispin) = itt
                  l_sorted=.FALSE.
               ENDIF
            ENDDO
         ENDDO

         ! now generate a map-to the global g-vectors
         lapw_gf%g2map(:)=0
         do n2=1,lapw_gf%nv2(ispin)
            lapw_gf%g2map(n2)=lapw_gf%global2Dmap(lapw_gf%k1p(n2,ispin),lapw_gf%k2p(n2,ispin))
         enddo
         IF (any(lapw_gf%g2map(:lapw_gf%nv2(ispin))==0)) THEN
            DO n2=1,lapw_gf%nv2(ispin)
                 write(*,*) n2,lapw_gf%g2map(n2),lapw_gf%k1p(n2,ispin),lapw_gf%k2p(n2,ispin)
            ENDDO
            CALL juDFT_error("Mapping failed in gf_apws")
         ENDIF

!
!     now add third component of g-vector
!
         n=0
         DO n2=1, lapw_gf%nv2(ispin)
            DO j3 = -mk3,mk3
               s(1) = bk(1)+lapw_gf%k1p(n2,ispin)+(2*ispin-3)/2.0       &
     &              *nococonv%qss(1)
               s(2) = bk(2)+lapw_gf%k2p(n2,ispin)+(2*ispin-3)/2.0       &
     &              *nococonv%qss(2)
               s(3) = bk(3) + j3 + (2*ispin - 3)/2.0*nococonv%qss(3)
               s3=dot_product(s,matmul(s,cell%bbmat))
               IF ((.NOT.lapw_gf%l_cylinder).AND.s3>rk2) CYCLE
               n = n + 1
!     new stop added ! no reduction of gmax anymore!
               IF (n>lapw_gf%nvd)                                       &
     &             CALL juDFT_error('Too many basis functions !',calledby="gf_apws")
               lapw%gvec(1,n,ispin) = lapw_gf%k1p(n2,ispin)
               lapw%gvec(2,n,ispin) = lapw_gf%k2p(n2,ispin)
               lapw%gvec(3,n,ispin) = j3
               lapw%kp(n,ispin) = n2
               lapw%rk(n,ispin) = SQRT(s3)
            ENDDO
         ENDDO
         IF (.NOT.lapw_gf%l_cylinder)  lapw%nv(ispin)=n
!--->    sort by shell-metzner
         m = n
   81    m = m/2
         IF (m<=0) GO TO 131
         ki = n - m
         j = 1
   91    i = j
  101    l = i + m
         IF (lapw%rk(i,ispin)>lapw%rk(l,ispin)) GO TO 121
  111    j = j + 1
         IF (j>ki) GO TO 81
         GO TO 91
  121    t = lapw%rk(i,ispin)
         lapw%rk(i,ispin) = lapw%rk(l,ispin)
         lapw%rk(l,ispin) = t
         itt = lapw%kp(i,ispin)
         lapw%kp(i,ispin)=lapw%kp(l,ispin)
         lapw%kp(l,ispin)=itt
         itt = lapw%gvec(1,i,ispin)
         lapw%gvec(1,i,ispin) = lapw%gvec(1,l,ispin)
         lapw%gvec(1,l,ispin) = itt
         itt = lapw%gvec(2,i,ispin)
         lapw%gvec(2,i,ispin) = lapw%gvec(2,l,ispin)
         lapw%gvec(2,l,ispin) = itt
         itt = lapw%gvec(3,i,ispin)
         lapw%gvec(3,i,ispin) = lapw%gvec(3,l,ispin)
         lapw%gvec(3,l,ispin) = itt
         i = i - m
         IF (i<1) GO TO 111
         GO TO 101
  131    CONTINUE

         !Now check for  largest k+g in sphere
         IF (lapw_gf%l_cylinder) THEN
            DO n=lapw%nv(ispin),1,-1
                IF (lapw%rk(n,jspin)<lapw_gf%rkmax) exit
            ENDDO
            lapw_gf%nv_sphere(ispin)=n
            write(oUnit,*) "Basis setup:"
            write(oUnit,*) "Cylinder:,",lapw%nv(ispin)," sphere:",lapw_gf%nv_sphere(ispin)
         ELSE
            lapw_gf%nv_sphere(ispin)=lapw%nv(ispin)
         ENDIF

         IF ((.NOT. noco%l_ss) .AND. (jspins==2) ) THEN
             lapw%nv(jspins-(jspin-1)) =  lapw%nv(jspin)
            DO i = 1, lapw%nv(jspin)
               lapw%rk(i,jspins-(jspin-1)) = lapw%rk(i,jspin)
               lapw%gvec(:,i,jspins-(jspin-1)) = lapw%gvec(:,i,jspin)
               lapw%kp(i,jspins-(jspin-1)) = lapw%kp(i,jspin)
            ENDDO
             lapw_gf%nv2(jspins-(jspin-1)) =  lapw_gf%nv2(jspin)
             lapw_gf%nv_sphere(jspins-(jspin-1)) =  lapw_gf%nv_sphere(jspin)
             DO i = 1, lapw_gf%nv2(jspin)
                lapw_gf%rkp(i,jspins-(jspin-1))=lapw_gf%rkp(i,jspin)
                lapw_gf%k1p(i,jspins-(jspin-1))=lapw_gf%k1p(i,jspin)
                lapw_gf%k2p(i,jspins-(jspin-1))=lapw_gf%k2p(i,jspin)
             ENDDO
         ENDIF


                                !spin loop
      ENDDO

      !<-- complete the modern t_lapw: vk/gk, k1..k3 copies, LO setup
      DO ispin = 1,jspins
         DO k = 1,lapw%nv(ispin)
            lapw%vk(:,k,ispin) = lapw%bkpt + lapw%gvec(:,k,ispin)        &
     &           + (ispin-1.5)*nococonv%qss
            lapw%gk(:,k,ispin) = MATMUL(TRANSPOSE(cell%bmat),            &
     &           lapw%vk(:,k,ispin))/MAX(lapw%rk(k,ispin),1.0e-30)
         ENDDO
      ENDDO
      IF (ALLOCATED(lapw%k1)) DEALLOCATE(lapw%k1,lapw%k2,lapw%k3)
      ALLOCATE(lapw%k1(SIZE(lapw%gvec,2),jspins))
      ALLOCATE(lapw%k2(SIZE(lapw%gvec,2),jspins))
      ALLOCATE(lapw%k3(SIZE(lapw%gvec,2),jspins))
      lapw%k1 = lapw%gvec(1,:,:)
      lapw%k2 = lapw%gvec(2,:,:)
      lapw%k3 = lapw%gvec(3,:,:)
      lapw%num_local_cols = lapw%nv

      IF (ANY(atoms%nlo > 0))                                            &
     &     CALL lapw%init_lo(atoms,input,sym,noco,nococonv,cell)

      lapw%nv_tot = lapw%nv(1)
      lapw%nmat = lapw%nv(1) + atoms%nlotot
      IF (noco%l_noco) lapw%nv_tot = lapw%nv_tot + lapw%nv(2)
      IF (noco%l_noco) lapw%nmat = lapw%nv_tot + 2*atoms%nlotot
      !>

      !<--Set the GF-dimensions correctly for noco
      !   (the *_sphere counts are plane-wave only, LOs excluded)
      lapw_gf%nmat = lapw%nv(jspin)
      lapw_gf%nmat_sphere = lapw_gf%nv_sphere(jspin)
      IF (noco%L_noco) lapw_gf%nmat = lapw%nv(1)+lapw%nv(2)
      IF (noco%L_noco) lapw_gf%nmat_sphere = lapw_gf%nv_sphere(1)+lapw_gf%nv_sphere(2)

      IF (noco%l_noco) THEN
         lapw_gf%nv_tot = lapw%nv(1)+lapw%nv(2)
         lapw_gf%nv_tot_sphere = lapw_gf%nv_sphere(1)+lapw_gf%nv_sphere(2)
         lapw_gf%nv2_tot = lapw_gf%nv2(1)+lapw_gf%nv2(2)
      ELSE
         lapw_gf%nv_tot = lapw%nv(jspin)
         lapw_gf%nv_tot_sphere = lapw_gf%nv_sphere(jspin)
         lapw_gf%nv2_tot = lapw_gf%nv2(jspin)
      ENDIF
      !>


      RETURN
      END SUBROUTINE gf_apws
      !<-- S:gf_apws_dim(qss,kpts,embinp,cell,lapw_gf)

      SUBROUTINE gf_apws_DIM(jspins,                                    &
     &     qss,kpts,embinp,cell,rkmax,rmtlmax,                          &
     &     lapw_gf,layer)
!******************************************
!     Calculates the max-no of lapw's
!     D. Wortmann
!******************************************
      USE m_gf_types
      USE m_boxdim
      IMPLICIT NONE
      !<--Arguments

      INTEGER,INTENT(IN)            :: jspins
      REAL   ,INTENT(IN)            :: qss(3),rmtlmax(2)
      REAL   ,INTENT(IN)            :: rkmax
      TYPE(t_kpts),INTENT(IN)       :: kpts
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_embinp),INTENT(INOUT)  :: embinp
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      INTEGER,INTENT(IN)            :: layer

      !>
      !<--Locals

      INTEGER :: mk1,mk2,mk3
      REAL    :: s(3),rk2,s2,s3
      INTEGER :: ispin,nk,j1,j2,j3,n,n2,nn
      INTEGER,ALLOCATABLE :: k1p(:),k2p(:)
      REAL                :: gkmax,gkmax2
      !>

      lapw_gf%rkmax = rkmax
      embinp%napw(layer)=ABS(embinp%napw(layer))

      !<--Determine maximal boxsize!

      IF (embinp%napw(layer)/=0) THEN
         lapw_gf%l_cylinder=.TRUE.
      ELSE
         lapw_gf%l_cylinder=.FALSE.
      ENDIF
      CALL boxdim(cell%bmat,s(1),s(2),s(3))
      mk1 = int( lapw_gf%rkmax/s(1) ) + 4
      mk2 = int( lapw_gf%rkmax/s(2) ) + 4
      IF (lapw_gf%l_cylinder) THEN
         mk3 = embinp%napw(layer)
      ELSE
         mk3=INT( lapw_gf%rkmax/s(3) ) + 4
      ENDIF
      ALLOCATE(k1p((2*mk1+1)*(2*mk2+1)),k2p((2*mk1+1)*(2*mk2+1)))

      !>

      !<--loop over all kpts and spins

      lapw_gf%nvd=0
      lapw_gf%nv2d=0
      gkmax2 = 0.0
      gkmax = 0.0
                   !Only needed for spin-spiral
      DO ispin=1,2
         DO nk=1,kpts%nkpt
            !<--Calculate 2D-Basis functions!

            rk2 = lapw_gf%rkmax*lapw_gf%rkmax
            n2= 0
            s(3)=0
            DO  j1 = -mk1,mk1
               s(1) = kpts%bk(1,nk) + j1 + (2*ispin - 3)/2.0*qss(1)
               DO  j2 = -mk2,mk2
                  s(2) = kpts%bk(2,nk) + j2 + (2*ispin - 3)/2.0*qss(2)
                  s2=dot_product(s,matmul(s,cell%bbmat))
                  IF (s2<=rk2) THEN
                     gkmax2=max(gkmax2,s2)
                     n2=n2+1
                     k1p(n2)=j1
                     k2p(n2)=j2
                  ENDIF
               ENDDO
            ENDDO
            lapw_gf%nv2d=max(n2,lapw_gf%nv2d)

            !>
            IF (lapw_gf%l_cylinder) THEN
               lapw_gf%nvd=max(lapw_gf%nvd,lapw_gf%nv2d*(2*embinp%napw(layer)+1))
            ELSE
               !<--Calculate the 3d-Basis

               n=0
               DO nn=1, n2
                  DO j3 = -mk3,mk3
                     s(1) = kpts%bk(1,nk)+k1p(nn)+(2*ispin-3)/2.*qss(1)
                     s(2) = kpts%bk(2,nk)+k2p(nn)+(2*ispin-3)/2.*qss(2)
                     s(3) = kpts%bk(3,nk)+j3 +(2*ispin-3)/2. *qss(3)
                     s3 = dot_PRODUCT(s,MATMUL(s,cell%bbmat))
                     IF (s3>rk2) CYCLE
                     gkmax=max(s3,gkmax)
                                                                 !Set na
                     embinp%napw(layer)=max(j3,embinp%napw(layer))
                                                 !maximum if in sc-mode
                     n = n + 1
                  ENDDO
               ENDDO
               lapw_gf%nvd=max(n,lapw_gf%nvd)

               !>
            ENDIF
         ENDDO
      ENDDO



      !>
      !>
      DEALLOCATE(k1p,k2p)
      WRITE(oUnit,*) "Basis setup:"
      WRITE(oUnit,*) "Maximal dimensions:"
      WRITE(oUnit,*) "nvd:",lapw_gf%nvd
      WRITE(oUnit,*) "nv2d:",lapw_gf%nv2d
      WRITE(oUnit,*) "Matching at MT:"
      IF (gkmax == 0.0)                                                 &
     &     gkmax = sqrt(gkmax2**2+(cell%bmat(3,3)*embinp%napw(layer))**2)
      WRITE(oUnit,*) "Maximal |k+g|=",gkmax
      WRITE(oUnit,*) "Maximal r*l  =",rmtlmax(2)
      WRITE(oUnit,*) "Minimal r*l  =",rmtlmax(1)

      CALL gf_apws_globalmap(cell,lapw_gf)

      END SUBROUTINE

      !>

      !<-- S:gf_apws_globalmap(cell,lapw_gf)

      SUBROUTINE gf_apws_globalmap(cell,                                &
     &     lapw_gf)
!******************************************
!   Creates the mapping list of the g-vectors needed to store the
!   embedding potentials
!******************************************
      USE m_gf_types
      USE m_boxdim
      USE m_gf_math,ONLY:sort
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      !>
      !<--Locals

      INTEGER :: mk1,mk2
      REAL    :: s(3),rk2,s_min
      INTEGER :: j1,j2,n,n2
      INTEGER,ALLOCATABLE :: k1p(:),k2p(:),sindex(:)
      REAL   ,ALLOCATABLE :: kl(:)

      !>


      !<-- Determine maximal boxsize!

      CALL boxdim(cell%bmat,s(1),s(2),s(3))
      mk1 = int( lapw_gf%rkmax/s(1) ) + 4
      mk2 = int( lapw_gf%rkmax/s(2) ) + 4

      ALLOCATE(sindex((2*mk1+1)*(2*mk2+1)),kl((2*mk1+1)*(2*mk2+1))      &
     &     ,k1p((2*mk1+1)*(2*mk2+1)),k2p((2*mk1+1)*(2*mk2+1)))

      !>

      !<-- Calculate 2D-g vectors!
      rk2 = lapw_gf%rkmax*lapw_gf%rkmax
      n2 = 0
      s(3)= 0
      DO  j1 = -mk1,mk1
         DO  j2 = -mk2,mk2
            !Try all combinations with maximal k-vectors
            s(1:2) = (/j1,j2/)
            s_min  = dot_PRODUCT(s,MATMUL(s,cell%bbmat))
            s(1:2) = (/j1+0.5,j2+0.0/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1+0.5,j2+0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1+0.0,j2+0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1-0.5,j2+0.0/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1-0.5,j2-0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1+0.0,j2-0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1+0.5,j2-0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            s(1:2) = (/j1-0.5,j2+0.5/)
            s_min  = MIN(s_min,dot_PRODUCT(s,MATMUL(s,cell%bbmat)))
            IF (s_min <= rk2) THEN
               !We found a g-vector
               n2 = n2+1
               k1p(n2) = j1
               k2p(n2) = j2
               kl(n2) = sqrt(1.0*(j1**2+j2**2))
            ENDIF
         ENDDO
      ENDDO
      !>
      !<-- Sort and create map (use sort from gf_tools)
      CALL sort(kl(:n2),sindex(:n2))

      ALLOCATE(lapw_gf%global2Dlist(2,n2))
      ALLOCATE(lapw_gf%g2map(n2))

      lapw_gf%global2Dlist(1,:) = k1p(sindex(:n2))
      lapw_gf%global2Dlist(2,:) = k2p(sindex(:n2))
      ALLOCATE(lapw_gf%global2Dmap(-mk1:mk1,-mk2:mk2))
      lapw_gf%global2Dmap = 0
      DO n = 1,n2
         lapw_gf%global2Dmap(lapw_gf%global2Dlist(1,n),lapw_gf%global2Dlist(2,n))&
     &        = n
      ENDDO
      !>

      DO N = 1,N2
         WRITE(oUnit,*) N,lapw_gf%global2Dlist(:,n)
         IF (lapw_gf%global2Dmap(lapw_gf%global2Dlist(1,n),lapw_gf%global2Dlist(2&
     &        ,n)) /= n) CALL juDFT_error                                 &
     &        ("Internal error in gf_apws_globalmap")
      ENDDO
      lapw_gf%g_max = (/MAXVAL(lapw_gf%global2Dlist(1,:))                &
     &     ,MAXVAL(lapw_gf%global2Dlist(2,:))/)
      DEALLOCATE(k1p,k2p,kl,sindex)

      END SUBROUTINE

      !>

      END
