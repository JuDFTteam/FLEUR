!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_lapwdim
CONTAINS
  SUBROUTINE lapw_dim(kpts,cell,input,noco,nococonv,forcetheo,atoms,nbasfcn,juPhon)
    !
    !*********************************************************************
    !     determines dimensions of the lapw basis set with |k+G|<rkmax.
    !  Generalization of the old apws_dim routine
    !*********************************************************************
    USE m_boxdim
    USE m_types_fleurinput
    USE m_types_forcetheo_extended
    IMPLICIT NONE
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_nococonv),INTENT(IN)  :: nococonv
    TYPE(t_noco),INTENT(IN)      :: noco
     
    CLASS(t_forcetheo),INTENT(IN):: forcetheo
    TYPE(t_atoms),INTENT(IN)     :: atoms
    INTEGER, INTENT(OUT)         :: nbasfcn

    TYPE(t_juPhon), INTENT(IN)   :: juPhon

    !local variable for init
    INTEGER               :: nvd,nv2d,addx,addy,addz
    TYPE(t_lapw) :: lapw

    INTEGER j1,j2,j3,mk1,mk2,mk3,iofile,ksfft,q,nk,nv,nv2
    INTEGER ispin,nvh(2),nv2h(2)

    REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3),gmaxp,qss(3)
    REAL,ALLOCATABLE:: q_vectors(:,:)
    REAL            :: bkpt(3)
    ! ..
    !
    !------->          ABBREVIATIONS
    !

    !   iofile      : device number for in and output
    !   gmax        : cut-off wavevector for charge density
    !   rkmax       : cut-off for |g+k|
    !   gmaxp       : gmaxp = gmax/rkmax, ideal: gmaxp=2
    !   arltv(i)    : length of reciprical lattice vector along
    !                 direction (i)
    !
    !---> Determine rkmax box of size mk1, mk2, mk3,
    !     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
    !
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

    !     (add 1+1 due to integer rounding, strange k_vector in BZ)
    mk1 = int(input%rkmax/arltv1) + 2
    mk2 = int(input%rkmax/arltv2) + 2
    mk3 = int(input%rkmax/arltv3) + 2

    rkm = input%rkmax
    rk2 = rkm*rkm

    !Determine the q-vector(s) to use
    IF (juPhon%l_dfpt) THEN
       ALLOCATE(q_vectors(3,SIZE(juPhon%qvec,2)))
       q_vectors=juPhon%qvec
    ELSE
       SELECT TYPE(forcetheo)
       TYPE IS (t_forcetheo_ssdisp)
          ALLOCATE(q_vectors(3,SIZE(forcetheo%qvec,2)))
          q_vectors=forcetheo%qvec
       TYPE IS (t_forcetheo_dmi)
          ALLOCATE(q_vectors(3,SIZE(forcetheo%qvec,2)))
          q_vectors=forcetheo%qvec
       TYPE IS (t_forcetheo_jij)
          ALLOCATE(q_vectors(3,SIZE(forcetheo%qvec,2)))
          q_vectors=forcetheo%qvec
       CLASS IS (t_forcetheo) ! DEFAULT
          ALLOCATE(q_vectors(3,1))
          q_vectors(:,1)=nococonv%qss
       END SELECT
    END IF

    if (any(abs(nococonv%qss-q_vectors(:,1))>1E-4)) CALL judft_warn("q-vector for self-consistency should be first in list for force-theorem")


    nvd = 0 ; nv2d = 0
    DO q=1,SIZE(q_vectors,2)
       qss=q_vectors(:,q)
       DO nk=1,kpts%nkpt
          bkpt=kpts%bk(:,nk)
          !---> obtain vectors
          !---> in a spin-spiral calculation different basis sets are used for
          !---> the two spin directions, because the cutoff radius is defined
          !---> by |G + k +/- qss/2| < rkmax.
          DO ispin = 1,2
             addX = abs(NINT((bkpt(1) + (2*ispin - 3)/2.0*qss(1))/arltv1))
             addY = abs(NINT((bkpt(2) + (2*ispin - 3)/2.0*qss(2))/arltv2))
             addZ = abs(NINT((bkpt(3) + (2*ispin - 3)/2.0*qss(3))/arltv2))
             nv = 0
             nv2 = 0
             DO j1 = -mk1-addX,mk1+addX
                s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*qss(1)
                DO j2 = -mk2-addY,mk2+addY
                   s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*qss(2)
                   !--->          nv2 for films
                   s(3) = 0.0
                   !r2 = dotirp(s,s,cell%bbmat)
                   r2 = dot_product(matmul(s,cell%bbmat),s)
                   IF (r2.LE.rk2) nv2 = nv2 + 1
                   DO j3 = -mk3-addz,mk3+addz
                      s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
                      !r2 = dotirp(s,s,cell%bbmat)
                      r2 = dot_product(matmul(s,cell%bbmat),s)
                      IF (r2.LE.rk2) THEN
                         nv = nv + 1
                      END IF
                   END DO
                END DO
             END DO
        
             nvh(ispin)  = nv
             nv2h(ispin) = nv2
          END DO
          nvd=MAX(nvd,MAX(nvh(1),nvh(2)))
          nv2d=MAX(nv2d,MAX(nv2h(1),nv2h(2)))

       ENDDO !k-loop
    ENDDO !q-loop

    nbasfcn = nvd + atoms%nat*atoms%nlod*(2*atoms%llod+1)
    IF (noco%l_noco) nbasfcn = 2*nbasfcn
    call lapw%init_dim(nvd,nv2d,nbasfcn)

  END SUBROUTINE lapw_dim

  
END MODULE m_lapwdim
