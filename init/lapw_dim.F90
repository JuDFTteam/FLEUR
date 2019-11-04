!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_lapwdim
CONTAINS
  SUBROUTINE lapw_dim(kpts,cell,input,noco,oneD,forcetheo,atoms)
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
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_oneD),INTENT(IN)      :: oneD
    CLASS(t_forcetheo),INTENT(IN):: forcetheo
    TYPE(t_atoms),INTENT(IN)     :: atoms

    !local variable for init
    INTEGER               :: nvd,nv2d,nbasfcn
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
       q_vectors(:,1)=noco%qss
    END SELECT

    if (any(abs(noco%qss-q_vectors(:,1))>1E-4)) CALL judft_warn("q-vector for self-consistency should be first in list for force-theorem")


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
             nv = 0
             nv2 = 0
             DO j1 = -mk1,mk1
                s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*qss(1)
                DO j2 = -mk2,mk2
                   s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*qss(2)
                   !--->          nv2 for films
                   s(3) = 0.0
                   !r2 = dotirp(s,s,cell%bbmat)
                   r2 = dot_product(matmul(s,cell%bbmat),s)
                   IF (r2.LE.rk2) nv2 = nv2 + 1
                   DO j3 = -mk3,mk3
                      s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
                      !r2 = dotirp(s,s,cell%bbmat)
                      r2 = dot_product(matmul(s,cell%bbmat),s)
                      IF (r2.LE.rk2) THEN
                         nv = nv + 1
                      END IF
                   END DO
                END DO
             END DO
             !-odim
             IF (oneD%odd%d1) THEN
                nv2 = 0
                s(1) = 0.0
                s(2) = 0.0
                DO j3 = -mk3,mk3
                   s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*qss(3)
                   !r2 = dotirp(s,s,cell%bbmat)
                   r2 = dot_product(matmul(s,cell%bbmat),s)

                   IF (r2.LE.rk2) THEN
                      nv2 = nv2 + 1
                   END IF
                END DO
             END IF
             !+odim
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

  SUBROUTINE lapw_fft_dim(cell,input,noco,stars)
    !
    !*********************************************************************
    !     determines dimensions of the lapw basis set with |k+G|<rkmax.
    !  Generalization of the old apws_dim routine
    !*********************************************************************
    USE m_boxdim
    USE m_ifft,     ONLY : ifft235
    USE m_types

    IMPLICIT NONE
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_stars),INTENT(INOUT)  :: stars

    INTEGER j1,j2,j3,mk1,mk2,mk3,iofile,ksfft,q,nk,nv,nv2
    INTEGER ispin,nvh(2),nv2h(2)

    REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3),gmaxp
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



    !---> Determine the dimensions kq1d, kq2d, kq3d
    !     of the dimension of the charge density fft-box
    !     needed for the fast calculation of pw density
    !     (add 1 due to integer rounding,
    !      factor 2 due to positive domain)
    !
    gmaxp = 2.0
    CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)
    !
    mk1 = int(gmaxp*input%rkmax/arltv1) + 1
    mk2 = int(gmaxp*input%rkmax/arltv2) + 1
    mk3 = int(gmaxp*input%rkmax/arltv3) + 1

    !---> add + 1 in spin spiral calculation, to make sure that all G's are
    !---> still within the FFT-box after being shifted by the spin spiral
    !---> q-vector.
    IF (noco%l_ss) THEN
       mk1 = mk1 + 1
       mk2 = mk2 + 1
       mk3 = mk3 + 1
    ENDIF
    !
    stars%kq1_fft = 2*mk1
    stars%kq2_fft = 2*mk2
    stars%kq3_fft = 2*mk3
    !
    !---> fft's are usually fastest for low primes
    !     (restrict kqid to: kqid=  (2**P) * (3**Q) * (5**R)
    !
    ksfft = 1
    !     ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
    !                      0  PROGRAM, RADIX-2 ONLY
    !                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5

    stars%kq1_fft = ifft235(6,ksfft,stars%kq1_fft,gmaxp)
    stars%kq2_fft = ifft235(6,ksfft,stars%kq2_fft,gmaxp)
    stars%kq3_fft = ifft235(6,ksfft,stars%kq3_fft,gmaxp)


  END SUBROUTINE lapw_fft_dim
END MODULE m_lapwdim
