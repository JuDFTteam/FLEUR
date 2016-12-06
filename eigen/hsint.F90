!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsint
CONTAINS
  SUBROUTINE hsint(input,noco,jij,stars, vpw,lapw,jspin,&
       n_size,n_rank,bkpt,cell,atoms,l_real,hamOvlp)
    !*********************************************************************
    !     initializes and sets up the hamiltonian and overlap matrices
    !     for the interstitial. only the lower triangle of the hermitian
    !     matrices are stored in compact real mode such that if h(i,j),
    !     i.ge.j, is hermitian and a is real, then
    !       a(i,j)=real( h(i,j) )  and  a(j,i)=aimag( h(i,j) )
    !                    m. weinert  1986
    !
    ! For the eigenvector parallelization each pe calculates an equal share
    ! of columns labeled nc. Then the starting element of a columns nc is
    !
    !     ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2
    !
    ! and, if a non-collinear matrix has to be set up, the starting column
    ! for the second spin-direction is
    !
    !     nc = int( 1. + (nv - n_rank - 1)/n_size ) + 1 .
    !
    ! For this direction, the outer loop starts at
    !
    ! istart = n_rank + (nc - 1)*n_size - nv .                        gb99
    !
    ! for a lo-calculation nv has to be replaced by nv+nlotot         gb01
    !
    !*********************************************************************
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_jij),INTENT(IN)        :: jij
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(INOUT)    :: lapw
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n_size,n_rank,jspin
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (INOUT) :: vpw(stars%n3d)
    REAL,    INTENT (IN)    :: bkpt(3) 
    LOGICAL,INTENT(IN)      :: l_real
    !     ..
    !     .. Local Scalars ..
    COMPLEX th,ts,phase
    REAL b1(3),b2(3),r2
    INTEGER i,i1,i2,i3,ii,in,j,ig3,ispin,l
    INTEGER istart,nc

    COMPLEX ust1,vp1
    COMPLEX, ALLOCATABLE :: vpw1(:)  ! for J constants
    !     ..
    ! ..
    if (l_real) THEN
       hamOvlp%a_r=0.0
       hamOvlp%b_r=0.0
    ELSE
       hamOvlp%a_c=0.0
       hamOvlp%b_c=0.0
    ENDIF
    ust1 = stars%ustep(1)
    ispin = jspin
    lapw%nmat = lapw%nv(ispin)

    !---> pk non-collinear
    IF (noco%l_noco) THEN
       !---> determine spin-up spin-up part of Hamiltonian- and overlapp-matrix
       !---> reload V_11
       READ (25) (vpw(ig3),ig3=1,stars%ng3)

       !--- J const
       IF( jij%l_J) THEN
          ALLOCATE ( vpw1(stars%n3d) )
          READ (25) (vpw1(ig3),ig3=1,stars%ng3)
       ENDIF
       !--- J const

       lapw%nmat = lapw%nv(1) + lapw%nv(2)
       ispin = 1

       !--- J const
       IF (jij%l_J) THEN
          DO i = 1,stars%ng3
             vpw(i) = (vpw(i) + vpw1(i))/2.
          END DO
       ENDIF
       !--- J const

       vp1 = REAL(vpw(1))
    ENDIF
    !---> pk non-collinear

    vp1 = vpw(1)
    !---> loop over (k+g')
    ii = 0
    DO  i = n_rank+1, lapw%nv(ispin), n_size
       !--->    loop over (k+g)
       DO  j = 1,i - 1
          ii = ii + 1
          !-->     determine index and phase factor
          i1 = lapw%k1(i,ispin) - lapw%k1(j,ispin)
          i2 = lapw%k2(i,ispin) - lapw%k2(j,ispin)
          i3 = lapw%k3(i,ispin) - lapw%k3(j,ispin)
          in = stars%ig(i1,i2,i3)
          IF (in.EQ.0) CYCLE
          phase = stars%rgphs(i1,i2,i3)
          !+APW_LO
          IF (input%l_useapw) THEN
             b1(1) = bkpt(1)+lapw%k1(i,ispin) ; b2(1) = bkpt(1)+lapw%k1(j,ispin)
             b1(2) = bkpt(2)+lapw%k2(i,ispin) ; b2(2) = bkpt(2)+lapw%k2(j,ispin)
             b1(3) = bkpt(3)+lapw%k3(i,ispin) ; b2(3) = bkpt(3)+lapw%k3(j,ispin)
             r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)   

             th = phase*(0.5*r2*stars%ustep(in)+vpw(in))
          ELSE
             th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,ispin)**2)*stars%ustep(in) + vpw(in))
          ENDIF
          !-APW_LO
          !--->    determine matrix element and store
          ts = phase*stars%ustep(in)
          if (l_real) THEN
          hamOvlp%a_r(ii) = REAL(th)
          hamOvlp%b_r(ii) = REAL(ts)
else
          hamOvlp%a_c(ii) = th
          hamOvlp%b_c(ii) = ts
endif
       ENDDO
       !--->    diagonal term (g-g'=0 always first star)
       ii = ii + 1
       if (l_real) THEN
       hamOvlp%a_r(ii) = 0.5*lapw%rk(i,ispin)*lapw%rk(i,ispin)*REAL(ust1) + REAL(vp1)
       hamOvlp%b_r(ii) = REAL(ust1)
else
       hamOvlp%a_c(ii) = 0.5*lapw%rk(i,ispin)*lapw%rk(i,ispin)*ust1 + vp1
       hamOvlp%b_c(ii) = ust1
endif
    ENDDO

    !---> pk non-collinear
    IF (noco%l_noco) THEN
       !+gb99
       nc = INT( 1. + (lapw%nv(1)+atoms%nlotot - n_rank - 1)/n_size )
       istart = n_rank + nc*n_size - (lapw%nv(1)+atoms%nlotot)
       !      ii = (nv(1)+nlotot+1)*(nv(1)+nlotot+2)/2 - 1
       ii = nc*(n_rank-n_size+1) + n_size*(nc+1)*nc/2 + lapw%nv(1)+atoms%nlotot
       !-gb99
       ispin = 2
       !---> determine spin-down spin-down part of Hamiltonian- and ovlp-matrix
       !---> reload V_22

       !--- J constants 
       IF(.NOT.jij%l_J) THEN
          READ (25) (vpw(ig3),ig3=1,stars%ng3)
          vp1 = REAL(vpw(1))
       ENDIF
       !--- J constants

       !---> loop over (k+g')
       DO i = istart+1, lapw%nv(ispin), n_size
          nc = nc + 1
          !--->    loop over (k+g)
          DO j = 1,i - 1
             !-gb99      ii = (nv(1)+i-1)*(nv(1)+i)/2 + nv(1) + j
             ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2 + lapw%nv(1)+atoms%nlotot + j
             !--->       determine index and phase factor
             i1 = lapw%k1(i,ispin) - lapw%k1(j,ispin)
             i2 = lapw%k2(i,ispin) - lapw%k2(j,ispin)
             i3 = lapw%k3(i,ispin) - lapw%k3(j,ispin)
             in = stars%ig(i1,i2,i3)
             IF (in.EQ.0) THEN
                WRITE (*,*) 'HSINT: G-G'' not in star i,j= ',i,j
             ELSE
                phase = stars%rgphs(i1,i2,i3)
                !+APW_LO
                IF (input%l_useapw) THEN
                   b1(1) = bkpt(1)+lapw%k1(i,ispin) ; b2(1) = bkpt(1)+lapw%k1(j,ispin)
                   b1(2) = bkpt(2)+lapw%k2(i,ispin) ; b2(2) = bkpt(2)+lapw%k2(j,ispin)
                   b1(3) = bkpt(3)+lapw%k3(i,ispin) ; b2(3) = bkpt(3)+lapw%k3(j,ispin)
                   r2 = DOT_PRODUCT(MATMUL(b2,cell%bbmat),b1)   

                   th = phase*( 0.5*r2*stars%ustep(in) + vpw(in) )
                ELSE
                   th = phase* (0.25* (lapw%rk(i,ispin)**2+lapw%rk(j,ispin)**2)*stars%ustep(in) + vpw(in))
                ENDIF
                !-APW_LO
                ts = phase*stars%ustep(in)
                hamOvlp%a_c(ii) = th
                hamOvlp%b_c(ii) = ts
             ENDIF
          ENDDO
          !--->    diagonal term (g-g'=0 always first star)
          !-gb99   ii = (nv(1)+i)*(nv(1)+i+1)/2
          ii = ii + 1
          hamOvlp%a_c(ii) = 0.5*lapw%rk(i,ispin)*lapw%rk(i,ispin)*ust1 + vp1
          hamOvlp%b_c(ii) = ust1
       ENDDO

       !---> determine spin-down spin-up part of Hamiltonian- and ovlp-matrix
       !---> reload real part of V_21
       READ (25) (vpw(ig3),ig3=1,stars%ng3)
       nc = INT( 1. + (lapw%nv(1)+atoms%nlotot - n_rank - 1)/n_size )
       !
       !---> loop over (k+g')
       !
       DO i = istart+1, lapw%nv(2), n_size
          nc = nc + 1
          !--->    loop over (k+g)
          DO j = 1,lapw%nv(1)
             !-gb99      ii = (nv(1)+i-1)*(nv(1)+i)/2 + j
             ii = (nc-1)*( n_rank - n_size + 1 ) + n_size*(nc-1)*nc/2 + j
             !--->       determine index and phase factor
             i1 = lapw%k1(i,2) - lapw%k1(j,1)
             i2 = lapw%k2(i,2) - lapw%k2(j,1)
             i3 = lapw%k3(i,2) - lapw%k3(j,1)
             in = stars%ig(i1,i2,i3)
             IF (in.EQ.0) THEN
                WRITE (*,*) 'HSINT: G-G'' not in star i,j= ',i,j
             ELSE
                hamOvlp%a_c(ii) = stars%rgphs(i1,i2,i3)*vpw(in) 
                !--- J constants 
                IF(jij%l_J) THEN
                   hamOvlp%a_c(ii) = 0
                ENDIF
                !--- J constants

             ENDIF
          ENDDO
       ENDDO
       !---> pk non-collinear
    ENDIF

    IF (jij%l_J) DEALLOCATE (vpw1)

    RETURN
  END SUBROUTINE hsint
END MODULE m_hsint
