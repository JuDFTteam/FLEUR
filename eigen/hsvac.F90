!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsvac
  USE m_juDFT
CONTAINS
  !-----------------------------------------------------------
  !Calculate the vacuum contribution to the Hamiltonian and
  !Overlap matrix
  !-----------------------------------------------------------
  SUBROUTINE hsvac(&
       vacuum,stars,DIMENSION, atoms, jsp,input,vxy,vz,evac,cell,&
       bkpt,lapw,sym, noco,jij, n_size,n_rank,nv2,l_real,hmat,smat)
 

    USE m_vacfun
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)  :: DIMENSION
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_vacuum),INTENT(IN)     :: vacuum
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_jij),INTENT(IN)        :: jij
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_lapwmat),INTENT(INOUT) :: hmat,smat
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp   ,n_size,n_rank
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (INOUT) :: vxy(:,:,:,:)!(vacuum%nmzxyd,stars%ng2-1,2,:)
    INTEGER, INTENT (OUT):: nv2(DIMENSION%jspd)
    REAL,    INTENT (INOUT) :: vz(vacuum%nmzd,2,4)
    REAL,    INTENT (IN) :: evac(2,DIMENSION%jspd)
    REAL,    INTENT (IN) :: bkpt(3)

    LOGICAL,INTENT(IN)    :: l_real
    !     ..
    !     .. Local Scalars ..
    COMPLEX hij,sij,apw_lo,c_1
    REAL d2,gz,sign,th,wronk
    INTEGER i,i2,ii,jj,ik,j,jk,k,jspin,ipot,ii0,sb,i0
    INTEGER ivac,irec,imz,igvm2,igvm2i
    INTEGER jspin1,jspin2,jmax,jsp_start,jsp_end
    INTEGER i_start,nc,nc_0
    !     ..
    !     .. Local Arrays ..
    INTEGER kvac1(DIMENSION%nv2d,DIMENSION%jspd),kvac2(DIMENSION%nv2d,DIMENSION%jspd)
    INTEGER map2(DIMENSION%nvd,DIMENSION%jspd)
    COMPLEX tddv(DIMENSION%nv2d,DIMENSION%nv2d),tduv(DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX tudv(DIMENSION%nv2d,DIMENSION%nv2d),tuuv(DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX vxy_help(stars%ng2-1)
    COMPLEX a(DIMENSION%nvd,DIMENSION%jspd),b(DIMENSION%nvd,DIMENSION%jspd)
    REAL ddnv(DIMENSION%nv2d,DIMENSION%jspd),dudz(DIMENSION%nv2d,DIMENSION%jspd)
    REAL duz(DIMENSION%nv2d,DIMENSION%jspd), udz(DIMENSION%nv2d,DIMENSION%jspd)
    REAL uz(DIMENSION%nv2d,DIMENSION%jspd)
    ! l_J auxiliary potential array
    COMPLEX, ALLOCATABLE :: vxy1(:,:,:)
    !     ..

    d2 = SQRT(cell%omtil/cell%area)

    IF (jij%l_J) ALLOCATE (vxy1(vacuum%nmzxyd,stars%ng2-1,2))

    !--->    set up mapping function from 3d-->2d lapws

    DO jspin = 1,input%jspins
       nv2(jspin) = 0
       k_loop:DO  k = 1,lapw%nv(jspin)
          DO  j = 1,nv2(jspin)
             IF (lapw%k1(k,jspin).EQ.kvac1(j,jspin)&
                  .AND. lapw%k2(k,jspin).EQ.kvac2(j,jspin)) THEN
                map2(k,jspin) = j
                CYCLE k_loop
             END IF
          ENDDO
          nv2(jspin) = nv2(jspin) + 1
          IF (nv2(jspin)>DIMENSION%nv2d)  CALL juDFT_error("hsvac:dimension%nv2d",calledby ="hsvac")
          kvac1(nv2(jspin),jspin) = lapw%k1(k,jspin)
          kvac2(nv2(jspin),jspin) = lapw%k2(k,jspin)
          map2(k,jspin) = nv2(jspin)
       ENDDO k_loop
    ENDDO
    !--->    loop over the two vacuua (1: upper; 2: lower)
    DO ivac = 1,2
       sign = 3. - 2.*ivac !+/- 1
       DO sb = 1,MERGE(3,1,noco%l_noco) !loop over different (spin)-blocks of the potential
          !--->       get the wavefunctions and set up the tuuv, etc matrices          
          CALL lapw%spinblock(noco%l_noco,sb,jsp,jspin1,jspin2,ii,jj,ipot)
          jspin=jsp
          CALL vacfun(&
               vacuum,DIMENSION,stars,&
               jsp,input,noco,sb,&
               sym, cell,ivac,evac(1,1),bkpt,vxy(:,:,ivac,ipot),vz(:,:,:),kvac1,kvac2,nv2,&
               tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
          !
          !--->       generate a and b coeffficients
          !
          IF (noco%l_noco) THEN
             DO jspin = 1,input%jspins
                DO k = 1,lapw%nv(jspin)
                   gz = sign*cell%bmat(3,3)*lapw%k3(k,jspin)
                   i2 = map2(k,jspin)
                   th = gz*cell%z1
                   c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronk)
                   a(k,jspin) = - c_1 * CMPLX(dudz(i2,jspin), gz*udz(i2,jspin) )
                   b(k,jspin) =   c_1 * CMPLX(duz(i2,jspin), gz* uz(i2,jspin) )
                ENDDO
             ENDDO
          ELSE
             DO k = 1,lapw%nv(jsp)
                gz = sign*cell%bmat(3,3)*lapw%k3(k,jsp)
                i2 = map2(k,jsp)
                th = gz*cell%z1
                c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronk)
                a(k,1) = - c_1 * CMPLX(dudz(i2,jsp), gz*udz(i2,jsp) )
                b(k,1) =   c_1 * CMPLX(duz(i2,jsp), gz* uz(i2,jsp) )
             ENDDO
          ENDIF
          !--->       update hamiltonian and overlap matrices
          IF (sb<3) THEN
             DO  i0 = 1,smat%matsize2
                i=smat%local_rk_map(i0,jspin)
                ik = map2(i,jspin)
                DO j = 1,i - 1 !TODO check noco case
                   !--->             overlap: only  (g-g') parallel=0       '
                   IF (map2(j,jspin).EQ.ik) THEN
                      sij = CONJG(a(i,jspin))*a(j,jspin) + &
                           CONJG(b(i,jspin))*b(j,jspin)*ddnv(ik,jspin1)
                      !+APW_LO
                      IF (input%l_useapw) THEN
                         apw_lo = CONJG(a(i,jspin)*  uz(ik,jspin1) + b(i,jspin)* udz(ik,jspin1) ) &
                              * (a(j,jspin)* duz(ik,jspin1) + b(j,jspin)*dudz(ik,jspin1) )&
                              +      (a(j,jspin)*  uz(ik,jspin1) + b(j,jspin)* udz(ik,jspin1) ) &
                              * CONJG(a(i,jspin)* duz(ik,jspin1) + b(i,jspin)*dudz(ik,jspin1) )
                         !            IF (i.lt.10) write (3,'(2i4,2f20.10)') i,j,apw_lo
                         IF (l_real) THEN
                            hmat%data_r(jj+j,ii+i0) = hmat%data_r(jj+j,ii+i0) + 0.25 * REAL(apw_lo) 
                         ELSE 
                            hmat%data_c(jj+j,ii+i0) = hmat%data_c(jj+j,ii+i0) + 0.25 * apw_lo
                         ENDIF
                      ENDIF
                      !Overlapp Matrix
                      IF (l_real) THEN
                         smat%data_r(jj+j,ii+i0) = smat%data_r(jj+j,ii+i0) + REAL(sij)
                      ELSE 
                         smat%data_c(jj+j,ii+i0) = smat%data_c(jj+j,ii+i0) + sij
                      ENDIF
                   END IF
                ENDDO
                !Diagonal term of Overlapp matrix, Hamiltonian later
                sij = CONJG(a(i,jspin))*a(i,jspin) + CONJG(b(i,jspin))*b(i,jspin)*ddnv(ik,jspin1)
                IF (l_real) THEN
                   smat%data_r(jj+j,ii+i0) = smat%data_r(jj+j,ii+i0) + REAL(sij)
                ELSE
                   smat%data_c(jj+j,ii+i0) = smat%data_c(jj+j,ii+i0) + sij
                ENDIF
             ENDDO
          ENDIF

          !--->    hamiltonian update
          DO i0=1,hmat%matsize2
             i=hmat%local_rk_map(i0,sb)
             ik = map2(i,jspin1)
             DO j = 1,MERGE(lapw%nv(1),i,sb==3)
                jk = map2(j,jspin2)
                hij = CONJG(a(i,jspin1))* (tuuv(ik,jk)*a(j,jspin2) +tudv(ik,jk)*b(j,jspin2))&
                     + CONJG(b(i,jspin1))* (tddv(ik,jk)*b(j,jspin2) +tduv(ik,jk)*a(j,jspin2))
                IF (l_real) THEN
                   hmat%data_r(jj+j,ii+i0) = hmat%data_r(jj+j,ii+i0) + REAL(hij)
                ELSE
                   hmat%data_c(jj+j,ii+i0) = hmat%data_c(jj+j,ii+i0) + hij
                ENDIF
             ENDDO
          ENDDO

          !--->    end of loop over different parts of the potential matrix
       ENDDO

       !---> end of loop over vacua
    ENDDO


  END SUBROUTINE hsvac
END MODULE m_hsvac
