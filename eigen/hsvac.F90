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
       vacuum,stars, mpi,jsp,input,v,evac,cell,&
       lapw,sym, noco,hmat,smat)
 

    USE m_vacfun
    USE m_types
    IMPLICIT NONE
    
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_vacuum),INTENT(IN)     :: vacuum
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_potden),INTENT(IN)     :: v
    CLASS(t_mat),INTENT(INOUT)    :: hmat(:,:),smat(:,:)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: evac(2,input%jspins)
    !     ..
    !     .. Local Scalars ..
    COMPLEX hij,sij,apw_lo,c_1
    REAL d2,gz,sign,th,wronk
    INTEGER i,i2,ii,jj,ik,j,jk,k,jspin,ii0,i0
    INTEGER ivac,irec,imz,igvm2,igvm2i,s1,s2
    INTEGER jspin1,jspin2,jmax,jsp_start,jsp_end
    INTEGER i_start,nc,nc_0
    !     ..
    !     .. Local Arrays ..
    INTEGER:: nv2(input%jspins)
    INTEGER kvac(2,lapw%dim_nv2d(),input%jspins)
    INTEGER map2(lapw%dim_nvd(),input%jspins)
    COMPLEX tddv(lapw%dim_nv2d(),lapw%dim_nv2d()),tduv(lapw%dim_nv2d(),lapw%dim_nv2d())
    COMPLEX tudv(lapw%dim_nv2d(),lapw%dim_nv2d()),tuuv(lapw%dim_nv2d(),lapw%dim_nv2d())
    COMPLEX vxy_help(stars%ng2-1)
    COMPLEX a(lapw%dim_nvd(),input%jspins),b(lapw%dim_nvd(),input%jspins)
    REAL ddnv(lapw%dim_nv2d(),input%jspins),dudz(lapw%dim_nv2d(),input%jspins)
    REAL duz(lapw%dim_nv2d(),input%jspins), udz(lapw%dim_nv2d(),input%jspins)
    REAL uz(lapw%dim_nv2d(),input%jspins)
    !     ..


    d2 = SQRT(cell%omtil/cell%area)

    !--->    set up mapping function from 3d-->2d lapws

    DO jspin = 1,input%jspins
       nv2(jspin) = 0
       k_loop:DO  k = 1,lapw%nv(jspin)
          DO  j = 1,nv2(jspin)
             IF (all(lapw%gvec(1:2,k,jspin)==kvac(1:2,j,jspin))) THEN
                map2(k,jspin) = j
                CYCLE k_loop
             END IF
          ENDDO
          nv2(jspin) = nv2(jspin) + 1
          IF (nv2(jspin)>lapw%dim_nv2d())  CALL juDFT_error("hsvac:lapw%dim_nv2d()",calledby ="hsvac")
          kvac(1:2,nv2(jspin),jspin) = lapw%gvec(1:2,k,jspin)
          map2(k,jspin) = nv2(jspin)
       ENDDO k_loop
    ENDDO
    !--->    loop over the two vacuua (1: upper; 2: lower)
    DO ivac = 1,2
       sign = 3. - 2.*ivac !+/- 1
       DO jspin1=MERGE(1,jsp,noco%l_noco),MERGE(2,jsp,noco%l_noco) !loop over global spin
          s1=MIN(SIZE(hmat,1),jspin1) !in colinear case s1=1
          DO jspin2=MERGE(1,jsp,noco%l_noco),MERGE(2,jsp,noco%l_noco) !loop over global spin
             s2=MIN(SIZE(hmat,1),jspin2) !in colinear case s2=1
          !--->       get the wavefunctions and set up the tuuv, etc matrices          
             CALL vacfun(&
                  vacuum,stars,&
                  input,noco,jspin1,jspin2,&
                  sym, cell,ivac,evac,lapw%bkpt,v%vacxy,v%vacz,kvac,nv2,&
                  tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
          !
          !--->       generate a and b coeffficients
          !
             DO jspin = MIN(jspin1,jspin2),MAX(jspin1,jspin2)
                DO k = 1,lapw%nv(jspin)
                   gz = sign*cell%bmat(3,3)*lapw%k3(k,jspin)
                   i2 = map2(k,jspin)
                   th = gz*cell%z1
                   c_1 = CMPLX( COS(th), SIN(th) )/ (d2*wronk)
                   a(k,jspin) = - c_1 * CMPLX(dudz(i2,jspin), gz*udz(i2,jspin) )
                   b(k,jspin) =   c_1 * CMPLX(duz(i2,jspin), gz* uz(i2,jspin) )
                ENDDO
             ENDDO
          !--->       update hamiltonian and overlap matrices
          IF (jspin1==jspin2) THEN
             DO  i = mpi%n_rank+1,lapw%nv(jspin2),mpi%n_size
                i0=(i-1)/mpi%n_size+1 !local column index
                ik = map2(i,jspin2)
                DO j = 1,i - 1 !TODO check noco case
                   !--->             overlap: only  (g-g') parallel=0       '
                   IF (map2(j,jspin1).EQ.ik) THEN
                      sij = CONJG(a(i,jspin2))*a(j,jspin2) + &
                           CONJG(b(i,jspin2))*b(j,jspin2)*ddnv(ik,jspin2)
                      !+APW_LO
                      IF (input%l_useapw) THEN
                         apw_lo = CONJG(a(i,jspin1)*  uz(ik,jspin1) + b(i,jspin1)* udz(ik,jspin1) ) &
                              * (a(j,jspin1)* duz(ik,jspin1) + b(j,jspin1)*dudz(ik,jspin1) )&
                              +      (a(j,jspin1)*  uz(ik,jspin1) + b(j,jspin1)* udz(ik,jspin1) ) &
                              * CONJG(a(i,jspin1)* duz(ik,jspin1) + b(i,jspin1)*dudz(ik,jspin1) )
                         !            IF (i.lt.10) write (3,'(2i4,2f20.10)') i,j,apw_lo
                         IF (hmat(1,1)%l_real) THEN
                            hmat(s1,s2)%data_r(j,i0) = hmat(s1,s2)%data_r(j,i0) + 0.25 * REAL(apw_lo) 
                         ELSE 
                            hmat(s1,s2)%data_c(j,i0) = hmat(s1,s2)%data_c(j,i0) + 0.25 * apw_lo
                         ENDIF
                      ENDIF
                      !Overlapp Matrix
                      IF (hmat(1,1)%l_real) THEN
                         smat(s1,s2)%data_r(j,i0) = smat(s1,s2)%data_r(j,i0) + REAL(sij)
                      ELSE 
                         smat(s1,s2)%data_c(j,i0) = smat(s1,s2)%data_c(j,i0) + sij
                      ENDIF
                   END IF
                ENDDO
                !Diagonal term of Overlapp matrix, Hamiltonian later
                sij = CONJG(a(i,jspin2))*a(i,jspin2) + CONJG(b(i,jspin2))*b(i,jspin2)*ddnv(ik,jspin2)
                IF (hmat(1,1)%l_real) THEN
                   smat(s2,s1)%data_r(j,i0) = smat(s1,s2)%data_r(j,i0) + REAL(sij)
                ELSE
                   smat(s2,s1)%data_c(j,i0) = smat(s1,s2)%data_c(j,i0) + sij
                ENDIF
             ENDDO
          ENDIF

          !--->    hamiltonian update
          DO  i = mpi%n_rank+1,lapw%nv(jspin1),mpi%n_size
             i0=(i-1)/mpi%n_size+1 !local column index
             ik = map2(i,jspin1)
             DO j = 1,MERGE(i,lapw%nv(jspin2),jspin1==jspin2)
                jk = map2(j,jspin2)
                hij = CONJG(a(i,jspin1))* (tuuv(ik,jk)*a(j,jspin2) +tudv(ik,jk)*b(j,jspin2))&
                     + CONJG(b(i,jspin1))* (tddv(ik,jk)*b(j,jspin2) +tduv(ik,jk)*a(j,jspin2))
                IF (hmat(1,1)%l_real) THEN
                   hmat(s2,s1)%data_r(j,i0) = hmat(s2,s1)%data_r(j,i0) + REAL(hij)
                ELSE
                   hmat(s2,s1)%data_c(j,i0) = hmat(s2,s1)%data_c(j,i0) + hij
                ENDIF
             ENDDO
          ENDDO

          !--->    end of loop over different parts of the potential matrix
       ENDDO

       !---> end of loop over vacua
    ENDDO
 ENDDO

  END SUBROUTINE hsvac
END MODULE m_hsvac
