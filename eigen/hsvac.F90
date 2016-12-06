MODULE m_hsvac
  USE m_juDFT
CONTAINS
  SUBROUTINE hsvac(&
       vacuum,stars,DIMENSION, atoms, jsp,input,vxy,vz,evac,cell,&
       bkpt,lapw,sym, noco,jij, n_size,n_rank,nv2,l_real,hamOvlp)
    !*********************************************************************
    !     adds in the vacuum contributions to the the hamiltonian and
    !     overlap matrices. as written, each k-point calculates the
    !     vacuum functions again since there is only a single vacuum
    !     parameter per vacuum.
    !                m. weinert
    !*********************************************************************
    !     modified by R.Podloucky for speeding up and microtaskin
    !*********************************************************************

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
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp   ,n_size,n_rank
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (INOUT) :: vxy(vacuum%nmzxyd,stars%n2d-1,2)
    INTEGER, INTENT (OUT):: nv2(DIMENSION%jspd)
    REAL,    INTENT (INOUT) :: vz(vacuum%nmzd,2,4)
    REAL,    INTENT (IN) :: evac(2,DIMENSION%jspd)
    REAL,    INTENT (IN) :: bkpt(3)

    LOGICAL,INTENT(IN)    :: l_real
    !     ..
    !     .. Local Scalars ..
    COMPLEX hij,sij,apw_lo,c_1
    REAL d2,gz,sign,th,wronk,fac1
    INTEGER i,i2,ii,ik,j,jk,k,jspin,ipot,npot,ii0
    INTEGER ivac,irec,imz,igvm2,igvm2i
    INTEGER jspin1,jspin2,jmax,jsp_start,jsp_end
    INTEGER i_start,nc,nc_0
    !     ..
    !     .. Local Arrays ..
    INTEGER kvac1(DIMENSION%nv2d,DIMENSION%jspd),kvac2(DIMENSION%nv2d,DIMENSION%jspd)
    INTEGER map2(DIMENSION%nvd,DIMENSION%jspd)
    COMPLEX tddv(DIMENSION%nv2d,DIMENSION%nv2d),tduv(DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX tudv(DIMENSION%nv2d,DIMENSION%nv2d),tuuv(DIMENSION%nv2d,DIMENSION%nv2d)
    COMPLEX vxy_help(stars%n2d-1)
    COMPLEX a(DIMENSION%nvd,DIMENSION%jspd),b(DIMENSION%nvd,DIMENSION%jspd)
    REAL ddnv(DIMENSION%nv2d,DIMENSION%jspd),dudz(DIMENSION%nv2d,DIMENSION%jspd)
    REAL duz(DIMENSION%nv2d,DIMENSION%jspd), udz(DIMENSION%nv2d,DIMENSION%jspd)
    REAL uz(DIMENSION%nv2d,DIMENSION%jspd)
    ! l_J auxiliary potential array
    COMPLEX, ALLOCATABLE :: vxy1(:,:,:)
    !     ..

    d2 = SQRT(cell%omtil/cell%area)

    IF (jij%l_J) ALLOCATE (vxy1(vacuum%nmzxyd,stars%n2d-1,2))

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
       sign = 3. - 2.*ivac
       npot = 1
       !---> pk non-collinear
       IF (noco%l_noco) THEN
          !--->       if the two vacuua are equivalent, the potential file has to
          !--->       be backspaced, because the potential is the same at both
          !--->       surfaces of the film
          IF ((ivac.EQ.2) .AND. (vacuum%nvac.EQ.1)) THEN
             DO irec = 1,4
                BACKSPACE (25)
             ENDDO
          ENDIF
          !--->       load the non-warping part of the potential
          READ (25)((vz(imz,ivac,ipot),imz=1,vacuum%nmzd),ipot=1,4)
          npot = 3
          ! for J-coeff. we average the up-up and down-down parts and off-diagonal elements of the
          ! potential matrix to zero
          IF (jij%l_J) THEN
             vz(:,ivac,1) = (vz(:,ivac,1) + vz(:,ivac,2))/2.
             vz(:,ivac,2) =  vz(:,ivac,1)
             vz(:,ivac,3) = 0.0
             vz(:,ivac,4) = 0.0
          END IF
       ENDIF
       !---> pk non-collinear

       DO ipot = 1,npot
          !--->       get the wavefunctions and set up the tuuv, etc matrices
          IF (noco%l_noco) THEN
             IF (.NOT.jij%l_J) THEN
                READ (25)((vxy(imz,igvm2,ivac), imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
             END IF
             ! l_J we want to average the diagonal elements of the potential matrix
             IF (jij%l_J .AND. ipot.EQ.1) THEN
                READ (25)((vxy(imz,igvm2,ivac), imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
                READ (25)((vxy1(imz,igvm2,ivac), imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
                vxy(:,:,ivac) = (vxy(:,:,ivac)+vxy1(:,:,ivac))/2.
             END IF

             IF (jij%l_J .AND. ipot.EQ.3) THEN
                READ (25)((vxy(imz,igvm2,ivac), imz=1,vacuum%nmzxy),igvm2=1,stars%ng2-1)
             END IF

             IF (vacuum%nvac.EQ.1 .AND. ivac.EQ.2 .AND.(.NOT.sym%zrfs) ) THEN
                !--->          In this case (inversion, no z-reflection and thus no
                !--->          2d-inversion) the coeff. of a 2d-star (containing G) of
                !--->          the second vacuum is equal to the invers 2d-star
                !--->          (containing -G) of the first vacuum.
                ! l_J no need to do this symmetrizaion twice, if the potential vxy is the same
                IF (.NOT.jij%l_J .OR. (jij%l_J .AND. ipot.EQ.1)) THEN
                   DO imz = 1,vacuum%nmzxy
                      DO igvm2 = 2,stars%ng2
                         !--->                   find the index of the invers 2d-star
                         igvm2i = stars%ig2(stars%ig(-stars%kv2(1,igvm2),-stars%kv2(2,igvm2),0))
                         vxy_help(igvm2-1) = vxy(imz,igvm2i-1,2)
                      ENDDO
                      DO igvm2 = 2,stars%ng2
                         vxy(imz,igvm2-1,2) = vxy_help(igvm2-1)
                      ENDDO
                   ENDDO
                END IF ! jij%l_J
             ENDIF ! ivac-vacuum%nvac
             ! l_J we want the off-diagonal potential matrix elements to be zero
             IF (jij%l_J .AND. ipot.EQ.3) vxy(:,:,ivac)=CMPLX(0.,0.)
          ENDIF
          CALL vacfun(&
               vacuum,DIMENSION,stars,&
               jsp,input,noco,ipot,&
               sym, cell,ivac,evac(1,1),bkpt,vxy(1,1,ivac),vz,kvac1,kvac2,nv2,&
               tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv,wronk)
          fac1 = 1.0 / (d2*wronk)
          !
          !--->       generate a and b coeffficients
          !
          IF (noco%l_noco) THEN
             DO jspin = 1,input%jspins
                DO k = 1,lapw%nv(jspin)
                   gz = sign*cell%bmat(3,3)*lapw%k3(k,jspin)
                   i2 = map2(k,jspin)
                   th = gz*cell%z1
                   c_1 = fac1 * CMPLX( COS(th), SIN(th) )
                   a(k,jspin) = - c_1 * CMPLX(dudz(i2,jspin), gz*udz(i2,jspin) )
                   b(k,jspin) =   c_1 * CMPLX(duz(i2,jspin), gz* uz(i2,jspin) )
                ENDDO
             ENDDO
          ELSE
             DO k = 1,lapw%nv(jsp)
                gz = sign*cell%bmat(3,3)*lapw%k3(k,jsp)
                i2 = map2(k,jsp)
                th = gz*cell%z1
                c_1 = fac1 * CMPLX( COS(th), SIN(th) )
                a(k,1) = - c_1 * CMPLX(dudz(i2,jsp), gz*udz(i2,jsp) )
                b(k,1) =   c_1 * CMPLX(duz(i2,jsp), gz* uz(i2,jsp) )
             ENDDO
          ENDIF
          !--->       update hamiltonian and overlap matrices
          IF (ipot.EQ.1 .OR. ipot.EQ.2) THEN
             jspin = ipot
             !+gb||
             IF (ipot.EQ.1) THEN
                nc = 0
                i_start = n_rank
             ELSE
                nc = nc + atoms%nlotot
                nc_0 = nc
                i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size) 
             ENDIF
             !-gb||
             DO  i = i_start+1, lapw%nv(jspin), n_size
                ik = map2(i,jspin)
                nc = nc + 1
                IF (ipot.EQ.1) THEN
                   jspin = 1
                   ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                ELSEIF (ipot.EQ.2) THEN
                   jspin = 2
                   ii0=nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)+ lapw%nv(1)+atoms%nlotot
                ENDIF
                jspin1 = jsp
                IF (noco%l_noco) jspin1 = jspin
                DO j = 1,i - 1
                   ii = ii0 + j
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
                            hamOvlp%a_r(ii) = hamOvlp%a_r(ii) + 0.25 * REAL(apw_lo) 
                         ELSE 
                            hamOvlp%a_c(ii) = hamOvlp%a_c(ii) + 0.25 * apw_lo
                         ENDIF
                      ENDIF
                      !+APW_LO
                      IF (l_real) THEN
                         hamOvlp%b_r(ii) = hamOvlp%b_r(ii) + REAL(sij)
                      ELSE 
                         hamOvlp%b_c(ii) = hamOvlp%b_c(ii) + sij
                      ENDIF
                   END IF
                ENDDO
                ii = ii0 + i
                sij = CONJG(a(i,jspin))*a(i,jspin) + CONJG(b(i,jspin))*b(i,jspin)*ddnv(ik,jspin1)
                IF (l_real) THEN

                   hamOvlp%b_r(ii) = hamOvlp%b_r(ii) + REAL(sij)
                ELSE
                   hamOvlp%b_c(ii) = hamOvlp%b_c(ii) + sij
                ENDIF
             ENDDO
          ENDIF

          !--->    hamiltonian update
          IF (ipot.EQ.1) THEN
             jspin1 = 1
             jspin2 = 1
             nc = 0
             i_start = n_rank
          ELSEIF (ipot.EQ.2) THEN
             jspin1 = 2
             jspin2 = 2
             nc = nc_0
             i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size) 
          ELSEIF (ipot.EQ.3) THEN
             jspin1 = 2
             jspin2 = 1
             nc = nc_0
             i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size) 
          ENDIF
          DO i = i_start+1, lapw%nv(jspin1), n_size
             ik = map2(i,jspin1)
             nc = nc + 1
             IF (ipot.EQ.1) THEN
                ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                jmax = i
             ELSEIF (ipot.EQ.2) THEN
                ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1) + lapw%nv(1)+atoms%nlotot
                jmax = i
             ELSEIF (ipot.EQ.3) THEN
                ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                jmax = lapw%nv(jspin2)
             ENDIF
             DO j = 1,jmax
                ii = ii0 + j
                jk = map2(j,jspin2)
                hij = CONJG(a(i,jspin1))* (tuuv(ik,jk)*a(j,jspin2) +tudv(ik,jk)*b(j,jspin2))&
                     + CONJG(b(i,jspin1))* (tddv(ik,jk)*b(j,jspin2) +tduv(ik,jk)*a(j,jspin2))
                IF (l_real) THEN

                   hamOvlp%a_r(ii) = hamOvlp%a_r(ii) + REAL(hij)
                ELSE
                   hamOvlp%a_c(ii) = hamOvlp%a_c(ii) + hij
                ENDIF
             ENDDO
          ENDDO

          !--->    end of loop over different parts of the potential matrix
       ENDDO

       !---> end of loop over vacua
    ENDDO

    IF (jij%l_J) DEALLOCATE (vxy1)

  END SUBROUTINE hsvac
END MODULE m_hsvac
