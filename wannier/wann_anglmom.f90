!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wann_anglmom
  !***********************************************************************
  !     Compute matrix elements of angular momentum operator 
  !     in the muffin-tin spheres.
  !
  !     Frank Freimuth
  !***********************************************************************
CONTAINS
  SUBROUTINE wann_anglmom(atoms,usdus,jspin,acof,bcof,ccof, mmn)
    USE m_types
    IMPLICIT NONE
    !     .. scalar arguments ..
    TYPE(t_atoms),INTENT(in)::atoms
    TYPE(t_usdus),INTENT(in)::usdus
    INTEGER,INTENT(IN)      ::jspin
    !     .. array arguments ..
    COMPLEX, INTENT (in)  :: ccof(-atoms%llod:,:,:,:) !ccof(-llod:llod,noccbd,atoms%nlod,natd)
    COMPLEX, INTENT (in)  :: acof(:,0:,:)!acof(noccbd,0:lmd,natd)
    COMPLEX, INTENT (in)  :: bcof(:,0:,:)!bcof(noccbd,0:lmd,natd)
    COMPLEX, INTENT (inout) :: mmn(:,:,:)!mmn(3,noccbd,noccbd)
    !     .. local scalars ..
    LOGICAL :: l_select
    INTEGER :: i,j,l,lo,lop,m,natom,nn,ntyp
    INTEGER :: nt1,nt2,lm,n,ll1,indat
    COMPLEX :: suma_z,sumb_z
    COMPLEX :: suma_p,sumb_p
    COMPLEX :: suma_m,sumb_m
    COMPLEX :: suma_x,sumb_x
    COMPLEX :: suma_y,sumb_y
    REAL    :: lplus,lminus
    !     ..
    !     .. local arrays ..
    COMPLEX, ALLOCATABLE :: qlo_z(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: qlo_p(:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: qlo_m(:,:,:,:,:)

    COMPLEX, ALLOCATABLE :: qaclo_z(:,:,:,:),qbclo_z(:,:,:,:)
    COMPLEX, ALLOCATABLE :: qaclo_p(:,:,:,:),qbclo_p(:,:,:,:)
    COMPLEX, ALLOCATABLE :: qaclo_m(:,:,:,:),qbclo_m(:,:,:,:)
    !     ..
    !     .. intrinsic functions ..
    INTRINSIC conjg

    ALLOCATE (qlo_z(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%nlod,atoms%ntype) &
           ,qaclo_z(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype),&
           qbclo_z(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype) )

    ALLOCATE (qlo_p(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%nlod,atoms%ntype) &
            ,qaclo_p(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype),&
            qbclo_p(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype) )

    ALLOCATE (qlo_m(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%nlod,atoms%ntype)&
               ,qaclo_m(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype),&
                qbclo_m(SIZE(acof,1),SIZE(acof,1),atoms%nlod,atoms%ntype) )

    INQUIRE(file='select_anglmom',exist=l_select)
    WRITE(*,*)'select_anglmom: ',l_select
    IF(l_select) THEN
       OPEN(866,file='select_anglmom')
       READ(866,*)indat
       CLOSE(866)
       WRITE(*,*)'anglmom for atom=',indat
       WRITE(*,*)atoms%ntype
       WRITE(*,*)atoms%neq(indat)
    ENDIF

    !-----> lapw-lapw-Terms
    DO i = 1,SIZE(acof,1)            
       DO j = 1,SIZE(acof,1)
          nt1 = 1
          DO n = 1,atoms%ntype
             nt2 = nt1 + atoms%neq(n) - 1
             DO l = 0,atoms%lmax(n)
                suma_z = CMPLX(0.,0.); sumb_z = CMPLX(0.,0.)
                suma_m = CMPLX(0.,0.); sumb_m = CMPLX(0.,0.)
                suma_p = CMPLX(0.,0.); sumb_p = CMPLX(0.,0.)
                IF(l_select .AND. (n.NE.indat)) CYCLE
                ll1 = l* (l+1)
                DO m = -l,l
                   lm = ll1 + m
                   lplus=SQRT(REAL( (l-m)*(l+m+1) ) )
                   lminus=SQRT(REAL( (l+m)*(l-m+1) ) )
                   DO natom = nt1,nt2
                      suma_z = suma_z + acof(i,lm,natom)*&
                                                     CONJG(acof(j,lm,natom))*REAL(m)
                      sumb_z = sumb_z + bcof(i,lm,natom)*&
                                                     CONJG(bcof(j,lm,natom))*REAL(m)
                      IF(m+1.LE.l)THEN
                         suma_p = suma_p + acof(i,lm,natom)*&
                                                        CONJG(acof(j,lm+1,natom))*lplus
                         sumb_p = sumb_p + bcof(i,lm,natom)*&
                                                        CONJG(bcof(j,lm+1,natom))*lplus
                      ENDIF
                      IF(m-1.GE.-l)THEN
                         suma_m = suma_m + acof(i,lm,natom)*&
                                                        CONJG(acof(j,lm-1,natom))*lminus
                         sumb_m = sumb_m + bcof(i,lm,natom)*&
                                                        CONJG(bcof(j,lm-1,natom))*lminus
                      ENDIF
                   ENDDO
                ENDDO
                mmn(3,j,i) = mmn(3,j,i) + (suma_z+sumb_z*usdus%ddn(l,n,jspin))

                suma_x=0.5*(suma_p+suma_m)
                sumb_x=0.5*(sumb_p+sumb_m)
                mmn(1,j,i) = mmn(1,j,i) + (suma_x+sumb_x*usdus%ddn(l,n,jspin))

                suma_y=CMPLX(0.0,-0.5)*(suma_p-suma_m)
                sumb_y=CMPLX(0.0,-0.5)*(sumb_p-sumb_m)
                mmn(2,j,i) = mmn(2,j,i) + (suma_y+sumb_y*usdus%ddn(l,n,jspin))
             ENDDO ! l
             nt1 = nt1 + atoms%neq(n)
          ENDDO ! n
       ENDDO ! j
    ENDDO ! i


    !---> Terms involving local orbitals.
    qlo_z = 0.0; qlo_p = 0.0; qlo_m = 0.0
    qaclo_z = 0.0; qaclo_p = 0.0; qaclo_m = 0.0
    qbclo_z = 0.0; qbclo_p = 0.0; qbclo_m = 0.0

    natom = 0
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          IF(l_select .AND. (ntyp.NE.indat)) CYCLE
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                lplus=SQRT(REAL( (l-m)*(l+m+1) ) )
                lminus=SQRT(REAL( (l+m)*(l-m+1) ) )
                DO i = 1,SIZE(acof,1)
                   DO j = 1,SIZE(acof,1)
                      qbclo_z(j,i,lo,ntyp) = qbclo_z(j,i,lo,ntyp) + (&
                                    bcof(i,lm,natom) * CONJG(ccof(m,j,lo,natom)) +&
                                    ccof(m,i,lo,natom)*CONJG(bcof(j,lm,natom)) )*REAL(m)

                      qaclo_z(j,i,lo,ntyp) = qaclo_z(j,i,lo,ntyp) + (&
                                    acof(i,lm,natom) * CONJG(ccof(m,j,lo,natom)) +&
                                    ccof(m,i,lo,natom)*CONJG(acof(j,lm,natom)) )*REAL(m)
                      IF(m+1.LE.l)THEN
                         qbclo_p(j,i,lo,ntyp) = qbclo_p(j,i,lo,ntyp) + (&
                                         bcof(i,lm,natom) * CONJG(ccof(m+1,j,lo,natom)) +&
                                         ccof(m,i,lo,natom)*CONJG(bcof(j,lm+1,natom)) )*lplus

                         qaclo_p(j,i,lo,ntyp) = qaclo_p(j,i,lo,ntyp) + (&
                                         acof(i,lm,natom) * CONJG(ccof(m+1,j,lo,natom)) +&
                                         ccof(m,i,lo,natom)*CONJG(acof(j,lm+1,natom)) )*lplus
                      ENDIF
                      IF(m-1.GE.-l)THEN
                         qbclo_m(j,i,lo,ntyp) = qbclo_m(j,i,lo,ntyp) + (&
                                         bcof(i,lm,natom) * CONJG(ccof(m-1,j,lo,natom)) +&
                                         ccof(m,i,lo,natom)*CONJG(bcof(j,lm-1,natom)) )*lminus

                         qaclo_m(j,i,lo,ntyp) = qaclo_m(j,i,lo,ntyp) + (&
                                         acof(i,lm,natom) * CONJG(ccof(m-1,j,lo,natom)) +&
                                         ccof(m,i,lo,natom)*CONJG(acof(j,lm-1,natom)) )*lminus
                      ENDIF

                   ENDDO !j
                ENDDO !i
             ENDDO !m
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      lplus=SQRT(REAL( (l-m)*(l+m+1) ) )
                      lminus=SQRT(REAL( (l+m)*(l-m+1) ) )
                      DO i = 1,SIZE(acof,1)
                         DO j = 1,SIZE(acof,1)
                            qlo_z(j,i,lop,lo,ntyp) = qlo_z(j,i,lop,lo,ntyp) + &
                                                     CONJG(ccof(m,j,lop,natom))&
                                                                *ccof(m,i,lo,natom)*REAL(m)
                            IF(m+1.LE.l)THEN
                               qlo_p(j,i,lop,lo,ntyp) = &
                                                       qlo_p(j,i,lop,lo,ntyp) + &
                                                        CONJG(ccof(m+1,j,lop,natom))&
                                                             *ccof(m,i,lo,natom)*lplus

                            ENDIF
                            IF(m-1.GE.-l)THEN
                               qlo_m(j,i,lop,lo,ntyp) = &
                                                       qlo_m(j,i,lop,lo,ntyp) + &
                                                        CONJG(ccof(m-1,j,lop,natom))&
                                                             *ccof(m,i,lo,natom)*lminus
                            ENDIF
                         ENDDO ! j
                      ENDDO ! i
                   ENDDO ! m
                ENDIF
             ENDDO ! lop
          ENDDO ! lo
       ENDDO ! nn
    ENDDO ! ntyp
    !---> perform summation of the coefficients with the integrals
    !---> of the radial basis functions
    DO ntyp = 1,atoms%ntype
       IF(l_select .AND. (ntyp.NE.indat) ) CYCLE
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          DO j = 1,SIZE(acof,1)
             DO i = 1,SIZE(acof,1)
                mmn(3,i,j)= mmn(3,i,j)  + &
                                    qaclo_z(i,j,lo,ntyp)*usdus%uulon(lo,ntyp,jspin) +&
                                    qbclo_z(i,j,lo,ntyp)*usdus%dulon(lo,ntyp,jspin)  

                suma_p=qaclo_p(i,j,lo,ntyp)*usdus%uulon(lo,ntyp,jspin) +&
                                     qbclo_p(i,j,lo,ntyp)*usdus%dulon(lo,ntyp,jspin)

                suma_m=qaclo_m(i,j,lo,ntyp)*usdus%uulon(lo,ntyp,jspin) +&
                                     qbclo_m(i,j,lo,ntyp)*usdus%dulon(lo,ntyp,jspin)

                suma_x=            0.5*(suma_p+suma_m)
                suma_y=CMPLX(0.0,-0.5)*(suma_p-suma_m)

                mmn(1,i,j)= mmn(1,i,j)  + suma_x
                mmn(2,i,j)= mmn(2,i,j)  + suma_y 

             ENDDO !i
          ENDDO !j 
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO j = 1,SIZE(acof,1)
                   DO i = 1,SIZE(acof,1)
                      mmn(3,i,j) = mmn(3,i,j)  + &
                                          qlo_z(i,j,lop,lo,ntyp)*usdus%uloulopn(lop,lo,ntyp,jspin)
                      suma_p=qlo_p(i,j,lop,lo,ntyp)*usdus%uloulopn(lop,lo,ntyp,jspin)
                      suma_m=qlo_m(i,j,lop,lo,ntyp)*usdus%uloulopn(lop,lo,ntyp,jspin)
                      mmn(1,i,j) = mmn(1,i,j) + 0.5*(suma_p+suma_m)
                      mmn(2,i,j) = mmn(2,i,j) + &
                                            CMPLX(0.0,-0.5)*(suma_p-suma_m)
                   ENDDO ! i
                ENDDO ! j
             ENDIF
          ENDDO !lop
       ENDDO !lo 
    ENDDO !ntyp 
    DEALLOCATE ( qlo_z,qaclo_z,qbclo_z )
    DEALLOCATE ( qlo_m,qaclo_m,qbclo_m )
    DEALLOCATE ( qlo_p,qaclo_p,qbclo_p )

  END SUBROUTINE wann_anglmom
END MODULE m_wann_anglmom
