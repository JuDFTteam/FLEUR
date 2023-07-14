!--------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP not_used
#else
#define CPP_OMP $OMP
#endif
MODULE m_hlomat
   use m_matmul_dgemm
  IMPLICIT NONE
  !***********************************************************************
  ! updates the hamiltonian  matrix with the contributions from the local
  ! orbitals.
  ! p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE hlomat(input,atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
       ntyp,na,fjgj,alo1,blo1,clo1, igSpinPr,igSpin,chi,hmat,l_fullj,l_ham,lapwq,fjgjq)

    USE m_hsmt_ab
    USE m_types
!    USE m_types_mpimat
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN),TARGET   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: fmpi
    TYPE(t_usdus),INTENT(IN)  :: ud
    TYPE(t_tlmplm),INTENT(IN) :: tlmplm
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_nococonv),INTENT(IN),TARGET   :: nococonv
    TYPE(t_fjgj),INTENT(IN),TARGET   :: fjgj


    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp
    INTEGER, INTENT (IN) :: ilSpinPr,ilSpin !spin for usdus and tlmplm
    INTEGER, INTENT (IN) :: igSpin,igSpinPr
    COMPLEX, INTENT (IN) :: chi
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN) :: alo1(:,:),blo1(:,:),clo1(:,:)

    CLASS(t_mat),INTENT (INOUT) :: hmat
    LOGICAL, INTENT(IN) :: l_fullj, l_ham

    TYPE(t_lapw), OPTIONAL, INTENT(IN),TARGET :: lapwq
    TYPE(t_fjgj), OPTIONAL, INTENT(IN),TARGET :: fjgjq
    !     ..
    ! Local Scalars
      COMPLEX :: prod,axx,bxx,cxx
      INTEGER :: invsfct,l,lm,lmp,lo,lolo,lolop,lop,lp,i,lo_lmax
      INTEGER :: mp,nkvec,nkvecp,lmplm,loplo,kp,m,mlo,mlolo,mlolo_new,lolop_new
      INTEGER :: locol,lorow,n,k,ab_size,ab_size_Pr,s
      LOGICAL :: l_samelapw

      ! Local Arrays
      COMPLEX, ALLOCATABLE :: abCoeffs(:,:), ax(:,:), bx(:,:), cx(:,:)
      COMPLEX, ALLOCATABLE :: abclo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: abCoeffsPr(:,:), axPr(:,:), bxPr(:,:), cxPr(:,:)
      COMPLEX, ALLOCATABLE :: abcloPr(:,:,:,:)

      TYPE(t_lapw) ,POINTER:: lapwPr
      TYPE(t_fjgj) ,POINTER:: fjgjPr

      l_samelapw = .FALSE.
      IF (.NOT.PRESENT(lapwq)) l_samelapw = .TRUE.
      IF (.NOT.l_samelapw) THEN
         lapwPr => lapwq
         fjgjPr => fjgjq
      ELSE
         lapwPr => lapw
         fjgjPr => fjgj
      END IF

      ! Synthesize a and b
      
      lo_lmax=maxval(atoms%llo)

      ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      ALLOCATE(ax(2*lo_lmax+1,MAXVAL(lapw%nv)),bx(2*lo_lmax+1,MAXVAL(lapw%nv)),cx(2*lo_lmax+1,MAXVAL(lapw%nv)))      
      ALLOCATE(abCoeffsPr(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapwPr%nv)))
      ALLOCATE(axPr(MAXVAL(lapwPr%nv),2*lo_lmax+1),bxPr(MAXVAL(lapwPr%nv),2*lo_lmax+1),cxPr(MAXVAL(lapwPr%nv),2*lo_lmax+1))
      ALLOCATE(abcloPr(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))

      !$acc data create(abcoeffs,abclo,abcoeffsPr,abcloPr)
      !$acc data copyin(alo1,blo1,clo1,fjgjPr,fjgjpr%fj,fjgjpr%gj)
      CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpinPr,igSpinPr,ntyp,na,cell,lapwPr,fjgjPr,abCoeffsPr(:,:),ab_size_Pr,.TRUE.,abcloPr,alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr))

      !we need the "unprimed" abcoeffs
      IF (ilSpin==ilSpinPr.AND.igSpinPr==igSpin.AND.l_samelapw) THEN
         !$acc kernels present(abcoeffs,abcoeffsPr,abclo,abcloPr)
         if (l_fullj) abcoeffs=abcoeffsPr !TODO automatic alloc on GPU????
         abclo=abcloPr
         !$acc end kernels
      ELSE     
         ALLOCATE(abCoeffs(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapw%nv)))       
         CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,ntyp,na,cell,lapw,fjgj,abCoeffs(:,:),ab_size,.TRUE.,abclo,alo1(:,ilSpin),blo1(:,ilSpin),clo1(:,ilSpin))
      END IF
   

      mlo=0;mlolo=0;mlolo_new=0
      DO m=1,ntyp-1
         mlo=mlo+atoms%nlo(m)
         mlolo=mlolo+atoms%nlo(m)*(atoms%nlo(m)+1)/2
         mlolo_new=mlolo_new+atoms%nlo(m)**2
      END DO

      IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
         ! If this atom is the first of two atoms related by inversion, the
         ! contributions to the overlap matrix of both atoms are added simultaneously.
         ! Where it is made use of the fact, that the sum of these contributions
         ! is twice the real part of the contribution of each atom. Note, that
         ! in this case there are twice as many (2*(2*l+1)) k-vectors.
         ! (compare abccoflo and comments there).
         IF (sym%invsat(na) == 0) invsfct = 1
         IF (sym%invsat(na) == 1) invsfct = 2
         CALL timestart("LAPW-LO")
         ! Calculate the hamiltonian matrix elements with the regular
         ! LAPW basis-functions   
         DO lo = 1,atoms%nlo(ntyp)
            l = atoms%llo(lo,ntyp)
            s = tlmplm%h_loc2_nonsph(ntyp) 
            !$acc data create(axpr,bxpr,cxpr)
            call blas_matmul(maxval(lapw%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_loc_LO(0:2*s-1,l*l:,ntyp,ilSpinPr,ilSpin),axPr,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
            call blas_matmul(maxval(lapw%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_loc_LO(0:2*s-1,s+l*l:,ntyp,ilSpinPr,ilSpin),bxPr,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
            call blas_matmul(maxval(lapw%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_LO(0:2*s-1,-l:,lo+mlo,ilSpinPr,ilSpin),cxPr,cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
       
            !Calculate the a,b,c coefs
            !DO m = -l,l
            !   lm = l* (l+1) + m
            !   s = tlmplm%h_loc2_nonsph(ntyp) 
            !   axPr(:,l+1+m) = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_loc_LO(0:2*s-1,lm,ntyp,ilSpinPr,ilSpin))
            !   bxPr(:,l+1+m) = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_loc_LO(0:2*s-1,s+lm,ntyp,ilSpinPr,ilSpin))
            !   cxPr(:,l+1+m) = matmul(transpose(conjg(abCoeffsPr(0:2*s-1,:))),tlmplm%h_LO(0:2*s-1,m,lo+mlo,ilSpinPr,ilSpin))
            !ENDDO  

          


            !LAPW LO contributions
            !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abclo,axpr,bxpr,cxpr)&
            !$acc & copyin(lapw,lapw%nv,lapw%index_lo,fmpi,fmpi%n_size,fmpi%n_rank)&
            !$acc & default(none)
            DO nkvec = 1,invsfct*(2*l+1)
               locol= lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
               IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
                  locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
                  IF (hmat%l_real) THEN
                    DO m=-l,l
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           hmat%data_r(kp,locol) = hmat%data_r(kp,locol) &
                                               & + REAL(chi) * invsfct * (&
                                               & abclo(1,m,nkvec,lo) *  axPr(kp,l+1+m) + &
                                               & abclo(2,m,nkvec,lo) *  bxPr(kp,l+1+m) + &
                                               & abclo(3,m,nkvec,lo) *  cxPr(kp,l+1+m) )
                        ENDDO   
                    END DO
                  ELSE
                    DO m=-l,l
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           hmat%data_c(kp,locol) = hmat%data_c(kp,locol) &
                                               & + chi * invsfct * ( &
                                               & abclo(1,m,nkvec,lo) *  axPr(kp,(l+1+m)) + &
                                               & abclo(2,m,nkvec,lo) *  bxPr(kp,(l+1+m)) + &
                                               & abclo(3,m,nkvec,lo) *  cxPr(kp,(l+1+m)) )
                        ENDDO   
                     END DO
                  END IF
                     ! Jump to the last matrix element of the current row
               END IF
            END DO
            !$acc end kernels
            !$acc end data
         ENDDO
         CALL timestop("LAPW-LO")
         IF (l_fullj) THEN
            CALL timestart("LO-LAPW")
            DO lo = 1,atoms%nlo(ntyp)
               l = atoms%llo(lo,ntyp)
               s = tlmplm%h_loc2_nonsph(ntyp) 
               !$acc data create(axpr,bxpr,cxpr)
               call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_loc_LO(0:2*s-1,l*l:,ntyp,ilSpinPr,ilSpin),abCoeffs,ax,cmplx(1.0,0.0),cmplx(0.0,0.0),'N','T')
               call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_loc_LO(0:2*s-1,s+l*l:,ntyp,ilSpinPr,ilSpin),abCoeffs,bx,cmplx(1.0,0.0),cmplx(0.0,0.0),'N','T')
               call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_LO2(0:2*s-1,-l:,lo+mlo,ilSpinPr,ilSpin),abCoeffs,cx,cmplx(1.0,0.0),cmplx(0.0,0.0),'N','T')

               !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abcloPr,ax,bx,cx)&
               !$acc & copyin(lapwPr,lapwPr%nv,lapwPr%index_lo,fmpi,fmpi%n_size,fmpi%n_rank)&
               !$acc & default(none)
               DO nkvec = 1,invsfct*(2*l+1)
                  lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lo,na)+nkvec
                  IF (hmat%l_real) THEN
                     DO m=-l,l
                        DO kp = 1,lapw%nv(igSpin)
                           hmat%data_r(lorow,kp) = hmat%data_r(lorow,kp) &
                                               & + REAL(chi) * invsfct * (&
                                               & CONJG(abcloPr(1,m,nkvec,lo)) *  ax(l+1+m,kp) + &
                                               & CONJG(abcloPr(2,m,nkvec,lo)) *  bx(l+1+m,kp) + &
                                               & CONJG(abcloPr(3,m,nkvec,lo)) *  cx(l+1+m,kp) )
                        END DO   
                     END DO
                  ELSE
                     DO m=-l,l
                        DO kp = 1,lapw%nv(igSpin)
                           hmat%data_c(lorow,kp) = hmat%data_c(lorow,kp) &
                                               & + chi * invsfct * ( &
                                               & CONJG(abcloPr(1,m,nkvec,lo)) *  ax((l+1+m),kp) + &
                                               & CONJG(abcloPr(2,m,nkvec,lo)) *  bx((l+1+m),kp) + &
                                               & CONJG(abcloPr(3,m,nkvec,lo)) *  cx((l+1+m),kp) )
                        END DO   
                     END DO
                  END IF
               END DO
               !$acc end kernels
               !$acc end data
            END DO
            CALL timestop("LO-LAPW")
         END IF
         CALL timestart("LO-LO")
         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abcoeffs,abclo,abcoeffsPr,abcloPr) &
         !$acc & copyin(atoms,lapw,lapwPr,tlmplm,lapw%nv(:),lapwPr%nv(:))&
         !$acc & copyin(tlmplm%tuloulo_newer(:,:,:,:,ntyp,ilSpinPr,ilSpin))&
         !$acc & copyin(tlmplm%h_loc_LO(:,:,ntyp,ilSpinPr,ilSpin),tlmplm%h_LO(:,:,:,ilSpinPr,ilSpin),tlmplm%h_LO2(:,:,:,ilSpinPr,ilSpin),tlmplm%h_loc2_nonsph)&
         !$acc & copyin(lapw%index_lo(:,na),lapwPr%index_lo(:,na),tlmplm%h_loc2,atoms%llo(:,ntyp),atoms%nlo(ntyp))&
         !$acc & copyin(fmpi, fmpi%n_size, fmpi%n_rank)&
         !$acc & default(none)
         DO lo = 1,atoms%nlo(ntyp)
            l = atoms%llo(lo,ntyp)
            ! Calculate the hamiltonian matrix elements with other local
            ! orbitals at the same atom and with itself
            DO nkvec = 1,invsfct* (2*l+1)
               locol = lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
               IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
                  locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
                  ! Calculate the Hamiltonian matrix elements with different
                  ! local orbitals at the same atom
                  DO lop = 1, MERGE(lo,atoms%nlo(ntyp),igSpinPr==igSpin.AND..NOT.l_fullj)
                     !IF (lop==lo) CYCLE No sepcial treatment needed for same LO
                     lp = atoms%llo(lop,ntyp)
                     DO nkvecp = 1,invsfct* (2*lp+1)
                        lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lop,na)+nkvecp
                        if (lorow>lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec) cycle
                        DO m = -l,l
                           lm = l*(l+1) + m
                           DO mp = -lp,lp
                              lmp = lp* (lp+1) + mp
                              s = tlmplm%h_loc2_nonsph(ntyp)
                              ! This is the product abclo^* tlmplm abclo in the 1,2,3 index (a,b,c coeff)
                              axx = tlmplm%h_loc_LO(lmp,lm,ntyp,ilSpinPr,ilSpin)*abclo(1,m,nkvec,lo)+ &
                                    tlmplm%h_loc_LO(lmp,lm+s,ntyp,ilSpinPr,ilSpin)*abclo(2,m,nkvec,lo) + &
                                    tlmplm%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin)*abclo(3,m,nkvec,lo)
                              bxx = tlmplm%h_loc_LO(lmp+s,lm,ntyp,ilSpinPr,ilSpin)*abclo(1,m,nkvec,lo) + &
                                    tlmplm%h_loc_LO(lmp+s,lm+s,ntyp,ilSpinPr,ilSpin)*abclo(2,m,nkvec,lo)  + &
                                    tlmplm%h_LO(lmp+s,m,lo+mlo,ilSpinPr,ilSpin) * abclo(3,m,nkvec,lo)
                              !cxx = tlmplm%tulou(lm,mp,lop+mlo,ilSpinPr,ilSpin) * abclo(1,m,nkvec,lo)+ &
                              !      tlmplm%tulod(lm,mp,lop+mlo,ilSpinPr,ilSpin)* abclo(2,m,nkvec,lo) + &
                              !      tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,ilSpinPr,ilSpin) * abclo(3,m,nkvec,lo)
                              cxx = tlmplm%h_LO2(lm,mp,lop+mlo,ilSpinPr,ilSpin) * abclo(1,m,nkvec,lo)+ &
                                    tlmplm%h_LO2(lm+s,mp,lop+mlo,ilSpinPr,ilSpin)* abclo(2,m,nkvec,lo) + &
                                    tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,ilSpinPr,ilSpin) * abclo(3,m,nkvec,lo)
                              prod= CONJG(abcloPr(1,mp,nkvecp,lop)) * axx + &
                                    CONJG(abcloPr(2,mp,nkvecp,lop)) * bxx + &
                                    CONJG(abcloPr(3,mp,nkvecp,lop)) * cxx
                            
                              IF (hmat%l_real) THEN
                                 hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + REAL(chi) * invsfct * REAL(prod)
                              ELSE
                                 hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + chi * prod
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO                  
               END IF !If this lo to be calculated by fmpi rank
            END DO
         END DO ! end of lo = 1,atoms%nlo loop
         !$acc end kernels
         CALL timestop("LO-LO")
      END IF

      !$acc end data
      !$acc end data
   END SUBROUTINE hlomat
END MODULE m_hlomat
