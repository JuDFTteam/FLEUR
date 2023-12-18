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
MODULE m_hslomat
   use m_matmul_dgemm
   USE m_types
   USE m_hsmt_fjgj
   IMPLICIT NONE
  !***********************************************************************
  ! updates the hamiltonian  matrix with the contributions from the local
  ! orbitals.
  ! p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE hslomat(atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
       ntyp,na,fjgj,alo1,blo1,clo1, igSpinPr,igSpin,chi,l_fullj,hmat,smat,lapwq,fjgjq)

    USE m_hsmt_ab
    IMPLICIT NONE
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
    CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat
    LOGICAL, INTENT(IN) :: l_fullj

    TYPE(t_lapw), OPTIONAL, INTENT(IN),TARGET :: lapwq
    TYPE(t_fjgj), OPTIONAL, INTENT(IN),TARGET :: fjgjq
    !     ..
    ! Local Scalars
      INTEGER :: lo_lmax,ab_size,ab_size_Pr

      ! Local Arrays
      COMPLEX, ALLOCATABLE :: abCoeffs(:,:)
      COMPLEX, ALLOCATABLE :: abclo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: abCoeffsPr(:,:)
      COMPLEX, ALLOCATABLE :: abcloPr(:,:,:,:)

      TYPE(t_lapw) ,POINTER:: lapwPr
      TYPE(t_fjgj) ,POINTER:: fjgjPr

      
      IF (present(lapwq)) THEN
         lapwPr => lapwq
         fjgjPr => fjgjq
      ELSE
         lapwPr => lapw
         fjgjPr => fjgj
      END IF

      ! Synthesize a and b
      
      lo_lmax=maxval(atoms%llo)

      ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      ALLOCATE(abcloPr(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod))
      ALLOCATE(abCoeffsPr(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapwPr%nv)))
   
      !$acc data create(abcoeffs,abclo,abcoeffsPr,abcloPr)
      !$acc data copyin(alo1,blo1,clo1,fjgjPr,fjgjpr%fj,fjgjpr%gj,tlmplm,tlmplm%h_loc_LO,tlmplm%h_lo)
      CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpinPr,igSpinPr,ntyp,na,cell,lapwPr,fjgjPr,abCoeffsPr(:,:),ab_size_Pr,.TRUE.,abcloPr,alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr))

      !we need the "unprimed" abcoeffs
      IF (ilSpin==ilSpinPr.AND.igSpinPr==igSpin.AND..not.present(lapwq)) THEN
         !$acc kernels present(abcoeffs,abcoeffsPr,abclo,abcloPr)
         if (l_fullj) abcoeffs=abcoeffsPr !TODO automatic alloc on GPU????
         abclo=abcloPr
         !$acc end kernels
      ELSE     
         ALLOCATE(abCoeffs(0:2*atoms%lnonsph(ntyp)*(atoms%lnonsph(ntyp)+2)+1,MAXVAL(lapw%nv)))       
         CALL hsmt_ab(sym,atoms,noco,nococonv,ilSpin,igSpin,ntyp,na,cell,lapw,fjgj,abCoeffs(:,:),ab_size,.TRUE.,abclo,alo1(:,ilSpin),blo1(:,ilSpin),clo1(:,ilSpin))
      END IF
   
      

      IF ((sym%invsat(na) == 0) .OR. (sym%invsat(na) == 1)) THEN
         ! If this atom is the first of two atoms related by inversion, the
         ! contributions to the overlap matrix of both atoms are added simultaneously.
         ! Where it is made use of the fact, that the sum of these contributions
         ! is twice the real part of the contribution of each atom. Note, that
         ! in this case there are twice as many (2*(2*l+1)) k-vectors.
         ! (compare abccoflo and comments there).
         call hslomat_LAPW_LO(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abCoeffsPr,abclo,chi,hmat,smat)
 
         IF (l_fullj) THEN
            if (fmpi%n_size>1) call judft_bug("Full matrix not implemented in MPI case")
            call hslomat_LO_LAPW(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abCoeffs,abcloPr,chi,hmat,smat)
         END IF

         call hslomat_LO_LO(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abclo,abcloPr,l_fullj,chi,hmat,smat)
      END IF

      !$acc end data
      !$acc end data
   END SUBROUTINE hslomat

   subroutine hslomat_LAPW_LO(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abCoeffsPr,abclo,chi,hmat,smat)
      IMPLICIT NONE
      INTEGER,INTENT(IN)        :: na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_lapw),INTENT(IN)  :: lapw,lapwPr
      TYPE(t_tlmplm),INTENT(IN):: tlmplm
      TYPE(t_usdus),INTENT(IN) :: ud
      TYPE(t_mpi),INTENT(IN)   :: fmpi
      COMPLEX,INTENT(IN)        :: abCoeffsPr(0:,:),abclo(:,-atoms%llod:,:,:),chi
      CLASS(t_mat),INTENT(INOUT):: hmat
      CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat
   
   
      INTEGER :: lo,l,s,nkvec,locol,m,kp,invsfct,mlo,lm
      COMPLEX, ALLOCATABLE :: abcxPr(:,:,:)
      


      ALLOCATE(abcxPr(3,MAXVAL(lapwPr%nv),2*maxval(atoms%llo)+1))
      !$acc data create(abcxpr)
      invsfct=sym%invsat(na)+1
      mlo=sum(atoms%nlo(:ntyp-1))
      CALL timestart("LAPW-LO")
      !
      ! Calculate the matrix elements with the regular
      ! LAPW basis-functions   
      !
      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         s = tlmplm%h_loc2_nonsph(ntyp) 
         
        
         call blas_matmul(maxval(lapwPr%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_loc_LO(0:2*s-1,l*l:,ntyp,ilSpinPr,ilSpin),abcxPr(1,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
         call blas_matmul(maxval(lapwPr%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_loc_LO(0:2*s-1,s+l*l:,ntyp,ilSpinPr,ilSpin),abcxPr(2,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
         call blas_matmul(maxval(lapwPr%nv),2*l+1,2*s,abCoeffsPr,tlmplm%h_LO(0:2*s-1,-l:,lo+mlo,ilSpinPr,ilSpin),abcxPr(3,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'C')
       
         !LAPW LO contributions
         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abclo,abcxpr)&
         !$acc & copyin(lapw,lapw%nv,lapw%index_lo,fmpi,fmpi%n_size,fmpi%n_rank,lo,na,igSpin)
         DO nkvec = 1,invsfct*(2*l+1)
            locol= lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
            IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
               locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
               IF (hmat%l_real) THEN
                 DO m=-l,l
                     DO kp = 1,lapwPr%nv(igSpinPr)
                        hmat%data_r(kp,locol) = hmat%data_r(kp,locol) +  invsfct * dot_product(conjg(abclo(:,m,nkvec,lo)),abcxPr(:,kp,l+1+m))
                     ENDDO   
                 END DO
               ELSE
                 DO m=-l,l
                     DO kp = 1,lapwPr%nv(igSpinPr)
                        hmat%data_c(kp,locol) = hmat%data_c(kp,locol) + chi * dot_product(conjg(abclo(:,m,nkvec,lo)),abcxPr(:,kp,l+1+m))
                     ENDDO   
                     !call blas_matmul(1,3,lapw%nv(igSpinPr),abclo(:,m,nkvec,lo),abcxPr(:,:,l+1+m),hmat%data_c(:,lo),cmplx(1.,0.),cmplx(1.,0.),'C')
                  END DO
               END IF                 
            END IF
         END DO
         !$acc end kernels            
         !Calculate the same for the overlap
         IF (ilspin==ilSpinPr.and.present(smat)) THEN
            !$acc kernels present(abCoeffsPr,ud,ud%ddn,ud%uulon,ud%dulon)copyin(l,ntyp,ilspin)
            DO m = -l,l
               lm = l* (l+1) + m
               s = size(abCoeffsPr,1)/2
               abcxPr(1,:,l+1+m) = conjg(abCoeffsPr(lm,:)) !only "a-a" term
               abcxPr(2,:,l+1+m) = conjg(abCoeffsPr(s+lm,:))*ud%ddn(l,ntyp,ilspin) !only "b-b" term
               abcxPr(3,:,l+1+m) = conjg(abCoeffsPr(lm,:))*ud%uulon(lo,ntyp,ilspin) + conjg(abCoeffsPr(s+lm,:))*ud%dulon(lo,ntyp,ilspin) 
            enddo   
            !$acc end kernels
            !$acc kernels present(smat,smat%data_c,smat%data_r,abclo,abcxpr)&
            !$acc & copyin(lapw,lapw%nv,lapw%index_lo,fmpi,fmpi%n_size,fmpi%n_rank,lo,na,igSpin)
            DO nkvec = 1,invsfct*(2*l+1)
               locol= lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec ! This is the column of the matrix
               IF (MOD(locol-1,fmpi%n_size) == fmpi%n_rank) THEN ! Only this MPI rank calculates this column
                  locol=(locol-1)/fmpi%n_size+1 ! This is the column in local storage
                  IF (smat%l_real) THEN
                  DO m=-l,l
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           smat%data_r(kp,locol) = smat%data_r(kp,locol) + invsfct *  dot_product(conjg(abclo(:,m,nkvec,lo)),abcxPr(:,kp,l+1+m))
                        ENDDO   
                  END DO
                  ELSE
                  DO m=-l,l
                        DO kp = 1,lapwPr%nv(igSpinPr)
                           smat%data_c(kp,locol) = smat%data_c(kp,locol) + chi *  dot_product(conjg(abclo(:,m,nkvec,lo)),abcxPr(:,kp,l+1+m))
                        ENDDO   
                     END DO
                  END IF
               END IF
            END DO
            !$acc end kernels    
         ENDIF   
         
      ENDDO !loop over lo
      !$acc end data
      CALL timestop("LAPW-LO")
   END subroutine


   subroutine hslomat_LO_LAPW(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abCoeffs,abcloPr,chi,hmat,smat)
      IMPLICIT NONE
      INTEGER,INTENT(IN)        :: na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_lapw),INTENT(IN)  :: lapw,lapwPr
      TYPE(t_tlmplm),INTENT(IN):: tlmplm
      TYPE(t_usdus),INTENT(IN) :: ud
      TYPE(t_mpi),INTENT(IN)   :: fmpi
      COMPLEX,INTENT(IN)        :: abCoeffs(0:,:),abcloPr(:,-atoms%llod:,:,:),chi
      CLASS(t_mat),INTENT(INOUT):: hmat
      CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat
   
   
      INTEGER :: lo,l,s,nkvec,m,kp,invsfct,mlo,lm,lorow
      COMPLEX, ALLOCATABLE :: abcx(:,:,:)
      
      mlo=sum(atoms%nlo(:ntyp-1))

      ALLOCATE(abcx(3,MAXVAL(lapwPr%nv),2*maxval(atoms%llo)+1))
      !$acc data create(abcx)
      invsfct=sym%invsat(na)+1
      CALL timestart("LO-LAPW")
      DO lo = 1,atoms%nlo(ntyp)
         l = atoms%llo(lo,ntyp)
         s = tlmplm%h_loc2_nonsph(ntyp) 
         call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_loc_LO(0:2*s-1,l*l:,ntyp,ilSpinPr,ilSpin),abCoeffs,abcx(1,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'T','N')
         call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_loc_LO(0:2*s-1,s+l*l:,ntyp,ilSpinPr,ilSpin),abCoeffs,abcx(2,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'T','N')
         call blas_matmul(2*l+1,maxval(lapw%nv),2*s,tlmplm%h_LO2(0:2*s-1,-l:,lo+mlo,ilSpinPr,ilSpin),abCoeffs,abcx(3,:,:),cmplx(1.0,0.0),cmplx(0.0,0.0),'T','N')

         !$acc kernels present(hmat,hmat%data_c,hmat%data_r,abcloPr,abcx)&
         !$acc & copyin(lapwPr,lapwPr%nv,lapwPr%index_lo,fmpi,fmpi%n_size,fmpi%n_rank,invsfct,igSpinPr,lo,na,l)
         DO nkvec = 1,invsfct*(2*l+1)
            lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lo,na)+nkvec
            IF (hmat%l_real) THEN
               DO m=-l,l
                  DO kp = 1,lapw%nv(igSpin)
                     hmat%data_r(lorow,kp) = hmat%data_r(lorow,kp) + invsfct * dot_product(abcloPr(:,m,nkvec,lo),abcx(:,l+1+m,kp))
                  END DO   
               END DO
            ELSE
               DO m=-l,l
                  DO kp = 1,lapw%nv(igSpin)
                     hmat%data_c(lorow,kp) = hmat%data_c(lorow,kp)  + chi * dot_product(abcloPr(:,m,nkvec,lo),abcx(:,l+1+m,kp))
                  END DO   
               END DO
            END IF
         END DO
         !$acc end kernels
         !Calculate the same for the overlap
         IF (ilspin==ilspinPr.and.present(smat)) THEN
            !$acc kernels present(abCoeffs,ud,ud%ddn,ud%uulon,ud%dulon)copyin(l,ntyp,ilspin)
            DO m = -l,l
               lm = l* (l+1) + m
               s = size(abCoeffs,1)/2
               abcx(1,:,l+1+m) = abCoeffs(lm,:) !only "a-a" term
               abcx(2,:,l+1+m) = abCoeffs(s+lm,:)*ud%ddn(l,ntyp,ilspin) !only "b-b" term
               abcx(3,:,l+1+m) = abCoeffs(lm,:)*ud%uulon(lo,ntyp,ilspin) + conjg(abCoeffs(s+lm,:))*ud%dulon(lo,ntyp,ilspin) 
            enddo   
            !$acc end kernels
            !$acc kernels present(smat,smat%data_c,smat%data_r,abcloPr,abcx)&
            !$acc & copyin(lapw,lapw%nv,lapw%index_lo,fmpi,fmpi%n_size,fmpi%n_rank,lo,na,igSpin)
            DO nkvec = 1,invsfct*(2*l+1)
               lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lo,na)+nkvec
               IF (smat%l_real) THEN
                  DO m=-l,l
                     DO kp = 1,lapwPr%nv(igSpinPr)
                        smat%data_r(lorow,kp) = smat%data_r(lorow,kp) + invsfct *  dot_product(abcloPr(:,m,nkvec,lo),abcx(:,kp,l+1+m))
                     ENDDO   
                  END DO
               ELSE
                  DO m=-l,l
                     DO kp = 1,lapwPr%nv(igSpinPr)
                        smat%data_c(lorow,kp) = smat%data_c(lorow,kp) + chi *  dot_product(abcloPr(:,m,nkvec,lo),abcx(:,kp,l+1+m))
                     ENDDO   
                  END DO
               END IF  
            END DO
            !$acc end kernels    
         ENDIF  
         
      END DO !loop over lo
      !$acc end data
      CALL timestop("LO-LAPW")
   END subroutine


   subroutine hslomat_LO_LO(na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr,atoms,sym,lapw,lapwPr,tlmplm,ud,fmpi,abclo,abcloPr,l_full_matrix,chi,hmat,smat)
      IMPLICIT NONE
      INTEGER,INTENT(IN)        :: na,ntyp,ilspinPr,ilSpin,igSpin,igSpinPr
      TYPE(t_atoms),INTENT(IN) :: atoms
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_lapw),INTENT(IN)  :: lapw,lapwPr
      TYPE(t_tlmplm),INTENT(IN):: tlmplm
      TYPE(t_usdus),INTENT(IN) :: ud
      TYPE(t_mpi),INTENT(IN)   :: fmpi
      COMPLEX,INTENT(IN)        :: abclo(:,-atoms%llod:,:,:),abcloPr(:,-atoms%llod:,:,:),chi
      LOGICAL,INTENT(IN)        :: l_full_matrix
      CLASS(t_mat),INTENT(INOUT):: hmat
      CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat
   
   
      INTEGER  :: lo,l,s,nkvec,nkvecp,m,mp,lmp,lp,invsfct,mlo,lm,lorow,locol,lop
      COMPLEX  :: abcxx(3)
      mlo=sum(atoms%nlo(:ntyp-1))
      invsfct=sym%invsat(na)+1
      CALL timestart("LO-LO")
      !$acc kernels present(smat,smat%data_c,smat%data_r,hmat,hmat%data_c,hmat%data_r,abclo,abcloPr) &
      !$acc & copyin(atoms,lapw,lapwPr,tlmplm,lapw%nv(:),lapwPr%nv(:))&
      !$acc & copyin(tlmplm%tuloulo_newer(:,:,:,:,ntyp,ilSpinPr,ilSpin))&
      !$acc & copyin(tlmplm%h_loc_LO(:,:,ntyp,ilSpinPr,ilSpin),tlmplm%h_LO(:,:,:,ilSpinPr,ilSpin),tlmplm%h_LO2(:,:,:,ilSpinPr,ilSpin),tlmplm%h_loc2_nonsph)&
      !$acc & copyin(lapw%index_lo(:,na),lapwPr%index_lo(:,na),tlmplm%h_loc2,atoms%llo(:,ntyp),atoms%nlo(ntyp))&
      !$acc & copyin(fmpi, fmpi%n_size, fmpi%n_rank,ud,ud%ddn,ud%dulon,ud%uulon)&
      !$acc & create(abcxx)&
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
               DO lop = 1, MERGE(lo,atoms%nlo(ntyp),igSpinPr==igSpin.AND..NOT.l_full_matrix)
                  !IF (lop==lo) CYCLE No sepcial treatment needed for same LO
                  lp = atoms%llo(lop,ntyp)
                  DO nkvecp = 1,invsfct* (2*lp+1)
                     lorow = lapwPr%nv(igSpinPr)+lapwPr%index_lo(lop,na)+nkvecp
                     if (lorow>lapw%nv(igSpin)+lapw%index_lo(lo,na)+nkvec.AND..NOT.l_full_matrix) cycle
                     DO m = -l,l
                        lm = l*(l+1) + m
                        DO mp = -lp,lp
                           lmp = lp* (lp+1) + mp
                           s = tlmplm%h_loc2_nonsph(ntyp)
                           ! This is the product abclo^* tlmplm abclo in the 1,2,3 index (a,b,c coeff)
                           abcxx(1) = tlmplm%h_loc_LO(lmp,lm,ntyp,ilSpinPr,ilSpin)*abclo(1,m,nkvec,lo)+ &
                                 tlmplm%h_loc_LO(lmp,lm+s,ntyp,ilSpinPr,ilSpin)*abclo(2,m,nkvec,lo) + &
                                 tlmplm%h_LO(lmp,m,lo+mlo,ilSpinPr,ilSpin)*abclo(3,m,nkvec,lo)
                           abcxx(2) = tlmplm%h_loc_LO(lmp+s,lm,ntyp,ilSpinPr,ilSpin)*abclo(1,m,nkvec,lo) + &
                                 tlmplm%h_loc_LO(lmp+s,lm+s,ntyp,ilSpinPr,ilSpin)*abclo(2,m,nkvec,lo)  + &
                                 tlmplm%h_LO(lmp+s,m,lo+mlo,ilSpinPr,ilSpin) * abclo(3,m,nkvec,lo)
                           abcxx(3) = tlmplm%h_LO2(lm,mp,lop+mlo,ilSpinPr,ilSpin) * abclo(1,m,nkvec,lo)+ &
                                 tlmplm%h_LO2(lm+s,mp,lop+mlo,ilSpinPr,ilSpin)* abclo(2,m,nkvec,lo) + &
                                 tlmplm%tuloulo_newer(mp,m,lop,lo,ntyp,ilSpinPr,ilSpin) * abclo(3,m,nkvec,lo)
                         
                           IF (hmat%l_real) THEN
                              hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + REAL(chi) * invsfct * dot_product(abcloPr(:,mp,nkvecp,lop),abcxx)
                           ELSE
                              hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + chi * dot_product(abcloPr(:,mp,nkvecp,lop),abcxx)
                           END IF
                        enddo   
                        if (ilspin==ilspinPr .AND. l==lp .and. present(smat)) THEN !Same for overlap matrix
                              ! This is the product abclo^* <u|u> * abclo in the 1,2,3 index (a,b,c coeff)
                           abcxx(1) = abclo(1,m,nkvec,lo)   +  0    + ud%uulon(lo,ntyp,ilspin) * abclo(3,m,nkvec,lo)
                           abcxx(2) = 0+ ud%ddn(l,ntyp,ilspin)*abclo(2,m,nkvec,lo)  + ud%dulon(lo,ntyp,ilspin) * abclo(3,m,nkvec,lo)
                           abcxx(3) = ud%uulon(lop,ntyp,ilspin)* abclo(1,m,nkvec,lo)+ &
                                 ud%dulon(lop,ntyp,ilspin)* abclo(2,m,nkvec,lo) + &
                                 ud%uloulopn(lop,lo,ntyp,ilspin)* abclo(3,m,nkvec,lo)
                        
                           IF (smat%l_real) THEN
                              smat%data_r(lorow,locol) = smat%data_r(lorow,locol) +  invsfct * dot_product(abcloPr(:,m,nkvecp,lop),abcxx)
                           ELSE
                              smat%data_c(lorow,locol) = smat%data_c(lorow,locol) + chi * dot_product(abcloPr(:,m,nkvecp,lop),abcxx)
                           END IF
                     
                        ENDIF   
                     END DO
                  END DO
               END DO                  
            END IF !If this lo to be calculated by fmpi rank
         END DO
      END DO ! end of lo = 1,atoms%nlo loop
      !$acc end kernels
      CALL timestop("LO-LO")
   

   end subroutine


END MODULE m_hslomat
