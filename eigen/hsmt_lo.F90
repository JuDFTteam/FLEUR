!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef _OPENACC
#define CPP_OMP !no OMP
#define CPP_ACC $acc
#else
#define CPP_OMP $OMP
#define CPP_ACC !no ACC
#endif
MODULE m_hsmt_lo
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC hsmt_lo
CONTAINS
  SUBROUTINE hsmt_lo(Input,Atoms,Sym,Cell,fmpi,Noco,nococonv,Lapw,Ud,Tlmplm,FjGj,N,Chi,ilSpinPr,ilSpin,igSpinPr,igSpin,Hmat,set0,l_fullj,l_ham,Smat,lapwq,fjgjq)
    USE m_hlomat
    USE m_slomat
    USE m_setabc1lo
    USE m_types_mpimat
    USE m_types
    USE m_hsmt_fjgj
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)      :: fmpi
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_nococonv),INTENT(IN) :: nococonv
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_usdus),INTENT(IN)    :: ud
    TYPE(t_tlmplm),INTENT(IN)   :: tlmplm
    TYPE(t_fjgj),INTENT(IN)     :: fjgj
    LOGICAL,INTENT(IN)          :: l_fullj, l_ham, set0  !if true, initialize the LO-part of the matrices with zeros
    TYPE(t_lapw),OPTIONAL,INTENT(IN) :: lapwq
    TYPE(t_fjgj), OPTIONAL, INTENT(IN) :: fjgjq

    CLASS(t_mat),INTENT(INOUT)::hmat
    CLASS(t_mat),INTENT(INOUT),OPTIONAL::smat

    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)   :: n
    INTEGER, INTENT (IN) :: ilSpinPr,ilSpin,igSpinPr,igSpin !spins
    COMPLEX, INTENT(IN)  :: chi

    !     ..
    !     .. Local Scalars ..
    INTEGER na,nn,usp
    INTEGER l,nkvec,kp
    !     ..
    !     .. Local Arrays ..
    REAL alo1(atoms%nlod,input%jspins),blo1(atoms%nlod,input%jspins),clo1(atoms%nlod,input%jspins)
    CALL timestart("LO setup")

    IF (set0) THEN
       SELECT TYPE (hmat)
       TYPE IS (t_mpimat)
          l = hmat%global_size2
       CLASS DEFAULT
          l = hmat%matsize2
       END SELECT

       !CPP_OMP PARALLEL DEFAULT(none) &
       !CPP_OMP SHARED(fmpi,l,lapw,hmat,smat,igSpin) &
       !CPP_OMP PRIVATE(nkvec,kp)
       !CPP_OMP DO
       !CPP_ACC kernels present(hmat,hmat%data_r,hmat%data_c)copyin(fmpi,lapw,lapw%nv)
       DO  nkvec =  fmpi%n_rank+1, l, fmpi%n_size
          IF( nkvec > lapw%nv(igSpin)) THEN
             kp=(nkvec-1)/fmpi%n_size+1
             IF (hmat%l_real) THEN
                hmat%data_r(:,kp) = 0.0
             ELSE
                hmat%data_c(:,kp) = CMPLX(0.0,0.0)
             ENDIF
          ENDIF
       ENDDO
       !CPP_ACC end kernels
       !CPP_OMP END DO
       IF ( present(smat)) THEN
          !CPP_OMP DO
          !CPP_ACC kernels present(smat,smat%data_r,smat%data_c)copyin(fmpi,lapw,lapw%nv)
          DO  nkvec =  fmpi%n_rank+1, l, fmpi%n_size
             IF( nkvec > lapw%nv(igSpin)) THEN
                kp=(nkvec-1)/fmpi%n_size+1
                IF (smat%l_real) THEN
                   smat%data_r(:,kp) = 0.0
                ELSE
                   smat%data_c(:,kp) = CMPLX(0.0,0.0)
                ENDIF
             ENDIF
          ENDDO
          !CPP_ACC end kernels
          !CPP_OMP END DO
       ENDIF
       !CPP_OMP END PARALLEL
    ENDIF

    na = atoms%firstAtom(n) - 1
    DO nn = 1,atoms%neq(n)
       na = na + 1
       IF ((sym%invsat(na).EQ.0) .OR. (sym%invsat(na).EQ.1)) THEN


          IF (atoms%nlo(n).GE.1) THEN


             !--->          set up the a,b and c  coefficients
             !--->          for the local orbitals, if necessary.
             !--->          actually, these are the fj,gj equivalents
             DO usp=min(ilSpinPr,ilSpin),max(ilSpinPr,ilSpin)
               CALL setabc1lo(atoms,n,ud,usp,alo1,blo1,clo1)
             enddo

             !--->       add the local orbital contribution to the overlap and
             !--->       hamiltonian matrix, if they are used for this atom.

               IF (ilSpinPr==ilSpin) THEN
                  IF (.NOT.PRESENT(smat)) THEN
                     IF (.NOT.PRESENT(lapwq)) CALL judft_error("Bug in hsmt_lo, called without smat")
                  ELSE
                     IF (PRESENT(lapwq)) THEN
                        CALL slomat(input,atoms,sym,fmpi,lapw,cell,nococonv,n,na,&
                           ilSpinPr,ud, alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr),fjgj,&
                           igSpinPr,igSpin,chi,smat,l_fullj,lapwq,fjgjq)
                     ELSE
                        CALL slomat(input,atoms,sym,fmpi,lapw,cell,nococonv,n,na,&
                           ilSpinPr,ud, alo1(:,ilSpinPr),blo1(:,ilSpinPr),clo1(:,ilSpinPr),fjgj,&
                           igSpinPr,igSpin,chi,smat,l_fullj)
                     END IF
                  END IF
               END IF
               CALL timestart("hlomat")
               IF (PRESENT(lapwq)) THEN
                  CALL hlomat(input,atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
                     n,na,fjgj,alo1,blo1,clo1,igSpinPr,igSpin,chi,hmat,l_fullj,l_ham,lapwq,fjgjq)
               ELSE
                  CALL hlomat(input,atoms,fmpi,lapw,ud,tlmplm,sym,cell,noco,nococonv,ilSpinPr,ilSpin,&
                     n,na,fjgj,alo1,blo1,clo1,igSpinPr,igSpin,chi,hmat,l_fullj,l_ham)
               END IF
               CALL timestop("hlomat")
            END IF
         END IF
         ! End loop over equivalent atoms
      END DO
      CALL timestop("LO setup")

      RETURN
   END SUBROUTINE hsmt_lo

END MODULE m_hsmt_lo
