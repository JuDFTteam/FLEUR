!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_tlmplm

   implicit none
CONTAINS
   SUBROUTINE dfpt_tlmplm(atoms,sym,sphhar,input,noco,enpara,hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,conj_V,iDtype_col)
      !! Get the (lm) matrix elements for the perturbed potential, which differs slightly from the base
      !! case of tlmplm for V/H.
      USE m_types
      USE m_local_hamiltonian, ONLY: add_nonsph, extract_nonsph

      IMPLICIT NONE

      TYPE(t_mpi),INTENT(IN)       :: fmpi
      TYPE(t_enpara),INTENT(IN) :: enpara
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)  :: noco
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_potden),INTENT(IN)    :: vTot
      TYPE(t_tlmplm),INTENT(INOUT) :: tdV1
      TYPE(t_hub1inp),INTENT(IN)   :: hub1inp
      TYPE(t_hub1data),INTENT(INOUT)::hub1data

      TYPE(t_potden), INTENT(IN) :: v1real, v1imag

      LOGICAL, INTENT(IN) :: conj_V

      INTEGER, INTENT(IN), OPTIONAL :: iDtype_col

      INTEGER :: iSpinV1, iSpinPr, iSpin, iPart, n, nlims(2), j1, j2
      COMPLEX :: one

      REAL, ALLOCATABLE :: vr1(:, :)

        ALLOCATE( vr1(SIZE(v1real%mt,1),0:SIZE(v1real%mt,2)-1))

        CALL timestart("tlmplm")
        CALL tdV1%init(atoms,input%jspins,.FALSE.)

        nlims(1) = 1
        nlims(2) = atoms%ntype
        IF (PRESENT(iDtype_col)) nlims = [iDtype_col,iDtype_col]

        !$OMP PARALLEL DO DEFAULT(NONE)&
        !$OMP PRIVATE(n,one,iSpinV1,iSpinPr,iSpin,iPart,vr1,j1,j2)&
        !$OMP SHARED(noco,atoms,sym,sphhar,enpara,tdV1,vTot,v1real,v1imag,conj_V,nlims)&
        !$OMP SHARED(fmpi,input,hub1inp,hub1data)
        DO n = nlims(1), nlims(2)
            CALL tdV1%radfun(n)%generate_radial_functions(atoms,input,enpara,fmpi,vTot,n,hub1data)
            DO iSpinV1 = 1, MERGE(4, input%jspins, any(noco%l_unrestrictMT))
                iSpinPr = 1; iSpin = 1
                IF (iSpinV1.EQ.2.OR.iSpinV1.EQ.3) iSpinPr = 2
                IF (iSpinV1.EQ.2.OR.iSpinV1.EQ.4) iSpin   = 2
                DO iPart = 1, 2
                    IF (.NOT.conj_V) THEN
                       IF (iPart.EQ.1) one = CMPLX(1.0, 0.0)
                       IF (iPart.EQ.2) one = CMPLX(0.0, 1.0)
                       IF (iPart.EQ.1) vr1 = v1real%mt(:, :, n, iSpinV1)
                       IF (iPart.EQ.2) vr1 = v1imag%mt(:, :, n, iSpinV1)
                    ELSE
                       IF (iPart.EQ.1) one = CMPLX(1.0, 0.0)
                       IF (iPart.EQ.2) one = CMPLX(0.0,-1.0)
                       IF (iSpinV1==1.OR.iSpinV1==2) THEN
                          IF (iPart.EQ.1) vr1 = v1real%mt(:, :, n, iSpinV1)
                          IF (iPart.EQ.2) vr1 = v1imag%mt(:, :, n, iSpinV1)
                       ELSE IF (iSpinV1==3) THEN
                          IF (iPart.EQ.1) vr1 = v1real%mt(:, :, n, 4)
                          IF (iPart.EQ.2) vr1 = v1imag%mt(:, :, n, 4)
                       ELSE
                          IF (iPart.EQ.1) vr1 = v1real%mt(:, :, n, 3)
                          IF (iPart.EQ.2) vr1 = v1imag%mt(:, :, n, 3)
                       END IF
                    END IF
                    CALL add_nonsph(tdV1,n,atoms,sym,sphhar,input,hub1inp,vr1,0,iSpinPr,iSpin,one)
                END DO
            END DO
            DO j1 = 1, SIZE(tdV1%h,4)
               DO j2 = 1, SIZE(tdV1%h,5)
                  CALL extract_nonsph(tdV1,atoms,n,j1,j2)
               END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
        CALL timestop("tlmplm")

    END SUBROUTINE dfpt_tlmplm
END MODULE m_dfpt_tlmplm
