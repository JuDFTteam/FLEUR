!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_tlmplm

CONTAINS
   SUBROUTINE dfpt_tlmplm(atoms,sym,sphhar,input,noco,enpara,hub1inp,hub1data,vTot,fmpi,tdV1,v1real,v1imag,conj_V,iDtype_col)
      !! Get the (lm) matrix elements for the perturbed potential, which differs slightly from the base
      !! case of tlmplm for V/H.
      USE m_types
      USE m_tlmplm

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

      INTEGER :: iSpinV1, iSpinPr, iSpin, iPart, n, offs, nlims(2)
      COMPLEX :: one

      REAL, ALLOCATABLE :: vr1(:, :)

      TYPE(t_usdus)    :: uddummy
      TYPE(t_potden)   :: vxdummy
      TYPE(t_nococonv) :: nococonvdummy

        ALLOCATE( vr1(SIZE(v1real%mt,1),0:SIZE(v1real%mt,2)-1))

        call uddummy%init(atoms,input%jspins)
        CALL timestart("tlmplm")
        CALL tdV1%init(atoms,input%jspins,.FALSE.)

        nlims(1) = 1
        nlims(2) = atoms%ntype
        IF (PRESENT(iDtype_col)) nlims = [iDtype_col,iDtype_col]

        !$OMP PARALLEL DO DEFAULT(NONE)&
        !$OMP PRIVATE(n,one,iSpinV1,iSpinPr,iSpin,vr1,offs)&
        !$OMP SHARED(noco,nococonvdummy,atoms,sym,sphhar,enpara,tdV1,uddummy,vTot,vxdummy,v1real,v1imag,conj_V,nlims)&
        !$OMP SHARED(fmpi,input,hub1inp,hub1data)
        DO n = nlims(1), nlims(2)
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
                    CALL tlmplm(n, sphhar, atoms, sym, enpara, nococonvdummy, iSpinPr, iSpin, iSpinV1, fmpi, &
                              & vTot, vxdummy, input, hub1inp, hub1data, tdV1, uddummy, 0.0, one, .TRUE., vr1)
                END DO
            END DO

            offs = tdV1%h_loc2_nonsph(n)
            tdV1%h_loc_nonsph(0:offs-1,0:offs-1,n,:,:)    = tdV1%h_loc(0:offs-1,0:offs-1,n,:,:)
            tdV1%h_loc_nonsph(offs:offs+offs-1,0:offs-1,n,:,:)  = tdV1%h_loc(tdV1%h_loc2(n):offs+tdV1%h_loc2(n)-1,0:offs-1,n,:,:)
            tdV1%h_loc_nonsph(0:offs-1,offs:offs+offs-1,n,:,:)  = tdV1%h_loc(0:offs-1,tdV1%h_loc2(n):offs+tdV1%h_loc2(n)-1,n,:,:)
            tdV1%h_loc_nonsph(offs:offs+offs-1,offs:offs+offs-1,n,:,:)= tdV1%h_loc(tdV1%h_loc2(n):offs+tdV1%h_loc2(n)-1,tdV1%h_loc2(n):offs+tdV1%h_loc2(n)-1,n,:,:)
        END DO
        !$OMP END PARALLEL DO
        CALL timestop("tlmplm")

    END SUBROUTINE dfpt_tlmplm
END MODULE m_dfpt_tlmplm
