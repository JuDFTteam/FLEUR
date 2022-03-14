!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_tlmplm

CONTAINS
    SUBROUTINE dfpt_tlmplm(atoms,sym,sphhar,input,noco,nococonv,enpara,hub1inp,hub1data,inden,vTot,vx,fmpi,td,ud,l_all_l,v1)
        USE m_types
        USE m_tlmplm

        IMPLICIT NONE

        TYPE(t_mpi),INTENT(IN)       :: fmpi
        TYPE(t_enpara),INTENT(IN) :: enpara
        TYPE(t_input),INTENT(IN)     :: input
        TYPE(t_noco),INTENT(IN)  :: noco
        TYPE(t_nococonv),INTENT(IN)  :: nococonv
        TYPE(t_sym),INTENT(IN)       :: sym
        TYPE(t_sphhar),INTENT(IN)    :: sphhar
        TYPE(t_atoms),INTENT(IN)     :: atoms
        TYPE(t_potden),INTENT(IN)    :: inDen
        TYPE(t_potden),INTENT(IN)    :: vTot,vx
        TYPE(t_tlmplm),INTENT(INOUT) :: td
        TYPE(t_usdus),INTENT(INOUT)  :: ud
        TYPE(t_hub1inp),INTENT(IN)   :: hub1inp
        TYPE(t_hub1data),INTENT(INOUT)::hub1data

        LOGICAL, INTENT(IN) :: l_all_l

        TYPE(t_potden), OPTIONAL, INTENT(IN) :: v1

        INTEGER :: iSpinV1, iSpinPr, iSpin, iPart, n
        COMPLEX :: one

        CALL timestart("tlmplm")
        CALL td%init(atoms,input%jspins,.FALSE.)
        !$OMP PARALLEL DO DEFAULT(NONE)&
        !$OMP PRIVATE(n,one,iSpinV1,iSpinPr,iSpin)&
        !$OMP SHARED(noco,nococonv,atoms,sym,sphhar,enpara,td,ud,vTot,vx,v1)&
        !$OMP SHARED(fmpi,input,hub1inp,hub1data,l_all_l)
        DO  n = 1,atoms%ntype
            DO iSpinV1 = 1, MERGE(4, input%jspins, any(noco%l_unrestrictMT))
                IF (iSpinV1.EQ.1) iSpinPr = 1; iSpin = 1
                IF (iSpinV1.EQ.2) iSpinPr = 2; iSpin = 2
                IF (iSpinV1.EQ.3) iSpinPr = 2; iSpin = 1
                IF (iSpinV1.EQ.4) iSpinPr = 1; iSpin = 2
                DO iPart = 1, 2
                    IF (iPart.EQ.1) one = CMPLX(1.0, 0.0)
                    IF (iPart.EQ.2) one = CMPLX(0.0, 1.0)
                    CALL tlmplm(n, sphhar, atoms, sym, enpara, nococonv, iSpinPr, iSpin, iSpinV1, fmpi, &
                              & vTot, vx, input, hub1inp, hub1data, td, ud, 0.0, 0, one, v1)
                END DO
            END DO
        END DO
        !$OMP END PARALLEL DO
        CALL timestop("tlmplm")

    END SUBROUTINE dfpt_tlmplm
END MODULE m_dfpt_tlmplm
