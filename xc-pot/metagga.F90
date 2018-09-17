!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_metagga

    TYPE t_realspace_potden
        REAL, ALLOCATABLE  :: is(:,:), mt(:,:)
    END TYPE t_realspace_potden

    PUBLIC  :: calc_EnergyDen
    PRIVATE :: calc_EnergyDen_auxillary_weights, t_realspace_potden, subtract_RS, multiply_RS

    INTERFACE OPERATOR (-)
        PROCEDURE subtract_RS
    END INTERFACE OPERATOR (-)

    INTERFACE OPERATOR (*)
        PROCEDURE multiply_RS
END INTERFACE OPERATOR (*)
CONTAINS
    SUBROUTINE calc_kinEnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, den, atoms, enpara, stars,&
                                 vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, kinEnergyDen)
#ifdef CPP_LIBXC 
        USE m_types_setup
        USE m_types_potden
        USE m_types_kpts
        USE m_types_mpi
        USE m_types_enpara
        USE m_types_misc
        USE m_types_regionCharges
        USE m_types_dos
        USE m_types_cdnval
        USE m_cdnval


        IMPLICIT NONE

        INTEGER,                  INTENT(in)           :: eig_id
        TYPE(t_mpi),              INTENT(in)           :: mpi
        TYPE(t_kpts),             INTENT(in)           :: kpts
        TYPE(t_noco),             INTENT(in)           :: noco
        TYPE(t_input),            INTENT(in)           :: input
        TYPE(t_banddos),          INTENT(in)           :: banddos
        TYPE(t_cell),             INTENT(in)           :: cell
        TYPE(t_potden),           INTENT(in)           :: den
        TYPE(t_atoms),            INTENT(in)           :: atoms 
        TYPE(t_enpara),           INTENT(in)           :: enpara
        TYPE(t_stars),            INTENT(in)           :: stars 
        TYPE(t_vacuum),           INTENT(in)           :: vacuum 
        TYPE(t_dimension),        INTENT(in)           :: DIMENSION
        TYPE(t_sphhar),           INTENT(in)           :: sphhar 
        TYPE(t_sym),              INTENT(in)           :: sym 
        TYPE(t_potden),           INTENT(in)           :: vTot
        TYPE(t_oneD),             INTENT(in)           :: oneD
        TYPE(t_results),          INTENT(in)           :: results
        TYPE(t_realspace_potden), INTENT(inout)        :: kinEnergyDen

        ! local vars

        TYPE(t_potden)                    :: EnergyDen
        TYPE(t_realspace_potden)          :: den_RS, EnergyDen_RS, vTot_RS


        CALL calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
                            vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, EnergyDen)
        
        CALL transform_to_grid(input, noco, sym, stars, cell, den, atoms, sphhar, EnergyDen, vTot, den_RS, EnergyDen_RS, vTot_RS)

        kinEnergyDen = EnergyDen_RS - vTot_RS * den_RS
#else
        USE m_juDFT_stop
        CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
    END SUBROUTINE calc_kinEnergyDen

    SUBROUTINE calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
                              vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, EnergyDen)
        ! calculates the energy density
        ! EnergyDen = \sum_i n_i(r) \varepsilon_i
        ! where n_i(r) is the one-particle density
        ! and \varepsilon_i are the eigenenergies
         
        
        USE m_types_setup
        USE m_types_potden
        USE m_types_kpts
        USE m_types_mpi
        USE m_types_enpara
        USE m_types_misc
        USE m_types_regionCharges
        USE m_types_dos
        USE m_types_cdnval
        USE m_cdnval

        IMPLICIT NONE

        INTEGER,           INTENT(in)           :: eig_id
        TYPE(t_mpi),       INTENT(in)           :: mpi
        TYPE(t_kpts),      INTENT(in)           :: kpts
        TYPE(t_noco),      INTENT(in)           :: noco
        TYPE(t_input),     INTENT(in)           :: input
        TYPE(t_banddos),   INTENT(in)           :: banddos
        TYPE(t_cell),      INTENT(in)           :: cell
        TYPE(t_atoms),     INTENT(in)           :: atoms 
        TYPE(t_enpara),    INTENT(in)           :: enpara
        TYPE(t_stars),     INTENT(in)           :: stars 
        TYPE(t_vacuum),    INTENT(in)           :: vacuum 
        TYPE(t_dimension), INTENT(in)           :: DIMENSION
        TYPE(t_sphhar),    INTENT(in)           :: sphhar 
        TYPE(t_sym),       INTENT(in)           :: sym 
        TYPE(t_potden),    INTENT(in)           :: vTot
        TYPE(t_oneD),      INTENT(in)           :: oneD
        TYPE(t_results),   INTENT(in)           :: results
        TYPE(t_potden),    INTENT(inout)        :: EnergyDen


        ! local
        INTEGER                         :: jspin
 
        TYPE(t_regionCharges)           :: regCharges
        TYPE(t_dos)                     :: dos
        TYPE(t_moments)                 :: moments
        TYPE(t_results)                 :: tmp_results
        TYPE(t_cdnvalJob)               :: cdnvalJob
        TYPE(t_potden)                  :: aux_den, real_den

        CALL regCharges%init(input, atoms)
        CALL dos%init(input,        atoms, DIMENSION, kpts, vacuum)
        CALL moments%init(input,    atoms)
        tmp_results = results

        DO jspin = 1,input%jspins
            CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin)

            ! replace brillouin weights with auxillary weights
            CALL calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, cdnvalJob%weights)

            CALL cdnval(eig_id, mpi, kpts, jspin, noco, input, banddos, cell, atoms, &
                        enpara, stars, vacuum, DIMENSION, sphhar, sym, vTot, oneD, cdnvalJob, &
                        EnergyDen, regCharges, dos, tmp_results, moments)
        ENDDO

    END SUBROUTINE calc_EnergyDen

    SUBROUTINE calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, f_ik)
        USE m_types_kpts
        USE m_eig66_io
        IMPLICIT NONE
        ! calculates new (auxillary-)weights as
        ! f_iks = w_iks * E_iks
        !, where  f_iks are the new (auxillary-)weights
        ! w_iks are the weights used in brillouin zone integration
        ! E_iks are the eigen energies


        INTEGER,      INTENT(in)        :: eig_id
        INTEGER,      INTENT(in)        :: jspin
        TYPE(t_kpts), INTENT(in)        :: kpts
        REAL,         INTENT(inout)     :: f_ik(:,:) ! f_ik(band_idx, kpt_idx)

        ! local vars
        REAL                       :: w_i(SIZE(f_ik,dim=1)), e_i(SIZE(f_ik,dim=1))
        INTEGER                    :: ikpt

        DO ikpt = 1,kpts%nkpt
            CALL read_eig(eig_id,ikpt,jspin, eig=e_i, w_iks=w_i)

            f_ik(:,ikpt) = e_i * w_i
        ENDDO
    END SUBROUTINE calc_EnergyDen_auxillary_weights

    SUBROUTINE transform_to_grid(input, noco, sym, stars, cell, den, atoms, sphhar, EnergyDen, vTot, den_RS, EnergyDen_RS, vTot_RS)
        USE m_types_potden
        USE m_types_setup
        USE m_types_xcpot_libxc
        USE m_types_xcpot
        USE m_juDFT_stop
        USE m_pw_tofrom_grid
        USE m_mt_tofrom_grid

        IMPLICIT NONE 
        
        TYPE(t_potden),           INTENT(in)        :: den, EnergyDen, vTot 
        TYPE(t_input),            INTENT(in)        :: input
        TYPE(t_noco),             INTENT(in)        :: noco
        TYPE(t_sym),              INTENT(in)        :: sym
        TYPE(t_stars),            INTENT(in)        :: stars
        TYPE(t_cell),             INTENT(in)        :: cell
        TYPE(t_atoms),            INTENT(in)        :: atoms
        TYPE(t_sphhar),           INTENT(in)        :: sphhar
        TYPE(t_realspace_potden), INTENT(out)       :: den_RS, EnergyDen_RS, vTot_RS ! could be changed to a real-space type

        !local vars
        TYPE(t_xcpot_libxc) ::aux_xcpot
        TYPE(t_gradients)   :: tmp_grad
        INTEGER, PARAMETER  :: id_corr = 9, id_exch = 1
        INTEGER             :: nsp, n



        !make some auxillary xcpot, that is not a GGA (we don't need gradients)
        CALL aux_xcpot%init(input%jspins, id_exch, id_corr, id_exch, id_corr)
        IF(aux_xcpot%is_gga()) CALL juDFT_error("aux_xcpot must not be GGA", &
                                                hint="choose id_corr and id_exch correctly")

        ! interstitial part
        CALL init_pw_grid(aux_xcpot,stars,sym,cell)

        CALL pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, den%pw,       tmp_grad, den_RS%is)
        CALL pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, EnergyDen%pw, tmp_grad, EnergyDen_RS%is)
        CALL pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, vTot%pw,      tmp_grad, vTot_RS%is)

        CALL finish_pw_grid()

        ! muffin tins
        nsp=(atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
        CALL init_mt_grid(nsp,input%jspins,atoms,sphhar,aux_xcpot,sym)

        DO n = 1,atoms%ntype
            CALL mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, den%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, den_RS%mt)
            CALL mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, EnergyDen%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, EnergyDen_RS%mt)
            CALL mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, vTot%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, vTot_RS%mt)
        ENDDO

        CALL finish_mt_grid()
    END SUBROUTINE transform_to_grid

    FUNCTION subtract_RS(rs1, rs2) RESULT(rs_out)
        IMPLICIT NONE

        TYPE(t_realspace_potden), INTENT(in)  :: rs1, rs2
        TYPE(t_realspace_potden)              :: rs_out

        WRITE (*,*) "MT subtraction not implemented"

        rs_out%is = rs1%is - rs2%is 
    END FUNCTION subtract_RS 

    FUNCTION multiply_RS(rs1, rs2) RESULT(rs_out)
        IMPLICIT NONE
        
        TYPE(t_realspace_potden), INTENT(in)  :: rs1, rs2
        TYPE(t_realspace_potden)              :: rs_out

        WRITE (*,*) "MT multiplication not implemented"

        rs_out%is = rs1%is * rs2%is 
    END FUNCTION multiply_RS

END MODULE m_metagga
