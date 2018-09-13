!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_metagga

    type t_realspace_potden
        real, allocatable  :: is(:,:), mt(:,:)
    end type t_realspace_potden

    public  :: calc_EnergyDen
    private :: calc_EnergyDen_auxillary_weights, t_realspace_potden, subtract_RS, multiply_RS

    interface operator (-)
        procedure subtract_RS
    end interface operator (-)

    interface operator (*)
        procedure multiply_RS
end interface operator (*)
contains
    subroutine calc_kinEnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, den, atoms, enpara, stars,&
                                 vacuum, dimension, sphhar, sym, vTot, oneD, results, kinEnergyDen)
#ifdef CPP_LIBXC 
        use m_types_setup
        use m_types_potden
        use m_types_kpts
        use m_types_mpi
        use m_types_enpara
        use m_types_misc
        use m_types_regionCharges
        use m_types_dos
        use m_types_cdnval
        use m_cdnval


        implicit none

        integer,                  intent(in)           :: eig_id
        type(t_mpi),              intent(in)           :: mpi
        type(t_kpts),             intent(in)           :: kpts
        type(t_noco),             intent(in)           :: noco
        type(t_input),            intent(in)           :: input
        type(t_banddos),          intent(in)           :: banddos
        type(t_cell),             intent(in)           :: cell
        type(t_potden),           intent(in)           :: den
        type(t_atoms),            intent(in)           :: atoms 
        type(t_enpara),           intent(in)           :: enpara
        type(t_stars),            intent(in)           :: stars 
        type(t_vacuum),           intent(in)           :: vacuum 
        type(t_dimension),        intent(in)           :: dimension
        type(t_sphhar),           intent(in)           :: sphhar 
        type(t_sym),              intent(in)           :: sym 
        type(t_potden),           intent(in)           :: vTot
        type(t_oneD),             intent(in)           :: oneD
        type(t_results),          intent(in)           :: results
        type(t_realspace_potden), intent(inout)        :: kinEnergyDen

        ! local vars

        type(t_potden)                    :: EnergyDen
        type(t_realspace_potden)          :: den_RS, EnergyDen_RS, vTot_RS


        call calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
                            vacuum, dimension, sphhar, sym, vTot, oneD, results, EnergyDen)
        
        call transform_to_grid(input, noco, sym, stars, cell, den, atoms, sphhar, EnergyDen, vTot, den_RS, EnergyDen_RS, vTot_RS)

        kinEnergyDen = EnergyDen_RS - vTot_RS * den_RS
#else
        use m_juDFT_stop
        call juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
    end subroutine calc_kinEnergyDen

    subroutine calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
                              vacuum, dimension, sphhar, sym, vTot, oneD, results, EnergyDen)
        ! calculates the energy density
        ! EnergyDen = \sum_i n_i(r) \varepsilon_i
        ! where n_i(r) is the one-particle density
        ! and \varepsilon_i are the eigenenergies
         
        
        use m_types_setup
        use m_types_potden
        use m_types_kpts
        use m_types_mpi
        use m_types_enpara
        use m_types_misc
        use m_types_regionCharges
        use m_types_dos
        use m_types_cdnval
        use m_cdnval

        implicit none

        integer,           intent(in)           :: eig_id
        type(t_mpi),       intent(in)           :: mpi
        type(t_kpts),      intent(in)           :: kpts
        type(t_noco),      intent(in)           :: noco
        type(t_input),     intent(in)           :: input
        type(t_banddos),   intent(in)           :: banddos
        type(t_cell),      intent(in)           :: cell
        type(t_atoms),     intent(in)           :: atoms 
        type(t_enpara),    intent(in)           :: enpara
        type(t_stars),     intent(in)           :: stars 
        type(t_vacuum),    intent(in)           :: vacuum 
        type(t_dimension), intent(in)           :: dimension
        type(t_sphhar),    intent(in)           :: sphhar 
        type(t_sym),       intent(in)           :: sym 
        type(t_potden),    intent(in)           :: vTot
        type(t_oneD),      intent(in)           :: oneD
        type(t_results),   intent(in)           :: results
        type(t_potden),    intent(inout)        :: EnergyDen


        ! local
        integer                         :: jspin
 
        type(t_regionCharges)           :: regCharges
        type(t_dos)                     :: dos
        type(t_moments)                 :: moments
        type(t_results)                 :: tmp_results
        type(t_cdnvalJob)               :: cdnvalJob
        type(t_potden)                  :: aux_den, real_den

        call regCharges%init(input, atoms)
        call dos%init(input,        atoms, dimension, kpts, vacuum)
        call moments%init(input,    atoms)
        tmp_results = results

        do jspin = 1,input%jspins
            call cdnvalJob%init(mpi,input,kpts,noco,results,jspin)

            ! replace brillouin weights with auxillary weights
            call calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, cdnvalJob%weights)

            call cdnval(eig_id, mpi, kpts, jspin, noco, input, banddos, cell, atoms, &
                        enpara, stars, vacuum, dimension, sphhar, sym, vTot, oneD, cdnvalJob, &
                        EnergyDen, regCharges, dos, tmp_results, moments)
        enddo

    end subroutine calc_EnergyDen

    subroutine calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, f_ik)
        use m_types_kpts
        use m_eig66_io
        implicit none
        ! calculates new (auxillary-)weights as
        ! f_iks = w_iks * E_iks
        !, where  f_iks are the new (auxillary-)weights
        ! w_iks are the weights used in brillouin zone integration
        ! E_iks are the eigen energies


        integer,      intent(in)        :: eig_id
        integer,      intent(in)        :: jspin
        type(t_kpts), intent(in)        :: kpts
        real,         intent(inout)     :: f_ik(:,:) ! f_ik(band_idx, kpt_idx)

        ! local vars
        real                       :: w_i(size(f_ik,dim=1)), e_i(size(f_ik,dim=1))
        integer                    :: ikpt

        do ikpt = 1,kpts%nkpt
            call read_eig(eig_id,ikpt,jspin, eig=e_i, w_iks=w_i)

            f_ik(:,ikpt) = e_i * w_i
        enddo
    end subroutine calc_EnergyDen_auxillary_weights

    subroutine transform_to_grid(input, noco, sym, stars, cell, den, atoms, sphhar, EnergyDen, vTot, den_RS, EnergyDen_RS, vTot_RS)
        use m_types_potden
        use m_types_setup
        use m_types_xcpot_libxc
        use m_types_xcpot
        use m_juDFT_stop
        use m_pw_tofrom_grid
        use m_mt_tofrom_grid

        implicit none 
        
        type(t_potden),           intent(in)        :: den, EnergyDen, vTot 
        type(t_input),            intent(in)        :: input
        type(t_noco),             intent(in)        :: noco
        type(t_sym),              intent(in)        :: sym
        type(t_stars),            intent(in)        :: stars
        type(t_cell),             intent(in)        :: cell
        type(t_atoms),            intent(in)        :: atoms
        type(t_sphhar),           intent(in)        :: sphhar
        type(t_realspace_potden), intent(out)       :: den_RS, EnergyDen_RS, vTot_RS ! could be changed to a real-space type

        !local vars
        type(t_xcpot_libxc) ::aux_xcpot
        type(t_gradients)   :: tmp_grad
        integer, parameter  :: id_corr = 9, id_exch = 1
        integer             :: nsp, n



        !make some auxillary xcpot, that is not a GGA (we don't need gradients)
        call aux_xcpot%init(input%jspins, id_exch, id_corr)
        if(aux_xcpot%is_gga()) call juDFT_error("aux_xcpot must not be GGA", &
                                                hint="choose id_corr and id_exch correctly")

        ! interstitial part
        call init_pw_grid(aux_xcpot,stars,sym,cell)

        call pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, den%pw,       tmp_grad, den_RS%is)
        call pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, EnergyDen%pw, tmp_grad, EnergyDen_RS%is)
        call pw_to_grid(aux_xcpot, input%jspins, noco%l_noco, stars, cell, vTot%pw,      tmp_grad, vTot_RS%is)

        call finish_pw_grid()

        ! muffin tins
        nsp=(atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
        call init_mt_grid(nsp,input%jspins,atoms,sphhar,aux_xcpot,sym)

        do n = 1,atoms%ntype
            call mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, den%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, den_RS%mt)
            call mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, EnergyDen%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, EnergyDen_RS%mt)
            call mt_to_grid(aux_xcpot, input%jspins, atoms,    sphhar, vTot%mt(:,0:,n,:), &
                            nsp,       n,            tmp_grad, vTot_RS%mt)
        enddo

        call finish_mt_grid()
    end subroutine transform_to_grid

    function subtract_RS(rs1, rs2) result(rs_out)
        implicit none

        type(t_realspace_potden), intent(in)  :: rs1, rs2
        type(t_realspace_potden)              :: rs_out

        write (*,*) "MT subtraction not implemented"

        rs_out%is = rs1%is - rs2%is 
    end function subtract_RS 

    function multiply_RS(rs1, rs2) result(rs_out)
        implicit none
        
        type(t_realspace_potden), intent(in)  :: rs1, rs2
        type(t_realspace_potden)              :: rs_out

        write (*,*) "MT multiplication not implemented"

        rs_out%is = rs1%is * rs2%is 
    end function multiply_RS

end module m_metagga
