module m_opc_setup

    use m_slater
    use m_types

    implicit none

    contains

    subroutine opc_setup(input, atoms, fmpi, v, den, jspin, corrections)

        type(t_input),  intent(in)  :: input
        type(t_atoms),  intent(in)  :: atoms
        type(t_mpi),    intent(in)  :: fmpi
        type(t_potden), intent(in)  :: v, den
        integer,        intent(in)  :: jspin
        real, allocatable, intent(out) :: corrections(:)

        complex :: mmpmat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
        real, allocatable :: f(:,:)

        integer :: i_opc, l, atomType, index, m
        real :: orb_z

        allocate(corrections(atoms%n_opc), source=0.0)

        call slater(input, jspin, atoms, v%mt(:,0,:,jspin), l_write=fmpi%irank==0, slater_parameters=f)
        do i_opc = 1, atoms%n_opc
            atomType = atoms%lda_opc(i_opc)%atomType
            l = atoms%lda_opc(i_opc)%l

            !Calculate the expectation value of the orbital moment
            !from den%mmpmat (denisty matrices for OPC are behind LDA+U and LDA+HIA)

            index = atoms%n_u + atoms%n_hia + i_opc
            mmpmat = den%mmpmat(:,:,index, jspin)

            orb_z = 0.0
            do m = - l, l
                orb_z = orb_z + m * real(mmpmat(m,m))
            enddo
            !calculate the racah parameter as B = (9F2 - 5F4)/441
            corrections(i_opc) = - (9*f(1,i_opc) - 5*f(2,i_opc))/441.0 * orb_z
        enddo
    
    end subroutine

end module m_opc_setup