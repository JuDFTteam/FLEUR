!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_secvar
    use m_types_mat
    use m_types_mpi
    use m_types_lapw
    use m_types_atoms
    use m_judft
    IMPLICIT NONE
    PRIVATE
    TYPE :: t_secvar
        class(t_mat),allocatable :: mat(:,:) !The matrix representation of the operator
        integer :: ikpt,eig_id !k-point and eigenvalue id
        integer :: ne_first, ne_second !Number of eigenvalues in the first and second variation
        integer :: nmat_first, nmat_second !Size of the matrix in the first and second variation
        TYPE(t_mpi),   POINTER :: fmpi  => NULL()
        TYPE(t_lapw),  POINTER :: lapw  => NULL()
        TYPE(t_atoms), POINTER :: atoms => NULL()
        class(t_mat), allocatable :: eigvec !Eigenvectors of the second variation problem
        type(t_mat), allocatable :: zmat(:) !Matrix to store the eigenvectors in the original LAPW basis
        real,allocatable :: eigval(:) !Eigenvalues of the second variation problem
        real, allocatable :: eig(:,:) !Eigenvalues of the original problem for each spin channel
    contains
        procedure :: initialize
        procedure :: add_diagonal_elements 
        procedure :: diagonalize
        procedure :: store_eigvec
    END TYPE 
    public :: t_secvar
    
    contains
    
    SUBROUTINE initialize(this, l_noco, ikpt, eig_id, input, fmpi, lapw, atoms)
        use m_eig66_io
        use m_types_input
        CLASS(t_secvar), INTENT(INOUT) :: this
        LOGICAL, INTENT(IN) :: l_noco
        INTEGER, INTENT(IN) :: ikpt, eig_id
        TYPE(t_input),         INTENT(IN) :: input
        TYPE(t_mpi),   TARGET, INTENT(IN) :: fmpi
        TYPE(t_lapw),  TARGET, INTENT(IN) :: lapw
        TYPE(t_atoms), TARGET, INTENT(IN) :: atoms
        
        INTEGER :: i, j, jsp, nspin 
        
        ! Implementation of initialization for SOC second variation
        this%ikpt = ikpt
        this%eig_id = eig_id
        this%fmpi  => fmpi
        this%lapw  => lapw
        this%atoms => atoms

        this%ne_first=input%neig
        this%ne_second=input%neig*2

        this%nmat_first=lapw%nmat
        nspin=merge(2,1,l_noco)
        this%nmat_second=nspin*this%lapw%nmat

        allocate(t_mat::this%mat(1:nspin,1:nspin))
        do i=1,nspin
           do j=1,nspin
              call this%mat(i,j)%init(merge(.false.,input%l_real,l_noco),this%ne_first,this%ne_first)
           end do
        end do

        !Read the eigenvalues of the original problem for each spin channel, needed for the diagonalization in second variation
        
        allocate(this%zmat(input%jspins))
        allocate(this%eig(this%ne_first,input%jspins))
        DO jsp=1,input%jspins
            call this%zmat(jsp)%init(input%l_real,this%nmat_first,this%ne_first) 
           CALL read_eig(eig_id,this%ikpt,jsp, neig=this%ne_first,eig=this%eig(:,jsp),zmat=this%zmat(jsp))
        end do

    END SUBROUTINE initialize

    SUBROUTINE add_diagonal_elements(this)
        use m_types_mat
        use m_types_mpimat
        CLASS(t_secvar), INTENT(INOUT) :: this
        INTEGER :: i,jsp,jsp_in

        DO jsp=1,size(this%mat,1)
        jsp_in = min(size(this%eig,2),jsp)
        Associate(mat=> this%mat(jsp,jsp))
        select type(mat)
        type is (t_mpimat)
            CALL judft_bug("add_diagonal_elements not implemented for t_mpimat")
        type is(t_mat)
            if (mat%l_real) then  
                DO i=1,SIZE(mat%data_r,1)
                    mat%data_r(i,i) = mat%data_r(i,i) + this%eig(i,jsp_in)
                END DO
            else
                DO i=1,SIZE(mat%data_c,1)
                    mat%data_c(i,i) = mat%data_c(i,i) + cmplx(this%eig(i,jsp_in),0.0)
                END DO
            end if
        end select
        end associate
    enddo
    END SUBROUTINE add_diagonal_elements


    subroutine diagonalize(this)
        use m_eigen_diag
        CLASS(t_secvar), INTENT(INOUT) :: this
        class(t_mat),allocatable :: hmat

        allocate(this%eigval(this%ne_second))
        allocate(t_mat::this%eigvec)
        call this%eigvec%init(.false., this%ne_second, this%ne_second)
        CALL timestart("Matrix redistribution")
        ALLOCATE (hmat, mold=this%mat(1, 1))
        call hmat%init(.false., this%ne_second, this%ne_second)
        call hmat%copy(this%mat(1,1),1,1)
        call hmat%copy(this%mat(1,2),1,this%ne_first+1)
        call hmat%copy(this%mat(2,1),this%ne_first+1,1)
        call hmat%copy(this%mat(2,2),this%ne_first+1,this%ne_first+1)
        call timestop("Matrix redistribution")

        CALL timestart("Diagonalization in second variation")
        call eigen_diag_std(hmat, this%ne_second, this%eigval, this%eigvec)
        call timestop("Diagonalization in second variation")

    END SUBROUTINE diagonalize


    subroutine store_eigvec(this)
        use m_types_mat
        use m_eig66_io
        CLASS(t_secvar), INTENT(INOUT) :: this
        type(t_mat) :: eig, backtransformed
        integer :: jsp

        call eig%init(.false., this%ne_first, this%ne_second)
        call backtransformed%init(.false., this%nmat_first, this%ne_second)
        DO jsp=1,size(this%zmat)
            !Split the eigenvectors into a spin-up and down part
            call eig%copy(this%eigvec, 1, 1, m1=merge(1, this%ne_first+1, jsp==1), m2=1)
            eig%data_c=conjg(eig%data_c) !Conjugate the eigenvectors, since they are stored in the adjoint form in this%eigvec
            !Now transform the eigenvectors back to the original LAPW basis
            call this%zmat(jsp)%multiply(eig, backtransformed)
            !Store the backtransformed eigenvectors in the original matrix format
            call write_eig(this%eig_id, this%ikpt, jsp, neig=this%ne_second, neig_total=this%ne_second, eig=this%eigval, zmat=backtransformed)
        end do

    END SUBROUTINE store_eigvec
END MODULE m_types_secvar