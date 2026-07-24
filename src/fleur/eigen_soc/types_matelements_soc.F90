!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_matelements_soc
    USE m_types_matelements
    USE m_types_mat
    USE m_types_rsoc
    USE m_types_atoms
    USE m_types_input
    USE m_types_mpi
    USE m_types_enpara
    USE m_types_lapw
    USE m_types_cell
    USE m_types_usdus
    USE m_types_noco
    USE m_types_nococonv
    USE m_types_sym
    USE m_types_potden
    USE m_judft
    IMPLICIT NONE

    TYPE, EXTENDS(t_matelements) :: t_matelements_soc
        TYPE(t_atoms),  POINTER :: atoms  => NULL()
        TYPE(t_noco),   POINTER :: noco   => NULL()
        TYPE(t_input),  POINTER :: input  => NULL()
        TYPE(t_sym),    POINTER :: sym    => NULL()
        TYPE(t_cell),   POINTER :: cell   => NULL()
        TYPE(t_enpara), POINTER :: enpara => NULL()
        TYPE(t_lapw),   POINTER :: lapw   => NULL()
        TYPE(t_potden), POINTER :: vtot   => NULL()
        TYPE(t_rsoc),   POINTER :: rsoc   => NULL()
        TYPE(t_mpi),    POINTER :: fmpi   => NULL()
        TYPE(t_nococonv),POINTER :: nococonv => NULL()
    CONTAINS
        PROCEDURE :: init => init
        PROCEDURE :: calc_matrix_elements => calc_matrix_elements
    END TYPE t_matelements_soc

CONTAINS

    subroutine init(this, atoms, noco, input, sym, cell, enpara, lapw, vtot, rsoc,  fmpi, nococonv)
        CLASS(t_matelements_soc), INTENT(INOUT) :: this
        TYPE(t_atoms),  TARGET, INTENT(IN) :: atoms
        TYPE(t_noco),   TARGET, INTENT(IN) :: noco
        TYPE(t_input),  TARGET, INTENT(IN) :: input
        TYPE(t_sym),    TARGET, INTENT(IN) :: sym
        TYPE(t_cell),   TARGET, INTENT(IN) :: cell
        TYPE(t_enpara), TARGET, INTENT(IN) :: enpara
        TYPE(t_lapw),   TARGET, INTENT(IN) :: lapw
        TYPE(t_potden), TARGET, INTENT(IN) :: vtot
        TYPE(t_rsoc),   TARGET, INTENT(IN) :: rsoc
        TYPE(t_mpi),    TARGET, INTENT(IN) :: fmpi
        TYPE(t_nococonv), TARGET, INTENT(IN) :: nococonv

        !SOC is a spinor operator coupling the two spin channels, while the
        !first-variation wave functions are pure spin states; the result
        !matrix hence has a 2x2 spinor structure
        this%spinoroperator = .TRUE.
        this%spinorwavefcts = .FALSE.

        this%atoms  => atoms
        this%noco   => noco
        this%input  => input
        this%sym    => sym
        this%cell   => cell
        this%enpara => enpara
        this%lapw   => lapw
        this%vtot   => vtot
        this%rsoc   => rsoc
        this%fmpi   => fmpi
        this%nococonv => nococonv
    end subroutine init

    subroutine calc_matrix_elements(this, zmat, abc, radfun, usdus)
        use m_types_abc
        use m_types_radfun
        use m_types_nococonv

        CLASS(t_matelements_soc), INTENT(INOUT) :: this
        TYPE(t_mat),    INTENT(IN) :: zMat(:)   !unused, SOC works on the abc coefficients only
        TYPE(t_abc),    INTENT(IN) :: abc(:,:)  !(2,ntype)
        TYPE(t_radfun), INTENT(IN) :: radfun(:) !unused, the radial integrals are precomputed in rsoc
        TYPE(t_usdus),  INTENT(IN) :: usdus     !unused

        INTEGER :: num_bands
        integer :: n, l, m, lm, ll1, jcof, icof
        integer :: i, j, j0, i1, j1, na, lm1, m1
        complex :: cof_lm

        if (.not.allocated(this%mat)) then
            call judft_bug("calc_matrix_elements: The result matrix is not allocated.")
        end if
        if (size(this%mat,1) /= 2 .or. size(this%mat,2) /= 2 ) then
            call judft_bug("calc_matrix_elements: The matrix must be a 2x2 spinor matrix.")
        end if
        if (size(abc,1) /= 2 .or. size(abc,2) /= this%atoms%ntype) then
            call judft_bug("calc_matrix_elements: The abc coefficients must have shape (2,ntype).")
        end if
        num_bands = size(abc(1,1)%cof,1)
        !matsize1 equals the global number of bands also for the row-cyclic
        !distributed t_mpimat blocks (full rows are local there)
        if (this%mat(1,1)%matsize1 /= num_bands) then
            call judft_bug("calc_matrix_elements: The matrix size does not match the abc coefficients.")
        end if

        DO n = 1,this%atoms%ntype
         DO j1=1,2
            DO i1=1,2
               !Distribute the column (ket) band index row-cyclically over the
               !sub_comm group; for n_size==1/n_rank==0 this is the full 1..num_bands
               !loop with j0==j, i.e. identical to the serial case.
               DO j = this%fmpi%n_rank+1, num_bands, this%fmpi%n_size
                  j0 = (j-1)/this%fmpi%n_size + 1
                  DO na = 1, this%atoms%neq(n)
                     DO l = 1,this%atoms%lmax(n)
                        ll1 = l*(l+1)
                        DO m = -l,l
                           lm = ll1 + m
                           DO jcof=1,abc(1,n)%n_r(l)
                              cof_lm = CMPLX(0.,0.)
                              DO m1 = -l,l
                                 lm1 = ll1 + m1
                                 cof_lm = cof_lm + this%rsoc%soangl(l,m,i1,l,m1,j1)*conjg(abc(j1,n)%cof(j,lm1,jcof,na))
                              ENDDO
                              DO icof=1,abc(1,n)%n_r(l)
                                 DO i=1, num_bands
                                    this%mat(i1,j1)%data_c(i,j0) = this%mat(i1,j1)%data_c(i,j0) +&
                                    abc(i1,n)%cof(i,lm,icof,na)*this%rsoc%rso(icof,jcof,n,l,i1,j1)*cof_lm
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO  ! na
               ENDDO     ! j-states
            ENDDO   !i1 spin
         ENDDO !!j1 spin
      enddo !n atom types

    END SUBROUTINE calc_matrix_elements
END MODULE m_types_matelements_soc  
