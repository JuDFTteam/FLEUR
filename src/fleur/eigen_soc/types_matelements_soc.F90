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
        PROCEDURE :: add_matrix_elements => add_matrix_elements
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

    subroutine add_matrix_elements(this, zmat, mat_inout)
        use m_types_abc
        use m_types_radfun
        use m_types_nococonv

        CLASS(t_matelements_soc), INTENT(INOUT) :: this
        TYPE(t_mat),  INTENT(IN) :: zMat(:)
        CLASS(t_mat), INTENT(INOUT),optional,target :: mat_inout(:,:)

        TYPE(t_mat), pointer::mat(:,:)
      class(t_mat), pointer :: mat_block

        type(t_abc) :: abc(2)
        type(t_radfun) :: radfun
        type(t_usdus) :: usdus
        INTEGER :: num_bands
        integer :: n, l, m, lm, ll1, jcof, icof
      integer :: i, j, i1, j1, na, lm1, m1, i1_int, j1_int
        complex :: cof_lm


        if (present(mat_inout)) then
            mat => mat_inout
        else
            mat => this%mat
        end if   

        if (size(mat,1) /= 2 .or. size(mat,2) /= 2 ) then
            call judft_bug("add_matrix_elements: The matrix must be a 2x2 spinor matrix.")
        end if
        num_bands = mat(1,1)%matsize1

        
        DO n = 1,this%atoms%ntype
         call radfun%generate_radial_functions(this%atoms, this%input, this%enpara, this%fmpi, this%vtot, n,usdus_out=usdus)
         !calculate the abc coefficients for this atom type and both spins
         call abc(1)%init(this%input, this%atoms, radfun%n_r, num_bands, n)
         call abc(1)%calc_abc(this%input, this%atoms, this%sym, this%cell, this%lapw, num_bands, usdus, this%noco, this%nococonv, 1, n, zMat(1))

         if (this%input%jspins == 2) then
            call abc(2)%init(this%input, this%atoms, radfun%n_r, num_bands, n)
            call abc(2)%calc_abc(this%input, this%atoms, this%sym, this%cell, this%lapw, num_bands, usdus, this%noco, this%nococonv, 2, n, zMat(2))
         else
            abc(2) = abc(1)
         end if
         DO j1=1,2
            j1_int = min(this%input%jspins,j1)
            DO i1=1,2
               i1_int = min(this%input%jspins,i1)
               DO j = 1,num_bands
                  DO na = 1, this%atoms%neq(n)
                     DO l = 1,this%atoms%lmax(n)
                        ll1 = l*(l+1)
                        DO m = -l,l
                           lm = ll1 + m
                           DO jcof=1,radfun%n_r(l)
                              cof_lm = CMPLX(0.,0.)
                              DO m1 = -l,l
                                 lm1 = ll1 + m1
                                 cof_lm = cof_lm + this%rsoc%soangl(l,m,i1,l,m1,j1)*conjg(abc(j1)%cof(j,lm1,jcof,na))
                              ENDDO
                              DO icof=1,radfun%n_r(l)
                                 if (abs(radfun%integral(icof,jcof,l,i1_int,j1_int))<1E-15) cycle
                                 DO i=1, num_bands
                                    mat(i1,j1)%data_c(i,j) = mat(i1,j1)%data_c(i,j) +&
                                    abc(i1)%cof(i,lm,icof,na)*this%rsoc%rso(icof,jcof,n,l,i1,j1)*cof_lm*radfun%integral(icof,jcof,l,i1_int,j1_int)
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

    END SUBROUTINE add_matrix_elements
END MODULE m_types_matelements_soc  
