!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


module m_types_dfpt

    use m_judft
    implicit none
    
    private
    public :: t_dfpt

    type,abstract :: t_dfpt

        integer :: iDtypeLoop = 1 ! perform scf loop over atom 
        logical :: l_phonon = .true.  ! phonon type perturbation  
        real, allocatable :: qVectors(:,:)

        contains 

        procedure :: init => init_scf
        procedure :: perform_scf 
        procedure(init_child), deferred :: init_child 
        procedure(postprocessing_scf), deferred :: postprocessing_scf 
        procedure(postprocessing_qpoint), deferred :: postprocessing_qpoint 
        procedure(q_indepent_properties), deferred :: q_indepent_properties 
                
    end type 


    interface 
        subroutine init_child(this,fi,nqpts)
            import t_dfpt
            use m_types
            class(t_dfpt),intent(inout)    :: this
            type(t_fleurinput), intent(in) :: fi 
            integer, intent(in)            :: nqpts
        end subroutine init_child
    end interface 

    interface 
        subroutine q_indepent_properties(this,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
            import t_dfpt
            use m_types
            class(t_dfpt),intent(inout)   :: this       
            type(t_fleurinput), intent(in):: fi 
            type(t_mpi), intent(in)       :: fmpi
            type(t_stars),intent(in)      :: stars
            type(t_sphhar),intent(in)     :: sphhar
            class(t_xcpot),intent(in)     :: xcpot
            type(t_nococonv),intent(in)   :: nococonv
            type(t_hybdat),intent(inout)  :: hybdat
            type(t_potden),intent(in)     :: rho 
            type(t_potden),intent(in)     :: vTot
            type(t_potden), intent(in)    :: grRho3(3), grVtot3(3), grVC3(3), grVext3(3),grgrVext3x3(3,3)
        end subroutine q_indepent_properties
    end interface 


    interface 
        subroutine postprocessing_scf(this,fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,rho,vTot,grRho3,grVext3,grVc3,den1,vTot1,den1Im,vTot1Im,vC1)
            import t_dfpt
            use m_types
            class(t_dfpt),intent(inout)     :: this
            type(t_fleurinput), intent(in)  :: fi 
            type(t_stars),intent(in)        :: stars
            type(t_stars),intent(in)        :: starsq
            type(t_sphhar),intent(in)       :: sphhar
            class(t_xcpot),intent(in)       :: xcpot
            type(t_nococonv),intent(in)     :: nococonv
            type(t_hybdat),intent(inout)    :: hybdat
            type(t_mpi), intent(in)         :: fmpi
            type(t_kpts), intent(in)        :: qpts
            integer, allocatable,intent(in) :: q_list(:)
            integer, intent(in)             :: iQ,iDtype,iDir,eig_id,dfpt_eig_id,dfpt_eig_id2
            type(t_enpara),intent(inout)    :: enpara
            type(t_results),intent(inout)   :: results, results1
            logical,intent(in)              :: l_real
            type(t_potden),intent(in)       :: rho,vTot,grRho3(3),grVext3,grVC3,den1,vTot1,den1Im,vTot1Im,vC1         
        end subroutine postprocessing_scf
    end interface 


    interface 
        subroutine postprocessing_qpoint(this,fi,fmpi,qpts,iQ)
            import t_dfpt
            use m_types
            class(t_dfpt),intent(inout)   :: this         
            type(t_fleurinput),intent(in) :: fi
            type(t_mpi),intent(in)        :: fmpi
            type(t_kpts),intent(in)       :: qpts
            integer,intent(in)            :: iQ
        end subroutine postprocessing_qpoint
    end interface 

 
    subroutine init_scf(this,l_phonon,qvec)

        use m_types_fleurinput

        type(t_dfpt), intent(inout) :: this
        logical, intent(in) :: l_phonon 
        real, allocatable, intent(in) :: qvec(:,:)
        integer :: nqpts

        ! set up necessary input that seperates the different types of scf calculations 

        this%l_phonon = l_phonon
        if (this%l_phonon) this%iDtypeLoop = fi%atoms%nat

        allocate(this%qVectors,mold=qvec)
        this%qVectors = qvec

        nqpts = size(qvec,2)

        this%init_child(fi,nqpts)

    end subroutine


    subroutine perform_scf(this,fi,fmpi,stars,sphhar,xcpot,forcetheo,enpara,nococonv,hybdat,rho,vTot,vxc,&
                            results,resultsq, results1, eig_id,q_eig_id, & 
                            dfpt_eig_id,dfpt_eig_id2,l_minusq,resultsqm,results1m,qm_eig_id, dfpt_eigm_id, dfpt_eigm_id2)
        use m_eigen 
        use m_fermie 
        use m_dfpt_sternheimer
        use m_dfpt_generate_gradient
        use m_types

        type(t_dfpt), intent(inout) :: this 
        type(t_fleurinput), intent(in)  :: fi 
        type(t_mpi), intent(in)         :: fmpi
        type(t_stars),intent(in)      :: stars
        type(t_sphhar),intent(in)     :: sphhar
        class(t_xcpot),intent(in)     :: xcpot
        type(t_forcetheo),intent(in)  :: forcetheo
        type(t_enpara),intent(inout)     :: enpara
        type(t_nococonv),intent(in)   :: nococonv
        type(t_hybdat),intent(inout)     :: hybdat
        type(t_potden),intent(in)     :: rho 
        type(t_potden),intent(in)     :: vTot
        type(t_potden),intent(in)     :: vxc
        type(t_results),intent(inout)       :: results
        type(t_results),intent(inout)       :: resultsq
        type(t_results),intent(inout)       :: resultsqm
        type(t_results),intent(inout)       :: results1 , results1m
        integer, intent(in) :: eig_id, q_eig_id, dfpt_eig_id,dfpt_eig_id2,qm_eig_id,dfpt_eigm_id,dfpt_eigm_id2
        logical, intent(in) :: l_minusq

        ! gradient variables 
        type(t_potden)   :: grRho3(3),grVtot3(3),grVc3(3),grVext3(3),grgrVext3x3(3,3)

        ! scf variables 
        type(t_potden) :: den1, den1Im, vTot1, vTot1Im, vC1,vC1Im , vTot1m, vTot1mIm
        type(t_kpts) :: kqpts, kmqpts , qpts 
        type(t_stars)  :: starsq, starsmq

        integer, allocatable :: q_list(:)

#ifdef CPP_MPI
        integer :: ierr
#endif 

        l_real = fi%sym%invs.and.(.not.fi%noco%l_soc).and.(.not.fi%noco%l_noco).and.fi%atoms%n_hia==0


        ! calculate the gradients as this is currently a necessary input for the sternheimer in any perturbation 

        ! Generate the gradients of the density and the various potentials, that will be used at different points in the programm.
        ! The density gradient is calculated by numerical differentiation, while the potential gradients are constructed (from the
        ! density gradient) by a Weinert construction, just like the potentials are from the density.
        ! This is done to ensure good continuity.
        call timestart("Gradient generation")
        call dfpt_generate_gradient(fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)
        call timestop("Gradient generation")
        call dfpt_generate_gradient() 

        ! create a kpts type that contains the necessary q-vectors 
        qpts = fi%kpts
        deallocate(qpts%bk)
        allocate(qpts%bk,mold=this%%qVectors) ! this is not nice. Maybe change the expected type in dfpt_sternheimer/make_stars
        qpts%bk(:, :size(this%qVectors,2)) = this%qVectors
        allocate(q_list(size(this%qVectors,2)))  
        q_list = (/(iArray, iArray=1,SIZE(this%qVectors,2), 1)/)

        q_start = fi%juPhon%startq
        q_stop = merge(fi%juPhon%stopq,size(q_list),fi%juPhon%stopq/=0)

        call this%q_indepent_properties(fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVC3,grVext3,grgrVext3x3)

        do iQ = q_start, q_stop
            call timestart("q-Point")

            l_gamma =.FALSE.
            if  (sum(abs(qpts%bk(:,q_list(iQ))))== 0.0) then
                l_gamma =.TRUE.
            end if 
            
            kqpts = fi%kpts
            do ikpt = 1, fi%kpts%nkpt
                kqpts%bk(:, ikpt) = kqpts%bk(:, ikpt) + qpts%bk(:,q_list(iQ))
            end do 

            if (l_minusq) then 
                kmqpts = fi%kpts
                do ikpt = 1, fi%kpts%nkpt
                    kmqpts%bk(:, ikpt) = kmqpts%bk(:, ikpt) - qpts%bk(:,q_list(iQ))
                end do 
            end if 

            call timestart("Eigenstuff at k+q")
            ! Get the eigenstuff at k+q
            call resultsq%reset_results(fi%input)

            call eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv,  &
                    hybdat, 1, q_eig_id, resultsq, rho, vTot, vxc, hub1data, &
                    qpts%bk(:,q_list(iQ)))

            ! Fermi level and occupancies
            call timestart("determination of fermi energy")
            call fermie(q_eig_id, fmpi, kqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, resultsq)
            call timestop("determination of fermi energy")

#ifdef CPP_MPI
            call MPI_BCAST(resultsq%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
            call MPI_BCAST(resultsq%w_iks, size(resultsq%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif
            call timestop("Eigenstuff at k+q")


            if (l_minusq) then
                call timestart("Eigenstuff at k-q")
                ! Solve eigenvalue equation on k-q grid
                call resultsqm%reset_results(fi%input)

                call eigen(fi, fmpi, stars, sphhar, xcpot, forcetheo, enpara, nococonv,  &
                        hybdat, 1, qm_eig_id, resultsqm, rho, vTot, vxc, hub1data, &
                        -qpts%bk(:,q_list(iQ)))

                ! Fermi level and occupancies
                call timestart("determination of fermi energy")
                call fermie(qm_eig_id, fmpi, kmqpts, fi%input, fi%noco, enpara%epara_min, fi%cell, resultsqm)
                call timestop("determination of fermi energy")

#ifdef CPP_MPI
                call MPI_BCAST(resultsqm%ef, 1, MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
                call MPI_BCAST(resultsqm%w_iks, SIZE(resultsqm%w_iks), MPI_DOUBLE_PRECISION, 0, fmpi%mpi_comm, ierr)
#endif

                call timestop("Eigenstuff at k-q")
            end if


            do iDtype = 1 , this%iDtypeLoop ! change this to local variable that is called by a get function
                call timestart("Typeloop")
                    do iDir = 1 , 3 
                        call timestart("Dirloop")
                        dfpt_tag = ''
                        write(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir

                        if (fmpi%irank==0) then 
                            write(*,*) 'Starting calculation for:'
                            write(*,*) ' q         = ', fi%juPhon%qvec(:,q_list(iQ))
                            write(*,*) ' atom      = ', iDtype
                            write(*,*) ' direction = ', iDir
                        end if 

                        ! reset arrays 
                        call starsq%reset_stars()
                        call den1%reset_dfpt()
                        call den1Im%reset_dfpt()
                        call vTot1%reset_dfpt()
                        call vTot1Im%reset_dfpt()
                        call vC1%reset_dfpt()
                        call vC1Im%reset_dfpt()
                        call results1%reset_results(fi%input)
                        if (l_minusq) then 
                            call starsmq%reset_stars()
                            call vTot1m%reset_dfpt()
                            call vTot1mIm%reset_dfpt()
                            call results1m%reset_results(fi%input)
                            ! I am unsure why there is no vC1m 
                        end if 
                            
                        if (fmpi%irank==0) write(*,*) '-------------------------'
                        ! Legacy Comment: Dont dare to delete it
                        ! This is where the magic happens. The Sternheimer equation is solved
                        ! iteratively, providing the scf part of dfpt calculations.
                        if (l_minusq) then 
                            call timestart("Sternheimer with -q")
                            call dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, &
                                                rho, vTot, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                                dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im, &
                                                starsmq, resultsqm, dfpt_eigm_id, dfpt_eigm_id2, qm_eig_id, results1m, vTot1m, vTot1mIm)
                            call timestop("Sternheimer with -q")
                        else
                            call timestart("Sternheimer")
                            call dfpt_sternheimer(fi, xcpot, sphhar, stars, starsq, nococonv, qpts, fmpi, results, resultsq, enpara, hybdat, &
                                                rho, vTot, grRho3(iDir), grVtot3(iDir), grVext3(iDir), q_list(iQ), iDtype, iDir, &
                                                dfpt_tag, eig_id, l_real, results1, dfpt_eig_id, dfpt_eig_id2, q_eig_id, &
                                                den1, vTot1, den1Im, vTot1Im, vC1, vC1Im)
                            call timestop("Sternheimer")
                        end if 

                        
                        call this%postprocessing_scf(fi,stars,starsq,sphhar,xcpot,nococonv,hybdat,fmpi,qpts,q_list,iQ,iDtype,iDir,eig_id,dfpt_eig_id, &
                                          dfpt_eig_id2,enpara,results,results1,l_real,rho,vTot,grRho3,grVext3,grVc3,den1,vTot1,den1Im,vTot1Im,vC1)

                    end do ! iDir
                call timestop("Typeloop")
            end do ! iDtype
            call this%postprocessing_qpoint(fi,fmpi,qpts,iQ)
        end do ! iQ

    end subroutine perform_scf
    

end module m_types_dfpt