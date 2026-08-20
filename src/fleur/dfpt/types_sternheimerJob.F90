!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------



module m_types_sternheimerJob

    use m_juDFT 
    
    implicit none 

    private 
    public :: t_sternheimerJob


    type t_sternheimerJob

        ! all will have the size nJob  
        integer,allocatable :: iJobList(:)     
        integer,allocatable :: iQList(:) 
        integer,allocatable :: iDirList(:)
        integer,allocatable :: iDtypeList(:)
        logical,allocatable :: needs_postprocessing(:)
        logical,allocatable :: needs_eigen(:)

        ! List for task the MPI process has to perform 
        integer,allocatable :: JobListMPI(:)     


        logical  :: l_phonon = .false.
        logical  :: l_BEC = .false.
        logical  :: l_efield = .false.
        logical  :: l_bfield = .false.
        logical  :: l_efield_screened = .false.
        logical  :: l_IBScorrection = .false. ! Calculate incomplete Basis-Set-Corrections: Pulay + Surface 
        
        contains 

        procedure :: init => init_sternheimerJob

    end type t_sternheimerJob



    contains 



    subroutine init_sternheimerJob(this,fi,l_phonon,l_BEC,l_efield,l_bfield)

        use m_types_fleurinput
        
        class(t_sternheimerJob),intent(inout) :: this
        type(t_fleurinput),intent(in) :: fi 
        logical,optional,intent(in)   :: l_phonon 
        logical,optional,intent(in)   :: l_BEC 
        logical,optional,intent(in)   :: l_efield 
        logical,optional,intent(in)   :: l_bfield 


        integer :: jobSize , iJob , iQ , iDtype, iDir, iQ_last
        integer :: q_start, q_stop
        ! we currently work with logicals to circumvent circular dependence
        ! it would be nicer to work with select type(t_dfpt_scf) 
        logical :: l_ph = .false. 
        logical :: l_borncharges  = .false. 
        logical :: l_efield_pert  = .false. 
        logical :: l_bfield_pert  = .false. 

        if (present(l_phonon)) l_ph = l_phonon
        if (present(l_BEC)) l_borncharges = l_BEC
        if (present(l_efield)) l_efield_pert = l_efield
        if (present(l_bfield)) l_bfield_pert = l_bfield

        if (allocated(this%iJobList)) deallocate(this%iJobList)
        if (allocated(this%iQList)) deallocate(this%iQList)
        if (allocated(this%iDirList)) deallocate(this%iDirList)
        if (allocated(this%iDtypeList)) deallocate(this%iDtypeList)
        if (allocated(this%needs_postprocessing)) deallocate(this%needs_postprocessing)
        if (allocated(this%needs_eigen)) deallocate(this%needs_eigen)



        this%l_phonon      = .false.
        this%l_efield      = .false.
        this%l_BEC         = .false.
        this%l_bfield      = .false.

        if (l_ph) then 
            this%l_phonon = .true.
            this%l_IBScorrection = .true. 

            ! Introduce the option to calculate a fraction of the input 
            ! Useful for restarting a calculation
            q_start = fi%dfpt%startq
            q_stop = merge(fi%dfpt%stopq,size(fi%dfpt%qvec,2),fi%dfpt%stopq/=0)

            jobSize = 3 * fi%atoms%nat * (q_stop - q_start + 1 ) 

            allocate(this%iJobList(jobSize))
            allocate(this%iQList(jobSize))
            allocate(this%iDirList(jobSize))
            allocate(this%iDtypeList(jobSize))
            allocate(this%needs_postprocessing(jobSize))
            allocate(this%needs_eigen(jobSize))

            q_start = fi%dfpt%startq
            q_stop = merge(fi%dfpt%stopq,size(fi%dfpt%qvec,2),fi%dfpt%stopq/=0)


            iJob = 1 
            iQ_last = 1 
            this%needs_postprocessing(:) = .false. 
            this%needs_eigen(:) = .false. 
            ! set first iJob to true 
            this%needs_eigen(1) = .true. 


            do iQ = q_start , q_stop
                do iDtype = 1 , fi%atoms%nat
                    do iDir = 1 , 3 
                        this%iJobList(iJob)   = iJob
                        this%iQList(iJob)     = iQ
                        this%iDtypeList(iJob) = iDtype 
                        this%iDirList(iJob)   = iDir 
                        if (iQ .ne. iQ_last) this%needs_eigen(iJob) = .true.
                        iJob = iJob + 1
                        iQ_last = iQ  
                    end do ! iDir 
                end do !iDtype
                this%needs_postprocessing(iJob-1) = .true.
            end do !iQ 
        end if 

        if (l_efield_pert) then 

            this%l_efield = .true. 

            jobSize = 3 ! all 3 cartesian directions 

            allocate(this%iJobList(jobSize))
            allocate(this%iQList(jobSize))
            allocate(this%iDirList(jobSize))
            allocate(this%iDtypeList(jobSize))
            allocate(this%needs_postprocessing(jobSize))
            allocate(this%needs_eigen(jobSize))

            iJob = 1 
            this%needs_postprocessing(:) = .false. 
            this%needs_eigen(:) = .true. 
            
            do iDir = 1 , 3 
                this%iJobList(iJob)   = iJob
                this%iQList(iJob)     = iDir
                this%iDtypeList(iJob) = 1 
                this%iDirList(iJob)   = iDir 
                iJob = iJob + 1 
            end do !iDir 
        end if 

        if (l_bfield_pert) then 

            this%l_bfield = .true. 

            jobSize = 1 ! all 3 cartesian directions 

            allocate(this%iJobList(jobSize))
            allocate(this%iQList(jobSize))
            allocate(this%iDirList(jobSize))
            allocate(this%iDtypeList(jobSize))
            allocate(this%needs_postprocessing(jobSize))
            allocate(this%needs_eigen(jobSize))

            iJob = 1 
            this%needs_postprocessing(:) = .false. 
            this%needs_eigen(:) = .true. 
            
            do iDir = 1 , 1 !WIP
                this%iJobList(iJob)   = iJob
                this%iQList(iJob)     = iDir
                this%iDtypeList(iJob) = 1 
                this%iDirList(iJob)   = iDir 
                iJob = iJob + 1 
            end do !iDir 
        end if 

        if (l_borncharges) then 
            this%l_BEC           = .true. 
            this%l_IBScorrection = .true. 

            jobSize = 3 * 3 * fi%atoms%nat      

            allocate(this%iJobList(jobSize))
            allocate(this%iQList(jobSize))
            allocate(this%iDirList(jobSize))
            allocate(this%iDtypeList(jobSize))
            allocate(this%needs_postprocessing(jobSize))
            allocate(this%needs_eigen(jobSize))


            iJob = 1 
            iQ_last = 1 
            this%needs_postprocessing(:) = .false. 
            this%needs_eigen(:) = .false. 
            this%needs_eigen(1) = .true. 

            do iQ = 1 , 3 
                do iDtype = 1 , fi%atoms%nat
                    do iDir = 1 , 3 
                        this%iJobList(iJob)   = iJob
                        this%iQList(iJob)     = iQ
                        this%iDtypeList(iJob) = iDtype 
                        this%iDirList(iJob)   = iDir 
                        if (iQ .ne. iQ_last)  this%needs_eigen(iJob) = .true. 
                        iJob = iJob + 1
                        iQ_last = iQ
                    end do !iDir
                end do !iDtype 
            end do !iQ
        end if 
 

        ! This is a placeholder for future mpi parallelisation 
        ! Please be careful when distributing the jobs
        ! Carefully think about the distribution of the q points
        ! As the postprocessing_q_point needs the full information of 
        ! Sternheimer scf for that specific q point. This could lead to long
        ! waiting times for some iranks. 
        allocate(this%JobListMPI(jobSize))
        do iJob = 1 , size(this%JobListMPI)
            this%JobListMPI(iJob) =  iJob
        end do 

    end subroutine init_sternheimerJob 

end module m_types_sternheimerJob