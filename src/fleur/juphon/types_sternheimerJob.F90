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

        logical  :: l_phonon = .false.
        logical  :: l_BEC = .false.
        logical  :: l_efield = .false.
        logical  :: l_efield_screened = .false.
        logical  :: l_IBScorrection = .false. ! Calculate incomplete Basis-Set-Corrections: Pulay + Surface 
        
        contains 

        procedure :: init => init_sternheimerJob

    end type t_sternheimerJob



    contains 



    subroutine init_sternheimerJob(this,fi,l_phonon,l_BEC,l_efield)

        use m_types
        
        class(t_sternheimerJob),intent(inout) :: this
        type(t_fleurinput),intent(in) :: fi 
        logical,optional,intent(in)   :: l_phonon 
        logical,optional,intent(in)   :: l_BEC 
        logical,optional,intent(in)   :: l_efield 


        integer :: jobSize , iJob , iQ , iDtype, iDir
        ! we currently work with logicals to circumvent circular dependence
        ! it would be nicer to work with select type(t_dfpt) 
        logical :: l_ph = .false. 
        logical :: l_borncharges  = .false. 
        logical :: l_efield_pert  = .false. 

        if (present(l_phonon)) l_ph = l_phonon
        if (present(l_BEC)) l_borncharges = l_BEC
        if (present(l_efield)) l_efield_pert = l_efield

        if (allocated(this%iJobList)) deallocate(this%iJobList)
        if (allocated(this%iQList)) deallocate(this%iQList)
        if (allocated(this%iDirList)) deallocate(this%iDirList)
        if (allocated(this%iDtypeList)) deallocate(this%iDtypeList)
        if (allocated(this%needs_postprocessing)) deallocate(this%needs_postprocessing)



        this%l_phonon      = .false.
        this%l_efield      = .false.
        this%l_BEC         = .false.

        if (l_ph) then 
            this%l_phonon = .true.
            this%l_IBScorrection = .true. 

            jobSize = 3 * fi%atoms%nat * size(fi%juPhon%qvec,2)      

            allocate(this%iJobList(jobSize))
            allocate(this%iQList(jobSize))
            allocate(this%iDirList(jobSize))
            allocate(this%iDtypeList(jobSize))
            allocate(this%needs_postprocessing(jobSize))

            iJob = 1 
            this%needs_postprocessing(:) = .false. 
            
            do iQ = 1 , size(fi%juPhon%qvec,2)      
                do iDtype = 1 , fi%atoms%nat
                    do iDir = 1 , 3 
                        this%iJobList(iJob)   = iJob
                        this%iQList(iJob)     = iQ
                        this%iDtypeList(iJob) = iDtype 
                        this%iDirList(iJob)   = iDir 
                        iJob = iJob + 1 
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

            iJob = 1 
            this%needs_postprocessing(:) = .false. 
            
            do iDir = 1 , 3 
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

            iJob = 1 
            this%needs_postprocessing(:) = .false. 
            
            do iQ = 1 , 3 
                do iDtype = 1 , fi%atoms%nat
                    do iDir = 1 , 3 
                        this%iJobList(iJob)   = iJob
                        this%iQList(iJob)     = iQ
                        this%iDtypeList(iJob) = iDtype 
                        this%iDirList(iJob)   = iDir 
                        iJob = iJob + 1
                    end do !iDir
                end do !iDtype 
            end do !iQ
        end if 
 

        do iJob = 1 , jobSize
            write(114,*) "jobNumber" , iJob
            write(114,*) "iQ" , this%iQList(iJob)
            write(114,*) "iDtype" , this%iDtypeList(iJob)
            write(114,*) "iDir" , this%iDirList(iJob)
            write(114,*) "postprocessing?" , this%needs_postprocessing(iJob)

        end do 

    end subroutine init_sternheimerJob 

end module m_types_sternheimerJob