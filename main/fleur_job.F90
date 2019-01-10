!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_jobs
    USE m_juDFT
    IMPLICIT NONE
    PRIVATE
    CHARACTER(LEN=30),PARAMETER:: NOT_A_JOBFILE=".__NOT__A__JOBFILE__"

    TYPE t_job
        INTEGER           :: PE_requested
        INTEGER           :: mpi_comm
        CHARACTER(LEN=20) :: directory
    END TYPE

    PUBLIC:: t_job,fleur_job_arguments,fleur_job_init,fleur_job_distribute,fleur_job_execute

CONTAINS
    SUBROUTINE fleur_job_single(jobs)
        TYPE(t_job),ALLOCATABLE,INTENT(OUT)::jobs(:)

        ALLOCATE(jobs(1))
        jobs(1)%PE_requested=0
        jobs(1)%directory="."
    END SUBROUTINE

    SUBROUTINE read_jobfile(jobs,file)
        TYPE(t_job),ALLOCATABLE,INTENT(OUT)::jobs(:)
        CHARACTER(LEN=*)::file

        LOGICAL:: l_file
        INTEGER:: njobs,i
        INQUIRE(FILE=file,EXIST=l_file)
        IF (l_file) THEN
            OPEN(99,FILE=file,STATUS="old")
        ELSE
            WRITE(*,*) "job input file not found"
            WRITE(*,*) "You specified an invalid filename:",file
            STOP "JOB FILE MISSING"
        ENDIF
        !Count the number of lines in job-file
        njobs=0
        DO
            READ(99,*,END=100)
            njobs=njobs+1
        ENDDO
100     REWIND(99)
        ALLOCATE(jobs(njobs))
        DO i=1,njobs
            READ(99,*) jobs(i)%PE_REQUESTED,jobs(i)%directory
        ENDDO
        CLOSE(99)
    END SUBROUTINE

    SUBROUTINE jobs_fromcommandline(jobs,no_jobs)
        TYPE(t_job),ALLOCATABLE,INTENT(INOUT)::jobs(:)
        INTEGER,INTENT(INOUT):: no_jobs

        TYPE(t_job),ALLOCATABLE ::jobs_tmp(:)
        INTEGER:: i
        CHARACTER(LEN=30)::str

        IF(ALLOCATED(jobs)) THEN
            no_jobs=size(jobs)+no_jobs
            ALLOCATE(jobs_tmp(size(jobs)))
            jobs_tmp=jobs
            DEALLOCATE(jobs)
            ALLOCATE(jobs(no_jobs))
            jobs(:size(jobs_tmp))=jobs_tmp
            no_jobs=size(jobs_tmp)+1
            DEALLOCATE(jobs_tmp)
        ELSE
            ALLOCATE(jobs(no_jobs))
            no_jobs=1
        ENDIF

        DO i=1,command_argument_count()
            CALL get_command_argument(i,str)
            IF(adjustl(str)=="-j") THEN
                CALL get_command_argument(i+1,str)
                IF (index(str,":")>1) THEN
                    READ(str(:index(str,":")-1),*) jobs(no_jobs)%PE_requested
                    jobs(no_jobs)%directory=str(index(str,":")+1:)
                    no_jobs=no_jobs+1
                ELSE
                    PRINT *,"Illegal job-description"
                    PRINT *,"You specified:",str
                    STOP "ILLEGAL DESCRIPTION"
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE

    SUBROUTINE jobs_on_commandine(jobfile,no_jobs)
        INTEGER,INTENT(OUT)::no_jobs
        CHARACTER(LEN=*)::jobfile

        INTEGER i
        CHARACTER(LEN=30)::str
        jobfile=NOT_A_JOBFILE
        no_jobs=0
        DO i=1,command_argument_count()
            CALL get_command_argument(i,str)
            IF(adjustl(str)=="-j") THEN
                no_jobs=no_jobs+1
            ENDIF
            IF (adjustl(str)=="-f") THEN
                CALL get_command_argument(i+1,jobfile)
            ENDIF
        ENDDO
    END SUBROUTINE

    SUBROUTINE fleur_job_arguments(jobs)
        TYPE(t_job),ALLOCATABLE,INTENT(OUT)::jobs(:)

        CHARACTER(LEN=30):: file
        INTEGER          :: no_jobs_commandline

        CALL jobs_on_commandine(file,no_jobs_commandline)
        IF (file.NE.NOT_A_JOBFILE)  &
            CALL read_jobfile(jobs,file)
        IF (no_jobs_commandline>0) &
            CALL jobs_fromcommandline(jobs,no_jobs_commandline)
        IF (.NOT.allocated(jobs)) &
            CALL fleur_job_single(jobs)
    END SUBROUTINE

    SUBROUTINE fleur_job_init()
      USE m_fleur_help
      use m_judft
        INTEGER:: irank=0
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
        INTEGER ierr(3), i
        CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,i,ierr)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
        IF(irank.EQ.0) THEN
           !$    IF (i<MPI_THREAD_SERIALIZED) THEN
           !$       WRITE(*,*) ""
           !$       WRITE(*,*) "Linked MPI version does not support OpenMP."
           !$       WRITE(*,*) ""
           !$       WRITE(*,*) "To solve this problem please do one of:"
           !$       WRITE(*,*) "   1. Link an adequate MPI version."
           !$       WRITE(*,*) "   2. Use fleur without MPI."
           !$       WRITE(*,*) "   3. Compile and use fleur without OpenMP."
           !$       WRITE(*,*) ""
           !$       CALL juDFT_error("MPI not usable with OpenMP")
           !$    END IF
           !Select the io-mode from the command-line
        END IF
#endif
        IF (irank==0) THEN
           CALL fleur_help()
        END IF
    END SUBROUTINE

    SUBROUTINE fleur_job_execute(jobs)
        USE m_fleur
        TYPE(t_job),INTENT(IN) ::jobs(:)

        INTEGER:: njob=1
        INTEGER:: irank=0

#ifdef CPP_MPI
        INTEGER:: ierr
      INCLUDE 'mpif.h'
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

        !find the number of the job for this PE
        DO njob=1,size(jobs)
            IF (jobs(njob)%mpi_comm==MPI_UNDEFINED) CYCLE
            CALL MPI_COMM_RANK(jobs(njob)%mpi_comm,irank,ierr)
            IF (irank.NE.MPI_UNDEFINED) EXIT
        ENDDO
#endif
        if (njob>size(jobs)) THEN
            print *, "GLOBAL-PE:",irank," does nothing"
            return
        endif
        !change directory
        CALL chdir(jobs(njob)%directory)
        !Call FLEUR

        CALL fleur_execute(jobs(njob)%mpi_comm)

    END SUBROUTINE

    SUBROUTINE fleur_job_distribute(jobs)
        TYPE(t_job),INTENT(INOUT)::jobs(:)
#ifdef CPP_MPI
      INCLUDE 'mpif.h'
        INTEGER:: i,free_pe,isize,irank,min_pe,new_comm,ierr

        CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,ierr)

        if (irank==0) print *,"Running on ",isize," PE"
        !First determine if there are PE that should be distributed automatically
        IF (any(jobs%pe_requested==0)) THEN
            i=count(jobs%pe_requested==0)
            free_pe=isize-sum(jobs%pe_requested)
            if (irank==0) print *,i," jobs are distributed on ",free_pe," unassigned PE"
            i=free_pe/i

            IF (i<1) THEN
                if (irank==0) PRINT *,"Not enough PE after automatic assignment of jobs"
                STOP "NOT enough PE"
            ELSE
                WHERE (jobs%pe_requested==0) jobs%pe_requested=i
            ENDIF
        ENDIF
        free_pe=isize-sum(jobs%pe_requested)
        IF (free_pe<0) THEN
            if (irank==0) PRINT *,"Not enough PE for assignment of jobs"
            STOP "NOT enough PE"
        ENDIF
        IF (free_pe>0.and.irank==0)    PRINT *,"WARNING, there are unused PE"

        !Now create the groups
        DO i=1,size(jobs)
            min_pe=sum(jobs(:i-1)%PE_requested)
            IF ((irank.GE.min_pe).AND.(irank<min_pe+jobs(i)%PE_requested)) EXIT
        ENDDO
        jobs%mpi_comm=MPI_UNDEFINED
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,i,irank,new_comm,ierr)
        IF (i.LE.size(jobs)) THEN
            PRINT* ,"PE:",irank," works on job ",i," in ",jobs(i)%directory
            jobs(i)%mpi_comm=new_comm
        ENDIF

#else
        IF (size(jobs)>1) THEN
            PRINT*, "Cannot run multiple jobs without MPI"
            STOP "NO MPI"
        ENDIF
        IF (sum(jobs%pe_requested)>1) THEN
            PRINT*, "You cannot request a multiple PE job without MPI"
            STOP "NO MPI"
        ENDIF
        jobs(1)%mpi_comm=1
#endif
    END SUBROUTINE
END MODULE

PROGRAM fleurjob
    USE m_fleur_jobs
    USE m_juDFT
    IMPLICIT NONE
    TYPE(t_job),ALLOCATABLE::jobs(:)
    CALL judft_init()
    CALL fleur_job_init()
    CALL fleur_job_arguments(jobs)
    CALL fleur_job_distribute(jobs)
    CALL fleur_job_execute(jobs)
END
