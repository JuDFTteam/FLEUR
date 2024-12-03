!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_elph_linewidth

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT
    USE m_types 
    USE m_constants


    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_ph_linewidth(fi,fmpi,qpts,results,resultsq,results1,eigenVals,gmat,iQ,nbasfcnq_min,nuWindow,ph_linewidth)
        ! This subroutine calculates the phonon linewdith
        ! Currently implemented is the linewidth calcualtion with smearing 

        USE m_dosbin
        USE m_smooth
        USE m_dfpt_fermie, ONLY : sfermi
        USE m_dfpt_tetra_double
        !USE m_intgr, ONLY : intgz0
        USE m_smooth
        USE m_dfpt_tetra_single

        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_mpi), INTENT(IN) :: fmpi
        TYPE(t_kpts), INTENT(IN) :: qpts
        TYPE(t_results), INTENT(IN) :: results,resultsq, results1
        REAL,ALLOCATABLE,INTENT(IN) :: eigenVals(:)
        COMPLEX,ALLOCATABLE,INTENT(INOUT) :: gmat(:,:,:,:,:)
        INTEGER, INTENT(IN) :: iQ
        INTEGER, INTENT(IN) :: nbasfcnq_min
        INTEGER, INTENT(IN) :: nuWindow(2,2)
        REAL, ALLOCATABLE, INTENT(OUT)    :: ph_linewidth(:)

        INTEGER :: iNupr,nu, ispin , gridPoint , nk , nk_i ,  nZero , iMode , noccbd, ind, indPr
        REAL :: emin, emax , x ,xq , allowed
        REAL, ALLOCATABLE :: eGrid(:),linewidth(:,:) , kInt(:,:)
        COMPLEX,ALLOCATABLE:: kInt_gmat(:,:,:,:) !(nu',nu,jsp,normal_mode)
        REAL,ALLOCATABLE :: gmatCollect(:,:,:,:,:)

        REAL , ALLOCATABLE :: gauss(:) 
        REAL :: intOut
#ifdef CPP_MPI
        INTEGER :: ierr
#endif 


        gmat(:,:,:,:,:) = ABS(gmat(:,:,:,:,:))**2 
        

        SELECT CASE(fi%juphon%i_integration)

            ! 1 = Full calculation of linewidth with one delta function via Gaussian Smearing 
            ! 2 = Double delta approximation of linewidth calculated with Gaussian Smearing
            ! 3 = Double delta approximation of linewidth with Tetra Method 


            !   For the different integration methods we need different sizes of the gmat 
            !   For the linewidth in its original formulation we need to multiply the fermi functions
            !   In order to not get combinatinos of nu -> unocc , nu' --> unocc we have to set these terms 0
            !   We always scatter from a ground state property into a perturbed state 
            !   For the tetragonal method it could happen that at a tetra we have two occupied and two unoccupied at nu
            !   In order to get the correct interpolation we also need the contribution at the unoccupied states from nu'

            CASE(1)
                
                ! Cut out all contributions coming from nu-> unoccuppied 
                !DO ispin = 1 , fi%input%jspins
                !    DO nk = 1 , fi%kpts%nkpt
                !        noccbd  = COUNT(results%w_iks(:,nk,ispin)*2.0/fi%input%jspins>1.e-8)
                !        DO nu = 1 , size(gmat,2)
                !            allowed = 1. 
                !            IF (nu .GT. noccbd) allowed = 0.  
                !            gmat(:,nu,nk,ispin,:) = allowed * gmat(:,nu,nk,ispin,:) 
                !        END DO  ! nu 
                !    END DO ! nk 
                !END DO !ispin 

                ! mutliply with fermi function 
                IF (fi%juphon%i_integration == 1 ) THEN 
                    DO ispin = 1 , fi%input%jspins
                        DO nk_i = 1 , size(fmpi%k_list)
                            nk = fmpi%k_list(nk_i)
                            !noccbd  = COUNT(results%w_iks(:,nk,ispin)*2.0/fi%input%jspins>1.e-8)
                            DO nu = nuWindow(1,1) , nuWindow(1,2)  
                                ind =  nu - nuWindow(1,1) + 1 
                                x = (results%eig(nu,nk,ispin)-results%ef)/fi%input%tkb
                                DO iNupr = nuWindow(2,1) , nuwindow(2,2)
                                    indPr = iNupr- nuWindow(2,1) + 1 
                                    xq = (resultsq%eig(iNupr,nk,ispin)-results%ef)/fi%input%tkb
                                    gmat(indPr,ind,nk_i,ispin,:) = gmat(indPr,ind,nk_i,ispin,:)*(sfermi(x) - sfermi(xq))
                                END DO ! iNupr 
                            END DO  ! nu 
                        END DO ! nk 
                    END DO 
                END IF 

 
                eMin = - 8 * fi%input%tkb
                eMax =   8 * fi%input%tkb
                ALLOCATE(linewidth(fi%banddos%ndos_points,fi%input%jspins))  
                ALLOCATE(eGrid(fi%banddos%ndos_points))  
                ALLOCATE(ph_linewidth(3*fi%atoms%ntype))
                ALLOCATE(gauss(fi%banddos%ndos_points))
                gauss = 0.0 
                linewidth = 0.0 
                ph_linewidth = 0.0 
                eGrid = 0.0 

                
                DO gridPoint=1,fi%banddos%ndos_points
                    eGrid(gridPoint)=emin+(emax-emin)/(fi%banddos%ndos_points-1.0)*(gridPoint-1.0)
                    IF (ABS(eGrid(gridPoint)) .LT. 1e-8) nZero = gridPoint 
                    !gauss(gridPoint) = 1/SQRT(tpi_const*fi%juphon%smearingGauss**2) *EXP( -0.5*eGrid(gridPoint)**2/(fi%juphon%smearingGauss)**2)
                END DO

                DO iMode = 1 , 3*fi%atoms%ntype
                    IF (eigenVals(iMode) .GE. 0.0 ) THEN 
                        
                        ! If omega becomes negative the deltra distribution is never satisfied as IM(eig,eigq) = 0.0 
                        CALL dos_bin_transport(fmpi,fi%input%jspins,fi%kpts%wtkpt,eGrid,results%eig(nuWindow(1,1):nuWindow(1,2),:,:)  &
                        &                      ,resultsq%eig(nuWindow(2,1):nuWindow(2,2),:,:), REAL(gmat(:,:,:,:,iMode)), linewidth, -SQRT(eigenVals(iMode)))


#ifdef CPP_MPI
                        CALL MPI_ALLREDUCE(MPI_IN_PLACE,linewidth,size(linewidth),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
                        CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#endif 
                        IF (fmpi%irank == 0 ) THEN 
                            DO ispin = 1 , fi%input%jspins
                                CALL smooth(eGrid,linewidth(:,ispin),fi%juPhon%smearingGauss,size(eGrid))
                            END DO 


                            DO ispin = 1 , fi%input%jspins

                                ph_linewidth(iMode) =  ph_linewidth(iMode) +  tpi_const * linewidth(nZero,ispin)
                            END DO 
                        END IF 
                    ELSE
                        IF (fmpi%irank ==0 ) THEN 
                            write(*,*) '-------------------------'
                            write(*,*) 'linewidth: Eigenvalue imaginary --> Phonon linewidth set to zero'
                            write(*,*) '-------------------------'
                        END IF 
                    END IF 
                END DO 
            

            CASE(2)

                ! Cut out all contributions coming from nu-> unoccuppied 
                DO ispin = 1 , fi%input%jspins
                    DO nk = 1 , fi%kpts%nkpt
                        noccbd  = COUNT(results%w_iks(:,nk,ispin)*2.0/fi%input%jspins>1.e-8)
                        DO nu = 1 , size(gmat,2)
                            allowed = 1. 
                            IF (nu .GT. noccbd) allowed = 0.  
                            gmat(:,nu,nk,ispin,:) = allowed * gmat(:,nu,nk,ispin,:) 
                        END DO  ! nu 
                    END DO ! nk 
                END DO !ispin

                eMin = - 4 * fi%input%tkb
                eMax =   4 * fi%input%tkb
                ALLOCATE(linewidth(fi%banddos%ndos_points,fi%input%jspins))  
                ALLOCATE(eGrid(fi%banddos%ndos_points))  
                ALLOCATE(ph_linewidth(3*fi%atoms%ntype))
                ALLOCATE(gauss(fi%banddos%ndos_points))
                gauss = 0.0 
                linewidth = 0.0 
                ph_linewidth = 0.0 
                eGrid = 0.0


                DO gridPoint=1,fi%banddos%ndos_points
                    eGrid(gridPoint)=emin+(emax-emin)/(fi%banddos%ndos_points-1.0)*(gridPoint-1.0)
                    IF (ABS(eGrid(gridPoint)) .LT. 1e-8) nZero = gridPoint 
                END DO

                DO iMode = 1 , 3*fi%atoms%ntype
                    IF (eigenVals(iMode) .GE. 0.0 ) THEN 
                        ! Gaussian Smearing is already done in routine, in case of double delta distribution
                        ! Warning still needs to be tested 
                        CALL dos_bin_double(fi%input%jspins,fi%kpts%wtkpt,eGrid,results%eig(:size(gmat,2),:,:)  &
                        &     ,resultsq%eig(:nbasfcnq_min,:,:), REAL(gmat(:nbasfcnq_min,:,:,:,iMode)), fi%juphon%smearingGauss ,linewidth, results%ef)

                        DO ispin = 1 , fi%input%jspins
                            ! factor two for spin deg. is calculated in dos_bin 
                            ph_linewidth(iMode) =  ph_linewidth(iMode) +  tpi_const * SQRT(eigenVals(iMode)) * linewidth(nZero,ispin)
                        END DO 
                    ELSE
                        write(*,*) '-------------------------'
                        write(*,*) 'linewidth: Eigenvalue imaginary --> Phonon linewidth set to zero'
                        write(*,*) '-------------------------'
                    END IF 
                END DO 
            


            CASE(3)

                CALL timestart("k-Integration el-ph")
                CALL dfpt_tetra_double(fi,results,resultsq, results1, gmat,nbasfcnq_min,kInt_gmat)
                CALL timestop("k-Integration el-ph")

                DEALLOCATE(gmat)
                ALLOCATE(ph_linewidth( 3*fi%atoms%ntype))
                
                ph_linewidth = 0.0 
                DO iMode = 1 , 3*fi%atoms%ntype
                    DO ispin = 1 , fi%input%jspins
                        DO nu = 1 , size(kInt_gmat,2)
                            DO iNupr = 1 , nbasfcnq_min
                                ph_linewidth(iMode) =  ph_linewidth(iMode) + tpi_const* 2/fi%input%jspins * SQRT(eigenVals(iMode))*REAL(kInt_gmat(iNupr,nu,ispin,iMode))
                            END DO 
                        END DO 
                    END DO 
                END DO 
                
            CASE(4)
                ALLOCATE(linewidth(3*fi%atoms%ntype,fi%input%jspins))
                ALLOCATE(ph_linewidth(3*fi%atoms%ntype))
                ALLOCATE(gmatCollect(size(gmat,1),size(gmat,2),fi%kpts%nkpt,size(gmat,4),size(gmat,5)))

                ph_linewidth = 0.0 
                linewidth = 0.0 
                gmatCollect = 0.0
                
                
                DO ispin = 1 , fi%input%jspins
                    DO nk_i = 1 , size(fmpi%k_list)
                        nk = fmpi%k_list(nk_i)
                        DO nu = nuWindow(1,1) , nuWindow(1,2)  
                            ind =  nu - nuWindow(1,1) + 1 
                            x = (results%eig(nu,nk,ispin)-results%ef)/fi%input%tkb
                            DO iNupr = nuWindow(2,1) , nuwindow(2,2)
                                indPr = iNupr- nuWindow(2,1) + 1 
                                xq = (resultsq%eig(iNupr,nk,ispin)-results%ef)/fi%input%tkb
                                gmat(indPr,ind,nk_i,ispin,:) = gmat(indPr,ind,nk_i,ispin,:)*(sfermi(x) - sfermi(xq))
                            END DO ! iNupr 
                        END DO  ! nu 
                    END DO ! nk 
                END DO 


                DO nk_i = 1 , size(fmpi%k_list)
                    nk = fmpi%k_list(nk_i)
                    gmatCollect(:,:,nk,:,:) = REAL(gmat(:,:,nk_i,:,:))
                END DO 
#ifdef CPP_MPI
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,gmatCollect,size(gmatCollect),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
                CALL MPI_BARRIER(fmpi%MPI_COMM,ierr)
#else 
                gmatCollect = REAL(gmat)         
#endif 
                !DEALLOCATE(gmat)
                IF (fmpi%irank==0) THEN 
                    CALL timestart("k-Integration el-ph")
                    CALL dfpt_tetra_single(fi,results,resultsq,results1,gmatCollect,eigenVals,nuWindow,linewidth)
                    CALL timestop("k-Integration el-ph")
                END IF 


                DO iMode = 1 , 3*fi%atoms%ntype
                    DO ispin = 1 , fi%input%jspins
                        ph_linewidth(iMode) = ph_linewidth(iMode) + tpi_const * linewidth(iMode,ispin)
                    END DO 
                END DO 



        END SELECT

        IF (fmpi%irank == 0 ) THEN 
            IF ( iQ .EQ. fi%juPhon%startq ) THEN 
                open( 110, file="linewidth", status='replace', action='write', form='formatted')
                write(110,*) "q-Point", qpts%bk(:,iQ) 
                write(*,*) '-------------------------'
                write(*,*) "Linewidth q-Point", qpts%bk(:,iQ) 
                write(110,*) ph_linewidth(:) 
                write(*,*) ph_linewidth(:) 
                write(*,*) '-------------------------'
            
            ELSE
                open( 110, file="linewidth", status="old", position="append", action="write")
                write(110,*) "q-Point", qpts%bk(:,iQ) 
                write(*,*) '-------------------------'
                write(*,*) "Linewidth q-Point", qpts%bk(:,iQ) 
                write(110,*) ph_linewidth(:) 
                write(*,*) ph_linewidth(:) 
                write(*,*) '-------------------------'
            END IF 
        END IF 

                    
        close(110)



    END SUBROUTINE dfpt_ph_linewidth




END MODULE m_dfpt_elph_linewidth