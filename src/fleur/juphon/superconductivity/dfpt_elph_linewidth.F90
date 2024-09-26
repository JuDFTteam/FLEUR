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
    SUBROUTINE dfpt_ph_linewidth(fi,qpts,results,resultsq,results1,eigenVals,gmat,iQ,nbasfcnq_min,ph_linewidth)
        ! This subroutine calculates the phonon linewdith
        ! Currently implemented is the linewidth calcualtion with smearing 

        USE m_dosbin
        USE m_smooth
        USE m_dfpt_fermie, ONLY : sfermi
        USE m_dfpt_tetra
        USE m_intgr, ONLY : intgz0

        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_kpts), INTENT(IN) :: qpts
        TYPE(t_results), INTENT(IN) :: results,resultsq, results1
        REAL,ALLOCATABLE,INTENT(IN) :: eigenVals(:)
        COMPLEX,ALLOCATABLE,INTENT(INOUT) :: gmat(:,:,:,:,:)
        INTEGER, INTENT(IN) :: iQ
        INTEGER, INTENT(IN) :: nbasfcnq_min
        REAL, ALLOCATABLE, INTENT(OUT)    :: ph_linewidth(:)

        INTEGER :: iNupr,nu, ispin , gridPoint , nk , nZero , iMode
        REAL :: emin, emax , x ,xq
        REAL, ALLOCATABLE :: eGrid(:),linewidth(:,:)
        COMPLEX,ALLOCATABLE:: kInt_gmat(:,:,:,:) !(nu',nu,jsp,normal_mode)

        REAL , ALLOCATABLE :: gauss(:) 
        REAL :: smearing  , intOut



        smearing = 1e-7
        gmat(:,:,:,:,:) = ABS(gmat(:,:,:,:,:))**2 


        SELECT CASE(fi%input%bz_integration)

            CASE(BZINT_METHOD_HIST,BZINT_METHOD_GAUSS)
                

                DO ispin = 1 , fi%input%jspins
                    DO nk = 1 , fi%kpts%nkpt
                        DO nu = 1 , size(gmat,2)
                            x = (results%eig(nu,nk,ispin)-results%ef)/fi%input%tkb
                            DO iNupr = 1 , size(gmat,1)
                                xq = (resultsq%eig(iNupr,nk,ispin)-results%ef)/fi%input%tkb
                                gmat(iNupr,nu,nk,ispin,:) = gmat(iNupr,nu,nk,ispin,:)*(sfermi(x) - sfermi(xq))
                            END DO ! iNupr 
                        END DO  ! nu 
                    END DO ! nk 
                END DO 

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
                    gauss(gridPoint) = 1/SQRT(tpi_const*smearing ) *EXP( -eGrid(gridPoint)**2/(2*smearing))
                END DO

                DO iMode = 1 , 3*fi%atoms%ntype
                    IF (eigenVals(iMode) .GE. 0.0 ) THEN 
                        ! If omega becomes negative the deltra distribution is never satisfied as IM(eig,eigq) = 0.0 
                        CALL dos_bin_transport(fi%input%jspins,fi%kpts%wtkpt,eGrid,results%eig(:size(gmat,2),:,:) , resultsq%eig(:nbasfcnq_min,:,:), REAL(gmat(:,:,:,:,iMode)), linewidth, -SQRT(eigenVals(iMode)),spinDeg=.FALSE.)
                        !CALL smooth(eGrid,linewidth(:,1),fi%banddos%sig_dos,size(eGrid))
                        DO ispin = 1 , fi%input%jspins
                            CALL intgz0(gauss*linewidth(:,ispin), eGrid(2)-eGrid(1) , size(gauss) , intOut,.FALSE.)
                            !ph_linewidth(iMode) =  tpi_const * SQRT(eigenVals(iMode))/fi%kpts%nkpt*REAL(linewidth(nZero,1))
                            ph_linewidth(iMode) =  ph_linewidth(iMode) +  tpi_const/fi%input%jspins * SQRT(eigenVals(iMode))/fi%kpts%nkpt*intOut
                        END DO 
                    ELSE
                        CALL juDFT_warn('linewidth: Eigenvalue imaginary --> Phonon linewidth set to zero')
                    END IF 
                END DO 
            
            CASE(BZINT_METHOD_TETRA)

                CALL timestart("k-Integration el-ph")
                CALL dfpt_tetra_int(fi,results,resultsq, results1, gmat,nbasfcnq_min,kInt_gmat)
                CALL timestop("k-Integration el-ph")

                DEALLOCATE(gmat)
                ALLOCATE(ph_linewidth( 3*fi%atoms%ntype))
                
                ph_linewidth = 0.0 
                DO iMode = 1 , 3*fi%atoms%ntype
                    DO ispin = 1 , fi%input%jspins
                        DO nu = 1 , size(kInt_gmat,2)
                            DO iNupr = 1 , nbasfcnq_min
                                ph_linewidth(iMode) =  ph_linewidth(iMode) + tpi_const/fi%input%jspins * SQRT(eigenVals(iMode))/fi%kpts%nkpt*REAL(kInt_gmat(iNupr,nu,ispin,iMode))
                            END DO 
                        END DO 
                    END DO 
                END DO 
                


        END SELECT

        IF ( iQ .EQ. fi%juPhon%startq ) THEN 
            open( 110, file="linewidth", status='replace', action='write', form='formatted')
            write(110,*) "q-Point", qpts%bk(:,iQ) 
            write(*,*) "Linewidth q-Point", qpts%bk(:,iQ) 
            write(110,*) ph_linewidth(:) 
            write(*,*) ph_linewidth(:) 
        
        ELSE
            open( 110, file="linewidth", status="old", position="append", action="write")
            write(110,*) "q-Point", qpts%bk(:,iQ) 
            write(*,*) "Linewidth q-Point", qpts%bk(:,iQ) 
            write(110,*) ph_linewidth(:) 
            write(*,*) ph_linewidth(:) 
        END IF 

                    
        close(110)



    END SUBROUTINE dfpt_ph_linewidth




END MODULE m_dfpt_elph_linewidth