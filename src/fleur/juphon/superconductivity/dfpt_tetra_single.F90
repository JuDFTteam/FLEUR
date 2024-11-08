!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_tetra_single
  

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT
    USE m_types
    USE m_constants
    IMPLICIT NONE 
    
    CONTAINS 

    SUBROUTINE dfpt_tetra_single(fi,results,resultsq,results1,gmat,shift,nuWindow)!,kInt)

        USE m_npy
        !
        ! This routine calculates /int F(k,k+q) \delta(\epsilon_k - \epsilon_{k+q} - \omega)
        ! At the corners of the tetrahedons we put \tilde{\epsilon} = \epsilon_k - \epsilon_{k+q}
        ! Thereby reducing the problem to the hypersurface \omega crossing the hypersurface of \epsilon_k - \epsilon_{k+q}
        ! Implementation done as in: "Journal of Physics C: Solid State Physics 12.15 (1979): 2991"

        TYPE(t_fleurinput), INTENT(IN) :: fi 
        TYPE(t_results),    INTENT(IN) :: results 
        TYPE(t_results),    INTENT(IN) :: resultsq 
        TYPE(t_results),    INTENT(IN) :: results1 
        REAL, ALLOCATABLE,  INTENT(IN) :: gmat(:,:,:,:,:) ! iNupr,nu,kpts,sigma , iMode
        REAL,               INTENT(IN) :: shift(:) ! eigenvalues dynmat
        INTEGER,            INTENT(IN) :: nuWindow(2,2)
        !REAL, ALLOCATABLE,  INTENT(OUT):: kInt

        REAL, ALLOCATABLE :: eig_nondeg(:,:,:) 
        REAL :: eig(4), tmpMat(4), valArea
        INTEGER :: ind, indPr, comInd, ispin, itet, nu, iNupr, i, j, icorn, jcorn, ncorners , iMode


        ALLOCATE(eig_nondeg(size(gmat,1)*size(gmat,2),fi%kpts%nkpt,fi%input%jspins))

        eig_nondeg = 10000000 


        ! Film has tetra 3 corners
        ! Bulk tetra has 4 corners (if layered system 2 corners will be degenerate)
        ncorners=SIZE(fi%kpts%ntetra,1)

        CALL timestart("Tetrahedon Degeneracy Test k")
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO itet = 1 , fi%kpts%ntet
                DO nu = nuWindow(1,1), nuWindow(1,2) 
                    ind =  nu - nuWindow(1,1) + 1 
                    DO iNupr = nuWindow(2,1), nuWindow(2,2)
                        indPr = iNupr- nuWindow(2,1) + 1
                        comInd = (ind-1) * (nuWindow(2,2) - nuWindow(2,1)+1) + (indPr)  
                        write(2100,*)  "commInd" , comInd , "indPr" , indPr , "ind" , ind , "indMax" , nuWindow(2,2) - nuWindow(2,2)+1
                        DO i=1, ncorners !corners
                            icorn = fi%kpts%ntetra(i,itet)
                            eig_nondeg(comInd,icorn,ispin) = results%eig(nu,icorn,ispin) - resultsq%eig(iNupr,icorn,ispin) 
                            DO j = i+1,ncorners !corner
                                jcorn = fi%kpts%ntetra(j,itet)
                                eig_nondeg(comInd,jcorn,ispin) = results%eig(nu,jcorn,ispin) - resultsq%eig(iNupr,jcorn,ispin) 
                                IF (abs(eig_nondeg(comInd,icorn,ispin)-eig_nondeg(comInd,jcorn,ispin)).LT.fi%juPhon%eDiffcut) THEN 
                                    eig_nondeg(comInd,icorn,ispin) = eig_nondeg(comInd,icorn,ispin) + i*fi%juPhon%eDiffcut*itet 
                                    eig_nondeg(comInd,jcorn,ispin) = eig_nondeg(comInd,jcorn,ispin) - i*fi%juPhon%eDiffcut*itet  
                                END IF     
                            END DO !j
                        END DO !i
                    END DO !iNupr
                END DO !nu 
            END DO !itet 
        END DO !ispin 
        call save_npy("eig_nondeg.npy", eig_nondeg)
        CALL timestop("Tetrahedon Degeneracy Test k")

        CALL timestart("Tetra int")
        DO iMode = 1 , size(gmat,5)
            IF (shift(iMode) .LT. 0) THEN
                WRITE(*,*) '-------------------------'
                WRITE(*,*) 'linewidth: Eigenvalue imaginary --> Phonon linewidth set to zero'
                WRITE(*,*) '-------------------------' 
                CYCLE 
            END IF 
            DO ispin = 1, MERGE(1,fi%input%jspins,fi%noco%l_noco)
                DO itet = 1 , fi%kpts%ntet
                    !DO nu = nuWindow(1,1), nuWindow(1,2) 
                    DO ind = 1, size(gmat,2)
                        !ind =  nu - nuWindow(1,1) + 1 
                        !DO iNupr = nuWindow(2,1), nuWindow(2,2)
                        DO indPr = 1, size(gmat,1)
                            !indPr = iNupr- nuWindow(2,1) + 1
                            comInd = (ind-1) * (nuWindow(2,2) - nuWindow(2,1)+1) + (indPr)  
                            write(2101,*)  "commInd" , comInd , "indPr" , indPr , "ind" , ind , "indMax" , nuWindow(2,2) - nuWindow(2,2)+1
                            DO i = 1 , ncorners
                                icorn = fi%kpts%ntetra(i,itet)
                                eig(i) = eig_nondeg(comInd,icorn,ispin)
                                tmpMat(i) = gmat(indPr,ind,icorn,ispin,iMode)  
                            END DO !corners
                            CALL tetra_area(eig,tmpMat,shift(iMode),valArea)
                        END DO !iNuPr
                    END DO !nu
                END DO  !itet
            END DO !ispin 
        END DO !iMode 
        CALL timestop("Tetra int")

    END SUBROUTINE dfpt_tetra_single


    SUBROUTINE tetra_area(eig,tmpMat,shift,valArea) 

        REAL, INTENT(INOUT) :: eig(4)
        REAL, INTENT(INOUT) :: tmpMat(4)
        REAL, INTENT(IN)    :: shift
        REAL, INTENT(OUT)   :: valArea 


        REAL :: interSect(4), tmpVec(3)
        REAL :: prefactor,scaling 

        INTEGER :: ind 

        valArea = 0.0

        CALL sorting(eig,tmpMat)

        ! Case 1 

        IF ( (eig(1) .LT. shift) .AND. (eig(2) .GE. shift)  ) THEN 
            ! tmpVec contains f_{2,1} , f_{3,1} , f_{4,1} 
            tmpVec(1) = (shift - eig(1)) / (eig(2)- eig(1))
            tmpVec(1) = (shift - eig(1)) / (eig(3)- eig(1))
            tmpVec(1) = (shift - eig(1)) / (eig(3)- eig(1))

            scaling = tmpVec(1)*tmpVec(2)*tmpVec(3)
            
            prefactor = 3* scaling / (  shift - eig(1)  )

            interSect(1) = 1/3 * (  (1-tmpVec(1)) + (1-tmpVec(2)) + (1-tmpVec(3)) )
            interSect(2) = 1/3 * tmpVec(2)
            interSect(3) = 1/3 * tmpVec(3)
            interSect(4) = 1/3 * tmpVec(4)


            DO ind = 1 , 4 
                valArea = valArea + scaling* interSect(ind)*tmpMat(ind) 
            END DO 

        ELSE IF ( (eig(2) .LT. shift) .AND. (eig(3) .GE. shift)  ) THEN 


        ELSE IF ( (eig(3) .LT. shift) .AND. (eig(4) .GE. shift)  ) THEN 
                      

        ELSE 
            valArea = 0.0 
        END IF 

    END SUBROUTINE tetra_area 

    SUBROUTINE sorting(arr1,arr2)


        REAL, INTENT(INOUT) :: arr1(4)
        REAL, INTENT(INOUT) :: arr2(4)

        REAL :: tmp
        INTEGER :: i,j

        DO i = 1,size(arr1)
            DO j = i+1,size(arr1)
                IF (arr1(i) .GT. arr1(j)) THEN
                    tmp = arr1(j)
                    arr1(j) = arr1(i)
                    arr1(i) = tmp 
                
                    tmp = arr2(j)
                    arr2(j) = arr2(i)
                    arr2(i) = tmp 
                END IF
            END DO !j  
        END DO !i
    
    END SUBROUTINE sorting

    

END MODULE m_dfpt_tetra_single