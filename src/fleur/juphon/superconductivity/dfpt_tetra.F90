!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_tetra
  

#ifdef CPP_MPI
    USE mpi
#endif
    USE m_juDFT
    USE m_types
    USE m_constants
    IMPLICIT NONE 

CONTAINS 
    SUBROUTINE dfpt_tetra_int(fi,results,resultsq,results1,gmat,nbasfcnq_min,k_Int)
        !
        ! subroutine that calculates k integration of type
        ! \sum_k F(k) \delta(\Omega - \varepsilon_{k,\nu'}) \delta( E - \varepsilon_{k,\nu}) 
        ! for Bulk systems for film we need triangular method
        ! Linear Interpolation of the matrix element at the middle of the intersection line
        ! Method implemented as  "P.B. Allen, phys. stat. sol. (b) 120,629 (1983)" 
        ! 
        ! Here no NOCO spin logic is implemented
        USE m_types
        USE m_types_kpts
        USE m_types_juPhon
        USE m_constants
        USE m_npy
        TYPE(t_fleurinput), INTENT(IN) :: fi
        !TYPE(t_mpi),INTENT(IN)         :: fmpi
        TYPE(t_results), INTENT(IN)    :: results
        TYPE(t_results), INTENT(IN)    :: resultsq
        TYPE(t_results), INTENT(IN)    :: results1
        COMPLEX, ALLOCATABLE,INTENT(INOUT) :: gmat(:,:,:,:,:)  !(nu',nu,kpoints,jsp,normal_mode)\
        INTEGER, INTENT(IN) :: nbasfcnq_min
        COMPLEX, ALLOCATABLE, INTENT(INOUT) :: k_int(:,:,:,:)
        
        INTEGER :: noccbd_max, iMode ,icase, ncorners, ispin ,itet , nu , i , icorn , j , jcorn , iNupr
        REAL, ALLOCATABLE :: eigk(:,:,:),eigq(:,:,:),eigk_nondeg(:,:,:),eigq_nondeg(:,:,:),eigkVal(:),eigqVal(:)
        COMPLEX :: area 
        COMPLEX,ALLOCATABLE :: tmp_gmat(:)
        REAL :: efermi 


        ALLOCATE(eigk(size(gmat,2),fi%kpts%nkpt,fi%input%jspins)) ! Not sure about noco logic  
        ALLOCATE(eigq(nbasfcnq_min,fi%kpts%nkpt,fi%input%jspins))
        ALLOCATE(eigk_nondeg,mold=eigk)
        ALLOCATE(eigq_nondeg,mold=eigq)

        eigk(:,:,:) = results%eig(:size(gmat,2),:,:)
        eigq(:,:,:) = resultsq%eig(:nbasfcnq_min,:,:)
        

        ALLOCATE(k_int(nbasfcnq_min,size(gmat,2),size(gmat,4),size(gmat,5) ) )
        k_int = 0 

        ! Consider renormalization of fermi energy
        efermi = results%ef + results1%ef


        ! Film has tetra 3 corners
        ! Bulk tetra has 4 corners (if layered system 2 corners will be degenerate)
        ncorners=SIZE(fi%kpts%ntetra,1)
        ALLOCATE(eigkVal(ncorners))
        ALLOCATE(eigqVal(ncorners))
        ALLOCATE(tmp_gmat(ncorners))

    
        CALL save_npy("gmat_in_call.npy",gmat)
        eigk_nondeg=0
        eigq_nondeg=0

        !! Initial thoughts:
        !! We have to adjust fi%kpts%voltet as we store a normalized volume in t_kpts
        !! voltet_norm = voltet * fi%kpts%ntet / volbz 
        !! We ned voltet not voltet_norm 
        !! However, going from sum to integral we need 
        !! sum_k --> V/(tpi_const)**3 int 
        !! fi%kpts%voltet accounts for this factor, comment remains for future insight
        !ALLOCATE(voltetra,mold=fi%kpts%voltet)
        !volbz = tpi_const**3/fi%cell%omtil ! 
        !
        !voltetra = fi%kpts%voltet * volbz / fi%kpts%ntet


        eigk_nondeg(:,:,:)=eigk(:,:,:)
        eigq_nondeg(:,:,:)=eigq(:,:,:)

        
        !
        ! This part is from /dos/tetra_dos.F90
        ! care for degeneracy shift the edges 
        ! If no tetra has degeneracy then k+q also contains no dengenerate corners
        
        
        CALL timestart("Tetrahedon Degeneracy Test k")
        ! for the states nu (occupied states (mostly))
        ! we need some unoccupied states at some corners 
        ! e.g. e1, e2 , e3 occupied and e4 is unoccupied
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO itet = 1 , fi%kpts%ntet
                DO nu = 1, size(eigk,1) 
                    DO i=1, ncorners !corners
                        icorn = fi%kpts%ntetra(i,itet)
                        DO j = i+1,ncorners !corner
                            jcorn = fi%kpts%ntetra(j,itet)
                            IF (abs(eigk_nondeg(nu,icorn,ispin)-eigk_nondeg(nu,jcorn,ispin)).LT.fi%juPhon%eDiffcut) THEN 
                                eigk_nondeg(nu,icorn,ispin) = eigk_nondeg(nu,icorn,ispin) + i*fi%juPhon%eDiffcut*itet ! maybe just rewrite this as fi%juPhon%edifCut only 
                                eigk_nondeg(nu,jcorn,ispin) = eigk_nondeg(nu,jcorn,ispin) - i*fi%juPhon%eDiffcut*itet  
                            END IF     
                        END DO !j
                    END DO !i
                END DO !nu 
            END DO !itet 
        END DO !ispin 
        CALL timestop("Tetrahedon Degeneracy Test k")

        
        ! for the states iNupr (occupied and unoccupied nu')
        CALL timestart("Tetrahedon Degeneracy Test k+q")
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO itet = 1 , fi%kpts%ntet
                DO iNupr = 1, nbasfcnq_min
                    DO i=1,ncorners !corners
                        icorn = fi%kpts%ntetra(i,itet)
                        DO j = i+1,ncorners !corner
                            jcorn = fi%kpts%ntetra(j,itet)
                            IF (abs(eigq_nondeg(iNupr,icorn,ispin)-eigq_nondeg(iNupr,jcorn,ispin)).LT.fi%juPhon%eDiffcut) THEN
                                eigq_nondeg(iNupr,icorn,ispin) = eigq_nondeg(iNupr,icorn,ispin) + i*fi%juPhon%eDiffcut*itet
                                eigq_nondeg(iNupr,jcorn,ispin) = eigq_nondeg(iNupr,jcorn,ispin) - i*fi%juPhon%eDiffcut*itet   
                            END IF     
                        END DO !j
                    END DO !i
                END DO !iNpur
            END DO !itet 
        END DO !ispin 
        CALL timestop("Tetrahedon Degeneracy Test k+q")

        CALL timestart("Area of Intersection")
        DO iMode = 1 , size(gmat,5)
            DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
                DO itet = 1 , fi%kpts%ntet
                    DO nu = 1 ,size(gmat,2)
                        DO iNupr = 1,nbasfcnq_min
                            DO i=1,ncorners !corners
                                icorn = fi%kpts%ntetra(i,itet)
                                eigkVal(i)  = eigk_nondeg(nu,icorn,ispin)
                                eigqVal(i) = eigq_nondeg(iNupr,icorn,ispin)
                                tmp_gmat(i) = gmat(iNupr,nu,icorn,ispin,iMode) ! we give the nu' nu element for the k points at the tetra corners
                            END DO !corners
                            call tetra_area(eigkVal,eigqVal,results%ef,fi%kpts%voltet(itet),tmp_gmat,area,icase) !results%ef voltetra(itet)
                            k_int(iNupr,nu,ispin,iMode) = k_int(iNupr,nu,ispin,iMode) +  area 
                        END DO !iNupr
                    END DO !nu
                END DO !itet
            END DO !iSpin
        END DO !iMode
        CALL timestop("Area of Intersection")
        

        !This is needed if not fi%kpts%voltet is taken
        !k_int(:,:,:,:) = fi%cell%omtil/(tpi_const)**3 * k_int(:,:,:,:)
        !Consider the weight of one tetra
        k_int(:,:,:,:) = k_int(:,:,:,:) / fi%kpts%ntet
        
    END SUBROUTINE dfpt_tetra_int

    SUBROUTINE tetra_area(eigk,eigq,efermi,voltet,gmat,area,icase)
        ! In notation of the Paper "P.B. Allen, phys. stat. sol. (b) 120,629 (1983)" 
        ! tmp_mat --> p_1,p_2,p_3 etc.
        ! tmp_arr --> a_1,a_2,a_3 etc.
        ! prefac --> a,b,c
        REAL, INTENT(INOUT) :: eigk(:)
        REAL, INTENT(INOUT) :: eigq(:)
        REAL, INTENT(IN)    :: efermi
        REAL, INTENT(IN)    :: voltet
        COMPLEX, INTENT(INOUT) :: gmat(:)
        COMPLEX, INTENT(OUT) :: area
        INTEGER, INTENT(OUT) :: icase
        
        REAL                :: tmp_arr(3)
        REAL                :: prefac
        INTEGER :: i,j, case
        REAL    :: intersection_val , f , intersection_val2 ,f2 
        COMPLEX :: interpol_mat, interpol_mat2, tmp_mat(3)!, area 


        ! area can be complex from matrix element 


        !!!
        !!! sort the energies
        !!! e(1) < e(2) < e(3) < e(4)
        !!! also sort eigq and gmat that pairs still match 
        call sorting(eigk,r_arr2=eigq,c_arr2=gmat)

        
        IF ( eigk(1) .LT. efermi  .AND. efermi .LT. eigk(2) ) THEN 
            !case=1
            prefac= ( efermi - eigq(1) ) / ( efermi - eigk(1) )

            ! temporary arrray for energy quotient
            tmp_arr(1) = (eigq(2) - eigq(1)) / (eigk(2) - eigk(1)) 
            tmp_arr(2) = (eigq(3) - eigq(1)) / (eigk(3) - eigk(1)) 
            tmp_arr(3) = (eigq(4) - eigq(1)) / (eigk(4) - eigk(1)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (gmat(2) -gmat(1)) / (eigk(2) - eigk(1))
            tmp_mat(2)= (gmat(3) - gmat(1)) / (eigk(3) - eigk(1)) 
            tmp_mat(3)= (gmat(4) - gmat(1)) / (eigk(4) - eigk(1)) 



            call surface_intersection(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)

            f =  (efermi - eigk(1))/ ( (eigk(2) - eigk(1)) * (eigk(3) - eigk(1)) * (eigk(4) - eigk(1)) )

            interpol_mat  = gmat(1) + 0.5 * (efermi - eigk(1)) * interpol_mat 


            area = 6*voltet*f*intersection_val * interpol_mat 
            
            icase = 1 

        ELSE IF ( eigk(2) .LT. efermi  .AND. efermi .LT. eigk(3) ) THEN 
            !case=2 
            !!!
            !!! This contribution I_0
            !!!
            prefac= ( efermi - eigq(1) ) / ( efermi - eigk(1) )

            tmp_arr(1) = (eigq(2) - eigq(1)) / (eigk(2) - eigk(1)) 
            tmp_arr(2) = (eigq(3) - eigq(1)) / (eigk(3) - eigk(1)) 
            tmp_arr(3) = (eigq(4) - eigq(1)) / (eigk(4) - eigk(1)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (gmat(2) -gmat(1)) / (eigk(2) - eigk(1))
            tmp_mat(2)= (gmat(3) - gmat(1)) / (eigk(3) - eigk(1)) 
            tmp_mat(3)= (gmat(4) - gmat(1)) / (eigk(4) - eigk(1)) 

            call surface_intersection(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)

            f =  (efermi - eigk(1))/ ( (eigk(2) - eigk(1)) * (eigk(3) - eigk(1)) * (eigk(4) - eigk(1)) )

            interpol_mat  = gmat(1) + 0.5 * (efermi - eigk(1)) * interpol_mat 

        
            !!!
            !!! This is contribution I_1
            !!!
            prefac= ( efermi - eigq(2) ) / ( efermi - eigk(2) )

            tmp_arr(1) = (eigq(1) - eigq(2)) / (eigk(1) - eigk(2)) 
            tmp_arr(2) = (eigq(3) - eigq(2)) / (eigk(3) - eigk(2)) 
            tmp_arr(3) = (eigq(4) - eigq(2)) / (eigk(4) - eigk(2)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (gmat(1) - gmat(2)) / (eigk(1) - eigk(2))
            tmp_mat(2)= (gmat(3) - gmat(2)) / (eigk(3) - eigk(2)) 
            tmp_mat(3)= (gmat(4) - gmat(2)) / (eigk(4) - eigk(2)) 

            call surface_intersection(tmp_arr,prefac,intersection_val2,tmp_mat,interpol_mat2)

            f2= (efermi - eigk(2))/ ( (eigk(2) - eigk(1)) * (eigk(3) - eigk(2)) * (eigk(4) - eigk(2)) )

            interpol_mat2  = gmat(2) + 0.5 * (efermi - eigk(2)) * interpol_mat2 

            area = 6 * voltet * (f*intersection_val * interpol_mat - f2*intersection_val2 * interpol_mat2)

            icase = 2
        ELSE IF ( eigk(3) .LT. efermi  .AND. efermi .LT. eigk(4) ) THEN 
            !case=3 
            prefac= ( efermi - eigq(4) ) / ( efermi - eigk(4) )

            tmp_arr(1) = (eigq(1) - eigq(4)) / (eigk(1) - eigk(4)) 
            tmp_arr(2) = (eigq(2) - eigq(4)) / (eigk(2) - eigk(4)) 
            tmp_arr(3) = (eigq(3) - eigq(4)) / (eigk(3) - eigk(4)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (gmat(1) - gmat(4)) / (eigk(1) - eigk(4))
            tmp_mat(2)= (gmat(2) - gmat(4)) / (eigk(2) - eigk(4)) 
            tmp_mat(3)= (gmat(3) - gmat(4)) / (eigk(3) - eigk(4)) 

            call surface_intersection(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)

            f= (eigk(4)  - efermi )/ ( (eigk(4) - eigk(1)) * (eigk(4) - eigk(2)) * (eigk(4) - eigk(3)) )

            interpol_mat = gmat(4) + 0.5* (efermi - eigk(4)) *interpol_mat

            area = 6 * voltet*f*intersection_val * interpol_mat

            icase = 3 
        ELSE 
            area = 0 
        END IF

    END SUBROUTINE tetra_area

    SUBROUTINE surface_intersection(arr_sort,prefac,intersection_val, arr_mat,interpol_mat) 
        
        !!! Calculates the intersection line between the k+q and k hypersurface
        !!! Interpolates the matrix element at this intersection line 
        REAL,INTENT(INOUT) :: arr_sort(3)
        REAL,INTENT(IN)    :: prefac
        REAL,INTENT(INOUT)   :: intersection_val
        COMPLEX, INTENT(INOUT):: arr_mat(3)
        COMPLEX, INTENT(OUT)  :: interpol_mat

        REAL    :: tmp 
        INTEGER :: i, j 

        CALL sorting(arr_sort,c_arr2=arr_mat)

        IF ( arr_sort(1) .LT. prefac .AND. prefac .LT. arr_sort(2) ) THEN
            intersection_val = (prefac - arr_sort(1)) / ((arr_sort(2) - arr_sort(1)) * (arr_sort(3) - arr_sort(1)))
            
            interpol_mat = 2* arr_mat(1) + (prefac - arr_sort(1))/(arr_sort(2)-arr_sort(1)) * ( arr_mat(2) - arr_mat(1)) &
                            +(prefac - arr_sort(1))/(arr_sort(3)-arr_sort(1)) * (arr_mat(3)-arr_mat(1)) 
        ELSE IF ( arr_sort(2) .LT. prefac .AND. prefac .LT. arr_sort(3) ) THEN
            intersection_val = (arr_sort(3) - prefac ) / ((arr_sort(3) - arr_sort(2)) * (arr_sort(3) - arr_sort(1)))

            interpol_mat = 2* arr_mat(3) + (arr_sort(3) - prefac )/( arr_sort(3)-arr_sort(2)) * ( arr_mat(2) - arr_mat(3)) &
                            +(arr_sort(3) - prefac)/(arr_sort(3)-arr_sort(1)) * (arr_mat(1)-arr_mat(3)) 
        ELSE
            intersection_val = 0
            interpol_mat = 0
            
        END IF 

    END SUBROUTINE surface_intersection

    SUBROUTINE sorting(arr1,r_arr2,c_arr2)

        !
        ! This subroutine sorts arr form least to highest value
        ! If arr2 is present, then arr2 is sorted that the pairs between arr1 and arr2 remain
        ! after arr1 is sorrted

        REAL, INTENT(INOUT) :: arr1(:)
        REAL, OPTIONAL, INTENT(INOUT) :: r_arr2(:)
        COMPLEX, OPTIONAL, INTENT(INOUT) :: c_arr2(:)

        REAL :: tmp
        COMPLEX :: c_tmp
        INTEGER :: i,j

        DO i = 1,size(arr1)
            DO j = i+1,size(arr1)
                IF (arr1(i) .GT. arr1(j)) THEN
                    tmp = arr1(j)
                    arr1(j) = arr1(i)
                    arr1(i) = tmp 
                
                    IF (PRESENT(r_arr2)) THEN
                        tmp = r_arr2(j)
                        r_arr2(j) = r_arr2(i)
                        r_arr2(i) = tmp 
                    END If

                    IF (PRESENT(c_arr2)) THEN
                        c_tmp = c_arr2(j)
                        c_arr2(j) = c_arr2(i)
                        c_arr2(i) = c_tmp 
                    END IF 
                END IF
            END DO !j  
        END DO !i
    
    END SUBROUTINE sorting









END MODULE m_dfpt_tetra