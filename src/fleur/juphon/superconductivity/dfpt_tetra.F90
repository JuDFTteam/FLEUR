!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


MODULE m_dfpt_tetra
    IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_k_int(fi,fmpi , results, eig, eigq , noccbd, nbasfcnq, el_ph_mat )
        !
        ! subroutine that calculates k integration of type
        ! \sum_k F(k) \delta(\Omega - \varepsilon_{k,\nu'}) \delta( E - \varepsilon_{k,\nu}) 
        ! for Bulk systems for film we need triangular method
        ! Linear Interpolation of the matrix element at the middle of the intersection line
        ! Method implemented as  "P.B. Allen, phys. stat. sol. (b) 120,629 (1983)" 
        ! 
        USE m_types
        USE m_types_kpts
        USE m_types_juPhon
        USE m_constants
        USE m_npy
        TYPE(t_fleurinput), INTENT(IN) :: fi
        TYPE(t_mpi),INTENT(IN)         :: fmpi
        TYPE(t_results), INTENT(IN)    :: results
        REAL, INTENT(IN) :: eig(:,:,:), eigq(:,:,:) 
        INTEGER, INTENT(IN) :: noccbd(:,:), nbasfcnq(:,:)
        COMPLEX, INTENT(IN) :: el_ph_mat(:,:,:,:)
        REAL, ALLOCATABLE :: eig_nondeg(:,:,:) , eig_nondegq(:,:,:)
        !REAL, ALLOCATABLE :: eig_q_sort(:,:,:,:) !nu,iNupr,kpt, spin
        TYPE(t_kpts)  :: kpts_local
        INTEGER :: ispin,itet,nu,iNupr,i,j,icorn,jcorn
        REAL    :: tmp_e,tmp_eq
        LOGICAL :: l_timeReversalCheck
        INTEGER :: ncorners

        REAL, ALLOCATABLE :: eigval(:), eigvalq(:),vol_tetra(:)
        COMPLEX, ALLOCATABLE :: tmp_el_ph(:)

        integer :: swapped

        REAL :: volbz
        COMPLEX :: k_int !kint k_resolved with gmat can be complex 
        COMPLEX, ALLOCATABLE :: k_resolved(:,:,:)


        !ALLOCATE(eig_nondeg,mold=eig)   !what is neigd 
        !ALLOCATE(eig_nondegq,mold=eigq)
        
        ALLOCATE(eig_nondeg(MINVAL(noccbd),size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco))) 
        ALLOCATE(eig_nondegq(MINVAL(nbasfcnq),size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco))) 
        
        ALLOCATE(k_resolved(MINVAL(noccbd),MINVAL(nbasfcnq),MERGE(1,fi%input%jspins,fi%noco%l_noco)))

        !ALLOCATE(eig_q_sort(MINVAL(noccbd),MINVAL(nbasfcnq),size(fmpi%k_list),MERGE(1,fi%input%jspins,fi%noco%l_noco)))
        
        print *, "Fermi in Tetra"
        print *, results%ef

        !l_timeReversalCheck = .FALSE.
        !IF(.NOT.fi%banddos%band.AND..NOT.fi%banddos%dos) THEN
        !    IF(fi%noco%l_soc.OR.fi%noco%l_ss) l_timeReversalCheck = .TRUE.
        !END IF

        !CALL kpts_local%init(fi%sym, fi%input%film, fi%hybinp%l_hybrid .or. fi%input%l_rdmft, l_timeReversalCheck)
        !CALL kpts_local%initTetra(fi%input, fi%cell, fi%sym, fi%noco%l_soc .OR. fi%noco%l_ss)
  
        !print *, "I should have initialized the new k point type"

        !DO itet=1 , kpts_local%ntet
        !   DO i=1,4 !corners
        !      icorn = kpts_local%ntetra(i,itet)
        !      write(2221,*),  icorn , kpts_local%bk(:,icorn), "local kpoints"
        !   END DO 
        !END DO 

        ! Film has tetra 3 corners
        ! Bulk tetra has 4 corners (if layered system 2 corners will be degenerate)
        ncorners=SIZE(fi%kpts%ntetra,1)
        ALLOCATE(eigval(ncorners))
        ALLOCATE(eigvalq(ncorners))
        ALLOCATE(tmp_el_ph(ncorners))

        eig_nondeg=0
        eig_nondegq=0


        !! We have to adjust fi%kpts%voltet as we store a normalized voluem in t_kpts
        !! voltet_norm = voltet * fi%kpts%ntet / volbz 
        !! We ned voltet not voltet_norm 
        ALLOCATE(vol_tetra,mold=fi%kpts%voltet)
        volbz = tpi_const**3/fi%cell%vol
        
        vol_tetra = fi%kpts%voltet * volbz / fi%kpts%ntet


        eig_nondeg(:,:,:)=eig(:MINVAL(noccbd),:,:)
        eig_nondegq(:,:,:)=eigq(:MINVAL(nbasfcnq),:,:)

        !print *, ncorners ,"Size "
        call save_npy("tetra_in_routine.npy",fi%kpts%ntetra)

        DO itet=1 , fi%kpts%ntet
            write(0002,*), fi%kpts%voltet(itet)
           DO i=1,ncorners !corners
              icorn = fi%kpts%ntetra(i,itet)
              write(2222,*), icorn,  fi%kpts%bk(:,icorn), "fi kpts"
           END DO 
        END DO 

        !STOP

        !
        ! This part is from /dos/tetra_dos.F90
        ! care for degeneracy shift the edges 
        ! If no tetra has degeneracy then k+q also contains no dengenerate corners
        
        
        CALL timestart("Tetrahedon Degenercy Test k")
        ! for the states nu (occupied states)
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO itet = 1 , fi%kpts%ntet
                DO nu = 1, MINVAL(noccbd) 
                    DO i=1, ncorners !corners
                        icorn = fi%kpts%ntetra(i,itet)
                        DO j = i+1,ncorners !corner
                            jcorn = fi%kpts%ntetra(j,itet)
                            IF (abs(eig_nondeg(nu,icorn,ispin)-eig_nondeg(nu,jcorn,ispin)).LT.fi%juPhon%eDiffcut) THEN 
                                eig_nondeg(nu,icorn,ispin) = eig_nondeg(nu,icorn,ispin) + i*fi%juPhon%eDiffcut*itet
                                eig_nondeg(nu,jcorn,ispin) = eig_nondeg(nu,jcorn,ispin) - i*fi%juPhon%eDiffcut*itet  
                            END IF     
                        END DO !j
                    END DO !i
                END DO !nu 
            END DO !itet 
        END DO !ispin 
        call save_npy("eigen_nondeg.npy",eig_nondeg)
        CALL timestop("Tetrahedon Degenercy Test k")

        
        ! for the states iNupr (occupied and unoccupied nu')
        CALL timestart("Tetrahedon Degenercy Test k+q")
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO itet = 1 , fi%kpts%ntet
                DO iNupr = 1,  MINVAL(nbasfcnq)
                    DO i=1,ncorners !corners
                        icorn = fi%kpts%ntetra(i,itet)
                        DO j = i+1,ncorners !corner
                            jcorn = fi%kpts%ntetra(j,itet)
                            IF (abs(eig_nondegq(iNupr,icorn,ispin)-eig_nondegq(iNupr,jcorn,ispin)).LT.fi%juPhon%eDiffcut) THEN
                                eig_nondegq(iNupr,icorn,ispin) = eig_nondegq(iNupr,icorn,ispin) + i*fi%juPhon%eDiffcut*itet
                                eig_nondegq(iNupr,jcorn,ispin) = eig_nondegq(iNupr,jcorn,ispin) - i*fi%juPhon%eDiffcut*itet   
                            END IF     
                        END DO !j
                    END DO !i
                END DO !iNpur
            END DO !itet 
        END DO !ispin 
        call save_npy("eigen_nondegq.npy",eig_nondegq)
        CALL timestop("Tetrahedon Degenercy Test k+q")
    


        CALL timestart("Area of Intersection")
        DO ispin = 1 , MERGE(1,fi%input%jspins,fi%noco%l_noco)
            DO nu = 1 , 1!MINVAL(noccbd)  ! TODO think the do loop has to be adjusted to accurately sum the right states
                DO iNupr=65, MINVAL(nbasfcnq)
                    k_int=0.0
                    DO itet = 1 , fi%kpts%ntet
                        DO i=1,ncorners !corners
                            icorn = fi%kpts%ntetra(i,itet)
                            eigval(i)  = eig_nondeg(nu,icorn,ispin)
                            eigvalq(i) = eig_nondegq(iNupr,icorn,ispin)
                            tmp_el_ph(i) = el_ph_mat(nu,iNupr,icorn,ispin) ! we give the nu' nu element for the k points at the tetra corners
                        END DO !i
                        call tetra_area(eigval,eigvalq,0.00,vol_tetra(itet),tmp_el_ph,k_int) !results%ef
                        if (k_int .NE. 0) THEN 
                            write(3500,*), "Outer stuff"
                            write(3500,*), k_int, itet , nu , iNupr
                            write(3500,*), eigval(1) , eigval(2) , eigval(3) ,eigval(4)
                            write(3500,*), eigvalq(1) , eigvalq(2) , eigvalq(3) ,eigvalq(4)
                        end if 
                    END DO !itet 
                    k_resolved(nu,iNupr,ispin) = k_int
                END DO !iNupr 
            END DO !nu
        END DO !ispin 
        CALL timestop("Area of Intersection")
        
        call save_npy("k_resolved.npy",k_resolved)
        !call save_npy("eig_sorted.npy",eig_nondeg)
        !call save_npy("eigq_sorted.npy",eig_q_sort)

        
    END SUBROUTINE dfpt_k_int



    SUBROUTINE tetra_area(eig_k, eig_q ,e_fermi,voltet, el_ph_mat ,k_int)

        REAL, INTENT(INOUT) :: eig_k(:)
        REAL, INTENT(INOUT) :: eig_q(:)
        REAL, INTENT(IN)    :: e_fermi
        REAL, INTENT(IN)    :: voltet
        COMPLEX, INTENT(INOUT) :: el_ph_mat(:)
        COMPLEX, INTENT(INOUT) :: k_int

        
        REAL                :: tmp_arr(3)
        REAL                :: prefac
        INTEGER :: i,j, case
        REAL    :: tmp_e, tmp_eq, intersection_val , f , intersection_val2 ,f2 
        COMPLEX :: interpol_mat, interpol_mat2, tmp_mat(3), area 

        ! area can be complex from matrix element 
    
        !!!
        !!! sort the energies
        !!! e(1) < e(2) < e(3) < e(4)
        call sorting(eig_k,r_arr2=eig_q)


        IF ( eig_k(1) .LT. e_fermi  .AND. e_fermi .LT. eig_k(2) ) THEN 
            write(4000,*), "CASE 1"
            !case=1
            prefac= ( e_fermi - eig_q(1) ) / ( e_fermi - eig_k(1) )
            ! temporary arrray for energy quotient
            tmp_arr(1) = (eig_q(2) - eig_q(1)) / (eig_k(2) - eig_k(1)) 
            tmp_arr(2) = (eig_q(3) - eig_q(1)) / (eig_k(3) - eig_k(1)) 
            tmp_arr(3) = (eig_q(4) - eig_q(1)) / (eig_k(4) - eig_k(1)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (el_ph_mat(2) -el_ph_mat(1)) / (eig_k(2) - eig_k(1))
            tmp_mat(2)= (el_ph_mat(3) - el_ph_mat(1)) / (eig_k(3) - eig_k(1)) 
            tmp_mat(3)= (el_ph_mat(4) - el_ph_mat(1)) / (eig_k(4) - eig_k(1)) 


            call intersection_surfaces(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)
            
            f =  (e_fermi - eig_k(1))/ ( (eig_k(2) - eig_k(1)) * (eig_k(3) - eig_k(1)) * (eig_k(4) - eig_k(1)) )

            interpol_mat  = el_ph_mat(1) + 0.5 * (e_fermi - eig_k(1)) * interpol_mat 


            area = 6*voltet*f*intersection_val * interpol_mat 
        

        ELSE IF ( eig_k(2) .LT. e_fermi  .AND. e_fermi .LT. eig_k(3) ) THEN 
            !case=2 
            !!!
            !!! This contribution I_0
            !!!
            write(4000,*), "CASE 2"
            prefac= ( e_fermi - eig_q(1) ) / ( e_fermi - eig_k(1) )
            tmp_arr(1) = (eig_q(2) - eig_q(1)) / (eig_k(2) - eig_k(1)) 
            tmp_arr(2) = (eig_q(3) - eig_q(1)) / (eig_k(3) - eig_k(1)) 
            tmp_arr(3) = (eig_q(4) - eig_q(1)) / (eig_k(4) - eig_k(1)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (el_ph_mat(2) -el_ph_mat(1)) / (eig_k(2) - eig_k(1))
            tmp_mat(2)= (el_ph_mat(3) - el_ph_mat(1)) / (eig_k(3) - eig_k(1)) 
            tmp_mat(3)= (el_ph_mat(4) - el_ph_mat(1)) / (eig_k(4) - eig_k(1)) 

            call intersection_surfaces(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)

            f =  (e_fermi - eig_k(1))/ ( (eig_k(2) - eig_k(1)) * (eig_k(3) - eig_k(1)) * (eig_k(4) - eig_k(1)) )

            interpol_mat  = el_ph_mat(1) + 0.5 * (e_fermi - eig_k(1)) * interpol_mat 


            !!!
            !!! This is contribution I_1
            !!!
            prefac= ( e_fermi - eig_q(2) ) / ( e_fermi - eig_k(2) )
            tmp_arr(1) = (eig_q(1) - eig_q(2)) / (eig_k(1) - eig_k(2)) 
            tmp_arr(2) = (eig_q(3) - eig_q(2)) / (eig_k(3) - eig_k(2)) 
            tmp_arr(3) = (eig_q(4) - eig_q(2)) / (eig_k(4) - eig_k(2)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (el_ph_mat(1) -el_ph_mat(2)) / (eig_k(1) - eig_k(2))
            tmp_mat(2)= (el_ph_mat(3) - el_ph_mat(2)) / (eig_k(3) - eig_k(2)) 
            tmp_mat(3)= (el_ph_mat(4) - el_ph_mat(2)) / (eig_k(4) - eig_k(2)) 

            call intersection_surfaces(tmp_arr,prefac,intersection_val2,tmp_mat,interpol_mat2)

            f2= (e_fermi - eig_k(2))/ ( (eig_k(2) - eig_k(1)) * (eig_k(3) - eig_k(2)) * (eig_k(4) - eig_k(2)) )

            interpol_mat2  = el_ph_mat(2) + 0.5 * (e_fermi - eig_k(2)) * interpol_mat2 

            area = 6 * voltet * (f*intersection_val * interpol_mat - f2*intersection_val2 * interpol_mat2)


        ELSE IF ( eig_k(3) .LT. e_fermi  .AND. e_fermi .LT. eig_k(4) ) THEN 
            !case=3 
            write(4000,*), "CASE 3"
            prefac= ( e_fermi - eig_q(4) ) / ( e_fermi - eig_k(4) )
            tmp_arr(1) = (eig_q(1) - eig_q(4)) / (eig_k(1) - eig_k(4)) 
            tmp_arr(2) = (eig_q(2) - eig_q(4)) / (eig_k(2) - eig_k(4)) 
            tmp_arr(3) = (eig_q(3) - eig_q(4)) / (eig_k(3) - eig_k(4)) 

            ! temporary array for matrix element interpolation 
            tmp_mat(1)= (el_ph_mat(1) -el_ph_mat(4)) / (eig_k(1) - eig_k(4))
            tmp_mat(2)= (el_ph_mat(2) - el_ph_mat(4)) / (eig_k(2) - eig_k(4)) 
            tmp_mat(3)= (el_ph_mat(3) - el_ph_mat(4)) / (eig_k(3) - eig_k(4)) 

            call intersection_surfaces(tmp_arr,prefac,intersection_val,tmp_mat,interpol_mat)

            f= (eig_k(4)  - e_fermi )/ ( (eig_k(4) - eig_k(1)) * (eig_k(4) - eig_k(2)) * (eig_k(4) - eig_k(3)) )

            interpol_mat = el_ph_mat(3) + 0.5* (e_fermi - eig_k(3)) *interpol_mat

            area = 6 * voltet*f*intersection_val * interpol_mat

        ELSE 
            area = 0 
        END IF

        k_int = k_int + area

    END SUBROUTINE tetra_area


    SUBROUTINE intersection_surfaces(arr_sort,prefac,intersection_val, arr_mat,interpol_mat) 
        
        !!! Calculates the intersection line between the k+q and k hypersurface
        !!! Interpolates the matrix element at this intersection line 
        REAL,INTENT(INOUT) :: arr_sort(3)
        REAL,INTENT(IN)    :: prefac
        REAL,INTENT(OUT)   :: intersection_val
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
        END IF 



    END SUBROUTINE intersection_surfaces

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