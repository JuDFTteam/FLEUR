!--------------------------------------------------------------------------------
! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_potden_offset
    USE m_juDFT
#ifdef CPP_MPI
    use mpi
#endif

    CONTAINS 
    
    SUBROUTINE dfpt_potden_offset(ispin, fmpi, stars, cell, atoms, potden, potdenImag,l_correct,l_den,offset_out)
        USE m_types
        USE m_constants
        USE m_convol
        USE m_intgr, ONLY: intgr3

        IMPLICIT NONE

        INTEGER,           INTENT(IN)       :: ispin
        TYPE(t_mpi),       INTENT(IN)       :: fmpi
        TYPE(t_stars),     INTENT(IN)       :: stars        
        TYPE(t_cell),      INTENT(IN)       :: cell
        TYPE(t_atoms),     INTENT(IN)       :: atoms
        TYPE(t_potden),    INTENT(INOUT)    :: potden
        TYPE(t_potden),    INTENT(INOUT)    :: potdenImag
        LOGICAL,           INTENT(IN)       :: l_correct
        LOGICAL,           INTENT(IN)       :: l_den
        COMPLEX,OPTIONAL,  INTENT(OUT)      :: offset_out

        COMPLEX, ALLOCATABLE                :: dummy_w(:) , dummy_const(:) , dummy_const_w(:)
        COMPLEX                             :: int_ir, int_tot, int_mt_tot , vol_ir  , offset , sum_ints
        REAL, ALLOCATABLE                   :: mt_r(:,:) , mt_Im(:,:) , int_mt_r(:) , int_mt_Im(:) 
        INTEGER                             :: iType , f_write , count
        REAL                                :: vol_mts
        LOGICAL                             :: l_do

#ifdef CPP_MPI
        integer:: ierr
#endif

        ALLOCATE(dummy_w,mold = potden%pw(:,ispin))        
        ALLOCATE(dummy_const,mold = potden%pw(:,ispin)) 
        ALLOCATE(dummy_const_w,mold = potden%pw(:,ispin))        
        ALLOCATE(mt_r( atoms%jmtd,atoms%ntype))
        ALLOCATE(mt_Im(atoms%jmtd,atoms%ntype))
        ALLOCATE(int_mt_r( atoms%ntype))
        ALLOCATE(int_mt_Im( atoms%ntype))


        dummy_w = 0.0 
        dummy_const = 0.0 
        dummy_const_w = 0.0 
        
        mt_r = 0.0
        mt_Im = 0.0 
        int_mt_r = 0.0 
        int_mt_Im = 0.0 
        int_mt_tot = 0.0 

        vol_mts = 0.0 

        f_write = 4000
        IF (l_den) f_write = 4500 


        count = 0 
        l_do=.TRUE.

        dummy_const(1) = 1.0 

        CALL timestart("DFPT potden offset correction")

        IF ( fmpi%irank==0 .AND. norm2(stars%center) < 1e-8 ) THEN 
            ! IR integral
            ! Only G=0 survives
            CALL convol(stars, dummy_w,potden%pw(:,ispin))
            int_ir = dummy_w(1) * cell%omtil

            ! Convolution for the constant shift we will add onto G=0
            CALL convol(stars, dummy_const_w, dummy_const(:))
            vol_ir = dummy_const_w(1) * cell%omtil

            ! MT part
            ! Only l=0 survives
            DO iType = 1, atoms%ntype
                IF (.NOT. l_den ) THEN
                    mt_r(1:atoms%jri(iType),iType) = atoms%rmsh( 1:atoms%jri(iType),iType) **2 * potden%mt(1:atoms%jri(iType),0,iType,ispin)
                    mt_Im(1:atoms%jri(iType),iType) = atoms%rmsh( 1:atoms%jri(iType),iType) **2 * potdenImag%mt(1:atoms%jri(iType),0,iType,ispin)
                ELSE 
                    mt_r(1:atoms%jri(iType),iType) = potden%mt(1:atoms%jri(iType),0,iType,ispin) 
                    mt_Im(1:atoms%jri(iType),iType) = potdenImag%mt(1:atoms%jri(iType),0,iType,ispin)
                END IF 
                CALL intgr3(mt_r(:,iType),atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),int_mt_r(iType))
                CALL intgr3(mt_Im(:,iType),atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),int_mt_Im(iType))
                
                vol_mts = vol_mts + atoms%volmts(iType)

                int_mt_tot = int_mt_tot + sfp_const* CMPLX(int_mt_r(iType),int_mt_Im(iType))
                write(f_write+2,*)  "Atom Type iDtype", iType , "Value", sfp_const*CMPLX(int_mt_r(iType),int_mt_Im(iType))
            END DO 

            sum_ints =  int_ir + int_mt_tot       
            offset =  (int_ir + int_mt_tot)/ (cell%omtil)
            
            IF ( l_correct) THEN
                potden%pw(1,ispin) = potden%pw(1,ispin) - offset
                ! this only works for densities atm 
                IF (l_den) THEN
                    DO iType = 1, atoms%ntype 
                        potden%mt(1:atoms%jri(iType),0,iType,ispin) = potden%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                           - sfp_const * REAL(offset) * atoms%rmsh( 1:atoms%jri(iType),iType) **2
                        potdenImag%mt(1:atoms%jri(iType),0,iType,ispin) = potdenImag%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                               - sfp_const* AIMAG(offset) * atoms%rmsh( 1:atoms%jri(iType),iType) **2
                    END DO 
                ELSE
                    DO iType = 1, atoms%ntype 
                        potden%mt(1:atoms%jri(iType),0,iType,ispin) = potden%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                           - sfp_const * REAL(offset) 
                        potdenImag%mt(1:atoms%jri(iType),0,iType,ispin) = potdenImag%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                               - sfp_const* AIMAG(offset)
                    END DO 
                END IF
            END IF 
            
            IF ( l_correct ) THEN 
                write(f_write,*)  "Offset that was corrected for", offset
                write(f_write,*)  "Sum", int_mt_tot + int_ir
            ELSE 
                write(f_write+1,*) "Sum", int_mt_tot + int_ir
                write(f_write+1,*) "Int ir", int_ir 
                write(f_write+1,*) "Int mt", int_mt_tot
            END IF 
        END IF 

        IF (PRESENT(offset_out)) offset_out = offset
        IF (l_correct) THEN
            call potden%distribute(fmpi%mpi_comm)
            call potdenImag%distribute(fmpi%mpi_comm)
        END IF 
!#ifdef CPP_MPI
!        CALL MPI_BCAST( potden%pw, size(potden%pw), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )        
!#endif

        CALL timestop("DFPT potden offset correction")

    END SUBROUTINE dfpt_potden_offset

    SUBROUTINE dfpt_surface_offset(ispin, fmpi, starsq,  stars, cell, atoms, potden1, potden1Imag, potden, grpotden, iDtype , l_correct,l_den, offset_out)
        USE m_types
        USE m_constants
        USE m_convol
        USE m_intgr, ONLY: intgr3

        IMPLICIT NONE

        INTEGER,           INTENT(IN)       :: ispin
        TYPE(t_mpi),       INTENT(IN)       :: fmpi
        TYPE(t_stars),     INTENT(IN)       :: starsq,stars        
        TYPE(t_cell),      INTENT(IN)       :: cell
        TYPE(t_atoms),     INTENT(IN)       :: atoms
        TYPE(t_potden),    INTENT(INOUT)    :: potden1 , potden1Imag
        TYPE(t_potden),    INTENT(IN)    :: potden , grpotden
        INTEGER,           INTENT(IN)       :: iDtype
        LOGICAL,           INTENT(IN)       :: l_correct
        LOGICAL,           INTENT(IN)       :: l_den
        COMPLEX,OPTIONAL,  INTENT(OUT)      :: offset_out

        COMPLEX, ALLOCATABLE                :: dummy_w(:) , dummy_const(:) , dummy_const_w(:)
        COMPLEX                             :: int_ir, int_tot, int_mt_tot , vol_ir  , offset , sum_ints
        REAL, ALLOCATABLE                   :: mt_r(:,:) , mt_Im(:,:) , int_mt_r(:) , int_mt_Im(:) 
        INTEGER                             :: iType , f_write , count
        REAL                                :: vol_mts
        LOGICAL                             :: l_do
#ifdef CPP_MPI
        integer:: ierr
#endif

        ALLOCATE(dummy_w,mold = potden1%pw(:,ispin))        
        ALLOCATE(dummy_const,mold = potden1%pw(:,ispin)) 
        ALLOCATE(dummy_const_w,mold = potden1%pw(:,ispin))        
        ALLOCATE(mt_r( atoms%jmtd,atoms%ntype))
        ALLOCATE(mt_Im(atoms%jmtd,atoms%ntype))
        ALLOCATE(int_mt_r( atoms%ntype))
        ALLOCATE(int_mt_Im( atoms%ntype))


        dummy_w = 0.0 
        dummy_const = 0.0 
        dummy_const_w = 0.0 
        count = 0     
        mt_r = 0.0
        mt_Im = 0.0 
        int_mt_r = 0.0 
        int_mt_Im = 0.0 
        int_mt_tot = 0.0 
        sum_ints = 0.0 
        vol_mts = 0.0 

        f_write = 4500
        IF (l_den) f_write = 4550 


        l_do = .TRUE.


        dummy_const(1) = 1.0 

        CALL timestart("DFPT vgen offset correction")

        IF ( fmpi%irank==0 .AND. norm2(starsq%center) < 1e-8 ) THEN 
            ! IR integral
            ! Only G=0 survives
            CALL dfpt_surface_convol(1,stars, starsq, potden1,potden,dummy_w)
            int_ir = dummy_w(1) * cell%omtil

            ! Convolution for the constant shift we will add onto G=0
            CALL convol(stars, dummy_const_w, dummy_const(:))
            vol_ir = dummy_const_w(1) * cell%omtil

            ! MT part
            ! Only l=0 survives
            DO iType = 1, atoms%ntype
                IF (.NOT. l_den ) THEN
                    mt_r(1:atoms%jri(iType),iType) = atoms%rmsh( 1:atoms%jri(iType),iType) **2 * potden1%mt(1:atoms%jri(iType),0,iType,ispin)
                    mt_Im(1:atoms%jri(iType),iType) = atoms%rmsh( 1:atoms%jri(iType),iType) **2 * potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin)
                    IF (iType == iDtype) THEN
                        mt_r(1:atoms%jri(iType),iType) = mt_r(1:atoms%jri(iType),iType) + atoms%rmsh( 1:atoms%jri(iType),iType) **2 * grpotden%mt(1:atoms%jri(iType),0,iType,ispin)
                        !mt_Im(1:atoms%jri(iType),iType) =  mt_Im(1:atoms%jri(iType),iType) + atoms%rmsh( 1:atoms%jri(iType),iType) **2 * %mt(1:atoms%jri(iType),0,iType,ispin)
                    END IF
                ELSE 
                    mt_r(1:atoms%jri(iType),iType) = potden1%mt(1:atoms%jri(iType),0,iType,ispin)
                    mt_Im(1:atoms%jri(iType),iType) = potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin)
                    IF (iType == iDtype) THEN
                        mt_r(1:atoms%jri(iType),iType) = mt_r(1:atoms%jri(iType),iType) + grpotden%mt(1:atoms%jri(iType),0,iType,ispin)
                    END IF 
                END IF 
                CALL intgr3(mt_r(:,iType),atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),int_mt_r(iType))
                CALL intgr3(mt_Im(:,iType),atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),int_mt_Im(iType))

                
                vol_mts = vol_mts + atoms%volmts(iType)
                

                int_mt_tot = int_mt_tot + sfp_const * CMPLX(int_mt_r(iType),int_mt_Im(iType)) 
                write(f_write+2,*)  "Atom Type iDtype", iType , "Value", sfp_const*CMPLX(int_mt_r(iType),int_mt_Im(iType))
            END DO 

            sum_ints = int_ir + int_mt_tot
            offset =  (int_ir + int_mt_tot)/ cell%omtil 
            
            IF ( l_correct) THEN
                potden1%pw(1,ispin) = potden1%pw(1,ispin) - offset
                IF (l_den) THEN
                    DO iType = 1, atoms%ntype 
                        potden1%mt(1:atoms%jri(iType),0,iType,ispin) = potden1%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                           - sfp_const * REAL(offset) * atoms%rmsh( 1:atoms%jri(iType),iType) **2
                        potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin) = potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                               - sfp_const* AIMAG(offset) * atoms%rmsh( 1:atoms%jri(iType),iType) **2
                    END DO 
                ELSE
                    DO iType = 1, atoms%ntype 
                        potden1%mt(1:atoms%jri(iType),0,iType,ispin) = potden1%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                           - sfp_const * REAL(offset) 
                        potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin) = potden1Imag%mt(1:atoms%jri(iType),0,iType,ispin) &
                        &                                               - sfp_const* AIMAG(offset)
                    END DO 
                END IF
            END IF 
            
            IF ( l_correct ) THEN 
                write(f_write,*)  "Offset that was corrected for", offset
                write(f_write,*)  "Sum", int_mt_tot + int_ir
            ELSE 
                !write(f_write+1,*) "Offset", offset , 
                write(f_write+1,*)"Sum", int_mt_tot + int_ir
                write(f_write+1,*) "Int ir", int_ir 
                write(f_write+1,*) "Int mt", int_mt_tot
            END IF 
        END IF 

        IF (PRESENT(offset_out)) offset_out = offset
        IF (l_correct) THEN 
            call potden1%distribute(fmpi%mpi_comm)
            call potden1Imag%distribute(fmpi%mpi_comm)
        END IF 
!#ifdef CPP_MPI
!        CALL MPI_BCAST( potden%pw, size(potden%pw), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )        
!#endif

        CALL timestop("DFPT vgen offset correction")

    END SUBROUTINE dfpt_surface_offset



    SUBROUTINE dfpt_surface_convol(ispin,stars,starsq,vTot1,vTot,pw_w)
        USE m_types
        USE m_constants
        USE m_fft3d

        IMPLICIT NONE 
        
        INTEGER, INTENT(IN)         :: ispin 
        TYPE(t_stars),INTENT(IN)    :: stars , starsq
        TYPE(t_potden), INTENT(IN)  :: vTot1,vTot
        COMPLEX,ALLOCATABLE,INTENT(INOUT) :: pw_w(:)
        

        INTEGER                         :: i, ifft3
        REAL,    ALLOCATABLE :: fftwork(:), vre(:), v1re(:), v1im(:)
        COMPLEX, ALLOCATABLE :: v1full(:)

       

        ifft3 = 27*stars%mx1*stars%mx2*stars%mx3 !TODO: What if starsq/=stars in that regard?

        ALLOCATE(vre(ifft3),v1re(ifft3),v1im(ifft3),v1full(ifft3),fftwork(ifft3))
    
        vre = 0.0; v1re = 0.0; v1im = 0.0
        CALL fft3d(v1re, v1im, vTot1%pw(:, ispin), starsq, +1)
        CALL fft3d(vre, fftwork, vTot%pw(:, ispin), stars, +1)
        v1full = CMPLX(v1re,v1im) * stars%ufft +  vre * starsq%ufft1!-q
        v1re =  REAL(v1full)
        v1im = AIMAG(v1full)
        CALL fft3d(v1re, v1im, pw_w, starsq, -1)
    

    END SUBROUTINE dfpt_surface_convol

END MODULE m_dfpt_potden_offset