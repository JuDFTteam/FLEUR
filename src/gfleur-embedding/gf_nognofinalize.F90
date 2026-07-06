!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_nognofinalize 
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_nognofinalize(g,layer,lapw,lapw_gf,single_en)
!**********************************************                         
!     Perform the final matrix multiplications                          
!     with the matrix of the eigenvectors ("uhumatrix").                
!     This finally turns g into the energy sum over                     
!     the Green's functions. For all equilibrium energies               
!     the "free" Green's function has to be added, which                
!     is given by the energy denominator sandwiched between             
!     the matrix of the eigenvectors.                                   
!     May 2007, Frank Freimuth                                          
!**********************************************                         
      use m_juDFT 
      USE m_gf_energies,ONLY:gf_z,direction,gf_weightz,gf_noen 
      USE m_gf_spectrum
      USE m_gf_types
      IMPLICIT NONE 
      TYPE(t_lapw),INTENT(IN)::lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      COMPLEX,INTENT(INOUT) :: g(:,:) 
      INTEGER,INTENT(IN)    :: layer
      INTEGER,INTENT(IN),OPTIONAL::single_en
      COMPLEX,ALLOCATABLE   :: gtemp(:,:),uhu(:,:)
      COMPLEX,ALLOCATABLE   :: enedenom(:)

      INTEGER matsize,en,j,en_min,en_max
      COMPLEX z 
                                                                        
      matsize=size(g,1) 
      IF(size(uhumatrix,1)<matsize)CALL juDFT_error("uhumatrix")
      ALLOCATE( uhu(matsize,matsize) )
      IF(size(uhumatrix,1)==matsize) THEN
            uhu=uhumatrix
      ELSE
            CALL juDFT_error("BUG in sphere code",calledby="gf_nognofinalize")
            uhu(:lapw_gf%nv_tot_sphere,:)=uhumatrix(:lapw_gf%nv_tot_sphere,:matsize)
            if (matsize>lapw_gf%nv_tot_sphere) then !copy LO'part
                uhu(lapw_gf%nv_tot_sphere+1:,:)=uhumatrix(lapw%nv_tot+1:,:matsize)
            endif
      endif


      ALLOCATE( gtemp(matsize,matsize) )

      ALLOCATE( enedenom(matsize)) 
      enedenom=cmplx(0.0,0.0) 

      en_min=1
      en_max=gf_noen()

      if (present(single_en)) THEN
         en_min=single_en
         if (single_en==0) en_min=1
         en_max=single_en
      endif

      DO en=en_min,en_max
                                   !ignore the non-equilibrium energies

         IF(direction(en)/=0)CYCLE
         z=gf_z(en,layer)
         IF (present(single_en)) THEN
           enedenom(:)=enedenom(:)+1.0/(z-cmplx(uhueigval(:),0.0))
         ELSE
           enedenom(:)=enedenom(:)+gf_weightz(en)/(z-cmplx(uhueigval(:),0.0))
         ENDIF
      ENDDO

      DO j=1,matsize 
         g(j,j)=g(j,j)+enedenom(j) 
      ENDDO 


      CALL zgemm('N','C',matsize,matsize,matsize,              &
     &      cmplx(1.0,0.0),g,matsize,uhu,                         &
     &      matsize,cmplx(0.0,0.0),gtemp,matsize)                       
      CALL zgemm('N','N',matsize,matsize,matsize,              &
     &      cmplx(1.0,0.0),uhu,matsize,gtemp,                     &
     &      matsize,cmplx(0.0,0.0),g,matsize)                           
                                                                        
      END SUBROUTINE gf_nognofinalize 
      END                                           
