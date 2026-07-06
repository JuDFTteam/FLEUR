!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_ab_coef 
      IMPLICIT NONE
      private
      complex,allocatable,save :: ab_coef(:,:,:,:,:)
      integer,save             :: ab_spin
      public:: gf_ab_coef_calc,gf_ab_coef_delete,gf_ab_coef_matrix,gf_ab_coef_vector
      CONTAINS 
      SUBROUTINE gf_ab_coef_calc(l_noco,jspin,bk,sym,el,vr0,atoms,cell,lapw)
      !Calculate the ab_coefs and store in a local module variable
      USE m_gf_types
	  !Arguments
      logical,intent(in)       :: l_noco
      INTEGER,INTENT(IN)       :: jspin
      REAL,     INTENT(IN)     :: bk(3)
      REAL,INTENT(IN)          :: vr0(:,:,:)
      REAL,INTENT(IN)          :: el(0:,:,:)
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_lapw),INTENT(IN)  :: lapw
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_atoms),INTENT(IN) :: atoms
	  !Locals
      if (l_noco) THEN
           ab_spin=2
      else
      	   ab_spin=1
      endif
	  ALLOCATE(ab_coef(maxval(lapw%nv_sphere),0:MAXVAL(atoms%lmax)*(MAXVAL(atoms%lmax)+2),atoms%nat,2,ab_spin))

	  if (ab_spin==1) then
	      call priv_ab_coef(jspin,bk,sym,el(:,:,jspin),vr0(:,:,jspin),atoms,cell,lapw,  &
	                        ab_coef(:,:,:,:,1))
	  else
	      call priv_ab_coef(1,bk,sym,el(:,:,1),vr0(:,:,1),atoms,cell,lapw,  &
	                        ab_coef(:,:,:,:,1))
	      call priv_ab_coef(2,bk,sym,el(:,:,2),vr0(:,:,2),atoms,cell,lapw,  &
	                        ab_coef(:,:,:,:,2))
	  endif
      END SUBROUTINE
	  subroutine gf_ab_coef_delete()
	  !simple subroutine to clean memory
	  if (allocated(ab_coef)) deallocate(ab_coef)
	  end subroutine

	  function gf_ab_coef_matrix(lmpmax,minat,maxat,jspin)
	  !returns the full matrix of ab-coefficients for the current spin
	  integer,intent(in):: jspin,lmpmax,minat,maxat
	  complex           :: gf_ab_coef_matrix(size(ab_coef,1),0:lmpmax-1,maxat-minat+1,size(ab_coef,4))

	  if (ab_spin==1) then
	  	gf_ab_coef_matrix=ab_coef(:,0:lmpmax-1,minat:maxat,:,1)
	  else
	  	gf_ab_coef_matrix=ab_coef(:,0:lmpmax-1,minat:maxat,:,jspin)
      endif
      end function

      function gf_ab_coef_vector(lmp,nt,i,jspin)
      !returns only the vector of ab_coef for all g
      integer,intent(in) :: lmp,nt,i,jspin
      complex           :: gf_ab_coef_vector(size(ab_coef,1))
      !locals
      integer :: jsp

      jsp=jspin
      if (ab_spin==1) jsp=1

      gf_ab_coef_vector=ab_coef(:,lmp,nt,i,jsp)
      end function

      SUBROUTINE priv_ab_coef(jspin,bk,sym,                               &
     &     el,vr0,atoms,cell,lapw,ab)                                   
!*********************************************************************  
!     Returns the LAPW a+b coefs for the given k-point.                 
!     WARNING: array ab is pretty large                                 
!                                           Daniel Wortmann             
!*********************************************************************  
      USE m_fleur_interface,ONLY:fleur_ylm,fleur_dsphbs,fleur_sphbes 
      USE m_constants ,ONLY: pimach 
      USE m_gf_types 
      USE m_fleur_interface,ONLY:fleur_radfun 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN)       :: jspin 
                                              !k-point                  
      REAL,     INTENT(IN)     :: bk(3) 
                                              ! sph. potential          
      REAL,INTENT(IN)          :: vr0(:,:) 
                                              ! enparas                 
      REAL,INTENT(IN)          :: el(0:,:) 
      TYPE(t_sym),INTENT(IN)   :: sym 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_cell),INTENT(IN)  :: cell 
      TYPE(t_atoms),INTENT(IN) :: atoms 
      COMPLEX, INTENT(OUT)     :: ab(:,0:,:,:) 
                                                                        
!     locals                                                            
                                             !loops                     
      INTEGER :: n,lm,l,m,nat,ntyp,na,nap 
                                             !constants                 
      REAL    :: tpi 
                                             !k-vectors                 
      REAL    :: fkr(3),fkp(3),fk(3) 
      COMPLEX :: ylm((size(el,1))**2) 
                                         !i**l                          
      COMPLEX :: ipower(0:size(el,1)-1) 
      REAL    :: fj(0:size(el,1)-1),dfj(0:size(el,1)-1) 
      REAL    :: r1,const,df 
      COMPLEX :: wronk,phase 
                                                                        
                                !for radfun                             
      INTEGER :: noded,nodeu 
                                !not used                               
      REAL    :: rdum 
                                        !not used afterwards            
      REAL,ALLOCATABLE :: f(:,:),g(:,:) 
                                                      !MT-lapws         
      REAL,ALLOCATABLE :: us(:),dus(:),duds(:),uds(:) 
                                                                        
      ALLOCATE(f(maxval(atoms%jri),2),g(maxval(atoms%jri),2)) 
      ALLOCATE(us(0:maxval(atoms%lmax)),dus(0:maxval(atoms%lmax)),      &
     &     duds(0:maxval(atoms%lmax)),uds(0:maxval(atoms%lmax)))
                                                                        
!  Some constants                                                       
      tpi=2*pimach() 
      const=4*pimach()/SQRT(cell%omtil) 
      DO l=0,maxval(atoms%lmax)
         ipower(l)=CMPLX(0.0,1.0)**l 
      ENDDO 
                                                                        
      nat=0 
      DO ntyp=1,atoms%ntype 
!                                                                       
!   The radial functions                                                
!                                                                       
         DO l=0,atoms%lmax(ntyp)
            CALL fleur_radfun(                                          &
     &           l,el(l,ntyp),vr0(:,ntyp),ntyp,atoms,f,g,us(l)          &
     &           ,dus(l),uds(l),duds(l),rdum,nodeu,noded,rdum)          
         ENDDO 
         DO na=1,atoms%neq(ntyp) 
            nat=nat+1 
            DO n=1,lapw%nv_sphere(jspin)
               fk(1)=bk(1)+lapw%k%k1(n,jspin) 
               fk(2)=bk(2)+lapw%k%k2(n,jspin) 
               fk(3)=bk(3)+lapw%k%k3(n,jspin) 
               phase =EXP(CMPLX(0.0, tpi* (fk(1)*atoms%taual(1,nat)+    &
     &              fk(2)*atoms%taual(2,nat)+                           &
     &              fk(3)*atoms%taual(3,nat))))                         
               nap=sym%invtab(atoms%ngopr(nat)) 
               fkr=MATMUL(fk,real(sym%mrot(:,:,nap))) 
               fkp=MATMUL(fkr,cell%bmat) 
               !CALL cotra3(fkr,fkp,cell%bmat)                          
!     ----> generate spherical harmonics                                
               CALL fleur_ylm(                                          &
     &              maxval(atoms%lmax),fkp,                             &
     &              ylm)                                                
               r1 = atoms%rmt(ntyp)*lapw%k%rk(n,jspin) 
               CALL fleur_sphbes(                                       &
     &              maxval(atoms%lmax),r1,                              &
     &              fj)                                                 
               CALL fleur_dsphbs(                                       &
     &              maxval(atoms%lmax),r1,fj,                           &
     &              dfj)                                                
               DO  l = 0,atoms%lmax(ntyp)
                  df = lapw%k%rk(n,jspin)*dfj(l) 
                  wronk=uds(l)*dus(l)-                                  &
     &                 us(l)*duds(l)                                    
                                    !This disabled code has to be used f
                  IF (.FALSE.) THEN 
                     fj(l) = 1.0*const * fj(l)/us(l) 
                     dfj(l) = 0.0d0 
                  ELSE 
                     dfj(l)=const*(dus(l)*fj(l)-df*us(l))               &
     &                    /wronk                                        
                     fj(l)=const*(df*uds(l)-fj(l)*duds(l))              &
     &                    /wronk                                        
                  ENDIF 
                  DO m=-l,l 
                     lm = l*(l+1) + m 
                     ab(n,lm,nat,1)=ipower(l)*phase*                    &
     &                    CONJG(ylm(lm+1))*fj(l)                        
                     ab(n,lm,nat,2)=ipower(l)*phase*                    &
     &                    CONJG(ylm(lm+1))*dfj(l)                       
                  ENDDO 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
      DEALLOCATE(f,g,us,dus,duds,uds) 
      END SUBROUTINE 
                                                                        
      END                                           
