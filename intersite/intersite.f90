MODULE m_intersite

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_intersite
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module calculates the intersite green's function (WIP)
   !
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   !------------------------------------------------------------------------------
   USE m_juDFT
   USE m_types

   CONTAINS

   SUBROUTINE intersite_coeffs(kpt,atoms,sym,cell,input,ispin,noccbd,tetweights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)

      USE m_constants
      USE m_sel_sites

      IMPLICIT NONE

      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_cell),           INTENT(IN)     :: cell
      TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
      TYPE(t_usdus),          INTENT(IN)     :: usdus
      TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs
      TYPE(t_input),          INTENT(IN)     :: input

      INTEGER,                INTENT(IN)     :: ispin
      INTEGER,                INTENT(IN)     :: noccbd

      REAL,                   INTENT(IN)     :: kpt(3)
      REAL,                   INTENT(IN)     :: wtkpt
      REAL,                   INTENT(IN)     :: tetweights(:,:)
      INTEGER,                INTENT(IN)     :: ind(:,:)
      REAL,                   INTENT(IN)     :: eig(noccbd)

      INTEGER  i_gf,ntype,n,np,i_np,n_nearest,nn,l,lp,lm,lmp,m,mp,j,i
      REAL     wk, fac
      LOGICAL  l_zero
      COMPLEX  phase,cil,tmp

      INTEGER :: at(27*atoms%nat)
      REAL    :: dist(3,27*atoms%nat)
      wk = wtkpt/greensfCoeffs%del

      !Loop through the atoms
      DO nType = 1, atoms%ntype 
         n = SUM(atoms%neq(:ntype-1))
         DO nn = 1, atoms%neq(ntype)
            fac = 1.0/atoms%neq(ntype)
            n = n+1
            !Find the nearest neighbours to the current atom 
            CALL n_neighbours(n,atoms,cell,2,n_nearest,dist,at)
            DO i_np = 1, n_nearest
               np = at(i_np)
               !Bloch phase:
               phase = exp(ImagUnit*dot_product(kpt,dist(:,i_np)))

               DO i = 1, noccbd
                  l_zero = .true.
                  IF(input%tria) THEN
                     !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                     IF(ANY(tetweights(ind(i,1):ind(i,2),i).NE.0.0)) l_zero = .false.
                  ELSE
                     !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
                     j = NINT((eig(i)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
                     IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) ) l_zero = .false.
                  END IF

                  IF(l_zero) CYCLE

                  DO l = 0, lmaxU_const
                     DO lp = 0, lmaxU_const
                        cil = ImagUnit**(l-lp)
                        DO m = -l,l
                           lm = l*(l+1)+m
                           DO mp = -lp,lp
                              lmp = lp*(lp+1)+mp
                              IF(input%tria) THEN
                                 !We need to differentiate the weights with respect to energy (can maybe be done analytically)
                                 DO j = ind(i,1), ind(i,2)
                                    tmp = cil * phase * conjg(eigVecCoeffs%acof(i,lm,n,ispin))*eigVecCoeffs%acof(i,lmp,np,ispin)
                                    greensfCoeffs%uu_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%uu_int(j,n,np,m,mp,ispin) - fac * tetweights(j,i) * pi_const * REAL(tmp)
                                    
                                    tmp = cil * phase * conjg(eigVecCoeffs%bcof(i,lm,n,ispin))*eigVecCoeffs%bcof(i,lmp,np,ispin)                                                                    
                                    greensfCoeffs%dd_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%dd_int(j,n,np,m,mp,ispin) - fac * tetweights(j,i) * pi_const * REAL(tmp)
                                    
                                    tmp = cil * phase * conjg(eigVecCoeffs%acof(i,lm,n,ispin))*eigVecCoeffs%bcof(i,lmp,np,ispin)                                                                
                                    greensfCoeffs%du_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%du_int(j,n,np,m,mp,ispin) - fac * tetweights(j,i) * pi_const * REAL(tmp)
                                                                                                   
                                    tmp = cil * phase * conjg(eigVecCoeffs%bcof(i,lm,n,ispin))*eigVecCoeffs%acof(i,lmp,np,ispin)
                                    greensfCoeffs%ud_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%ud_int(j,n,np,m,mp,ispin) - fac * tetweights(j,i) * pi_const * REAL(tmp)
                                 ENDDO
                              ELSE
                                 tmp = cil * phase * conjg(eigVecCoeffs%acof(i,lm,n,ispin))*eigVecCoeffs%acof(i,lmp,np,ispin)
                                 greensfCoeffs%uu_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%uu_int(j,n,np,m,mp,ispin) - fac * wk * pi_const * REAL(tmp)
                                 
                                 tmp = cil * phase * conjg(eigVecCoeffs%bcof(i,lm,n,ispin))*eigVecCoeffs%bcof(i,lmp,np,ispin)                                                                    
                                 greensfCoeffs%dd_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%dd_int(j,n,np,m,mp,ispin) - fac * wk * pi_const * REAL(tmp)
                                 
                                 tmp = cil * phase * conjg(eigVecCoeffs%acof(i,lm,n,ispin))*eigVecCoeffs%bcof(i,lmp,np,ispin)                                                                
                                 greensfCoeffs%du_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%du_int(j,n,np,m,mp,ispin) - fac * wk * pi_const * REAL(tmp)
                                                                                                
                                 tmp = cil * phase * conjg(eigVecCoeffs%bcof(i,lm,n,ispin))*eigVecCoeffs%acof(i,lmp,np,ispin)
                                 greensfCoeffs%ud_int(j,n,np,lm,lmp,ispin) = greensfCoeffs%ud_int(j,n,np,m,mp,ispin) - fac * wk * pi_const * REAL(tmp)
                              END IF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE intersite_coeffs

   !SUBROUTINE prepare_intersite(atoms,vr,jspins,gf,ef,greensfCoeffs,sym,lmax)
!
   !   !sets up the information for the intersite green's function in the green's function type 
   !   !and calculates the radial function R_l, which is the KKR-radial function
!
   !   USE m_radsra
   !   USE m_constants
!
   !   IMPLICIT NONE
!
   !   TYPE(t_atoms),          INTENT(IN)     :: atoms
   !   TYPE(t_greensf),        INTENT(INOUT)  :: gf
   !   TYPE(t_greensfCoeffs),  INTENT(IN)     :: greensfCoeffs
   !   TYPE(t_sym),            INTENT(IN)     :: sym
!
   !   INTEGER,                INTENT(IN)     :: jspins
   !   INTEGER,                INTENT(IN)     :: lmax
!
   !   REAL,                   INTENT(IN)     :: ef
   !   REAL,                   INTENT(IN)     :: vr(atoms%jmtd,atoms%ntype,jspins)
!
   !   INTEGER  n,jspin,j,l,nodeu
   !   REAL     e,c
!
   !   REAL :: u_mt(lmax),du_mt(lmax),alpha(lmax)
   !   
   !   IF(gf%l_onsite.OR..NOT.greensfCoeffs%l_intersite) CALL juDFT_error("Green's function not initialized for Intersite calculations",calledby="prepare_intersite")
!
   !   !
   !   !Give the Green#s function the information about the energy grid
   !   !
   !   gf%ne    = greensfCoeffs%ne
   !   gf%e_bot = greensfCoeffs%e_bot
   !   gf%e_top = greensfCoeffs%e_top
   !   gf%del   = greensfCoeffs%del 
   !   gf%sigma = greensfCoeffs%sigma 
   !   !
   !   !Write the coefficients to the corresponding arrays
   !   !
   !   CALL move_ALLOC(greensfCoeffs%uu_int,gf%uu_int)
   !   CALL move_ALLOC(greensfCoeffs%dd_int,gf%dd_int)
   !   CALL move_ALLOC(greensfCoeffs%du_int,gf%du_int)
   !   CALL move_ALLOC(greensfCoeffs%ud_int,gf%ud_int)
   !   !
   !   !Calculate the KKR-radial functions
   !   !
   !   DO n = 1, atoms%nat 
   !      DO jspin = 1, jspins
   !         DO j = 1, gf%ne 
   !            e = j*gf%del+gf%e_bot
   !            IF(e.LT.0) CYCLE
   !            DO l = 0, lmax
   !               !
   !               !Calculate the radial function R_l^\alpha(r,E)
   !               !
   !               c = c_light(1.0)
   !               !This calculates the normed radial function
   !               CALL radsra(e,l,vr(:,n,jspin),atoms%rmsh(1,n),atoms%dx(n),atoms%jri(n),&
   !                           c, u_mt(l),du_mt(l),nodeu,gf%R(j,:,1,l,jspin),gf%R(j,:,2,l,jspin))
   !            ENDDO
   !            !Now we calculate the factor resulting from the different boundary condition
   !            !CALL bound(e,atoms%rmsh(atoms%jri(n),n),alpha(:),u_mt(:),du_mt(:)/u_mt(:),lmax)
   !            !Now rescale R_l accordingly
   !            DO l = 0, lmax
   !               gf%R(j,:,:,l,jspin) = alpha(l) * gf%R(j,:,:,l,jspin) 
   !            ENDDO
   !         ENDDO
   !      ENDDO
   !   ENDDO
!
   !END SUBROUTINE prepare_intersite
!
   !SUBROUTINE calc_intersite(gf,jr,jrp,atoms,vr,enpara,sym,jspins,lmax,ef,gmmpMat)
   !      
   !      !calculates the intersite green's function from the previously calculated coefficients for r and r'
   !         
   !      USE m_types
   !      USE m_radfun
   !      USE m_constants
   !      USE m_kkintgr
   !      
   !      IMPLICIT NONE
!
   !      TYPE(t_greensf),       INTENT(IN)     :: gf 
   !      TYPE(t_atoms),          INTENT(IN)     :: atoms
   !      TYPE(t_enpara),         INTENT(IN)     :: enpara
   !      TYPE(t_sym),            INTENT(IN)     :: sym
!
   !      INTEGER,                INTENT(IN)     :: jspins
   !      INTEGER,                INTENT(IN)     :: jr,jrp
!
   !      REAL,                   INTENT(IN)     :: ef
   !      REAL,                   INTENT(IN)     :: vr(atoms%jmtd,atoms%ntype,jspins)
   !      COMPLEX,                INTENT(OUT)    :: gmmpMat(gf%ne,atoms%nat,atoms%nat,0:lmax**2+lmax,0:lmax**2+lmax,input%jspins)
!
!
   !      INTEGER        i_gf,n,jspin,l,lp,lm,lmp,m,mp,j,l_dim
   !      INTEGER        nodeu,noded
   !      REAL           e,c,wronk
   !      TYPE(t_usdus)  usdus
!
   !      !Arrays for LAPW radial functions
   !      REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
   !      REAL,    ALLOCATABLE :: im_gmmpMat(:,:,:,:,:)
   !      
   !      ldim = lmax**2 + lmax
   !      
   !      ALLOCATE (im_gmmpMat(gf%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
   !      im_gmmpMat = 0.0
   !      ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspins) )
   !      ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jspins) )
   !      CALL usdus%init(atoms,jspins)
!
   !      DO i_gf = 1, gf%n_gf
   !         n = gf%atomType(i_gf)
!
   !         !EQ.47
   !         !
   !         !  R_l^\alpha(r,E)ImG_{L,L'}^{\alpha,\alpha'}(E)R_l'^\alpha'(r',E) = R_l^\alpha(r,E)\sqrt(E)\Theta(E)\delta_{L,L'}^{\alpha,\alpha'}R_l'^\alpha'(r',E) + ImG_{L,L'}^{\alpha,\alpha'}(r,r',E)
   !         !
!
   !         !
   !         !Calculate the radial functions needed
   !         !
   !         CALL timestart("Inter-Site: Radial functions")
   !         DO jspin = 1, jspins
   !            DO l = 0, lmax
   !               !
   !               !Calculate the LAPW radial function u_l and \dot{u}_l
   !               !
   !               CALL radfun(l,n,jspin,enpara%el0(l,n,jspin),vr(:,n,jspin),atoms,f(:,:,l,jspin),g(:,:,l,jspin),usdus, nodeu,noded,wronk)
   !            ENDDO
   !         ENDDO
   !         CALL timestop("Inter-Site: Radial functions")
!
   !         CALL timestart("Inter-Site: Imaginary part")
   !         DO jspin = 1, jspins
   !            DO l = 0, lmaxU_const
   !               DO lp = 0, lmaxU_const
   !                  DO j = 1, gf%ne
   !                     !
   !                     !ImG_{L,L'}^{\alpha,\alpha'}(r,r',E) = EQ.44
   !                     !
   !                     DO m = -l,l
   !                        lm = l*(l+1)+m
   !                        DO mp = -lp,lp
   !                           lmp = lp*(lp+1)+mp
   !                           im_gmmpMat(j,i_gf,lm,lmp,jspin) = im_gmmpMat(j,i_gf,lm,lmp,jspin) + &
   !                                                            gf%uu_int(:,i_gf,lm,lmp,jspin) * (f(jr,1,l,jspin)*f(jrp,1,lp,jspin)+f(jr,2,l,jspin)*f(jrp,2,lp,jspin)) +&
   !                                                            gf%dd_int(:,i_gf,lm,lmp,jspin) * (g(jr,1,l,jspin)*g(jrp,1,lp,jspin)+g(jr,2,l,jspin)*g(jrp,2,lp,jspin)) +&
   !                                                            gf%du_int(:,i_gf,lm,lmp,jspin) * (f(jr,1,l,jspin)*g(jrp,1,lp,jspin)+f(jr,2,l,jspin)*g(jrp,2,lp,jspin)) +&
   !                                                            gf%ud_int(:,i_gf,lm,lmp,jspin) * (g(jr,1,l,jspin)*f(jrp,1,lp,jspin)+g(jr,2,l,jspin)*f(jrp,2,lp,jspin))
!
   !                        ENDDO
   !                     ENDDO
   !                     !
   !                     !R_l^\alpha(r,E)\sqrt(E)\Theta(E)\delta_{L,L'}^{\alpha,\alpha'}R_l'^\alpha'(r',E)
   !                     !
   !                     e = j*gf%del+gf%e_bot
   !                     IF(e.GT.0.AND.l.EQ.lp) THEN
   !                        DO m = -l, l
   !                           lm = l*(l+1) + m
   !                           im_gmmpMat(j,i_gf,lm,lm,jspin) = im_gmmpMat(j,i_gf,lm,lm,jspin) + sqrt(e) * &
   !                                                            (gf%R(j,jr,1,l,jspin)*gf%R(j,jrp,1,l,jspin)+gf%R(j,jr,2,l,jspin)*gf%R(j,jrp,2,l,jspin))
   !                        ENDDO
   !                     ENDIF 
   !                  ENDDO
   !               ENDDO
   !            ENDDO
   !         ENDDO
   !         CALL timestop("Inter-Site: Imaginary part")
   !         !
   !         !taking care of spin degeneracy in non-magnetic case
   !         !
   !         IF(jspins.EQ.1) gf%im_gmmpMat(:,i_gf,:,:,1) = 2.0 * gf%im_gmmpMat(:,i_gf,:,:,1)
   !         !
   !         ! Calculate DOS
   !         !
!
   !         !
   !         ! Determine energy cutoff
   !         !
!
   !         !
   !         !EQ.50: R_l^\alpha(r,z)ImG_{L,L'}^{\alpha,\alpha'}(z)R_l'^\alpha'(r',z) = 1/pi int_E_b^E_c dE' 1/(E'-z) R_l^\alpha(r,E)ImG_{L,L'}^{\alpha,\alpha'}(E)R_l'^\alpha'(r',E)
   !         !
   !         !CALL timestart("On-Site: Kramer-Kronigs-Integration")
   !         !DO jspin = 1, jspins
   !         !   DO l = 0, lmaxU_const
   !         !      DO lp = 0, lmaxU_const
   !         !         DO m = -l,l
   !         !            lm = l*(l+1)+m
   !         !            DO mp = -lp,lp
   !         !               lmp = lp*(lp+1)+mp
   !         !               CALL kkintgr_complex(gf%nz,gf%e(:),gf%ne,gf%sigma,gf%del,gf%e_bot,&
   !         !                           im_gmmpMat(:,i_gf,lm,lmp,jspin),gmmpMat(:,i_gf,lm,lmp,jspin))
!!
   !         !            ENDDO
   !         !         ENDDO
   !         !      ENDDO
   !         !   ENDDO
   !         !ENDDO
   !         !CALL timestop("On-Site: Kramer-Kronigs-Integration")
   !      ENDDO
!
   !   END SUBROUTINE calc_intersite

END MODULE m_intersite