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

   SUBROUTINE intersite_coeffs(kpt,atoms,sym,cell,input,ispin,nbands,tetweights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)

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
      INTEGER,                INTENT(IN)     :: nbands

      REAL,                   INTENT(IN)     :: kpt(3)
      REAL,                   INTENT(IN)     :: wtkpt
      REAL,                   INTENT(IN)     :: tetweights(:,:)
      INTEGER,                INTENT(IN)     :: ind(:,:)
      REAL,                   INTENT(IN)     :: eig(nbands)

      INTEGER  i_gf,ntype,n,np,i_np,n_nearest,nn,l,lp,lm,lmp,m,mp,j,ib
      REAL     fac,weight
      LOGICAL  l_zero
      COMPLEX  phase,cil,tmp

      INTEGER :: at(27*atoms%nat)
      REAL    :: dist(3,27*atoms%nat)

      !Loop through the atoms
      DO nType = 1, atoms%ntype 
         DO n = SUM(atoms%neq(:ntype-1)) + 1, SUM(atoms%neq(:ntype))
            fac = 1.0/atoms%neq(ntype)
            !Find the nearest neighbours to the current atom 
            CALL n_neighbours(n,atoms,cell,2,n_nearest,dist,at)
            DO i_np = 1, n_nearest
               np = at(i_np)
               !Bloch phase:
               phase = exp(ImagUnit*dot_product(kpt,dist(:,i_np)))

               DO ib = 1, nbands
                  l_zero = .true.
                  IF(input%tria) THEN
                     !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
                     IF(ANY(tetweights(ind(ib,1):ind(ib,2),ib).NE.0.0)) l_zero = .false.
                  ELSE
                     !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
                     j = NINT((eig(ib)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
                     IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) ) l_zero = .false.
                  END IF

                  IF(l_zero) CYCLE

                  DO ie = MERGE(ind(ib,1),j,input%tria), MERGE(ind(ib,2),j,input%tria)

                     weight = MERGE(tetweights(ie,ib),wtkpt/greensfCoeffs%del,input%tria)

                     DO l = 0, lmaxU_const
                        DO lp = 0, lmaxU_const
                           cil = ImagUnit**(l-lp)
                           DO m = -l,l
                              lm = l*(l+1)+m
                              DO mp = -lp,lp
                                 lmp = lp*(lp+1)+mp
                                 tmp = cil * phase * conjg(eigVecCoeffs%acof(ib,lm,n,ispin))*eigVecCoeffs%acof(ib,lmp,np,ispin)
                                 greensfCoeffs%uu_int(ie,n,np,lm,lmp,ispin) = greensfCoeffs%uu_int(ie,n,np,m,mp,ispin) - pi_const * weight * REAL(tmp)
                                 
                                 tmp = cil * phase * conjg(eigVecCoeffs%bcof(ib,lm,n,ispin))*eigVecCoeffs%bcof(ib,lmp,np,ispin)                                                                    
                                 greensfCoeffs%dd_int(ie,n,np,lm,lmp,ispin) = greensfCoeffs%dd_int(ie,n,np,m,mp,ispin) - pi_const * weight * REAL(tmp)
                                 
                                 tmp = cil * phase * conjg(eigVecCoeffs%acof(ib,lm,n,ispin))*eigVecCoeffs%bcof(ib,lmp,np,ispin)                                                                
                                 greensfCoeffs%du_int(ie,n,np,lm,lmp,ispin) = greensfCoeffs%du_int(ie,n,np,m,mp,ispin) - pi_const * weight * REAL(tmp)
                                                                                                
                                 tmp = cil * phase * conjg(eigVecCoeffs%bcof(ib,lm,n,ispin))*eigVecCoeffs%acof(ib,lmp,np,ispin)
                                 greensfCoeffs%ud_int(ie,n,np,lm,lmp,ispin) = greensfCoeffs%ud_int(ie,n,np,m,mp,ispin) - pi_const * weight * REAL(tmp)
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE intersite_coeffs
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