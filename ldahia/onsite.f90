MODULE m_onsite

!------------------------------------------------------------------------------
!
! MODULE: m_onsite
!
!> @author
!> Henning Janßen
!
! DESCRIPTION: 
!>  This module contains the functions to calculate the imaginary part of the 
!>  onsite GF with and without radial dependence
!>  Further we can transform this imaginary part to obtain the onsite GF
!>  using the Kramer Kronig Transformation
!
!------------------------------------------------------------------------------

USE m_juDFT
USE m_types
USE m_constants

CONTAINS

SUBROUTINE onsite_coeffs(atoms,input,ispin,nbands,tetweights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)

   !This Subroutine calculates the contribution to the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !of the current k-Point (it is called in cdnval) inside the MT-sphere 
   !and sums over the Brillouin-Zone using the histogram method or linear tetrahedron method
   !It is essentially the l-density of states in a (m,mp) matrix with an additional factor - pi

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
   TYPE(t_usdus),         INTENT(IN)    :: usdus
   TYPE(t_greensfCoeffs), INTENT(INOUT) :: greensfCoeffs
   TYPE(t_input),         INTENT(IN)    :: input 

   !-Scalar Arguments 
   INTEGER,               INTENT(IN)    :: ispin  !Current spin index
   INTEGER,               INTENT(IN)    :: nbands !Number of bands to be considered
   REAL,                  INTENT(IN)    :: wtkpt  !Weight of the current k-point

   !-Array Arguments
   REAL,                  INTENT(IN)    :: tetweights(greensfCoeffs%ne,nbands) !Precalculated tetrahedron weights for the current k-point
   INTEGER,               INTENT(IN)    :: ind(nbands,2)                       !Gives the range where the tetrahedron weights are non-zero
   REAL,                  INTENT(IN)    :: eig(nbands)                         !Eigenvalues for the current k-point

   !-Local Scalars
   LOGICAL l_zero
   INTEGER i_gf,ib,ie,j,nType,natom,l,m,mp,lm,lmp,ilo,ilop
   REAL    fac,weight



   !Loop through the gf elements to be calculated
   DO i_gf = 1, atoms%n_gf

      l     = atoms%onsiteGF(i_gf)%l
      nType = atoms%onsiteGF(i_gf)%atomType
      fac = 1.0  / atoms%neq(nType)

      !Loop through equivalent atoms
      DO natom = SUM(atoms%neq(:nType-1)), SUM(atoms%neq(:nType))

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(natom,l,nType,ispin,wtkpt,nbands,i_gf,fac,input) &
         !$OMP SHARED(atoms,eigVecCoeffs,usdus,greensfCoeffs,eig,tetweights,ind) &
         !$OMP PRIVATE(ib,ie,j,m,mp,lm,lmp,ilo,ilop,l_zero,weight)

         !$OMP DO
         !Loop through bands
         DO ib = 1, nbands

            !Check wether there is a non-zero weight for the energy window
            l_zero = .true.
            IF(input%tria) THEN
               !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
               IF(ANY(tetweights(ind(ib,1):ind(ib,2),ib).NE.0.0)) l_zero = .false.
            ELSE
               !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
               j = NINT((eig(ib)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
               IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) )         l_zero = .false.
            END IF

            IF(l_zero) CYCLE

            !Choose the relevant energy points depending on the bz-integration method
            DO ie = MERGE(ind(ib,1),j,input%tria), MERGE(ind(ib,2),j,input%tria)
               !weight for the bz-integration including spin-degeneracy
               weight = 2.0/input%jspins * fac * MERGE(tetweights(ie,ib),wtkpt/greensfCoeffs%del,input%tria)    
               !
               ! contribution from states
               !    
               DO m = -l, l
                  lm = l*(l+1)+m
                  DO mp = -l,l
                     lmp = l*(l+1)+mp
                        IF(input%onsite_sphavg) THEN
                           !CHANGE: We drop the dependency of u/u_dot on spin and average over the two 
                           greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                                     REAL(conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                                                          conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                                                          SUM(usdus%ddn(l,nType,:)/input%jspins))
                        ELSE
                           greensfCoeffs%uu(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                                  conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                           greensfCoeffs%dd(ie,i_gf,m,mp,ispin) = greensfCoeffs%dd(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                                  conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                           greensfCoeffs%ud(ie,i_gf,m,mp,ispin) = greensfCoeffs%ud(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                                  conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                           greensfCoeffs%du(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                                  conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                        END IF
                  ENDDO
               ENDDO
               !
               ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
               !
               DO ilo = 1, atoms%nlo(nType)
                  IF(atoms%llo(ilo,nType).EQ.l) THEN
                     !There is a local orbital here
                     !We just add it to the projdos in the spherical average case (What to do in the radial dependence case??)
                     DO m = -l, l
                        lm = l*(l+1)+m
                        DO mp = -l, l
                           lmp = l*(l+1)+mp
                           IF(input%onsite_sphavg) THEN
                              greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) &
                                                                        - pi_const *(  usdus%uulon(ilo,nType,ispin) * (&
                                                                        conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                                                        conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) )&
                                                                        + usdus%dulon(ilo,nType,ispin) * (&
                                                                        conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                                                        conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)))

                              DO ilop = 1, atoms%nlo(nType)
                                 IF (atoms%llo(ilop,nType).EQ.l) THEN

                                    greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) &
                                                                              - pi_const * usdus%uloulopn(ilo,ilop,nType,ispin) *&
                                                                              conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,ib,ilo,natom,ispin)

                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO  
            ENDDO !ie
         ENDDO !ib
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO !natom
   ENDDO !i_gf

END SUBROUTINE onsite_coeffs

SUBROUTINE calc_onsite(atoms,enpara,vr,jspins,greensfCoeffs,gOnsite,sym,ef,beta,l_sphavg,onsite_exc_split,hub1)

   USE m_kkintgr
   USE m_radfun
   USE m_gfcalc
   USE m_cfmat
   USE m_hubbard1_io

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_enpara),         INTENT(IN)     :: enpara
   TYPE(t_greensfCoeffs),  INTENT(IN)     :: greensfCoeffs
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_hub1ham),        INTENT(INOUT)  :: hub1

   INTEGER,                INTENT(IN)     :: jspins
   REAL,                   INTENT(IN)     :: ef
   REAL,                   INTENT(IN)     :: beta
   REAL,                   INTENT(IN)     :: vr(atoms%jmtd,atoms%ntype,jspins)
   REAL,                   INTENT(OUT)     :: onsite_exc_split

   LOGICAL,                INTENT(IN)     :: l_sphavg

   TYPE(t_usdus) usdus
   INTEGER i_gf,i,l,m,mp,jr,noded,nodeu,n,j,jspin,i_hia,it,is,isi,natom,n_op
   INTEGER e_cut(2)
   REAL wronk,fac
   CHARACTER(len=30) :: filename

   REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)
   COMPLEX mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_gf,jspins), sym_diag(-lmaxU_const:lmaxU_const)
   REAL,    ALLOCATABLE :: im(:,:,:,:,:)
   REAL, ALLOCATABLE :: e(:)
   COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
   COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspins) )
   ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jspins) )

   CALL usdus%init(atoms,jspins)

   IF(greensfCoeffs%sigma.NE.0.0) THEN
      !construct an energy grid for smoothing (required by the function in m_smooth)
      ALLOCATE (e(greensfCoeffs%ne))
      DO i = 1, greensfCoeffs%ne
         e(i) = greensfCoeffs%del * (i-1) + greensfCoeffs%e_bot
      ENDDO
   END IF

   ALLOCATE( im(MAXVAL(gOnsite%nr(:)),greensfCoeffs%ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins) )

   DO i_gf = 1, atoms%n_gf
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType

      im = 0.0
      DO jspin = 1, jspins
         !The functions f and g can probably be taken from another call of the routine
         IF(gOnsite%nr(i_gf).NE.1) THEN
            CALL radfun(l,n,jspin,enpara%el0(l,n,jspin),vr(:,n,jspin),atoms,&
                           f(:,:,l,jspin),g(:,:,l,jspin),usdus, nodeu,noded,wronk)
         ENDIF

         DO m = -l, l
            DO mp = -l,l
               !calculate the radial dependence
               DO jr = 1, gOnsite%nr(i_gf)
                  IF(gOnsite%nr(i_gf).NE.1) THEN
                     im(jr,:,m,mp,jspin) = im(jr,:,m,mp,jspin) + &
                                          greensfCoeffs%uu(:,i_gf,m,mp,jspin) * (f(jr,1,l,jspin)*f(jr,1,l,jspin)+f(jr,2,l,jspin)*f(jr,2,l,jspin)) +&
                                          greensfCoeffs%dd(:,i_gf,m,mp,jspin) * (g(jr,1,l,jspin)*g(jr,1,l,jspin)+g(jr,2,l,jspin)*g(jr,2,l,jspin)) +&
                                          greensfCoeffs%du(:,i_gf,m,mp,jspin) * (f(jr,1,l,jspin)*g(jr,1,l,jspin)+f(jr,2,l,jspin)*g(jr,2,l,jspin))
                  ELSE
                     im(1,:,m,mp,jspin) = greensfCoeffs%projdos(:,i_gf,m,mp,jspin)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      !Enforcing that the projected density of states follows the local symmetries
      !
     
      DO jr = 1, gOnsite%nr(i_gf)
         DO j = 1, greensfCoeffs%ne
            DO jspin = 1, jspins
               n_tmp = 0.0
               n_tmp(-l:l,-l:l) = im(jr,j,-l:l,-l:l,jspin)
               im(jr,j,-l:l,-l:l,jspin) = 0.0
               n_op = 0
               DO natom = SUM(atoms%neq(:n-1)) + 1, SUM(atoms%neq(:n-1)) + atoms%neq(n)
                  DO it = 1, sym%invarind(natom)
                     is = sym%invarop(natom,it)
                     isi = sym%invtab(is)
                     d_tmp(:,:) = cmplx(0.0,0.0)
                     DO m = -l,l
                        DO mp = -l,l
                           d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                        ENDDO
                     ENDDO
                     DO m = -l,l
                        sym_diag(m) = d_tmp(m,m)
                        d_tmp(m,m) = 0.0
                     ENDDO
                     !Exclude all symmetries that would prevent splitting of the levels
                     IF(ANY(d_tmp(:,:).NE.0.0)) CYCLE
                     n_op = n_op + 1
                     DO m = -l,l
                        d_tmp(m,m) = sym_diag(m)
                     ENDDO
                     nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                     n1_tmp =  matmul( nr_tmp, d_tmp )
                     DO m = -l,l
                        DO mp = -l,l
                           im(jr,j,m,mp,jspin) = im(jr,j,m,mp,jspin) + conjg(n1_tmp(m,mp))
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               IF(n_op.EQ.0) CALL juDFT_error("Symmetry error",calledby="calc_onsite")
               fac = 1.0/(REAL(n_op))
               im(jr,j,:,:,jspin) = im(jr,j,:,:,jspin) * fac
            ENDDO
         ENDDO            
      ENDDO
      !
      !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 
      !
      CALL greensf_cutoff(im,atoms,gOnsite%nr(i_gf),l,n,jspins,greensfCoeffs%ne,greensfCoeffs%del,greensfCoeffs%e_bot,greensfCoeffs%e_top,l_sphavg,ef,onsite_exc_split,e_cut)

      !DO i_hia = 1, atoms%n_hia
      !   IF(hub1%lda_u(i_hia)%l.NE.l.OR.hub1%lda_u(i_hia)%atomType.NE.n) CYCLE
      !   CALL cfcontrib(-1/pi_const * im(1,:,:,:,:),l,n,e_cut(1),e_cut(2),greensfCoeffs%ne,greensfCoeffs%del,jspins,hub1%ccfmat(i_hia,:,:)) 
      !ENDDO
      

      CALL timestart("On-Site: Kramer-Kronigs-Integration")
      DO jspin = 1, jspins
         DO m= -l,l
            DO mp= -l,l
               DO jr = 1, gOnsite%nr(i_gf)
                  !G(E+i delta)
                  CALL kkintgr(im(jr,:,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                              gOnsite%gmmpMat(jr,:,i_gf,m,mp,jspin,1),gOnsite%e,.false.,gOnsite%mode,gOnsite%nz,1)
                  !G(E-i delta)
                  CALL kkintgr(im(jr,:,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                              gOnsite%gmmpMat(jr,:,i_gf,m,mp,jspin,2),gOnsite%e,.true.,gOnsite%mode,gOnsite%nz,1)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL timestop("On-Site: Kramer-Kronigs-Integration")
      !CALL occmtx(gOnsite,i_gf,atoms,sym,jspins,beta*hartree_to_ev_const,ef,mmpMat(:,:,i_gf,:))
      !CALL write_gf("gdft",gOnsite,i_gf)
      !CALL ldosmtx("g",gOnsite,i_gf,atoms,sym,jspins)
   ENDDO

END SUBROUTINE calc_onsite


SUBROUTINE greensf_cutoff(im,atoms,nr,l,n,jspins,ne,del,e_bot,e_top,l_sphavg,ef,onsite_exc_split,e_cut)
   !This Subroutine determines the cutoff energy for the kramers-kronig-integration
   !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
   !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

   USE m_intgr
   USE m_kkintgr
   
   IMPLICIT NONE

   REAL,                INTENT(INOUT)  :: im(nr,ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
   TYPE(t_atoms),       INTENT(IN)     :: atoms

   INTEGER,             INTENT(IN)     :: nr
   INTEGER,             INTENT(IN)     :: l
   INTEGER,             INTENT(IN)     :: n
   INTEGER,             INTENT(IN)     :: jspins
   INTEGER,             INTENT(IN)     :: ne

   REAL,                INTENT(IN)     :: del
   REAL,                INTENT(IN)     :: e_bot
   REAL,                INTENT(IN)     :: e_top
   REAL,                INTENT(IN)     :: ef
   REAL,                INTENT(OUT)    :: onsite_exc_split

   INTEGER,             INTENT(OUT)    :: e_cut(2)

   LOGICAL,             INTENT(IN)     :: l_sphavg


   REAL, ALLOCATABLE :: fDOS(:,:)
   REAL :: projDOS(ne+1,jspins,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   REAL, ALLOCATABLE :: tmp(:,:)

   INTEGER i,m,mp,j,n_c,kkintgr_cut,jr,ispin,cut

   REAL integral, m_up, m_dwn,norm_up,norm_dwn
   REAL a,b, imag, n_states
   LOGICAL l_write,l_checkprojDOS

   ALLOCATE(fDOS(ne,jspins))
   ALLOCATE(tmp(ne,jspins))

   l_write=.false.
   l_checkprojDOS = .false.

   !For debugging: print out the proj. density of states integrated over energy
   IF (l_checkprojDOS) THEN
      projDOS = 0.0
      DO ispin = 1, jspins
         DO m = -l , l
            DO mp = -l , l
               DO j = 1, INT((ef-e_bot)/del) + 1
                  IF(l_sphavg) THEN
                     imag = im(1,j,m,mp,ispin)
                  ELSE
                     CALL intgr3(im(:,j,m,mp,ispin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)
                  END IF
                  projDOS(j,ispin,m,mp) = projDOS(j,ispin,m,mp) + imag
                  IF(m.EQ.mp.AND.imag.GT.0) CALL juDFT_error("SIGN ERROR")
               ENDDO
               CALL trapz(projDOS(1:ne,ispin,m,mp), del, ne, projDOS(ne+1,ispin,m,mp))
            ENDDO
         ENDDO
      ENDDO
      WRITE(*,"(7f14.8)") -1/pi_const * projDOS(ne+1,1,:,:)
      WRITE(*,"(7f14.8)") -1/pi_const * projDOS(ne+1,2,:,:)
   ENDIF


   fDOS = 0.0

   !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS 

   !n_f(e) = -1/pi * TR[Im(G_f(e))]
   DO ispin = 1, jspins
      DO m = -l , l
         DO j = 1, ne
            IF(l_sphavg) THEN
               imag = im(1,j,m,m,ispin)
            ELSE
               CALL intgr3(im(:,j,m,m,ispin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)
            END IF
            fDOS(j,ispin) = fDOS(j,ispin) + imag
         ENDDO
      ENDDO
   ENDDO
   fDOS(:,:) = -1/pi_const*fDOS(:,:)


   IF(l_write) THEN
      
      OPEN(1337,file="fDOS_up.txt",action="write",status="replace")

      DO i = INT(ne/2.0), ne
         WRITE(1337,*) ((i-1)*del + e_bot-ef)*hartree_to_ev_const, fDOS(i,1)/hartree_to_ev_const
      ENDDO

      CLOSE(unit = 1337)

      IF(jspins.EQ.2) THEN
         OPEN(1337,file="fDOS_dwn.txt",action="write",status="replace")

         DO i = INT(ne/2.0), ne
            WRITE(1337,*) ((i-1)*del + e_bot-ef)*hartree_to_ev_const, -fDOS(i,2)/hartree_to_ev_const
         ENDDO

         CLOSE(unit = 1337)
      ENDIF
   ENDIF
   
   IF(jspins.EQ.2) fDOS(:,1) = fDOS(:,1) + fDOS(:,2)

   CALL trapz(fDOS(1:ne,1), del, ne, integral)

   n_states = 2*(2*l+1)
   
   !WRITE(*,*) "Integral over fDOS: ", integral

   kkintgr_cut = ne

   IF(integral.LT.n_states-0.5) THEN
      ! If the integral is to small we stop here to avoid problems
      CALL juDFT_error("fDOS-integral too small: make sure numbands is big enough", calledby="greensf_cutoff")
      
   ELSE IF((integral.GT.n_states).AND.((integral-n_states).GT.0.001)) THEN
      !IF the integral is bigger than 14, search for the cutoff using the bisection method   

      a = e_bot
      b = e_top

      DO

         n_c = INT(((a+b)/2.0-e_bot)/del)+1
         CALL trapz(fDOS(1:n_c,1),del,n_c,integral)

         IF((ABS(integral-n_states).LT.0.001).OR.(ABS(a-b)/2.0.LT.del)) THEN

            kkintgr_cut = INT(((a+b)/2.0-e_bot)/del)+1
            EXIT

         ELSE IF((integral-n_states).LT.0) THEN
            a = (a+b)/2.0
         ELSE IF((integral-n_states).GT.0) THEN
            b = (a+b)/2.0
         END IF

      ENDDO

      CALL trapz(fDOS(1:kkintgr_cut,1),del,kkintgr_cut,integral)

      !WRITE(*,*) "CALCULATED CUTOFF: ", kkintgr_cut
      !WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
   END IF

   e_cut(1) = 1
   e_cut(2) = kkintgr_cut
   !If we are in the magnetic case we want to calculate the effectiv exchange interaction j0 from the Gf
   !For that we can use the difference in the center of gravity of the up and down bands
   !This is calculated here
   cut = 2500
   IF(jspins.EQ.2) THEN
      fDOS(:,1) = fDOS(:,1)-fDOS(:,2)
      !multiply fDOS(E)*E
      tmp = 0.0
      DO j = 1, cut
         tmp(j,:) = fDOS(j,:) *((j-1) * del + e_bot)
      ENDDO
      CALL trapz(tmp(1:cut,1),del,cut,m_up)
      CALL trapz(tmp(1:cut,2),del,cut,m_dwn)
      CALL trapz(fDOS(1:cut,1),del,cut,norm_up)
      CALL trapz(fDOS(1:cut,2),del,cut,norm_dwn)
      onsite_exc_split = m_dwn/norm_dwn - m_up/norm_up
   ENDIF

   !Now we set the imaginary part of the greens function to zero above this cutoff
   DO ispin = 1, jspins
      DO i = kkintgr_cut+1, ne
         DO m = -l, l
            DO mp = -l, l
               IF(l_sphavg) THEN
                  im(1,i,m,mp,ispin) = 0.0
               ELSE
                  DO jr = 1, nr
                     im(jr,i,m,mp,ispin) = 0.0
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE greensf_cutoff

END MODULE m_onsite