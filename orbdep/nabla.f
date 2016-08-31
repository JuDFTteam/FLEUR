      MODULE m_nabla
      use m_juDFT
        
      CONTAINS

      SUBROUTINE nabla(
     >                 ispecies,number_of_j1,grid_size,delta_x,
     >                 nstd,ntypd,j1,l1,lmax,ms,ri,psi,phi,dphi,
     <                 psi_phi)
!----------------------------------------------------------------
!
!     ispecies     ... number of the atom (itype)
!     number_of_j1 ... number of core wavefunction
!     grid_size    ... rumber of radial gridpoints (jri)
!     delta_x      ... logarithmic increment of grid (dx)
!     ri(grid_size)... radial mesh
!     nstd     ... number of corelevels (dimension)
!     ntypd    ... number of atoms (dimension)
!     j1,l1    ... quantum numbers of the core wavefunction
!     ms       ... -1/2 or +1/2 for spin 1 or 2
!     lmax     ... parameter, 3 (s,p,d,f)
!     psi(r)   ... core wavefunction
!     phi(r,l) ... valence wavefunction
!    dphi(r,l) ... radial derivative of valence wavefunction
!     
!----------------------------------------------------------------
       USE m_clebsch
       USE m_intgr, ONLY : intgr3
       IMPLICIT NONE  

       INTEGER, INTENT(IN) :: ispecies, number_of_j1, grid_size
       INTEGER, INTENT(IN) :: l1, lmax, nstd, ntypd
       REAL,    INTENT(IN) :: delta_x, j1, ms
       REAL,    INTENT(IN) :: psi(grid_size), ri(grid_size)
       REAL,    INTENT(IN) :: phi(grid_size,0:lmax)
       REAL,    INTENT(IN) ::dphi(grid_size,0:lmax)
       COMPLEX, INTENT(OUT):: psi_phi(nstd,(lmax+1)**2,3*ntypd)

       INTEGER :: m1, l2, m2, index, alloc_error, lmn1, lmn2
       REAL  :: result, result1, total_result, spin
       REAL  :: mu, cnst_one_over_sqrt_two, cnst_zero
       COMPLEX :: cnst_i
       REAL, DIMENSION(:), POINTER :: f
       spin = 0.50
       cnst_one_over_sqrt_two = 1.0/sqrt(2.0)
       cnst_i = cmplx(0.0,1.0)
       cnst_zero = 0.0

       NULLIFY(f)

       IF ( ASSOCIATED(f) ) THEN
        WRITE(6,*)'nabla: f association status:',ASSOCIATED(f)
        STOP
       ENDIF

       ALLOCATE ( f(grid_size),STAT=alloc_error )
       IF (alloc_error /= 0)  CALL juDFT_error("Couldn't allocate f",
     +      calledby ="nabla")
      
       lmn1 = 2 * (number_of_j1 - 1)  * l1
       mu = -j1

       DO WHILE (mu <= j1)
        lmn1 = lmn1 + 1
        m1 = INT(mu - ms) 
        lmn2 = 0
        DO l2 = 0, lmax
         DO m2 = -l2, l2
            lmn2 = lmn2 + 1
            IF(l1 == l2 + 1)THEN  
             total_result = 0.00
             result  = 0.00
             result1 = 0.00
!     
! (l+1)/srqt[(2l+1)(2l+3)] < phi_core | (d phi_valence / dr ) >
!
             f(:) = psi(:) * dphi(:,l2) ! assumed to be already multiplied with * ri(:) * ri(:)
      
             CALL intgr3(f,ri,delta_x,grid_size,result)
      
             result = result * (l2 + 1.00) /
     +                         sqrt((2.00*l2 +1.00)*(2.*l2+3.00))
      
!
!  - l(l+1)/srqt[(2l+1)(2l+3)] < phi_core | (1/r) phi_valence >
!
             f(:) = psi(:) * phi(:,l2) ! assumed to be already multiplied with 1G* ri(:)
             CALL intgr3(f,ri,delta_x,grid_size,result1)
             result1 = - result1 * l2 * (l2 + 1.0) /              
     +                   sqrt((2.0 * l2 + 1.00) * (2. * l2 + 3.00)) 
!
! Sum up and decorate with Clebsch-Gordon coefficients
!
             result  = result + result1
             total_result = result*clebsch(real(l1),spin,mu-ms,ms,j1,mu)
      
             index = (ispecies - 1) * 3 + 1
             psi_phi(lmn1,lmn2,index)= cgc(l2,1,l1,m2,1,m1) *       ! left polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
      
             index = index + 1
             psi_phi(lmn1,lmn2,index)= cgc(l2,1,l1,m2,-1,m1) *      ! right polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
      
             index = index + 1
             psi_phi(lmn1,lmn2,index)= cgc(l2,1,l1,m2,0,m1) *        ! z-polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
      
            ELSEIF(l1== l2-1)THEN  
      
              result  =  cnst_zero 
              result1 =  cnst_zero
!
! l/srqt[(2l-1)(2l+1)] < phi_core | (d phi_valence / dr ) >
!
              f(:) = psi(:)* dphi(:,l2)   * ri(:) * ri(:)
              CALL intgr3(f,ri,delta_x,grid_size,result)
              result = result * l2 / sqrt((2.0*l2 - 1.0)*(2.0*l2 + 1.0)) 
!
!   l(l+1)/srqt[(2l-1)(2l+1)] < phi_core | (1/r) phi_valence >
!
              f(:) = psi(:)* phi(:,l2) * ri(:)
              CALL intgr3(f,ri,delta_x,grid_size,result1)
              result1 = result1 * l2 * (l2 + 1.0) /
     +                  sqrt((2.00 * l2 - 1.00) * (2.00 * l2 + 1.00)) 
!
! Sum up and decorate with Clebsch-Gordon coefficients
!
              result = result + result1
              total_result= result*clebsch(real(l1),spin,mu-ms,ms,j1,mu)
      
      
              index = (ispecies - 1) * 3 + 1
              psi_phi(lmn1,lmn2,index) = cgc(l2,1,l1,m2,1,m1) *          ! left polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
              index = index + 1
              psi_phi(lmn1,lmn2,index) = cgc(l2,1,l1,m2,-1,m1) *         ! right polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
              index = index + 1
              psi_phi(lmn1,lmn2,index) = cgc(l2,1,l1,m2,0,m1) *          ! z-polarization
     +                                 total_result / cgc(l2,1,l1,0,0,0)
            ENDIF         
          ENDDO
        ENDDO
        mu = mu + 1.00
       ENDDO
       DEALLOCATE(f,STAT=alloc_error)
       IF (alloc_error /= 0)  CALL juDFT_error("Couldn't deallocate f"
     +      ,calledby ="nabla")
      END SUBROUTINE nabla

      FUNCTION cgc(l1,l2,l3,m1,m2,m3)

      IMPLICIT NONE  
      INTEGER :: l1, l2, l3, m1, m2, m3
      REAL  :: two_l1p1, two_l1p2, l1pm3, l1pm3p1, l1mm3p1, l1mm3, cgc 

      IF (m3 /= m1 + m2) THEN
       cgc = 0.0
       RETURN
      END IF 
!     gb  m3 = m1 + m2
      two_l1p1 = 2 * l1 + 1
      two_l1p2 = 2 * l1 + 2
      l1pm3 = l1 + m3
      l1pm3p1 = l1 + m3 + 1
      l1mm3p1 = l1 - m3 + 1
      l1mm3 = l1 - m3 
      cgc = 0.0 
      IF (l3 == l1 + 1) THEN
          IF (m2 == 1) then
           cgc = sqrt( (l1pm3 * l1pm3p1) / (two_l1p1 * two_l1p2))   
          ELSEIF (m2 == 0) THEN
           cgc = sqrt( (l1mm3p1 * l1pm3p1) / (two_l1p1 * (l1 + 1)))
          ELSEIF (m2 == -1) THEN
           cgc = sqrt( (l1mm3 * l1mm3p1) / (two_l1p1 * two_l1p2))   
          END IF
      ELSE IF(l3 == l1 -1) THEN
          IF (m2 == 1) then
           cgc = sqrt( (l1mm3 * l1mm3p1) / (2.d0 * l1 * two_l1p1))   
          ELSEIF (m2 == 0) THEN
           cgc = -sqrt( (l1mm3 * l1pm3) / (l1 * (two_l1p1)))   
          ELSEIF (m2 == -1) THEN
           cgc = sqrt( (l1pm3p1 * l1pm3) / (2.0d0 * l1 * two_l1p1))   
          END IF
      END IF
      END FUNCTION cgc


      END MODULE m_nabla
