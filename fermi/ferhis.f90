MODULE m_ferhis
CONTAINS
  SUBROUTINE ferhis(input,kpts,mpi,results, index,idxeig,idxkpt,idxjsp,n,&
       nstef,ws,spindg,weight, e,ne,we, noco,jij,cell)
    !***********************************************************************
    !
    !     This subroutine determines the fermi energy and the sum of the
    !     single particle eigenvalues by histogram method.
    !
    !
    !     Theory :   zelec(nwd) = sum{ sum{ sum{ we * f(e) } } }
    !                             sp    e    k
    !
    !
    !                seigv = sum{ sum{ sum{ w(k) * f(e) * e }
    !                         sp   e    k
    !
    !                the weights w(k) are normalized: sum{w(k)} = 1
    !                                                  k                -6
    !                         a) 1                           for kT < 10
    !                we    = {                           -1             -6
    !                         b){ 1+exp(e(k,nu) -ef)/kt) }   for kt >=10
    !
    !                in case a) we choose the Fermi energy the highest
    !                           valence state
    !                        b) we choose as Fermi energy the arithmetric
    !                           mean between the highest occupied and lowest
    !                           unoccupied state if this energy difference
    !                           Delta E <= kT, otherwise as a).
    !
    !                                      stefan bl"ugel , kfa , oct 1991
    !
    !               free energy and extrapolation T -> 0  implemented
    !                         (see M.J.Gillan, J.Phys.: Condens. Matter 1,
    !                          (1989) 689-711 )
    !
    !                                      peter richard, kfa, jun. 1995
    !
    !               adapted to flapw7
    !
    !                                      philipp kurz, kfa, oct. 1995
    !               entropy calculation changed
    !                     
    !                                      r.pentcheva, kfa, may  1996
    !
    !***********************************************************************
    USE m_efnewton
    USE m_types
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_mpi),INTENT(IN)          :: mpi
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_kpts),INTENT(IN)         :: kpts
    TYPE(t_noco),INTENT(IN),OPTIONAL         :: noco
    TYPE(t_jij),INTENT(IN),OPTIONAL          :: jij
    TYPE(t_cell),INTENT(IN),OPTIONAL         :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)  ::  n ,nstef
    REAL,INTENT(IN)     ::  spindg,ws,weight
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: idxeig(:)!(dimension%neigd*kpts%nkptd*dimension%jspd)
    INTEGER, INTENT (IN) :: idxjsp(:)!(dimension%neigd*kpts%nkptd*dimension%jspd)
    INTEGER, INTENT (IN) :: idxkpt(:)!(dimension%neigd*kpts%nkptd*dimension%jspd)
    INTEGER, INTENT (IN) ::  INDEX(:)!(dimension%neigd*kpts%nkptd*dimension%jspd)
    INTEGER, INTENT (IN) ::     ne(:,:)!(kpts%nkptd,dimension%jspd)
    REAL,    INTENT (IN) ::      e(:)!(kpts%nkptd*dimension%neigd*dimension%jspd)
    REAL,    INTENT (INOUT) ::  we(:)!(kpts%nkptd*dimension%neigd*dimension%jspd)

    !--- J constants
    !--- J constants

    !     ..
    !     .. Local Scalars ..
    REAL,PARAMETER:: del=1.e-6
    REAL :: efermi,emax,emin,entropy,fermikn,gap,&
              wfermi,wvals,w_below_emin,w_near_ef,tkb
    INTEGER ink,inkem,j,js,k,kpt,nocc,nocst,i,nspins

    !     .. Local Arrays ..      
    REAL :: qc(3)

    !     ..
    !***********************************************************************
    !------->          ABBREVIATIONS
    !
    !     eig        : array of eigenvalues within all energy-windows
    !     wtkpt      : list of the weights of each k-point (from inp-file)
    !     e          : linear list of the eigenvalues within the highest
    !                  energy-window
    !     we         : list of weights of the eigenvalues in e
    !     w          : array of weights (output, needed to calculate the
    !                  new charge-density)
    !     zelec      : number of electrons in a window
    !     spindg     : spindegeneracy (2 in nonmagnetic calculations)
    !     seigv      : weighted sum of the occupied valence eigenvalues
    !     seigsc     : weighted sum of the semi-core eigenvalues
    !     seigscv    : sum of seigv and seigsc
    !     ts         : entropy contribution to the free energy
    !     tkb        : value of temperature (kt) broadening around fermi
    !                  energy in htr units
    !     ef         : fermi energy determined to obtain charge neutrality
    !     wfermi     : weight of the occupation number for states at the
    !                  fermi energy.
    !     fd         : fermi dirac distribution
    !     fermikn    : fermi dirac distribution for the k-th point 
    !                  and n-th state
    !**********************************************************************
    !     ..
    nspins=input%jspins
    if (noco%l_noco) nspins=1
    tkb=input%tkb !might be modified if we have an insulator
    IF ( mpi%irank == 0 ) THEN
       WRITE (6,FMT='(/)')
       WRITE (6,FMT='(''FERHIS:  Fermi-Energy by histogram:'')')
    END IF

    efermi = results%ef
    IF (nstef.LT.n) THEN
       gap = e(INDEX(nstef+1)) - results%ef
       IF ( mpi%irank == 0 ) WRITE (6,FMT=8050) gap
    END IF
    IF ( mpi%irank == 0 ) THEN
       WRITE ( 6,FMT=8010) spindg* (ws-weight)
       WRITE (16,FMT=8010) spindg* (ws-weight)
    END IF
    !
    !---> DETERMINE OCCUPATION AT THE FERMI LEVEL
    !
    wfermi = ws - weight
    !                                          -6
    !======> DETERMINE FERMI ENERGY for kT >= 10
    !
    !
    IF (tkb.GE.del) THEN
       !
       !---> TEMPERATURE BROADENING
       !
       IF (nstef.LT.n) THEN
          !
          !--->    STATES ABOVE EF AVAILABLE           
          !
          results%ef = 0.5* (e(INDEX(nstef+1))+results%ef)
          emax = results%ef + 8.0*tkb
          emin = results%ef - 8.0*tkb
          w_near_ef = 0.0
          w_below_emin = 0.0
          inkem = 0
          ink_loop: DO ink = 1,n

             IF (e(INDEX(ink)).LT.emin) THEN
                inkem = ink
                w_below_emin = w_below_emin + we(INDEX(ink))
             ELSE IF (e(INDEX(ink)).GT.emax) THEN
                EXIT ink_loop
             END IF

          ENDDO ink_loop
          IF (ink>n) THEN
             IF ( mpi%irank == 0 ) THEN
                WRITE (6,*) 'CAUTION!!!  All calculated eigenvalues ', 'are below ef + 8kt.'
                WRITE (16,*) 'CAUTION!!!  All calculated eigenvalues ', 'are below ef + 8kt.'
             END IF
          ENDIF

          w_near_ef = weight - w_below_emin

          IF (w_near_ef.GT.del) THEN
             !
             !--->       STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE:
             !--->            ADJUST FERMI-ENERGY BY NEWTON-METHOD
             !
             nocst = ink - 1
             CALL ef_newton(n,mpi%irank, inkem,nocst,index,tkb,e, w_near_ef,results%ef,we)
             !
             IF ( mpi%irank == 0 ) THEN
                WRITE (16,FMT=8030) results%ef,spindg*weight, spindg*w_below_emin,spindg* (w_below_emin+w_near_ef)
                WRITE (6,FMT=8030) results%ef,spindg*weight, spindg*w_below_emin,spindg* (w_below_emin+w_near_ef)
             END IF

          ELSE
             !
             !--->       NO STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE
             !
             IF ( mpi%irank == 0 ) WRITE (6,FMT=8020)
             nocst = nstef
             we(INDEX(nocst)) = we(INDEX(nocst)) - wfermi
             results%ef = efermi
             tkb = 0.0
          END IF
       ELSE
          !
          !--->    NO STATES ABOVE EF AVAILABLE
          !
          tkb = 0.0
          nocst = nstef
          we(INDEX(nocst)) = we(INDEX(nocst)) - wfermi
       END IF

    ELSE
       !
       !---> NO TEMPERATURE BROADENING IF tkb < del
       !
       nocst = nstef
       we(INDEX(nocst)) = we(INDEX(nocst)) - wfermi
    END IF
    !
    !      write(6,*) nocst,'    nocst in ferhis'
    !      do  ink = 1,nocst
    !         write(6,*) ink,index(ink),we(index(ink)),
    !     +      '    ink,index(ink),we(index(ink)): weights for eigenvalues'
    !      end do
    !
    !
    !=======>   DETERMINE OCCUPATION NUMBER AND WEIGHT OF EIGENVALUES
    !                     FOR EACH K_POINT
    !
    results%w_iks(:,:,:) = 0.0

    IF ( mpi%irank == 0 ) WRITE (6,FMT=8080) nocst
    DO i=1,nocst
       results%w_iks(idxeig(INDEX(i)),idxkpt(INDEX(i)),idxjsp(INDEX(i))) = we(INDEX(i))
    ENDDO
    !
    !======>   CHECK SUM OF VALENCE WEIGHTS
    !

    wvals = 0.0
    DO  js = 1,nspins
       DO  k = 1,kpts%nkpt
          wvals = wvals + SUM(results%w_iks(:ne(k,js),k,js))
       ENDDO
    ENDDO

    IF ( mpi%irank == 0 ) WRITE (6,FMT=8070) wvals
    !
    !
    !=======>   DETERMINE ENTROPY
    !
    !
    ! --->   formula for entropy:
    !
    !        entropy = - two * sum wtkpt(kpt) * 
    !                          kpt
    !                       { sum ( f(e(kpt,nu))*log(f(e(kpt,nu)))
    !                          n
    !                              +(1-f(e(kpt,nu)))*log(1-f(e(kpt,nu))) )  }
    !
    !        here we have   w(n,kpt,js)= spindg*wghtkp(kpt)*f(e(kpt,n))
    !
    entropy = 0.0
    DO js = 1,nspins
       DO kpt = 1 , kpts%nkpt
          DO nocc=1,ne(kpt,js) 
             fermikn = results%w_iks(nocc,kpt,js)/kpts%wtkpt(kpt)
             IF ( fermikn .GT. 0.0 .AND. fermikn .LT. 1.0 ) &
                  entropy = entropy + kpts%wtkpt(kpt) * ( fermikn * LOG( fermikn) + ( 1.0 - fermikn) * LOG( 1.0 - fermikn) )
          END DO
       END DO
    ENDDO
    entropy = -spindg*entropy
    results%ts = tkb*entropy
    IF ( mpi%irank == 0 ) WRITE (6,FMT=8060) entropy,entropy*3.0553e-6 !: boltzmann constant in htr/k



    !
    !=======>   DETERMINE SINGLE PARTICLE ENERGY
    !
    !

    results%seigv = spindg*DOT_PRODUCT(e(INDEX(:nocst)),we(INDEX(:nocst)))
    IF ( mpi%irank == 0 ) WRITE (6,FMT=8040) results%seigv

    !--- J constants
    IF (PRESENT(jij)) THEN
       IF (jij%l_J) THEN
          IF (jij%l_disp) THEN
             qc=MATMUL(noco%qss,cell%bmat)
             WRITE (114,FMT=1001) noco%qss(1),noco%qss(2),noco%qss(3),SQRT(dot_product(qc,qc)),results%seigv
          ELSE
             WRITE (114,FMT=1002) noco%qss(1),noco%qss(2),noco%qss(3),results%seigv
          ENDIF
       ENDIF
1001   FORMAT (4(f14.10,1x),f20.10)
1002   FORMAT (3(f14.10,1x),f20.10)
    ENDIF
    !--- J constants 

    !
    ! 7.12.95 r.pentcheva   seigscv = seigsc + seigv   will be
    ! calculated in fermie
    !
8000 FORMAT (/,10x,'==>efrmhi: not enough wavefunctions.',i10,2e20.10)
8010 FORMAT (10x,'charge neutrality (T=0)     :',f11.6,'    (zero if ',&
         &       'the highest occ. eigenvalue is "entirely" occupied)')
8020 FORMAT (/,10x,'no eigenvalues within 8 tkb of ef',&
         &       ' reverts to the t=0 k method.')
8030 FORMAT (/,5x,'-->  new fermi energy            :',f11.6,' htr',&
         &       /,10x,'valence charge              :',f11.6,' e ',/,10x,&
         &       'actual charge blw ef-8kt    :',f11.6,' e ',/,10x,&
         &       'actual charge blw ef+8kt    :',f11.6,' e ')
8040 FORMAT (/,10x,'sum of val. single particle energies: ',f20.10,&
         &       ' htr',/)
8050 FORMAT (/,10x,'bandgap                     :',f11.6,' htr')
8060 FORMAT (10x,'entropy         :',f11.6,' *kb htr/K =',&
         &       f10.5,' htr/K')
8070 FORMAT (10x,'sum of the valence weights  :',f11.6)
8080 FORMAT (10x,'number of occ. states       :',i10)

  END SUBROUTINE ferhis
END MODULE m_ferhis
