!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_ferhis
CONTAINS
  SUBROUTINE ferhis(input,kpts,fmpi, index,idxeig,idxkpt,idxjsp,nspins,n,&
                    nstef,ws,spindg,weight, e,ne,we, noco,cell,ef,seigv,w_iks,results,l_output)
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
    USE m_types
    USE m_constants
    USE m_efnewton
    USE m_xmlOutput

    IMPLICIT NONE

    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_mpi),INTENT(IN)          :: fmpi
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_kpts),INTENT(IN)         :: kpts
    TYPE(t_noco),INTENT(IN),OPTIONAL         :: noco
    TYPE(t_cell),INTENT(IN),OPTIONAL         :: cell
    !     ..
    !     .. Scalar Arguments ..
    INTEGER,INTENT(IN)  ::  nspins,n ,nstef
    REAL,INTENT(IN)     ::  spindg,ws,weight
    REAL,INTENT(INOUT)  ::  ef,seigv
    REAL,INTENT(OUT)    ::  w_iks(:,:,:)
    !     .. Scalar Arguments ..
    LOGICAL, INTENT(IN) :: l_output
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: idxeig(:)!(input%neig*kpts%nkpt*dimension%jspd)
    INTEGER, INTENT (IN) :: idxjsp(:)!(input%neig*kpts%nkpt*dimension%jspd)
    INTEGER, INTENT (IN) :: idxkpt(:)!(input%neig*kpts%nkpt*dimension%jspd)
    INTEGER, INTENT (IN) ::  INDEX(:)!(input%neig*kpts%nkpt*dimension%jspd)
    INTEGER, INTENT (IN) ::     ne(:,:)!(kpts%nkpt,dimension%jspd)
    REAL,    INTENT (IN) ::      e(:)!(kpts%nkpt*input%neig*dimension%jspd)
    REAL,    INTENT (INOUT) ::  we(:)!(kpts%nkpt*input%neig*dimension%jspd)

    !--- J constants
    !--- J constants

    !     ..
    !     .. Local Scalars ..
    REAL,PARAMETER:: del=1.e-6
    REAL :: efermi,emax,emin,entropy,fermikn,gap,&
              wfermi,wvals,w_below_emin,w_near_ef,tkb, seigvTemp
    INTEGER ink,inkem,j,js,k,kpt,nocc,nocst,i

    !     .. Local Arrays ..      
    REAL :: qc(3)
    CHARACTER(LEN=20)    :: attributes(2)

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

    tkb=input%tkb !might be modified if we have an insulator
    IF ( fmpi%irank == 0 .and. l_output) THEN
       WRITE (oUnit,FMT='(/)')
       WRITE (oUnit,FMT='(''FERHIS:  Fermi-Energy by histogram:'')')
    END IF

    efermi = ef
    IF (nstef.LT.n) THEN
       gap = e(INDEX(nstef+1)) - ef
       results%bandgap = gap*hartree_to_ev_const
       IF ( fmpi%irank == 0 .and. l_output) THEN
          attributes = ''
          WRITE(attributes(1),'(f20.10)') gap*hartree_to_ev_const
          WRITE(attributes(2),'(a)') 'eV'
          CALL writeXMLElement('bandgap',(/'value','units'/),attributes)
          WRITE (oUnit,FMT=8050) gap
       END IF
    END IF
    IF ( fmpi%irank == 0 .and. l_output) THEN
       WRITE (oUnit,FMT=8010) spindg* (ws-weight)
    END IF
    !
    !---> DETERMINE OCCUPATION AT THE FERMI LEVEL
    !
    wfermi = ws - weight
    !                                          -6
    !======> DETERMINE FERMI ENERGY for kT >= 10
    !
    !
    IF ((tkb.GE.del).AND.(nstef.NE.0)) THEN
       !
       !---> TEMPERATURE BROADENING
       !
       IF (nstef.LT.n) THEN
          !
          !--->    STATES ABOVE EF AVAILABLE           
          !
          ef = 0.5* (e(INDEX(nstef+1))+ef)
          emax = ef + 8.0*tkb
          emin = ef - 8.0*tkb
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
             IF ( fmpi%irank == 0 .and. l_output) THEN
                WRITE (oUnit,*) 'CAUTION!!!  All calculated eigenvalues ', 'are below ef + 8kt.'
             END IF
          ENDIF

          w_near_ef = weight - w_below_emin

          IF (w_near_ef.GT.del) THEN
             !
             !--->       STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE:
             !--->            ADJUST FERMI-ENERGY BY NEWTON-METHOD
             !
             nocst = ink - 1
             CALL ef_newton(n,fmpi%irank, inkem,nocst,index,tkb,e, w_near_ef,ef,we)
             !
             IF ( fmpi%irank == 0 .and. l_output) THEN
                WRITE (oUnit,FMT=8030) ef,spindg*weight, spindg*w_below_emin,spindg* (w_below_emin+w_near_ef)
             END IF

          ELSE
             !
             !--->       NO STATES BETWEEN EF-8kt AND EF+8kt AVAILABLE
             !
             IF ( fmpi%irank == 0 .and. l_output) WRITE (oUnit,FMT=8020)
             nocst = nstef
             we(INDEX(nocst)) = we(INDEX(nocst)) - wfermi
             ef = efermi
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

    ELSE IF (nstef.NE.0) THEN
       !
       !---> NO TEMPERATURE BROADENING IF tkb < del
       !
       nocst = nstef
       we(INDEX(nocst)) = we(INDEX(nocst)) - wfermi
    ELSE
       ! zero occupation
       nocst = nstef
    END IF
    !
    !      write(oUnit,*) nocst,'    nocst in ferhis'
    !      do  ink = 1,nocst
    !         write(oUnit,*) ink,index(ink),we(index(ink)),
    !     +      '    ink,index(ink),we(index(ink)): weights for eigenvalues'
    !      end do
    !
    !
    !=======>   DETERMINE OCCUPATION NUMBER AND WEIGHT OF EIGENVALUES
    !                     FOR EACH K_POINT
    !
    w_iks(:,:,:) = 0.0

    IF ( fmpi%irank == 0 .and. l_output) WRITE (oUnit,FMT=8080) nocst
    DO i=1,nocst
       w_iks(idxeig(INDEX(i)),idxkpt(INDEX(i)),idxjsp(INDEX(i))) = we(INDEX(i))
    ENDDO
    !
    !======>   CHECK SUM OF VALENCE WEIGHTS
    !

    wvals = 0.0
    DO  js = 1,nspins
       DO  k = 1,kpts%nkpt
          wvals = wvals + SUM(w_iks(:ne(k,js),k,js))
       ENDDO
    ENDDO

    IF ( fmpi%irank == 0 .and. l_output) WRITE (oUnit,FMT=8070) wvals
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
             fermikn = w_iks(nocc,kpt,js)/kpts%wtkpt(kpt)
             IF ( fermikn .GT. 0.0 .AND. fermikn .LT. 1.0 ) &
                  entropy = entropy + kpts%wtkpt(kpt) * ( fermikn * LOG( fermikn) + ( 1.0 - fermikn) * LOG( 1.0 - fermikn) )
          END DO
       END DO
    ENDDO
    entropy = -spindg*entropy
    results%ts = tkb*entropy
    results%tkb_loc = tkb
    IF ( fmpi%irank == 0 .and. l_output) WRITE (oUnit,FMT=8060) entropy,entropy*3.0553e-6 !: boltzmann constant in htr/k



    !
    !=======>   DETERMINE SINGLE PARTICLE ENERGY
    !
    !

    seigv = seigv+spindg*DOT_PRODUCT(e(INDEX(:nocst)),we(INDEX(:nocst)))
    seigvTemp = seigv
    IF (noco%l_soc .AND. (.NOT. noco%l_noco)) THEN
       seigvTemp = seigvTemp / 2.0
    END IF
    IF (fmpi%irank == 0 .and. l_output) THEN
       attributes = ''
       WRITE(attributes(1),'(f20.10)') seigvTemp
       WRITE(attributes(2),'(a)') 'Htr'
       CALL writeXMLElement('sumValenceSingleParticleEnergies',(/'value','units'/),attributes)
       WRITE (oUnit,FMT=8040) seigvTemp
    END IF

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
8070 FORMAT (10x,'sum of the valence weights  :',f12.6)
8080 FORMAT (10x,'number of occ. states       :',i10)

  END SUBROUTINE ferhis
END MODULE m_ferhis
