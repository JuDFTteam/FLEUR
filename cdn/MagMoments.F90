!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_magmoments
    !! Module to encapsulate routines to calculate and output the spin- and orbital Moments
    contains
    SUBROUTINE spinMoments(input,atoms,noco,nococonv,moments,den)
        !! Calculate and output the spin moments
        !! Either the moments have to be given as calculated in cdnval (moments-argument)
        !! or the moments are calculated from the density 
        USE m_types
        USE m_constants
        USE m_xmlOutput
        USE m_intgr, ONLY : intgr3
        TYPE(t_input), INTENT(IN)           :: input
        TYPE(t_atoms), INTENT(IN)           :: atoms
        TYPE(t_noco), INTENT(IN)            :: noco
        TYPE(t_nococonv), INTENT(IN)        :: nococonv
        TYPE(t_moments),INTENT(IN),OPTIONAL :: moments
        TYPE(t_potden),INTENT(IN),OPTIONAL  :: den
     
        INTEGER                       :: iType
        REAL                          :: up,down,local_m(0:3),global_m(0:3)
        COMPLEX                       :: off_diag
        LOGICAL                       :: l_offdiag
     
        if (input%jspins==1) return ! No spin-pol calculation!
        if (.not.present(moments).and..not.present(den)) return !No data provided
        if (present(moments)) call priv_print_spin_density_at_nucleus(input,atoms,moments)
        

        !Header of output
        write(ounit,*)
        write(oUnit,'(a,19x,a,19x,a,19x,a)') "Spin Magn. mom.  |","Global Frame","  | ","Local Frame"
        write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
        IF(.NOT.(l_noco.or.l_soc)) THEN
           write(oUnit,'(a,5x,a,2(" | ",5(a,5x)))') "Atom "," m    ","mx   ","my   ","mz   ","alpha","beta ","mx   ","my   ","mz   ","alpha","beta "
        ELSE
           write(oUnit,'(a,5x,a,2(" | ",5(a,5x)))') "Atom ","|m|   ","mx   ","my   ","mz   ","alpha","beta ","mx   ","my   ","mz   ","alpha","beta "
        END IF   
        CALL openXMLElement('magneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))


        DO iType = 1, atoms%ntype
           !Collect the data
           if (present(moments)) THEN
              up=moments%chmom(iType,1)
              down=moments%chmom(iType,2)
           elseif(present(den)) THEN
              CALL intgr3(den%mt(:,0,iType,1),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),up)
              CALL intgr3(den%mt(:,0,iType,2),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),down)
              up=up*sfp_const;down=down*sfp_const
           endif   
           off_diag=0.0
           l_offdiag=.false.
           if (present(moments).and.noco%l_mperp) THEN
               off_diag=moments%qa21(itype)
               l_offdiag=.true.
           elseif (present(den)) THEN
              if (size(den%mt,4)==4) THEN
                 CALL intgr3(den%mt(:,0,iType,3),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),off1)
                 CALL intgr3(den%mt(:,0,iType,4),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),off2)
                 off_diag=cmplx(off1,off2)*sfp_const
                 l_offdiag=.true.
              endif
           endif   
           !Determine local vector and angles
           local_m=nococonv%denmat_to_mag(up,down,off_diag)
           !Rotate to global frame
           if (noco%l_noco) then
                call nococonv%rotdenmat(itype, up, down, off_diag, toGlobal=.true.)
           else
                call nococonv%rotdenmat(nococonv%phi,nococonv%theta,up, down, off_diag, toGlobal=.true.)
           endif
           global_m=nococonv%denmat_to_mag(up,down,off_diag)

           call priv_print_mt_moment(itype,noco%l_noco,l_offdiag,noco%l_soc,local_m(1:3),global_m(1:3),"smm")
        enddo   
     
        CALL closeXMLElement('magneticMomentsInMTSpheres')
        write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
        write(oUnit,*)
        if (noco%l_soc) write(ounit,*) "SOC calculation: Global frame aligned to lattice"
     
    END SUBROUTINE 

    SUBROUTINE orbMoments(input,atoms,noco,nococonv,moments)
        !! Calculate and output the orbital moments
        !! Either the moments have to be given as calculated in cdnval (moments-argument)
        !! or the moments are calculated from the density 
        USE m_types
        USE m_constants
        USE m_xmlOutput
        
        USE m_intgr, ONLY : intgr3
        TYPE(t_input), INTENT(IN)           :: input
        TYPE(t_atoms), INTENT(IN)           :: atoms
        TYPE(t_noco), INTENT(IN)            :: noco
        TYPE(t_nococonv), INTENT(IN)        :: nococonv
        TYPE(t_moments),INTENT(IN)          :: moments
     
        INTEGER                       :: iType
        REAL                          :: local_m(0:3),global_m(3),theta,phi
        if (input%jspins==1) return ! No spin-pol calculation!
        

        !Header of output
        write(ounit,*)
        write(oUnit,'(a,19x,a,19x,a,19x,a)') "Orbital  moments |","Global Frame","  | ","Local Frame"
        write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
        write(oUnit,'(a,5x,a,2(" | ",5(a,5x)))') "Atom ","|m|   ","mx   ","my   ","mz   ","alpha","beta ","mx   ","my   ","mz   ","alpha","beta "
           
        CALL openXMLElement('orbitalMomentsInMTSpheres',(/'units'/),(/'muBohr'/))


        DO iType = 1, atoms%ntype
           !Collect the data
            global_m(:)=moments%clmom(:,iType,1)+moments%clmom(:,iType,2)
            if (noco%l_noco) THEN
                theta = nococonv%beta(iType)
                phi   = nococonv%alph(iType)
            else
                if (noco%l_soc) THEN    
                    theta = nococonv%theta
                    phi   = nococonv%phi
                else
                    theta=0.0
                    phi=0.0
                endif        
            endif
            
           !Rotate to local frame
           local_m(1:3)=global_m
           local_m(0)=0 
           call nococonv%rot_magvec(phi,theta,local_m, toGlobal=.false.)
           
           call priv_print_mt_moment(itype,.true.,.true.,.true.,local_m(1:3),global_m,"omm")
        enddo   
     
        CALL closeXMLElement('orbitalMomentsInMTSpheres')
        write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
        write(oUnit,*)
        
    END SUBROUTINE 
     
    subroutine priv_print_mt_moment(itype,l_noco,l_offdiag,l_soc,magmomL,magmom,grepstring)
        USE m_xmlOutput
        USE m_constants
        USE m_types_noco
        USE m_polangle
        integer,intent(in):: itype 
        logical,intent(in):: l_noco,l_offdiag,l_soc
        real,intent(in)   :: magmom(3)
        real,intent(in)   :: magmomL(3)
        character(len=*),INTENT(IN) :: grepstring
        
        
     
        character(len=15):: label 
        character(len=30):: attributes(2)      
        real :: alphal,betal,alpha,beta

        !Calculate angles
        CALL pol_angle(magmomL(1),magmomL(2),magmomL(3),betal,alphal,.true.)
        CALL pol_angle(magmom(1),magmom(2),magmom(3),beta,alpha,.true.)
                
        if (l_offdiag) THEN
           write(oUnit,'(i6,1x,f9.6,2(" | ",5(f9.6,1x)),"  mm ")') itype,sqrt(dot_product(magmom(1:3),magmom(1:3))),magmom,alpha,beta,magmomL,alphaL,betaL
        elseif (l_noco.or.l_soc) THEN
           write(oUnit,'(i6,1x,f9.6," | ",5(f9.6,1x)," |    ---       ---    ",f9.6,"    ---          ---           ",a)') &
            itype,sqrt(dot_product(magmom(1:3),magmom(1:3))),magmom,alpha,beta,magmomL(3),grepstring
        else
           write(oUnit,'(i6,1x,f9.6," |    ---       ---       ---       ---         ---   |    ---       ---       ---       ---          ---           ",a)') &
           itype,magmom(3),grepstring
        endif

        WRITE(attributes(1),'(i0)') iType
        WRITE(attributes(2),'(3(f9.5,1x))') magmom(1),magmom(2),magmom(3)
        label="globalMagMoment"
        CALL writeXMLElementFormPoly(label,(/'atomType','vec     '/),&
                        attributes,reshape((/8,3,6,30/),(/2,2/)))
        WRITE(attributes(2),'(3(f9.5,1x))') magmomL(1),magmomL(2),magmomL(3)
        label="localMagMoment "
        CALL writeXMLElementFormPoly(label,(/'atomType','vec     '/),&
                        attributes,reshape((/8,3,6,30/),(/2,2/)))
                               
     
    end subroutine
    subroutine priv_print_spin_density_at_nucleus(input,atoms,moments)
        USE m_types
        USE m_constants
        TYPE(t_input), INTENT(IN)           :: input
        TYPE(t_atoms), INTENT(IN)           :: atoms
        TYPE(t_moments),INTENT(IN)          :: moments
        INTEGER                       :: iType
        REAL                          :: sval,stot,scor
        WRITE (oUnit,FMT=8000)
        DO iType = 1,atoms%ntype
         sval = moments%svdn(iType,1) - moments%svdn(iType,input%jspins)
         stot = moments%stdn(iType,1) - moments%stdn(iType,input%jspins)
         scor = stot - sval
         WRITE (oUnit,FMT=8010) iType,stot,sval,scor,moments%svdn(iType,1),moments%stdn(iType,1)
        END DO
        8000 FORMAT (/,/,10x,'spin density at the nucleus:',/,10x,'type',t25,&
                 'total',t42,'valence',t65,'core',t90,&
                 'majority valence and total density',/)
        8010 FORMAT (i13,2x,3e20.8,5x,2e20.8)
    end subroutine
end module    
