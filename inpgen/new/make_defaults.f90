!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_defaults
  USE m_juDFT
!---------------------------------------------------------------------
!  Check muffin tin radii and determine a reasonable choice for MTRs.
!  Derive also other parameters for the input file, to provide some
!  help in the out-file.                                        gb`02
!---------------------------------------------------------------------
CONTAINS

  subroutine make_FleurInputSchema()
     USE iso_c_binding
     implicit none
    integer :: errorStatus
    logical :: l_exist
    interface
       function dropInputSchema() bind(C, name="dropInputSchema")
         use iso_c_binding
         INTEGER(c_int) dropInputSchema
       end function dropInputSchema
    end interface
    
    inquire(file="FleurInputSchema",exist=l_exist)
    if (l_exist) return
    errorStatus = 0
    errorStatus = dropInputSchema()
    IF(errorStatus.NE.0) THEN
       call judft_error('Error: Cannot print out FleurInputSchema.xsd')
    END IF
  end subroutine make_FleurInputSchema


  
  SUBROUTINE make_defaults(atoms,sym,vacuum,input,stars,&
&                   xcpot,noco,hybrid)
      USE m_types
      TYPE(t_atoms),INTENT(IN)    ::atoms
      TYPE(t_sym),INTENT(IN)      ::sym
      
      TYPE(t_vacuum),INTENT(INOUT)::vacuum
      TYPE(t_input),INTENT(INOUT) ::input
      TYPE(t_stars),INTENT(INOUT) ::stars
      TYPE(t_xcpot),INTENT(INOUT) ::xcpot
      TYPE(t_noco),INTENT(INOUT)  ::noco
      TYPE(t_hybrid),INTENT(INOUT)::hybrid
      
      
      IMPLICIT NONE

      !
      !input
      !
      input%delgau = input%tkb
      IF (noco%l_noco) input%jspins = 2
       
      IF ( ANY(atoms%nlo(:).NE.0) ) THEN
        input%ellow = -1.8
      ELSE
        input%ellow = -0.8  
      ENDIF
      IF (input%film) THEN
         input%elup = 0.5
      ELSE
         input%elup = 1.0
      ENDIF 
      IF (l_hyb) THEN
         input%ellow = input%ellow -  2.0
         input%elup  = input%elup  + 10.0
         input%gw_neigd = max( nint(input%zelec)*10, 60 )
         input%minDistance = 1.0e-5
         input%ctail = .false.
      ELSE
        input%gw_neigd = 0
      END IF

      input%rkmax   = real(NINT(input%rkmax   * 10  ) / 10.)
      IF (noco%l_ss) input%ctail = .FALSE.
 
      !
      ! stars
      !
      stars%gmax     = merge(stars%gmax,3.0*input%rkmax,stars%gmax>0)
      stars%gmax     = real(NINT(stars%gmax    * 10  ) / 10.)
      stars%gmaxInit = stars%gmax

      !
      !xcpot
      !
      xcpot%gmaxxc  = merge(xcpot%gmaxxc,3.0*input%rkmax,xcpot%gmaxxc>0)
      xcpot%gmaxxc  = real(NINT(xcpot%gmaxxc  * 10  ) / 10.)
      

      !
      !vacuum
      !
      IF (.not.input%film) THEN
         vacuum%dvac = a3(3) 
      Else
         vacuum%dvac = REAL(NINT(vacuum%dvac*100)/100.)
      ENDIF
      vacuum%nvac = 2
      IF (sym%zrfs.OR.sym%invs) vacuum%nvac = 1
      IF (oneD%odd%d1) vacuum%nvac = 1
   
      !
      !noco
      !
      ALLOCATE(noco%l_relax(atoms%ntype),noco%b_con(2,atoms%ntype))
      ALLOCATE(noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype))
      noco%qss = MERGE(noco%qss,[0.0,0.0,0.0],noco%l_ss)
      noco%l_relax(:) = .FALSE.
      noco%alphInit(:) = 0.0
      noco%alph(:) = 0.0
      noco%beta(:) = 0.0
      noco%b_con(:,:) = 0.0
      
      !
      !hybrid
      !
      hybrid%gcutm1       = input%rkmax - 0.5
      ALLOCATE(hybrid%lcutwf(atoms%ntype))
      ALLOCATE(hybrid%lcutm1(atoms%ntype))
      ALLOCATE(hybrid%select1(4,atoms%ntype))
      hybrid%lcutwf      = atoms%lmax - atoms%lmax / 10
      hybrid%lcutm1      = 4
      hybrid%select1(1,:) = 4
      hybrid%select1(2,:) = 0
      hybrid%select1(3,:) = 4
      hybrid%select1(4,:) = 2
      hybrid%l_hybrid = l_hyb
      hybrid%gcutm1 = real(NINT(hybrid%gcutm1 * 10  ) / 10.)
    

    END SUBROUTINE make_defaults
  END MODULE m_make_defaults
