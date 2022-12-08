!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdninf
CONTAINS
  SUBROUTINE cdninf(input,sym,noco,atoms,vacuum,&
                    cell,kpts,eigdos)
    !***********************************************************************
    !     this subroutine calculates the charge distribution of each state
    !     and writes this information to the out file. If dos or vacdos
    !     are .true. it also write the necessary information for dos or
    !     bandstructure plots to the file dosinp and vacdos respectivly
    !***********************************************************************
    !       changed this subroutine slightly for parallisation of dosinp&
    !       vacdos output (argument z replaced by ksym,jsym, removed sympsi
    !       call)                                        d.wortmann 5.99
    !
    !******** ABBREVIATIONS ************************************************
    !     qal      : l-like charge of each state
    !     qvac     : vacuum charge of each state
    !     qvlay    : charge in layers (z-ranges) in the vacuum of each state
    !     starcoeff: T if star coefficients have been calculated
    !     qstars   : star coefficients for layers (z-ranges) in vacuum
    !
    !***********************************************************************
    USE m_types
    USE m_types_dos
    USE m_types_vacdos
    USE m_types_eigdos
    USE m_constants
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_vacuum),INTENT(IN)      :: vacuum
    TYPE(t_noco),INTENT(IN)        :: noco
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_cell),INTENT(IN)        :: cell
    TYPE(t_kpts),INTENT(IN)        :: kpts
    TYPE(t_atoms),INTENT(IN)       :: atoms
    CLASS(t_eigdos_list),INTENT(IN)   :: eigdos(:)
    CLASS(t_eigdos),pointer :: dos
    Type(t_vacdos),pointer  :: vacdos
    !     ..
    !     .. Local Scalars ..
    REAL qalmax,qishlp,qvacmt,qvact
    INTEGER i,iband,ilay,iqispc,iqvacpc,ityp,itypqmax,ivac,l,lqmax
    INTEGER ikpt,jspin
    !     ..
    !     .. Local Arrays ..
    INTEGER iqalpc(0:3,atoms%ntype),max_l_type(2)
    CHARACTER chstat(1:4)
    !     ..
    !     .. Data statements ..
    DATA chstat/'s','p','d','f'/
    !     ..

    dos=>eigdos(1)%p

    vacdos=>null()
    if (size(eigdos)>1) THEN
      associate(vd=>eigdos(2)%p)
      select type(vd)
      type is (t_vacdos)
        vacdos=>vd
      end select  
      end associate
    endif  

    select type(dos)
    type is (t_dos)

    DO jspin=1,input%jspins
      DO ikpt=1,kpts%nkpt

    IF (input%film) THEN
       WRITE (oUnit,FMT=8000) (kpts%bk(i,ikpt),i=1,2)
8000   FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,&
            &          'int',t22,'vac',t28,'spheres(s,p,d,f)')
    ELSE
       WRITE (oUnit,FMT=8010) (kpts%bk(i,ikpt),i=1,3)
8010   FORMAT (/,3x,'q(atom,l): k=',3f10.5,/,/,t8,'e',t13,'max',t18,&
            &          'int',t24,'spheres(s,p,d,f)')
    END IF
8020 FORMAT (1x,3e20.12,i6,e20.12)

    DO iband = 1,count(dos%eig(:,ikpt,jspin)<1E99)
      if (associated(vacdos)) THEN
        qvact=sum(vacdos%qvac(iband,:,ikpt,jspin))
      else
        qvact = 0
      endif
       iqvacpc = NINT(qvact*100.0)
       !qvacmt = qvact
       QVACMT=0.0
       iqalpc(0:3,:) = NINT(dos%qal(0:3,:,iband,ikpt,jspin)*100.0)
       DO l=0,3
         qvacmt=qvacmt+dot_product(dos%qal(l,:,iband,ikpt,jspin),atoms%neq)
       ENDDO
       max_l_type=maxloc(dos%qal(0:3,:,iband,ikpt,jspin))
       qishlp = 1.0 - qvacmt
       IF (noco%l_noco) qishlp = dos%qis(iband,ikpt,jspin)
       iqispc = NINT(qishlp*100.0)

       IF (input%film) THEN
          WRITE (oUnit,FMT=8040) dos%eig(iband,ikpt,jspin),chstat(max_l_type(1)),max_l_type(2),&
               &        iqispc,iqvacpc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
8040      FORMAT (f10.4,2x,a1,i2,2x,2i3, (t26,6 (4i3,1x)))
       ELSE
          WRITE (oUnit,FMT=8080) dos%eig(iband,ikpt,jspin),chstat(max_l_type(1)),max_l_type(2),&
               &        iqispc, ((iqalpc(l,ityp),l=0,3),ityp=1,atoms%ntype)
8080      FORMAT (f10.4,2x,a1,i2,2x,i3, (t26,6 (4i3,1x)))
       END IF
    END DO
  ENDDO
ENDDO
end select
  END SUBROUTINE cdninf
END MODULE m_cdninf
