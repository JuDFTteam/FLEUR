!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relaxio
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_relax,write_relax
CONTAINS
  SUBROUTINE write_relax(positions,forces,energies,displace)
    REAL,INTENT(in):: positions(:,:,:)
    REAL,INTENT(in):: forces(:,:,:)
    REAL,INTENT(in):: energies(:)
    REAL,INTENT(in):: displace(:,:)

    INTEGER :: no_steps,n,ntype,step
    No_steps=SIZE(positions,3)
    ntype=SIZE(positions,2)
    IF (ntype.NE.SIZE(forces,2).OR.ntype<SIZE(displace,2).OR.&
         no_steps.NE.SIZE(forces,3).OR.no_steps.NE.SIZE(energies)) &
         CALL judft_error("BUG in relax_io")
    OPEN(765,file="relax.inp",status="replace")
    WRITE(765,*) "<relaxation>"
    !write current set of displacements
    WRITE(765,*) "  <displacements>"
    DO n=1,SIZE(displace,2)
       WRITE(765,"(a,i0,a,3(f15.10,1x),a)") &
            '    <displace   na="',n,'">',displace(:,n),'</displace>'
    END DO
    WRITE(765,"(a)") '  </displacements>'
    
    !Write all known old positions,forces and energies
    WRITE(765,*) "  <relaxation-history>"
    DO step=1,no_steps
       WRITE(765,"(a,f20.10,a)") '    <step energy="',energies(step),'">'
       DO n=1,ntype
          WRITE(765,"(a,i0,a,3(f15.10,1x),a)") &
               '      <pos   n="',n,'">',positions(:,n,step),'</pos>'
          WRITE(765,"(a,i0,a,3(f15.10,1x),a)") &
               '      <force n="',n,'">',forces(:,n,step),'</force>'
       END DO
       WRITE(765,"(a)") '    </step>'
    ENDDO
    WRITE(765,*) "  </relaxation-history>"
    WRITE(765,*) "</relaxation>"
    CLOSE(765)
  END SUBROUTINE write_relax
  
  SUBROUTINE read_relax(positions,forces,energies)
    USE m_xmlIntWrapFort 
    USE m_calculator
    REAL,INTENT(OUT),ALLOCATABLE:: positions(:,:,:)
    REAL,INTENT(OUT),ALLOCATABLE:: forces(:,:,:)
    REAL,INTENT(OUT),ALLOCATABLE:: energies(:)
    
    INTEGER:: no_steps
    INTEGER:: ntype,step,n
    CHARACTER(len=50):: path,p,str
    no_steps=xmlGetNumberOfNodes('/fleurInput/relaxation/relaxation-history/step')
    IF (no_steps==0) THEN
       ALLOCATE(positions(0,0,0),forces(0,0,0),energies(0))
       RETURN
    ELSE
       ntype=xmlGetNumberOfNodes('/fleurInput/relaxation/relaxation-history/step[1]/pos')
    END IF
    ALLOCATE(positions(3,ntype,no_steps))
    ALLOCATE(forces(3,ntype,no_steps))
    ALLOCATE(energies(no_steps))
    DO step=1,no_steps
       WRITE(path,"(a,i0,a)") '/fleurInput/relaxation/relaxation-history/step[',step,']'
       energies(step)=evaluateFirstOnly(xmlGetAttributeValue(TRIM(path)//"/energy"))
       DO n=1,ntype
          WRITE(p,"(a,a,i0,a)") TRIM(path),"/pos[",n,"]"
          str=xmlGetAttributeValue(p)
          positions(:,n,step)=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
          WRITE(p,"(a,a,i0,a)") TRIM(path),"/force[",n,"]"
          str=xmlGetAttributeValue(p)
          forces(:,n,step)=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
       ENDDO
    END DO
  END SUBROUTINE read_relax

  SUBROUTINE read_displacements(cell,atoms)
    USE m_xmlIntWrapFort 
    USE m_calculator
    USE m_types
    TYPE(t_cell),INTENT(INOUT) :: cell
    TYPE(t_atoms),INTENT(INOUT):: atoms

    
    INTEGER:: n,na
    REAL   :: d(3)
    CHARACTER(len=50):: path,str
    
    IF (xmlGetNumberOfNodes('/fleurInput/relaxation/displacements').NE.0) THEN
       !read displacements and apply to positions
       DO n=1,xmlGetNumberOfNodes('/fleurInput/relaxation/displacements/displace')
          WRITE(path,"(a,i0,a)") '/fleurInput/relaxation/displacements/displace[',n,']'
          str=xmlGetAttributeValue(path)
          d=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
          na=INT(evaluateFirstOnly(xmlGetAttributeValue(TRIM(path)//"/@na")))
          IF (na>SIZE(atoms%taual,2)) CALL judft_error&
               ("Wrong number of displacements. na can not be larger than number of atoms in setup")
          atoms%taual(:,na)=atoms%taual(:,na)+d
          atoms%pos(:,na) = matmul(cell%amat,atoms%taual(:,na))
       END DO
    END IF
  END SUBROUTINE read_displacements
END MODULE m_relaxio
