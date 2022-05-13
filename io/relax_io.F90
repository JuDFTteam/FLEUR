!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relaxio
  !This module handles IO to the relax.xml file
  !The writing is done directly to relax.xml
  !The reading uses the libxml interface to inp.xml. Hence the relax.xml has to be included here.
  USE m_judft
  USE m_constants
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_relax,write_relax,apply_displacements,read_displacements,rotate_to_all_sites
CONTAINS
  SUBROUTINE write_relax(positions,forces,energies,displace)
    REAL,INTENT(in):: positions(:,:,:)
    REAL,INTENT(in):: forces(:,:,:)
    REAL,INTENT(in):: energies(:)
    REAL,INTENT(in):: displace(:,:)

    INTEGER :: no_steps,n,ntype,step
    No_steps=SIZE(positions,3)
    ntype=SIZE(positions,2)
    IF (ntype.NE.SIZE(forces,2).OR.ntype.NE.SIZE(displace,2).OR.&
         no_steps.NE.SIZE(forces,3).OR.no_steps.NE.SIZE(energies))THEN
       CALL judft_error("BUG in relax_io")
    ENDIF
    OPEN(765,file="relax.xml",status="replace")
    WRITE(765,*) "<!-- Attention, absolute coordinates used here -->"
    WRITE(765,*) "<relaxation>"
    !write current set of displacements
    WRITE(765,*) "  <displacements>"
    DO n=1,SIZE(displace,2)
       WRITE(765,"(a,3(f15.10,1x),a)") &
            '    <displace>',displace(:,n),'</displace>'
    END DO
    WRITE(765,"(a)") '  </displacements>'

    !Write all known old positions,forces and energies
    WRITE(765,*) "  <relaxation-history>"
    DO step=1,no_steps
       WRITE(765,"(a,f20.10,a)") '    <step energy="',energies(step),'">'
       DO n=1,ntype
          WRITE(765,"(a,6(f15.10,1x),a)") &
               '      <posforce>',positions(:,n,step),forces(:,n,step),'</posforce>'
       END DO
       WRITE(765,"(a)") '    </step>'
    ENDDO
    WRITE(765,*) "  </relaxation-history>"
    WRITE(765,*) "</relaxation>"
    CLOSE(765)
  END SUBROUTINE write_relax

  SUBROUTINE read_relax(positions,forces,energies)
    USE m_types_xml
    USE m_calculator
    REAL,INTENT(INOUT),ALLOCATABLE:: positions(:,:,:)
    REAL,INTENT(INOUT),ALLOCATABLE:: forces(:,:,:)
    REAL,INTENT(INOUT),ALLOCATABLE:: energies(:)

    REAL,ALLOCATABLE::rtmp(:,:,:)
    INTEGER:: no_steps
    INTEGER:: ntype,step,n
    CHARACTER(len=100):: path,p,str

    TYPE(t_xml)::xml
    call xml%init()
    no_steps=xml%GetNumberOfNodes('/fleurInput/relaxation/relaxation-history/step')
    ntype=SIZE(positions,2)
    IF (no_steps==0) THEN
       IF (.NOT.ALLOCATED(positions)) ALLOCATE(positions(0,0,0),forces(0,0,0),energies(0))
       RETURN
    END IF
    IF (ALLOCATED(positions)) THEN
       !Assume that we got already a set of positions, forces, energy and extend that list
       rtmp=positions
       DEALLOCATE(positions)
       ALLOCATE(positions(3,ntype,no_steps+SIZE(rtmp,3)))
       positions(:,:,no_steps+1:)=rtmp
       rtmp=forces
       DEALLOCATE(forces)
       ALLOCATE(forces(3,ntype,no_steps+SIZE(rtmp,3)))
       forces(:,:,no_steps+1:)=rtmp
       rtmp(1,1,:)=energies
       DEALLOCATE(energies)
       ALLOCATE(energies(no_steps+SIZE(rtmp,3)))
       energies(no_steps+1:)=rtmp(1,1,:)
    ELSE
       ALLOCATE(positions(3,ntype,no_steps))
       ALLOCATE(forces(3,ntype,no_steps))
       ALLOCATE(energies(no_steps))
    END IF
    DO step=1,no_steps
       WRITE(path,"(a,i0,a)") '/fleurInput/relaxation/relaxation-history/step[',step,']'
       energies(step)=evaluateFirstOnly(xml%GetAttributeValue(TRIM(path)//"/@energy"))
       DO n=1,ntype
          WRITE(p,"(a,a,i0,a)") TRIM(path),"/posforce[",n,"]"
          str=xml%GetAttributeValue(p)
          positions(:,n,step)=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
          Forces(:,n,step)=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
       ENDDO
    END DO
    call xml%FreeResources()
  END SUBROUTINE read_relax


  SUBROUTINE read_displacements(atoms,disp)
    USE m_types_xml
    USE m_calculator
    USE m_types
    TYPE(t_atoms),INTENT(in)::atoms
    REAL,INTENT(out)::disp(:,:)
    CHARACTER(len=50):: path,str
    INTEGER :: n

    TYPE(t_xml)::xml

    disp=0.0
    IF (xml%GetNumberOfNodes('/fleurInput/relaxation/displacements')==0) RETURN
    !read displacements and apply to positions
    IF (atoms%ntype.NE.xml%GetNumberOfNodes('/fleurInput/relaxation/displacements/displace')) CALL judft_error("Wrong number of displacements in relaxation")
    DO n=1,atoms%ntype
       WRITE(path,"(a,i0,a)") '/fleurInput/relaxation/displacements/displace[',n,']'
       str=xml%GetAttributeValue(path)
       disp(:,n)=(/evaluateFirst(str),evaluateFirst(str),evaluateFirst(str)/)
    END DO
  END SUBROUTINE read_displacements

  SUBROUTINE apply_displacements(cell,input,vacuum ,sym,noco,atoms,gfinp)
    USE m_types
    USE m_chkmt
    USE m_constants
    USE m_mapatom
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_cell),INTENT(IN)    :: cell
     
    TYPE(t_sym),INTENT(INOUT)  :: sym
    TYPE(t_noco),INTENT(IN)    :: noco
    type(t_gfinp), intent(in)  :: gfinp
    character(len=:), allocatable :: error_output

    TYPE(t_atoms),INTENT(INOUT):: atoms


    INTEGER:: n,indx(2)
    REAL   :: disp(3,atoms%ntype),disp_all(3,atoms%nat),taual0(3,atoms%nat),overlap(0:atoms%ntype,atoms%ntype)

    CALL read_displacements(atoms,disp)
    !change displacements to internal coordinates
    disp=MATMUL(cell%bmat,disp)/tpi_const
    IF (ALL(ABS(disp)<1E-8)) RETURN
    !Fist make sure original MT spheres do not overlap
    CALL chkmt(atoms,input,vacuum,cell ,.TRUE.,overlap=overlap)
    IF(ANY(overlap>0.0)) CALL judft_error("Overlapping MT-spheres in relaxation before even applying any displacement",hint="You messed up your setup")

    taual0=atoms%taual !Store original positions

    !Now check for overlapping mt-spheres
    overlap=1.0
    DO WHILE(ANY(overlap>1E-10))
       atoms%taual=taual0
       CALL rotate_to_all_sites(disp,atoms,cell,sym,disp_all)
       atoms%taual=taual0+disp_all
       atoms%pos=MATMUL(cell%amat,atoms%taual)
       CALL chkmt(atoms,input,vacuum,cell ,.TRUE.,overlap=overlap)
       IF (ANY(overlap>0.0)) THEN
          IF (ANY(overlap(0,:)>1E-10)) CALL judft_error("Atom spills out into vacuum during relaxation")
          indx=MAXLOC(overlap(1:,:)) !the two indices with the most overlap
          !Try only 90% of displacement
          disp(:,indx(1))=disp(:,indx(1))*0.9
          disp(:,indx(2))=disp(:,indx(2))*0.9
          WRITE(*,*) "Attention: Overlapping MT-spheres. Reduced displacement by 10%"
          WRITE(*,*) indx,overlap(indx(1),indx(2))
          WRITE(oUnit,'(a,2(i0,1x),f12.8)') "Attention, overlapping MT-spheres: ",indx,overlap(indx(1),indx(2))
          ! WRITE(error_output, '(3a,f12.8,a)') "Overlapping MT-spheres during relaxation: ", atoms%label(sum(atoms%neq(:indx(1)-1))+1),&
          ! &atoms%label(sum(atoms%neq(:indx(2)-1))+1), overlap(indx(1),indx(2)),&
          ! &NEW_LINE('A')//"Treat as an error: writing rescaled displacements to relax.xml is not implemented"
          error_output = "Overlapping MT-spheres during relaxation: " // NEW_LINE('A') // &
                         strip(atoms%label(sum(atoms%neq(:indx(1)-1))+1)) // " " // &
                         strip(atoms%label(sum(atoms%neq(:indx(2)-1))+1)) // &
                         " olap: "//float2str(overlap(indx(1),indx(2)))//NEW_LINE('A')//&
                         "Treat as an error: writing rescaled displacements to relax.xml is not implemented"
          CALL judft_error(error_output)
       END IF
    END DO

    CALL mapatom(sym,atoms,cell,input,noco,gfinp)

    WRITE(oUnit,*) "Atomic positions including displacements:"
    DO n=1,atoms%nat
       WRITE(oUnit,"(i4,6(1x,f12.5))") n,atoms%taual(:,n),atoms%pos(:,n)
    ENDDO

  END SUBROUTINE apply_displacements

  SUBROUTINE rotate_to_all_sites(disp,atoms,cell,sym,disp_all)
    USE m_types
    REAL,INTENT(in)          :: disp(:,:)
    TYPE(t_atoms),INTENT(in) :: atoms
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_sym),INTENT(IN)   :: sym
    REAL,INTENT(out)         :: disp_all(:,:)

    INTEGER:: iType,iAtom,jop, startAtom
    REAL   :: tau0(3),tau0_rot(3),tau_rot(3)


    DO iType = 1, atoms%ntype
       startAtom = SUM(atoms%neq(:iType-1))+1
       tau0=atoms%taual(:,startAtom)
       DO iAtom = startAtom, SUM(atoms%neq(:iType))
          jop = sym%invtab(sym%ngopr(iAtom))
          tau0_rot=MATMUL(1.*sym%mrot(:,:,jop),tau0)+sym%tau(:,jop) !translation will cancel, included for clarity
          tau_rot=MATMUL(1.*sym%mrot(:,:,jop),tau0+disp(:,iType))+sym%tau(:,jop)
          disp_all(:,iAtom)=tau_rot-tau0_rot
       END DO
    END DO
  END SUBROUTINE rotate_to_all_sites
END MODULE m_relaxio
