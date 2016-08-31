!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_xsf_io
  USE m_types
  !-----------------------------------------------
  ! DESC:subroutines to write xsf-files for xcrysden
  !                 Daniel Wortmann, (06-01-26)
  !-----------------------------------------------
  ! Bohr radius a0, http://physics.nist.gov/cgi-bin/cuu/Value?eqbohrrada0
  REAL, PRIVATE, PARAMETER :: a0 = 0.52917720859 ! in Angstroem

CONTAINS
  !<-- S:S: xsf_WRITE_atoms(fileno,film,amat,neq(:ntype),zatom(:ntype),pos)
  SUBROUTINE xsf_WRITE_atoms(fileno,atoms,film,od,amat)
    !-----------------------------------------------
    !     Writes the crystal dimensions&atomic positions
    !           (last modified: 2004-00-00) D. Wortmann
    !-----------------------------------------------
    IMPLICIT NONE
    !<--Arguments
    INTEGER,INTENT(IN)     :: fileno
    TYPE(t_atoms),INTENT(IN):: atoms
    LOGICAL,INTENT(IN)      :: film
    LOGICAL,INTENT(IN)     :: od
    REAL,INTENT(IN)        :: amat(3,3)
    !>
    !<-- Locals
    INTEGER             :: n,nn,na
    !>
    IF (film) THEN
       IF (od) THEN
          WRITE(fileno,*) "POLYMERE"
       ELSE
          WRITE(fileno,*) "SLAB"
       ENDIF
    ELSE
       WRITE(fileno,*) "CRYSTAL"
    ENDIF

    WRITE(fileno,*) "PRIMVEC"
    ! Write in atomic units
    WRITE(fileno,'(3(f0.7,1x))') amat(:,1)*a0
    WRITE(fileno,'(3(f0.7,1x))') amat(:,2)*a0
    WRITE(fileno,'(3(f0.7,1x))') amat(:,3)*a0

    WRITE(fileno,*) "PRIMCOORD"
    WRITE(fileno,*) SUM(atoms%neq)," 1"
    na = 1
    DO n = 1,SIZE(atoms%neq)
       DO nn = 1,atoms%neq(n)
          WRITE(fileno,'(i4,2x,3(f0.7,1x))') NINT(atoms%zatom(n)),&
               &            atoms%pos(:,na)*a0
          na=na+1
       ENDDO
    ENDDO
    WRITE(fileno,*)
  END SUBROUTINE xsf_WRITE_atoms
  !> 
  !<-- S: xsf_write_header(fileno,twodim,desc,vec1,vec2,vec3,zero,grid)
  SUBROUTINE xsf_WRITE_header(fileno,twodim,desc,vec1,vec2,vec3,zero&
       &     ,grid)
    !-----------------------------------------------
    !  writes the beginning of a gid-datablock
    !           (last modified: 2004-00-00) D. Wortmann
    !-----------------------------------------------
    IMPLICIT NONE
    !<--Arguments
    INTEGER,INTENT(IN)     :: fileno,grid(:)
    LOGICAL,INTENT(IN)     :: twodim
    REAL   ,INTENT(IN)     :: vec1(:),vec2(:),vec3(:),zero(:)
    CHARACTER(LEN =*),INTENT(IN) :: desc 
    !>

    IF (twodim) THEN
       WRITE(fileno,*) "BEGIN_BLOCK_DATAGRID_2D"
       WRITE(fileno,*) desc
       WRITE(fileno,*) "BEGIN_DATAGRID_2D_A"
       WRITE(fileno,'(3i7)') grid(1:2)
       WRITE(fileno,'(3(f12.7,1x))') zero*a0
       WRITE(fileno,'(3(f12.7,1x))') vec1*a0
       WRITE(fileno,'(3(f12.7,1x))') vec2*a0
    ELSE
       WRITE(fileno,*) "BEGIN_BLOCK_DATAGRID_3D"
       WRITE(fileno,*) desc
       WRITE(fileno,*) "BEGIN_DATAGRID_3D_A"
       WRITE(fileno,'(3i7)') grid(1:3)
       WRITE(fileno,'(3(f12.7,1x))') zero*a0
       WRITE(fileno,'(3(f12.7,1x))') vec1*a0
       WRITE(fileno,'(3(f12.7,1x))') vec2*a0
       WRITE(fileno,'(3(f12.7,1x))') vec3*a0
    ENDIF
  END SUBROUTINE xsf_WRITE_header
  !> 
  !<-- S: xsf_write_newblock(fileno,twodim,vec1,vec2,vec3,zero,grid)
  SUBROUTINE xsf_WRITE_newblock(fileno,twodim,vec1,vec2&
       &     ,vec3,zero,grid)
    !-----------------------------------------------
    !  writes the beginning of a new gid-datablock for second spin
    !           (last modified: 2004-00-00) D. Wortmann
    !-----------------------------------------------
    IMPLICIT NONE
    !<--Arguments
    INTEGER,INTENT(IN)     :: fileno,grid(:)
    LOGICAL,INTENT(IN)     :: twodim
    REAL   ,INTENT(IN)     :: vec1(:),vec2(:),vec3(:),zero(:)
    !>

    IF (twodim) THEN
       WRITE(fileno,*) "END_DATAGRID_2D"
       WRITE(fileno,*) "BEGIN_DATAGRID_2D_B"
       WRITE(fileno,'(3i7)') grid(1:2)
       WRITE(fileno,'(3(f12.7,1x))') zero*a0
       WRITE(fileno,'(3(f12.7,1x))') vec1*a0
       WRITE(fileno,'(3(f12.7,1x))') vec2*a0
    ELSE
       WRITE(fileno,*) "END_DATAGRID_3D"
       WRITE(fileno,*) "BEGIN_DATAGRID_3D_B"
       WRITE(fileno,'(3i7)') grid(1:3)
       WRITE(fileno,'(3(f12.7,1x))') zero*a0
       WRITE(fileno,'(3(f12.7,1x))') vec1*a0
       WRITE(fileno,'(3(f12.7,1x))') vec2*a0
       WRITE(fileno,'(3(f12.7,1x))') vec3*a0
    ENDIF
  END SUBROUTINE xsf_WRITE_newblock
  !> 
  !<-- S: xsf_write_endblock(fileno,twodim)
  SUBROUTINE xsf_write_endblock(fileno,twodim)
    !-----------------------------------------------
    !
    !           (last modified: 2004-00-00) D. Wortmann
    !-----------------------------------------------
    IMPLICIT NONE
    !<--Arguments
    INTEGER,INTENT(IN)     :: fileno
    LOGICAL,INTENT(IN)     :: twodim
    !>

    IF (twodim) THEN
       WRITE(fileno,*) "END_DATAGRID_2D"                      
       WRITE(fileno,*) "END_BLOCK_DATAGRID_2D" 
    ELSE
       WRITE(fileno,*) "END_DATAGRID_3D"                      
       WRITE(fileno,*) "END_BLOCK_DATAGRID_3D" 
    ENDIF
  END SUBROUTINE xsf_write_endblock
  !> 

  SUBROUTINE xsf_WRITE_force(fileno,atoms,film,od,amat,force)
    !-----------------------------------------------
    !     Writes the crystal dimensions&force positions
    !           (last modified: 2004-00-00) D. Wortmann
    !-----------------------------------------------
    IMPLICIT NONE
    !<--Arguments
    INTEGER,INTENT(IN)       :: fileno
    TYPE(t_atoms),INTENT(IN) :: atoms
    LOGICAL,INTENT(IN)       :: film
    LOGICAL,INTENT(IN)       :: od
    REAL,INTENT(IN)          :: amat(3,3)
    INTEGER,INTENT(IN)     :: force ! number of atoms + force vectors
    !>
    !<-- Locals
    INTEGER             :: n,nn,na
    !>
    IF (film) THEN
       IF (od) THEN
          WRITE(fileno,*) "POLYMERE"
       ELSE
          WRITE(fileno,*) "SLAB"
       ENDIF
    ELSE
       WRITE(fileno,*) "CRYSTAL"
    ENDIF

    WRITE(fileno,*) "PRIMVEC"
    WRITE(fileno,'(3(f0.7,1x))') amat(:,1)*a0
    WRITE(fileno,'(3(f0.7,1x))') amat(:,2)*a0
    WRITE(fileno,'(3(f0.7,1x))') amat(:,3)*a0

    WRITE(fileno,*) "PRIMCOORD"
    WRITE(fileno,*) force," 1"
    na = 1
    DO n = 1,SIZE(atoms%neq)
       DO nn = 1,atoms%neq(n)
          WRITE(fileno,'(i4,2x,3(f0.7,1x))') NINT(atoms%zatom(n)),&
               &            atoms%pos(:,na)*a0
          na=na+1
       ENDDO
    ENDDO
    WRITE(fileno,*)
  END SUBROUTINE xsf_WRITE_force
  !> 
  !-----------------------------------------------
END MODULE m_xsf_io

 
