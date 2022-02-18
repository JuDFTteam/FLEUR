!--------------------------------------------------------------------------------
! Copyright (c) 2020 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_magMultipoles
  USE m_types
  USE m_juDFT
  USE m_constants
  IMPLICIT NONE
CONTAINS
  SUBROUTINE magMultipoles(sym,stars, atoms,cell, sphhar, vacuum, input, noco,nococonv,outden)
    USE m_plot
    USE m_divergence
    USE m_mpmom
    TYPE(t_input),INTENT(IN)                    :: input
    TYPE(t_atoms), INTENT(IN)                   :: atoms
    TYPE(t_sphhar), INTENT(IN)                  :: sphhar
    TYPE(t_sym), INTENT(IN)                     :: sym
    TYPE(t_noco), INTENT(IN)                    :: noco
    TYPE(t_nococonv), INTENT(IN)                :: nococonv
    TYPE(t_stars),INTENT(IN)                    :: stars
    TYPE(t_cell),INTENT(IN)                     :: cell
    TYPE(t_vacuum),INTENT(IN)                   :: vacuum
    TYPE(t_potden), INTENT(in)                  :: outden


    TYPE(t_potden) :: cden,m_den(3),div
    COMPLEX:: qlmo(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    INTEGER:: n,l,m


    CALL cden%init(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,0,input%jspins,.FALSE.,.FALSE.,1001,&
         vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
    DO n=1,3
       CALL m_den(n)%init(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,0,input%jspins,.FALSE.,.FALSE.,1001,&
            vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
    ENDDO
    CALL div%init(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,0,input%jspins,.FALSE.,.FALSE.,1001,&
    vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
    allocate(div%pw_w,mold=div%pw)

    !Generate magnetization out of density matrix
    CALL matrixsplit(sym,stars, atoms, sphhar, vacuum, input, noco,nococonv, 1.0, &
         outden, cden, m_den(1), m_den(2), m_den(3))
    !Calcalate divergence
    CALL divergence(input,stars,atoms,sphhar,vacuum,sym,cell,noco,m_den,div)
    qlmo = 0.0
    CALL mt_moments( input, atoms, sym,sphhar, div%mt(:,:,:,1), POTDEN_TYPE_POTCOUL,qlmo,.FALSE.)

    WRITE(oUnit,*) "Magnetic Multipoles:"
    DO n=1,atoms%ntype
       WRITE(oUnit,*) "Atom type:",n
       DO l=0,4
          WRITE(oUnit,"(10(2f12.7,3x))") (qlmo(m,l,n),m=-l,l)
       ENDDO
    ENDDO
  END SUBROUTINE magMultipoles
END MODULE m_magMultipoles
