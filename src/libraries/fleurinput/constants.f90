!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_constants
  use m_types_fleurinput_base
  IMPLICIT NONE

  TYPE,EXTENDS(t_fleurinput_base)::t_constants
   CONTAINS
     PROCEDURE:: read_xml=>read_xml_constants
     PROCEDURE:: mpi_bc=>mpi_bc_constants
  END TYPE t_constants

  REAL                        :: warp_factor=1.0  !should be set from input later
  INTEGER,          PARAMETER :: oUnit = 50 !This is the unit number for the out file
  INTEGER,          PARAMETER :: noState_const = 0
  INTEGER,          PARAMETER :: coreState_const = 1
  INTEGER,          PARAMETER :: valenceState_const = 2
  INTEGER,          PARAMETER :: lmaxU_const = 3
  COMPLEX,          PARAMETER :: cmplx_0=(0.0,0.0)
  COMPLEX,          PARAMETER :: cmplx_1=(1.0,0.0)
  REAL,             PARAMETER :: pi_const=3.1415926535897932
  REAL,             PARAMETER :: tpi_const=2.*3.1415926535897932
  REAL,             PARAMETER :: fpi_const=4.*3.1415926535897932
  REAL,             PARAMETER :: sfp_const=SQRT(4.*3.1415926535897932)
  real,             parameter :: sqrt_2=sqrt(2.0)
  COMPLEX,          PARAMETER :: ImagUnit=(0.0,1.0)
  REAL,             PARAMETER :: hartree_to_ev_const=27.211386245988 ! value from  https://physics.nist.gov/cgi-bin/cuu/Value?hrev
                                                                     ! 27.211386245988(53)
  REAL,             PARAMETER :: bohr_to_angstrom_const=0.529177210903 ! value from https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
                                                                       ! 0.529177210903(80)
  REAL,             PARAMETER :: eVac0Default_const = -0.25
  CHARACTER(len=9), PARAMETER :: version_const = 'fleur 41'
  CHARACTER(len=49), PARAMETER :: version_const_MaX = '     MaX-Release 8.0          (www.max-centre.eu)'
  CHARACTER(len=4), PARAMETER :: inputFileVersion_const = '0.39'
  ! outputFileVersion_const is defined in juDFT/xmlOutput.f90
  ! For version number update also update version numbers in files io/xml/FleurInputSchema.xsd, io/xml/FleurOutputSchema.xsd, fleurInput/types_xml, and update the files
  ! io/xml/inputSchema.h.backup, io/xml/outputSchema.h.backup, dropInputSchema.c, dropOutputSchema.c, inputSchema_old.h, outputSchema_old.h.
  REAL, PARAMETER             :: boltzmann_const = 3.1668114e-6 ! value is given in Hartree/Kelvin
  REAL, PARAMETER             :: r0_const = 0.0000236216 ! r_0 as in approximation of the nuclear radius R = r_0*A^(1/3), value in Bohr radii.

  INTEGER, PARAMETER :: POTDEN_TYPE_OTHER        = 0    ! POTDEN_TYPE <= 0 ==> undefined
  INTEGER, PARAMETER :: POTDEN_TYPE_POTTOT       = 1    ! 0 < POTDEN_TYPE <= 1000 ==> potential
  INTEGER, PARAMETER :: POTDEN_TYPE_POTCOUL      = 2
  INTEGER, PARAMETER :: POTDEN_TYPE_POTX         = 3
  INTEGER, PARAMETER :: POTDEN_TYPE_POTYUK       = 4
  INTEGER, PARAMETER :: POTDEN_TYPE_EnergyDen    = 5
  INTEGER, PARAMETER :: POTDEN_TYPE_CRYSTALFIELD = 6
  INTEGER, PARAMETER :: POTDEN_TYPE_DEN          = 1001 ! 1000 < POTDEN_TYPE ==> density

  INTEGER, PARAMETER :: KPTS_KIND_UNSPECIFIED = 0
  INTEGER, PARAMETER :: KPTS_KIND_MESH        = 1
  INTEGER, PARAMETER :: KPTS_KIND_PATH        = 2
  INTEGER, PARAMETER :: KPTS_KIND_TRIA_BULK   = 3
  INTEGER, PARAMETER :: KPTS_KIND_TRIA        = 4
  INTEGER, PARAMETER :: KPTS_KIND_SPEX_MESH   = 5

  CHARACTER(LEN=11),DIMENSION(0:5),PARAMETER :: kptsKindString_consts = &
     (/'unspecified','mesh       ','path       ','tria-bulk  ','tria       ','SPEX mesh  ' /)

  INTEGER, PARAMETER :: BZINT_METHOD_HIST  = 0
  INTEGER, PARAMETER :: BZINT_METHOD_GAUSS = 1
  INTEGER, PARAMETER :: BZINT_METHOD_TRIA  = 2
  INTEGER, PARAMETER :: BZINT_METHOD_TETRA = 3

   ! 11 == plot not yet implemented
   INTEGER, PARAMETER :: PLOT_INPDEN=1
   INTEGER, PARAMETER :: PLOT_INPDEN_N_CORE=11
   INTEGER, PARAMETER :: PLOT_POT_TOT=2
   INTEGER, PARAMETER :: PLOT_POT_EXT=11
   INTEGER, PARAMETER :: PLOT_POT_COU=3
   INTEGER, PARAMETER :: PLOT_POT_VXC=4
   INTEGER, PARAMETER :: PLOT_OUTDEN_Y_CORE=5
   INTEGER, PARAMETER :: PLOT_MIXDEN_Y_CORE=6
   INTEGER, PARAMETER :: PLOT_MIXDEN_N_CORE=11

  CHARACTER(2),DIMENSION(0:103),PARAMETER :: namat_const=(/&
       'va',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',&
       'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',&
       ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',&
       'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',&
       'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',&
       'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',&
       'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',&
       'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',&
       'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/)

  CHARACTER(7),DIMENSION(29),PARAMETER :: coreStateList_const=(/&
       '(1s1/2)','(2s1/2)','(2p1/2)','(2p3/2)','(3s1/2)',&
       '(3p1/2)','(3p3/2)','(4s1/2)','(3d3/2)','(3d5/2)',&
       '(4p1/2)','(4p3/2)','(5s1/2)','(4d3/2)','(4d5/2)',&
       '(5p1/2)','(5p3/2)','(6s1/2)','(4f5/2)','(4f7/2)',&
       '(5d3/2)','(5d5/2)','(6p1/2)','(6p3/2)','(7s1/2)',&
       '(5f5/2)','(5f7/2)','(6d3/2)','(6d5/2)' /)

  INTEGER,DIMENSION(29),PARAMETER :: coreStateNumElecsList_const=(/& ! This is the number of electrons per spin
       1, 1, 1, 2, 1, 1, 2, 1, 2, 3, 1, 2, 1, 2,&
       3, 1, 2, 1, 3, 4, 2, 3, 1, 2, 1, 3, 4, 2, 3/)

  INTEGER,DIMENSION(29),PARAMETER :: coreStateNprncList_const=(/&
       1, 2, 2, 2, 3, 3, 3, 4, 3, 3, 4, 4, 5, 4, 4,&
       5, 5, 6, 4, 4, 5, 5, 6, 6, 7, 5, 5, 6, 6/)
  INTEGER,DIMENSION(29),PARAMETER :: coreStateKappaList_const=(/&
       -1,-1, 1,-2,-1, 1,-2,-1, 2,-3, 1,-2,-1, 2,-3,&
       1,-2,-1, 3,-4, 2,-3, 1,-2,-1, 3,-4, 2,-3/)

  ! Nonrelativistic notation for core electron states
  CHARACTER(4),DIMENSION(18),PARAMETER :: nr_coreStateList_const=(/&
       '(1s)','(2s)','(2p)','(3s)','(3p)','(4s)','(3d)','(4p)','(5s)','(4d)','(5p)','(6s)','(4f)','(5d)','(6p)','(7s)','(5f)','(6d)' /)
 
  INTEGER,DIMENSION(29),PARAMETER :: coreStateToNRList_const=(/&
       1, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 14, 14, 15, 15, 16, 17, 17, 18, 18 /)
  

  REAL, DIMENSION(118), PARAMETER :: atomicMasses_const = (/  1.01,   4.00,   6.94,   9.01,  10.81,  12.01,  14.01,  16.00,  19.00, & ! up to fluorine
                                                             20.18,  22.99,  24.31,  26.98,  28.09,  30.97,  32.06,  35.45,  39.95, & ! up to argon
                                                             39.10,  40.08,  44.96,  47.87,  50.94,  52.00,  54.94,  55.85,  58.93, & ! up to cobalt
                                                             58.69,  63.55,  65.38,  69.72,  72.63,  74.92,  78.97,  79.90,  83.80, & ! up to krypton
                                                             85.47,  87.62,  88.91,  91.22,  92.91,  95.95,  97.40, 101.07, 102.91, & ! up to ruthenium
                                                            106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29, & ! up to xenon
                                                            132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 146.00, 150.36, 151.96, & ! up to europium
                                                            157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, & ! up to hafnium
                                                            180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, & ! up to thallium
                                                            207.20, 208.98, 209.98, 210.00, 222.00, 223.00, 226.00, 227.00, 232.04, & ! up to thorium
                                                            231.04, 238.03, 237.00, 244.00, 243.00, 247.00, 247.00, 251.00, 252.00, & ! up to einsteinium
                                                            257.00, 258.00, 259.00, 262.00, 267.00, 269.00, 270.00, 272.00, 273.00, & ! up to hassium
                                                            277.00, 281.00, 281.00, 285.00, 286.00, 289.00, 288.00, 293.00, 294.00, 294.00/)

  real, parameter :: massInElectronMasses = 1836.15 ! Ratio of proton to electron mass 

  CHARACTER(4),DIMENSION(6),PARAMETER :: nobleGasConfigList_const=(/'[He]','[Ne]','[Ar]','[Kr]','[Xe]','[Rn]'/)

  INTEGER,DIMENSION(6),PARAMETER :: nobleGasNumStatesList_const=(/1, 4, 7, 12, 17, 24/)

  ! juPhon stuff:
  complex, dimension(3, 3), parameter  :: Tmatrix0 = reshape( &
                                              [(0.7071067811865475, 0.),(0.,-0.7071067811865475),(0., 0.),&
                                               (0., 0.),(0., 0.),(1., 0.), &
                                               (-0.7071067811865475,0.),(0.,-0.7071067811865475),(0., 0.)],[3,3])!< Klüppelberg PhD thesis 4.18
  complex, dimension(3, 3), parameter :: c_im =   reshape([ &
                                                (1.4472025091165353, 0.),(0., 1.4472025091165353),(0., 0.),&
                                                 (0., 0.), (0., 0.),(2.0466534158929770, 0.),&
                                                (-1.4472025091165353, 0.),(0., 1.4472025091165353),(0., 0.)&


                                              & ], [3, 3] )   !< Klüppelberg PhD thesis 4.28
  

  integer, parameter, dimension(3)    :: dirvecx = [1, 0, 0]
  integer, parameter, dimension(3)    :: dirvecy = [0, 1, 0]
  integer, parameter, dimension(3)    :: dirvecz = [0, 0, 1]
  integer, parameter, dimension(3, 3) :: id3x3   = reshape([dirvecx, dirvecy, dirvecz], [3, 3])

  real,private,parameter:: sqrt32=sqrt(3.0/2.0)
  complex, dimension(3,3,5) :: mat2ord  !initialization below due to problems with nvhpc
  
  DATA mat2ord / &
       ! Tensor 1
       (sqrt32, 0.0),  (0.0, sqrt32),  (0.0, 0.0), &
       (0.0, sqrt32),  (-sqrt32, 0.0), (0.0, 0.0), &
       (0.0, 0.0),            (0.0, 0.0),            (0.0, 0.0), &
       ! Tensor 2
       (0.0, 0.0),            (0.0, 0.0),            (sqrt32, 0.0), &
       (0.0, 0.0),            (0.0, 0.0),            (0.0, sqrt32), &
       (sqrt32, 0.0),  (0.0, sqrt32),  (0.0, 0.0), &
       ! Tensor 3
       (-1.0, 0.0),           (0.0, 0.0),            (0.0, 0.0), &
       (0.0, 0.0),            (-1.0, 0.0),           (0.0, 0.0), &
       (0.0, 0.0),            (0.0, 0.0),            (2.0, 0.0), &
       ! Tensor 4
       (0.0, 0.0),            (0.0, 0.0),            (-sqrt32, 0.0), &
       (0.0, 0.0),            (0.0, 0.0),            (0.0, sqrt32), &
       (-sqrt32, 0.0), (0.0, sqrt32),  (0.0, 0.0), &
       ! Tensor 5
       (sqrt32, 0.0),  (0.0, -sqrt32), (0.0, 0.0), &
       (0.0, -sqrt32), (-sqrt32, 0.0), (0.0, 0.0), &
       (0.0, 0.0),            (0.0, 0.0),            (0.0, 0.0) /
                                                  
CONTAINS

  REAL PURE FUNCTION pimach()
    IMPLICIT NONE
    !  This subprogram supplies the value of the constant PI correct to
    !  machine precision where

    !  PI=3.1415926535897932384626433832795028841971693993751058209749446

    pimach = 3.1415926535897932
  END FUNCTION pimach

  REAL ELEMENTAL FUNCTION c_light(fac)
    IMPLICIT NONE
    !  This subprogram supplies the value of c according to
    !  NIST standard 13.1.99
    !  Hartree and Rydbergs changed by fac = 1.0 or 2.0

    REAL, INTENT (IN) :: fac
    c_light = 137.0359895e0 * fac * warp_factor
    !c_light = 1e6*fac
  END FUNCTION c_light

  SUBROUTINE  read_xml_constants(this,xml)
    USE m_types_xml
    CLASS(t_constants),INTENT(INout)::this
    TYPE(t_xml),INTENT(inout)   ::xml

    IF (xml%GetNumberOfNodes('/fleurInput/calculationSetup/expertModes/@warp_factor')==1)&
    warp_factor=evaluateFirstOnly(xml%GetAttributeValue('/fleurInput/calculationSetup/expertModes/@warp_factor'))
  END SUBROUTINE read_xml_constants
  SUBROUTINE mpi_bc_constants(this,mpi_comm,irank)
    USE m_mpi_bc_tool
    CLASS(t_constants),INTENT(INOUT)::this
    INTEGER,INTENT(IN):: mpi_comm
    INTEGER,INTENT(IN),OPTIONAL::irank
    INTEGER ::rank
    IF (PRESENT(irank)) THEN
       rank=irank
    ELSE
       rank=0
    END IF
    CALL mpi_bc(warp_factor,rank,mpi_comm)
  END SUBROUTINE mpi_bc_constants
END MODULE m_constants
