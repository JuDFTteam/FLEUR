      MODULE m_constants
      IMPLICIT NONE

      REAL,PARAMETER:: pi_const=3.1415926535897932
      REAL,PARAMETER:: tpi_const=2.*3.1415926535897932
      REAL,PARAMETER:: fpi_const=4.*3.1415926535897932
      REAL,PARAMETER:: sfp_const=sqrt(4.*3.1415926535897932)
      REAL,PARAMETER:: hartree_to_ev_const=27.21138505
      CHARACTER(2),DIMENSION(0:103),PARAMETER:: namat_const=(/
     &      'va',' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',         
     &     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     &     ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     &     'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     &     'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',
     &     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     &     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',
     &     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw'/)

      CONTAINS
!------------------------------------------------------------------------
      REAL PURE FUNCTION pimach()
!
!     This subprogram supplies the value of the constant PI correct to
!     machine precision where
!
!     PI=3.1415926535897932384626433832795028841971693993751058209749446
!
      pimach = 3.1415926535897932
      END FUNCTION pimach
!------------------------------------------------------------------------
      REAL ELEMENTAL FUNCTION c_light(fac)
!
!     This subprogram supplies the value of c according to
!     NIST standard 13.1.99 
!     Hartree and Rydbergs changed by fac = 1.0 or 2.0
!
      REAL, INTENT (IN) :: fac
      c_light = 137.0359895e0 * fac 

      RETURN
      END FUNCTION c_light
      END MODULE m_constants

