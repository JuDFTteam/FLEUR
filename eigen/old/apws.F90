MODULE m_apws
  use m_juDFT
  !*********************************************************************
  !     determines the lapw list such that |k+G|<rkmax.
  !     bk(i) is the nk k-point given in internal (i.e. b1,b2,b3) units.
  !        m. weinert  1986
  !     unit 29 removed gb 2004
  !*********************************************************************
  !     modified for explicit use of z-reflection symmetry in seclr4.f
  !        g. bihlmayer '96
  !     subroutine boxdim added to treat non-orthogonal lattice vectors
  !        s.bluegel, IFF, 18.Nov.97
  !*********************************************************************
CONTAINS

  !REMOVED and REPLACED by types_lapw

END MODULE m_apws
