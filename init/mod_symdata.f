      MODULE m_symdata
!-------------------------------------------------------------------------+
! Contains all arrays for the determination of the 2D symmetry elements   !
! from a given name of plane-group                                        !
!                                                                   gb`02 !
!-------------------------------------------------------------------------+
      
      INTEGER ord2(25)    ! Number of 2D symmetry operations
      LOGICAL l_c2(25)    ! whether plane group contains the c_2
      REAL    tau2(2,3)   ! translations for the generators
      INTEGER gen2(2,2,9) ! rotation matrices for the generators
      INTEGER spg2(3,25)  ! generators for 2d space groups
      INTEGER gnt2(3,25)  ! translations for 2d space groups
      CHARACTER(len=4) :: namgr2(25) ! names of 2d space groups
      CHARACTER(len=4) :: nammap(20) ! names as in the inp-file
!
! All 20 names that can be given in the inp-file
!
      DATA nammap/'p1  ','p2  ','pmy ','pgy ','cmy ','pmm ','pmg ',
     +            'pgg ','cmm ','p4  ','p4m ','p4g ','p3  ','p3m1',
     +            'p31m','p6  ','p6m ','pm  ','pg  ','cm  '/
!
! Number of 2D symmetry operations for all 17+8 plane groups
!
      DATA ord2/1,2,2,2,2,4,4,4,4,4,8,8,3,6,6,6,12,2,2,2,3,6,6,6,12/
!
! Whether the plane group contains the c_2 or not:
!
      DATA l_c2/.false.,.true.,3*.false.,7*.true.,3*.false.,
     +           2*.true.,6*.false.,2*.true./
!
! Names for  all 17+8 plane groups as given in the International tables
!
      DATA namgr2/'p1  ','p2  ','pm  ','pg  ','cm  ','p2mm','p2mg',
     +            'p2gg','cmm ','p4  ','p4mm','p4gm','p3  ','p3m1',
     +            'p31m','p6  ','p6mm','pm  ','pg  ','cm  ','p3  ',
     +            'p3m1','p31m','p6  ','p6mm'/
!
! Generators: c_3 is for hx3 (angle = 120), c_3' for hex (angle = 60)
!             m_d and m_t are 'diagonal' mirror planes for p3m1 and p31m
!
      DATA gen2/ 1, 0, 0, 1,  -1, 0, 0,-1, -1, 0, 0, 1,   !   E, c_2, m_x
     +           0,-1, 1, 0,   0, 1,-1,-1,  0,-1,-1, 0,   ! c_4, c_3, m_d
     +           0, 1, 1, 0,   1, 0, 0,-1,  0,-1, 1,-1/   ! m_t, m_y, c_3'
!
! and translation vectors for the generators:
!
      DATA tau2/ 0.5, 0.0,  0.0, 0.5,  0.5, 0.5/
!
! Generators for the groups (excluding the identity):
!
      DATA spg2/ 0, 0, 0,  2, 0, 0,  3, 0, 0,  3, 0, 0,   ! p1,  p2,   pm,   pg
     +           6, 0, 0,  2, 3, 0,  2, 3, 0,  2, 3, 0,   ! cm,  p2mm, p2mg, p2gg
     +           2, 6, 0,  2, 4, 0,  2, 3, 4,  2, 3, 4,   ! cmm, p4,   p4mm, p4gm
     +           5, 0, 0,  5, 6, 0,  5, 7, 0,  2, 5, 0,   ! p3,  p3m1, p31m, p6
     +           2, 5, 6,  8, 0, 0,  8, 0, 0,  7, 0, 0,   ! p6mm,pm,   pg,   cm
     +           9, 0, 0,  9, 7, 0,  9, 6, 0,  2, 9, 0,   ! p3,  p3m1, p31m, p6
     +           2, 9, 6/
!
! translation vectors for non-symmorphic groups :
!
      DATA gnt2/ 0, 0, 0,  0, 0, 0,  0, 0, 0,  2, 0, 0,   ! pg
     +           0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 3, 0,   ! p2mg, p2gg
     +           0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 3, 0,   ! p4gm
     +           0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0,   !
     +           0, 0, 0,  0, 0, 0,  1, 0, 0,  0, 0, 0,   ! pg
     +           0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0,   !
     +           0, 0, 0/
!
      END MODULE m_symdata
