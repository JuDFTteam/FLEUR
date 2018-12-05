!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_kpts
  INTEGER,PARAMETER:: kpts_by_number=1
  INTEGER,PARAMETER:: kpts_by_mesh  =2
  INTEGER,PARAMETER:: kpts_by_list  =3
  
  TYPE t_kpts
     INTEGER :: specificationType 
                                  
     !no
     INTEGER               :: nkpt
     INTEGER               :: ntet
     REAL                  :: posScale
     LOGICAL               :: l_gamma
     !(3,nkpt) k-vectors internal units
     REAL,ALLOCATABLE      :: bk(:,:)
     !(nkpts) weights
     REAL,ALLOCATABLE      :: wtkpt(:)
     INTEGER               :: nkptf !<k-vectors in full BZ
     INTEGER               :: nkpt3(3)
     REAL                  :: kPointDensity(3) ! only used if k point set is defined as density
     REAL   ,ALLOCATABLE   :: bkf(:,:)
     INTEGER,ALLOCATABLE   :: bkp(:)
     INTEGER,ALLOCATABLE   :: bksym(:)
     INTEGER                       :: numSpecialPoints
     INTEGER, ALLOCATABLE          :: specialPointIndices(:)
     CHARACTER(LEN=50),ALLOCATABLE :: specialPointNames(:)
     REAL   ,ALLOCATABLE           :: specialPoints(:,:)
     INTEGER,ALLOCATABLE           :: ntetra(:,:)
     REAL   ,ALLOCATABLE           :: voltet(:)
     REAL   ,ALLOCATABLE           :: sc_list(:,:) !list for all information about folding of bandstructure (need for unfoldBandKPTS)((k(x,y,z),K(x,y,z),m(g1,g2,g3)),(nkpt),k_original(x,y,z))
  ENDTYPE t_kpts


  SUBROUTINE read_xml_kpts()

      kpts%ntet = 1
      kpts%numSpecialPoints = 1
      ALLOCATE(kpts%specialPoints(3,kpts%numSpecialPoints))
    END SUBROUTINE read_xml_kpts

  SUBROUTINE init_kpts()


     IF ( banddos%dos .AND. banddos%ndir == -3 ) THEN
          WRITE(*,*) 'Recalculating k point grid to cover the full BZ.'
          CALL gen_bz(kpts,sym)
          kpts%nkpt = kpts%nkptf
          DEALLOCATE(kpts%bk,kpts%wtkpt)
          ALLOCATE(kpts%bk(3,kpts%nkptf),kpts%wtkpt(kpts%nkptf))
          kpts%bk(:,:) = kpts%bkf(:,:)
          IF (kpts%nkpt3(1)*kpts%nkpt3(2)*kpts%nkpt3(3).NE.kpts%nkptf) THEN
             IF(kpts%l_gamma) THEN
                kpts%wtkpt = 1.0 / (kpts%nkptf-1)
                DO i = 1, kpts%nkptf
                   IF(ALL(kpts%bk(:,i).EQ.0.0)) kpts%wtkpt(i) = 0.0
                END DO
             ELSE
                CALL juDFT_error("nkptf does not match product of nkpt3(i).",calledby="fleur_init")
             END IF
          ELSE
             kpts%wtkpt = 1.0 / kpts%nkptf
          END IF
       END IF
     END SUBROUTINE init_kpts
END MODULE m_types_kpts
