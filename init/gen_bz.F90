      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! gen_bz generates the (whole) Brillouin zone from the          !
      ! (irreducible) k-points given in the kpts file.                !
      !                                                               !
      !                                     M.Betzinger (09/07)       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
      MODULE m_gen_bz

      CONTAINS

      SUBROUTINE gen_bz( kpts,sym)

      !     bk       ::    irreducible k-points
      !     nkpt     ::    number of irr. k-points
      !     bkf    ::    all k-points
      !     nkptf      ::    number of all k-points
      !     bkp       ::    k-point parent
      !     bksym     ::    symmetry operation, that connects the parent
      !                     k-point with the current one

      USE m_util  ,ONLY: modulo1
      USE m_types

      IMPLICIT NONE
      TYPE(t_kpts),INTENT(INOUT) :: kpts
      TYPE(t_sym),INTENT(IN)     :: sym

!     - local scalars -
      INTEGER                 ::  ic,iop,ikpt,ikpt1
      REAL                    ::  rdum
      LOGICAL                 ::  ldum
      
!     - local arrays - 
      INTEGER,ALLOCATABLE     ::  iarr(:)
      REAL                    ::  wt(kpts%nkpt)
      REAL                    ::  rrot(3,3,sym%nsym),rotkpt(3)
      REAL,ALLOCATABLE        ::  rarr1(:,:)

      ALLOCATE( kpts%bkf(3,sym%nsym*kpts%nkpt),kpts%bkp(sym%nsym*kpts%nkpt),kpts%bksym(sym%nsym*kpts%nkpt) )
      
      ! symmetrie operation in reciprocal space
      DO iop=1,sym%nsym
        IF( iop .le. sym%nop ) THEN
          rrot(:,:,iop) = transpose( sym%mrot(:,:,sym%invtab(iop)) )
        ELSE
          rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
        END IF
      END DO


      ! read in number of kpts
      OPEN (41,file='kpts',form='formatted',status='old')
      READ (41,*) ikpt, rdum ! jump over first line
      DO ikpt=1,kpts%nkpt
        READ (41,*) rdum,rdum,rdum,wt(ikpt)
      END DO
      CLOSE (41)
      kpts%nkptf = nint( sum( wt ) )
#ifdef CPP_DEBUG
      WRITE(*,*) 'whole BZ consists of',kpts%nkptf,' k-points'
#endif


      !apply symmetrie operations to all k-points of IBZ
      kpts%bkf = 0
      kpts%bkf(:,1) = 0
      kpts%bkp(1)      = 1
      kpts%bksym(1)    = 1
      ic          = 1
      
      DO iop=1,sym%nsym
        DO ikpt=2,kpts%nkpt
          ldum = .false.
          rotkpt = matmul ( rrot(:,:,iop), kpts%bk(:,ikpt) )
          !transform back into IBZ
          rotkpt = modulo1( rotkpt,kpts%nkpt3 )
          DO ikpt1=1,ic
            IF (maxval( abs( kpts%bkf(:,ikpt1) - rotkpt )) .le. 1E-06)THEN
              ldum = .true.
              EXIT
            END IF
          END DO
          
          IF( .not. ldum ) THEN
            ic = ic + 1
            kpts%bkf(:,ic) = rotkpt
            kpts%bkp(ic)      = ikpt
            kpts%bksym(ic)    = iop
          END IF
        END DO
      END DO
      
      IF ( kpts%nkptf /= ic ) STOP 'gen_bz: error kpts/symmetrie '

      !reallocate bkf,bkp,bksym
      ALLOCATE (iarr(kpts%nkptf))
      iarr = kpts%bkp(:kpts%nkptf)
      DEALLOCATE( kpts%bkp )
      ALLOCATE ( kpts%bkp(kpts%nkptf) )
      kpts%bkp = iarr
      iarr= kpts%bksym(:kpts%nkptf)
      DEALLOCATE ( kpts%bksym )
      ALLOCATE ( kpts%bksym(kpts%nkptf) )
      kpts%bksym = iarr
      DEALLOCATE( iarr )
      ALLOCATE ( rarr1(3,kpts%nkptf) )
      rarr1 = kpts%bkf(:,:kpts%nkptf)
      DEALLOCATE ( kpts%bkf )
      ALLOCATE ( kpts%bkf(3,kpts%nkptf) )
      kpts%bkf = rarr1
      
     

#ifdef CPP_DEBUG
       !test
      DO ikpt=2,kpts%nkptf
        IF ( kpts%bkp(ikpt) .ne. ikpt) THEN
          rotkpt = matmul( rrot(:,:,kpts%bksym(ikpt)),kpts%bkf(:,kpts%bkp(ikpt)))
          rotkpt = modulo1(rotkpt,kpts%nkpt3)
          IF( maxval( abs( rotkpt - kpts%bkf(:,ikpt) ) ) .gt. 1e-8) THEN
            STOP 'gen_bz: failure kpts%bksym,bkp'
          END IF
        END IF
      END DO

      WRITE(*,*)'whole BZ',kpts%nkptf
      DO ikpt=1,kpts%nkptf
        WRITE(*,*) kpts%bkf(:,ikpt)
      END DO
#endif
      END SUBROUTINE gen_bz

      END MODULE m_gen_bz
