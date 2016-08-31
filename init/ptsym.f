      MODULE m_ptsym
      use m_juDFT
!********************************************************************
!     determines the point group symmetry for each representative
!     atom and then check whether there are several with the same
!     local symmetry.
!
!     input:  nops      number of operations in space group
!             mrot      rotation matrices in INTERNAL coordinates
!             tau       non-primitive translations, in INTERNAL coord.
!
!             ntype     number of atom types
!             neq       number of equivalent atoms of each type
!             pos       atomic positions in INTERNAL coord.
!
!     output: nsymt     number of symmetry kinds
!             typsym    symmetry kind for each atom type
!             nrot      number of operations for each symmetry kind
!             locops    mapping of operations to space group list
!*********************************************************************
      CONTAINS
      SUBROUTINE ptsym(
     >                 ntype,natd,neq,pos,nops,mrot,tau,lmax,
     <                 nsymt,typsym,nrot,locops)

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: ntype,neq(ntype),natd
      INTEGER, INTENT (IN) :: nops,mrot(3,3,nops),lmax(ntype)
      REAL,    INTENT (IN) :: tau(3,nops),pos(3,natd)

      INTEGER, INTENT(OUT) :: locops(nops,natd),nrot(natd)
      INTEGER, INTENT(OUT) :: nsymt,typsym(natd)

      REAL, PARAMETER :: eps=1.e-7

      INTEGER :: iop,irot,n,na,nn,nsym
      INTEGER :: indsym(natd),indsym1(natd)
      REAL    :: v(3),sv(3)

!--->    loop over representative atoms
      na = 1
      DO n = 1,ntype
          v(:) = pos(:,na)

!--->     loop over space group operations to see which belong
!--->     to point group: sv = {R|t_R}v - v = Rv + t_R - v
          iop = 0
          DO irot = 1,nops
            sv = matmul( real(mrot(:,:,irot)) , v ) + tau(:,irot) - v
!--->       check whether sv is a lattice vector ( sv integer)
            IF ( ANY( ABS( sv - ANINT(sv) ) > eps ) ) CYCLE

!--->       this operation belongs to the point group
            iop = iop + 1
            locops(iop,na) = irot
          ENDDO

          nrot(na) = iop
          na = na + neq(n)
      ENDDO

!--->    check that the number of operations in local groups are correct
      na = 1
      DO n = 1, ntype
          IF ( neq(n)*nrot(na) .NE. nops ) THEN
            WRITE (6,'(/a,i3)') ' symmetry is incorrect for atom',na
            WRITE (6,'(" neq=",i3,", nrot=",i3,", nops=",i3)')            
     &               neq(n),nrot(na),nops
            CALL juDFT_error("symmetry is incorrect for some atomp"
     +           ,calledby ="ptsym")
          ENDIF
          na = na + neq(n)
      ENDDO

!--->    now determine unique symmetry kinds

      nsymt     = 1
      typsym(1) = 1
      indsym(1) = 1
      indsym1(1)= 1

      na = 1
      atom_loop: DO n = 1, ntype
        IF (na > 1) THEN 

          symm_loop: DO nsym = 1, nsymt
            IF ( nrot(na) .NE. nrot(indsym(nsym)) ) CYCLE

            DO irot=1,nrot(na)
               IF(locops(irot,na).NE.locops(irot,indsym(nsym))) THEN
                  CYCLE symm_loop  ! try next symmetry type
               ENDIF
            ENDDO
            IF ( lmax(n).NE.lmax(indsym1(nsym)) ) CYCLE

!--->    same symmetry as a previous one:
            typsym(na) = nsym
            na = na + 1
            equi : DO nn = 2, neq(n)
               typsym(na) = nsym
               na = na + 1
            ENDDO equi
            CYCLE atom_loop   ! go to next atom

          ENDDO symm_loop

!--->       new symmetry kind
          nsymt = nsymt + 1
          typsym(na) = nsymt
          indsym(nsymt) = na
          indsym1(nsymt) = n

        ENDIF
        na = na + 1
        equi_loop : DO nn = 2, neq(n) 
           typsym(na) = nsymt
           na = na + 1
        ENDDO equi_loop

      ENDDO atom_loop

!--->    pack locops array
      DO n = 2, nsymt
         nrot(n) = nrot(indsym(n))
         DO irot = 1,nrot(n)
            locops(irot,n) = locops(irot,indsym(n))
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ptsym
      END MODULE m_ptsym
