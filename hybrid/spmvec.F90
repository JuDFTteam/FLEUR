      MODULE m_spmvec
      
      CONTAINS
       
      SUBROUTINE spmvec(&
     &           atoms,lcutm,maxlcutm,nindxm,maxindxm,&
     &           nbasm,nbasp,hybrid,ikpt,kpts,cell,&
     &           coulomb_mt1,coulomb_mt2,coulomb_mt3,&
     &           coulomb_mtir,vecin,&
     &           vecout)
      
      USE m_wrapper
      USE m_constants
      USE m_types
      IMPLICIT NONE
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      
      ! - scalars -
      INTEGER, INTENT(IN) ::  ikpt 
      INTEGER, INTENT(IN) ::  nbasm,nbasp 
      INTEGER, INTENT(IN) ::  maxindxm,maxlcutm
            
      ! - arrays -
      INTEGER, INTENT(IN) ::  lcutm(atoms%ntype)
      INTEGER, INTENT(IN) ::  nindxm(0:maxlcutm,atoms%ntype)
      
      
      REAL   , INTENT(IN) ::  coulomb_mt1(maxindxm-1,maxindxm-1,&
     &                                    0:maxlcutm,atoms%ntype)
#ifdef CPP_INVERSION      
      REAL   , INTENT(IN) ::  coulomb_mt2(maxindxm-1,-maxlcutm:maxlcutm,&
     &                                    0:maxlcutm+1,atoms%nat)
      REAL   , INTENT(IN) ::  coulomb_mt3(maxindxm-1,atoms%nat,nat)
#ifdef CPP_IRCOULOMBAPPROX
      REAL   , INTENT(IN) ::  coulomb_mtir(:,:)
#else
      REAL   , INTENT(IN) ::  coulomb_mtir(:)
#endif
      REAL   , INTENT(IN) ::  vecin (nbasm)
      REAL   , INTENT(OUT)::  vecout(nbasm)
#else
      COMPLEX, INTENT(IN) ::  coulomb_mt2(maxindxm-1,-maxlcutm:maxlcutm,&
     &                                    0:maxlcutm+1,atoms%nat)
      COMPLEX, INTENT(IN) ::  coulomb_mt3(maxindxm-1,atoms%nat,atoms%nat)
#ifdef CPP_IRCOULOMBAPPROX
      COMPLEX, INTENT(IN) ::  coulomb_mtir(:,:)
#else
      COMPLEX, INTENT(IN) ::  coulomb_mtir(:)
#endif
      COMPLEX, INTENT(IN) ::  vecin (nbasm)
      COMPLEX, INTENT(OUT)::  vecout(nbasm)
#endif
      
      ! - local scalars -
      INTEGER             ::  itype,ieq,iatom,ishift
      INTEGER             ::  itype1,ieq1,iatom1,ishift1
      INTEGER             ::  indx0,indx1,indx2,indx3,indx4
      INTEGER             ::  i,ibasm,igptm
      INTEGER             ::  l
      INTEGER             ::  n,m
      
      REAL                ::  gnorm
      
      COMPLEX, PARAMETER  ::  img = (0d0,1d0)
      ! - local arrays -
      
#ifdef CPP_INVERSION
      REAL                ::  vecinhlp (nbasm)
      REAL   ,ALLOCATABLE ::  coulhlp(:,:)
#else
      REAL                ::  vecr(maxindxm-1),veci(maxindxm-1)
      COMPLEX             ::  vecinhlp (nbasm)
      COMPLEX,ALLOCATABLE ::  coulhlp(:,:)
#endif
 
  
      
      vecinhlp = vecin
      
      CALL reorder(nbasm,nbasp,atoms,lcutm, maxlcutm,nindxm,1, vecinhlp)
  

      ibasm = 0
      iatom = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l = 0,lcutm(itype)
            DO m = -l,l
              ibasm = ibasm + nindxm(l,itype) - 1
            END DO
          END DO
        END DO
      END DO
      
      
      ! compute vecout for the indices from 0:ibasm
      iatom = 0
      indx1 = 0; indx2 = 0; indx3 = ibasm
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l = 0,lcutm(itype)
            DO m = -l,l
              indx1 = indx1 + 1
              indx2 = indx2 + nindxm(l,itype) - 1
              indx3 = indx3 + 1

              
              vecout(indx1:indx2) = matmul(coulomb_mt1(:nindxm(l,itype)-1,:nindxm(l,itype)-1, l,itype),&
     &                vecinhlp(indx1:indx2))

! #ifdef CPP_INVERSION
!               CALL DGEMV('N',nindxm(l,itype)-1,nindxm(l,itype)-1,1d0,coulomb_mt1(:,:,l,itype),maxindxm-1,
!      &                    vecinhlp(indx1:),1,0d0,vecout(indx1:),1)
! #else
!               CALL DSYMV('L',nindxm(l,itype)-1,1d0,coulomb_mt1(:,:,l,itype),maxindxm-1,REAL(vecinhlp(indx1:)),1,0d0,vecr,1)
!               CALL DSYMV('L',nindxm(l,itype)-1,1d0,coulomb_mt1(:,:,l,itype),maxindxm-1,IMAG(vecinhlp(indx1:)),1,0d0,veci,1)
!              
!               vecout(indx1:indx2) = vecr(1:nindxm(l,itype)-1) + img * veci(1:nindxm(l,itype)-1) 
! #endif  
              
              
              vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:nindxm(l,itype)-1,m,l,iatom) * vecinhlp(indx3) 
            
              indx1 = indx2
            END DO

            
          END DO
        END DO
      END DO
     
      IF( indx2 .ne. ibasm ) &
     &  STOP 'spmvec: error counting basis functions'
      
      IF( ikpt .eq. 1) THEN       
        iatom = 0
        indx0 = 0
        DO itype = 1,atoms%ntype
          ishift = sum( (/ ((2*l+1)*(nindxm(l,itype)-1), l=0,lcutm(itype) )/) )
          DO ieq = 1,atoms%neq(itype)
            iatom = iatom + 1
            l = 0
            m = 0
            
            indx1 = indx0 + 1 
            indx2 = indx1 + nindxm(l,itype) - 2
            
            iatom1 = 0
            indx3  = ibasm
            DO itype1 = 1,atoms%ntype
              ishift1 = (lcutm(itype1)+1)**2
              DO ieq1 = 1,atoms%neq(itype1)
                iatom1 = iatom1 + 1
                indx4  = indx3 + (ieq1 - 1)*ishift1 + 1
                IF( iatom .eq. iatom1 ) CYCLE
                
                vecout(indx1:indx2) = vecout(indx1:indx2)&
     &                              + coulomb_mt3(:nindxm(l,itype)-1, iatom1,iatom) * vecinhlp(indx4)

              END DO
              indx3 = indx3 + atoms%neq(itype1)*ishift1
            END DO
            
            IF( indx3 .ne. nbasp ) STOP 'spmvec: error counting index indx3'
     
            vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:nindxm(l,itype)-1,0, maxlcutm+1,iatom) * vecinhlp(indx3+1)
     
            indx0 = indx0 + ishift
          END DO

        END DO
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm
      
#ifdef CPP_IRCOULOMBAPPROX
        
      indx0 = sum( (/ ( ((2*l+1)*atoms%neq(itype),l=0,lcutm(itype)), itype=1,atoms%ntype ) /) )+ hybrid%ngptm
      indx1 = sum( (/ ( ((2*l+1)*atoms%neq(itype),l=0,lcutm(itype)), itype=1,atoms%ntype ) /) )

#ifdef CPP_INVERSION
      CALL DGEMV('N',indx1,indx0,1d0,coulomb_mtir,(maxlcutm+1)**2*atoms,&
     &          vecinhlp(ibasm+1:),1,0d0,vecout(ibasm+1:),1)
      
      CALL DGEMV('T',indx1,hybrid,1d0,coulomb_mtir(:indx1,indx1+1:),&
     &          indx1,vecinhlp(ibasm+1:),1,0d0,vecout(ibasm+indx1+1:),1)
#else    
      CALL ZGEMV('N',indx1,indx0,(1d0,0d0),coulomb_mtir,&
     &          (maxlcutm+1)**2*atoms,vecinhlp(ibasm+1:),&
     &          1,(0d0,0d0),vecout(ibasm+1:),1)
      
      CALL ZGEMV('C',indx1,hybrid,(1d0,0d0),coulomb_mtir(:indx1,indx1+1:)&
     &          ,indx1,vecinhlp(ibasm+1:),1,(0d0,0d0),&
     &          vecout(ibasm+indx1+1:),1)
#endif

!       vecout(ibasm+1:ibasm+indx1) = matmul( coulomb_mtir(:indx1,:indx0),vecinhlp(ibasm+1:ibasm+indx0) )
!       vecout(ibasm+indx1+1:ibasm+indx0) = matmul( conjg(transpose(coulomb_mtir(:indx1,indx1+1:indx0))),
!      &                                            vecinhlp(ibasm+1:ibasm+indx1) )

      
      indx0 = ibasm + indx1
      IF( indx0 .ne. nbasp ) STOP 'spmvec: error indx0'
      DO i = 1,hybrid%ngptm 
        indx0 = indx0 + 1
        igptm = hybrid%pgptm(i)
        gnorm = sqrt(sum(matmul(kpts%bk(:) + hybrid%gptm(:,igptm),cell%bmat)**2))
        IF( gnorm .eq. 0 ) CYCLE
        vecout(indx0) = vecout(indx0) + fpi*vecinhlp(indx0)/gnorm
      END DO
      
#else

      indx1 = sum( (/ ( ((2*l+1)*atoms%neq(itype),l=0,lcutm(itype)),&
     &                                      itype=1,atoms%ntype ) /) )+ hybrid%ngptm(ikpt)
#ifdef CPP_INVERSION
      CALL dspmv('U',indx1,1d0,coulomb_mtir,vecinhlp(ibasm+1),1,0d0, vecout(ibasm+1),1)
#else
      call zhpmv('U',indx1,(1d0,0d0),coulomb_mtir,vecinhlp(ibasm+1), 1,(0d0,0d0),vecout(ibasm+1),1)
#endif

#endif
      iatom = 0
      indx1 = ibasm; indx2 = 0; indx3 = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO l = 0,lcutm(itype)
            n = nindxm(l,itype)
            DO m = -l,l
              indx1 = indx1 + 1
              indx2 = indx2 + 1
              indx3 = indx3 + n - 1
              
              vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt2(:n-1,m,l,iatom), vecinhlp(indx2:indx3))
              indx2 = indx3
            END DO
             
          END DO
        END DO
      END DO
      
      
      IF( ikpt .eq. 1 ) THEN
        iatom = 0
        indx0 = 0
        DO itype = 1,atoms%ntype
          ishift = sum( (/ ((2*l+1)*(nindxm(l,itype)-1), l=0,lcutm(itype) )/) )
          DO ieq = 1,atoms%neq(itype)
            iatom = iatom + 1
            indx1 = indx0 + 1
            indx2 = indx1 + nindxm(0,itype) - 2
            vecout(nbasp+1) = vecout(nbasp+1) + dotprod(coulomb_mt2(:nindxm(0,itype)-1,0, maxlcutm+1,iatom), vecinhlp(indx1:indx2))
        
           indx0 = indx0 + ishift
          END DO
        END DO
         
        iatom = 0
        indx0 = ibasm
        DO itype = 1,atoms%ntype
          ishift = (lcutm(itype)+1)**2
          DO ieq = 1,atoms%neq(itype)
            iatom = iatom + 1
            indx1 = indx0 + 1
            
            iatom1 = 0
            indx2  = 0
            DO itype1 = 1,atoms%ntype
              ishift1 = sum( (/ ( (2*l+1)*(nindxm(l,itype1)-1), l=0,lcutm(itype1)) /))
              DO ieq1 = 1,atoms%neq(itype1)
                iatom1 = iatom1 + 1
                IF( iatom1 .eq. iatom ) CYCLE
                
                indx3  = indx2 + (ieq1 - 1)*ishift1 + 1
                indx4  = indx3 + nindxm(0,itype1) - 2

                vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt3(:nindxm(0,itype1)-1,iatom,iatom1), vecinhlp(indx3:indx4))

              END DO
              indx2 = indx2 + atoms%neq(itype1)*ishift1
            END DO
            indx0 = indx0 + ishift
          END DO
        END DO
        IF( indx0 .ne. nbasp ) STOP 'spmvec: error index counting (indx0)'
      END IF
      
      CALL reorder(nbasm,nbasp,atoms,lcutm,&
     &             maxlcutm,nindxm,2,&
     &             vecout)
      

      END SUBROUTINE spmvec
      
      
      SUBROUTINE reorder(nbasm,nbasp,atoms,lcutm, maxlcutm,nindxm,imode, vec)
        USE m_types
      IMPLICIT NONE
      TYPE(t_atoms),INTENT(IN)   :: atoms
      
      ! - scalars -
      INTEGER, INTENT(IN)   ::  maxlcutm
      INTEGER, INTENT(IN)   ::  nbasm,nbasp
      INTEGER, INTENT(IN)   ::  imode
      
        
      ! - arrays -
      INTEGER, INTENT(IN)   ::  lcutm(atoms%ntype)
      INTEGER, INTENT(IN)   ::  nindxm(0:maxlcutm,atoms%ntype)
#ifdef CPP_INVERSION
      REAL   , INTENT(INOUT)::  vec(nbasm)
#else
      COMPLEX, INTENT(INOUT)::  vec(nbasm)
#endif
      ! - local scalars -
      INTEGER               ::  itype,ieq
      INTEGER               ::  indx1,indx2
      INTEGER               ::  l
      INTEGER               ::  n,m
      ! - local arrays -
#ifdef CPP_INVERSION
      REAL                  ::  vechlp(nbasm)
#else
      COMPLEX               ::  vechlp(nbasm)
#endif


      IF( imode .ne. 1 .and. imode .ne. 2 ) STOP 'reorder: imode equals neither 1 nor 2'

      vechlp = vec
      
      IF( imode .eq. 1 ) THEN
        indx1 = 0
        indx2 = 0
        DO itype = 1,atoms%ntype
          DO ieq = 1,atoms%neq(itype)
            DO l = 0,lcutm(itype)
              DO m = -l,l
                DO n = 1,nindxm(l,itype)-1
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  vec(indx1) = vechlp(indx2)
                END DO
                indx2 = indx2 + 1
              END DO
            END DO
          END DO
        END DO
      
        indx2 = 0
        DO itype = 1,atoms%ntype
          DO ieq = 1,atoms%neq(itype)
            DO l = 0,lcutm(itype)
              DO m = -l,l
                indx1 = indx1 + 1
                indx2 = indx2 + nindxm(l,itype)
                vec(indx1) = vechlp(indx2)
              END DO
            END DO
          END DO
        END DO        
      ELSE IF (imode .eq. 2) THEN
        indx1 = 0
        indx2 = 0
        DO itype = 1,atoms%ntype
          DO ieq = 1,atoms%neq(itype)
            DO l = 0,lcutm(itype)
              DO m = -l,l
                DO n = 1,nindxm(l,itype)-1
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  vec(indx2) = vechlp(indx1)
                END DO
                indx2 = indx2 + 1
              END DO
            END DO
          END DO
        END DO
        
        indx2 = 0
        DO itype = 1,atoms%ntype
          DO ieq = 1,atoms%neq(itype)
            DO l = 0,lcutm(itype)
              DO m = -l,l
                indx1 = indx1 + 1
                indx2 = indx2 + nindxm(l,itype)
                vec(indx2) = vechlp(indx1)
              END DO
            END DO
          END DO
        END DO        
      END IF
      !IR must not be rearranged

      END SUBROUTINE reorder
      
      END MODULE
