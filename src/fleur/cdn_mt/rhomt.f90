MODULE m_rhomt
CONTAINS
  SUBROUTINE rhomt(atoms,we,ne,eigVecCoeffs,denCoeffs,ispin)
    !     ***************************************************************
    !     perform the sum over m (for each l) and bands to set up the
    !     coefficient of spherical charge densities in subroutine
    !     cdnval                                   c.l.fu
    !     ***************************************************************

    USE m_types
    IMPLICIT NONE

    INTEGER,              INTENT(IN)    :: ne, ispin
    TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
    REAL,                 INTENT(IN)    :: we(:)!(nobd)
    TYPE(t_atoms),        INTENT(IN)    :: atoms
    TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

    INTEGER i,l,lm ,n,na,natom,m

    natom = 0
    DO n = 1,atoms%ntype
       DO na = 1,atoms%neq(n)
          natom = natom + 1
          DO l = 0,atoms%lmax(n)
             !     -----> sum over m
             DO m = -l,l
                lm = l* (l+1) + m
                !     -----> sum over occupied bands
                DO i = 1,ne
                   denCoeffs%uu(l,n,ispin) = denCoeffs%uu(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%abcof(i,lm,0,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,0,natom,ispin)))
                   denCoeffs%dd(l,n,ispin) = denCoeffs%dd(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%abcof(i,lm,1,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ispin)))
                   denCoeffs%du(l,n,ispin) = denCoeffs%du(l,n,ispin) +&
                      we(i) * REAL(eigVecCoeffs%abcof(i,lm,0,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ispin)))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

!   YAMASHITA_ADD   

    natom = 0
    DO n = 1,atoms%ntype
       DO na = 1,atoms%neq(n)
          natom = natom + 1
          DO l = 0,atoms%lmax(n)
             !     -----> sum over m
             DO m = -l,l
                lm = l* (l+1) + m

!     -----> sum over occupied bands
                DO i = 1,ne 

       !            denCoeffs%uu(l,n,ispin) = denCoeffs%uu(l,n,ispin) +&
       !               we(i) * REAL(eigVecCoeffs%abcof(i,lm,0,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,0,natom,ispin)))
       !            denCoeffs%dd(l,n,ispin) = denCoeffs%dd(l,n,ispin) +&
       !               we(i) * REAL(eigVecCoeffs%abcof(i,lm,1,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ispin)))
       !            denCoeffs%du(l,n,ispin) = denCoeffs%du(l,n,ispin) +&
       !               we(i) * REAL(eigVecCoeffs%abcof(i,lm,0,natom,ispin)*CONJG(eigVecCoeffs%abcof(i,lm,1,natom,ispin)))

!     --------- x and y components

                ! F term

!                eigVecCoeffs%kexyabcof(i,l,m,0,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


!                eigVecCoeffs%kexyabcof(i,l,m,1,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


!                eigVecCoeffs%kexyabcof(i,l,m,2,natom,ispin)=&
!      &           eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


!                eigVecCoeffs%kexyabcof(i,l,m,3,natom,ispin)=&
!      &           eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)                  
                !   F dot term

!                eigVecCoeffs%kexyabcof(i,l,m,4,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)                  

!                eigVecCoeffs%kexyabcof(i,l,m,5,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)                   

!                eigVecCoeffs%kexyabcof(i,l,m,6,natom,ispin)=&
!      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)                  

!                eigVecCoeffs%kexyabcof(i,l,m,7,natom,ispin)=&
!      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)                     


                !   G term

!                eigVecCoeffs%kexyabcof(i,l,m,8,natom,ispin)=(-l)*&
!      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)              

!                eigVecCoeffs%kexyabcof(i,l,m,9,natom,ispin)=(-l)*&
!      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


!                eigVecCoeffs%kexyabcof(i,l,m,10,natom,ispin)=(l+1)*&
!      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)                  

!                eigVecCoeffs%kexyabcof(i,l,m,11,natom,ispin)=(l+1)*&
!      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)                  

                ! G dot term

!                eigVecCoeffs%kexyabcof(i,l,m,12,natom,ispin)=(-l)*&
!      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin) 

!                eigVecCoeffs%kexyabcof(i,l,m,13,natom,ispin)=(-l)*&
!      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)                  
        
!                eigVecCoeffs%kexyabcof(i,l,m,14,natom,ispin)=(l+1)*&
!      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)     
        
!                eigVecCoeffs%kexyabcof(i,l,m,15,natom,ispin)=(l+1)*&
!      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)  

!     --------- z component


      !    J term

!                eigVecCoeffs%kezabcof(i,l,m,0,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

!                eigVecCoeffs%kezabcof(i,l,m,1,natom,ispin)=&
!      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


      !    J dot term

!                eigVecCoeffs%kezabcof(i,l,m,2,natom,ispin)=&
!      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

!                eigVecCoeffs%kezabcof(i,l,m,3,natom,ispin)=&
!      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
    


      !    K term

!                eigVecCoeffs%kezabcof(i,l,m,4,natom,ispin)=(-l)*&
!      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

!                eigVecCoeffs%kezabcof(i,l,m,5,natom,ispin)=(1+l)*&
!      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)



      !    K dot term 

!                eigVecCoeffs%kezabcof(i,l,m,6,natom,ispin)=(-l)*&
!      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

!                eigVecCoeffs%kezabcof(i,l,m,7,natom,ispin)=(1+l)*&
!      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin) 

                

                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE rhomt
END MODULE m_rhomt
