    SUBROUTINE kedmt(atoms,we,ne,eigVecCoeffs,denCoeffs,ispin,&
                    &input,vTot,enpara,fmpi,usdus,hub1data)

    USE m_genMTBasis
    USE m_types
    USE m_radfun
    USE m_constants
    use m_types_xcpot
    USE m_calcDenCoeffs

    IMPLICIT NONE

    INTEGER,              INTENT(IN)    :: ne, ispin
    TYPE(t_eigVecCoeffs), INTENT(INOUT) :: eigVecCoeffs ! Should be renamed
    REAL,                 INTENT(INOUT)    :: we(:)!(nobd)
    TYPE(t_atoms),        INTENT(INOUT)    :: atoms
    TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

    TYPE(t_potden),        INTENT(INOUT)    :: vTot
    TYPE(t_enpara),        INTENT(INOUT)    :: enpara
    TYPE(t_mpi),           INTENT(INOUT)    :: fmpi
    TYPE(t_input),         INTENT(INOUT)    :: input
    TYPE (t_usdus),        INTENT(INOUT)    :: usdus
    TYPE(t_hub1data),       OPTIONAL, INTENT(INOUT) :: hub1data

    INTEGER i,l,lm ,n,na,natom,m,imesh,itype

    REAL,    ALLOCATABLE  :: f(:,:,:,:),g(:,:,:,:),flo(:,:,:,:) ! radial functions

    ! temporary stored
    REAL,           ALLOCATABLE              :: df(:)
    REAL,           ALLOCATABLE              :: dg(:)
 
    ! should be stored in type
    REAL,           ALLOCATABLE              :: dfr(:,:,:,:,:)
    REAL,           ALLOCATABLE              :: dgr(:,:,:,:,:)



    ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd,input%jspins)) 
    ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd,input%jspins))
    ALLOCATE (flo(atoms%jmtd,2,atoms%nlod,input%jspins)) 

    ALLOCATE( df(atoms%jmtd) )
    ALLOCATE( dg(atoms%jmtd) )
    df=0
    dg=0
    
    ALLOCATE( dfr(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,2) )
    ALLOCATE( dgr(atoms%jmtd,2,0:atoms%lmaxd,input%jspins,2) )
    dfr=0
    dgr=0

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

                ! F term
                !eigvec%kexsbvof -> kxeigveccoef%abcof(*,*,0..15,...)

                eigVecCoeffs%kexabcof(i,lm,0,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,0,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


                eigVecCoeffs%kexabcof(i,lm,1,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,1,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


                eigVecCoeffs%kexabcof(i,lm,2,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,2,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


                eigVecCoeffs%kexabcof(i,lm,3,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,3,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


                !   F dot term

                eigVecCoeffs%kexabcof(i,lm,4,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,4,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,5,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,5,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,6,natom,ispin)=&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,6,natom,ispin)=&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,7,natom,ispin)=&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,7,natom,ispin)=&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                !   G term

                eigVecCoeffs%kexabcof(i,lm,8,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,8,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,9,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,9,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


                eigVecCoeffs%kexabcof(i,lm,10,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,10,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,11,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,11,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

                ! G dot term

                eigVecCoeffs%kexabcof(i,lm,12,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,12,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,13,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,13,natom,ispin)=(-l)*&
      &            eigVecCoeffs%KEDU(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,14,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,14,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kexabcof(i,lm,15,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)
                eigVecCoeffs%keyabcof(i,lm,15,natom,ispin)=(l+1)*&
      &            eigVecCoeffs%KEDD(lm,-1)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)


!     --------- z component


      !    J term

                eigVecCoeffs%kezabcof(i,lm,0,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

                eigVecCoeffs%kezabcof(i,lm,1,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)


      !    J dot term

                eigVecCoeffs%kezabcof(i,lm,2,natom,ispin)=&
      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kezabcof(i,lm,3,natom,ispin)=&
      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)



      !    K term

                eigVecCoeffs%kezabcof(i,lm,4,natom,ispin)=(-l)*&
      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)

                eigVecCoeffs%kezabcof(i,lm,5,natom,ispin)=(1+l)*&
      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,0,natom,ispin)



      !    K dot term

                eigVecCoeffs%kezabcof(i,lm,6,natom,ispin)=(-l)*&
      &           eigVecCoeffs%KEDU(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)

                eigVecCoeffs%kezabcof(i,lm,7,natom,ispin)=(1+l)*&
      &           eigVecCoeffs%KEDD(lm,0)*eigVecCoeffs%abcof(i,lm,1,natom,ispin)



                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO


    DO iType = 1, atoms%ntype

      CALL genMTBasis(atoms,enpara,vTot,fmpi,iType,ispin,usdus,&
              &f(:,:,0:,ispin),g(:,:,0:,ispin),flo(:,:,:,ispin),&
              &hub1data=hub1data)

    ENDDO

     natom = 0
     DO n = 1,atoms%ntype
        itype=n
      DO na = 1,atoms%neq(n)
       natom = natom + 1
       DO l = 0,atoms%lmax(n)
  
       call derivative_loc_ked( f(:,1,l,ispin), itype, atoms, df )
       dfr(:,1,l,ispin,1)=df(:)
       call derivative_loc_ked( f(:,2,l,ispin), itype, atoms, df )
       dfr(:,2,l,ispin,1)=df(:)
       call derivative_loc_ked( g(:,1,l,ispin), itype, atoms, dg )
       dgr(:,1,l,ispin,1)=dg(:)
       call derivative_loc_ked( g(:,2,l,ispin), itype, atoms, dg )
       dgr(:,2,l,ispin,1)=dg(:)


       DO imesh = 1, atoms%jmtd ! not atoms%jri(itype)

        dfr(imesh,1,l,ispin,2)=&
               &f(imesh,1,l,ispin)/atoms%rmsh(imesh, itype)
        dfr(imesh,2,l,ispin,2)&
               &=f(imesh,2,l,ispin)/atoms%rmsh(imesh, itype)
        dgr(imesh,1,l,ispin,2)=&
               &g(imesh,1,l,ispin)/atoms%rmsh(imesh, itype)
        dgr(imesh,2,l,ispin,2)=&
               &g(imesh,2,l,ispin)/atoms%rmsh(imesh, itype)


        !WRITE(1994,"(4I5,F15.9)") ispin,natom,l,imesh,dfr(imesh,1,l,ispin,1)
        !WRITE(1995,"(4I5,F15.9)") ispin,natom,l,imesh,dfr(imesh,1,l,ispin,2)

       ENDDO

       ENDDO
      ENDDO
     ENDDO

     DEALLOCATE(f)
     DEALLOCATE(g)
     DEALLOCATE(flo)
     DEALLOCATE(df)
     DEALLOCATE(dg)


    ! atoms%lmaxd

    END SUBROUTINE


    SUBROUTINE derivative_loc_ked(f, itype, atoms, df)
    
    USE m_types


    integer,       intent(in)  :: itype
    type(t_atoms), intent(inout)  :: atoms                          
    real,          intent(inout)  :: f(atoms%jri(itype))
                                    
    real,          intent(out) :: df(atoms%jri(itype))
    real                       :: h, r, d21, d32, d43, d31, d42, d41, df1, df2, s
    real                       :: y0, y1, y2
    integer                    :: i, n

    n = atoms%jri(itype)
    h = atoms%dx(itype)
    r = atoms%rmsh(1, itype)

      ! use Lagrange interpolation of 3rd order (and averaging) for points 3 to n
    d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)    
    d31 = d21 + d32      ; d42 = d32 + d43    
    d41 = d31 + d43
                  
    df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1d0/d21 - 1d0/d31 - 1d0/d41) * f(1)&
            &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)

    df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1d0/d21 - 1d0/d32 - 1d0/d42) * f(2)&     
            &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)

    df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +&
            &  ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(3) + d31*d32 / (d41*d42*d43)*f(4)

      do i = 3, n - 2
      d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
      d31 = d42 ; d42 = d42 * exp(h)
      d41 = d41 * exp(h)
      df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1d0/d21 - 1d0/d32 - 1d0/d42) * f(i) + &
              &d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) *f(i+2)
              
      df(i) = ( df1 + df2 ) / 2
      
      df1   = d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42)* f(i) +&
              &(1d0/d31 + 1d0/d32 - 1d0/d43)*f(i+1)+d31*d32 / (d41*d42*d43) * f(i+2)
      enddo

      df(n-1) = df1
      df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -&      
              &d41*d42 / (d31*d32*d43) * f(n-1) + ( 1d0/d41 + 1d0/d42 + 1d0/d43 ) * f(n)
               
      ! for first two points use Lagrange interpolation of second order for log(f(i))
      ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
      s = sign(1d0,f(1))

      if(sign(1d0,f(2)) /= s .or. sign(1d0,f(3))  /= s .or. any(abs(f(:3)) < 1e0)) then
                       
      d21   = r * (exp(h)-1)
      d32   = d21 * exp(h)
      d31   = d21 + d32
      s     = df(3) / (d31*d32) - f(1) / (d21*d31**2) + f(2) / (d21*d32**2) - f(3) /&
      &(d31**2*d32) - f(3) / (d31*d32**2)

      df(1) = - (d21+d31) / (d21*d31) * f(1) + d31 / (d21*d32) * f(2) - d21 /&
      &(d31*d32) * f(3) + d21*d31 * s
      
      df(2) = - (d21-d32) / (d21*d32) * f(2) - d32 / (d21*d31) * f(1)&
      &+ d21 / (d31*d32) * f(3) - d21*d32 * s

      else

         y0    = log(abs(f(1)))         
         y1    = log(abs(f(2)))                           
         y2    = log(abs(f(3)))                           
         df(1) = ( - 3*y0/2 + 2*y1 - y2/2 ) * f(1) / (h*r)                                    
         df(2) = (y2-y0)/2* f(2) / (h*r*exp(h))
      endif

    END SUBROUTINE        
