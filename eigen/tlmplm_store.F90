!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_tlmplm_store
! Instead of storing data to tmat&tmas files this module stores into module variables
! used to transfer the results from tlmplm&density matrix in case of lda+u from eigen
! into force_a21
!      D.W 2014
    USE m_types
    IMPLICIT NONE
    PRIVATE
    TYPE(t_tlmplm)            :: td_stored
    COMPLEX,ALLOCATABLE       :: vs_mmp_stored(:,:,:,:)
    PUBLIC write_tlmplm, read_tlmplm, read_tlmplm_vs_mmp
CONTAINS
    SUBROUTINE write_tlmplm(td,vs_mmp,ldau,ispin,jspin,jspins)
        TYPE(t_tlmplm),INTENT(IN) :: td
        COMPLEX,INTENT(IN)        :: vs_mmp(:,:,:,:)
        LOGICAL,INTENT(IN)        :: ldau
        INTEGER,INTENT(IN)        :: ispin,jspin,jspins

        INTEGER:: llod,lmplmd,ntypd,lmd,mlot_d,mlolot_d

        IF (.NOT.allocated(td_stored%tuu)) THEN
            lmplmd=size(td%tuu,1)-1
            ntypd=size(td%tuu,2)
            ALLOCATE(td_stored%tuu(0:lmplmd,ntypd,jspins))
            ALLOCATE(td_stored%tud(0:lmplmd,ntypd,jspins))
            ALLOCATE(td_stored%tdd(0:lmplmd,ntypd,jspins))
            ALLOCATE(td_stored%tdu(0:lmplmd,ntypd,jspins))
            lmd=size(td%tdulo,1)-1
            llod=size(td%tdulo,2)!not actual llod but -llod:llod
            mlot_d = size(td%tdulo,3)
            mlolot_d = size(td%tuloulo,3)
            ALLOCATE(td_stored%tdulo(0:lmd,llod,mlot_d,jspins))
            ALLOCATE(td_stored%tuulo(0:lmd,llod,mlot_d,jspins))
            ALLOCATE(td_stored%tuloulo(llod,llod,mlolot_d,jspins))
            ALLOCATE(td_stored%ind(0:lmd,0:lmd,ntypd,jspins))
            IF (ldau) &
                ALLOCATE(vs_mmp_stored(  &
                size(vs_mmp,1),size(vs_mmp,2),size(vs_mmp,3),jspins))
        ENDIF

        td_stored%tuu(:,:,jspin) = td%tuu(:,:,ispin)
        td_stored%tud(:,:,jspin) = td%tud(:,:,ispin)
        td_stored%tdu(:,:,jspin) = td%tdu(:,:,ispin)
        td_stored%tdd(:,:,jspin) = td%tdd(:,:,ispin)

        td_stored%tdulo(:,:,:,jspin)   = td%tdulo(:,:,:,ispin)
        td_stored%tuulo(:,:,:,jspin)   = td%tuulo(:,:,:,ispin)
        td_stored%tuloulo(:,:,:,jspin) = td%tuloulo(:,:,:,ispin)
        td_stored%ind(:,:,:,jspin)     = td%ind(:,:,:,ispin)
        IF (ldau) vs_mmp_stored(:,:,:,jspin)=vs_mmp(:,:,:,ispin)


    END SUBROUTINE write_tlmplm

    SUBROUTINE read_tlmplm(n,jspin,nlo,tuu,tud,tdu,tdd,ind,tuulo,tuloulo,tdulo)
        COMPLEX,INTENT(OUT)::tuu(:),tdd(:),tud(:),tdu(:)
        INTEGER,INTENT(OUT)::ind(:,:)
        COMPLEX,INTENT(OUT)::tuulo(:,:,:),tdulo(:,:,:),tuloulo(:,:,:)
        INTEGER,INTENT(IN) :: n,jspin,nlo(:)

        INTEGER:: mlo,mlolo
        tuu=td_stored%tuu(:size(tuu,1),n,jspin)
        tud=td_stored%tud(:size(tuu,1),n,jspin)
        tdu=td_stored%tdu(:size(tuu,1),n,jspin)
        tdd=td_stored%tdd(:size(tuu,1),n,jspin)
        ind=td_stored%ind(:size(ind,1),:size(ind,2),n,jspin)

        IF (nlo(n)>0) THEN
            mlo=sum(nlo(:n-1))+1
            mlolo=dot_product(nlo(:n-1),nlo(:n-1)+1)/2+1
            tuulo(:,:,mlo:(mlo+nlo(n)-1))=&
               td_stored%tuulo(:size(tuulo,1),:size(tuulo,2),mlo:mlo+nlo(n)-1,jspin)
            tdulo(:,:,mlo:(mlo+nlo(n)-1))=&
               td_stored%tdulo(:size(tuulo,1),:size(tuulo,2),mlo:mlo+nlo(n)-1,jspin)
            tuloulo(:,:,mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1)=&
               td_stored%tuloulo(:size(tuloulo,1),:size(tuloulo,2),mlolo:mlolo+nlo(n)*(nlo(n)+1)/2-1,jspin)
        ENDIF
    
    END SUBROUTINE read_tlmplm

    SUBROUTINE read_tlmplm_vs_mmp(jspin,n_u,vs_mmp)

       INTEGER, INTENT(IN)  :: jspin, n_u
       COMPLEX, INTENT(OUT) :: vs_mmp(:,:,:)

       IF(n_u.GT.0) THEN
          vs_mmp(:,:,:) = vs_mmp_stored(:,:,:,jspin)
       END IF

    END SUBROUTINE read_tlmplm_vs_mmp

END MODULE m_tlmplm_store
