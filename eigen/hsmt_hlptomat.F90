module m_hsmt_hlptomat
#include "juDFT_env.h"
    implicit none
    contains
    subroutine hsmt_hlptomat(nlotot,nv,sub_comm,chi11,chi21,chi22,aahlp,aa,bbhlp,bb)
    !hsmt_hlptomat: aa/bbhlp - to -aa/bb matrix
    !Rotate the aahlp&bbhlp arrays from the local spin-frame into the global frame
    !and add the data to the aa&bb arrays, call mingeselle in distributed case

#ifdef CPP_MPI
      USE m_mingeselle
#endif
        implicit none
        integer, intent(in)                :: nlotot,nv(:),sub_comm
        complex, intent(in)                :: chi11,chi21,chi22
        complex, allocatable,intent(inout) :: aahlp(:)
        complex,             intent(inout) :: aa(:)
        complex, optional,intent(inout)    :: bb(:),bbhlp(:)

        integer :: ii,ij,ki,kj,n_rank,n_size

        REAL :: aa_r(1),bb_r(1) !dummy arguments for mingeselle
#ifdef CPP_MPI
#include "mpif.h"
        CALL MPI_COMM_RANK(sub_comm,n_rank,ki)
        CALL MPI_COMM_SIZE(sub_comm,n_size,ki)
#else
        n_size=1
#endif
        IF (n_size==1) THEN
            DO ki =  1, nv(1)+nlotot
                !--->        spin-up spin-up part
                ii = (ki-1)*(ki)/2
                ij = (ki-1)*(ki)/2
                aa(ij+1:ij+ki)=aa(ij+1:ij+ki)+chi11*aahlp(ii+1:ii+ki)
                if (present(bb)) bb(ij+1:ij+ki)=bb(ij+1:ij+ki)+chi11*bbhlp(ii+1:ii+ki)
                !--->        spin-down spin-down part
                ij = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2+nv(1)+nlotot
                aa(ij+1:ij+ki)=aa(ij+1:ij+ki)+chi22*aahlp(ii+1:ii+ki)
                if (present(bb)) bb(ij+1:ij+ki)=bb(ij+1:ij+ki)+chi22*bbhlp(ii+1:ii+ki)
                !--->        spin-down spin-up part, lower triangle
                ij = (nv(1)+nlotot+ki-1)*(nv(1)+nlotot+ki)/2
                aa(ij+1:ij+ki)=aa(ij+1:ij+ki)+chi21*aahlp(ii+1:ii+ki)
                if (present(bb)) bb(ij+1:ij+ki)=bb(ij+1:ij+ki)+chi21*bbhlp(ii+1:ii+ki)
                !--->        spin-down spin-up part, upper triangle.
                DO kj = 1,ki-1
                    ij = (nv(1)+nlotot+kj-1)*(nv(1)+nlotot+kj)/2 + ki
                    aa(ij) = aa(ij) + conjg(aahlp(ii+kj))*chi21
                    if (present(bb)) bb(ij) = bb(ij) + conjg(bbhlp(ii+kj))*chi21
                ENDDO
            ENDDO
        ELSE
            aa(:size(aahlp)) = aa(:size(aahlp))+aahlp*chi11
            aahlp = conjg(aahlp)*chi21
            IF (present(bb).and.nlotot>1) THEN
              !CALL juDFT_error("noco+LO and EVP is broken")
              bb(:size(aahlp)) = bb(:size(aahlp))+bbhlp*chi11
              bbhlp = conjg(bbhlp)*chi21
            ENDIF
            CALL mingeselle(SUB_COMM,n_size,n_rank,nv,&
                aahlp,.false.,aa_r,aa)
            IF (present(bb).and.nlotot>1) CALL mingeselle(SUB_COMM,n_size,n_rank,nv,&
                bbhlp,.false.,bb_r,bb)
        ENDIF

    end subroutine

    
end module m_hsmt_hlptomat
