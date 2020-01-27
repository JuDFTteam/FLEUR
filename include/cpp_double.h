#ifdef CPP_DOUBLE
# define CPP_LAPACK_DOUBLE
# define CPP_BLAS_DOUBLE
#endif

#ifdef CPP_LAPACK_DOUBLE
# define CPP_LAPACK_cgesv  zgesv
# define CPP_LAPACK_chetrf zhetrf
# define CPP_LAPACK_chetri zhetri
# define CPP_LAPACK_sgeev  dgeev
# define CPP_LAPACK_sgetrf dgetrf
# define CPP_LAPACK_sgetri dgetri
# define CPP_LAPACK_chegvd zhegvd
# define CPP_LAPACK_chegvx zhegvx
# define CPP_LAPACK_cheevx zheevx
# define CPP_LAPACK_cgetrs zgetrs
# define CPP_LAPACK_ssygvx dsygvx
# define CPP_LAPACK_ctrtrs ztrtrs
# define CPP_LAPACK_cpotrf zpotrf
# define CPP_LAPACK_spotrf dpotrf
# define CPP_LAPACK_ssygst dsygst
# define CPP_LAPACK_ssyevr dsyevr
# define CPP_LAPACK_strtrs dtrtrs
# define CPP_LAPACK_chegst zhegst
# define CPP_LAPACK_cheevr zheevr
# define CPP_LAPACK_slamch dlamch
# define CPP_LAPACK_spptrf dpptrf
# define CPP_LAPACK_sspevx dspevx
# define CPP_LAPACK_sspev  dspev 
# define CPP_LAPACK_sspgst dspgst
# define CPP_LAPACK_ssygv  dsygv
# define CPP_LAPACK_chegv  zhegv
# define CPP_LAPACK_stptrs dtptrs
# define CPP_LAPACK_chpevx zhpevx
# define CPP_LAPACK_cheev  zheev
# define CPP_LAPACK_ssyev  dsyev
# define CPP_LAPACK_cheevd zheevd
# define CPP_LAPACK_chpgst zhpgst
# define CPP_LAPACK_chpgv  zhpgv
# define CPP_LAPACK_cpptrf zpptrf
# define CPP_LAPACK_ctptrs ztptrs
# define CPP_LAPACK_cgeevx zgeevx
# define CPP_LAPACK_cgehrd zgehrd
# define CPP_LAPACK_cunghr zunghr
# define CPP_LAPACK_chseqr zhseqr
# define CPP_LAPACK_chsein zhsein
# define CPP_LAPACK_sgesv  dgesv
# define CPP_LAPACK_sgesv  dgesv
# define CPP_LAPACK_sygvd  dygvd
#else
# define CPP_LAPACK_chegvd chegvd
# define CPP_LAPACK_cheevx cheevx
# define CPP_LAPACK_chegvx chegvx
# define CPP_LAPACK_cgetrs cgetrs
# define CPP_LAPACK_ssygvx ssygvx
# define CPP_LAPACK_ctrtrs ctrtrs
# define CPP_LAPACK_cpotrf cpotrf
# define CPP_LAPACK_spotrf spotrf
# define CPP_LAPACK_ssygst ssygst
# define CPP_LAPACK_ssyevr ssyevr
# define CPP_LAPACK_strtrs strtrs
# define CPP_LAPACK_chegst chegst
# define CPP_LAPACK_cheevr cheevr
# define CPP_LAPACK_slamch slamch
# define CPP_LAPACK_spptrf spptrf
# define CPP_LAPACK_sspevx sspevx
# define CPP_LAPACK_sspev  sspev 
# define CPP_LAPACK_sspgst sspgst
# define CPP_LAPACK_ssygv  ssygv
# define CPP_LAPACK_chegv  chegv
# define CPP_LAPACK_stptrs stptrs
# define CPP_LAPACK_chpevx chpevx
# define CPP_LAPACK_cheev  cheev
# define CPP_LAPACK_ssyev  ssyev
# define CPP_LAPACK_cheevd cheevd
# define CPP_LAPACK_chpgst chpgst
# define CPP_LAPACK_chpgv  chpgv
# define CPP_LAPACK_cpptrf cpptrf
# define CPP_LAPACK_ctptrs ctptrs
# define CPP_LAPACK_cgeevx cgeevx
# define CPP_LAPACK_cgehrd cgehrd
# define CPP_LAPACK_cunghr cunghr
# define CPP_LAPACK_chseqr chseqr
# define CPP_LAPACK_chsein chsein
# define CPP_LAPACK_sgesv  sgesv
# define CPP_LAPACK_sgesv  sgesv
#endif

#ifdef CPP_BLAS_DOUBLE
# define CPP_BLAS_cher   zher
# define CPP_BLAS_chpr   zhpr
# define CPP_BLAS_sgemm  dgemm
# define CPP_BLAS_cgemm  zgemm
# define CPP_BLAS_sasum	 dasum
# define CPP_BLAS_saxpy	 daxpy
# define CPP_BLAS_scnrm2 dznrm2
# define CPP_BLAS_scopy	 dcopy
# define CPP_BLAS_sdot	 ddot
# define CPP_BLAS_sgemv	 dgemv
# define CPP_BLAS_sspmv  dspmv
# define CPP_BLAS_sscal  dscal
# define CPP_BLAS_cscal  zscal
# define CPP_BLAS_ssum	 dsum
# define CPP_BLAS_caxpy	 zaxpy
# define CPP_BLAS_ccopy	 zcopy
# define CPP_BLAS_cdotc	 zdotc
# define CPP_BLAS_cdotu  zdotu
# define CPP_BLAS_cgemv	 zgemv
# define CPP_BLAS_chpmv  zhpmv
# define CPP_BLAS_cscal	 zscal
# define CPP_BLAS_csscal zdscal

#else
# define CPP_BLAS_cher   cher
# define CPP_BLAS_chpr   chpr
# define CPP_BLAS_sgemm  sgemm
# define CPP_BLAS_cgemm  cgemm
# define CPP_BLAS_sasum	 sasum
# define CPP_BLAS_saxpy	 saxpy
# define CPP_BLAS_scnrm2 scnrm2
# define CPP_BLAS_scopy	 scopy
# define CPP_BLAS_sdot	 sdot
# define CPP_BLAS_sgemv	 sgemv
# define CPP_BLAS_sspmv  sspmv
# define CPP_BLAS_sscal  sscal
# define CPP_BLAS_cscal  cscal
# define CPP_BLAS_ssum	 ssum
# define CPP_BLAS_caxpy	 caxpy
# define CPP_BLAS_ccopy	 ccopy
# define CPP_BLAS_cdotc	 cdotc
# define CPP_BLAS_cdotu  cdotu
# define CPP_BLAS_cgemv	 cgemv
# define CPP_BLAS_chpmv  chpmv
# define CPP_BLAS_cscal	 cscal
# define CPP_BLAS_csscal csscal
#endif

#ifdef CPP_DOUBLE
# define CPP_LAPACK_cgetri zgetri
# define CPP_LAPACK_cgetrf zgetrf
# define CPP_LAPACK_cgeev  zgeev
# define CPP_LAPACK_ssytri dsytri
# define CPP_LAPACK_ssytrf dsytrf
# define CPP_LAPACK_pdsygvx pdsygvx
# define CPP_LAPACK_pzhegvx pzhegvx
# define CPP_MPI_REAL      MPI_DOUBLE_PRECISION
# define CPP_MPI_COMPLEX   MPI_DOUBLE_COMPLEX
# define CPP_MPI_TYP_REAL    MPI_DOUBLE_PRECISION
# define CPP_MPI_TYP_COMPLEX MPI_DOUBLE_COMPLEX
#else
# define CPP_LAPACK_cgetri cgetri
# define CPP_LAPACK_cgetrf cgetrf
# define CPP_LAPACK_cgeev  cgeev
# define CPP_LAPACK_ssytri ssytri
# define CPP_LAPACK_ssytrf ssytrf
# define CPP_LAPACK_pdsygvx pssygvx
# define CPP_LAPACK_pzhegvx pchegvx
# define CPP_MPI_REAL      MPI_REAL
# define CPP_MPI_COMPLEX   MPI_COMPLEX
# define CPP_MPI_TYP_REAL    MPI_REAL
# define CPP_MPI_TYP_COMPLEX MPI_COMPLEX
#endif

