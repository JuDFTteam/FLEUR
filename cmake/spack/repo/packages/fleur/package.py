# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Fleur(CudaPackage, Package):
    """FLEUR (Full-potential Linearised augmented plane wave in EURope)
    is a code family for calculating groundstate as well as excited-state properties
    of solids within the context of density functional theory (DFT)."""

    homepage = "https://www.flapw.de/current"
    git = "https://iffgit.fz-juelich.de/fleur/fleur.git"

    license("MIT")
    phases = ["configure","build", "install"]
    version("develop", branch="develop")
    version("release", branch="release")
    version("7.0", commit="4db9bf591c4f6e01580d0671a7bf9ed4a566468b")
    version("5.1", tag="MaX-R5.1", commit="a482abd9511b16412c2222e2ac1b1a303acd454b", deprecated=True)
    version("5.0", tag="MaX-R5", commit="f2df362c3dad6ef39938807ea14e4ec4cb677723", deprecated=True)
    version("4.0", tag="MaX-R4", commit="ea0db7877451e6240124e960c5546318c9ab3953", deprecated=True)
    version("3.1", tag="MaX-R3.1", commit="f6288a0699604ad9e11efbfcde824b96db429404", deprecated=True)

    patch("elsi-config.patch",when="@7.0")

    variant("mpi", default=True, description="Enable MPI support")
    variant("hdf5", default=True, description="Enable HDF5 support")
    variant("scalapack", default=False, description="Enable SCALAPACK")
    variant("fft",default="internal",values=("internal", "mkl", "fftw"),description="Enable the use of Intel MKL FFT/FFTW provider")
    variant("elpa", default=False, description="Enable ELPA support")
    variant("elsi", default=False, description="Enable ELSI support")
    variant("magma", default=False, description="Enable Magma support")
    variant("libxc", default=False, description="Enable libxc support")
    variant("spfft", default=False, description="Enable spfft support")
    variant("wannier90", default=False, description="Enable wannier90 support")
    variant("openmp", default=True, description="Enable OpenMP support")
    variant("amd", default=False, description="Use some patch for AMD processors")
    variant("cuda",default=False,description="Use OpenACC on top of CUDA for NVIDIA GPUs")
    variant("cuda_arch",default=80 ,description="specify the CUDA architecture to build for")
    variant("build_type",default="RelWithDebInfo",description="The build type to build",values=("Debug", "Release", "RelWithDebInfo"))
    
    depends_on("cmake", type="build")
    depends_on("python@3:", type="build")
    depends_on("blas")
    depends_on("lapack")
    depends_on("libxml2")
    depends_on("mpi", when="+mpi")
    depends_on("intel-mkl", when="fft=mkl")
    depends_on("fftw-api", when="fft=fftw")
    depends_on("scalapack", when="+scalapack")
    depends_on("libxc", when="+libxc")
    depends_on("hdf5+hl+fortran", when="+hdf5")
    depends_on("magma+fortran", when="+magma")
    depends_on("wannier90", when="+wannier90")
    depends_on("spfft+fortran~openmp", when="+spfft~openmp")
    depends_on("spfft+fortran+openmp", when="+spfft+openmp")
    depends_on("elpa~openmp", when="+elpa~openmp")
    depends_on("elpa+openmp", when="+elpa+openmp")
    depends_on("elsi", when="+elsi")
    depends_on("cuda",when="+cuda")
    requires("%nvhpc",when="+cuda",msg="OpenACC on CUDA only with %nvhpc")

    conflicts("%intel@:16.0.4", msg="ifort version <16.0 will most probably not work correctly")
    conflicts("%gcc@:6.3.0", msg="gfortran version <v6.3 will most probably not work correctly")
    conflicts("%pgi@:18.4.0",msg="You need at least PGI version 18.4 but might still run into some problems.")
    conflicts("~hdf5", when="@7.0:" ,msg="Spack installs of versions >6.0 need hdf5")
    conflicts("~scalapack", when="+elpa", msg="ELPA requires scalapack support")
    conflicts("~scalapack", when="+elsi", msg="ELSI requires scalapack support")
    conflicts("@:5.0", when="fft=fftw", msg="FFTW interface is supported from Fleur v5.0")
    conflicts("@:5.0", when="+wannier90", msg="wannier90 is supported from Fleur v5.0")
    conflicts("@:4.0", when="+spfft", msg="SpFFT is supported from Fleur v4.0")
    conflicts("cuda_arch=none", when="+cuda",msg="CUDA architecture is required")

    def setup_build_environment(self, env):
        spec = self.spec

        if "+mpi" in spec:
            env.set("CC", spec["mpi"].mpicc, force=True)
            env.set("FC", spec["mpi"].mpifc, force=True)
            env.set("CXX", spec["mpi"].mpicxx, force=True)

    def configure(self,spec,prefix):
        sh = which("bash")

        #populate the args to the FLEUR configure.sh
        args = []
        #Some args are collected first
        link_opt=[]
        lib_opt=[]
        include_opt=[]
        
        #Blas+Lapack+XML2 is always require
        link_opt.append(spec["blas"].libs.link_flags)
        lib_opt.append(spec["blas"].prefix.lib)
        include_opt.append(spec["blas"].prefix.include)

        link_opt.append(spec["lapack"].libs.link_flags)
        lib_opt.append(spec["lapack"].prefix.lib)
        include_opt.append(spec["lapack"].prefix.include)

        link_opt.append(spec["libxml2"].libs.link_flags)
        lib_opt.append(spec["libxml2"].prefix.lib)
        include_opt.append(spec["libxml2"].prefix.include)
        include_opt.append(join_path(spec["libxml2"].prefix.include, "libxml2"))


        if "+cuda" in spec:
            link_opt.append(spec["cuda"].libs.link_flags)
            lib_opt.append(spec["cuda"].prefix.lib)
            args.append("-gpu")
            cuda_arch = spec.variants["cuda_arch"].value
            args.append(f"acc:cc{cuda_arch}")
        if "fft=mkl" in spec:
            link_opt.append(spec["intel-mkl"].libs.link_flags)
            lib_opt.append(spec["intel-mkl"].prefix.lib)
            include_opt.append(spec["intel-mkl"].prefix.include)
        if "fft=fftw" in spec:
            link_opt.append(spec["fftw-api"].libs.link_flags)
            lib_opt.append(spec["fftw-api"].prefix.lib)
            include_opt.append(spec["fftw-api"].prefix.include)
        if "+scalapack" in spec:
            link_opt.append(spec["scalapack"].libs.link_flags)
            lib_opt.append(spec["scalapack"].prefix.lib)
        if "+hdf5" in spec:
            link_opt.append(spec["hdf5"].libs.link_flags)
            lib_opt.append(spec["hdf5"].prefix.lib)
            include_opt.append(spec["hdf5"].prefix.include)
            args.append("-hdf5")
            args.append("true")
        else:    
            args.append("-hdf5")
            args.append("false")
        if "+magma" in spec:
            link_opt.append(spec["magma"].libs.link_flags)
            lib_opt.append(spec["magma"].prefix.lib)
            include_opt.append(spec["magma"].prefix.include)
        if "+wannier90" in spec:
            # Workaround: The library is not called wannier90.a/so
            #    for this reason spec['wannier90'].libs.link_flags fails!
            link_opt.append("-lwannier")
            lib_opt.append(spec["wannier90"].prefix.lib)
        if "+spfft" in spec:
            link_opt.append(spec["spfft"].libs.link_flags)
            # Workaround: The library is installed in /lib64 not /lib
            lib_opt.append(spec["spfft"].prefix.lib + "64")
            # Workaround: The library needs spfft.mod in include/spfft path
            include_opt.append(join_path(spec["spfft"].prefix.include, "spfft"))
        if "+elsi" in spec:
            link_opt.append(spec["elsi"].libs.link_flags)
            #workaround: additional dependencies
            link_opt.append("-lMatrixSwitch -lNTPoly -lOMM -lelpa -lfortjson")
           # Workaround: The library is installed in /lib64 not /lib
            lib_opt.append(spec["elsi"].prefix.lib)
            # Workaround: The library needs spfft.mod in include/spfft path
            include_opt.append(spec["elsi"].prefix.include)
        if "+elpa" in spec:
            link_opt.append(spec["elpa"].libs.link_flags)
            lib_opt.append(spec["elpa"].prefix.lib)
            # Workaround: The library needs elpa.mod in include/elpa_%VERS/modules
            include_opt.append(spec["elpa"].prefix.include)
            include_opt.append(spec["elpa"].headers.include_flags[2:])
            include_opt.append(
                join_path(spec["elpa"].headers.include_flags[2:], "modules")
            )
        if "+amd" in spec:
            args.append("-amd")
        
        #Now add collected options
        args.append("-link")
        args.append(" ".join(link_opt))
        args.append("-libdir")
        args.append(" ".join(lib_opt))
        args.append("-includedir")
        args.append(" ".join(include_opt))
        
        sh("configure.sh", *args)

    def build(self,spec,prefix):
        with working_dir("build"):
            make()
        

    def install(self, spec, prefix):
        with working_dir("build"):
            mkdirp(prefix.bin)
            if "+mpi" in spec:
                install("fleur_MPI", prefix.bin)
            else:
                install("fleur", prefix.bin)
            install("inpgen", prefix.bin)

    @run_after("build")
    @on_package_attributes(run_tests=True)
    def test(self):
        with working_dir("build"):
            sh = which("bash")
            sh("run_tests.sh")
        