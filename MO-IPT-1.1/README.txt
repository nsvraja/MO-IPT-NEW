                      A JNCASR-LA-SiGMA Software Distribution
                                MO-IPT package 
              Version 1.1 (Revision Placeholder) - April 2019
        Copyright 2019 Jawaharlal Nehru Centre for Advanced Scientific Research

                                
                This package contains code and sample data for generating
                Green's functions and self-energies of multi-orbital strongly
                correlated electron systems within dynamical mean field theory.
                The impurity solver is multi-orbital iterated perturbation
		theory. Detailed installation and operation instructions and
		other files may be found in the docs subdirectory. The code
		employs parallelization through message passing interface.
		For the latest version and other resources visit GitHub.
 
     ---------------------------------------------------------------------
                                       Description

     MO-IPT stands for multi-orbital iterated perturbation theory. The present code
implements the MO-IPT solver as described in Dasari et al [Dasari, N., Mondal, W., Zhang, P. et al. Eur. Phys. J. B (2016) 89: 202. https://doi.org/10.1140/epjb/e2016-70133-4] within
the dynamical mean field theory (DMFT) framework and may be used to obtain single-
-particle spectra and self-energies for strongly correlated model Hamiltonians as
well as real materials with multiple orbital degrees of freedom. The code is implemented
for zero as well as finite temperatures. The impurity solver uses the second-order
self-energy in an ansatz motivated by the continued fraction expansion of the self-
-energy. The ansatz has certain free parameters that are chosen to satisfy high
frequency and the atomic limits. Since the ansatz reproduces low frequency Fermi liquid
behaviour and the band limit by construction, the MO-IPT is expected to be a reasonable
interpolating approximation. Naturally, the MO-IPT cannot be expected to be accurate
close to phase transitions, etc, where exact methods such as QMC and NRG would be
far more accurate. For extensive benchmarking of the method with continuous time
Monte Carlo and other methods, please see Dasari et al [paper reference]. One of the main
limitations of the code is that the Hund's coupling is presently implemented only as
a density-density interaction. The main merits of the MO-IPT is that it is fast,
numerically inexpensive, can deal with many orbitals, and provides real frequency
results at zero and finite temperature. The main motivation of our implementation
is to carry out first principles calculations of strongly correlated materials, so
the integration with band structure results (from e.g WIEN2K) is also implemented
in this set of codes.

The basic single-orbital IPT code was developed by N.S.Vidhyadhiraja (nsvraja@gmail.com).
The multi-orbital extension and MPI wrapper for k-summation were done by Nagamalleswararao
Dasari (nagamalleswararao.d@gmail.com). The optimization and benchmarks were carried out
by Dasari, Peng () and Wasim (wasimr.mondal@gmail.com) with the assistance of Mark
Jarrell (jarrellphysics@gmail.com), Juana Moreno (moreno@phys.lsu.edu) and N.S.Vidhyadhiraja. 
The version 1.1 implements all convolutions through fast Fourier transforms if a uniform grid
is used.

Detailed operating instructions, testing procedure and physics description of the test data
can be found in doc/manual.pdf. Sketchy details of the prerequisites, setup and operation
are given in docs/instructions.txt and a more detailed  manual is docs/manual.pdf.
The sample data given in 'data' sub-directory contain README files which describe the
specific problem being solved.

    ---------------------------------------------------------------------

                                      Prerequisites
                  
MO-IPT-1.1 has been tested on the following system:
 
(1) Ubuntu 18.04
    Distributor ID:	Ubuntu
    Description:	Ubuntu 18.04.2 LTS
    Release:	18.04
    Codename:	bionic
    Hardware - Dual Intel Xeon Gold 6148 @ 2.4GHz - Total 40 cores 192GB RAM 6TB HDD
    Compiler - 
      COLLECT_GCC=/usr/bin/gfortran
      COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-linux-gnu/7/lto-wrapper
      OFFLOAD_TARGET_NAMES=nvptx-none
      OFFLOAD_TARGET_DEFAULT=1
      Target: x86_64-linux-gnu
      Configured with: ../src/configure -v --with-pkgversion='Ubuntu 7.3.0-27ubuntu1~18.04' --with-bugurl=file:///usr/share/doc/gcc-7/README.Bugs --enable-languages=c,ada,c++,go,brig,d,fortran,objc,obj-c++ --prefix=/usr --with-gcc-major-version-only --program-suffix=-7 --program-prefix=x86_64-linux-gnu- --enable-shared --enable-linker-build-id --libexecdir=/usr/lib --without-included-gettext --enable-threads=posix --libdir=/usr/lib --enable-nls --with-sysroot=/ --enable-clocale=gnu --enable-libstdcxx-debug --enable-libstdcxx-time=yes --with-default-libstdcxx-abi=new --enable-gnu-unique-object --disable-vtable-verify --enable-libmpx --enable-plugin --enable-default-pie --with-system-zlib --with-target-system-zlib --enable-objc-gc=auto --enable-multiarch --disable-werror --with-arch-32=i686 --with-abi=m64 --with-multilib-list=m32,m64,mx32 --enable-multilib --with-tune=generic --enable-offload-targets=nvptx-none --without-cuda-driver --enable-checking=release --build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu
    Thread model: posix
    gcc version 7.3.0 (Ubuntu 7.3.0-27ubuntu1~18.04)
    opempi-4.0.1

(2) Mac OS High Sierra Version 10.13.6
    Hardware - 1.4GHz Intel Core i5 4GB RAM 256HDD
    Compiler information: 
            COLLECT_GCC=/usr/local/bin/gfortran
            COLLECT_LTO_WRAPPER=/usr/local/libexec/gcc/x86_64-apple-darwin15.6.0/6.2.0/lto-wrapper
            Target: x86_64-apple-darwin15.6.0
            Configured with: ../gcc-6.2.0/configure --enable-languages=c++,fortran --with-gmp=/usr/local
            Thread model: posix
            gcc version 6.2.0 (GCC)
