set(CMAKE_Fortran_COMPILER "/opt/psi/Programming/intel/22.1/compiler/latest/linux/bin/intel64/ifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "20.2.5.20211109")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "i386-linux-gnu")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "i386-linux-gnu")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/opt/psi/Programming/intel/22.1/vpl/2022.0.0/include;/opt/psi/Programming/intel/22.1/tbb/2021.5.0/include;/opt/psi/Programming/intel/22.1/mpi/2021.5.0/include;/opt/psi/Programming/intel/22.1/mkl/2022.0.1/include;/opt/psi/Programming/intel/22.1/ipp/2021.5.1/include;/opt/psi/Programming/intel/22.1/ippcp/2021.5.0/include;/opt/psi/Programming/intel/22.1/dpl/2021.6.0/linux/include;/opt/psi/Programming/intel/22.1/dpcpp-ct/2022.0.0/include;/opt/psi/Programming/intel/22.1/dnnl/2022.0.1/cpu_dpcpp_gpu_dpcpp/lib;/opt/psi/Programming/intel/22.1/dev-utilities/2021.5.1/include;/opt/psi/Programming/intel/22.1/dal/2021.5.1/include;/opt/psi/Programming/intel/22.1/ccl/2021.5.0/include/cpu_gpu_dpcpp;/afs/psi.ch/sys/psi.merlin/Programming/intel/22.1/compiler/2022.0.1/linux/compiler/include/intel64;/afs/psi.ch/sys/psi.merlin/Programming/intel/22.1/compiler/2022.0.1/linux/compiler/include/icc;/afs/psi.ch/sys/psi.merlin/Programming/intel/22.1/compiler/2022.0.1/linux/compiler/include;/usr/local/include;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0/include-fixed;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/opt/psi/Programming/psi-python39/2021.11/lib;/opt/psi/Tools/emacs/29.1/lib64;/opt/psi/Tools/emacs/29.1/lib;/opt/psi/Tools/gnutls/3.5.19/lib64;/opt/psi/Tools/gnutls/3.5.19/lib;/opt/psi/Programming/Tk/8.6.4/lib64;/opt/psi/Programming/Tk/8.6.4/lib;/opt/psi/Programming/Tcl/8.6.4/lib64;/opt/psi/Programming/Tcl/8.6.4/lib;/opt/psi/Programming/intel/22.1/vpl/2022.0.0/lib;/opt/psi/Programming/intel/22.1/tbb/2021.5.0/lib/intel64/gcc4.8;/opt/psi/Programming/intel/22.1/mpi/2021.5.0/lib/release;/opt/psi/Programming/intel/22.1/mpi/2021.5.0/lib;/opt/psi/Programming/intel/22.1/mkl/2022.0.1/lib/intel64;/opt/psi/Programming/intel/22.1/ipp/2021.5.1/lib/intel64;/opt/psi/Programming/intel/22.1/ippcp/2021.5.0/lib/intel64;/opt/psi/Programming/intel/22.1/dnnl/2022.0.1/cpu_dpcpp_gpu_dpcpp/lib;/opt/psi/Programming/intel/22.1/dal/2021.5.1/lib/intel64;/opt/psi/Programming/intel/22.1/compiler/2022.0.1/linux/compiler/lib/intel64_lin;/opt/psi/Programming/intel/22.1/compiler/2022.0.1/linux/lib;/opt/psi/Programming/intel/22.1/clck/2021.5.0/lib/intel64;/opt/psi/Programming/intel/22.1/ccl/2021.5.0/lib/cpu_gpu_dpcpp;/afs/psi.ch/sys/psi.merlin/Programming/intel/22.1/compiler/2022.0.1/linux/compiler/lib/intel64_lin;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib/gcc/x86_64-pc-linux-gnu/9.3.0;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib/gcc;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib64;/lib64;/usr/lib64;/afs/psi.ch/sys/psi.merlin/Programming/gcc/9.3.0/lib;/lib;/usr/lib;/usr/lib/i386-linux-gnu")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
