cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

# set not CYGWIN_WIN32, to quench warning messages from CYGWIN cmake
# will do no harm if not CYGWIN
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# set project name and languanges that are involved
project(DALTON Fortran C CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# To link with fortran not c++
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/mpi
    ${CMAKE_SOURCE_DIR}/cmake/fortran-standard
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_BINARY_DIR}
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"                           OFF)
option(ENABLE_GPU            "Enable GPU acceleration"                          OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"                     OFF)
option(ENABLE_CUBLAS         "Enable the CUBLAS GPU library"                    OFF)
option(ENABLE_COLLAPSE       "Enable COLLAPSE clause for reorderings"           OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                              OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                             OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                       OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking"                  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int library"                           ON)
option(ENABLE_PELIB          "Enable Polarizable Embedding (PE) library"        ON)
option(ENABLE_QFITLIB        "Enable charge fitting library (QFITLIB)"          ON)
option(ALWAYS_RESET_EXTERNAL "Always remove builds stamps of external projects" ON)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                                 OFF)
option(ENABLE_XCFUN          "Enable XCFUN"                                     ON)
option(ENABLE_CRAY_WRAPPERS  "Enable cray wrappers for BLAS/LAPACK and MPI"     OFF)
option(ENABLE_LSEEK          "Enable lseek"                                     OFF)
option(ENABLE_VPOTDAMP       "Enable VPOTDAMP library"                          OFF)
option(ENABLE_QMMM_CUDA      "Enable QMMM CUDA library"                         OFF) # radovan: untested
option(ENABLE_CHEMSHELL      "Compile for ChemShell"                            OFF)
option(ENABLE_QCMATRIX       "Enable QcMatrix library"                          OFF)
option(ENABLE_OPENRSP        "Enable OPENRSP library (need QcMatrix library)"   OFF)
option(ENABLE_EFS            "Enable Echidna Fock Solver"                       OFF)
option(ENABLE_BUILTIN_BLAS   "Enable builtin BLAS implementation (slow)"        OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_BUILTIN_LAPACK "Enable builtin LAPACK implementation (slow)"      OFF)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_LARGE_TEST     "Enable Large-scale (long) test cases"             OFF)
option(ENABLE_DEVELOPER      "Enable Developer test cases and features"         OFF)
option(ENABLE_PYTHON         "Enable Python driver script"                      OFF)

set(EXTERNAL_LIBS)

include(ConfigVersion)
include(SourcesDALTON)
include(ConfigArchitecture)

include(ConfigOMP)
if(ENABLE_OPENMP)
    add_definitions(-DVAR_OMP)
    set(ENABLE_THREADED_MKL TRUE)
else()
    set(ENABLE_THREADED_MKL FALSE)
endif()

include(ConfigCompilerFlags)
include(ConfigExternal)

# math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

include(ConfigExplicitLibs)
include(ConfigMPI)
include(ConfigSafeGuards)
include(GenericMacros)
include(BinaryInfo)

# set code coverage
if(ENABLE_CODE_COVERAGE)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} gcov)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

if(ENABLE_GPU)
    message(STATUS "Enabling GPU acceleration (OpenACC/CUDA).")
    message(STATUS "If your compiler does not support OpenACC, this *will* fail.")
    if(NOT ENABLE_CRAY_WRAPPERS AND NOT ENABLE_CUBLAS)
        set(ENABLE_GPU OFF CACHE BOOL "Enable GPU acceleration" FORCE)
        set(ENABLE_OPENACC OFF CACHE BOOL "Enable OpenACC GPU acceleration" FORCE)
        set(ENABLE_CUDA OFF CACHE BOOL "Enable CUDA GPU acceleration" FORCE)
        message(STATUS "No GPU math library found, so GPU acceleration is disabled!")
        message(STATUS "GPU acceleration currently only works in combination with:")
        message(STATUS "'--cray' (Libsci_acc) OR '--cublas' (CUBLAS)")
    endif()
endif()

include_directories(
    ${CMAKE_SOURCE_DIR}/DALTON/include
    ${CMAKE_SOURCE_DIR}/DALTON/dft
    ${CMAKE_SOURCE_DIR}/DALTON/abacus
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${PROJECT_BINARY_DIR}
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

# set definitions
include(Definitions)
if(ENABLE_OPENRSP)
    include(DefinitionsOpenRSP)
endif()

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

# PCMSolver configuration
include(pcmsolver)

# configure INSTALL_WRKMEM and INSTALL_MMWORK
include(ConfigWORKMEM)

set(opt_lsdalton "0") #FALSE, ie. dalton
if(ENABLE_PYTHON)
configure_script(
    ${CMAKE_SOURCE_DIR}/dalton.py
    ${CMAKE_BINARY_DIR}/dalton
    )
else()
configure_script(
    ${CMAKE_SOURCE_DIR}/dalton.in
    ${CMAKE_BINARY_DIR}/dalton
    )
endif()
configure_script(
    ${CMAKE_SOURCE_DIR}/DALTON/dalton_config.in
    ${CMAKE_BINARY_DIR}/dalton_config.h
    )

# if BLAS and/or LAPACK not found, add own sources to the list of
# sources to compile
if(USE_BUILTIN_BLAS)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_BLAS_SOURCES}
        )
endif()
if(USE_BUILTIN_LAPACK)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_LAPACK_SOURCES}
        )
endif()

include(LibsDALTON)

include(ConfigTesting)

# copy basis/ to build directory
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})

# copy test scripts to build/test
# FIXME this will not work on windows
execute_process(COMMAND mkdir -p ${CMAKE_BINARY_DIR}/test)
execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/DALTON/test/runtest_v1.py        ${CMAKE_BINARY_DIR}/test)
execute_process(COMMAND cp ${CMAKE_SOURCE_DIR}/DALTON/test/runtest_dalton.py ${CMAKE_BINARY_DIR}/test)

# this controlls "make install" target
include(ConfigMakeInstall)

# give information about system, compiler flags, and size of static allocations
set(STATIC_MEM_INFO_BINARIES dalton)
include(ConfigInfo)
