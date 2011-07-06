cmake_minimum_required (VERSION 2.6)
project (READ_HYCOM)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_COMPILER "ifort")

# FFLAGS depend on the compiler
if (CMAKE_Fortran_COMPILER STREQUAL "ifort")
   set (CMAKE_Fortran_FLAGS_RELEASE "-convert big_endian")
endif (CMAKE_Fortran_COMPILER STREQUAL "ifort")

# build executables
set (NMPROGRAMS "fortran/read_hycom_zlevel.f90" "fortran/htn.f90")
set (EXECUTABLES "read_hycom.x" ${NMPROGRAMS})

add_executable ( ${EXECUTABLES} )
add_executable ( test.x "fortran/test.f90" )
