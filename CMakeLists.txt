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
set (NMPROGRAMS1 "fortran/read_hycom_zlevel.f90" "fortran/htn.f90")
set (NMPROGRAMS2 "fortran/read_hycom_2d.f90" "fortran/htn_2d.f90")
set (EXECUTABLES1 "read_hycom.x" ${NMPROGRAMS1})
set (EXECUTABLES2 "read_hycom_2d.x" ${NMPROGRAMS2})

add_executable ( ${EXECUTABLES1} )
add_executable ( ${EXECUTABLES2} )
add_executable ( test.x "fortran/test.f90" )
