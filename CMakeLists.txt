cmake_minimum_required(VERSION 3.0)

project(hamsom VERSION 0.1.0 LANGUAGES Fortran)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# various compilation switches

if(HAMSOM_USE_DEBUG)
   add_definitions(-DDEBUG)
endif(HAMSOM_USE_DEBUG)

if(HAMSOM_USE_PROFILING)
   add_definitions(-DPROF)
endif(HAMSOM_USE_PROFILING)

option(HAMSOM_METEO_NCEP "Use NCEP meteo forcing" OFF)
if(HAMSOM_METEO_NCEP)
    message("Using NCEP meteo forcing")
    add_definitions(-DNCEP)
else()
	message("Using COASTDAT meteo forcing")
    add_definitions(-DCOASTDAT) # need to link to NetCDF
endif()

add_definitions(-DSATLAS)
add_definitions(-DSASCII)
add_definitions(-DCARA)
add_definitions(-DNOFISH)

option(HAMSOM_USE_FABM "Include support for Framework for Aquatic Biogeochemical Models (fabm.net)" OFF)
if(HAMSOM_USE_FABM)
   message("Using FABM")
   add_definitions(-D_FABM_)
   set (HAMSOM_FABM hamsom_fabm.F90)
   set(FABM_FORCED_HOST hamsom)
   add_subdirectory(${PROJECT_SOURCE_DIR}/extern/fabm fabm EXCLUDE_FROM_ALL)
else()
   set (HAMSOM_BIO)
endif()

find_package(NetCDF REQUIRED)

option(HAMSOM_USE_PARALLEL "Use MPI" OFF)
if(HAMSOM_USE_PARALLEL)
   message("Using MPI")
   set (HAMSOM_PARALLEL parallel_fabm.F)
   find_package(MPI REQUIRED)
   add_definitions(-DMPI)
endif()

# Customize compiler flags
#[[---
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  add_definitions(-DGFORTRAN)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  add_definitions(-DIFORT)
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
    set(CMAKE_Fortran_STACK_SIZE 128000000)
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
endif()
---]]

# Compile library
add_library(hamsom
            ${HAMSOM_FABM}
            ${HAMSOM_BIO}
#            bio_ad4.F
            hydrodynamics_fabm.F
            ice_dynamics.F
            linear.F
            ncio.f90
	    ${HAMSOM_PARALLEL}
            setup_fabm.F
            subs.F
            trmice_cd2swr.F
            main_fabm.F
)
target_link_libraries(hamsom PRIVATE netcdf)
if(HAMSOM_USE_PARALLEL)
  target_link_libraries(hamsom PRIVATE MPI::MPI_Fortran)
endif()
if(HAMSOM_USE_FABM)
	target_link_libraries(hamsom PUBLIC fabm)
endif()

# Build HAMSOM executable.
add_executable(HAMSOM.out main_fabm.F)
target_include_directories(HAMSOM.out PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_link_libraries(HAMSOM.out PRIVATE hamsom)
if(HAMSOM_USE_PARALLEL)
  target_link_libraries(HAMSOM.out PRIVATE MPI::MPI_Fortran)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install HAMSOM in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install HAMSOM in")
endif()

# Install HAMSOM executable
# install(TARGETS HAMSOM.out DESTINATION bin)
