#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.9 FATAL_ERROR)
project(ADHD)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleADHD exampleADHD.cc ${sources} ${headers})
target_link_libraries(exampleADHD ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ADHD. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ADHD_SCRIPTS
  exampleADHD.in 
  exampleADHD.out 
  gui.mac
  init.mac 
  init_vis.mac 
  ionization_check_proton.mac
  ionization_check_He4.mac
  ionization_check_C12.mac 
  vis.mac
  #Particle Gun Macro
  long_run_proton.mac
  long_run_pbar.mac
  long_run_Deuteron.mac
  long_run_Dbar.mac
  long_run_He4.mac
  long_run_C12.mac
  ev16.mac
  antimatter_effect.mac
  #General Particle Source macros
  long_run_proton_GPS.mac
  long_run_pbar_GPS.mac
  long_run_Deuteron_GPS.mac
  long_run_Dbar_GPS.mac
  long_run_He4_GPS.mac
  long_run_C12_GPS.mac  
  long_run_electron_GPS.mac
  #
  long_run_Dbar_GPS_stat.mac
  long_run_pbar_GPS_stat.mac
  )

foreach(_script ${ADHD_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ADHD DEPENDS exampleADHD)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleADHD DESTINATION bin)

