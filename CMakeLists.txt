cmake_minimum_required(VERSION 3.14)

project(
    analysis-module
    VERSION 0.1.0
    DESCRIPTION "View hdf5 data recorded by the Data Recorder plugin"
    HOMEPAGE_URL "https://rtxi.org/"
    LANGUAGES C CXX
)

# These lines help with third-party tooling integration
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(RTXI_CMAKE_SCRIPTS "" CACHE PATH "In development mode, RTXI uses this to expose conan dependencies")
set(CMAKE_PREFIX_PATH "${RTXI_CMAKE_SCRIPTS}")
set(RTXI_PACKAGE_PATH "/usr/local/" CACHE PATH "Path hint to RTXI package information")
set(CMAKE_INSTALL_RPATH "${RTXI_PACKAGE_PATH}/lib")
list(APPEND CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
# Required to find qwt
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ---- find libraries ----
find_package(rtxi REQUIRED HINTS ${RTXI_PACKAGE_PATH})
find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets HINTS ${RTXI_CMAKE_SCRIPTS})
find_package(fmt REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(qwt REQUIRED)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


#################################################################################################
### You can modify within this region for finding additional libraries and linking them to your #
### custom plugin. Make sure to install them prior to configuration or else build will fail     #
### with linking and include errors!                                                            # 
#################################################################################################
add_library(
    analysis-module MODULE
    widget.cpp
    widget.hpp
)

# Consult library website for how to link them to your plugin using cmake
target_link_libraries(analysis-module PUBLIC 
    rtxi::rtxi rtxi::rtxidsp rtxi::rtxiplot rtxi::rtxififo Qt5::Core Qt5::Gui Qt5::Widgets dl fmt::fmt
    hdf5::hdf5_hl hdf5::hdf5 qwt::qwt
)

################################################################################################ 

# We need to tell cmake to use the c++ version used to compile the dependent library or else...
get_target_property(REQUIRED_COMPILE_FEATURE rtxi::rtxi INTERFACE_COMPILE_FEATURES)
target_compile_features(analysis-module PRIVATE ${REQUIRED_COMPILE_FEATURE})

install(
    TARGETS analysis-module
    DESTINATION ${RTXI_PACKAGE_PATH}/bin/rtxi_modules
)

