cmake_minimum_required (VERSION 2.8.11)

project (vandmc)

#Compile with c++11 or c++0x depending on the compiler version
if (CMAKE_COMPILER_IS_GNUCXX)
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.3" AND
            ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    elseif(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.6")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    endif(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.3" AND
            ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7")

    #Make compiler messages nice with colored tags.
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.9")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=auto")
	endif()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -fPIC")
endif()

#if user does not specify prefix we assign it to the exec directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   message(STATUS "Install Prefix not specified.")
   file(MAKE_DIRECTORY install)
   get_filename_component(INSTALL_DIR ${CMAKE_SOURCE_DIR}/install REALPATH)
   set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "Install Prefix" FORCE)
endif()
message(STATUS "Installing to ${CMAKE_INSTALL_PREFIX}")

#Define the default build type to be Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Build type, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
   message(STATUS "Build type not defined, using default: ${CMAKE_BUILD_TYPE}")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#Add additional CXX flags if we use the Debug option
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

#------------------------------------------------------------------------------

#Install options
option(BUILD_TOOLS "Build and install tool programs." OFF)
option(BUILD_RENDERER "Build and install vandmc detector renderer." OFF)

#------------------------------------------------------------------------------

#Find required packages.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

#Find simpleScan install (https://github.com/cthornsb/SimplePixieScan)
find_package (SimpleScan REQUIRED)
include_directories(${SimpleScan_INCLUDE_DIR})
mark_as_advanced(SimpleScan_SCAN_LIB SimpleScan_OPT_LIB)

#Find rcbuild executable.
find_program(RCBUILD_EXECUTABLE NAMES rcbuild PATHS ${SimpleScan_BIN_DIR})
if("${RCBUILD_EXECUTABLE}" STREQUAL "RCBUILD_EXECUTABLE-NOTFOUND")
	message(FATAL_ERROR "Could NOT find rcbuild (https://github.com/cthornsb/SimplePixieScan)")
else()
	message(STATUS "Found rcbuild: ${RCBUILD_EXECUTABLE}")
endif()

if(BUILD_RENDERER)
	#Find qt install.
	find_package (QT REQUIRED)
endif()

#Find ROOT install.
find_package (ROOT REQUIRED)
mark_as_advanced(FORCE GENREFLEX_EXECUTABLE ROOTCINT_EXECUTABLE ROOT_CONFIG_EXECUTABLE)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(TOP_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

set(DICTIONARY_PREFIX "VandmcDict" CACHE STRING "Prefix to root dictionary.")
set(DEF_STRUCT_FILE "${TOP_DIRECTORY}/dict/def.struct" CACHE PATH "Root structure definition file.")

add_subdirectory(dict)

include_directories(include)
add_subdirectory(source)

if(BUILD_RENDERER)
	add_subdirectory(tools/qtfiles)
endif()

if(BUILD_TOOLS)
	add_subdirectory(tools)
endif()

#Build/install the miscellaneous stuff
add_subdirectory(share)
