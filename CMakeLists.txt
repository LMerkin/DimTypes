# vim:ts=2:et
#=============================================================================#
#                        "CMakeLists.txt" (Top-Level):                        #
#=============================================================================#
#=============================================================================#
# CMake Settings:                                                             #
#=============================================================================#
# Export flags for convenience when using Vim's YCM plugin:
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.27)
CMAKE_POLICY(          VERSION 3.27)

set(TOOL_CHAIN "GCC" CACHE STRING "ToolChain, currently supported: GCC (default), CLang, NVHPC")
set_property(CACHE TOOL_CHAIN PROPERTY STRINGS GCC CLang NVHPC)
set(PROJ_NAME "DimTypes" CACHE STRING "Project name")
set(LIB_DIR "" CACHE STRING "Library output directory")
set(BIN_DIR "" CACHE STRING "Executable output directory")
option(UNCHECKED_MODE "Set UNCHECKED_MODE" OFF)

IF (NOT UNIX)
  MESSAGE(FATAL_ERROR "Non-UNIX platforms are not supported")
ENDIF()

#=============================================================================#
# Project Settings:                                                           #
#=============================================================================#
# NB: Also, CUDA_VER, CUDA_CC and CUDA_INCL may be passed from the CMake cmdl,
# but it is not an error if no CUDA is provided...

# NB: In this case, PROJECT_SOURCE_DIR = CMAKE_SOURCE_DIR
#              and  PROJECT_BINARY_DIR = CMAKE_BINARY_DIR ,
# so we use the "CMAKE_"-prefixed vars:
#
PROJECT(${PROJ_NAME} LANGUAGES CXX)

SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

IF (LIB_DIR)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}")
ENDIF()
IF (BIN_DIR)
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}")
ENDIF()

#=============================================================================#
# Compiler Settings:                                                          #
#=============================================================================#
MESSAGE(STATUS "Configuring for the ${TOOL_CHAIN} ToolChain...")

#-----------------------------------------------------------------------------#
IF("${TOOL_CHAIN}" STREQUAL "GCC")
#-----------------------------------------------------------------------------#
  # Other flags:
  SET(CXXFLAGS "${CXXFLAGS} \
    -march=native        -mtune=native        -std=gnu++23 \
    -Wall                -Wextra -Wshadow     -Wconversion \
    -Wold-style-cast     -Wcast-qual          -Wnon-virtual-dtor \
    -Woverloaded-virtual -Wctor-dtor-privacy  -Wnull-dereference \
    -Wlogical-op         -Wduplicated-cond    -Wduplicated-branches \
    -Wformat=2           -Wcast-align=strict")
  # Optimisation and debugging flags:
  # NB: Do NOT use for release build:
  # "-ffast-math" option -- it will result in incorrect treatment of "NaN"s!
  # -flto"               -- it intereacts in a strange way with header-only
  #                         code!
  SET(CXXFLAGS_RELEASE "-O3 -mavx2 -fno-semantic-interposition -fomit-frame-pointer ${CXXFLAGS}")
  SET(CXXFLAGS_RELWITHDEB "-O2 -g ${CXXFLAGS}")
  SET(CXXFLAGS_DEBUG "-O0 -g -fstack-protector-strong ${CXXFLAGS}")

#-----------------------------------------------------------------------------#
ELSEIF("${TOOL_CHAIN}" STREQUAL "CLang")
#-----------------------------------------------------------------------------#
  # Other flags:
  SET(CXXFLAGS "${CXXFLAGS} \
    -march=native        -mtune=native        -std=gnu++23 \
    -Weverything         -Wno-c++98-compat    -Wno-c++98-compat-pedantic \
    -Wno-padded          -Wno-weak-vtables    -Wextra      \
    -Wshadow-all         -Wconversion         -Wold-style-cast \
    -Wcast-qual          -Wnull-dereference   -Wlogical-op-parentheses   \
    -Wnon-virtual-dtor   -Woverloaded-virtual -Wdelete-non-virtual-dtor  \
    -Wformat=2           -Wcast-align         -Wextra-semi \
    -Wunreachable-code-aggressive             -Wrange-loop-analysis \
    -Wpessimizing-move   -Wredundant-move     -Wno-double-promotion \
    -Wno-float-equal     -Wno-switch-enum")
  # Optimisation and debugging flags:
  SET(CXXFLAGS_RELEASE "-O3 -mavx2 -fno-semantic-interposition -fomit-frame-pointer")
  SET(CXXFLAGS_RELWITHDEB "-O2 -g")
  SET(CXXFLAGS_DEBUG
      "-O0 -g -fstack-protector-strong -fsanitize=undefined,address \
       -fno-omit-frame-pointer -DDEBUG")

#-----------------------------------------------------------------------------#
ELSEIF("${TOOL_CHAIN}" STREQUAL "NVHPC")
#-----------------------------------------------------------------------------#
  # Other flags:
  SET(CXXFLAGS  "${CXXFLAGS} -march=native -mtune=native    -std=c++23 \
                -Wno-switch-enum           -Wno-float-equal \
                -Wno-missing-prototypes    -Wno-vla")
  # Optimisation and debugging flags:
  SET(CXXFLAGS_RELEASE "-O3 -mavx2")
  SET(CXXFLAGS_RELWITHDEB "-O2 -g")
  SET(CXXFLAGS_DEBUG "-O0 -g")

#-----------------------------------------------------------------------------#
ELSE()
  MESSAGE(FATAL_ERROR  "UnSupported ToolChain: ${TOOL_CHAIN}")
ENDIF()

#=============================================================================#
# Compiler-Independent Compilation Env:                                       #
#=============================================================================#
# Reset the standard flags for the selected mode, to prevent a mix-up:
SET(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS} ${CXXFLAGS_RELEASE}")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CXXFLAGS_RELWITHDEB}")
SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS} ${CXXFLAGS_DEBUG}")

ADD_DEFINITIONS(-DUNCHECKED_MODE=${UNCHECKED_MODE})
IF (UNCHECKED_MODE)
  ADD_DEFINITIONS(-DNDEBUG)
ENDIF()

#-----------------------------------------------------------------------------#
# Common Configs:                                                             #
#-----------------------------------------------------------------------------#
# Use CCache if available:
OPTION(USE_CCACHE "Use ccache" OFF)
IF (USE_CCACHE)
  FIND_PROGRAM(CCACHE_FOUND ccache)
  IF (CCACHE_FOUND)
    SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK    ccache)
  ENDIF()
ENDIF()

# During build, link executables with build rpath, not install rpath; they are
# then automatically re-linked on install. This is the default behaviour, any-
# way:
SET(CMAKE_SKIP_BUILD_RPATH         FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/Include")
LINK_DIRECTORIES   ("${CMAKE_BINARY_DIR}" "${LIB_DIR}")

#=============================================================================#
# Tests:                                                                      #
#=============================================================================#
SET(DIM_TESTS
  DimTest
  CEMathsTest
)

SET(SOURCES Include/DimTypes/DimTypes.hpp
            Include/DimTypes/Bits/CEMaths.hpp
            Include/DimTypes/Bits/Encodings.hpp
            Include/DimTypes/Bits/Macros.h)
add_library(${PROJ_NAME} INTERFACE ${SOURCES})

FOREACH (DimTest ${DIM_TESTS})
  ADD_EXECUTABLE(${DimTest} Tests/${DimTest}.cpp)
  target_link_libraries(${DimTest} ${PROJ_NAME})
ENDFOREACH()
