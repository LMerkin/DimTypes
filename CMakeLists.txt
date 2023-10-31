# vim:ts=2:et
#=============================================================================#
#                        "CMakeLists.txt" (Top-Level):                        #
#=============================================================================#
#=============================================================================#
# CMake Settings:                                                             #
#=============================================================================#
# Export flags for convenience when using Vim's YCM plugin:
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.23)
CMAKE_POLICY(          VERSION 3.23)

IF (NOT UNIX)
  MESSAGE(FATAL_ERROR "Non-UNIX platforms are not supported")
ENDIF()

#=============================================================================#
# Project Settings:                                                           #
#=============================================================================#
PROJECT(${PROJECT_NAME} LANGUAGES CXX)

SET (CMAKE_EXPORT_COMPILE_COMMANDS ON)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# NB: The following vars are normally provided via the command-line params:
IF (LIB_DIR)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${LIB_DIR}")
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}")
ENDIF()

IF (BIN_DIR)
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}")
ENDIF()

IF (NOT TOOL_CHAIN)
  SET(TOOL_CHAIN "NVHPC")
ENDIF()

IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "RelWithDebInfo")
ENDIF()

IF (NOT ENV_PREFIX)
  SET(ENV_PREFIX "/opt")
ENDIF()

#=============================================================================#
# Compiler Settings:                                                          #
#=============================================================================#
MESSAGE(STATUS "Configuring for the ${TOOL_CHAIN} ToolChain...")

#-----------------------------------------------------------------------------#
IF("${TOOL_CHAIN}" STREQUAL "GCC")
#-----------------------------------------------------------------------------#
  # Optimisation and debugging flags:
  IF ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    # NB: Do NOT use:
    # "-ffast-math" option -- it will result in incorrect treatment of "NaN"s!
    # -flto"               -- it intereacts in a strange way with header-only
    #                         code!
    SET(CXXFLAGS "-O3 -mavx2 -DNDEBUG -funroll-loops -fomit-frame-pointer")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    SET(CXXFLAGS "-O2 -g -funroll-loops")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    SET(CXXFLAGS "-O0 -g")
  ELSE()
    MESSAGE(WARNING "UnSupported BuildType: ${CMAKE_BUILD_TYPE}")
  ENDIF()
  # Other flags:
  # In G++ 12+, there are lots of bogus warnings regarding dangling pointers in
  # our macros, so we have disabled them:
  SET(CXXFLAGS "${CXXFLAGS} -march=native  -mtune=native \
                -Wall                    -Wextra          -Wcast-qual \
                -Winit-self              -Wwrite-strings  -Wmissing-braces \
                -Wdisabled-optimization  -Wconversion  \
                -Wno-stringop-truncation -Wno-dangling-pointer \
                -std=gnu++23 \
                -Wctor-dtor-privacy  -Wnon-virtual-dtor -Woverloaded-virtual \
                -Wno-class-memaccess -Wno-double-promotion")

  # CXX Flags for OpenACC: TODO...
#-----------------------------------------------------------------------------#
ELSEIF("${TOOL_CHAIN}" STREQUAL "CLang")
#-----------------------------------------------------------------------------#
  # Optimisation and debugging flags:
  IF ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(CXXFLAGS "-O3 -mavx2 -DNDEBUG -funroll-loops -fomit-frame-pointer")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    SET(CXXFLAGS "-O2 -g -funroll-loops")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"  )
    SET(CXXFLAGS "-O0 -g")
  ELSE()
    MESSAGE(WARNING "UnSupported BuildType: ${CMAKE_BUILD_TYPE}")
  ENDIF()
  # Other flags:
  SET(CXXFLAGS "${CXXFLAGS} -march=native  -mtune=native \
                -Weverything             -Wno-covered-switch-default  \
                -Wno-switch-enum         -Wno-exit-time-destructors   \
                -Wno-global-constructors -Wno-reserved-id-macro       \
                -Wno-padded              -Wno-unused-macros           \
                -Wno-gnu-zero-variadic-macro-arguments \
                -Wno-disabled-macro-expansion \
                -Wno-float-equal         -Wno-missing-prototypes \
                -Wno-vla                 -Wno-vla-extension      \
                -Wno-c99-extensions      -Wno-zero-length-array  \
                -Wno-packed              -Wno-undefined-func-template \
                -Wno-gnu-statement-expression \
                -std=gnu++23 \
                -Wno-c++98-compat        -Wno-c++98-compat-pedantic    \
                -Wno-weak-vtables        -Wno-double-promotion         \
                -Wno-weak-template-vtables    \
                -Wno-unsafe-buffer-usage -Wno-unused-const-variable    \
                -Wno-unused-template")

#-----------------------------------------------------------------------------#
ELSEIF("${TOOL_CHAIN}" STREQUAL "NVHPC")
#-----------------------------------------------------------------------------#
  # Optimisation and debugging flags:
  IF ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(CXXFLAGS "-O3 -mavx2 -DNDEBUG")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    SET(CXXFLAGS "-O2 -g")
  ELSEIF ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"  )
    SET(CXXFLAGS "-O0 -g")
  ELSE()
    MESSAGE(WARNING "UnSupported BuildType: ${CMAKE_BUILD_TYPE}")
  ENDIF()

  # Other flags:
  SET(CXXFLAGS "${CXXFLAGS} -march=native  -mtune=native  \
                -Wno-switch-enum         -Wno-float-equal \
                -Wno-missing-prototypes  -Wno-vla         \
                -std=c++23")
#-----------------------------------------------------------------------------#
ELSE()
  MESSAGE(FATAL_ERROR  "UnSupported UNIX ToolChain: ${TOOL_CHAIN}")
ENDIF()

#=============================================================================#
# Compiler-Independent Compilation Env:                                       #
#=============================================================================#
# Reset the standard flags for the selected mode, to prevent a mix-up:
SET(CMAKE_C_FLAGS_RELEASE          "")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO   "")
SET(CMAKE_C_FLAGS_DEBUG            "")

SET(CMAKE_CXX_FLAGS_RELEASE        "")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")
SET(CMAKE_CXX_FLAGS_DEBUG          "")

#-----------------------------------------------------------------------------#
# Special defs (mostly for Boost) and libs:                                   #
#-----------------------------------------------------------------------------#
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")

ADD_DEFINITIONS(-DUNCHECKED_MODE=${UNCHECKED_MODE})
# Ignore warnings in system headers:
ADD_DEFINITIONS(-isystem "${ENV_PREFIX}/include")

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
SET(CMAKE_INSTALL_RPATH            ${CMAKE_INSTALL_PREFIX}/lib)

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/Include"
                    "${ENV_PREFIX}/include")

LINK_DIRECTORIES   ("${PROJECT_BINARY_DIR}"
                    "${LIB_DIR}"
                    "${ENV_PREFIX}/lib")

#=============================================================================#
# Tests:                                                                      #
#=============================================================================#
SET(DIM_TESTS
  DimTest
  CEMathsTest
)

FOREACH (DimTest ${DIM_TESTS})
  ADD_EXECUTABLE(${DimTest} Tests/${DimTest}.cpp)
ENDFOREACH()
