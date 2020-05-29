
IF( CMAKE_COMPILER_IS_GNUCXX )
  # There are confirmed bugs in 4.8.1 and 4.8.2
  SET( GCC_MIN_VERSION 4.8.3)
  IF (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${GCC_MIN_VERSION})
    MESSAGE(FATAL_ERROR "CovidSim relies on c++11 standards that only exist in g++ version ${GCC_MIN_VERSION} or higher. Current version is ${CMAKE_CXX_COMPILER_VERSION}.")
  ENDIF()
ELSEIF( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  SET( CLANG_MIN_VERSION 3.5 )
  IF (${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS ${CLANG_MIN_VERSION})
    MESSAGE(FATAL_ERROR "CovidSim relies on c++11 standards that only exist in clang version ${CLANG_MIN_VERSION} or higher. Current version is ${CMAKE_CXX_COMPILER_VERSION}.")
  ENDIF()
ENDIF()

IF( APPLE )
  #Apple defines a macro named 'check' in AssertMacros.h unless this is used
  #to suppress the macro definition.
  ADD_DEFINITIONS( -D__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES=0 )
ENDIF()

#==================================================
# Default compiler flags, these can be modified under
# the advanced options using ccmake
#==================================================
IF( NOT DEFINED CMAKE_FLAGS_INIT )

  #===============================
  # Set the build type to release by default, but debug if the binary directory contains the name debug
  SET( BUILD_TYPE_STRING "Choose the type of build, options are: Debug Release RelWithDebInfo Memcheck." )

  IF( NOT CMAKE_BUILD_TYPE )
    IF( BIN_DIR_NAME MATCHES "DEBUG" )
      SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    ELSEIF( BIN_DIR_NAME MATCHES "RELEASE" OR BIN_DIR_NAME MATCHES "DEPLOY" )
      SET(CMAKE_BUILD_TYPE "Release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    ELSEIF( BIN_DIR_NAME MATCHES "RELWITHDEBINFO" )
      SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    ELSEIF( BIN_DIR_NAME MATCHES "MEMCHECK" )
      SET(CMAKE_BUILD_TYPE "Memcheck" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
    ELSE()
      #SET(CMAKE_BUILD_TYPE "Release" CACHE STRING ${BUILD_TYPE_STRING} FORCE)
      SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING ${BUILD_TYPE_STRING} FORCE) #Default to debug for now
    ENDIF()
  ENDIF()
  
  MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})
  
  #=============================

  SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG" CACHE STRING "C++ Release Flags" FORCE )

  #Compiler flags for the C++ compiler
  IF(CMAKE_COMPILER_IS_GNUCXX )

    SET( GNU_WARNING_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wunused-result -Winit-self -Wno-variadic-macros -Wno-vla -Wno-strict-overflow" )
    IF (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER 5 AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 6)
      SET( GNU_WARNING_FLAGS "${GNU_WARNING_FLAGS} -Wno-maybe-uninitialized")
    ENDIF()
    
    SET( CMAKE_CXX_FLAGS "${GNU_WARNING_FLAGS} -std=c++11 -fstrict-aliasing -Wstrict-aliasing -pedantic -Wnon-virtual-dtor" CACHE STRING "C++ Flags" FORCE)
    SET( CMAKE_C_FLAGS "${GNU_WARNING_FLAGS} -fstrict-aliasing -Wstrict-aliasing" CACHE STRING "C Flags" FORCE)
    IF( NOT CYGWIN )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pthread" CACHE STRING "C++ Flags" FORCE)
      SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -pthread" CACHE STRING "C Flags" FORCE)
    ELSE()
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U__STRICT_ANSI__ -Wa,-mbig-obj -Og" CACHE STRING "C++ Flags" FORCE)
    ENDIF()

    SET( CMAKE_CXX_FLAGS_DEBUG "-g -ftrapv -fbounds-check" CACHE STRING "C++ Debug Flags" FORCE )
    IF( NOT CYGWIN )
        SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0" CACHE STRING "C++ Debug Flags" FORCE )
    ELSE()
        SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og" CACHE STRING "C++ Debug Flags" FORCE )
    ENDIF()
    
    SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops" CACHE STRING "C++ Release Flags" FORCE )
    SET( CMAKE_CXX_FLAGS_MEMCHECK "-g -Os -fsanitize=address -fno-omit-frame-pointer" CACHE STRING "C++ Compiler Memory Check Flags" FORCE )
    SET( CMAKE_CXX_FLAGS_VECTORIZE "-O3 -ftree-vectorizer-verbose=7" CACHE STRING "C++ Release Flags" FORCE )

    SET( CMAKE_C_FLAGS_DEBUG "-g -O0 -ftrapv -fbounds-check" CACHE STRING "C Debug Flags" FORCE )
    SET( CMAKE_C_FLAGS_RELEASE "-O3 -funroll-loops" CACHE STRING "C Release Flags" FORCE )

    SET( CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "Flags used by the linker during the creation of dll's." FORCE )

    SET( GNU_NO_INLINE_FLAGS "-DALWAYS_INLINE=inline -fno-inline -fno-inline-functions -fno-inline-small-functions -fno-inline-functions-called-once -fno-default-inline -fno-implicit-inline-templates" )

    SET( CMAKE_EXE_LINKER_FLAGS_MEMCHECK "-fuse-ld=gold -Wl,--disable-new-dtags -Wl,--allow-shlib-undefined" CACHE STRING "Executable Link Flags For Memcheck" FORCE )
  ELSE()
    MESSAGE(FATAL_ERROR "Only works for GNU compilers")
  ENDIF()
  

  IF( MPI_CXX_FOUND )
    SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
  ENDIF()
  
ENDIF()

MARK_AS_ADVANCED( FORCE
                  CMAKE_CXX_FLAGS_DEBUG
                  CMAKE_CXX_FLAGS_RELEASE
                  CMAKE_CXX_FLAGS_MEMCHECK
                  CMAKE_EXE_LINKER_FLAGS_MEMCHECK
                  CMAKE_CXX_FLAGS_VECTORIZE 
                  CMAKE_CXX_FLAGS_ANALYZE
                  CMAKE_C_FLAGS
                  CMAKE_C_FLAGS_DEBUG
                  CMAKE_C_FLAGS_RELEASE
                )

IF( CMAKE_BUILD_TYPE )
  STRING( TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE )
  MARK_AS_ADVANCED( CLEAR CMAKE_CXX_FLAGS )
  MARK_AS_ADVANCED( CLEAR CMAKE_CXX_FLAGS_${BUILD_TYPE} )
ENDIF()

SET( CMAKE_FLAGS_INIT TRUE CACHE INTERNAL "Indicator that this is the first run of cmake" )

#==============================================================================
# Check that the compiler actually works with C++11 features

SET(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE}} ${CMAKE_EXE_LINKER_FLAGS_${BUILD_TYPE}}")

# try to compile a simple program to make sure the c++11 auto feature works
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <vector>
  int main()
  {
    std::vector<double> vec0{1, 2};
    auto vec = vec0;
    return 0;
  }
  "
  CPP11_AUTO_COMPILES)


# try to compile a simple program to make sure the c++11 shared pointer feature works
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <memory>
  int main()
  {
    std::shared_ptr<double> vec0;
    vec0 = std::make_shared<double>(2.0);
    return 0; 
  }
  "
  CPP11_SHAREDPTR_COMPILES)


UNSET(CMAKE_REQUIRED_FLAGS)

IF(NOT CPP11_AUTO_COMPILES OR NOT CPP11_SHAREDPTR_COMPILES)
  MESSAGE( "====================================================================" )
  MESSAGE( "Basic tests of c++11 cannot be compiled.")
  MESSAGE( "Please make sure your compiler supports all c++11 features." )
  MESSAGE( "" )
  MESSAGE( "See CMakeFiles/CMakeError.log for more details.")
  MESSAGE( "====================================================================" )
  MESSAGE( "" )
  MESSAGE( FATAL_ERROR "" )
ENDIF()

