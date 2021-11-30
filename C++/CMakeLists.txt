INCLUDE( CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

#Get the name of the binary directory
STRING( TOUPPER ${CMAKE_BINARY_DIR} BIN_DIR_NAME )

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)
PROJECT(CovidSim)

SET(CMAKE_CXX_STANDARD 11) # C++11...
SET(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
SET(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

#set( CMAKE_VERBOSE_MAKEFILE ON )

MESSAGE( "CMAKE_SOURCE_DIR:" ${CMAKE_SOURCE_DIR} )
MESSAGE( "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR} )

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeInclude )

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)

#Enable testing for the project
enable_testing()
 
#Bring the headers into the project
INCLUDE_DIRECTORIES(src)

INCLUDE(FindPkgConfig)
INCLUDE( CMakeInclude/FindNLOPT.cmake )
INCLUDE( CMakeInclude/CompileFlags.cmake )

SET(LIBS
    OptimizersLib
    ModelsLib
    LinearAlgebraLib
    ${NLOPT_LIBRARIES}
   )
   
ADD_SUBDIRECTORY( src )
#ADD_SUBDIRECTORY( test )

#set(LIB_NAME ${PROJECT_NAME}lib)
#add_library(${LIB_NAME} SHARED ${SOURCES})
 
#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")
 
add_executable(CovidOptFull src/CovidOptFull.cpp)
target_link_libraries(CovidOptFull ${LIBS})

add_executable(CovidOptPiecewise src/CovidOptPiecewise.cpp)
target_link_libraries(CovidOptPiecewise ${LIBS})

add_executable(CovidOptPiecewise4Layer src/CovidOptPiecewise4Layer.cpp)
target_link_libraries(CovidOptPiecewise4Layer ${LIBS})

add_executable(CovidOptGlobalBasis src/CovidOptGlobalBasis.cpp)
target_link_libraries(CovidOptGlobalBasis ${LIBS})
 
add_executable(CovidSim src/CovidSim.cpp)
target_link_libraries(CovidSim ${LIBS})

