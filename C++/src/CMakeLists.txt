INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( CovidLib_SRC
    Population.cpp
    Simulator.cpp
    Optimizer.cpp
  )

ADD_LIBRARY( CovidLib STATIC ${CovidLib_SRC} )

#ADD_SUBDIRECTORY( Grid )

