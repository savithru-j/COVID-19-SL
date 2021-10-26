INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( ModelsLib_SRC
    Population.cpp
    ObservedPopulation.cpp
    Simulator.cpp
    ModelParams.cpp
  )

ADD_LIBRARY( ModelsLib STATIC ${ModelsLib_SRC} )

#

