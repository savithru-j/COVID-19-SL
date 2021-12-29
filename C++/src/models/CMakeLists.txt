INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( ModelsLib_SRC
    Population.cpp
    Population2Layer.cpp
    Population4Layer.cpp
    ObservedPopulation.cpp
    Simulator.cpp
    ModelParams.cpp
    ModelParams2Layer.cpp
    ModelParams4Layer.cpp
  )

ADD_LIBRARY( ModelsLib STATIC ${ModelsLib_SRC} )

#

