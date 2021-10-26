INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( OptimizersLib_SRC
    Optimization.cpp
    OptimizerFull.cpp
    OptimizerPiecewise.cpp
    OptimizerGlobalBasis.cpp
  )

ADD_LIBRARY( OptimizersLib STATIC ${OptimizersLib_SRC} )

#

