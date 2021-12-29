INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( OptimizersLib_SRC
    Optimization.cpp
    OptimizerFull.cpp
    OptimizerPiecewise.cpp
    OptimizerPiecewise2Layer.cpp
    OptimizerPiecewise4Layer.cpp
    OptimizerGlobalBasis.cpp
  )

ADD_LIBRARY( OptimizersLib STATIC ${OptimizersLib_SRC} )

#

