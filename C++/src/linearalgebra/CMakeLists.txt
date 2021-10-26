INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( LinearAlgebraLib_SRC
	Matrix.cpp
  )

ADD_LIBRARY( LinearAlgebraLib STATIC ${LinearAlgebraLib_SRC} )

