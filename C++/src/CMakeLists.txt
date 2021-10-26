INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

ADD_SUBDIRECTORY( linearalgebra )
ADD_SUBDIRECTORY( models )
ADD_SUBDIRECTORY( optimizers )
ADD_SUBDIRECTORY( utils )

