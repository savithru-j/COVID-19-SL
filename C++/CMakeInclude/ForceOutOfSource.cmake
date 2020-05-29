##======================================================================
##
## Force build to be out of source
##
##======================================================================
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" sourcetree)
GET_FILENAME_COMPONENT(PARENTDIR ${CMAKE_SOURCE_DIR} PATH)
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PARENTDIR}" sourcesubtree)
IF(sourcetree OR sourcesubtree)
  MESSAGE( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
  MESSAGE( "!!  Please do not run cmake in the source tree. Run it in a build directory." )
  MESSAGE( "!!  CMake will not work until you DELETE the generated file and folder:" )
  MESSAGE( "!!" )
  MESSAGE( "!!  ${CMAKE_SOURCE_DIR}/CMakeFiles/" )
  MESSAGE( "!!  ${CMAKE_SOURCE_DIR}/CMakeCache.txt" )
  MESSAGE( "!!" )
  MESSAGE( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
  MESSAGE( "" )
  MESSAGE(FATAL_ERROR "")
ENDIF()
