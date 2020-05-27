 #
 # Find the NLopt includes and libraries
 #
 # NLopt Can be found at
 #     http://ab-initio.mit.edu/wiki/index.php/NLopt
 #

PKG_CHECK_MODULES(NLOPT QUIET libnlopt)

IF ( NLOPT_FOUND )
    SET(NLOPT_DEFINITIONS ${NLOPT_CFLAGS_OTHER})
ELSE()
    FIND_PATH(NLOPT_INCLUDE_DIRS nlopt.h
          HINTS $ENV{NLOPT_DIR}/include
          PATH_SUFFIXES nlopt )

    FIND_LIBRARY(NLOPT_LIBRARIES NAMES nlopt nlopt_cxx
             HINTS $ENV{NLOPT_DIR}/lib )
ENDIF()


INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(NLOPT  DEFAULT_MSG
                                  NLOPT_LIBRARIES NLOPT_INCLUDE_DIRS)

MARK_AS_ADVANCED(NLOPT_INCLUDE_DIRS NLOPT_LIBRARIES)

