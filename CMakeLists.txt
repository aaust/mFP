cmake_minimum_required(VERSION 2.6)
project(mFP)

set(MYFIT_SOURCES myFit.cc control.cc wave.cc event.cc likelihood.cc 3j.cc fitInfo.cc fitInfoDict.cc waveDict.cc)
set(MASSDEP_SOURCES massDep.cc bw.cc fitInfo.cc fitInfoDict.cc fitModel.cc wave.cc event.cc waveDict.cc)
set(PREDICT_SOURCES predict.cc bw.cc wave.cc event.cc fitInfo.cc fitInfoDict.cc waveDict.cc)
set(MOMENTS_SOURCES moments.cc 3j.cc wave.cc event.cc waveDict.cc)

# define path to macros and settings
LIST( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

FIND_PACKAGE(ROOT REQUIRED COMPONENTS Minuit2 MathMore)

# include directories 
INCLUDE_DIRECTORIES( 
	BEFORE
	${ROOT_INCLUDE_DIRS})


MACRO (ROOT_GENERATE_DICTIONARY LIBNAME INFILES INCLUDE_DIRS_IN LINKDEF_FILE OUTFILE)
  set (INCLUDE_DIRS)
  foreach (_current_FILE ${INCLUDE_DIRS_IN})
    set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})   
  endforeach (_current_FILE ${INCLUDE_DIRS_IN})

  SET(EXTRA_DICT_ARGS "")
  IF(ROOT_found_version GREATER 59900)
    STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1_rdict.pcm" OUTFILE_PCM "${OUTFILE}")
    STRING(REGEX REPLACE "^(.*)Dict\\.(.*)$" "\\1.rootmap" OUTFILE_RMF "${OUTFILE}")
    SET (OUTFILES ${OUTFILE} ${OUTFILE_PCM} ${OUTFILE_RMF})
    SET(EXTRA_DICT_ARGS -inlineInputHeader -rmf ${OUTFILE_RMF} -rml ${LIBNAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
  ELSE()
    STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" OUTFILE_H "${OUTFILE}")
    SET (OUTFILES ${OUTFILE} ${OUTFILE_H})
  ENDIF()
  
  ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
    COMMAND ${ROOT_CINT_EXECUTABLE}
    ARGS -f ${OUTFILE} ${EXTRA_DICT_ARGS} -c ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE}
    DEPENDS ${INFILES} ${LINKDEF_FILE})

  STRING(REGEX REPLACE "/" "" OUTFILE_NO_SLASH "${OUTFILE}")
  ADD_CUSTOM_TARGET(generate_${OUTFILE_NO_SLASH} DEPENDS ${OUTFILE})

  IF(ROOT_found_version GREATER 59900)
    # PCM files and rootmap-files should be installed
    INSTALL(FILES ${OUTFILE_RMF} ${OUTFILE_PCM}
            DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
  ENDIF()
ENDMACRO (ROOT_GENERATE_DICTIONARY)

ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/wave.h"
                         "${CMAKE_CURRENT_SOURCE_DIR}"
			 "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/waveDict.cc")
ROOT_GENERATE_DICTIONARY("${CMAKE_SHARED_LIBRARY_PREFIX}${PROJECT_NAME}"
                         "${CMAKE_CURRENT_SOURCE_DIR}/fitInfo.h"
                         "${CMAKE_CURRENT_SOURCE_DIR}"
			 "${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h"
			 "${CMAKE_CURRENT_BINARY_DIR}/fitInfoDict.cc")

add_executable(myFit ${MYFIT_SOURCES})
add_executable(massDep ${MASSDEP_SOURCES})
add_executable(predict ${PREDICT_SOURCES})
add_executable(moments ${MOMENTS_SOURCES})