set(ROOT_DICT_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

set(STRUCTURES_PREFIX "vandmcStructures")
set(STRUCTURES_SOURCE "${ROOT_DICT_BUILD_DIR}/${STRUCTURES_PREFIX}.cpp")
set(STRUCTURES_HEADER "${ROOT_DICT_BUILD_DIR}/${STRUCTURES_PREFIX}.hpp")
set(LINKDEF_FILE "${ROOT_DICT_BUILD_DIR}/LinkDef.h")

set(COPIED_STRUCTURES_HEADER "${TOP_DIRECTORY}/include/vandmcStructures.hpp")

set(CINT_OUTPUT_FILE "${ROOT_DICT_BUILD_DIR}/${DICTIONARY_PREFIX}.cpp")
set(CINT_PCM_FILE "${ROOT_DICT_BUILD_DIR}/${DICTIONARY_PREFIX}_rdict.pcm")

add_custom_target(GenerateDict ALL DEPENDS ${STRUCTURES_SOURCE} ${STRUCTURES_HEADER} ${COPIED_STRUCTURES_HEADER} ${LINKDEF_FILE} ${CINT_OUTPUT_FILE} ${CINT_PCM_FILE})

#Use rcbuild to generate data structure source files.
add_custom_command(OUTPUT ${STRUCTURES_SOURCE} ${STRUCTURES_HEADER} ${LINKDEF_FILE} COMMAND ${RCBUILD_EXECUTABLE} -i ${DEF_STRUCT_FILE} -d ${ROOT_DICT_BUILD_DIR} -P ${STRUCTURES_PREFIX} DEPENDS ${DEF_STRUCT_FILE} VERBATIM)

#Install the new Structures file.
install(FILES ${STRUCTURES_HEADER} DESTINATION include)

#Copy the new Structures file to the top-level include directory.
add_custom_command(OUTPUT ${COPIED_STRUCTURES_HEADER} COMMAND cp ${STRUCTURES_HEADER} ${COPIED_STRUCTURES_HEADER} DEPENDS ${STRUCTURES_HEADER} VERBATIM)

#Use rootcint to generate dictionary source file.
add_custom_command(OUTPUT ${CINT_OUTPUT_FILE} ${CINT_PCM_FILE} COMMAND ${ROOTCINT_EXECUTABLE} -f ${CINT_OUTPUT_FILE} -c ${STRUCTURES_HEADER} ${LINKDEF_FILE} DEPENDS ${STRUCTURES_HEADER} ${LINKDEF_FILE} VERBATIM)

#Set the scan sources that we will make a lib out of.
set(DictSources ${STRUCTURES_SOURCE} ${CINT_OUTPUT_FILE})

#Add the sources to the library.
add_library(DictObjects OBJECT ${DictSources})

#Create static library.
add_library(${DICTIONARY_PREFIX}Static STATIC $<TARGET_OBJECTS:DictObjects>)

# Build shared libs
if(${BUILD_SHARED})
	#Create shared library.
	add_library(${DICTIONARY_PREFIX} SHARED $<TARGET_OBJECTS:DictObjects>)
	install(TARGETS ${DICTIONARY_PREFIX} DESTINATION lib)
endif(${BUILD_SHARED})


