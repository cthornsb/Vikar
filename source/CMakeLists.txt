#Set the scan sources that we will make a lib out of.
set(CoreSources detectors.cpp geometry.cpp kindeux.cpp materials.cpp vandmc_core.cpp)

#Add the sources to the library.
add_library(CoreObjects OBJECT ${CoreSources})

#Generate a static library.
add_library(VandmcStatic STATIC $<TARGET_OBJECTS:CoreObjects>)

#Build simpleScan executable.
add_executable(vandmc vandmc.cpp)
target_link_libraries(vandmc VandmcStatic ${DICTIONARY_PREFIX}Static ${SimpleScan_SCAN_LIB} ${SimpleScan_OPT_LIB} ${ROOT_LIBRARIES})
install(TARGETS vandmc DESTINATION bin)

# Build shared libs
if(${BUILD_SHARED})
	add_library(Vandmc SHARED $<TARGET_OBJECTS:CoreObjects>)
	install(TARGETS Vandmc DESTINATION lib)
endif(${BUILD_SHARED})
