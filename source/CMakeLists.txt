#Set the scan sources that we will make a lib out of.
set(CoreSources detectors.cpp geometry.cpp kindeux.cpp materials.cpp vandmc_core.cpp)

#Add the sources to the library.
add_library(CoreObjects OBJECT ${CoreSources})

#Generate a static library.
add_library(VandmcStatic STATIC $<TARGET_OBJECTS:CoreObjects>)

#Build simpleScan executable.
add_executable(vandmc vandmc.cpp)
target_link_libraries(vandmc VandmcStatic ${DICTIONARY_PREFIX}Static ${ROOT_LIBRARIES})
install(TARGETS vandmc DESTINATION bin)
