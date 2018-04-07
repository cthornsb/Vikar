#Add the sources to the library.
add_library(ToolObjects OBJECT dataPack.cpp comConverter.cpp)

#Generate a static library.
add_library(VandmcToolStatic STATIC $<TARGET_OBJECTS:ToolObjects>)

#vandmc tools.

if(${BUILD_TOOLS_ELOSS})
	add_executable(eloss eloss.cpp)
	target_link_libraries(eloss VandmcStatic)
	install(TARGETS eloss DESTINATION bin)
endif()

if(${BUILD_TOOLS_KINDIST})
	add_executable(kindist Kindist.cpp)
	target_link_libraries(kindist VandmcStatic)
	install(TARGETS kindist DESTINATION bin)
endif()

if(${BUILD_TOOLS_RELATIVISTIC})
	add_executable(relativistic relativistic.cpp)
	target_link_libraries(relativistic VandmcStatic)
	install(TARGETS relativistic DESTINATION bin)
endif()

if(${BUILD_TOOLS_TESTSETUP})
	add_executable(testSetup testSetup.cpp)
	target_link_libraries(testSetup VandmcStatic VandmcToolStatic ${DICTIONARY_PREFIX}Static ${ROOT_LIBRARIES})
	install(TARGETS testSetup DESTINATION bin)
endif()

if(${BUILD_TOOLS_VANDLESOLIDANGLE})
	add_executable(vandleSolidAngle vandleSolidAngle.cpp)
	target_link_libraries(vandleSolidAngle VandmcStatic VandmcToolStatic ${DICTIONARY_PREFIX}Static ${ROOT_LIBRARIES})
	install(TARGETS vandleSolidAngle DESTINATION bin)
endif()

if(${BUILD_TOOLS_KINEMATICSGENERATOR})
	add_executable(kinematicsGenerator kinematicsGenerator.cpp)
	target_link_libraries(kinematicsGenerator VandmcStatic VandmcToolStatic ${ROOT_LIBRARIES})
	install(TARGETS kinematicsGenerator DESTINATION bin)
endif()

if(${BUILD_TOOLS_INTEGRATOR})
	add_executable(adIntegrator Integrator.cpp)
	target_link_libraries(adIntegrator VandmcStatic)
	install(TARGETS adIntegrator DESTINATION bin)
endif()

if(${BUILD_TOOLS_KINEMATICS})
	add_executable(kinematics Kinematics.cpp)
	target_link_libraries(kinematics VandmcStatic)
	install(TARGETS kinematics DESTINATION bin)
endif()

if(${BUILD_TOOLS_RANGE})
	add_executable(range range.cpp)
	target_link_libraries(range VandmcStatic)
	install(TARGETS range DESTINATION bin)
endif()

if(${BUILD_TOOLS_STRIPS})
	add_executable(strips strips.cpp)
	target_link_libraries(strips ${ROOT_LIBRARIES})
	install(TARGETS strips DESTINATION bin)
endif()

if(${BUILD_TOOLS_DETFILEMAKER})
	add_executable(detFileMaker detFileMaker.cpp)
	target_link_libraries(detFileMaker VandmcStatic)
	install(TARGETS detFileMaker DESTINATION bin)
endif()

# DEPRECATED 

#add_executable(angleConvert angleConvert.cpp)
#target_link_libraries(angleConvert VandmcStatic)

#add_executable(lightOutput light.cpp)
#target_link_libraries(lightOutput VandmcStatic)

#add_executable(testViewer testViewer.cpp)
#target_link_libraries(testViewer VandmcStatic)

#add_executable(vandmcFront vandmcFront.cpp)
#target_link_libraries(vandmcFront VandmcStatic)

#add_executable(root2raw root2raw.cpp)
#target_link_libraries(root2raw VandmcStatic)

