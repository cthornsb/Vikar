#vandmc tools.
add_executable(angleConvert angleConvert.cpp)
target_link_libraries(angleConvert VandmcStatic)

add_executable(eloss eloss.cpp)
target_link_libraries(eloss VandmcStatic)

add_executable(kindist Kindist.cpp)
target_link_libraries(kindist VandmcStatic)

#add_executable(lightOutput light.cpp)
#target_link_libraries(lightOutput VandmcStatic)

add_executable(relativistic relativistic.cpp)
target_link_libraries(relativistic VandmcStatic)

add_executable(testSetup testSetup.cpp)
target_link_libraries(testSetup VandmcStatic ${DICTIONARY_PREFIX}Static ${ROOT_LIBRARIES})

#add_executable(testViewer testViewer.cpp)
#target_link_libraries(testViewer VandmcStatic)

#add_executable(vandmcFront vandmcFront.cpp)
#target_link_libraries(vandmcFront VandmcStatic)

#add_executable(adIntegrator Integrator.cpp)
#target_link_libraries(adIntegrator VandmcStatic)

add_executable(kinematics Kinematics.cpp)
target_link_libraries(kinematics VandmcStatic)

add_executable(range range.cpp)
target_link_libraries(range VandmcStatic)

#add_executable(root2raw root2raw.cpp)
#target_link_libraries(root2raw VandmcStatic)

install(TARGETS angleConvert eloss kindist relativistic testSetup kinematics range DESTINATION bin)
