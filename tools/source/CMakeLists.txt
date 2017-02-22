#vandmc tools.
add_executable(angleConvert angleConvert.cpp)
target_link_libraries(angleConvert vandmcStatic)

add_executable(eloss eloss.cpp)
target_link_libraries(eloss vandmcStatic)

add_executable(kindist Kindist.cpp)
target_link_libraries(kindist vandmcStatic)

#add_executable(lightOutput light.cpp)
#target_link_libraries(lightOutput vandmcStatic)

add_executable(relativistic relativistic.cpp)
target_link_libraries(relativistic vandmcStatic)

add_executable(testSetup testSetup.cpp)
target_link_libraries(testSetup vandmcStatic ${DICTIONARY_PREFIX}Static ${ROOT_LIBRARIES})

#add_executable(testViewer testViewer.cpp)
#target_link_libraries(testViewer vandmcStatic)

#add_executable(vandmcFront vandmcFront.cpp)
#target_link_libraries(vandmcFront vandmcStatic)

#add_executable(adIntegrator Integrator.cpp)
#target_link_libraries(adIntegrator vandmcStatic)

add_executable(kinematics Kinematics.cpp)
target_link_libraries(kinematics vandmcStatic)

add_executable(range range.cpp)
target_link_libraries(range vandmcStatic)

#add_executable(root2raw root2raw.cpp)
#target_link_libraries(root2raw vandmcStatic)

install(TARGETS angleConvert eloss kindist relativistic testSetup kinematics range DESTINATION bin)
