set(RENDERER_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(RENDERER_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

set(QT_DEFINES "-DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED")
set(QT_INCPATH "-I/usr/include/qt4 -I/usr/include/qt4/QtGui")

set(CAMERA_HEADER ${RENDERER_SOURCE_DIR}/camera.h)
set(CAMERA_UI_HEADER ${RENDERER_SOURCE_DIR}/camera.ui)
set(MOC_FILENAME ${RENDERER_BUILD_DIR}/moc_camera.cpp)
set(UIC_FILENAME ${RENDERER_BUILD_DIR}/ui_camera.h)

#Use moc to generate the moc_camera source file.
add_custom_target(MocFile ALL DEPENDS ${MOC_FILENAME})
add_custom_command(OUTPUT ${MOC_FILENAME} COMMAND ${MOC_EXECUTABLE} ${QT_DEFINES} ${QT_INCPATH} ${CAMERA_HEADER} -o ${MOC_FILENAME} DEPENDS ${CAMERA_HEADER} VERBATIM)

#Use uic to generate the ui_camera header file.
add_custom_target(UicFile ALL DEPENDS ${UIC_FILENAME})
add_custom_command(OUTPUT ${UIC_FILENAME} COMMAND ${UIC_EXECUTABLE} ${CAMERA_UI_HEADER} -o ${UIC_FILENAME} DEPENDS ${CAMERA_UI_HEADER} VERBATIM)

#Build the renderer executable.
add_executable(renderer main.cpp camera.cpp ${MOC_FILENAME})
target_include_directories(renderer PUBLIC ${QT_INC_DIR} ${QT_GUI_INC_DIR} ${RENDERER_BUILD_DIR})
target_link_libraries(renderer VandmcStatic ${QT_GUI_LIB} ${QT_CORE_LIB})
install(TARGETS renderer DESTINATION bin)
