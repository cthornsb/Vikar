#Configure and install the module file
configure_file("modulefiles/vandmc" vandmc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/vandmc DESTINATION share/modulefiles/)
