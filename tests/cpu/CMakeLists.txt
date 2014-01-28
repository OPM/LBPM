# Add tests/executables that need to be compiled but not run
INSTALL_LBPM_EXE( lb2_Color_wia_mpi )


# Run the serial ConstrainedBubble inputs as a weekly test
CONFIGURE_FILE( ${LBPM_SOURCE_DIR}/example/ConstrainedBubble/Color.in ${CMAKE_CURRENT_BINARY_DIR}/Color.in COPYONLY )
CONFIGURE_FILE( ${LBPM_SOURCE_DIR}/example/ConstrainedBubble/Domain.in ${CMAKE_CURRENT_BINARY_DIR}/Domain.in COPYONLY )
ADD_LBPM_WEEKLY_TEST( lb2_Color_wia_mpi 1 )


