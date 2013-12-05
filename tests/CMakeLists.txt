# Copy files for the tests
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cylindertest ${CMAKE_CURRENT_BINARY_DIR}/cylindertest COPYONLY )


# Add the tests
ADD_LBPM_TEST( pmmc_cylinder )
ADD_LBPM_TEST( TestCylinderAreas )
ADD_LBPM_TEST( TestSphereCurvature )



