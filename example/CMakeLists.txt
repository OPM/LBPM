# Copy the examples to the install folder
INSTALL_EXAMPLE( Bubble )
INSTALL_EXAMPLE( ConstrainedBubble )
INSTALL_EXAMPLE( Piston )
INSTALL_EXAMPLE( Sph1896 )

# Create unit tests for each example
#   TEST_EXAMPLE( example exe N_procs )
TEST_EXAMPLE( Bubble            TestBubble        1  )
TEST_EXAMPLE( Sph1896           lb2_Color_wia_mpi 27 )
TEST_EXAMPLE( ConstrainedBubble lb2_Color_wia_mpi 1  )



