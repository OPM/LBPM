======================================
Color model -- Image Sequence Protocol
======================================

The image sequence protocol is designed to perform a set steady-state
simulations based on a sequence of 3D (8-bit) images provided by the user.
The images might be the output of a previous LBPM simulation, a sequence of
(segmented) experimental data, or data generated from a custom routine.
The image sequence protocol will apply the same set of flow conditions
to all images in the sequence. This means

* the image labels and any associated properties are the same
* the external boundary conditions are the same
* the physical simulation parameters are the same

The image sequence protocol does not set boundary conditions by default.
It is up to the user to determine the flow condition, with the understanding
that the same set of will be applied to each image in the sequence.

To enable the image sequence protocol, the following keys should be set
within the ``Color`` section of the input database

.. code-block:: bash

   protocol = "image sequence"
   image_sequence = "image1.raw", "image2.raw"
    
