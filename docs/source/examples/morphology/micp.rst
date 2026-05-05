************************************
MICP Morphology Simulator
************************************

To expand and improve the current morphological tools in LBPM, new
morphology pre-processors and simulators were added. These changes
encompass the implementation of state-of-the-art morphology algorithms,
as well as an additional MICP simulator in which the invading non-wetting
phase is injected through all faces of the sample.

.. code:: c

   Domain {
      Filename = "crop_bt_101x112x88_uint8.raw"
      N = 101, 112, 88            // domain size
      n = 101, 112, 88
      nproc = 1, 1, 1
      ReadValues = 0, 1, 2
      WriteValues = 0, 1, 2
      voxel_length = 1
      BC = 0
   }

   FM {
      Diameters = 1, 30, 1        // start, end, step
      SaveImage = true
      direction = "surround"      // axis of intrusion
      protocol = "micp"
   }
