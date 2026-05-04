======================================
VTK XML output format
======================================

Due to several limitations in the support for VisIt files in ParaView, we added support for VTK XML (.vti) files, which are natively supported by ParaView.

This output format can be enabled by setting format = "vtk" in the Visualization section of the database file. The fields to be written to output can be chosen using the same parameters as in the VisIt output.

.. code:: c

   Visualization {
      format = "vtk"
      save_phase_field = true
      save_pressure    = false
   }

The interval for writing output files is defined by the ``visualization_interval`` parameter in the ``Analysis`` section.

LBPM also generates and updates a ``LBM.pvd`` file, which allows ParaView to open all .vti files generated during the current simulation.
