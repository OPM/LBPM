# visit -debug 1 -valgrind engine_ser -cli -s /projects/JamesClure/LBPM-WIA/visit/testVisit.py


OpenDatabase("localhost:/projects/JamesClure/build/debug/tests/summary.LBM", 0)
#AddPlot("Mesh", "pointmesh", 1, 1)
AddPlot("Mesh", "trilist", 1, 1)
DrawPlots()
#quit()


