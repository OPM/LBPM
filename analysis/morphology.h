// Morphological opening routine
#include "common/Array.h"
#include "common/Domain.h"

double MorphOpen(DoubleArray &SignDist, char *id, std::shared_ptr<Domain> Dm, double VoidFraction);
double MorphGrow(DoubleArray &BoundaryDist, DoubleArray &Dist, Array<char> &id, std::shared_ptr<Domain> Dm, double TargetVol);