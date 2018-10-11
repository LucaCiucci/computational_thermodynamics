
#include <cstdlib>
#include "mathFunctions.h"

int sign(double x) {
	if (x == 0) return 0;
	if (x > 0) return 1;
	return -1;
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}