
#include <cstdlib>
#include "mathFunctions.h"
#include <math.h>

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

double abs(Vector3 vector)
{
	return sqrt( vector.X * vector.X + vector.Y * vector.Y + vector.Z * vector.Z );
}

double abs(Vector2 vector)
{
	return sqrt( vector.X * vector.X + vector.Y * vector.Y );
}