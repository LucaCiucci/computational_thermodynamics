#pragma once

#include "simulation.h"

constexpr double Epsilon = 0.000000001;

int sign(double x);/* {
	if (x == 0) return 0;
	if (x > 0) return 1;
	return -1;
}*/

double fRand(double fMin, double fMax);/*
{
	
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
	
}*/

double abs(Vector3);// TODO
double abs(Vector2);// TODO