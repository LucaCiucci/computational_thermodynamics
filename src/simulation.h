#pragma once
#ifndef SIMULATION_H
#define SIMULATION_H


// Math.h - STD math Library
#include <math.h>

#include <vector>

#include <string>

#include <iostream>

// TODO rearrange

struct Vector2
{
	// Default Constructor
	Vector2()
	{
		X = 0.0f;
		Y = 0.0f;
	}
	// Variable Set Constructor
	Vector2(double X_, double Y_)
	{
		X = X_;
		Y = Y_;
	}
	// Bool Equals Operator Overload
	bool operator==(const Vector2& other) const
	{
		return (this->X == other.X && this->Y == other.Y);
	}
	// Bool Not Equals Operator Overload
	bool operator!=(const Vector2& other) const
	{
		return !(this->X == other.X && this->Y == other.Y);
	}
	// Addition Operator Overload
	Vector2 operator+(const Vector2& right) const
	{
		return Vector2(this->X + right.X, this->Y + right.Y);
	}
	// Subtraction Operator Overload
	Vector2 operator-(const Vector2& right) const
	{
		return Vector2(this->X - right.X, this->Y - right.Y);
	}
	// double Multiplication Operator Overload
	Vector2 operator*(const double& other) const
	{
		return Vector2(this->X *other, this->Y * other);
	}

	// Scalar Product Operator Overload
	Vector2 operator*(const Vector2& right) const
	{
		return Vector2(this->X * right.X, this->Y * right.Y);
	}
	// Vectorial Product Operator Overload
	double operator%(const Vector2& right) const
	{
		return double(this->X*right.Y - this->Y*right.X);
	}

	// Positional Variables
	double X;
	double Y;
};

struct Vector3
{
	// Default Constructor
	Vector3()
	{
		X = 0.0f;
		Y = 0.0f;
		Z = 0.0f;
	}
	// Variable Set Constructor
	Vector3(double X_, double Y_, double Z_)
	{
		X = X_;
		Y = Y_;
		Z = Z_;
	}
	// Bool Equals Operator Overload
	bool operator==(const Vector3& other) const
	{
		return (this->X == other.X && this->Y == other.Y && this->Z == other.Z);
	}
	// Bool Not Equals Operator Overload
	bool operator!=(const Vector3& other) const
	{
		return !(this->X == other.X && this->Y == other.Y && this->Z == other.Z);
	}
	// Addition Operator Overload
	Vector3 operator+(const Vector3& right) const
	{
		return Vector3(this->X + right.X, this->Y + right.Y, this->Z + right.Z);
	}
	// Subtraction Operator Overload
	Vector3 operator-(const Vector3& right) const
	{
		return Vector3(this->X - right.X, this->Y - right.Y, this->Z - right.Z);
	}
	// double Multiplication Operator Overload
	Vector3 operator*(const double& other) const
	{
		return Vector3(this->X * other, this->Y * other, this->Z * other);
	}
	// double Division Operator Overload
	Vector3 operator/(const double& other) const
	{
		return Vector3(this->X / other, this->Y / other, this->Z / other);
	}
	
	// Scalar Product Operator Overload
	Vector3 operator*(const Vector3& right) const
	{
		return Vector3(this->X * right.X, this->Y * right.Y, this->Z * right.Z);
	}
	// Vectorial Product Operator Overload
	Vector3 operator%(const Vector3& right) const
	{
		return Vector3(this->Y*right.Z - this->Z*right.Y, this->Z*right.X - this->X*right.Z, this->X*right.Y - this->Y*right.X);
	}

	// Positional Variables
	double X;
	double Y;
	double Z;
};

struct Triangle3 {
	Vector3 v1, v2, v3;
};

struct Triangle2 {
	Vector2 v1, v2, v3;
};

struct Box {
	Vector3 min;
	Vector3 max;
};

struct Vertex {
	// Position Vector
	Vector3 Position;

	// Normal Vector
	Vector3 Normal;
};

struct Plane {
	double a, b, c, d;// ax+by+cz+d=0
};

struct Mesh {
	// Mesh Name
	std::string MeshName;
	// Vertex List
	std::vector<Vertex> Vertices;
	// Index List
	std::vector<unsigned int> Indices;

	//
	std::vector<Plane> plane;
};

struct Volume {
	Box box = { {0, 0, 0}, {0, 0, 0} };
	bool containsGas = true;
	bool containsSB = false;
	bool isActive = true;
	double deltaPotential = 0;
	Mesh mesh;
	/*it would e grate if we could assign a variable potential in the form
	a*x + b*y + c*z + d, but it would add more unnecessary complexity to the program and it will slow down the simulation a lot,
	maybe I'll do this in a future version
	*/
};

struct GasParticle {
	double potential = 0;
	Vector3 position;
	Vector3 speed;
};



struct SimulationSettings {
	double dt = 0.01;// max step dt
	double gasRadius = 0.001;// particles will collide when distance = 2 * gasRadius
	double sbRadius = 0.001;// radius of the softbody particles
	double finalTime = -1;// simulation end time, -1 to disable
	int maxStepNumber = 10;// maximum number of simulation steps
	double gasTemperature = 1;// average of points' kinetic energy
	int gasPointsNumber = 10;
    // TODO expand if necessary...
};

enum class EventType {
	gasGasCollision,
	gasMeshCollision,
	gasSbCollision,
	sbSbCollision,
	none
};

struct SimEvent {
	double relTime = 0.0;
	EventType eventType = EventType::none;
	int index1;// index of the first object (gas or sb)
	int index2;// index of the second object (gas or volume)
	int index3;// index of the subobjec (i.e. triangle index in the index2 volume)
	Vector3 position1;
	Vector3 position2;
};






//-------------------------------------------------------------------------------------------------

class Simulation { // TODO reorganize
public:
	//Simulation(int gasParticleNumber, double temperature);// TODO to complete

	bool addVolume(Volume);
	void ciao1(void) const { std::cout << "\nCiao dalla SIMULAZIONE!!!\n"; }// 1!!!!!!!!!!!!!!!!!!!!!!!1
	int getVolumesNumber(void) const;
	void setGasParticleNumber(int);
	int getGasParticleNumber(void) const;
	void setInitialTemperature(double);
	double getInitialTemperature(void) const;
	bool createGasParticles(void);
	bool performOneStep(void);

	//TODO delete these functions...
	void printVertex(int) const;
	void printPlanes(int) const;
	void printPoints(void) const;
	void test(void) const;

private:
	bool isInsideVolume(Vector3, int) const;// (point, volume index)
	bool isInside2dTriangle(Triangle2, Vector2) const;// TODO generalize for 3d triangle, too
	bool isInside3dTriangle(Triangle3, Vector3) const;
	Triangle3 triangle2dTo3d(Triangle2) const;
	Vector3 vector2dTo3d(Vector2) const;
	Box boundingBox(Volume) const;
	Box boundingBox(int) const;// (volume index)
	std::vector<int> volumesContainingVector(Vector3) const;
	double potential(Vector3) const;
	bool isInAllowedVolume(Vector3) const;
	bool simulateInTimeInterval(double);// (dt)
	SimEvent findFirstEvent(double) const;
	SimEvent findGasMeshEvent(double) const;
	SimEvent findFirstContainerCollision(int, double) const;
	SimEvent findFirstTriangleCollision(int, int, double) const;
	SimEvent findTriangleCollision(int, int, int, double) const;
	
	
	double temperature = 1;
	double GasParticleNumber = 10;
	SimulationSettings simSettings;
	std::vector<Volume> volumes;
	std::vector<GasParticle> gasParticles;
    
};


#endif // SIMULATION_H
