//#include "OBJ_Loader.h"
#include "simulation.h"
#include "mathFunctions.h"
//#include <random>

//-------------------------------------------------------------------------------------------------

int Simulation::getVolumesNumber(void) const { return volumes.size(); }

void Simulation::setGasParticleNumber(int _gasParticleNumber) { GasParticleNumber = _gasParticleNumber; }

int Simulation::getGasParticleNumber(void) const { return GasParticleNumber; }

void Simulation::setInitialTemperature(double _temperature) { temperature = _temperature; }

double Simulation::getInitialTemperature(void) const { return temperature; }

//-------------------------------------------------------------------------------------------------

bool Simulation::addVolume(Volume volume)
{
	Plane plane;
	Vector3 p1, p2, p3;
	Vector3 v1, v2,V;
	double Vm;
	for (int i = 0; i < volume.mesh.Indices.size(); i += 3)
	{
		p1 = volume.mesh.Vertices[volume.mesh.Indices[i]].Position;
		p2 = volume.mesh.Vertices[volume.mesh.Indices[i+1]].Position;
		p3 = volume.mesh.Vertices[volume.mesh.Indices[i+2]].Position;
		
		v1 = p2 - p1;
		v2 = p3 - p1;
		V = v1 % v2;// triangle area vector
		Vm = sqrt(V.X*V.X + V.Y*V.Y + V.Z*V.Z); // V's lenght

		plane.a = V.X / Vm;// coefficients for the equation of the plane containing the triangle
		plane.b = V.Y / Vm;// a*x + b*y + c*z + d = 0
		plane.c = V.Z / Vm;
		plane.d =
			- plane.a * p1.X
			- plane.b * p1.Y
			- plane.c * p1.Z;

		volume.mesh.plane.push_back(plane);
	}
	volume.box = boundingBox(volume);
	volumes.push_back(volume);
	return false;// TODO
}

bool Simulation::createGasParticles(void)// TODO avoid collisions!
{
	bool someVolumes = false;
	Box box;
	bool n1 = false;
	// find the box containing all the bounding boxes
	for (int i = 0; i < volumes.size(); i++)
	{
		if (volumes[i].containsGas) {
			if (!n1)
			{
					n1 = true;
					box = volumes[i].box;
					someVolumes = true;
			}
			box.min.X = fmin(box.min.X, volumes[i].box.min.X);
			box.min.Y = fmin(box.min.Y, volumes[i].box.min.Y);
			box.min.Z = fmin(box.min.Z, volumes[i].box.min.Z);

			box.max.X = fmax(box.max.X, volumes[i].box.max.X);
			box.max.Y = fmax(box.max.Y, volumes[i].box.max.Y);
			box.max.Z = fmax(box.max.Z, volumes[i].box.max.Z);
		}
	}
	if (!someVolumes) return false;// there are no volumes that should contain gas

	GasParticle newPoint;
	double rand;
	// creates points in a random position inside the right volumes
	for (int i = 0; i < GasParticleNumber; i++)
	{
		// tries random points until a contained point in the correct position is found
		do {
			newPoint.position.X = fRand(box.min.X, box.max.X);
			newPoint.position.Y = fRand(box.min.Y, box.max.Y);
			newPoint.position.Z = fRand(box.min.Z, box.max.Z);

			newPoint.potential = potential(newPoint.position);
		} while ((newPoint.potential >= 0) || !isInAllowedVolume(newPoint.position));// TODO qui c'è un errore!!!!!!!!!!!!! o forse è stato corretto o forse no?!?!?!?
		// assign random speed
		rand = fRand(- temperature * 2.0 / 3.0, temperature * 2.0 / 3.0); // energy along x direction (signed)
		newPoint.speed.X = (double)sign(rand) * (double)sqrt (2.0 * abs(rand));
		rand = fRand(-temperature * 2.0 / 3.0, temperature * 2.0 / 3.0);
		newPoint.speed.Y = (double)sign(rand) * (double)sqrt(2.0 * abs(rand));
		rand = fRand(-temperature * 2.0 / 3.0, temperature * 2.0 / 3.0);
		newPoint.speed.Z = (double)sign(rand) * (double)sqrt(2.0 * abs(rand));
		// add to the array of gas particles
		gasParticles.push_back(newPoint);
		/*
		std::cout << "eccomi:::::::" << isInAllowedVolume({ -0.997497, 0.997918, -0.613392 }) << std::endl;////////////////////
		std::cout << "eccolo-------" << isInAllowedVolume(newPoint.position) << std::endl;////////////////////
		*/
	}
	return true;
}

// run simulation for one step, then pause
bool Simulation::performOneStep(void)// TODO
{
	/*
	call the function that simulates from now to the now+dt time
	*/
	return false;
}

// TODO
// TODO perchè non va bene?!?!?!
/*Simulation::Simulation(int gasParticleNumber, double temperature)
{

}*/


//---------------------------------PRIVATE---------------------------------------------------------

// find the bounding for the provided volume
Box Simulation::boundingBox(Volume currVolume) const
{
	Box box;
	box.max = box.min = currVolume.mesh.Vertices[0].Position;
	for (int i = 1; i < currVolume.mesh.Vertices.size(); i++) {
		box.min.X = fmin(box.min.X, currVolume.mesh.Vertices[i].Position.X);
		box.min.Y = fmin(box.min.Y, currVolume.mesh.Vertices[i].Position.Y);
		box.min.Z = fmin(box.min.Z, currVolume.mesh.Vertices[i].Position.Z);

		box.max.X = fmax(box.max.X, currVolume.mesh.Vertices[i].Position.X);
		box.max.Y = fmax(box.max.Y, currVolume.mesh.Vertices[i].Position.Y);
		box.max.Z = fmax(box.max.Z, currVolume.mesh.Vertices[i].Position.Z);
		
	}
	return box;
}

//the same but accepting the index of the volume
Box Simulation::boundingBox(int volumeIndex) const
{
	Box _boundingBox = { {0, 0, 0,}, {0, 0, 0} };
	if (volumeIndex < volumes.size())
		_boundingBox = boundingBox(volumes[volumeIndex]);
	else
		std::cout << "OUT OF RANGE!" << std::endl;
	return _boundingBox;
}

// list the indexes of the box containing the point
std::vector<int> Simulation::volumesContainingVector(Vector3 position_vector) const
{
	std::vector<int> indexes;
	for (int i = 0; i < volumes.size(); i++)
	{
		if (isInsideVolume(position_vector, i))
			indexes.push_back(i);
	}
	return indexes;
}

//get potential at provided point
double Simulation::potential(Vector3 position_vector) const
{
	std::vector<int> list;
	double potential = 0;
	list = volumesContainingVector(position_vector);
	for (int i = 0; i < list.size(); i++)
	{
		potential += volumes[list[i]].deltaPotential;
	}
	return potential;
}

// is in a point where gas should be created
bool Simulation::isInAllowedVolume(Vector3 position_vector) const
{
	for (int i = 0; i < volumes.size(); i++)
		if (volumes[i].containsGas && isInsideVolume(position_vector, i))
			return true;
	return false;
}

// point is inside the provided volume
bool Simulation::isInsideVolume(Vector3 point, int volumeIndex) const // TODO optimize (avoid unnecessary vars, too)
{
	//example: isInsideVolume({ 1.0, 16.0, 1.0 }, 0)

	// find the box containing the volume
	Volume currVolume = volumes[volumeIndex];
	Box box = volumes[volumeIndex].box;
	/*if (box.max == box.min) {
		boundingBox(volumeIndex);
		box = volumes[volumeIndex].box;
	}*/
	// if the point is outside the box,
		// then it is outside the volume
	if (point.X < box.min.X) return false;
	if (point.Y < box.min.Y) return false;
	if (point.Z < box.min.Z) return false;

	if (point.X > box.max.X) return false;
	if (point.Y > box.max.Y) return false;
	if (point.Z > box.max.Z) return false;
	// else
	// RAY-CASTING
	// trace a ray ( for example {∀t∈ℝ0+ : x = Xp+t; y = Yp; z = Zp})
	/*
	The easiest way to do this, I suppose, is to look at the system from one side.
	check how many triangle the point is inside of, to do this we check if the point is above or below
	the plane containing the triangle (ignore planes perpendicular to the view plane), and check if the
	projection of the plane is inside the projection of the triangle (using area vectors)
	*/
	// for every triangle
	Plane currPlane;
	int counter = 0;
	double delta;
	Triangle2 triangle;
	for (int i = 0; i < volumes[volumeIndex].mesh.plane.size(); i++)
	{
		currPlane = volumes[volumeIndex].mesh.plane[i];// to make coding easier
		// this was the general idea -> // find intersection, if it is inside the triangle, increase counter (think about system with no solution or indetetrmined too! (solved) )
		// suppose we are looking from the point fo Z
		// if c = 0, this means the plane in perpendicular to our view plane, ne do not neeed to consider this triangle
		if (currPlane.c != 0)
		{
			// else, f(x,y,z) = ax+b+cz+d=0. the point is above the plane if (delta := ( sign(c) * sign(f(xp,yp,zp)) ) ) = +1, otherwise the point is below, if 0 the point is on the plane
			delta = sign(currPlane.c) * sign(currPlane.a * point.X + currPlane.b * point.Y + currPlane.c * point.Z + currPlane.d);
			triangle =
			{
				{ volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 0]].Position.X,volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 0]].Position.Y },
				{ volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 1]].Position.X,volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 1]].Position.Y },
				{ volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 2]].Position.X,volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * i + 2]].Position.Y }
			};
			// if is inside the triangle and delta > 0 (not equal, we consider points on the triangles as "outside" to simplify the creation of random points (but this is wrong?!)), then counter++
			if (isInside2dTriangle(triangle, { point.X, point.Y }) && (delta > 0)) counter++;
		}
	}
	// if counter % 2 = 1, point is inside volume
	if (counter % 2 == 1)
		return true;
	return false;
}

// used inside the isInsideVolume() funcion
bool Simulation::isInside2dTriangle(Triangle2 triangle, Vector2 point) const
{
	double triangleArea = abs( (triangle.v2 - triangle.v1) % (triangle.v3 - triangle.v1) );
	double subArea1 = abs((point - triangle.v1) % (point - triangle.v2));
	double subArea2 = abs((point - triangle.v2) % (point - triangle.v3));
	double subArea3 = abs((point - triangle.v3) % (point - triangle.v1));

	//if the areas are equal, then th point is inside the triangle
	if ((subArea1 + subArea2 + subArea3 - triangleArea) <= Epsilon) // approximation may cause areas not to coincide perfectly
		return true;
	return false;
}

// TODO
// simulate over time interval dt
bool Simulation::simulateInTimeInterval(double dt)// TODO
{
	
	// create a struct that identifies the first event and its informations (simEvent)
	SimEvent firstSimEvent;
	// call the funcion to find this event
	//firstSimEvent = findFirstEvent();
	// if there is an event, then perform first event.
	// apply space and (in case) speed variations, then if the time of the first event < dt,
	// then call this function in the new interval = dt - first event time
	
	
	return false;
}

// TODO
// TODO dt e NON simSettings.dt!!!!! anche nelle altre funzioni!!!!!!!
SimEvent Simulation::findFirstEvent(double dt) const
{
	SimEvent firstSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	/*// TODO
	newSimEvent = findGasMeshEvent();
	if (newSimEvent.relTime < firstSimEvent)
		firstSimEvent = newSimEvent;
	*/
	return SimEvent();
}

// TODO
SimEvent Simulation::findGasMeshEvent(double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	for (int i = 0; i < gasParticles.size(); i++)
	{
		/*
		newSimEvent = findFirstContainerCollision(i);
		if (newSimEvent.relTime < firstSimEvent)
			firstSimEvent = newSimEvent;
		*/
	}
	return SimEvent();
}


// TODO
SimEvent Simulation::findFirstContainerCollision(int gasPointIndex, double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	for (int i = 0; i < volumes.size(); i++)
	{
		/*
		newSimEvent = findFirstTriangleCollision(gasPointIndex, i);
		if (newSimEvent.relTime < firstSimEvent)
			firstSimEvent = newSimEvent;
		*/
	}
	return SimEvent();
}

// TODO
SimEvent Simulation::findFirstTriangleCollision(int gasPointIndex, int volumeIndex, double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	for (int i = 0; i < volumes[volumeIndex].mesh.plane.size(); i++) {
		/*
		newSimEvent = findTriangleCollision(gasPointIndex, volumesIndex, i);
		if (newSimEvent.relTime < firstSimEvent)
			firstSimEvent = newSimEvent;
		*/
	}
	return SimEvent();
}

SimEvent Simulation::findTriangleCollision(int gasPointIndex, int volumeIndex, int triangleIndex, double dt) const
{
	GasParticle particle = gasParticles[gasPointIndex];
	GasParticle particleF = particle;//particle at the final position, excluding collisions
	Plane plane = volumes[volumeIndex].mesh.plane[triangleIndex];
	// TODO add triangle object
	particleF.position = particleF.position + particleF.speed * dt;
	// if it crosses the plane
	if (
		sign( particle.position.X * plane.a + particle.position.Y * plane.b + particle.position.Z * plane.c + plane.d )
		!=
		sign(particleF.position.X * plane.a + particleF.position.Y * plane.b + particleF.position.Z * plane.c + plane.d)
		)
	{
		// then find when...
		double collisionRelativeTime = -( particle.position.X * plane.a + particle.position.Y * plane.b + particle.position.Z * plane.c + plane.d )
			/ ( particle.speed.X * plane.a + particle.speed.Y * plane.b + particle.speed.Z * plane.c );
		// if it appens before dt! (or it would be useless to calculate the remaining collision parameters)
		if (collisionRelativeTime <= dt)// TODO < o <=???
			return SimEvent();
		// ...and where
		Vector3 collisionPosition = particle.position + particle.speed * collisionRelativeTime;
		// then check if it is inside the triangle
		// TODO CONTINUA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// then calculate collion parameters
	}
	return SimEvent();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Simulation::printVertex(int volumeIndex) const
{
	if ( volumes.size()  > volumeIndex)
	{
		for (int i = 0; i < volumes[volumeIndex].mesh.Vertices.size(); i++) {
			std::cout << "P:\t" << volumes[volumeIndex].mesh.Vertices[i].Position.X << " \t" << volumes[volumeIndex].mesh.Vertices[i].Position.Y << " \t" << volumes[volumeIndex].mesh.Vertices[i].Position.Z << std::endl;
		}/*
		for (int i = 0; i < volumes[volumeIndex].mesh.Indices.size(); i++)
			std::cout << volumes[volumeIndex].mesh.Indices[i] << std::endl;*///TODO delete, only for test
	} else 
	{
		std::cout << "OUT OF RANGE!" << std::endl;
	}
	return;
}

void Simulation::printPlanes(int volumeIndex) const
{
	if (volumes.size() > volumeIndex)
	{
		std::cout << "\n\n\n";
		for (int i = 0; i < volumes[volumeIndex].mesh.plane.size(); i++) {
			std::cout
				<< "Plane:\t a = " << volumes[volumeIndex].mesh.plane[i].a
				<< "\t b = " << volumes[volumeIndex].mesh.plane[i].b
				<< "\t c = " << volumes[volumeIndex].mesh.plane[i].c
				<< "\t d = " << volumes[volumeIndex].mesh.plane[i].d
				<< std::endl;
		}/*
		for (int i = 0; i < volumes[volumeIndex].mesh.Indices.size(); i++)
			std::cout << volumes[volumeIndex].mesh.Indices[i] << std::endl;*///TODO delete, only for test
	}
	else
	{
		std::cout << "OUT OF RANGE!" << std::endl;
	}
	return;
}

void Simulation::printPoints(void) const
{
	std::cout << "\n\n\n";
	double energySum = 0, energy;
	for (int i = 0; i < gasParticles.size(); i++) {
		energy = gasParticles[i].speed.X * gasParticles[i].speed.X / 2 +
			gasParticles[i].speed.Y * gasParticles[i].speed.Y / 2 +
			gasParticles[i].speed.Z * gasParticles[i].speed.Z / 2;
		std::cout << "Point #" << i << "\tposition: x = "
			<< gasParticles[i].position.X << " y = "
			<< gasParticles[i].position.Y << " z = "
			<< gasParticles[i].position.Z << "       \tspeed: x = "
			<< gasParticles[i].speed.X << " y = "
			<< gasParticles[i].speed.Y << " z = "
			<< gasParticles[i].speed.Z << "              \tenergy: "
			<< energy << std::endl;
		energySum += energy;
	}
	std::cout << "\nAverage energy: " << energySum / gasParticles.size() << std::endl;
}

void Simulation::test(void) const
{
	std::cout << "\n\n" << isInsideVolume({ 1.01, 1.0, 1.0 }, 0) << std::endl;//non funzionaaaaaaaa (se capita su una linea (entro tolleranza epsilon su area non funziona))
	std::cout << isInsideVolume({ 1.0, 6.0, 1.0 }, 0) << std::endl;
	std::cout << isInsideVolume({ 1.0, 16.0, 1.0 }, 0) << std::endl;
	std::cout << isInsideVolume({ 1.0, -6.0, 1.0 }, 0) << std::endl;
	//std::cout           << isInsideVolume({ 0.01, 2.411, 0.012 }, 2) << std::endl;
	std::cout << isInAllowedVolume({ -0.997497, 0.997918, -0.613392 }) << std::endl;
	std::cout << isInsideVolume({ 1.9, 1.92, 1.93 }, 0) << std::endl;
	std::cout << "box x:: " << volumes[0].box.min.X << "box x:: " << volumes[0].box.min.Y << "box z:: " << volumes[0].box.min.Z << std::endl;
}