//#include "OBJ_Loader.h"
#include "simulation.h"
#include "mathFunctions.h"
#include <Windows.h>
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
		Vm = abs(V); // V's lenght

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
	simulateInTimeInterval(simSettings.dt);
	return false;
}

double Simulation::getDt(void) const { return simSettings.dt; }
bool Simulation::setDt(double dt) { simSettings.dt = dt; return true; }
;

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
	// TODO this commented part can be deleted, looks like the usage of isInside3dTriangle()
	// works fine, but check another time before deleting
	/*
	double triangleArea = abs( (triangle.v2 - triangle.v1) % (triangle.v3 - triangle.v1) );
	double subArea1 = abs((point - triangle.v1) % (point - triangle.v2));
	double subArea2 = abs((point - triangle.v2) % (point - triangle.v3));
	double subArea3 = abs((point - triangle.v3) % (point - triangle.v1));

	//if the areas are equal, then th point is inside the triangle
	if ((subArea1 + subArea2 + subArea3 - triangleArea) <= Epsilon) // approximation may cause areas not to coincide perfectly
		return true;
	return false;
	*/
	return isInside3dTriangle(triangle2dTo3d(triangle), vector2dTo3d(point));
}

bool Simulation::isInside3dTriangle(Triangle3 triangle, Vector3 point) const
{
	double triangleArea = abs((triangle.v2 - triangle.v1) % (triangle.v3 - triangle.v1));
	double subArea1 = abs((point - triangle.v1) % (point - triangle.v2));
	double subArea2 = abs((point - triangle.v2) % (point - triangle.v3));
	double subArea3 = abs((point - triangle.v3) % (point - triangle.v1));

	//if the areas are equal, then th point is inside the triangle
	if (abs(subArea1 + subArea2 + subArea3 - triangleArea) <= Epsilon * triangleArea) // approximation may cause areas not to coincide perfectly
		return true;
	return false;
}

// takes a 2d triangle and retun the corrisponding 3d triangle (with z = 0)
Triangle3 Simulation::triangle2dTo3d(Triangle2 triangle) const
{
	Triangle3 _triangle = { vector2dTo3d(triangle.v1), vector2dTo3d(triangle.v2), vector2dTo3d(triangle.v3) };
	return _triangle;// TODO temporary variable not necessary
}

Vector3 Simulation::vector2dTo3d(Vector2 vector) const
{
	return {vector.X, vector.Y, 0.0};
}

// TODO
// simulate over time interval dt
bool Simulation::simulateInTimeInterval(double dt)// TODO
{
	
	// create a struct that identifies the first event and its informations (simEvent)
	SimEvent firstSimEvent;
	// call the funcion to find this event
	firstSimEvent = findFirstEvent(dt);
	// if there is an event, then perform first event and apply space and (in case) speed variations
	performEvent(firstSimEvent);
	// then if the time of the first event < dt,
	if ((firstSimEvent.relTime > 0.0) && (firstSimEvent.relTime < dt))
		// then call this function in the new interval = dt - first event time
		simulateInTimeInterval(dt - firstSimEvent.relTime);	
	
	return false;
}

// TODO
SimEvent Simulation::findFirstEvent(double dt) const
{
	SimEvent firstSimEvent;
	firstSimEvent.relTime = dt * 2.0;
	// TODO
	SimEvent newSimEvent = findGasMeshEvent(dt);
	if ( (newSimEvent.relTime > 0.0) && ( newSimEvent.relTime < firstSimEvent.relTime ) )
		firstSimEvent = newSimEvent;
	/*newSimEvent = findGasGasEvent(dt);
	if ((newSimEvent.relTime > 0.0) && (newSimEvent.relTime < firstSimEvent.relTime))
		firstSimEvent = newSimEvent;*/
	

	if ((firstSimEvent.relTime > 0.0) && (firstSimEvent.relTime < dt))
		return firstSimEvent;
	return SimEvent();
}

// return info on the first collision gas-container, if it occurs
SimEvent Simulation::findGasMeshEvent(double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	for (int i = 0; i < gasParticles.size(); i++)
	{
		newSimEvent = findFirstContainerCollision(i, dt);
		if ( (newSimEvent.relTime > 0.0) && ( newSimEvent.relTime < firstSimEvent.relTime ) )
			firstSimEvent = newSimEvent;
	}
	if ((firstSimEvent.relTime > 0.0) && (firstSimEvent.relTime < dt))
		return firstSimEvent;
	return SimEvent();
}

//
SimEvent Simulation::findFirstContainerCollision(int gasPointIndex, double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = simSettings.dt * 2.0;
	for (int i = 0; i < volumes.size(); i++)
	{
		// if ti is an interesting volume only
		if (volumes[i].isActive)
		{
			newSimEvent = findFirstTriangleCollision(gasPointIndex, i, dt);
			if ( (newSimEvent.relTime > 0.0) && (newSimEvent.relTime < firstSimEvent.relTime))
				firstSimEvent = newSimEvent;
		}
	}
	if ((firstSimEvent.relTime > 0.0) && (firstSimEvent.relTime < dt))
		return firstSimEvent;
	return SimEvent();
}

// return info ont the event representing the collison between the given particle and given volume, if it occurs in the step.
SimEvent Simulation::findFirstTriangleCollision(int gasPointIndex, int volumeIndex, double dt) const
{
	SimEvent firstSimEvent, newSimEvent;
	firstSimEvent.relTime = dt * 2.0;
	for (int i = 0; i < volumes[volumeIndex].mesh.plane.size(); i++) {
		newSimEvent = findTriangleCollision(gasPointIndex, volumeIndex, i, dt);
		if ((newSimEvent.relTime > 0.0) && (newSimEvent.relTime < firstSimEvent.relTime))
			firstSimEvent = newSimEvent;
	}
	if ((firstSimEvent.relTime > 0.0) && (firstSimEvent.relTime < dt))
		return firstSimEvent;
	return SimEvent();
}

// TODO check
// return info ont the event representing the collison between the given particle and given triangle, if it occurs in the step.
SimEvent Simulation::findTriangleCollision(int gasPointIndex, int volumeIndex, int triangleIndex, double dt) const
{
	GasParticle particle = gasParticles[gasPointIndex];
	GasParticle particleF = particle;// particle at the final position, excluding collisions
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
		// if it appens before dt, then continue! (or it would be useless to calculate the remaining collision parameters (so return the defalult SimEvent()))
		if ( ( collisionRelativeTime > dt ) || ( collisionRelativeTime <= 0.0 ) )// TODO < o <=???// < because if the step end exactly when there is a collision wi would copute it in the nex step
			return SimEvent();
		// ...and where
		Vector3 collisionPosition = particle.position + particle.speed * collisionRelativeTime;

		Triangle3 triangle;
		triangle =
		{
			{
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 0]].Position.X,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 0]].Position.Y,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 0]].Position.Z
			},
			{
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 1]].Position.X,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 1]].Position.Y,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 1]].Position.Z
			},
			{
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 2]].Position.X,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 2]].Position.Y,
				volumes[volumeIndex].mesh.Vertices[volumes[volumeIndex].mesh.Indices[3 * triangleIndex + 2]].Position.Z
			}
		};
		// then check if it is inside the triangle
		if (isInside3dTriangle(triangle, collisionPosition))
		{
			SimEvent collision;
			collision.relTime = collisionRelativeTime;
			collision.eventType = EventType::gasMeshCollision;
			collision.index1 = gasPointIndex;// index of the particle
			collision.index2 = volumeIndex;// index of the volume
			collision.index3 = triangleIndex;// index of the triangle
			collision.position1 = collision.position2 = collisionPosition;

			return collision;
		}
		// then calculate collion parameters // TODO NOT NECESSARY!
	}
	// there is no collision
	return SimEvent();
}

bool Simulation::performEvent(SimEvent simEvent)
{
	EventType eventType = simEvent.eventType;
	switch (eventType)
	{
	case EventType::gasMeshCollision:
		//Beep(2000, 100);
		preformGasMeshCollision(simEvent);
		break;// TODO
	case EventType::gasGasCollision:
		break;// TODO
	case EventType::gasSbCollision:
		break;// TODO
	case EventType::sbSbCollision:
		break;// TODO
	case EventType::none:
		performNullEvent(simEvent);
		break;// TODO
	default:
		break;
	}
	//std::cout << "evento.........\n";
	return false;// TODO return true
}

bool Simulation::preformGasMeshCollision(SimEvent simEvent)
{
	double dt = simEvent.relTime;
	GasParticle particle = gasParticles[simEvent.index1];
	int volumeIndex = simEvent.index2;
	int triangleIndex = simEvent.index3;
	Plane plane = volumes[volumeIndex].mesh.plane[triangleIndex];

	gasPositionIncrement(simEvent.relTime, simEvent.index1);// increment position except for interested particle
	// perform the collision (hard code);
	// apply movement
	gasParticles[simEvent.index1].position = simEvent.position1;
	// calculate new speed
	Vector3 planeNormal = { plane.a, plane.b, plane.c };
	double speedN = particle.speed * planeNormal;// speed component on the triangle normal
	double externPotential = potential(simEvent.position1 + particle.speed * (dt * Epsilon));// potential on the other side
	double dPotential = externPotential - particle.potential;
	Vector3 vSpeedN = planeNormal * speedN;
	double cEnergyN = (0.5 * speedN * speedN);
	Vector3 otherComponent = particle.speed - vSpeedN;
	// if the particle does not have enought energy
	if (cEnergyN < dPotential)
	{
		// then reflection
		vSpeedN = -vSpeedN;// reflection!
		particle.speed = vSpeedN + otherComponent;
	}
	else
	{
		vSpeedN = vSpeedN / abs(vSpeedN) * sqrt(2.0 * (cEnergyN - dPotential));// particle loses energy if dPotential > 0
		particle.speed = vSpeedN + otherComponent;
	}
	particle.position = simEvent.position1 + particle.speed * (dt * Epsilon);
	gasParticles[simEvent.index1] = particle;

	// TODO SB position and speed increment
	return false;// TODO return true
}

bool Simulation::gasPositionIncrement(double dt, int exI = -1)
{
	for (int i = 0; i < exI; i++)
		gasParticles[i].position += gasParticles[i].speed * dt;
	for (int i = exI + 1; i < gasParticles.size(); i++)
		gasParticles[i].position += gasParticles[i].speed * dt;
	//std::cout << dt << std::endl;
	//std::cout << "ciaoooooooooooooooooooo" << std::endl;
	return true;
}

bool Simulation::performNullEvent(SimEvent simEvent)
{
	gasPositionIncrement(getDt());// increment position except for interested particle
	// TODO SB position and speed increment
	return false;
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
	Vector3 v1 = { 0.0,0.0,0.0 }, v2 = { 4.0,0.0,4.0 }, v3 = { 0.0,4.0,4.0 }, v = { 1.0,1.0,2.0 };
	std::cout << abs((v - v1) % (v - v2)) << std::endl;
	//std::cout << abs((v - v1) % (v - v3)) << std::endl;
	//std::cout << abs((v - v2) % (v - v3)) << std::endl;
	//std::cout << abs((v - v1) % (v - v2)) + abs((v - v2) % (v - v3)) + abs((v - v1) % (v - v3)) << std::endl;
	//std::cout << abs((v2 - v1) % (v3 - v1)) << std::endl;
	
	std::cout << isInside3dTriangle({ v1,v2,v3 }, v) << std::endl;
	
}

void Simulation::printxyz(void) const
{
	std::cout << gasParticles[0].position.X << ", "
		<< gasParticles[0].position.Y << ", "
		<< gasParticles[0].position.Z << std::endl;
}
