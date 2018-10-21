//#include <QCoreApplication>
//#include "prova.h"

#include <iostream>
#include <math.h>
#include <Windows.h>

#include "OBJ_Loader.h"
#include "simulation.h"
#include "mathFunctions.h"

//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------

//takes the file path and passes it to the loader
bool loadObjFile(objl::Loader*, std::string);
//takes	meshes from the loader and passes them to the simulation to create volumes
bool loadVolumes(objl::Loader, Simulation*, std::vector<double>);

//-------------------------------------------------------------------------------------------------

int main(/*int argc, char *argv[]*/)
{
    /*QCoreApplication a(argc, argv);

    return a.exec();*/

	objl::Loader loader;
	Simulation simulation;
	std::vector<double> potentials = { -100, -100, -100 };


	simulation.ciao1();
	std::cout << "ciao da questa build\t\n";
	if (!loadObjFile(&loader, "box_stack.obj"/*"a.obj"*/))// bool loadout =
	{
		std::cout << "OBJ loading FAIL!!!" << std::endl;
		return 0;
	}
	//std::cout << loadout << std::endl;
	//std::cout << loader.LoadedMeshes.size() << std::endl;
	//std::cout << loader.LoadedMeshes[0].MeshName << std::endl;
	loadVolumes(loader, &simulation, potentials);
	std::cout << "Number of Volumes: " << simulation.getVolumesNumber() << std::endl;
	std::cout << std::endl;
	//simulation.printVertex(0);// useful
	//simulation.printPlanes(0);// useful
	std::cout << "\n\n\n\n\n";
	simulation.setGasParticleNumber(100);
	simulation.setInitialTemperature(1);
	simulation.createGasParticles();
	simulation.test();
	simulation.printPoints();// useful
	simulation.setDt(0.003);
	for (int i = 0; i < 50000; i++) {
		simulation.performOneStep();
		//simulation.printPoints();
		//std::cout << simulation.getDt() << std::endl;
		//Sleep(50);
		simulation.printxyz(0);
		//Beep(2000, 50);
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------

bool loadObjFile(objl::Loader* loader, std::string fileName) {
	// load obj file, return true if loaded
	if (!loader->LoadFile(fileName))
		return false;
	/*
	//dsplay mashes names
	for (int i = 0; i < Loader.LoadedMeshes.size(); i++)
	{
		objl::Mesh curMesh = Loader.LoadedMeshes[i];
		std::cout << "mesh " << i << ", name:\t" << curMesh.MeshName << std::endl;
	}
	*/
	return true;
	//return true;
}

bool loadVolumes(objl::Loader loader, Simulation* simulation, std::vector<double> potentials)
{
	// TODO check for the validity of meshes and potentials

	//std::cout << "ciaoooooooooooooooooooooooo " << loader.LoadedMeshes.size()  << std::endl;//TODO delete

	for (int i = 0; i < fmin(loader.LoadedMeshes.size(), potentials.size()); i++)
	{
		Vector3 Position;
		Vector3 Normal;
		objl::Mesh currMesh;
		Mesh currVolumeMesh;
		Volume currVolume;

		//TODO to implement// ma cosa???????????? forse è già fatto
		
		// load a mesh
		currMesh = loader.LoadedMeshes[i];
		
		// create currMesh
		currVolumeMesh.MeshName = currMesh.MeshName;
		for (int j = 0; j < currMesh.Vertices.size(); j++)
		{
			/*
			currVolumeMesh.Vertices[j].Position.X = currMesh.Vertices[j].Position.X;
			currVolumeMesh.Vertices[j].Position.Y = currMesh.Vertices[j].Position.Y;
			currVolumeMesh.Vertices[j].Position.Z = currMesh.Vertices[j].Position.Z;

			currVolumeMesh.Vertices[j].Normal.X = currMesh.Vertices[j].Normal.X;
			currVolumeMesh.Vertices[j].Normal.Y = currMesh.Vertices[j].Normal.Y;
			currVolumeMesh.Vertices[j].Normal.Z = currMesh.Vertices[j].Normal.Z;
			*/
			Position.X = currMesh.Vertices[j].Position.X;
			Position.Y = currMesh.Vertices[j].Position.Y;
			Position.Z = currMesh.Vertices[j].Position.Z;

			Normal.X = currMesh.Vertices[j].Normal.X;
			Normal.Y = currMesh.Vertices[j].Normal.Y;
			Normal.Z = currMesh.Vertices[j].Normal.Z;

			currVolumeMesh.Vertices.push_back( { Position, Normal } );
		}
		currVolumeMesh.Indices = currMesh.Indices;
		
		// create currVolume
		//currVolume = { potentials[i], currVolumeMesh };// does not work with newest Volume
		currVolume.deltaPotential = potentials[i];
		currVolume.mesh = currVolumeMesh;
		
		// add currVolume to simulation volumes
		simulation->addVolume(currVolume);
	}
	return true;
}