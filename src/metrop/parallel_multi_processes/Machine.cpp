/*======================================================
 * Machine.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the class Machine to
 * assign machine geometry for local grids
 *======================================================*/


#include "Machine.h"


namespace wenchong
{

/*
 * constructor:
 * init Machine with command line arguments
 * and assign neighbours
 */
Machine::Machine(int argc, char* argv[])
{
	if (argc != 8)
	{
		std::cout << "Usage: ./exe L L np_x np_y nMeas nSweeps nTherms\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	// matrix should be of size LxL
	if (strcmp(argv[1], argv[2]) != 0)
	{
		std::cout << "Matrix width and length mismatch\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	Argv = argv;

	// get # of processes from command line
	nx = atoi(argv[3]);
	ny = atoi(argv[4]);

	// for the local lattice classes
	Measures = atoi(argv[5]);
	nSweeps = atoi(argv[6]);
	nThrow = atoi(argv[7]);

	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	// compute x, y coordinates in Machine geometry
	x = Rank % nx;
	y = Rank / nx;

	// assign Neighbour processor id with periodic boundary
	Neighbour[NORTH] = x + (y + 1) % ny * nx;
	Neighbour[SOUTH] = x + (y - 1 + ny) % ny * nx;
	Neighbour[EAST]  = (x + 1) % nx + y * nx;
	Neighbour[WEST]  = (x - 1 + nx) % nx + y * nx;
}


/*
 * Destructor
 */
Machine::~Machine()
{}

};


/*================ End of File ================*/