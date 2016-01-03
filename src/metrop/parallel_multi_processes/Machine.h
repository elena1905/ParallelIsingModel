/*=====================================================
 * Machine.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Machine to
 * assign machine geometry for local grids
 *=====================================================*/


#ifndef MACHINE_H_
#define MACHINE_H_


#include <iostream>
#include <stdarg.h>
#include <sstream>
#include <cstring>
#include "mpi.h"


// neighbour types
#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3


namespace wenchong
{

class Machine
{
public:
	Machine(int argc, char* argv[]);
	~Machine();

	int nProc;         // # of total processes
	int Rank;          // rank of process

	int nx;            // # of processes on x-axis
	int ny;            // # of processes on y-axis

	int x;             // x-coor of Machine geometry
	int y;             // y-coor of Machine geometry

	int Neighbour[4];  // neighbour types

	int Measures;      // # of measures
	int nSweeps;       // # of sweeps between two measures
	int nThrow;        // # of sweeps for thermalization

	char** Argv;       // argument values
};

};


#endif