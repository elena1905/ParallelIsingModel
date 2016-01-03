/*=====================================================
 * BaseLattice.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class BaseLattice
 * as an abstract base class for various algorithms
 * of the Ising Model
 *=====================================================*/


#ifndef BASELATTICE_H_
#define BASELATTICE_H_


#include <math.h>
#include <random>
#include <fstream>
#include <vector>
#include <stdexcept>
#include "mpi.h"

#include "Field.h"
#include "Communicator.h"


// the max delta time for evaluating the
// Rho(t) and Tau(t)
const int DELTA_TIME = 60;


namespace wenchong
{

class BaseLattice
{
public:
	BaseLattice(Machine* host, int init, double beta, unsigned int seed);
	virtual ~BaseLattice();

	virtual void update(int numSweeps) = 0;
	virtual void computeXt(void) = 0;
	void retrieveXt(const char* filenameXt);
	void computeRhoTau(void);
	void printLattice(void);

protected:
	int Size;       // size of the 2D lattice matrix
	int nRow;       // # of rows of the lattice matrix
	int nCol;       // # of cols of the lattice matrix

	int Measures;   // # of measures
	int nSweeps;    // # of sweeps between 2 measures
	double Beta;    // Beta = J/K(B)T
	double Factor;  // different for metrop & worm
	double Mean;    // Mean of the susceptibility
	double Var;     // variance of the susceptibility

	std::vector<double> Xt;  // to store susceptibility
	std::mt19937 Generator;  // random number generator
	std::uniform_real_distribution<double> UniformDist;

	Field* Spins;        // the 2D lattice spin matrix
	Communicator* Comms; // parallel data communication

	virtual void updateLattice(int evenOddFlag) = 0;
	void flipSpin(int row, int col);
	virtual bool isAccept(int row, int col) = 0; // accept-reject process
	int sumSpins(void);
};

};


#endif