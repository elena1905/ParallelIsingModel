/*=====================================================
 * main.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This file is for executing the Metropolis algorithm
 * with multi processes
 *=====================================================*/


#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "mpi.h"

#include "Metrop.h"
#include "Machine.h"


using namespace wenchong;
using namespace std;


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
	Machine* host = new Machine(argc, argv);

	
	//============ seeding the RNG in each process ============//
	unsigned int seed = 0;              // the seed to generate u~(0, 1)
	unsigned int seedOfSeed = 25938026; // the seed to generate random seed above
	mt19937 gen(seedOfSeed);
	uniform_int_distribution<int> randSeed(1000000, 99999999);

	// each process gets a fixed random seed,
	// then the seed will be passed into the Metrop class
	// to construct RNG and RN distribution type
	for (int i = 0; i < host->Rank + 1; i++)
	{
		seed = randSeed(gen);
	}

	
	//============ initialization ============//
	int init = 1;       // initial spin value
	double beta = log(1 + sqrt(2)) / 2; // beta = J/K(B)T
	
	// class Metrop encapsulates the Metropolis lattice system
	Metrop c = Metrop(host, init, beta, seed);
	
	
	//===== start counting programme execution wall time =====//
	double elapsedTime = MPI_Wtime();

	
	//============ thermalization ============//
	c.update(host->nThrow);
	
	
	//============ compute X(t) ============//
	c.computeXt();
	

	//===== compute programme execution wall time =====//
	elapsedTime = MPI_Wtime() - elapsedTime;
	cout << fixed << setprecision(6) << "P" << host->Rank
			 << ": Elapsed wall time: " << elapsedTime << " seconds\n\n";
	

	//===== get X(t) from file to compute Rho(t) and Tau(t) =====//
	//char* filenameXt = (char*)"xt.dat";
	//c.retrieveXt(filenameXt);
	

	//===== compute autocorrelation Rho(t) and Tau(t) =====//
	c.computeRhoTau();	


	delete host;
	
	MPI_Finalize();
	

	return 0;
}


/* =============== End of Programme =============== */