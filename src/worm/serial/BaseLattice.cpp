/*=====================================================
 * BaseLattice.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the non-pure
 * functions of the class BaseLattice
 *=====================================================*/


#include "BaseLattice.h"


namespace wenchong
{

/*
 * construct the 2D BaseLattice system
 */
BaseLattice::BaseLattice(Machine* host, int init, double beta, unsigned int seed) :
	Generator(seed),
	UniformDist(0.0, 1.0)
{
	if (init != 1 && init != -1)
		throw std::invalid_argument("BaseLattice::(): invalid initial spin");

	Comms = new Communicator();
	
	Spins = new Field(host);
	Spins->init(init);

	Size = Spins->nData;   // nRow * nCol
	nRow = Spins->nxLocal;
	nCol = Spins->nyLocal;

	Beta = beta;
	Measures = host->Measures;
	nSweeps = host->nSweeps;
	Mean = 0.0;
	Var = 0.0;
}


/*
 * destruct the BaseLattice system and release memory
 */
BaseLattice::~BaseLattice()
{
	delete Comms;
	delete Spins;
}


/*
 * flip the value of spin at given (row, col)
 */
void BaseLattice::flipSpin(int row, int col)
{
	if (row >= nRow || row < 0 || col >= nCol || col < 0)
		throw std::out_of_range("BaseLattice::(): row and col our of bounds");

	(*Spins)(row, col) *= -1;
}


/*
 * sum up all values of spin in the lattice
 */
int BaseLattice::sumSpins(void)
{
	return Spins->sumData();
}


/*
 * load X from file to Xt vector
 */
void BaseLattice::retrieveXt(const char* filenameXt)
{
	// only get data from file inside one process
	if (Spins->Host->Rank == ROOT)
	{
		int count = -1;
		double average = 0.0; // average spin in a sweep at time(meas) t
		std::ifstream ifsXt;

		ifsXt.open(filenameXt, std::ifstream::in);

		if (!ifsXt.is_open())
			throw std::invalid_argument("BaseLattice::(): Error opening file xt.dat!");
	
		// load X from file to Xt vector,
		// compute Mean and Var
		while (!ifsXt.eof())
		{
			count++;

			// computing X = magnetization ^ 2
			ifsXt >> average;

			// store X for Rho and Tau evaluation
			Xt.push_back(average);

			// calculate Mean and Variance of X
			Mean += average;
			Var += average * average;
		}
		
		// check if the number of data in the file
		// is equal to the number of measurements
		if (count != Measures)
		{
			throw std::invalid_argument("BaseLattice::(): data size not match");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	
		ifsXt.close();
	
		Mean /= (double)Measures;
		Var = Var / (double)Measures - Mean * Mean;
	}
}


/*
 * DO NOT call this method before computeXt() is called
 */
void BaseLattice::computeRhoTau(void)
{
	if (Spins->Host->Rank == ROOT)
	{
		if (Mean == 0.0 || Var == 0.0)
			throw std::invalid_argument("BaseLattice::(): uninitialised Mean and Var");

		double Rt = 0.0;  // autocovariance R(t)
		double Rho = 1.0; // autocorrelation Rho(t)
		double Tau = 0.5; // integrated autocorrelation time Tau(t)

		// filenames to store Rho and Tau
		char* filenameRho = (char*)"rho.dat";
		char* filenameTau = (char*)"tau.dat";

		// open files for writing
		std::ofstream ofsRho(filenameRho, std::ofstream::out);
		std::ofstream ofsTau(filenameTau, std::ofstream::out);

		if (!ofsRho.is_open() || !ofsTau.is_open())
			throw std::invalid_argument("BaseLattice::(): Error opening files!");
	
		// write the first values at t=0 to files
		ofsRho << "0 " << Rho << std::endl;
		ofsTau << "0 " << Tau << std::endl;
	
		// iterate through t from 1 to DELTA_TIME
		for (int t = 1; t <= DELTA_TIME; t++)
		{
			Rt = 0.0;
		
			// compute R(delta) = E[X(t)X(t+delta)] - u^2
			int i = 0;
			while ((i + t) < Measures)
			{
				Rt += Xt[i] * Xt[i + t];
				i += 1;
			}
			Rt = Rt / i - Mean * Mean;
		
			// Rho(t) = R(t) / (sigma^2)
			Rho = Rt / Var;
		
			// Tau(t) = 1/2 + sum(Rho(t))
			Tau += Rho;
		
			// write data to files
			ofsRho << t << " " << Rho << std::endl;
			ofsTau << t << " " << Tau << std::endl;
		}
	
		ofsRho.close();
		ofsTau.close();
	}
}


/*
 * print the spins of the BaseLattice matrix
 * for testing purpose
 */
void BaseLattice::printLattice(void)
{
	for (int i = 0; i < nRow; i++)
	{
		for (int j = 0; j < nCol; j++)
		{
			if ((*Spins)(i, j) < 0)
			{
				printf("%d ", (*Spins)(i, j));
			}
			else
			{
				printf(" %d ", (*Spins)(i, j));
			}

			// go to a new line when reaching end of row
			if ((j % Size) == (Size - 1))
			{
				printf("\n");
			}
		}
	}
}

};
