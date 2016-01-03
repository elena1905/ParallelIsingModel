/*=====================================================
 * Metrop.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the class Metrop
 * that implements the Metropolis algorithm
 *=====================================================*/


#include "Metrop.h"


namespace wenchong
{

/*
 * construct the 2D lattice system
 */
Metrop::Metrop(Machine* host, int init, double beta, unsigned int seed)
	 : BaseLattice(host, init, beta, seed)
{
	// pre-compute the exponetial factors for
	// different combinations of neighbour values
	// as defined below:
	// delta = -delta(E) / (K(B)T)
	// delta = -2 * sum(si * sj) * (J / K(B)T)
	// delta = -2 * Beta * [si * sum(sj)]
	// delta = Factor * COMBS#
	// e^delta = exp(delta)

	Factor = -2 * Beta;

	ExpoDelta[0] = exp(Factor * (double)COMBS0);
	ExpoDelta[1] = exp(Factor * (double)COMBS1);
	ExpoDelta[2] = exp(Factor * (double)COMBS2);
	ExpoDelta[3] = exp(Factor * (double)COMBS3);
	ExpoDelta[4] = exp(Factor * (double)COMBS4);
}


/*
 * update the local lattice with even-odd ordering
 */
void Metrop::update(int numSweeps)
{
	for (int i = 0; i < numSweeps; i ++)
	{
		// update even sites and exchange boundary data
		updateLattice(EVEN);
		Spins->packBuffer(EVEN);
		Comms->sendBoundaryData(Spins);

		// update odd sites and exchange boundary data
		updateLattice(ODD);
		Spins->packBuffer(ODD);
		Comms->sendBoundaryData(Spins);
	}
}


/*
 * a half sweep of even sites or odd sites,
 * a half sweep of each is a whole sweep
 */
void Metrop::updateLattice(int evenOddFlag)
{
	for (int i = 0; i < nRow; i++)
	{
		// decide the starting site to update
		int start = (i + evenOddFlag) % 2;

		for (int j = start; j < nCol; j += 2)
		{
			// accept-reject process:
			// if proposal is accepted, flip spin
			if (isAccept(i, j))
				flipSpin(i, j);
		}
	}
}


/*
 * metropolis accept-reject procces,
 * check if proposal is accepted
 */
bool Metrop::isAccept(int row, int col)
{
	int current = (*Spins)(row, col);
	int left    = (*Spins)(row, col - 1);
	int right   = (*Spins)(row, col + 1);
	int up      = (*Spins)(row - 1, col);
	int down    = (*Spins)(row + 1, col);

	// the exponential factors are pre-computed
	// by the rules below:
	// delta = -delta(E) / (K(B)T)
	// delta = -2 * sum(si * sj) * (J / K(B)T)
	// delta = -2 * Beta * [si * sum(sj)]
	// delta = Factor * COMBS#

	// e^delta = exp(delta), which is computed
	// at construction as ExpoDelta,
	// so compute the combination value of [si * sum(sj)]
	// below, and check which case the ExpoDelta should
	// be applied

	int comb = current * (left + right + up + down);
	double expo = 0.0;

	switch (comb)
	{
		case COMBS0:
			expo = ExpoDelta[0];
			break;
		case COMBS1:
			expo = ExpoDelta[1];
			break;
		case COMBS2:
			expo = ExpoDelta[2];
			break;
		case COMBS3:
			expo = ExpoDelta[3];
			break;
		case COMBS4:
			expo = ExpoDelta[4];
			break;
		default:
			break;
	}

	/*===== condition 1: exp >= 1 --> accept =====*/
	if (expo >= 1.0)
	{
		return true;
	}

	/*========== condition 2: exp < 1 ==========*/
	// - if u ~ U(0, 1) < exp --> accept
	// - otherwise --> don't accept
	//double u = (double)rand() / RAND_MAX;
	double u = UniformDist(Generator);

	if ((u - expo) < 0.0)
	{
		return true;
	}

	return false;
}


/*
 * compute susceptibility X=<m^2>
 */
void Metrop::computeXt(void)
{
	int max = Measures * nSweeps;  // number of total sweeps
	double average = 0.0; // average spin in a sweep at time(meas) t
	std::ofstream ofsXt;

	// only write output to file inside one process
	if (Spins->Host->Rank == ROOT)
	{
		char* filenameXt = (char*)"xt.dat";
		ofsXt.open(filenameXt, std::ofstream::out);

		if (!ofsXt.is_open())
			throw std::invalid_argument("Metrop::(): Error opening file xt.dat!");
	}
	
	// number of total sites of the global lattice
	double nGlobalSites = (double)(Spins->nxGlobal * Spins->nyGlobal);

	// comuting X=<m^2>
	for (int i = 0; i < max; i++)
	{
		update(1);
		
		// sum up the spins at time t at every nSweeps
		if ((i % nSweeps) == 0)
		{
			double localSum = (double)sumSpins();
			double globalSum = 0.0;

			// get global sum over spins
			Comms->computeGlobalSum(&localSum, &globalSum);

			// computing magnetization
			average = globalSum / nGlobalSites;

			// computing X = magnetization ^ 2
			average *= average;

			// store X for Rho and Tau evaluation
			Xt.push_back(average);

			// calculate Mean and Variance of X
			Mean += average;
			Var += average * average;
			
			// store X in file
			if (Spins->Host->Rank == ROOT)
				ofsXt << average << std::endl;	
		}
	}
	
	if (Spins->Host->Rank == ROOT)
		ofsXt.close();
	
	Mean /= (double)Measures;
	Var = Var / (double)Measures - Mean * Mean;
}

};
