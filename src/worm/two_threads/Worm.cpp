/*=====================================================
 * Metrop.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the class Worm that
 * implements the threaed Worm algorithm
 *=====================================================*/


#include "Worm.h"


namespace wenchong
{

/*
 * construct Worm
 */
Worm::Worm(Machine* host, int init, double beta, unsigned int seed)
	: BaseLattice(host, init, beta, seed)
{
	// init links
	nLinks = Size * 2;
	Links = new bool[nLinks];
	Locks = new omp_lock_t[nLinks];
	for (int i = 0; i < nLinks; i++)
	{
		Links[i] = false;
		omp_init_lock(&Locks[i]);
	}

	nKick = 0;
	TanhBeta = tanh(Beta);
	Factor = pow(nRow, 1.75) / (double)Size; // L^(7/4)

	for (int i = 0; i < 2; i++)
		Head[i] = Tail[i] = 0;

	// Head and Tail start from a random site
	randKick();
}


/*
 * destructor: release resources
 */
Worm::~Worm()
{
	for (int i = 0; i < nLinks; i++)
		omp_destroy_lock(&Locks[i]);

	delete [] Locks;
	delete [] Links;
}


/*
 * update the worm system
 */
void Worm::update(int numSweeps)
{
	for (int i = 0; i < numSweeps; i++)
	{
		// it's serial code, doesn't matter
		// even or odd ordering
		updateLattice(EVEN);
	}
}


/*
 * the worm crawls here, equivalent to a sweep,
 * using two threads for the Head and the Tail,
 * the initial Head/Tail position has been specified
 * in the constructor
 */
void Worm::updateLattice(int evenOddFlag)
{
	int halfSize = Size / 2;
	int newSite[2][2] = {0, 0, 0, 0};
	int* curSites[2] = {Head, Tail};
	nKick = 0;

	// there are two threads, so the iteration size
	// should be half of the original size
	for (int i = 0; i < halfSize; i++)
	{
		// create two threads for Head and Tail
		#pragma omp parallel for num_threads(2)
		for (int j = 0; j < 2; j++)
		{
			// generate random neighbour
			randNeighbour(curSites[j], newSite[j]);

			// decide the begin and end of the link,
			// if curSite is on the left of newSite, or
			// curSite is at the bottom of newSite,
			// curSite is the begin site, newSite is the end site
			int* beginSite = curSites[j];  // begin site of the link

			// otherwise, newSite is the begin site, and
			// curSite is the end site
			if ((newSite[j][X] + nRow - 1) % nRow == curSites[j][X] ||
				(newSite[j][Y] + 1) % nCol == curSites[j][Y])
			{
				beginSite = newSite[j];
			}

			// compute the beginSite's position, and the link ID
			// for positive Y direction starting from beginSite
			int beginPos = (beginSite[X] * nCol + beginSite[Y]) * 2;
			int linkID = beginPos + Y;

			// check for positive X direction
			if ((newSite[j][Y] + nRow - 1) % nRow == curSites[j][Y] ||
				(newSite[j][Y] + 1) % nCol == curSites[j][Y])
			{
				linkID = beginPos + X;
			}

			// lock the link with given linkID for
			// accept-reject process, if linkIDs of Head and Tail
			// are different, this piece of code will run in parallel,
			// otherwise, it will be locked and executed in serial
			omp_set_lock(&Locks[linkID]);
			if (isAccept(linkID))
			{
				// if link is off, turn it on,
				// otherwise, turn it on
				if (!Links[linkID])
					Links[linkID] = true;
				else
					Links[linkID] = false;

				// move curSite to newSite
				curSites[j][X] = newSite[j][X];
				curSites[j][Y] = newSite[j][Y];
			}
			omp_unset_lock(&Locks[linkID]);
		}

		// if Head and Tail meet
		if (Tail[X] == Head[X] && Tail[Y] == Head[Y])
		{
			// since the number of iterations is half of
			// the original size, nKick should be
			// incremented by 2
			nKick += 2;
			
			// generate RV u~[0, 1],
			// if u >= 0.5, kick the Head and Tail
			// together to a random site
			double u = UniformDist(Generator);
			if (u >= 0.5)
				randKick();
		}
	}
}


/*
 * generate rv(x, y) in the range ([0, nRow-1], [0, nCol-1]),
 * set Head and Tail to this random site
 */
void Worm::randKick(void)
{
	int x = rand() % nRow;
	int y = rand() % nCol;

	while (x == Head[X] && y == Head[Y])
	{
		x = rand() % nRow;
		y = rand() % nCol;
	}

	Head[X] = Tail[X] = x;
	Head[Y] = Tail[Y] = y;
}


/*
 * select the nearest neighbour randomly
 */
void Worm::randNeighbour(const int* curSite, int* newSite)
{
	newSite[X] = curSite[X];
	newSite[Y] = curSite[Y];

	int count = 0;
	double u = UniformDist(Generator);

	// each neighbour is selected with
	// equal probability
	if (u >= 0 && u < 0.25)
		newSite[X] = (curSite[X] + nRow - 1) % nRow;
	else if (u >= 0.25 && u < 0.5)
		newSite[Y] = (curSite[Y] + 1) % nCol;
	else if (u >= 0.5 && u < 0.75)
		newSite[X] = (curSite[X] + 1) % nRow;
	else
		newSite[Y] = (curSite[Y] + nCol - 1) % nCol;
}


/*
 * accept-reject process:
 * this is specially used by the threaded worm
 * accept with Prob = min(1, e^u * (2 * kl - 1)).
 * 
 * :. tanh(beta) = e^(-u),
 * .: Prob = min(1, tanh(beta)^(1 - 2 * kl)).
 */
bool Worm::isAccept(int linkID)
{
	/*
	 * condition 1: kl = 1,
	 * tanh(beta)^(1-2kl) = 2.41421356 > 1,
	 * --> accept
	 */
	if (Links[linkID])
		return true;

	/*
	 * condition 2: kl = 0,
	 * tanh(beta)^(1-2kl) = 0.41421356 < 1,
	 * - if u ~ U(0, 1) <= tanh(beta)^(1-2kl),
	 *   --> accept
	 * - otherwise,
	 *   --> don't accept
	 */
	double u = UniformDist(Generator);
	if (u <= TanhBeta)
		return true;

	return false;
}


/*
 * DO NOT use this function
 */
bool Worm::isAccept(int row, int col)
{
	return true;
}


/*
 * L^(7/4) / X = (nKick/N) * L^(7/4), where N = L*L
 */
void Worm::computeXt(void)
{
	int max = Measures * nSweeps;  // number of total sweeps
	double average = 0.0; // (nKick/N) * L^(7/4), N = L*L
	
	// open file for writing L^(7/4) / X
	char* filenameXt = (char*)"xt.dat";
	std::ofstream ofsXt(filenameXt, std::ofstream::out);
	
	if (!ofsXt.is_open())
	{
		throw std::invalid_argument("Worm::(): Error opening files!");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	for (int i = 0; i < max; i++)
	{
		update(1);
		
		// calculate Xt at time t after nSweeps sweeps
		if (((i + 1) % nSweeps) == 0)
		{
			// average = (nKick/N) * L^(7/4), where N = L*L
			// Factor = L^(7/4) / N
			average = Factor * (double)nKick;

			Xt.push_back(average);

			Mean += average;
			Var += average * average;
			
			// store average magnetization in file
			ofsXt << average << std::endl;
		}
	}
	
	ofsXt.close();
	
	Mean /= (double)Measures;
	Var = Var / (double)Measures - Mean * Mean;

	std::cout << "Mean: " << Mean << "\n\n";
}

};
