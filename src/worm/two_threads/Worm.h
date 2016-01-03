/*=====================================================
 * Worm.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Worm that
 * implements the Worm algorithm
 *=====================================================*/


#ifndef WORM_H_
#define WORM_H_


#include <iostream>
#include <omp.h>

#include "BaseLattice.h"


#define  X  0  // index for positive x direction
#define  Y  1  // index for positive y direction


#define  H  0
#define  T  1


namespace wenchong
{

class Worm : public BaseLattice
{
public:
	Worm(Machine* host, int init, double beta, unsigned int seed);
	~Worm();

	virtual void update(int numSweeps);
	virtual void computeXt(void);

private:
	int nLinks;    // size of matrix to store links
	int nKick;     // # of times that head and tail touch
	int Head[2];   // Head site, position fixed
	int Tail[2];   // Tail site that crawls
	int Neib[2][2];

	double TanhBeta; // tanh(beta) = e^(-u)

	// to store links between sites
	bool* Links;

	omp_lock_t Lock;

	void randKick(void);
	void randNeighbour(int* row, int* col);
	bool isLinkOn(const int* curSite, const int* newSite);
	bool isAdjacent(const int* curSite, const int* newSite);
	void turnOnLink(const int* curSite, const int* newSite);
	void turnOffLink(const int* curSite, const int* newSite);

	virtual void updateLattice(int evenOddFlag);
	virtual bool isAccept(int row, int col);
	void checkAccept(int* curSite, int* neibSite);
};

};


#endif