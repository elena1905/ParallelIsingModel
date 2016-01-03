/*=====================================================
 * Worm.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Worm that
 * implements the serial Worm algorithm
 *=====================================================*/


#ifndef WORM_H_
#define WORM_H_


#include <iostream>
#include <omp.h>

#include "BaseLattice.h"


#define  X  0  // index for positive x direction
#define  Y  1  // index for positive y direction


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

	double TanhBeta; // tanh(beta) = e^(-u)

	// to store links between sites
	bool* Links;

	void randKick(void);
	void randNeighbour(const int* curSite, int* newSite);

	virtual void updateLattice(int evenOddFlag);
	virtual bool isAccept(int row, int col);
	bool isAccept(int linkID);
};

};


#endif