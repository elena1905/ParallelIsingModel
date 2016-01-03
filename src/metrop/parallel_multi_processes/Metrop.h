/*=====================================================
 * Metrop.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Metrop that
 * implements the Metropolis algorithm
 *=====================================================*/


#ifndef METROP_H_
#define METROP_H_


#include "BaseLattice.h"


// valuses of combinations of [si * sum(sj)]
//const int COMBS[5] = {-4, -2, 0, 2, 4};
#define COMBS0 -4
#define COMBS1 -2
#define COMBS2  0
#define COMBS3  2
#define COMBS4  4


namespace wenchong
{

class Metrop : public BaseLattice
{
public:
	Metrop(Machine* host, int init, double beta, unsigned int seed);
	~Metrop() {}
	
	virtual void update(int numSweeps);
	virtual void computeXt(void);

private:
	double ExpoDelta[5];  // to store pre-computed factors
	virtual void updateLattice(int evenOddFlag);
	virtual bool isAccept(int row, int col);
};

};


#endif