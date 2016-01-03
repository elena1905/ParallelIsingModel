/*=====================================================
 * Communicator.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Communicator
 * for parallel data communications
 *=====================================================*/


#ifndef COMMUNICATOR_H_
#define COMMUNICATOR_H_


#include "mpi.h"

#include "Field.h"


namespace wenchong
{

class Communicator
{
public:
	Communicator();
	~Communicator();

	void sendBoundaryData(Field* f);
	void computeGlobalSum(double* localSum, double* globalSum);

private:
	int nSend;
	int nRecv;

	MPI_Request SendRequest[4];
	MPI_Request RecvRequest[4];

	MPI_Status SendStatus[4];
	MPI_Status RecvStatus[4];

	// to exchage data with neighbours
	void sendToEast(Field* f);
	void sendToWest(Field* f);
	void sendToNorth(Field* f);
	void sendToSouth(Field* f);
};

};


#endif