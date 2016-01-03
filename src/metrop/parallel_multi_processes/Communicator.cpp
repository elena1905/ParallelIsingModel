/*========================================================
 * Communicator.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the class Communicator
 * that manages parallel data communications
 *========================================================*/


#include "Communicator.h"


namespace wenchong
{

/*
 * constructor:
 * allocate resources to BoundaryComm structure
 */
Communicator::Communicator()
{
	nSend = 0;
	nRecv = 0;
}


/*
 * destructor:
 * free resources allocated to BoundaryComm structure
 */
Communicator::~Communicator()
{
	
}


/*
 * send and recv new boundary data
 */
void Communicator::sendBoundaryData(Field* f)
{
	sendToEast(f);
	sendToWest(f);
	sendToNorth(f);
	sendToSouth(f);

	// wait for local jobs to complete
	MPI_Waitall(nSend, SendRequest, SendStatus);
	MPI_Waitall(nRecv, RecvRequest, RecvStatus);

	MPI_Barrier(MPI_COMM_WORLD);

	nSend = 0;
	nRecv = 0;
}


/*
 * sum all the spin values in the lattice system globally
 */
void Communicator::computeGlobalSum(double* localSum, double* globalSum)
{
	MPI_Allreduce(localSum, globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
}


/*
 * send data to east
 * recv data from west
 */
void Communicator::sendToEast(Field* f)
{
	// send east boundary to east
	MPI_Isend(f->SendBuffer[EAST], f->nyBuffer, MPI_INT,
			  f->Host->Neighbour[EAST], 1000, MPI_COMM_WORLD,
			  SendRequest + nSend);
	nSend++;

	// recv west boundary from west
	MPI_Irecv(f->RecvBuffer[WEST], f->nyBuffer, MPI_INT,
			  f->Host->Neighbour[WEST], 1000, MPI_COMM_WORLD,
			  RecvRequest + nRecv);
	nRecv++;
}


/*
 * send data to west
 * recv data from east
 */
void Communicator::sendToWest(Field* f)
{
	// send west boundary data to west
	MPI_Isend(f->SendBuffer[WEST], f->nyBuffer, MPI_INT,
			  f->Host->Neighbour[WEST], 1001, MPI_COMM_WORLD,
			  SendRequest + nSend);
	nSend++;

	// recv east boundary from east
	MPI_Irecv(f->RecvBuffer[EAST], f->nyBuffer, MPI_INT,
			  f->Host->Neighbour[EAST], 1001, MPI_COMM_WORLD,
			  RecvRequest + nRecv);
	nRecv++;
}


/*
 * send data to north
 * recv data from south
 */
void Communicator::sendToNorth(Field* f)
{
	// send north boundary data to north
	MPI_Isend(f->SendBuffer[NORTH], f->nxBuffer, MPI_INT,
			  f->Host->Neighbour[NORTH], 1002, MPI_COMM_WORLD,
			  SendRequest + nSend);
	nSend++;

	// recv south boundary from south
	MPI_Irecv(f->RecvBuffer[SOUTH], f->nxBuffer, MPI_INT,
			  f->Host->Neighbour[SOUTH], 1002, MPI_COMM_WORLD,
			  RecvRequest + nRecv);
	nRecv++;
}


/*
 * send data to south
 * recv data from north
 */
void Communicator::sendToSouth(Field* f)
{
	// send south boundary data to south
	MPI_Isend(f->SendBuffer[SOUTH], f->nxBuffer, MPI_INT,
			  f->Host->Neighbour[SOUTH], 1003, MPI_COMM_WORLD,
			  SendRequest + nSend);
	nSend++;

	// recv north boundary data from norths
	MPI_Irecv(f->RecvBuffer[NORTH], f->nxBuffer, MPI_INT,
			  f->Host->Neighbour[NORTH], 1003, MPI_COMM_WORLD,
			  RecvRequest + nRecv);
	nRecv++;
}

};


/*================ End of File ================*/