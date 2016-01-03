/*=====================================================
 * Field.cpp
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This definition file implements the class Field for
 * local lattice storage and accessing
 *=====================================================*/


#include "Field.h"


namespace wenchong
{

/*
 * alloc resources to Field and init Field boundaries
 */
Field::Field(Machine* m)
{
	Host = m;

	// get # of points on x, y axises
	nxGlobal = atoi(Host->Argv[1]);
	nyGlobal = atoi(Host->Argv[2]);

	// compute nxLocal and nyLocal
	checkGrid();

	// active Data size
	nData = nxLocal * nyLocal;

	Data = new int[nData];

	// compute buffer sizes
	nxBuffer = nxLocal / 2 + 1;
	nyBuffer = nyLocal / 2 + 1;

	SendBuffer[NORTH] = new int[nxBuffer];
	SendBuffer[SOUTH] = new int[nxBuffer];
	SendBuffer[EAST]  = new int[nyBuffer];
	SendBuffer[WEST]  = new int[nyBuffer];

	RecvBuffer[NORTH] = new int[nxBuffer];
	RecvBuffer[SOUTH] = new int[nxBuffer];
	RecvBuffer[EAST]  = new int[nyBuffer];
	RecvBuffer[WEST]  = new int[nyBuffer];
}


/*
 * init data and boundary data
 */
void Field::init(int initVal)
{
	// init lattice data
	for (int i = 0; i < nData; i++)
		Data[i] = initVal;

	// init boundary data/receive buffer
	for (int i = 0; i < 4; i++)
	{
		if (i == NORTH || i == SOUTH)
		{
			for (int j = 0; j < nxBuffer; j++)
				RecvBuffer[i][j] = initVal;
		}
		else
		{
			for (int j = 0; j < nyBuffer; j++)
				RecvBuffer[i][j] = initVal;
		}
	}
}


/*
 * release resources allocated to Field
 */
Field::~Field()
{
	delete [] Data;

	for (int i = 0; i < 4; i++)
	{
		delete [] SendBuffer[i];
		delete [] RecvBuffer[i];
	}
}


/*
 * get the value of spin at given (row, col),
 * may need boundary data from receive buffers
 */
int& Field::operator() (int row, int col)
{
	if (row > nxLocal || row < -1 || col > nyLocal || col < -1)
	{
		throw std::out_of_range("Field::(): row and col our of bounds");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// get data from North buffer
	if (col == nyLocal)
		return RecvBuffer[NORTH][row / 2];

	// get data from South buffer
	if (col == -1)
		return RecvBuffer[SOUTH][row / 2];

	// get data from East buffer
	if (row == nxLocal)
		return RecvBuffer[EAST][col / 2];

	// get data from West buffer
	if (row == -1)
		return RecvBuffer[WEST][col / 2];

	// get data from local Data array
	return Data[row * nyLocal + col];
}


/*
 * sum up all values of spin in the lattice
 */
int Field::sumData(void)
{
	int sum = 0;

	for (int i = 0; i < nxLocal; i++)
		for (int j = 0; j < nyLocal; j++)
			sum += Data[i * nyLocal + j];

	return sum;
}


/*
 * pack boundary data into the send buffers
 */
void Field::packBuffer(int evenOddFlag)
{
	int row, col, start;

	start = evenOddFlag % 2;

	// pack data to West buffer
	for (int y = start; y < nyLocal; y += 2)
	{
		SendBuffer[WEST][y / 2] = Data[y];
	}

	// pack data to South buffer
	for (int x = start; x < nxLocal; x += 2)
	{
		SendBuffer[SOUTH][x / 2] = Data[x * nyLocal];
	}

	// pack data to East buffer
	row = nxLocal - 1;
	start = (row + evenOddFlag) % 2;
	for (int y = start; y < nyLocal; y += 2)
	{
		SendBuffer[EAST][y / 2] = Data[row * nyLocal + y];
	}

	// pack data to North buffer
	col = nyLocal - 1;
	start = (col + evenOddFlag) % 2;
	for (int x = start; x < nxLocal; x += 2)
	{
		SendBuffer[NORTH][x / 2] = Data[x * nyLocal + col];
	}
}


/*
 * check Grid and Machine size compatability and
 * compute offsets of Grid relative to global
 */
void Field::checkGrid(void)
{	
	// check machine size compatability
	if (Host->nx * Host->ny != Host->nProc)
	{
		std::cout << "Incompatible grid parallelisation: "
				  << Host->nx << " x " << Host->ny
				  << " on " << Host->nProc << " processes\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// check local geometry compatability
	nxLocal = nxGlobal / Host->nx;
	if (nxLocal * Host->nx != nxGlobal)
	{
		std::cout << "Grid mismatch in X-dir\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	nyLocal = nyGlobal / Host->ny;
	if (nyLocal * Host->ny != nyGlobal)
	{
		std::cout << "%ZGrid mismatch in Y-dir\n";
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// compute offsets relative to global coor system
	xOffset = Host->x * nxLocal;
	yOffset = Host->y * nyLocal;
}

};


/*================ End of File ================*/