/*=====================================================
 * Field.h
 *
 * @Author: Wenchong Chen
 *
 * Code for MSc Project
 *
 * This header file declares the class Field for
 * local lattice storage and accessing
 *=====================================================*/


#ifndef FIELD_H_
#define FIELD_H_


#include <iostream>
#include <stdexcept>
#include "mpi.h"

#include "Machine.h"


#define EVEN 0  // for even ordering
#define ODD  1  // for odd ordering


#define ROOT 0  // root processor


namespace wenchong
{

class Field
{
public:
	Field(Machine* m);
	~Field();

	void init(int initVal); // init Data array
	int sumData(void);      // sum data in Data array
	int& operator() (int row, int col);
	void packBuffer(int evenOddFlag); // pack boundary to send

	int nxGlobal;       // # of points on global x-axis
	int nyGlobal;       // # of points on global y-axis

	int nData;          // # of active/local data
	int nxLocal;        // # of points on local x-axis
	int nyLocal;        // # of points on local y-axis

	int xOffset;        // offset relative to global x-axis
	int yOffset;        // offset relative to global y-axis
	
	int nxBuffer;       // buffer size of x axis
	int nyBuffer;       // buffer size of y axis

	int* Data;          // data in the Field
	int* SendBuffer[4]; // buffer of boundary data to send
	int* RecvBuffer[4]; // buffer of boundary data to receive

	Machine* Host;      // host processor

private:
	void checkGrid(void); // check Grid and Machine compatability
};

};


#endif