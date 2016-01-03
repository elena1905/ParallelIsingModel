#===================================================
# Makefile
#
# @Author: Wenchong Chen
#
# Makefile for Metropolis algorithm of MSc Project
#===================================================


# variables
OBJS = main.o Machine.o Field.o Communicator.o BaseLattice.o Metrop.o
COMP = mpicxx -std=c++11 -O2 -lm -pg


# compile and link code
main: $(OBJS)
	$(COMP) -o main $(OBJS)

main.o: main.cpp Machine.h Field.h Communicator.h BaseLattice.h Metrop.h
	$(COMP) -c main.cpp

Machine.o: Machine.cpp Machine.h
	$(COMP) -c Machine.cpp

Field.o: Field.cpp Field.h Machine.h
	$(COMP) -c Field.cpp

Communicator.o: Communicator.cpp Communicator.h Field.h
	$(COMP) -c Communicator.cpp

BaseLattice.o: BaseLattice.cpp BaseLattice.h Field.h Communicator.h
	$(COMP) -c BaseLattice.cpp

Metrop.o: Metrop.cpp Metrop.h BaseLattice.h
	$(COMP) -c Metrop.cpp


# clean target
clean:
	rm -f main $(OBJS)
