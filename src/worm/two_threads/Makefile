#=====================================================
# Makefile
#
# @Author: Wenchong Chen
#
# Makefile for threaded Worm algorithm of MSc Project
#=====================================================


OBJS = main.o Machine.o Field.o Communicator.o BaseLattice.o Worm.o
COMP = mpicxx -std=c++11 -O2 -lm -pg -fopenmp


# compile and link code
main: $(OBJS)
	$(COMP) -o main $(OBJS)

main.o: main.cpp Machine.h Field.h Communicator.h BaseLattice.h Worm.h
	$(COMP) -c main.cpp

Machine.o: Machine.cpp Machine.h
	$(COMP) -c Machine.cpp

Field.o: Field.cpp Field.h Machine.h
	$(COMP) -c Field.cpp

Communicator.o: Communicator.cpp Communicator.h Field.h
	$(COMP) -c Communicator.cpp

BaseLattice.o: BaseLattice.cpp BaseLattice.h Field.h Communicator.h
	$(COMP) -c BaseLattice.cpp

Worm.o: Worm.cpp Worm.h BaseLattice.h
	$(COMP) -c Worm.cpp


# clean target
clean:
	rm -f main $(OBJS)
