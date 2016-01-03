#==================================================================
# README.TXT
#
# @Author: Wenchong Chen
# Copyright of this project and the code is owned and reserved by
# the author Wenchong Chen.
#
# Instructions for MSc Project
#==================================================================
#
# 1. Makefile is provided. In the command line,
#    type in 'make' to compile and link the code,
#    the executable file 'main' will be
#    generated.
#
# 2. The metrop code is written for multi processes,
#    use the following command to run the
#    executable file:
#
#    $mpirun -n nProc ./main L L nx_p ny_p nMeas nSweeps nTherms
#
#    where L is problem size,
#    nx_p is the # of processes on x-axis,
#    ny_y is the # of processes on y-axis,
#    nMeas is # of measurements,
#    nSweeps is # of sweeps between two measurements,
#    nTherms is # of thermalisation.
#
#    An example for the execution command is:
#
#    $mpirun -n 4 ./main 16 16 2 2 10000 2 55
#
# 3. The worm code is written for one process in the
#    parallel framework. For both the 1 thread version
#    and the two threads version, use the following
#    example to run the executable:
#    
#    $mpirun -n 1 ./main 16 16 1 1 10000 1 1000
#
#=======================================================
