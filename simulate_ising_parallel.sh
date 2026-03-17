#!/bin/bash

FOLDER_NAME=IsotropicIsingParallelTrajectories_D16 #IsotropicIsingTrajectories_D8
for tstar in `cat temp_sel.txt`
do
	echo T=$tstar
	mkdir -p $FOLDER_NAME/$tstar
	if test -f "$FOLDER_NAME/$tstar/ising_traj_simul.bitmat"
	then
		echo "Grid exists"
		continue
	fi
	g++ ParallelIsingSampler.cpp -g -w -o pis -O3 -I ./libraries_cpp -fopenmp -DTEMP_SIM=$tstar
	time ./pis
	mv ising_traj_simul.bitmat $FOLDER_NAME/$tstar/
done
