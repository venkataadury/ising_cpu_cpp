#!/bin/bash

FOLDER_NAME=AnisotropicIsingTrajectories_D8 #IsotropicIsingTrajectories_D8
for tstar in `cat temp_sel.txt`
do
	echo T=$tstar
	mkdir -p $FOLDER_NAME/$tstar
	if test -f "$FOLDER_NAME/$tstar/ising_traj.bitmat"
	then
		echo "Grid exists"
		continue
	fi
	g++ IsingSampler.cpp -g -o ising_sampler -O3 -I ~/cpp/libraries_cpp/ -w -DTEMP_SIM=$tstar
	time ./ising_sampler
	mv ising_traj.bitmat $FOLDER_NAME/$tstar/
done
