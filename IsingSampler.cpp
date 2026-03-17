#include <iostream>
#include <fstream>
#include "IsingClusterFlip.hpp"

#define GRID_DIM 8
#define EQUIL_STEPS 10000000 // 10M
#define SIM_STEPS 100000000 // 100M
#define WRITE_STEPS 50000 // 25K
#define X_COUPL 1.0
#define Y_COUPL 2.0
//#define BETA_SIM 1.0

int main()
{
    ClusterFlipIsingModel<GRID_DIM> grid(/*beta=*/100.0);
    grid.writeToStream(true);
    std::cout << "\n\n";
    grid.constructBondingGrids();
    grid.bondingGridToStream();
    std::cout << "\n\n";
    grid.updateClusterLabels();
    grid.clusterLabelsToStream();
    std::cout << "\n";
    grid.step();
    grid.writeToStream(true);
    
    
    const double BETA_SIM=1.0d/TEMP_SIM;
    ClusterFlipIsingModel<GRID_DIM> final_grid(/*beta=*/BETA_SIM,/*Jx=*/X_COUPL,/*Jy=*/Y_COUPL);
    for(long i=0;i<EQUIL_STEPS;i++) final_grid.step(); // Equilibriation
    std::cout << "Equilibriation completed\n";
    
    std::ofstream fout; fout.open("ising_traj.bitmat");
    for(long i=0;i<SIM_STEPS;i++)
    {
        final_grid.step();
        if(i%WRITE_STEPS==0)
        {
            fout << GRID_DIM << "\n";
            final_grid.writeToStream(false,fout);
            //std::cout << i << " of " << SIM_STEPS << "\n";
        }
    }
    fout.close();
    return 0;
}
