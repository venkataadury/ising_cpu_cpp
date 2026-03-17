#include <iostream>
#include <fstream>
#include "IsingClusterFlip.hpp"
#include <omp.h>

#define GRID_DIM 16
#define EQUIL_STEPS 5000000 // 10M
#define SIM_STEPS 20000000 // 100M
#define WRITE_STEPS 5000 // 25K
#define X_COUPL 1.0
#define Y_COUPL 1.0
#define N_REPLICAS 12
//#define BETA_SIM 1.0

int main()
{
    const double BETA_SIM=1.0d/TEMP_SIM;
    std::vector<ClusterFlipIsingModel<GRID_DIM>*> grids;
    for(int i=0;i<N_REPLICAS;i++) grids.push_back(new ClusterFlipIsingModel<GRID_DIM>(/*beta=*/BETA_SIM,/*Jx=*/X_COUPL,/*Jy=*/Y_COUPL));
    
#pragma omp parallel for
    for(int I=0;I<N_REPLICAS;I++)
    {
        for(long i=0;i<EQUIL_STEPS;i++) grids[I]->step(); // Equilibriation
    }
    std::cout << "Equilibriation completed\n";
    
    int N_SAMPLES=(int)(SIM_STEPS/WRITE_STEPS)+1;
    std::vector<bool**> buffer(N_SAMPLES);
    
    //volatile int loc=0;
    //volatile long last_step=0;
#pragma omp parallel for
    for(int I=0;I<N_REPLICAS;I++)
    {
        for(long i=1;i<=SIM_STEPS;i++)
        {
            grids[I]->step();
            if(i%WRITE_STEPS==0)// {last_step=i; buffer[loc++]=grids[I]->cloneSpinMatrix();}
            {
                int loct=(i/WRITE_STEPS);
                if((loct-I)%N_REPLICAS==0) buffer[loct]=grids[I]->cloneSpinMatrix();
            }
        }
    }
    
    int loc=(SIM_STEPS/WRITE_STEPS);
    std::cout << "Simulation ends with "<<loc<<" frames\n";
    
    std::ofstream fout; fout.open("ising_traj_simul.bitmat");
    for(int I=1;I<=loc;I++)
    {
        //std::cout << "Frame "<<I<<"\n";
        fout << GRID_DIM << "\n";
        for(int i=0;i<GRID_DIM;i++)
        {
            for(int j=0;j<GRID_DIM;j++)
            {
                fout << ((int)buffer[I][i][j])<<" ";
            }
            fout << "\n";
        }
    }
    fout.close();
    
    for(int I=1;I<loc;I++)
    {
        for(int i=0;i<GRID_DIM;i++) delete[] buffer[I][i];
        delete[] buffer[I];
    }
    
    /*
    std::ofstream fout; fout.open("ising_traj_simul.bitmat");
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
    */
    return 0;
}
