#ifndef SOLVER
#define SOLVER


#include "./grid.h"

class RimanSolver1D
{
private:
    Grid grid;

    bool check_grid(const Grid &grid) ;

public:
    RimanSolver1D(Grid* grid);
};


#endif