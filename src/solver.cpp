#include "solver.h"

RimanSolver1D::RimanSolver1D(Grid* grid)
{

}

RimanSolver1D::~RimanSolver1D()
{
}

bool RimanSolver1D::check_grid(const Grid &grid) 
{
    if(grid.cells_num_y > 1)
        throw std::runtime_error("Not 1d problem. Use another solver.");
}