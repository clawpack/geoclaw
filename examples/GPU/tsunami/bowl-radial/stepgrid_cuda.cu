#include "real.H"
#include "params.H"
#include "shallow_topo_riemann_solver.H"
#include "fused_solvers_headers.H"
#include "advance.H"
#include <cuda_runtime.h>

extern "C" void call_C_limited_riemann_update(
        const int cellsX, const int cellsY, const int ghostCells,
        const real startX, const real endX, const real startY, const real endY,
        const real dt,
        real* q, real* qNew, 
        real* coefficients,
        real* waveSpeedsX, real* waveSpeedsY,
        const int numStates, const int numCoefficients,
        real* cfls, const int ngrids, const int mcapa,
        const int id, const int dev_id) {

    // actually qNew holds the input old solution as well as new output solution
    // q is just temporary storage for intermediate solution between x-sweep and y-sweep

    cudaStream_t stream;

    get_cuda_stream(id, dev_id, &stream);

    pdeParam param(cellsX, cellsY, ghostCells, 
            numStates, NWAVES, numCoefficients,
            startX, endX, startY, endY, dt,
            q, qNew, 
            coefficients, 
            waveSpeedsX, waveSpeedsY,
            cfls, mcapa, id, dev_id);

    param.setOrderOfAccuracy(2);

    shallow_topo_horizontal swe_h;
    shallow_topo_vertical swe_v;
    
    limiter_MC phi;
    // limiter_VanLeer phi;

    limited_Riemann_Update(param, 
            swe_h, swe_v, 
            phi,stream);

}

