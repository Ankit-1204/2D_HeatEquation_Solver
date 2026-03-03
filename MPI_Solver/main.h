#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <string>
#include "grid.h"
#include "field.h"
#include "halo.h"
#include "operator.h"
#include "cg.h"

using namespace std;

// Example source (zero here). Replace if testing manufactured solution.
inline double default_source(double x,double y,double t){ return 0.0; }

// Simple CSV writer per rank for visualization
static void write_csv_local(const LocalGrid &lg, const Field &u, const string &prefix, int rank, int step) {
    string fname = prefix + "_rank_" + to_string(rank) + "_step_" + to_string(step) + ".csv";
    ofstream f(fname);
    f << "x,y,u\n";
    for (int j = 1; j <= lg.Ny_loc; ++j) {
        for (int i = 1; i <= lg.Nx_loc; ++i) {
            double x = lg.global_x(i);
            double y = lg.global_y(j);
            f << x << "," << y << "," << u.atc(i,j) << "\n";
        }
    }
    f.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, nprocs; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Problem parameters (edit as needed)
    int Nx = 100, Ny = 80;                      // global interior points
    double x0=0.0, x1=1.0, y0=0.0, y1=1.0;
    double alpha = 1.0;
    double dt = 1e-4;
    double t_final = 0.01;
    double theta = 0.5;                         // 1.0 -> BE, 0.5 -> CN
    bool use_jacobi = true;
    int out_every = 20;
    double tol = 1e-7;
    int maxit = 2000;

    int Px = (int)floor(sqrt((double)nprocs));
    while (Px > 1 && (nprocs % Px) != 0) --Px;
    if (Px == 0) Px = 1;
    int Py = nprocs / Px;
    if (rank == 0) cout << "Px="<<Px<<" Py="<<Py<<" procs="<<nprocs<<"\n";

    if (Nx % Px != 0 || Ny % Py != 0) {
        if (rank == 0) cerr << "ERROR: Require Nx divisible by Px and Ny divisible by Py for this simple partitioning.\n";
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int nx_loc = Nx / Px;
    int ny_loc = Ny / Py;

    int dims[2] = {Px, Py};
    int periods[2] = {0,0};
    MPI_Comm cart; MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart);

    int coords[2]; MPI_Cart_coords(cart, rank, 2, coords);
    int px = coords[0], py = coords[1];

    Neighbors nb;
    MPI_Cart_shift(cart, 0, 1, &nb.left, &nb.right); 
    MPI_Cart_shift(cart, 1, 1, &nb.down, &nb.up);   

    int gx = px * nx_loc; 
    int gy = py * ny_loc; 

    double hx = (x1 - x0) / (Nx + 1);
    double hy = (y1 - y0) / (Ny + 1);

    LocalGrid lg;
    lg.init(x0, x1, y0, y1, nx_loc, ny_loc, gx, gy);
    lg.hx = hx; lg.hy = hy;

    Field u_old(&lg), u_new(&lg), b(&lg);
    u_old.reset(&lg); u_new.reset(&lg); b.reset(&lg);

    // Initial condition: Gaussian
    double cx = 0.5*(x0+x1), cy = 0.5*(y0+y1);
    double sigma = 0.08;
    for (int j = 1; j <= lg.Ny_loc; ++j)
        for (int i = 1; i <= lg.Nx_loc; ++i) {
            double x = lg.global_x(i), y = lg.global_y(j);
            u_old.at(i,j) = exp(-((x-cx)*(x-cx) + (y-cy)*(y-cy)) / (2.0*sigma*sigma));
        }

    exchange_halo(cart, lg, u_old, nb);

    double h = hx; // assume hx==hy for gamma; generalize if needed
    double gamma = alpha * dt / (h*h);

    // function<void(const LocalGrid&, const Field&, Field&)> precond = nullptr;
    // if (use_jacobi) {
    //     precond = [&](const LocalGrid &L, const Field &r, Field &z) {
    //         double diag = 1.0 + theta * gamma * 4.0;
    //         double inv = 1.0 / diag;
    //         for (int j=1;j<=L.Ny_loc;++j) for (int i=1;i<=L.Nx_loc;++i) z.at(i,j) = inv * r.atc(i,j);
    //     };
    // }

    auto applyM = [&](const LocalGrid &L, const Field &x, Field &y) {
        apply_M_local(L, x, y, theta, gamma);
    };

    int nsteps = (int)ceil(t_final / dt);
    if (rank==0) cout << "nsteps="<<nsteps<<" dt="<<dt<<" theta="<<theta<<" gamma="<<gamma<<"\n";

    for (int step=0; step < nsteps; ++step) {
        double t = step * dt;
        double tnext = t + dt;

        exchange_halo(cart, lg, u_old, nb);

        // Build RHS b based on theta-method
        if (theta == 1.0) {
            // Backward Euler: b = u^n + dt * s^{n+1}
            for (int j=1;j<=lg.Ny_loc;++j) for (int i=1;i<=lg.Nx_loc;++i) {
                double x = lg.global_x(i), y = lg.global_y(j);
                b.at(i,j) = u_old.atc(i,j) + dt * default_source(x,y,tnext);
            }
        } else if (fabs(theta - 0.5) < 1e-12) {
            // Crank-Nicolson: b = (I - gamma/2 * A) u^n + dt*(s^{n+1}+s^n)/2
            Field Au(&lg); Au.reset(&lg);
            apply_A_local(lg, u_old, Au);
            for (int j=1;j<=lg.Ny_loc;++j) for (int i=1;i<=lg.Nx_loc;++i) {
                double x = lg.global_x(i), y = lg.global_y(j);
                b.at(i,j) = u_old.atc(i,j) - 0.5 * gamma * Au.atc(i,j) + dt * 0.5 * (default_source(x,y,t) + default_source(x,y,tnext));
            }
        } else {
            // General theta-method: b = (I - (1-theta)*gamma*A)u^n + dt*(theta s^{n+1} + (1-theta) s^n)
            Field Au(&lg); Au.reset(&lg);
            apply_A_local(lg, u_old, Au);
            for (int j=1;j<=lg.Ny_loc;++j) for (int i=1;i<=lg.Nx_loc;++i) {
                double x = lg.global_x(i), y = lg.global_y(j);
                b.at(i,j) = u_old.atc(i,j) - (1.0 - theta) * gamma * Au.atc(i,j) + dt * (theta * default_source(x,y,tnext) + (1.0-theta) * default_source(x,y,t));
            }
        }

        // initial guess = u_old
        copy_local(lg, u_old, u_new);
        exchange_halo(cart, lg, u_new, nb);

        // Solve M * u_new = b
        CGResult cres = pcg_solve_mpi(cart, lg, nb,
                                      [&](const LocalGrid &L, const Field &x, Field &y) {
                                          applyM(L, x, y); 
                                      },
                                       nullptr,
                                      b, u_new, maxit, tol);

        if (rank == 0 && step % out_every == 0) {
            cout << "step="<<step<<" iters="<<cres.iters<<" res="<<cres.resnorm<<"\n";
        }

        // move solution to u_old for next step (interior copy only)
        copy_local(lg, u_new, u_old);

        if (step % out_every == 0) write_csv_local(lg, u_old, "u", rank, step);
    }

    if (rank == 0) cout << "Done\n";
    MPI_Finalize();
    return 0;
}