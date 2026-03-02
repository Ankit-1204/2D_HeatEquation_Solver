#pragma once
#include <mpi.h>
#include <functional>
#include <cmath>
#include "grid.h"
#include "field.h"
#include "halo.h"

// Local linear algebra helpers (interior only)
inline double dot_local(const LocalGrid &lg, const Field &a, const Field &b) {
    double s = 0.0;
    for (int j = 1; j <= lg.Ny_loc; ++j)
        for (int i = 1; i <= lg.Nx_loc; ++i) s += a.atc(i,j) * b.atc(i,j);
    return s;
}
inline void axpy_local(const LocalGrid &lg, double a, const Field &x, Field &y) {
    for (int j = 1; j <= lg.Ny_loc; ++j)
        for (int i = 1; i <= lg.Nx_loc; ++i) y.at(i,j) += a * x.atc(i,j);
}
inline void copy_local(const LocalGrid &lg, const Field &src, Field &dst) {
    for (int j = 1; j <= lg.Ny_loc; ++j)
        for (int i = 1; i <= lg.Nx_loc; ++i) dst.at(i,j) = src.atc(i,j);
}

struct CGResult { int iters; double resnorm; };

// global reduction
inline double allreduce_sum(MPI_Comm comm, double local) {
    double g; MPI_Allreduce(&local, &g, 1, MPI_DOUBLE, MPI_SUM, comm);
    return g;
}

// Matrix-free PCG solver (MPI-aware).
// - comm: Cartesian communicator
// - lg: local grid
// - nb: neighbor ranks
// - applyM: function(localGrid, x, y) that computes y = M*x (reads ghosts of x)
// - precond: function(localGrid, r, z) computing z = M^{-1} r (or nullptr for none)
// - b, x: local Fields (ghosts should be managed by caller when required)
inline CGResult pcg_solve_mpi(MPI_Comm comm, const LocalGrid &lg, const Neighbors &nb,
                              std::function<void(const LocalGrid&, const Field&, Field&)> applyM,
                              std::function<void(const LocalGrid&, const Field&, Field&)> precond,
                              Field &b, Field &x, int maxit, double tol)
{
    int rank; MPI_Comm_rank(comm, &rank);

    Field r(&lg); r.reset(&lg);
    Field z(&lg); z.reset(&lg);
    Field p(&lg); p.reset(&lg);
    Field Ap(&lg); Ap.reset(&lg);
    Field tmp(&lg); tmp.reset(&lg);

    applyM(lg, x, tmp);
    for (int j=1;j<=lg.Ny_loc;++j) for (int i=1;i<=lg.Nx_loc;++i) r.at(i,j) = b.atc(i,j) - tmp.atc(i,j);

    if (precond) precond(lg, r, z); else copy_local(lg, r, z);
    copy_local(lg, z, p);

    double rz_local = dot_local(lg, r, z);
    double rz = allreduce_sum(comm, rz_local);
    double bnorm_local = dot_local(lg, b, b);
    double bnorm = sqrt(allreduce_sum(comm, bnorm_local));
    if (bnorm == 0.0) bnorm = 1.0;
    double res0 = sqrt(rz);
    if (res0 <= tol * bnorm) return {0, res0};

    double rz_old = rz;

    for (int iter = 1; iter <= maxit; ++iter) {
        // Ap = M * p   (ensure p halos are current)
        exchange_halo(comm, lg, p, nb);
        applyM(lg, p, Ap);

        double pAp_local = dot_local(lg, p, Ap);
        double pAp = allreduce_sum(comm, pAp_local);
        if (fabs(pAp) < 1e-24) {
            double rr_local = dot_local(lg, r, r);
            double rr = allreduce_sum(comm, rr_local);
            return {iter, sqrt(rr)};
        }

        double alpha = rz_old / pAp;

        axpy_local(lg, alpha, p, x); // x += alpha * p
        exchange_halo(comm, lg, x, nb);

        axpy_local(lg, -alpha, Ap, r); // r -= alpha*Ap

        double rr_local = dot_local(lg, r, r);
        double rr = allreduce_sum(comm, rr_local);
        double res = sqrt(rr);
        if (res <= tol * bnorm) return {iter, res};

        if (precond) precond(lg, r, z); else copy_local(lg, r, z);

        double rz_new_local = dot_local(lg, r, z);
        double rz_new = allreduce_sum(comm, rz_new_local);

        double beta = rz_new / rz_old;

        for (int j = 1; j <= lg.Ny_loc; ++j)
            for (int i = 1; i <= lg.Nx_loc; ++i)
                p.at(i,j) = z.atc(i,j) + beta * p.atc(i,j);

        rz_old = rz_new;
    }

    double rr_local = dot_local(lg, r, r);
    double rr = allreduce_sum(comm, rr_local);
    return {maxit, sqrt(rr)};
}