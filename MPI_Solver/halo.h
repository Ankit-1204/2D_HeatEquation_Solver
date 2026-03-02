#pragma once
#include <mpi.h>
#include <vector>
#include <functional>
#include "grid.h"
#include "field.h"

struct Neighbors { int left = MPI_PROC_NULL, right = MPI_PROC_NULL, down = MPI_PROC_NULL, up = MPI_PROC_NULL; };

inline double bc_left(double y)   { return 0.0; }
inline double bc_right(double y)  { return 0.0; }
inline double bc_bottom(double x) { return 0.0; }
inline double bc_top(double x)    { return 0.0; }

inline void exchange_halo(MPI_Comm cart_comm, const LocalGrid &lg, Field &u, const Neighbors &nb) {
    MPI_Status st;

    int nx = lg.Nx_loc;
    int ny = lg.Ny_loc;

    std::vector<double> send_left(ny), recv_left(ny), send_right(ny), recv_right(ny);
    for (int j = 1; j <= ny; ++j) {
        send_left[j-1]  = u.at(1, j); 
        send_right[j-1] = u.at(nx, j);   
    }

    // Exchange with left neighbor: send left column, receive neighbor's right interior column into recv_left
    if (nb.left != MPI_PROC_NULL) {
        MPI_Sendrecv(send_left.data(), ny, MPI_DOUBLE, nb.left, 0,
                     recv_left.data(), ny, MPI_DOUBLE, nb.left, 1,
                     cart_comm, &st);
    }
    // Exchange with right neighbor
    if (nb.right != MPI_PROC_NULL) {
        MPI_Sendrecv(send_right.data(), ny, MPI_DOUBLE, nb.right, 1,
                     recv_right.data(), ny, MPI_DOUBLE, nb.right, 0,
                     cart_comm, &st);
    }

    if (nb.left == MPI_PROC_NULL) {
        for (int j = 1; j <= ny; ++j) u.at(0, j) = bc_left(lg.global_y(j));
    } else {
        for (int j = 1; j <= ny; ++j) u.at(0, j) = recv_left[j-1];
    }
    if (nb.right == MPI_PROC_NULL) {
        for (int j = 1; j <= ny; ++j) u.at(nx+1, j) = bc_right(lg.global_y(j));
    } else {
        for (int j = 1; j <= ny; ++j) u.at(nx+1, j) = recv_right[j-1];
    }

    std::vector<double> send_down(nx), recv_down(nx), send_up(nx), recv_up(nx);
    for (int i = 1; i <= nx; ++i) {
        send_down[i-1] = u.at(i, 1);   
        send_up[i-1]   = u.at(i, ny);  
    }

    if (nb.down != MPI_PROC_NULL) {
        MPI_Sendrecv(send_down.data(), nx, MPI_DOUBLE, nb.down, 2,
                     recv_down.data(), nx, MPI_DOUBLE, nb.down, 3,
                     cart_comm, &st);
    }
    if (nb.up != MPI_PROC_NULL) {
        MPI_Sendrecv(send_up.data(), nx, MPI_DOUBLE, nb.up, 3,
                     recv_up.data(), nx, MPI_DOUBLE, nb.up, 2,
                     cart_comm, &st);
    }

    if (nb.down == MPI_PROC_NULL) {
        for (int i = 1; i <= nx; ++i) u.at(i, 0) = bc_bottom(lg.global_x(i));
    } else {
        for (int i = 1; i <= nx; ++i) u.at(i, 0) = recv_down[i-1];
    }
    if (nb.up == MPI_PROC_NULL) {
        for (int i = 1; i <= nx; ++i) u.at(i, ny+1) = bc_top(lg.global_x(i));
    } else {
        for (int i = 1; i <= nx; ++i) u.at(i, ny+1) = recv_up[i-1];
    }
}