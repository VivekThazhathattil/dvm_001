#ifndef SIMULATION_H
#define SIMULATION_H

#include "discrete_vortex_method.h"
#include <complex.h>
#include <stddef.h>

typedef struct SIMULATION_s{
  double U; // free-stream velocity
  double c; // cylinder radius
  double dt; // time step
  dvm_t *dvm;
  double x_min; // x-coordinate of the left-most point
  double x_max; // x-coordinate of the right-most point
  double y_min; // y-coordinate of the bottom-most point
  double y_max; // y-coordinate of the top-most point
  size_t nx; // number of grid points in the x-direction
  size_t ny; // number of grid points in the y-direction
  double *x; // x-coordinates of the grid points
  double *y; // y-coordinates of the grid points
  double **X; // x mesh grid
  double **Y; // y mesh grid
} sim_t;

sim_t* sim_create(double U, double c, double dt, double x_min, double x_max, double y_min, double y_max, size_t nx, size_t ny);
void sim_kill(sim_t *s);
void sim_create_meshgrid(sim_t *s);
double complex** sim_compute_complex_velocity_field(sim_t *s);
void write_vel_coords(double complex **u, sim_t *s);

#endif