#ifndef DISCRETE_VORTEX_METHOD_H
#define DISCRETE_VORTEX_METHOD_H

#include "complex_velocity.h"
#include "vortex.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef struct DVM_s{
  double U; // free-stream velocity
  double c; // cylinder radius
  size_t num_vorts; // number of vortices
  vortex_t *vs; // array of vortices
  double dt; // time step
} dvm_t;

dvm_t* dvm_create(double U, double c, size_t num_vorts);
void dvm_kill(dvm_t* dvm);
double complex dvm_get_complex_velocity(double complex z, dvm_t *dvm);
void dvm_append_vortex(dvm_t* dvm, vortex_t* v);
void dvm_debug_put_arbitrary_vortices(dvm_t *dvm, size_t num_vorts, \
    double x_min, double x_max, double y_min, double y_max);
void dvm_advect_vortices(dvm_t *dvm, double dt);
void dvm_create_nascent_vortices(dvm_t* dvm, double* theta);
double complex dvm_pos_from_angle(double c,double theta);
#endif
