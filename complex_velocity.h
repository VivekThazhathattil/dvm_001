#ifndef COMPLEX_VELOCITY_H
#define COMPLEX_VELOCITY_H

#include <complex.h>
#include <stddef.h>
#include "vortex.h"

double complex cv_uniform_cylinder(double complex z, double U, double c);
double complex cv_vortex_with_image(vortex_t *v, double complex z, double c);
void cv_sum_vortices_with_images(vortex_t *vs, size_t n, double complex z, double c, int self_idx, double complex *w_sum);
double complex complex_velocity(double complex z, double U, double c, vortex_t *vs, size_t n);

#endif

