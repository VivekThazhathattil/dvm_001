#ifndef VORTEX_H
#define VORTEX_H

#include <complex.h>
#include <stddef.h>

typedef struct VORTEX_s{
    double complex z;
    double gamma;
} vortex_t;

vortex_t vortex_create(double complex z, double gamma);
void vortex_set_position(vortex_t *v, double complex z);
void vortex_set_gamma(vortex_t *v, double gamma);
double complex vortex_get_position(vortex_t *v);
double vortex_get_gamma(vortex_t *v);

double complex vortex_get_velocity(vortex_t *v, double complex z, double c);

#endif

