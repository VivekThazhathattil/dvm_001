#ifndef VORTEX_H
#define VORTEX_H

#include <complex.h>
#include <stddef.h>

#define EPSILON 1e-10

#ifndef PI
#define PI 3.14159265358979323846
#endif  

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
void vortex_advect(vortex_t* v, double complex vel, double dt);
void vortex_set_nascent_vortex_strength(vortex_t *v, double Us, double dt);
double vortex_get_nascent_vortex_dist_from_centre(vortex_t *v, double Us);
void vortex_set_nascent_vortex_position(vortex_t *v, double o_p_mj, double theta_s);
vortex_t* vortex_create_nascent_vortex(double theta, double Us, double dt);

#endif

