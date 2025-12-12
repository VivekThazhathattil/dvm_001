#ifndef SEPARATION_POINT_H
#define SEPARATION_POINT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h> // For 'bool' type

#include "discrete_vortex_method.h" // For 'dvm_t' type

#ifndef PI
#define PI 3.14159265358979323846 // Definition of Pi
#endif

// Structure to hold all state variables and parameters
typedef struct {
    double c;            // Radius (R) of the cylinder
    double nu;           // Kinematic viscosity
    double Uinf;         // Freestream velocity
    double H_sep_thr;    // Shape-factor threshold for separation
    int N;               // Azimuthal resolution (number of points)
    double dphi;         // Angular step size
    double ds;           // Arc length step size

    double *phi;         // Angular position [rad]
    double *s;           // Arc length position [m]
    double *theta;      // Momentum thickness (current)
    double *theta_new;  // Momentum thickness (new)
    double *dtheta_ds;  // Spatial derivative of momentum thickness
    double *Ue;         // External flow velocity
    double *dUe_ds;     // Spatial derivative of Ue
    double *H;          // Shape factor
    double *Cf;         // Skin friction coefficient
    bool *sep;          // Separation flag
    double sep_angles_in_rad[2]; // Separation angles [rad]
    dvm_t *dvm_state;    // Discrete vortex model state
} sep_pt_t;

// Public function prototypes
sep_pt_t* sep_pt_module_init(double R, double nu, double Uinf, double H_sep_thr, int N, dvm_t *dvm);
void sep_pt_module_free(sep_pt_t *sep_pt);
void sep_pt_get_numerical_Ue(sep_pt_t *sep_pt);
void sep_pt_get_numerical_dUe_ds(sep_pt_t *sep_pt);
void sep_pt_compute_separation_points(sep_pt_t *sep_pt, double dt);
void sep_pt_get_sep_angles(sep_pt_t *sep_pt, double dt);

#endif // SEPARATION_POINT_H