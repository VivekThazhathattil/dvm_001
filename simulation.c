#include "simulation.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

sim_t* sim_create(double U, double c, double dt, double x_min, double x_max, double y_min, double y_max, size_t nx, size_t ny) {
    sim_t *s;
    s = (sim_t *) malloc(sizeof(sim_t));
    s->U = U;
    s->c = c;
    s->dt = dt;
    s->dvm = dvm_create(U, c, 0);
    s->x_min = x_min;
    s->x_max = x_max;
    s->y_min = y_min;
    s->y_max = y_max;
    s->nx = nx;
    s->ny = ny;
    s->x = NULL;
    s->y = NULL;
    s->X = NULL;
    s->Y = NULL;
    sim_create_meshgrid(s);
    return s;
}

void sim_kill(sim_t *s) {
    size_t i;
    dvm_kill(s->dvm);
    for(i = 0; i < s->nx; ++i) {
        free(s->X[i]);
    }
    free(s->X);
    for(i = 0; i < s->ny; ++i) {
        free(s->Y[i]);
    }
    free(s->Y);
    free(s);
}

void sim_create_meshgrid(sim_t *s) {
    double x, y;
    size_t i, j;

    s->x = (double *) malloc(s->nx * sizeof(double));
    s->y = (double *) malloc(s->ny * sizeof(double));

    for (i = 0; i < s->nx; ++i) {
        x = s->x_min + (s->x_max - s->x_min) * i / (s->nx - 1);
        s->x[i] = x;
    }

    for (i = 0; i < s->ny; ++i) {
        y = s->y_min + (s->y_max - s->y_min) * i / (s->ny - 1);
        s->y[i] = y;
    }

    s->X = (double **) malloc(s->ny * sizeof(double *));
    s->Y = (double **) malloc(s->ny * sizeof(double *));

    for (i = 0; i < s->ny; ++i) {
        s->X[i] = (double *) malloc(s->nx * sizeof(double));
        s->Y[i] = (double *) malloc(s->nx * sizeof(double));
        for (j = 0; j < s->nx; ++j) {
            s->X[i][j] = s->x[j];
            s->Y[i][j] = s->y[i];
        }
    }
}

double complex** sim_compute_complex_velocity_field(sim_t *s) {
    size_t i, j;
    double complex **velocity_field, z;
    velocity_field = (double complex **) malloc(s->ny * sizeof(double complex *));
    for (i = 0; i < s->ny; ++i) {
        velocity_field[i] = (double complex *) malloc(s->nx * sizeof(double complex));
    }
    for (i = 0; i < s->ny; ++i) {
        for (j = 0; j < s->nx; ++j) {
            z = s->X[i][j] + I * s->Y[i][j];
            velocity_field[i][j] = dvm_get_complex_velocity(z, s->dvm);
        }
    }  
    return velocity_field;
}

void write_vel_coords(double complex **u, sim_t *s) {
    size_t i, j;
    FILE *f;
    /* write u velocity, v velocity, X and Y coordinates to individual files as matrices in csv format */
    /* each file has ny rows and nx columns */
    f = fopen("unit_tests/u_velocity.csv", "w");
    for (i = 0; i < s->ny; ++i) {
        for (j = 0; j < s->nx; ++j) {
            fprintf(f, "%e, ", creal(u[i][j]));
        }
        fprintf(f, "\n");
    }
    fclose(f);

    f = fopen("unit_tests/v_velocity.csv", "w");
    for (i = 0; i < s->ny; ++i) {
        for (j = 0; j < s->nx; ++j) {
            fprintf(f, "%e, ", cimag(u[i][j]));
        }
        fprintf(f, "\n");
    }
    fclose(f);

    f = fopen("unit_tests/x_coordinates.csv", "w");
    for (i = 0; i < s->ny; ++i) {
        for (j = 0; j < s->nx; ++j) {
            fprintf(f, "%f, ", s->X[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);

    f = fopen("unit_tests/y_coordinates.csv", "w");
    for (i = 0; i < s->ny; ++i) {
        for (j = 0; j < s->nx; ++j) {
            fprintf(f, "%f, ", s->Y[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int main() {
    sim_t *s;
    double complex **u;

    s = sim_create(1.0, 1.0, 0.125, -4, 12, -4, 4, 200, 160);
    u = sim_compute_complex_velocity_field(s);
    write_vel_coords(u, s);    
    return 0;
}