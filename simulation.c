#include "simulation.h"
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

sim_config_t* sim_read_config(const char *filename) {
    char line[256];
    char key[64], value[64];

    sim_config_t *config = (sim_config_t *)malloc(sizeof(sim_config_t));

    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open config file %s\n", filename);
        free(config);
        return NULL;
    }
    
    while (fgets(line, sizeof(line), f) != NULL) {
        if (sscanf(line, "%63[^=]=%63s", key, value) == 2) {
            if (strcmp(key, "U") == 0) config->U = atof(value);
            else if (strcmp(key, "c") == 0) config->c = atof(value);
            else if (strcmp(key, "dt") == 0) config->dt = atof(value);
            else if (strcmp(key, "num_vorts") == 0) config->num_vorts = (size_t)atoi(value);
            else if (strcmp(key, "x_min") == 0) config->x_min = atof(value);
            else if (strcmp(key, "x_max") == 0) config->x_max = atof(value);
            else if (strcmp(key, "y_min") == 0) config->y_min = atof(value);
            else if (strcmp(key, "y_max") == 0) config->y_max = atof(value);
            else if (strcmp(key, "nx") == 0) config->nx = (size_t)atoi(value);
            else if (strcmp(key, "ny") == 0) config->ny = (size_t)atoi(value);
            else if (strcmp(key, "num_iters") == 0) config->niters = (size_t)atoi(value);
        }
    }
    fclose(f);
    return config;
}

sim_t* sim_create(double U, double c, double dt,\
     double x_min, double x_max, double y_min, double y_max,\
      size_t nx, size_t ny, size_t niters) {
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
    s->niters = niters;
    sim_create_meshgrid(s);
    return s;
}

void sim_kill(sim_t *s) {
    size_t i;
    dvm_kill(s->dvm);
    for(i = 0; i < s->ny; ++i) {
        free(s->X[i]);
    }
    free(s->X);
    for(i = 0; i < s->ny; ++i) {
        free(s->Y[i]);
    }
    free(s->Y);
    free(s->x);
    free(s->y);
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

void sim_write_double_field_to_csv(const char *filename, double **field, size_t ny, size_t nx) {
    FILE *f = fopen(filename, "w");
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            fprintf(f, "%f, ", field[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void sim_write_complex_field_real_part_to_csv(const char *filename, double complex **field, size_t ny, size_t nx) {
    FILE *f = fopen(filename, "w");
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            fprintf(f, "%e, ", creal(field[i][j]));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void sim_write_complex_field_imag_part_to_csv(const char *filename, double complex **field, size_t ny, size_t nx) {
    FILE *f = fopen(filename, "w");
    for (size_t i = 0; i < ny; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            fprintf(f, "%e, ", cimag(field[i][j]));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void sim_write_vel_coords(double complex **u, sim_t *s, int timestep) {
    char filename[256];
    sprintf(filename, "unit_tests/u_velocity_%d.csv", timestep);
    sim_write_complex_field_real_part_to_csv(filename, u, s->ny, s->nx);
    sprintf(filename, "unit_tests/v_velocity_%d.csv", timestep);
    sim_write_complex_field_imag_part_to_csv(filename, u, s->ny, s->nx);
    sprintf(filename, "unit_tests/x_coordinates_%d.csv", timestep);
    sim_write_double_field_to_csv(filename, s->X, s->ny, s->nx);
    sprintf(filename, "unit_tests/y_coordinates_%d.csv", timestep);
    sim_write_double_field_to_csv(filename, s->Y, s->ny, s->nx);
}

void sim_debug_put_arbitrary_vortices(dvm_t *dvm, size_t num_vorts, \
    double x_min, double x_max, double y_min, double y_max) {
    dvm_debug_put_arbitrary_vortices(dvm, num_vorts, x_min, x_max, y_min, y_max);
}

void sim_free_complex_field(double complex **field, size_t ny) {
    for (size_t i = 0; i < ny; ++i) {
        free(field[i]);
    }
    free(field);
}

void sim_advect_vortices(sim_t *s) {
    dvm_advect_vortices(s->dvm, s->dt);
}

int main() {
    sim_t *s;
    sep_pt_t *sep_pt;
    double complex **u;
    sim_config_t *config;

    int N_azim;
    double nu, H_sep_thr;

    N_azim = 400; // number of azimuthal points for the cylinder
    nu = 1.5e-5; // kinematic viscosity
    H_sep_thr = 2.7; // shape-factor threshold for separation

    config = sim_read_config("setup.config");
    if (config == NULL) return 1;

    s = sim_create(config->U, config->c, config->dt, config->x_min, config->x_max, 
                   config->y_min, config->y_max, config->nx, config->ny, config->niters);
    sim_debug_put_arbitrary_vortices(s->dvm, config->num_vorts, config->x_min, config->x_max, 
                                     config->y_min, config->y_max);

    sep_pt = sep_pt_module_init(s->c, nu, s->U, H_sep_thr, N_azim, s->dvm);

    for (int i = 0; i < s->niters; ++i) {
        #if DEBUG_OLD
            fprintf(stdout, "Simulation Timestep %d\n", i);
        #endif
        sep_pt_get_sep_angles(sep_pt, s->dt);
        //u = sim_compute_complex_velocity_field(s);
        //sim_write_vel_coords(u, s, i);
        //sim_free_complex_field(u, s->ny);
        sim_advect_vortices(s);
    }

    sim_kill(s);
    sep_pt_module_free(sep_pt);
    free(config);
    return 0;
}