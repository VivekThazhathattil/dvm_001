#include "separation_point.h"

/**
 * @brief Thwaites' integral closure relations to link the Thwaites parameter (lambda)
 * to the Shape Factor (H) and the skin friction parameter (l).
 */
static void thwaites_closure(double lambda, double *H_out, double *l_out) {
    double H, l;

    // Skin friction parameter 'l' closure
    l = 0.45 - 6.0 * lambda;

    // Shape factor 'H' closure (piecewise approximations)
    if (lambda >= 0.0 && lambda <= 0.1) {
        H = 2.61 - 3.75 * lambda + 5.24 * lambda * lambda;
    } else if (lambda < 0.0 && lambda >= -0.1) {
        H = 2.088 + 0.0731 * lambda + 0.14 * lambda * lambda;
    } else if (lambda > 0.1) {
        H = 2.61 - 3.75 * 0.1 + 5.24 * 0.1 * 0.1;
        // Approximation for larger positive lambda
        H -= 0.5 * (lambda - 0.1) / (1.0 + lambda);
        if (H < 1.4) H = 1.4;
    } else { // lambda < -0.1
        // Approximation for large negative lambda (separation/reversed flow)
        H = 3.5 + 10.0 * (-(lambda + 0.1));
        if (H > 10.0) H = 10.0;
    }

    *H_out = H;
    *l_out = l;
}

static void compute_derivatives_periodic(int N, const double *arr, double dx, double *darr_ds) {
    int i, ip, im;
    double dF;
    for (i = 0; i < N; ++i) {
        // Periodic boundary indices
        ip = (i + 1) % N;       // i+1
        im = (i - 1 + N) % N;   // i-1 (handle wrap-around for i=0)

        // Central difference scheme: (f(i+1) - f(i-1)) / (2 * dx)
        dF = (arr[ip] - arr[im]) / (2.0 * dx);
        darr_ds[i] = dF;
    }
}


/**
 * @brief Initializes and allocates memory for the sep_pt_t structure.
 */
sep_pt_t* sep_pt_module_init(double R, double nu, double Uinf, double H_sep_thr, int N, dvm_t *dvm) {
    int i;
    sep_pt_t *sep_pt;

    // Use R (radius) for variable c in the struct
    sep_pt = (sep_pt_t*)malloc(sizeof(sep_pt_t));
    if (!sep_pt) {
        fprintf(stderr, "ERR: sep_pt: sep_pt_module_init: Allocation failed for struct\n");
        return NULL;
    }
    
    sep_pt->c = R;
    sep_pt->nu = nu;
    sep_pt->Uinf = Uinf;
    sep_pt->H_sep_thr = H_sep_thr;
    sep_pt->N = N;
    sep_pt->dphi = 2.0 * PI / (double)N;
    sep_pt->ds = sep_pt->c * sep_pt->dphi;
    sep_pt->dvm_state = dvm;

    // Memory allocations
    sep_pt->phi = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->s = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->theta = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->theta_new = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->dtheta_ds = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->Ue = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->dUe_ds = (double*)malloc(sep_pt->N * sizeof(double)); 
    sep_pt->H = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->Cf = (double*)malloc(sep_pt->N * sizeof(double));
    sep_pt->sep = (bool*)malloc(sep_pt->N * sizeof(bool));

    // Check for allocation failures and free memory if any failed (CRITICAL FIX)
    if (!sep_pt->phi || !sep_pt->s || !sep_pt->theta || !sep_pt->theta_new || !sep_pt->dtheta_ds || !sep_pt->Ue || !sep_pt->dUe_ds || !sep_pt->H || !sep_pt->Cf || !sep_pt->sep) {
        fprintf(stderr, "ERR: sep_pt: Array allocation failed\n");
        sep_pt_module_free(sep_pt); // Use the free function to clean up
        return NULL;
    }

    for (i = 0; i < sep_pt->N; ++i) {
        sep_pt->phi[i] = i * sep_pt->dphi;
        sep_pt->s[i] = sep_pt->c * sep_pt->phi[i];
        sep_pt->theta[i] = 1e-4; // Initialize momentum thickness
    }

    // Initialize Ue and dUe_ds based on potential flow (used for first step)
    sep_pt_get_numerical_Ue(sep_pt);
    sep_pt_get_numerical_dUe_ds(sep_pt);
    compute_derivatives_periodic(sep_pt->N, sep_pt->theta, sep_pt->ds, sep_pt->dtheta_ds);

    return sep_pt;
}

/**
 * @brief Frees all dynamically allocated memory in the sep_pt_t structure.
 */
void sep_pt_module_free(sep_pt_t *sep_pt) {
    if (sep_pt) {
        free(sep_pt->phi);
        free(sep_pt->s);
        free(sep_pt->theta);
        free(sep_pt->theta_new);
        free(sep_pt->dtheta_ds);
        free(sep_pt->Ue);
        free(sep_pt->dUe_ds);
        free(sep_pt->H);
        free(sep_pt->Cf);
        free(sep_pt->sep);
        free(sep_pt);
    }
}

void sep_pt_get_numerical_Ue(sep_pt_t *sep_pt) {
    int i;
    double R;
    double complex z, V;

    R = sep_pt->c;

    for (i = 0; i < sep_pt->N; ++i) {
        z = R * cos(sep_pt->phi[i]) + I * R * sin(sep_pt->phi[i]);
        V = dvm_get_complex_velocity(z, sep_pt->dvm_state);
        sep_pt->Ue[i] = cabs(V); 
    }
}

void sep_pt_get_numerical_dUe_ds(sep_pt_t *sep_pt) {
    compute_derivatives_periodic(sep_pt->N, sep_pt->Ue, sep_pt->ds, sep_pt->dUe_ds);
}

void sep_pt_compute_separation_points(sep_pt_t *sep_pt, double dt) {
    int i;
    double Hloc, lloc;
    double lambda, Cfloc;
    double dtheta_dt;
    bool separated;

    if (sep_pt == NULL) return;
    
    // --- Logic Fix: Call Ue updates only once per step ---
    // Ue and dUe_ds are updated based on the (simplified) external flow assumption
    sep_pt_get_numerical_Ue(sep_pt);
    sep_pt_get_numerical_dUe_ds(sep_pt);

    // Compute the spatial derivative of theta
    compute_derivatives_periodic(sep_pt->N, sep_pt->theta, sep_pt->ds, sep_pt->dtheta_ds);

    /*--- Step 1: Compute H, Cf, and Separation Status ---*/
    for (i = 0; i < sep_pt->N; ++i) {
        lambda = (sep_pt->theta[i] * sep_pt->theta[i] / sep_pt->nu) * sep_pt->dUe_ds[i];
        
        thwaites_closure(lambda, &Hloc, &lloc);

        // Calculate Cf: Cf = 2 * l * nu / (Ue * theta)
        Cfloc = 0.0;
        if (sep_pt->Ue[i] > 1e-12 && sep_pt->theta[i] > 1e-12) {
            Cfloc = 2.0 * lloc * sep_pt->nu / (sep_pt->Ue[i] * sep_pt->theta[i]);
        } else {
            Cfloc = 0.0;
        }

        sep_pt->H[i] = Hloc;
        sep_pt->Cf[i] = Cfloc;

        // Check for separation
        #if DEBUG_OLD
            fprintf(stdout, "%f,",Cfloc);
        #endif
        separated = false;
        if (Cfloc < 0.0) separated = true;
        if (Hloc > sep_pt->H_sep_thr) separated = true;
        sep_pt->sep[i] = separated;
    }

    /*--- Step 2: Time Integration (Explicit Euler) ---*/
    for (i = 0; i < sep_pt->N; ++i) {
        /* Unsteady Momentum Integral Equation:
         * dtheta/dt = - Ue * dtheta/ds - ((2 + H)/2) * theta * dUe/ds + (Cf * Ue / 2)
         */
        dtheta_dt = - sep_pt->Ue[i] * sep_pt->dtheta_ds[i] 
                    - ((2.0 + sep_pt->H[i]) / 2.0) * sep_pt->theta[i] * sep_pt->dUe_ds[i] 
                    + 0.5 * sep_pt->Cf[i] * sep_pt->Ue[i];

        sep_pt->theta_new[i] = sep_pt->theta[i] + dt * dtheta_dt;

        // Keep theta positive and not unphysically small
        if (sep_pt->theta_new[i] < 1e-12) sep_pt->theta_new[i] = 1e-12;
    }

    /*--- Step 3: Rotate theta <- theta_new (Update state) ---*/
    for (i = 0; i < sep_pt->N; ++i) {
        sep_pt->theta[i] = sep_pt->theta_new[i];
    }
}

void sep_pt_find_sep_angles(sep_pt_t *sep_pt) {
    int i;
    bool curr_val, next_val;

    /*--- Assumption: flow cannot separate at the leading edge ---*/

    curr_val = sep_pt->sep[0];
    for(i = 1; i < sep_pt->N; ++i) {
        next_val = sep_pt->sep[i];
        if (next_val != curr_val) {
            sep_pt->sep_angles_in_rad[0] = sep_pt->phi[i];
            break;
        }
    }
    /*
    curr_val = sep_pt->sep[sep_pt->N-1];
    for(i = sep_pt->N - 2; i > 0; --i) {
        next_val = sep_pt->sep[i];
        if (next_val != curr_val) {
            sep_pt->sep_angles_in_rad[1] = sep_pt->phi[i];
            break;
        }
    }
    */

   /* 
   I'm taking the second separation point to be symmetrically opposite to the first 
   one because I'm not seeing the second separation point converging properly
   */
   sep_pt->sep_angles_in_rad[1] = -sep_pt->sep_angles_in_rad[0];
} 

void sep_pt_get_sep_angles(sep_pt_t *sep_pt, double dt) {
    sep_pt_compute_separation_points(sep_pt, dt);
    sep_pt_find_sep_angles(sep_pt);
    #if DEBUG_OLD
    /* Print difference of Ue and Ue analytical from potential flow solution*/
        for(int i = 0; i < sep_pt->N; ++i) {
            fprintf(stdout, "%f,", sep_pt->Ue[i] - 2*sep_pt->Uinf*sin(sep_pt->phi[i]));
        }
        fprintf(stdout, "\n");
    #endif
    #if DEBUG_OLD
    /* Print separation angles*/
        fprintf(stdout, "Separation Angles (in deg): %f, %f\n", sep_pt->sep_angles_in_rad[0] * 180.0 / PI, sep_pt->sep_angles_in_rad[1] * 180.0 / PI);
    #endif
    #if DEBUG
    /* Print separation angles without additional text*/
        fprintf(stdout, "%f, %f\n", sep_pt->sep_angles_in_rad[0] * 180.0 / PI, sep_pt->sep_angles_in_rad[1] * 180.0 / PI);
    #endif
    #if DEBUG_OLD
        fprintf(stdout, "%f\n", sep_pt->sep_angles_in_rad[0] * 180.0 / PI);
        fprintf(stdout, "%f\n", sep_pt->sep_angles_in_rad[0] * 180.0 / PI);
        fprintf(stdout, "\n");
    #endif
    #if DEBUG_OLD
    /* Print separation status*/
        for(int i = 0; i < sep_pt->N; ++i) {
            fprintf(stdout, "%d", sep_pt->sep[i]);
        }
        fprintf(stdout, "\n");
    #endif
    #if DEBUG_OLD
        for(int i = 0; i < sep_pt->N; ++i) {
            fprintf(stdout, "%f,", sep_pt->H[i]);
        }
    #endif
    #if DEBUG_OLD
    /* Print angles in degrees*/
        for(int i = 0; i < sep_pt->N; ++i) {
            fprintf(stdout, "%f,", sep_pt->phi[i] * 180 / PI);
        }
        fprintf(stdout, "\n");
    #endif
}