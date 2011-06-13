//
//  fit_functions.h
//  curvefit
//
//  Created by Claude Rogers on 6/8/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlin.h>

struct cfdata {
    char *model;
    int varlen;
    int datalen;
    double *x;
    double *y;
    double *m;
    double *d1;
    double *d2;
    double *d3;
};

void usage(char **argv);
int filelines(char *filepath);
void equation(struct cfdata *data, double *var);
void derivatives(struct cfdata *data, double *var);
void get_f(gsl_vector *fvect, struct cfdata *data, double *var);
void get_jac(double *jac, struct cfdata *data, double *var);
void solve_h(double *a, gsl_matrix *muImat, struct cfdata *data, double *g, double *h);
double get_rho(gsl_vector *fvect, gsl_vector *newfvect, struct cfdata *data,
               double mu, double *h, double *g);
void levenberg_marquardt(struct cfdata *data, double *var);