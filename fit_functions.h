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

enum Models {
    boltzmann,
    expdecay,
    gaussian,
    hill,
    ic50,
    mm,
    modsin
};

typedef struct {
  double *xdata;
  double *ydata;
  int size;
  int len;
} FileData;


typedef struct {
    enum Models model;
    int varlen;  // number of parameters for given model
    int datalen; // number of data points
    double *x;   // an array of the independent variable from input
    double *y;   // an array of the dependent variable from input
    double *m;   // calculated dependent variable based on xi and VAR
    double *d1;  // partial derivative of var[0] wrt y for xi
    double *d2;  // ...
    double *d3;  // ..., not used for mm model
    double *d4;  // only used in gaussian model
} CFData;

void usage(char **argv);

/*
 * Prints a help message
 */

void readfile(char *filepath, FileData *data);

/*
 * Counts the number of data points in the file
 */

void equation(CFData *data, double *var);

/*
 * Calculates values for m(x, VAR)
 */

void derivatives(CFData *data, double *var);

/*
 * Calculates values for the partial derivatives of the
 * vars.
 */
void get_f(gsl_vector *fvect, CFData *data, double *var);

/*
 * Calculates f(x). f(xi, VAR) = yi - m(xi, VAR). Minimizing
 * $\frac{1}{2}\sum[f(x, VAR)]^2$ is the goal of the program
 */
void get_jac(double *jac, CFData *data, double *var);

/*
 * Calculates the Jacobian
 */

void solve_h(double *a, gsl_matrix *muImat, CFData *data,
             double *g, double *h);

/*
 * Solves the linear system (A + muI)h = -g
 */


double get_rho(gsl_vector *fvect, gsl_vector *newfvect, CFData *data,
               double mu, double *h, double *g);

/*
 * Calculates an updating parameter
 */

double get_mu(int vlen, double *a);

/*
 * Calculates a scaling parameter
 */


void levenberg_marquardt(CFData *data, double *var);

/*
 * The main routine. Finds values for VAR than minimize the
 * function $\frac{1}{2}\Sigma[f(x, VAR)]^2$
 */

void output(char *filename, CFData *data, double *var);

/*
 * Creates a gnuplot script to plot results. TODO: add -o option
 * to allow users to specify a path to save these files.
 */

