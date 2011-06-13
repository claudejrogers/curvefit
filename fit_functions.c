//
//  fit_functions.c
//  curvefit
//
//  Created by Claude Rogers on 6/8/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include "fit_functions.h"

#define MAX_ITER 500
#define NORM_G 1e-10
#define NORM_H 1e-15

void usage(char **argv)
{
    printf("\nUsage: %s -f <filename> -m <model> [-x <var1> -y <var2> -z <var2>]\n"
           "\nFit experimental data to either an ic50 or Michaelis-Menten model\n\n"
           "INPUT:\nRequired. Provide the filepath to the input data.\n"
           "\t-f FILENAME\tFile must contain tab delimited x and y data\n"
           "MODEL:\nRequired. Enter model name.\n"
           "\t-m MODEL\tSupported models are \"ic50\" and \"mm\"\n"
           "INITIAL GUESS:\nOptional. Provide values for fit parameters. "
           "Default values are 1.0.\n\t-x VALUE\tValue for max-min or Vmax.\n"
           "\t-y VALUE\tValue for ic50 or Km.\n\t-z VALUE\tValue for Hill "
           "coefficient.\n\nOUTPUT:\nUpdated values for the variables are displayed"
           " for each iteration,\nalong with |f(x)|.\n\n", argv[0]);
}

int filelines(char *filepath)
{
    FILE *fp;
    fp = fopen(filepath, "r");
    int c, nl;
    nl = 0;
    while ((c = fgetc(fp)) != EOF)
        if (c == '\t')
            nl++;
    fclose(fp);
    return nl;
}

void equation(struct cfdata *data, double *var)
{
    int i;
    for (i = 0; i < data->datalen; i++) {
        if (!strcmp(data->model, "ic50"))
            data->m[i] = (1 - (var[0]/(1 + pow((var[1]/data->x[i]), var[2]))));
        else
            data->m[i] = ((var[0] * data->x[i]) / (var[1] + data->x[i]));
    }
}

void derivatives(struct cfdata *data, double *var)
{
    int i;
    for (i = 0; i < data->datalen; i++) {
        if (!strcmp(data->model, "ic50")) {
            data->d1[i] = (-(1/(1 + pow((var[1]/data->x[i]), var[2]))));
            data->d2[i] = ((var[0] * var[2] * pow((var[1]/data->x[i]), 
                                                  (var[2] - 1))) / 
                           (data->x[i] * pow((1 + (pow((var[1]/data->x[i]), 
                                                       var[2]))), 2)));
            data->d3[i] = ((var[0] * pow((var[1]/data->x[i]), var[2]) * 
                            log((var[1]/data->x[i]))) / 
                           (pow((1 + pow((var[1]/data->x[i]), var[2])), 2)));
        } else {
            data->d1[i] = (data->x[i] / (var[1] + data->x[i]));
            data->d2[i] = (-(var[0] * data->x[i]) / pow((var[1] + data->x[i]), 2.0));
        }
    }
}

void get_f(gsl_vector *fvect, struct cfdata *data, double *var)
{
    equation(data, var);
    double fv;
    int i;
    for (i = 0; i < data->datalen; i++) {
        fv = data->y[i] - data->m[i];
        gsl_vector_set(fvect, i, fv);
    }
}

void get_jac(double *jac, struct cfdata *data, double *var)
{
    derivatives(data, var);
    int i, j, k;
    i = 0;
    j = data->datalen;
    k = 2 * data->datalen;
    for (i = 0; i < data->datalen; i++, j++, k++) {
        if (!strcmp(data->model, "ic50")) {
            jac[i] = -data->d1[i];
            jac[j] = -data->d2[i];
            jac[k] = -data->d3[i];
        } else {
            jac[i] = -data->d1[i];
            jac[j] = -data->d2[i];
        }
    }
}

void solve_h(double *a, gsl_matrix *muImat, struct cfdata *data, double *g, double *h)
{
    int i, s, l, dl;
    l = data->varlen;
    dl = l * l;
    double *na, *ng;

    na = (double *) malloc(sizeof(double) * dl);
    ng = (double *) malloc(sizeof(double) * l);

    for (i = 0; i < dl; i++)
        na[i] = a[i];
    for (i = 0; i < l; i++)
        ng[i] = g[i];

    gsl_matrix_view nA = gsl_matrix_view_array(na, l, l);
    gsl_vector_view nG = gsl_vector_view_array(ng, l);
    gsl_vector_view nH = gsl_vector_view_array(h, l);
    gsl_permutation *p = gsl_permutation_alloc(l);

    gsl_matrix_add(&nA.matrix, muImat);
    gsl_vector_scale(&nG.vector, -1);
    gsl_linalg_LU_decomp(&nA.matrix, p, &s);
    gsl_linalg_LU_solve(&nA.matrix, p, &nG.vector, &nH.vector);

    free((void *) na);
    free((void *) ng);
    gsl_permutation_free(p);
}

double get_rho(gsl_vector *fvect, 
               gsl_vector *newfvect,
               struct cfdata *data,
               double mu, double *h, double *g)
{
    int i;
    double dF, dF1, dF2, dL;
    double *c;

    c = (double *) malloc(sizeof(double) * data->varlen);

    gsl_vector_view new_h = gsl_vector_view_array(h, data->varlen);

    gsl_blas_ddot(fvect, fvect, &dF1);
    gsl_blas_ddot(newfvect, newfvect, &dF2);

    dF = 0.5 * dF1 - 0.5 * dF2;

    for (i = 0; i < data->varlen; i++)
        c[i] = mu * h[i] - g[i];

    gsl_vector_view c_v = gsl_vector_view_array(c, data->varlen);
    gsl_blas_ddot(&new_h.vector, &c_v.vector, &dL);
    dL *= 0.5;

    free((void *) c);

    return dF/dL;
}

void levenberg_marquardt(struct cfdata *data, double *var)
{
    int i, k, v;
    double mu, rho, normH;
    double *new_var;
    double *j;
    double *a;
    double *g;
    double *h;

    new_var = (double *) malloc(sizeof(double)*data->varlen);
    j = (double *) malloc(sizeof(double)*data->datalen * data->varlen);
    a = (double *) malloc(sizeof(double)*data->varlen * data->varlen);
    g = (double *) malloc(sizeof(double)*data->varlen);
    h = (double *) malloc(sizeof(double)*data->varlen);

    k = 0;
    v = 2;

    gsl_vector *f = gsl_vector_alloc(data->datalen);
    gsl_vector *new_f = gsl_vector_alloc(data->datalen);
    gsl_matrix *muI = gsl_matrix_alloc(data->varlen, data->varlen);
    
    gsl_matrix_view JT = gsl_matrix_view_array(j, data->varlen, data->datalen);
    gsl_matrix_view A = gsl_matrix_view_array(a, data->varlen, data->varlen);

    gsl_vector_view G = gsl_vector_view_array(g, data->varlen);
    gsl_vector_view H = gsl_vector_view_array(h, data->varlen);

    get_f(f, data, var);
    get_jac(j, data, var);

    // Compute A = J^T dot J
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,
                   1.0, &JT.matrix, &JT.matrix,
                   0.0, &A.matrix);

    // compute g = J^T dot f
    gsl_blas_dgemv(CblasNoTrans, 1.0, &JT.matrix, f, 0.0, &G.vector);

    if (!strcmp(data->model, "ic50"))
        mu = 1e-3*(GSL_MAX_DBL(a[0], GSL_MAX_DBL(a[4], a[8])));
    else
        mu = 1e-3*(GSL_MAX_DBL(a[0], a[3]));

    while ((fabs(g[gsl_blas_idamax(&G.vector)]) >= NORM_G) && (k < MAX_ITER)) {
        k++;
        gsl_matrix_set_identity(muI);
        gsl_matrix_scale(muI, mu);

        // solve (A + muI)h = -g
        solve_h(a, muI, data, g, h);

        for (i = 0; i < data->varlen; i++) {
            new_var[i] = var[i] + h[i];
        }
        normH = gsl_blas_dnrm2(&H.vector);

        if (gsl_blas_dnrm2(&H.vector) <= NORM_H) {
            printf("var converged\n");
            k = MAX_ITER;
        }

        get_f(new_f, data, new_var);

        rho = get_rho(f, new_f, data, mu, h, g);

        if (rho > 0) {
            for (i = 0; i < data->varlen; i++)
                var[i] = new_var[i];
            gsl_vector_swap(f, new_f);
            get_jac(j, data, var);
            // Recalculate A
            gsl_blas_dgemm(CblasNoTrans,
                           CblasTrans, 1.0,
                           &JT.matrix,
                           &JT.matrix, 0.0,
                           &A.matrix);
            gsl_blas_dgemv(CblasNoTrans, 1.0, &JT.matrix, f, 0.0, &G.vector);
            mu *= GSL_MAX_DBL((1.0/3.0), (1 - pow((2*rho - 1), 3.0)));
            v = 2;
        } else {
            mu *= v;
            v *= 2;

        }
        printf("iter %3d: var =", k);
        for (i = 0; i < data->varlen; i++)
            printf(" %8.5f", var[i]);
        printf(", |f(x)| = %10.5g\n", gsl_blas_dnrm2(f));
    }

    printf("\nFitted Parameters:\n");
    if (!strcmp(data->model, "ic50")) {
        printf("\tdelta:\t%8.5f\n"
               "\t ic50:\t%8.5f\n"
               "\t Hill:\t%8.5f\n", var[0], var[1], var[2]);
    } else
        printf("\tVmax:\t%8.5f\n"
               "\t  Km:\t%8.5f\n", var[0], var[1]);

    free((void *) new_var);
    free((void *) j);
    free((void *) a);
    free((void *) g);
    free((void *) h);

    gsl_vector_free(f);
    gsl_vector_free(new_f);
    gsl_matrix_free(muI);
}