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
    printf("\nUsage: %s -f <filename> -m <model> [-x <var1> -y <var2> -z "
           "<var2>]\n"
           "\nFit experimental data to either an ic50 or Michaelis-Menten " 
           "model\n\n"
           "INPUT:\nRequired. Provide the filepath to the input data.\n"
           "\t-f FILENAME\tFile must contain tab delimited x and y data\n"
           "MODEL:\nRequired. Enter model name.\n"
           "\t-m MODEL\tSupported models are \"ic50\" and \"mm\"\n"
           "INITIAL GUESS:\nOptional. Provide values for fit parameters. "
           "Default values are 1.0.\n\t-x VALUE\tValue for max-min or Vmax.\n"
           "\t-y VALUE\tValue for ic50 or Km.\n\t-z VALUE\tValue for Hill "
           "coefficient.\n\nOUTPUT:\nUpdated values for the variables are "
           "displayed for each iteration,\nalong with |f(x)|.\n\n", argv[0]);
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
        switch (data->model) {
            case ic50:
                data->m[i] = (1 - (var[0]/(1 + pow((var[1]/data->x[i]), var[2]))));
                break;
            case mm:
                data->m[i] = ((var[0] * data->x[i]) / (var[1] + data->x[i]));
                break;
            default:
                exit(2);
        }
    }
}

void derivatives(struct cfdata *data, double *var)
{
    int i;
    for (i = 0; i < data->datalen; i++) {
        switch (data->model) {
            case ic50:
                data->d1[i] = (-(1/(1 + pow((var[1]/data->x[i]), var[2]))));
                data->d2[i] = ((var[0] * var[2] * pow((var[1]/data->x[i]), 
                                                      (var[2] - 1))) / 
                               (data->x[i] * pow((1 + (pow((var[1]/data->x[i]), 
                                                           var[2]))), 2)));
                data->d3[i] = ((var[0] * pow((var[1]/data->x[i]), var[2]) * 
                                log((var[1]/data->x[i]))) / 
                               (pow((1 + pow((var[1]/data->x[i]), var[2])), 2)));
                break;
            case mm:
                data->d1[i] = (data->x[i] / (var[1] + data->x[i]));
                data->d2[i] = (-(var[0] * data->x[i]) / pow((var[1] + data->x[i]), 
                                                            2.0));
                break;
            default:
                exit(3);
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
        switch (data->model) {
            case ic50:
                jac[i] = -data->d1[i];
                jac[j] = -data->d2[i];
                jac[k] = -data->d3[i];
                break;
            case mm:
                jac[i] = -data->d1[i];
                jac[j] = -data->d2[i];
                break;
            default:
                exit(4);
        }
    }
}

void solve_h(double *a, gsl_matrix *muImat, 
             struct cfdata *data, double *g, double *h)
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


double get_mu(int vlen, double *a)
{
    int i;
    int alen = vlen*vlen;
    double max = a[0];
    
    for (i=0; i < alen; i += (vlen + 1))
        max = (max < a[i]) ? a[i] : max;
    
    return 1e-3*max;
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
    
    mu = get_mu(data->varlen, a);

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
    switch (data->model) {
        case ic50:
            printf("\tdelta:\t%8.5f\n"
                   "\t ic50:\t%8.5f\n"
                   "\t Hill:\t%8.5f\n\n", var[0], var[1], var[2]);
            break;
        case mm:
            printf("\tVmax:\t%8.5f\n"
                   "\t  Km:\t%8.5f\n\n", var[0], var[1]);
            break;
        default:
            exit(5);
            
    }

    free((void *) new_var);
    free((void *) j);
    free((void *) a);
    free((void *) g);
    free((void *) h);

    gsl_vector_free(f);
    gsl_vector_free(new_f);
    gsl_matrix_free(muI);
}

void output(char *filename, struct cfdata *data, double *var)
{
    FILE *fdata, *fscript;
    int i;
    size_t len = strlen(filename);
    char data_tag[] = "_input.dat";
    char script_tag[] = "_script.plt";
    char plot_tag[] = "_plot.png";
    char *outdata;
    char *outscript;
    char *outplot;
    double minx, maxx, miny, maxy;
    
    outdata = (char *) malloc(sizeof(char) * (len + 11));
    outscript = (char *) malloc(sizeof(char) * (len + 12));
    outplot = (char *) malloc(sizeof(char) * (len + 10));
    
    for (i = 0; i < (len - 4); i++) {
        outdata[i] = outscript[i] = outplot[i] = filename[i];
    }
    outdata[i] = outscript[i] = outplot[i] = '\0';
    
    strncat(outdata, data_tag, strlen(data_tag));
    strncat(outscript, script_tag, strlen(script_tag));
    strncat(outplot, plot_tag, strlen(plot_tag));
    
    minx = maxx = data->x[0];
    miny = maxy = data->y[0];
    for (i = 0; i < data->datalen; i++) {
        double xi = data->x[i];
        double yi = data->y[i];
        minx = ((minx < xi) ? minx : xi);
        maxx = ((maxx > xi) ? maxx : xi);
        miny = ((miny < yi) ? miny : yi);
        maxy = ((maxy > yi) ? maxy : yi);
    }
    
    miny -= 0.2*maxy;
    maxy += 0.2*maxy;
    
    
    fdata = fopen(outdata, "w");
    for (i = 0; i < data->datalen; i++)
        fprintf(fdata, "%g\t%g\n", data->x[i], data->y[i]);
    fclose(fdata);
    
    fscript = fopen(outscript, "w");
    
    switch (data->model) {
        case ic50:
            fprintf(fscript, 
                    "set terminal png\n"
                    "set output \"%s\"\n"
                    "f(x) = 1 - (%f/(1 + (%f/x)**%f))\n"
                    "set xrange [%f:%f]\n"
                    "set yrange [%f:%f]\n"
                    "set log x\n"
                    "plot \"%s\" using 1:2 with points 4,"
                    " f(x) with lines 22\n", 
                    outplot, var[0], var[1], var[2], 
                    minx, maxx, miny, maxy, outdata);
            break;
        case mm:
            fprintf(fscript, 
                    "set terminal png\n"
                    "set output \"%s\"\n"
                    "f(x) = (%f * x) / (%f + x)\n"
                    "set xrange [%f:%f]\n"
                    "set yrange [%f:%f]\n"
                    "plot \"%s\" using 1:2 with points 4,"
                    " f(x) with lines 22\n", 
                    outplot, var[0], var[1], 
                    minx, maxx, miny, maxy, outdata);
            break;
        default:
            exit(6);
    }

    fclose(fscript);
    
    free((void *) outdata);
    free((void *) outscript);
    free((void *) outplot);
}
