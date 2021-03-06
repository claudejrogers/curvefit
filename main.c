//
//  main.c
//  curvefit
//
//  Created by Claude Rogers on 6/8/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fit_functions.h"

int main (int argc, char **argv)
{
    int i, c;
    char *filename = NULL;
    FileData fd;
    CFData values;
    double var[4];
    var[0] = var[1] = var[2] = var[3] = 1.0;

    opterr = 0;

    while ((c = getopt(argc, argv, "f:m:a:b:c:d:h")) != -1)
        switch (c)
    {
        case 'f':
            filename = optarg;
            break;
        case 'm':
            if (!strcmp(optarg, "boltzmann"))
                values.model = boltzmann;
            else if (!strcmp(optarg, "expdecay"))
                values.model = expdecay;
            else if (!strcmp(optarg, "gaussian"))
                values.model = gaussian;
            else if (!strcmp(optarg, "hill"))
                values.model = hill;
            else if (!strcmp(optarg, "ic50"))
                values.model = ic50;
            else if (!strcmp(optarg, "mm"))
                values.model = mm;
            else if (!strcmp(optarg, "modsin"))
                values.model = modsin;
            else {
                fprintf(stderr, "%s is not a supported model.\n", optarg);
                usage(argv);
                return 1;
            }
            break;
        case 'a':
            var[0] = (double)atof(optarg);
            break;
        case 'b':
            var[1] = (double)atof(optarg);
            break;
        case 'c':
            var[2] = (double)atof(optarg);
            break;
        case 'd':
            var[3] = (double)atof(optarg);
            break;
        case 'h':
            usage(argv);
            return 0;
        case '?':
            if (optopt == 'f' || optopt == 'm')
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            else if ((optopt != 'x' ||
                      optopt != 'y' ||
                      optopt != 'z') && isprint(optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            usage(argv);
            return 1;
        default:
            abort();
    }
    if (!filename) {
        fprintf(stderr, "Missing required argument\n");
        usage(argv);
        return 1;
    }

    readfile(filename, &fd);

    values.datalen = fd.len;

    switch (values.model) {
        case boltzmann:
        case gaussian:
            values.varlen = 4;
            break;
        case expdecay:
        case hill:
        case ic50:
        case modsin:
            values.varlen = 3;
            break;
        case mm:
            values.varlen = 2;
            break;
        default:
            fprintf(stderr, "Missing required argument\n");
            usage(argv);
            exit(7);
            break;
    }

    double *xv;
    double *yv;
    double *mv;
    double *d1v;
    double *d2v;
    double *d3v;
    double *d4v;

    xv =  (double *) malloc(sizeof(double) * fd.len);
    yv =  (double *) malloc(sizeof(double) * fd.len);
    mv =  (double *) malloc(sizeof(double) * fd.len);
    d1v = (double *) malloc(sizeof(double) * fd.len);
    d2v = (double *) malloc(sizeof(double) * fd.len);
    d3v = (double *) malloc(sizeof(double) * fd.len);
    d4v = (double *) malloc(sizeof(double) * fd.len);

    values.x = xv;
    values.y = yv;
    values.m = mv;
    values.d1 = d1v;
    values.d2 = d2v;
    values.d3 = d3v;
    values.d4 = d4v;

    for (i = 0; i < fd.len; i++){
        values.x[i] = fd.xdata[i];
        values.y[i] = fd.ydata[i];
    }

    free((void *) fd.xdata);
    free((void *) fd.ydata);

    levenberg_marquardt(&values, var);

    output(filename, &values, var);

    free((void *) xv);
    free((void *) yv);
    free((void *) mv);
    free((void *) d1v);
    free((void *) d2v);
    free((void *) d3v);
    free((void *) d4v);

    return 0;
}
