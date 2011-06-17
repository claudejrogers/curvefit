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
    FILE *fp;
    char *filename = NULL;
    struct cfdata values;
    double var[3];
    var[0] = var[1] = var[2] = 1.0;

    opterr = 0;

    while ((c = getopt(argc, argv, "f:m:x:y:z:h")) != -1)
        switch (c)
    {
        case 'f':
            filename = optarg;
            break;
        case 'm':
            if (!strcmp(optarg, "ic50"))
                values.model = ic50;
            else if (!strcmp(optarg, "mm"))
                values.model = mm;
            else {
                fprintf(stderr, "%s is not a supported model.\n", optarg);
                usage(argv);
                return 1;
            }
            break;
        case 'x':
            var[0] = (double)atof(optarg);
            break;
        case 'y':
            var[1] = (double)atof(optarg);
            break;
        case 'z':
            var[2] = (double)atof(optarg);
            break;
        case 'h':
            usage(argv);
            return 0;
        case '?':
            if (optopt == 'f' || optopt == 'm')
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
            else if ((optopt != 'x' || optopt != 'y' || optopt != 'z') && isprint(optopt))
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

    const unsigned int len = filelines(filename);

    fp = fopen(filename, "r");

    values.datalen = len;
    
    switch (values.model) {
        case ic50:
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

    xv = (double *) malloc(sizeof(double)*len);
    yv = (double *) malloc(sizeof(double)*len);
    mv = (double *) malloc(sizeof(double)*len);
    d1v = (double *) malloc(sizeof(double)*len);
    d2v = (double *) malloc(sizeof(double)*len);
    d3v = (double *) malloc(sizeof(double)*len);

    values.x = xv;
    values.y = yv;
    values.m = mv;
    values.d1 = d1v;
    values.d2 = d2v;
    values.d3 = d3v;

    float tempX, tempY;
    for (i = 0; fscanf(fp, "%f\t%f", &tempX, &tempY) != EOF; i++){
        values.x[i] = (double)tempX;
        values.y[i] = (double)tempY;
    }

    levenberg_marquardt(&values, var);

    fclose(fp);
    
    output(filename, &values, var);

    free((void *) xv);
    free((void *) yv);
    free((void *) mv);
    free((void *) d1v);
    free((void *) d2v);
    free((void *) d3v);

    return 0;
}
