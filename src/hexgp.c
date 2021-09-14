/** \file
 * Generate data for plotting hexagonal grids in gnuplot.
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "hexlib.h"

#define SQRT3 1.7320508075688772

typedef struct param 
{
    int N; // total number of pixels
    int R; // radius of the grid
    double p; // pitch
    char op;
} param_t;

param_t param_default = {
    .N = 3 * (1+3) * 3 + 1,
    .R = 3,
    .p = 10.0,
    .op = 'g'
};

void print_usage(const param_t *pm)
{
    printf("Usage:\n");
    printf("      operation: [-%c] -g: grid, -c: central and (q,r) label\n", pm->op);
    printf("      -p pixel pitch[%g] [mm]\n", pm->p);
    printf("      -R radius of the grid[%d]\n", pm->R);
}

int main(int argc, char **argv)
{
    int optC = 0;
    param_t pm;

    ssize_t i;
    int l, q, r;
    double d, x0, y0, x, y;
    double hexCoords[][2] = {
        {1.0, 0.0}, {0.5, SQRT3/2.0}, {-0.5, SQRT3/2.0},
        {-1.0, 0.0}, {-0.5, -SQRT3/2.0}, {0.5, -SQRT3/2.0}, {1.0, 0.0}};

    memcpy(&pm, &param_default, sizeof(pm));
    // parse switches
    while((optC = getopt(argc, argv, "cghp:R:")) != -1)
    {
        switch(optC)
        {
        case 'c':
            pm.op = optC;
            break;
        case 'g':
            pm.op = optC;
            break;
        case 'h':
            print_usage(&pm);
            return EXIT_FAILURE;
            break;
        case 'p':
            pm.p = atof(optarg);
            break;
        case 'R':
            pm.R = atoi(optarg);
            break;
        default:
            print_usage(&pm);
            return EXIT_FAILURE;
            break;
        }
    }
    argc -= optind;
    argv += optind;

    pm.N = 3*(1+pm.R)*pm.R + 1; // compute N from R
    d = pm.p / sqrt(3.0);

    for(l=0; l<pm.N; l++) {
        hex_l2qr(l, &q, &r);
        hex_qr2xy(pm.p, q, r, &x0, &y0);
        if(pm.op == 'g') {
            for(i=0; i<7; i++) {
                x = x0 + hexCoords[i][0] * d;
                y = y0 + hexCoords[i][1] * d;
                printf("%5d %5d %5d %10g %10g %10g %10g\n",
                       l, q, r, x0, y0, x, y);
            }
            printf("\n\n");
        }
        if(pm.op == 'c') {
            printf("%5d %5d %5d %10g %10g\n",
                   l, q, r, x0, y0);
        }
    }
    return EXIT_SUCCESS;
}
