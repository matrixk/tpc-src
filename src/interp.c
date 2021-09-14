#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "interp.h"

interp1_t *interp1_init(double *x, double *v, Int_t n, Int_t m)
/* n is the size of input arrays x and y.  m is the `order' of
 * interpolation -- number of points each interpolation would use (1D
 * case).
 */
{
    interp1_t *hdl;
    hdl = (interp1_t *)malloc(sizeof(interp1_t));
    hdl->xx = x;
    hdl->yy = v;
    hdl->n = n;
    hdl->mm = m;
    hdl->jsav = 0;
    hdl->cor = 0;
    hdl->dj = MIN(1, (Int_t)pow((double)n,0.25));
    return hdl;
}

int interp1_close(interp1_t *hdl)
{
    if(hdl)
        free(hdl);
    return 0;
}

/* Find the index jlo, where x is properly centered in [jlo..jlo+m-1] */
Int_t interp1_jlo(interp1_t *hdl, double x)
{
    Int_t jlo;
    jlo = hdl->cor ? interp1_hunt(hdl, x) : interp1_locate(hdl, x);
    return jlo;
}

Int_t interp1_locate(interp1_t *hdl, double x)
{
    Int_t ju, jm, jl;
    int ascnd; /* ascending? bool flag */

    if(hdl->n < 2 || hdl->mm < 2 || hdl->mm > hdl->n) {
        return -1;
    }
    ascnd = (hdl->xx[hdl->n - 1] >= hdl->xx[0]);
    jl = 0;
    ju = hdl->n - 1;
    while(ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if((x >= hdl->xx[jm]) == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    hdl->cor = abs(jl - hdl->jsav) > hdl->dj ? 0 : 1;
    hdl->jsav = jl;
    return MAX(0, MIN(hdl->n - hdl->mm, jl - ((hdl->mm - 2)>>1)));
}

Int_t interp1_hunt(interp1_t *hdl, double x)
{
    Int_t jl, jm, ju, inc=1;
    int ascnd;
    jl = hdl->jsav;
    if(hdl->n < 2 || hdl->mm < 2 || hdl->mm > hdl->n) {
        return -1;
    }
    ascnd = (hdl->xx[hdl->n - 1] >= hdl->xx[0]);

    if(jl < 0 || jl > hdl->n-1) {
        jl=0;
        ju = hdl->n - 1;
    } else {
        if((x >= hdl->xx[jl]) == ascnd) {
            for(;;) {
                ju = jl + inc;
                if(ju >= hdl->n - 1) { ju = hdl->n - 1; break;}
                else if((x < hdl->xx[ju]) == ascnd) break;
                else {
                    jl = ju;
                    inc += inc;
                }
            }
        } else {
            ju = jl;
            for(;;) {
                jl = jl - inc;
                if(jl <= 0) { jl = 0; break;}
                else if((x >= hdl->xx[jl]) == ascnd) break;
                else {
                    ju = jl;
                    inc += inc;
                }
            }
        }
    }
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1;
        if((x >= hdl->xx[jm]) == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    hdl->cor = abs(jl - hdl->jsav) > hdl->dj ? 0 : 1;
    hdl->jsav = jl;
    return MAX(0,MIN(hdl->n - hdl->mm, jl-((hdl->mm - 2)>>1)));
}

double interp1_linear(interp1_t *hdl, double x)
{
    Int_t j;

    j = interp1_jlo(hdl, x);
    if(hdl->xx[j] == hdl->xx[j+1]) {
        return hdl->yy[j];
    } else {
        return hdl->yy[j] + ((x - hdl->xx[j])/(hdl->xx[j+1] - hdl->xx[j]))
            * (hdl->yy[j+1] - hdl->yy[j]);
    }
}

interp2_t *interp2_init(double *x, double *y, double *v, Int_t nx, Int_t ny)
{
    interp2_t *hdl;

    hdl = (interp2_t *)malloc(sizeof(interp2_t));
    hdl->hxterp = interp1_init(x, x, nx, 2);
    hdl->hyterp = interp1_init(y, y, ny, 2);
    hdl->x = x;
    hdl->y = y;
    hdl->v = v;
    hdl->nx = nx;
    hdl->ny = ny;

    return hdl;
}

int interp2_close(interp2_t *hdl)
{
    if(hdl) {
        if(hdl->hxterp) interp1_close(hdl->hxterp);
        if(hdl->hyterp) interp1_close(hdl->hyterp);
        free(hdl);
    }
    return 0;
}

double interp2_bilinear(interp2_t *hdl, double x, double y)
{
    Int_t i, j;
    double yy, t, u;

    i = interp1_jlo(hdl->hxterp, x);
    j = interp1_jlo(hdl->hyterp, y);

    t = (x - hdl->hxterp->xx[i]) / (hdl->hxterp->xx[i+1] - hdl->hxterp->xx[i]);
    u = (y - hdl->hyterp->xx[j]) / (hdl->hyterp->xx[j+1] - hdl->hyterp->xx[j]);

    yy = (1.0 - t) * (1.0 - u) * hdl->v[j * hdl->nx + i]
        + t * (1.0 - u) * hdl->v[j * hdl->nx + i+1]
        + (1.0 - t) * u * hdl->v[(j+1) * hdl->nx + i]
        + t * u * hdl->v[(j+1) * hdl->nx + i + 1];
    return yy;
}

#ifdef INTERP_DEBUG_ENABLEMAIN
int main(int argc, char **argv) 
{
    double x[] = {-1.0, 0.0, 1.5, 3.0, 3.3, 3.8, 5.0, 6.0};
    double y[] = {0.0, 0.2, 0.5, 1.0, 1.2, 1.3, 4.5, 5.0};
    double v[64] = {3.0, 10.0, 1.0, 2.0,
                  -10.0, 2.0, 8.2, 3.3,
                  12.0, -3.0, 2.2, 12.5,
                  0.0, 10.0, -5.0, 0.0};
    double a, b;
    Int_t nx = 8, ny = 8, i, j;
    interp2_t *hdl;
    
    for(j=0; j<ny; j++) {
        for(i=0; i<nx; i++) {
            v[j * nx + i] = x[i] * y[j];
            printf("%g %g %g\n", x[i], y[j], v[j * nx + i]);
        }
        printf("\n");
    }
    printf("\n");

    hdl = interp2_init(x, y, v, nx, ny);
    for(a = -0.1; a<5.1; a+=0.1) {
        for(b = -1.1; b<6.1; b+= 0.1) {
            printf("%g %g %g\n", b, a, interp2_bilinear(hdl, b, a));
        }
        printf("\n");
    }

    interp2_close(hdl);
    
    return EXIT_SUCCESS;
}
#endif /* INTERP_DEBUG_ENABLEMAIN */
