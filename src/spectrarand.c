#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <fftw3.h>
#include "interp.h"
#include "spectrarand.h"

#ifndef LINE_MAX
#define LINE_MAX 4096
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

extern double rand0_1(void);

spectrarand_t *spectrarand_init_from_file(size_t n, double fs, const char *fname)
{
    ssize_t i;
    
    spectrarand_t *hdl=NULL;
    interp1_t *interp1hdl=NULL;

    size_t m;
    char line[LINE_MAX];
    double f, psd, df, *xx, *vv;
    FILE *fp;

    if((fp=fopen(fname, "r"))==NULL) {
        perror(fname);
        return NULL;
    }
    /* count number of lines with data */
    m = 0;
    while(fgets(line, LINE_MAX, fp)!=NULL) {
        if(line[0] == '#') continue; /* comment line */
        if(sscanf(line, "%lf%lf\n", &f, &psd) == 2)
            m++;
    }
    m++;
    xx = (double*)calloc(m, sizeof(double));
    vv = (double*)calloc(m, sizeof(double));

    rewind(fp);
    i = 1;
    while(fgets(line, LINE_MAX, fp)!=NULL) {
        if(line[0] == '#') continue; /* comment line */
        if(sscanf(line, "%lf%lf\n", &f, &psd) == 2) {
            xx[i] = f;
            vv[i] = psd;
            i++;
        }
    }
    fclose(fp);

    hdl = (spectrarand_t*)calloc(1, sizeof(spectrarand_t));
    hdl->as = (double*)calloc(n, sizeof(double));
    hdl->rs = (double*)calloc(n, sizeof(double));

    interp1hdl = interp1_init(xx, vv, m, 2);

    hdl->fs = fs;
    hdl->n = n;
    df = fs/(double)n;
    hdl->as[0] = 0.0; /* enforce 0 mean */
    for(i=1; i<n; i++) { /* note that due to Nyquist theorem, only i*df < fs/2 is useful */
        psd = interp1_linear(interp1hdl, i*df);
        hdl->as[i] = sqrt(psd * df);
    }

    hdl->fftwPlan = fftw_plan_r2r_1d(hdl->n, hdl->rs, hdl->rs, FFTW_HC2R, FFTW_ESTIMATE);
        
    interp1_close(interp1hdl);
    if(xx) free(xx);
    if(vv) free(vv);
    return hdl;
}

int spectrarand_close(spectrarand_t *hdl)
{
    if(hdl) {
        fftw_destroy_plan(hdl->fftwPlan);
        if(hdl->as)
            free(hdl->as);
        if(hdl->rs)
            free(hdl->rs);
        free(hdl);
    }
    fftw_cleanup();
    return 0;
}

double *spectrarand_rand(spectrarand_t *hdl)
{
    ssize_t i;
    double phi, re, im;

    hdl->rs[0] = hdl->as[0];

    for(i=1; i<=hdl->n/2; i++) {
        phi = rand0_1() * 2.0 * M_PI;
        re = cos(phi);
        im = sin(phi);

        hdl->rs[i] = hdl->as[i] * re;
        hdl->rs[hdl->n - i] = hdl->as[i] * im;
    }
    if(hdl->n % 2 == 0) /* When n is even, last one (real) is in the middle of the array. */
        hdl->rs[hdl->n/2] = hdl->as[hdl->n/2];
    
    fftw_execute(hdl->fftwPlan);
    return hdl->rs;
}

#ifdef SPECTRARAND_DEBUG_ENABLEMAIN
int main(int argc, char **argv)
{
    size_t n;
    ssize_t i;
    double fs;
    spectrarand_t *hdl;

    n = 65536;
    fs = 1.0e7;
    hdl = spectrarand_init_from_file(n, fs, argv[1]);
    for(i=0; i<n; i++) {
        printf("%24.16e %24.16e\n", fs/(double)n * i, hdl->as[i]);
    }
    printf("\n\n");
    spectrarand_rand(hdl);
    for(i=0; i<hdl->n; i++) {
        printf("%24.16e %24.16e\n", 1.0/fs*i, hdl->rs[i]);
    }
    printf("\n");
    for(i=0; i<hdl->n; i++) {
        printf("%24.16e %24.16e\n", 1.0/fs * hdl->n + 1.0/fs*i, hdl->rs[i]);
    }   
    spectrarand_close(hdl);
    return EXIT_SUCCESS;
}
#endif
