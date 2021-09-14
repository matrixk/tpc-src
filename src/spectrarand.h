/** \file 
 * Generate a random sequence that has a specified power spectral
 * density and uniformly distributed phase.
 */
#ifndef __SPECTRARAND_H__
#define __SPECTRARAND_H__

#include <fftw3.h>

typedef struct spectrarand_handle
{
    size_t n;    /**< number of sampling points */
    double fs;   /**< sampling frequency [Hz] */
    double *as ; /**< array holding amplitude spectra [V] */
    double *rs;  /**< random sequence */
    fftw_plan fftwPlan;
} spectrarand_t;

/** Initialize a spectrarand_t structure with a file containing a PSD.
 * The PSD is interpolated and sampled according to parameters n and
 * fs.
 * @param n number of sampling points.
 * @param fs sampling frequency [Hz].
 * @param fname PSD data file name.
 * @return allocated spectrarand_t *hdl.
 */
spectrarand_t *spectrarand_init_from_file(size_t n, double fs, const char *fname);

/** Release memory in spectrarand_t *hdl
 * @param hdl
 */
int spectrarand_close(spectrarand_t *hdl);

/** Generate random sequence.
 * @param hdl
 * @return hdl->rs
 */
double *spectrarand_rand(spectrarand_t *hdl);

#endif /* __SPECTRARAND_H__ */
