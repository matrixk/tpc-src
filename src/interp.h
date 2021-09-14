/** \file
 * Utilities for 1D and 2D interpolation.
 */
#ifndef __INTERP_H__
#define __INTERP_H__

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)<(y)?(y):(x))
#endif

#ifndef Int_t
typedef int Int_t; /**< could use ssize_t if desired. */
#endif

/** Structure holding data for 1D interpolation. */
typedef struct
{
    Int_t n, mm, jsav, cor, dj;
    double *xx, *yy;
} interp1_t;

/** Initialize for 1D interpolation.
 * @param n Size of input arrays x and v.
 * @param m `Order' of interpolation -- number of points each interpolation would use (1D case).
 * @return Pointer to an allocated interp1_t structure.
 */
interp1_t *interp1_init(double *x, double *v, Int_t n, Int_t m);
/** Destruct a interp1_t structure.
 * @param hdl Pointer to an interp1_t structure.
 */
int interp1_close(interp1_t *hdl);

/** Find the index jlo, where x is properly centered in [jlo..jlo+m-1] */
Int_t interp1_jlo(interp1_t *hdl, double x);
/** Locate jlo using bisection.  Internal use. */
Int_t interp1_locate(interp1_t *hdl, double x);
/** Hunt up and down for jlo of a correlated value.  Internal use. */
Int_t interp1_hunt(interp1_t *hdl, double x);

/** Piecewise linear interpolation.
 * @param x
 * @return interpolated value.
 */
double interp1_linear(interp1_t *hdl, double x);

/** Structure holding data for 2D interpolation.
 * x and y are values of (irregular) grid points, where v contains the
 * value at (x, y).  v is accessed as v[ y * nx + x ]
 */
typedef struct
{
    Int_t nx, ny;
    double *x, *y, *v;
    interp1_t *hxterp, *hyterp;
} interp2_t;

/** Initialize for 2D interpolation.
 * @param x x location of (irregular) grid points.
 * @param y y location of (irregular) grid points.
 * @param v value at (x, y).  v is accessed as v[ y * nx + x ].
 * @return Pointer to an allocated interp2_t structure.
 */
interp2_t *interp2_init(double *x, double *y, double *v, Int_t nx, Int_t ny);
/** Destruct a interp2_t structure.
 * @param hdl Pointer to an interp2_t structure.
 */
int interp2_close(interp2_t *hdl);
/** Bilinear interpolation. */
double interp2_bilinear(interp2_t *hdl, double x, double y);

#endif /* __INTERP_H__ */
