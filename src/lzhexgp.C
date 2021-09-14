/** \file
 * Save the projection of an lzH in hex grid format for gnuplot
 *
 * In root, use .L or .x to load or run this.
 */
#include <stdio.h>
#include <stdlib.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2I.h>

#include "hexlib.c"

/** Save the projection of an lzH in hex grid format for gnuplot.
 * @param[in] lzH
 * @param[in] fname output data file name.
 * @param[in] p pixel pitch
 */
Int_t lzhexgp(TH2I *lzH, const char *fname, Double_t p=1.0)
{
    FILE *fp;
    Int_t i, j, l, q, r, NX;
    TAxis *XAxis;
    TH1D *lH;
    Double_t x,y,v,err, x0, y0, d;

    const Double_t SQRT3=1.7320508075688772;
    Double_t hexCoords[][2] = {
        {1.0, 0.0}, {0.5, SQRT3/2.0}, {-0.5, SQRT3/2.0},
        {-1.0, 0.0}, {-0.5, -SQRT3/2.0}, {0.5, -SQRT3/2.0}, {1.0, 0.0}};

    if((fp = fopen(fname, "w"))==NULL) {
        perror(fname);
        return 0;
    }
    d = p / SQRT3;
    
    lH = lzH->ProjectionX();
    NX = lH->GetNbinsX();
    XAxis = lH->GetXaxis();
    for(i=1; i<=NX; i++) {
        x   = XAxis->GetBinLowEdge(i);
        v   = lH->GetBinContent(i);
        err = lH->GetBinError(i);
        l = i - 1;
        hex_l2qr(l, &q, &r);
        hex_qr2xy(p, q, r, &x0, &y0);
        for(j=0; j<7; j++) {
                x = x0 + hexCoords[j][0] * d;
                y = y0 + hexCoords[j][1] * d;
                fprintf(fp, "%5d %5d %5d %10g %10g %10g %10g %10g %10g\n",
                        l, q, r, x0, y0, x, y, v, err);
        }
        fprintf(fp, "\n\n");
    }
    fclose(fp);
    return 1;
}
