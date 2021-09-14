/** \file
 * Save a TH2 histogram's contents into a text .dat file for plotting in gnuplot.
 *
 * In root, use .L or .x to load or run this.
 */
#include <stdio.h>
#include <stdlib.h>
#include <TH2.h>

int h22GpDat(TH2 *h2, const char *fname)
{
    FILE *fp;
    Int_t i, j, NX, NY;
    TAxis *XAxis, *YAxis;
    Double_t x,y,v,err;

    if((fp = fopen(fname, "w"))==NULL) {
        perror(fname);
        return 0;
    }
    NX = h2->GetNbinsX();
    XAxis = h2->GetXaxis();
    NY = h2->GetNbinsY();
    YAxis = h2->GetYaxis();
    for(j=1; j<=NY; j++) {
        y = YAxis->GetBinLowEdge(j);
        for(i=1; i<=NX; i++) {
            x = XAxis->GetBinLowEdge(i);
            v = h2->GetBinContent(i, j);
            err = h2->GetBinError(i, j);
            fprintf(fp, "%g %g %24.16e %24.16e %d %d\n", x, y, v, err, i, j);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 1;
}
