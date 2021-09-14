/** \file
 * Generate pattern for Neural Network
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <getopt.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH2I.h>
#include <TH3I.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TUnuran.h>
#include <TUnuranMultiContDist.h>
#include <TMath.h>

#include <png.h>

#if defined(__MAKECINT__) || defined(__CINT__)
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<vector<int> >+;
#pragma link C++ class vector<vector<long> >+;
#endif

#ifndef MIN
#define MIN(a,b) ((a)>(b)?(b):(a))
#endif

/** Cholesky decomposition.  Adapted from NR3.
 * @param[in] a input array, must be square, symmetric and positive-definite, with n*n elements.
 * @param[in] n dimension of array a and L.
 * @param[out] L decomposition result, a lower triangular matrix.
 *               If (*L)=NULL is supplied, *L is allocated.
 * @return 1 when success, 0 if a is not positive-definite.
 */
int cholesky_decomp(const double *a, ssize_t n, double **L)
{
    ssize_t i, j, k;
    double sum, *el;
    if(*L == NULL) {
        (*L) = (double*)calloc(n*n, sizeof(double));
        if(*L == NULL) {
            fprintf(stderr, "%s(): *L allocation failure.\n", __FUNCTION__);
            return 0;
        }
    }
    el = *L;
    for(i=0; i<n*n; i++) {el[i] = a[i];}
    for(i=0; i<n; i++) {
        for(j=i; j<n; j++) {
            for(sum=el[i*n+j],k=i-1; k>=0; k--) sum -= el[i*n+k] * el[j*n+k];
            if(i == j) {
                if(sum <= 0.0) {
                    fprintf(stderr, "%s(): a is not positive-definite.\n", __FUNCTION__);
                    return 0;
                }
                el[i*n+i] = sqrt(sum);
            } else {
                el[j*n+i] = sum/el[i*n+i];
            }
        }
    }
    for(i=0; i<n; i++)
        for(j=0; j<i; j++)
            el[j*n+i] = 0.0;
    return 1;
}
/** n-dimmensional Gaussian (multivariate normal) distribution.
 * @param[in] n dimension.
 * @param[in] L Cholesky decomposed covariance matrix.
 * @param[in] mean vector of mean values.
 * @param[out] pt n-dimmensional Gaussian deviate output.  pt has to be pre-allocated.
 */
double *rand_gaussnd(TRandom *tr, ssize_t n, const double *L, const double *mean, double *pt)
{
#ifndef RAND_GAUSSND_NMAX
#define RAND_GAUSSND_NMAX 1024
#endif /* for efficiency */
    double spt[RAND_GAUSSND_NMAX];

    ssize_t i, j;
    double u, v, x, y, q;
    for(i=0; i<n; i++) { /* fill vector spt of independent normal deviates. */
        do{
            u = tr->Uniform();
            v = 1.7156*(tr->Uniform()-0.5);
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = x*x + y*(0.19600*y-0.25472*x);
        } while (q > 0.27597 && (q > 0.27846 || v*v > -4.*log(u)*u*u));
        spt[i] = v/u;
    }
    for(i=0; i<n; i++) {
        pt[i] = 0.0;
        for(j=0; j<=i; j++) pt[i] += L[i*n+j] * spt[j];
        pt[i] += mean[i];
    }
    return pt;
#undef RAND_GAUSSND_NMAX
}

typedef std::vector< std::vector<float> > vecvec;
// for xenon at 10bar
double Fano=0.14;
double Wi=24.8; // eV
double Dt=1.0; // transverse diffusion (sigma) [mm]
double Dl=2.0; // longitudinal diffusion (sigma) [mm]
int nbin=400, zbin=400;  // binning
double bmin=-200.0, bmax=200.0;
int pngscale=1;

int write_png(const char *fname, const char *title, const TH2I *xy, const TH2I *yz, const TH2I *zx)
{
    FILE *fp = NULL;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_bytep row = NULL;
    uint16_t *row16;
    int ret=1;

    ssize_t i, j;

    // Open png file for writing (binary mode)
    if((fp = fopen(fname, "wb")) == NULL) {
        fprintf(stderr, "Could not open file %s for writing\n", fname);
        ret = 0;
        goto out;
    }
    // Initialize png write structure
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if(png_ptr == NULL) {
        fprintf(stderr, "Could not allocate png write struct\n");
        ret = 0;
        goto out;
    }
    // Initialize png info structure
    info_ptr = png_create_info_struct(png_ptr);
    if(info_ptr == NULL) {
        fprintf(stderr, "Could not allocate info struct\n");
        ret = 0;
        goto out;
    }
    png_init_io(png_ptr, fp);
    // Write png header (16 bit color depth)
    png_set_IHDR(png_ptr, info_ptr, nbin, nbin,
                 16, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    // Set title
    if(title != NULL) {
        png_text title_text;
        title_text.compression = PNG_TEXT_COMPRESSION_NONE;
        title_text.key = (png_charp)"Title";
        title_text.text = (png_charp)title;
        png_set_text(png_ptr, info_ptr, &title_text, 1);
    }
    png_write_info(png_ptr, info_ptr);

    row = (png_bytep)malloc(sizeof(uint16_t)*3*nbin);
    row16 = (uint16_t*)row;
    png_set_swap(png_ptr); // for 16-bit endian

    for(i=0; i<nbin; i++) {
        for(j=0; j<nbin; j++) {
                 row16[j*3+0] = xy->GetBinContent(xy->GetBin(j+1, i+1)) * pngscale;
            if(i>=zbin) row16[j*3+1] = 0;
            else row16[j*3+1] = yz->GetBinContent(yz->GetBin(j+1, i+1)) * pngscale;
            if(j>=zbin) row16[j*3+2] = 0;
            else row16[j*3+2] = zx->GetBinContent(zx->GetBin(j+1, i+1)) * pngscale;
        }
        png_write_row(png_ptr, row);
    }
    // End write
    png_write_end(png_ptr, NULL);
out:
    if(fp != NULL) fclose(fp);
    if(info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if(png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    if(row != NULL) free(row);

    return ret;
}


int IoniImage(TTree *t1, const char *ofdir, int opng=-1)
{
    ssize_t i, j, k;
    char fname[PATH_MAX]={0}, ftmp[PATH_MAX], ftmp1[PATH_MAX];
    
    TFile *tfp=0;
    TTree *tp1=0;

    std::vector<int> *parentId=0;
    std::vector<double> *xp=0, *yp=0, *zp=0, *ed=0; // must be initialized to NULL
    int nIonTot;
    TH2I *xyH, *yzH, *zxH;
    int nIon;
    double mIon, sIon;
    double cov[9] = {Dt*Dt,   0.0,   0.0,
                       0.0, Dt*Dt,   0.0,
                       0.0,   0.0, Dl*Dl};
    double L[9], *Lp, sxyz[3], mean[3]={0.0, 0.0, 0.0};
    TRandom3 *tr = new TRandom3;
    Lp = L;
    cholesky_decomp(cov, 3, &Lp);

    xyH = new TH2I("xyH", "xy", nbin, bmin, bmax, nbin, bmin, bmax);
    yzH = new TH2I("yzH", "yz", nbin, bmin, bmax, zbin, bmin, bmax);
    zxH = new TH2I("zxH", "zx", zbin, bmin, bmax, nbin, bmin, bmax);
    
    t1->SetBranchAddress("parentid", &parentId);
    t1->SetBranchAddress("ed", &ed);
    t1->SetBranchAddress("xp", &xp);
    t1->SetBranchAddress("yp", &yp);
    t1->SetBranchAddress("zp", &zp);    

    if(opng<0) {
        tfp = new TFile(ofdir, "RECREATE");
        tfp->cd();
        tp1 = new TTree("p1", "IoniImage");
        tp1->Branch("nIonTot", &nIonTot, "nIonTot/I");
        tp1->Branch("xyH", "TH2I", &xyH);
        tp1->Branch("yzH", "TH2I", &yzH);
        tp1->Branch("zxH", "TH2I", &zxH);
    } else {
        strncpy(fname, ofdir, PATH_MAX-1); fname[PATH_MAX-1]='\0';
        i = strnlen(fname, PATH_MAX);
        if(fname[i-1]!='/') {
            fname[i]='/';
            fname[i+1] = '\0';
        }
        strncpy(ftmp, fname, PATH_MAX-1); ftmp[PATH_MAX-1]='\0';
        // check and create output dir
        struct stat sb;
        if(stat(fname, &sb)) { // failed
            if(errno == ENOENT) { // doesn't exist, make one
                if(mkdir(fname, 0755)) {
                    fprintf(stderr, "Unable to create directory %s\n", fname);
                    perror(fname);
                    return 0;
                }
            } else {
                fprintf(stderr, "Unable to access directory %s\n", fname);
                perror(fname);
                return 0;
            }
        } else {
            
        }
    }

    for(i=0; i<t1->GetEntries(); i++) {

        t1->GetEntry(i);
        xyH->Reset();
        yzH->Reset();
        zxH->Reset();

        //add randomness to the origin
        mean[0] = tr->Gaus(0, 0.1*bmax);
        mean[1] = tr->Gaus(0, 0.1*bmax);
        mean[2] = tr->Gaus(0, 0.1*bmax);
        
        nIonTot = 0;
        for(j=0; j<(ssize_t)ed->size(); j++) {
            mIon = (*ed)[j] * 1000.0 / Wi;
            sIon = std::sqrt(Fano * mIon);
            nIon = (int)tr->Gaus(mIon, sIon);
            if(nIon<0) nIon = 0;
            nIonTot += nIon;
            if(nIon == 0) continue;
            //handle diffusion
            for(k=0; k<nIon; k++) {
                rand_gaussnd(tr, 3, L, mean, sxyz);
                xyH->Fill((*xp)[j] + sxyz[0], (*yp)[j] + sxyz[1]);
                yzH->Fill((*yp)[j] + sxyz[1], (*zp)[j] + sxyz[2]);
                zxH->Fill((*zp)[j] + sxyz[2], (*xp)[j] + sxyz[0]);
            }
        }
        if(opng<0) {
            tp1->Fill();
        } else {
            sprintf(ftmp1, "p%07zd.png", i);
            strcpy(fname, ftmp);
            strcat(fname, ftmp1);
            write_png(fname, ftmp1, xyH, yzH, zxH);
        }
    }

    if(opng<0) {
        tp1->Write();
        tfp->CurrentFile()->Close();
        delete tfp;
    } else {
    }
    
    delete xyH;
    delete yzH;
    delete zxH;
    delete tr;

    return i;
}

#ifndef __CINT__

void print_usage(void)
{
    printf("Usage:\n");
    printf("      -n nbin[400] -z zbin[400] -i bmin[-200.0] -a bmax[200.0] binning[mm]\n");
    printf("      -w W[24.8] w-value[eV]\n");
    printf("      -f Fano[0.14] Fano factor\n");
    printf("      -p opng[-1] >=0 : output png\n");
    printf("      -s pngscale[1] multiply each RGB channel by pngscale\n");
    printf("      -t Dt[1.0] -l Dl[2.0] Transverse and longitudinal diffusion [mm]\n");
    printf("      -o outdir/\n");
    printf("      input.root\n");
}

/** Main entry if this file is compiled outside of root */
int main(int argc, char **argv)
{
    int optC = 0;
    int opng = -1;
    
    char *sRootFname;
    std::string ofdir, ofname;
    std::stringstream ss;
    

    TFile *tfs;
    TTree *t1s;

    // parse switches
    while((optC = getopt(argc, argv, "a:f:i:l:n:o:p:s:t:w:z:")) != -1)
    {
        switch(optC)
        {
        case 'a':
            bmax = atof(optarg);
            break;
        case 'f':
            Fano = atof(optarg);
            break;
        case 'i':
            bmin = atof(optarg);
            break;            
        case 'l':
            Dl = atof(optarg);
            break;
        case 'n':
            nbin = atoi(optarg);
            break;
        case 'o':
            ofdir = optarg;
            break;
        case 'p':
            opng = atoi(optarg);
            break;
        case 's':
            pngscale = atoi(optarg);
            break;
        case 't':
            Dt = atof(optarg);
            break;
        case 'w':
            Wi = atof(optarg);
            break;
        case 'z':
            zbin = atoi(optarg);
            break;
        default:
            print_usage();
            return EXIT_FAILURE;
            break;
        }
    }
    argc -= optind;
    argv += optind;
    if(argc<1 || argc>=2) {
        print_usage();
        return EXIT_FAILURE;
    }

    // signal
    sRootFname = argv[0];
    tfs = new TFile(sRootFname, "READ");
    if(tfs->IsZombie()) {
        std::cerr << "Error opening input file " << sRootFname << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "******************** Signal ************************" << std::endl;
    tfs->ls();
    t1s = (TTree*)(tfs->Get("t1"));
    t1s->ls();
    std::cout << "t1 has " << t1s->GetEntries() << " entries." << std::endl;

    IoniImage(t1s, ofdir.c_str(), opng);

    tfs->Close();
    delete tfs;

    return EXIT_SUCCESS;
}
#endif
