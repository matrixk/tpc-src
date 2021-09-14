#ifndef _LIB_H		
#define _LIB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TH1D.h>


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define ALPH 1.5
#define NDMX 50
#define MXDIM 10
#define TINY 1.0e-30
#define NR_END 1
#define FREE_ARG char*

int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?  (imaxarg1) : (imaxarg2))

int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))


long idum;      /* for ranno */
int ndim;       /* for fxn */

void rebin(double rc, int nd, double r[], double xin[], double xi[]);
double ran2(long *idum);


double *vec(long,long);

double myvegas(double regn[], int ndim, double (*fxn)(double [], double));
void vegas(double regn[], int ndim, double (*fxn)(double [], double), int init,
	unsigned long ncall, int itmx, int nprn, double *tgral, double *sd,
	double *chi2a);
	
	

using std::string;
using std::vector;

using namespace std;


void lib_matter(char * in, char *out){

	int length=strlen(in);
	int i=0;
	int j=0;
	for(i=0;i<length;i++){
		if(in[length-i-1]=='/') {
		j=length-i;
		break;
		}
	}

	for(i=j;i<length;i++){
		if(in[i]=='_'){
		    out[i-j]=0;
		}else{
			out[i-j]=in[i];
		}

	}


}


void lib_particle(char * in, char *out){

	int length=strlen(in);
	int i=0;
	int j=0;
	for(i=0;i<length;i++){
		if(in[length-i-1]=='/') {
		j=length-i;
		break;
		}
	}

	for(i=j;i<length;i++){
		if(in[i]=='_') {
		j=i+1;
		break;
		}
	}	
	
	for(i=j;i<length;i++){
		if(in[i]=='_'){
		    out[i-j]=0;
		}else{
			out[i-j]=in[i];
		}

	}
}

void lib_energy(char * in, char *out){

	int length=strlen(in);
	int i=0;
	int j=0;
	for(i=0;i<length;i++){
		if(in[length-i-1]=='/') {
		j=length-i;
		break;
		}
	}

int k=0;
for(k=0;k<2;k++){	
	for(i=j;i<length;i++){
		if(in[i]=='_') {
		j=i+1;
		break;
		}
	}	
}
	j+=strlen("E");
	
	for(i=j;i<length;i++){
		if(in[i]=='_'){
		    out[i-j]=0;
		}else{
			out[i-j]=in[i];
		}

	}
}


class EventData
{
public:

    EventData();
    ~EventData();
    
public:

    void Clear();
    
public:

    int m_eventId;
    int m_nbHits;
    int m_nbSteps;

    double m_totalEnergyDeposited;
    vector<int> *m_trackId;
    vector<int> *m_parentId;
    vector<string> *m_particleType;
    vector<string> *m_parentType;
    vector<string> *m_creatorProcess;
    vector<string> *m_depositingProcess;
    vector<double> *m_xp;               // position x,y,z
    vector<double> *m_yp;
    vector<double> *m_zp;
    vector<double> *m_xm;
    vector<double> *m_ym;
    vector<double> *m_zm;
    vector<double> *m_energyDeposited;
    vector<double> *m_kineticEnergy;
    vector<double> *m_time;

    vector<string> *m_primaryParticleType;
    double m_primaryParticleEnergy;
    double m_primaryXm;
    double m_primaryYm;
    double m_primaryZm;
    double m_primaryX;
    double m_primaryY;
    double m_primaryZ;
};


EventData::EventData()
{
    m_eventId = 0;
    m_nbHits = 0;
    m_nbSteps = 0;
    
    m_totalEnergyDeposited = 0.0;
    m_trackId = new vector<int>;
    m_parentId = new vector<int>;
    m_particleType = new vector<string>;
    m_parentType = new vector<string>;
    m_creatorProcess = new vector<string>;
    m_depositingProcess = new vector<string>;
    m_xp = new vector<double>;
    m_yp = new vector<double>;
    m_zp = new vector<double>;
    m_xm = new vector<double>;
    m_ym = new vector<double>;
    m_zm = new vector<double>;
    m_energyDeposited = new vector<double>;
    m_kineticEnergy = new vector<double>;
    m_time = new vector<double>;

    m_primaryParticleType = new vector<string>;
    m_primaryParticleEnergy = 0.0;
    m_primaryXm = 0.0;
    m_primaryYm = 0.0;
    m_primaryZm = 0.0;
    m_primaryX = 0.0;
    m_primaryY = 0.0;
    m_primaryZ = 0.0;
}

EventData::~EventData()
{
    delete m_trackId;
    delete m_parentId;
    delete m_particleType;
    delete m_parentType;
    delete m_creatorProcess;
    delete m_depositingProcess;
    delete m_xp;
    delete m_yp;
    delete m_zp;
    delete m_xm;
    delete m_ym;
    delete m_zm;
    delete m_energyDeposited;
    delete m_kineticEnergy;
    delete m_time;
    delete m_primaryParticleType;
}

void
EventData::Clear()
{
    m_eventId = 0;
    m_nbHits = 0;
    m_nbSteps = 0;

    m_totalEnergyDeposited = 0.0;
    m_trackId->clear();
    m_parentId->clear();
    m_particleType->clear();
    m_parentType->clear();
    m_creatorProcess->clear();
    m_depositingProcess->clear();
    m_xp->clear();
    m_yp->clear();
    m_zp->clear();
    m_xm->clear();
    m_ym->clear();
    m_zm->clear();
    m_energyDeposited->clear();
    m_kineticEnergy->clear();
    m_time->clear();

    m_primaryParticleType->clear();
    m_primaryParticleEnergy = 0.0;
    m_primaryXm = 0.0;
    m_primaryYm = 0.0;
    m_primaryZm = 0.0;
    m_primaryX = 0.0;
    m_primaryY = 0.0;
    m_primaryZ = 0.0;
}

void set_tree(TTree *t1,EventData *m_eventData){
    t1->SetBranchAddress("eventid", &m_eventData->m_eventId);
    t1->SetBranchAddress("nbhits", &m_eventData->m_nbHits);
    t1->SetBranchAddress("nbsteps", &m_eventData->m_nbSteps);

    t1->SetBranchAddress("etot", &m_eventData->m_totalEnergyDeposited);
    t1->SetBranchAddress("trackid",  &m_eventData->m_trackId);
    t1->SetBranchAddress("parentid",  &m_eventData->m_parentId);
    t1->SetBranchAddress("type",  &m_eventData->m_particleType);
    t1->SetBranchAddress("parenttype",  &m_eventData->m_parentType);
    t1->SetBranchAddress("creatproc",  &m_eventData->m_creatorProcess);
    t1->SetBranchAddress("edproc",  &m_eventData->m_depositingProcess);
    t1->SetBranchAddress("xp",  &m_eventData->m_xp);
    t1->SetBranchAddress("yp",  &m_eventData->m_yp);
    t1->SetBranchAddress("zp",  &m_eventData->m_zp);
    t1->SetBranchAddress("xm",  &m_eventData->m_xm);
    t1->SetBranchAddress("ym",  &m_eventData->m_ym);
    t1->SetBranchAddress("zm",  &m_eventData->m_zm);
    t1->SetBranchAddress("ed",  &m_eventData->m_energyDeposited);
    t1->SetBranchAddress("ek",  &m_eventData->m_kineticEnergy);
    t1->SetBranchAddress("time",  &m_eventData->m_time);

    t1->SetBranchAddress("type_pri",  &m_eventData->m_primaryParticleType);
    t1->SetBranchAddress("etot_pri", &m_eventData->m_primaryParticleEnergy);

    t1->SetBranchAddress("xm_pri", &m_eventData->m_primaryXm );
    t1->SetBranchAddress("ym_pri", &m_eventData->m_primaryYm);
    t1->SetBranchAddress("zm_pri", &m_eventData->m_primaryZm);

    t1->SetBranchAddress("xp_pri", &m_eventData->m_primaryX);
    t1->SetBranchAddress("yp_pri", &m_eventData->m_primaryY);
    t1->SetBranchAddress("zp_pri", &m_eventData->m_primaryZ);
}


void rebin(double rc, int nd, double r[], double xin[], double xi[])
{
	int i,k=0;
	double dr=0.0,xn=0.0,xo;

	for (i=1;i<nd;i++) {
		while (rc > dr) {
			dr += r[++k];
			xo=xn;
			xn=xi[k];
		}
		dr -= rc;
		xin[i]=xn-(xn-xo)*dr/r[k];
	}
	for (i=1;i<nd;i++) xi[i]=xin[i];
	xi[nd]=1.0;
}




double ran2(long *idum)
{
	
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;

  return 0;
}


double *vec(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
//	if (!v) ;//nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}



void vegas(double regn[], int ndim, double (*fxn)(double [], double), int init,
	unsigned long ncall, int itmx, int nprn, double *tgral, double *sd,
	double *chi2a)
{
	double ran2(long *idum);
	void rebin(double rc, int nd, double r[], double xin[], double xi[]);
	static int i,it,j,k,mds,nd,ndo,ng,npg,ia[MXDIM+1],kg[MXDIM+1];
	static double calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo;
	static double d[NDMX+1][MXDIM+1],di[NDMX+1][MXDIM+1],dt[MXDIM+1],
		dx[MXDIM+1], r[NDMX+1],x[MXDIM+1],xi[MXDIM+1][NDMX+1],xin[NDMX+1];
	static double schi,si,swgt;

	if (init <= 0) {
		mds=ndo=1;
		for (j=1;j<=ndim;j++) xi[j][1]=1.0;
	}
	if (init <= 1) si=swgt=schi=0.0;
	if (init <= 2) {
		nd=NDMX;
		ng=1;
		if (mds) {
			ng=(int)pow(ncall/2.0+0.25,1.0/ndim);
			mds=1;
			if ((2*ng-NDMX) >= 0) {
				mds = -1;
				npg=ng/NDMX+1;
				nd=ng/npg;
				ng=npg*nd;
			}
		}
		for (k=1,i=1;i<=ndim;i++) k *= ng;
		npg=IMAX(ncall/k,2);
		calls=npg*k;
		dxg=1.0/ng;
		for (dv2g=1,i=1;i<=ndim;i++) dv2g *= dxg;
		dv2g=(calls*dv2g)*(calls*dv2g)/npg/npg/(npg-1.0);
		xnd=nd;
		dxg *= xnd;
		xjac=1.0/calls;
		for (j=1;j<=ndim;j++) {
			dx[j]=regn[j+ndim]-regn[j];
			xjac *= dx[j];
		}
		if (nd != ndo) {
			for (i=1;i<=nd;i++) r[i]=1.0;
			for (j=1;j<=ndim;j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
			ndo=nd;
		}
		if (nprn >= 0) {
			printf("%s:  ndim= %3d  ncall= %8.0f\n",
				" Input parameters for vegas",ndim,calls);
			printf("%28s  it=%5d  itmx=%5d\n"," ",it,itmx);
			printf("%28s  nprn=%3d  ALPH=%5.2f\n"," ",nprn,ALPH);
			printf("%28s  mds=%3d  nd=%4d\n"," ",mds,nd);
			for (j=1;j<=ndim;j++) {
				printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
					" ",j,regn[j],j,regn[j+ndim]);
			}
		}
	}
	for (it=1;it<=itmx;it++) {
		ti=tsi=0.0;
		for (j=1;j<=ndim;j++) {
			kg[j]=1;
			for (i=1;i<=nd;i++) d[i][j]=di[i][j]=0.0;
		}
		for (;;) {
			fb=f2b=0.0;
			for (k=1;k<=npg;k++) {
				wgt=xjac;
				for (j=1;j<=ndim;j++) {
					xn=(kg[j]-ran2(&idum))*dxg+1.0;
					ia[j]=IMAX(IMIN((int)(xn),NDMX),1);
					if (ia[j] > 1) {
						xo=xi[j][ia[j]]-xi[j][ia[j]-1];
						rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
					} else {
						xo=xi[j][ia[j]];
						rc=(xn-ia[j])*xo;
					}
					x[j]=regn[j]+rc*dx[j];
					wgt *= xo*xnd;
				}
				f=wgt*(*fxn)(x,wgt);
				f2=f*f;
				fb += f;
				f2b += f2;
				for (j=1;j<=ndim;j++) {
					di[ia[j]][j] += f;
					if (mds >= 0) d[ia[j]][j] += f2;
				}
			}
			f2b=sqrt(f2b*npg);
			f2b=(f2b-fb)*(f2b+fb);
			if (f2b <= 0.0) f2b=TINY;
			ti += fb;
			tsi += f2b;
			if (mds < 0) {
				for (j=1;j<=ndim;j++) d[ia[j]][j] += f2b;
			}
			for (k=ndim;k>=1;k--) {
				kg[k] %= ng;
				if (++kg[k] != 1) break;
			}
			if (k < 1) break;
		}
		tsi *= dv2g;
		wgt=1.0/tsi;
		si += wgt*ti;
		schi += wgt*ti*ti;
		swgt += wgt;
		*tgral=si/swgt;
		*chi2a=(schi-si*(*tgral))/(it-0.9999);
		if (*chi2a < 0.0) *chi2a = 0.0;
		*sd=sqrt(1.0/swgt);
		tsi=sqrt(tsi);
		if (nprn >= 0) {
			printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
				" iteration no.",it,ti,tsi);
			printf("%s integral =%14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
				" all iterations:  ",*tgral,*sd,*chi2a);
			if (nprn) {
				for (j=1;j<=ndim;j++) {
					printf(" DATA FOR axis  %2d\n",j);
					printf("%6s%13s%11s%13s%11s%13s\n",
						"X","delta i","X","delta i","X","delta i");
					for (i=1+nprn/2;i<=nd;i += nprn+2) {
						printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
							xi[j][i],di[i][j],xi[j][i+1],
							di[i+1][j],xi[j][i+2],di[i+2][j]);
					}
				}
			}
		}
		for (j=1;j<=ndim;j++) {
			xo=d[1][j];
			xn=d[2][j];
			d[1][j]=(xo+xn)/2.0;
			dt[j]=d[1][j];
			for (i=2;i<nd;i++) {
				rc=xo+xn;
				xo=xn;
				xn=d[i+1][j];
				d[i][j] = (rc+xn)/3.0;
				dt[j] += d[i][j];
			}
			d[nd][j]=(xo+xn)/2.0;
			dt[j] += d[nd][j];
		}
		for (j=1;j<=ndim;j++) {
			rc=0.0;
			for (i=1;i<=nd;i++) {
				if (d[i][j] < TINY) d[i][j]=TINY;
				r[i]=pow((1.0-d[i][j]/dt[j])/
					(log(dt[j])-log(d[i][j])),ALPH);
				rc += r[i];
			}
			rebin(rc/xnd,nd,r,xin,xi[j]);
		}
	}
}

double (*vegas_function)(double []);

double trans_function(double pt[],double wgt)  //function parameter start from 1; wgt is not used;
{
	return vegas_function(pt+1);
}

class vegas_dimension{
public:
int number_of_dimension;
int	ncall;   //increase to get high precision 
int	itmax;    //increase to get high precision 
double *regn;
void set_dimension(int i){
number_of_dimension=i;
regn=vec(1,number_of_dimension*2+3);
}
void set_dim_start_stop(int dim,double start,double stop){
if(dim>=number_of_dimension) return;
regn[dim+1]=start;
regn[dim+1+number_of_dimension]=stop;
}

double (*afxn)(double [], double);


vegas_dimension(int dim){
set_dimension(dim);
afxn=trans_function;
		ncall=1500;  
		itmax=100;    
}
~vegas_dimension(){
free(regn);
}
};

double vegas(vegas_dimension* dv){


	int init,itmax,ncall,nprn;
	double avgi,chi2a,sd;	
	
		ncall=dv->ncall;   //increase to get high precision 
		itmax=dv->itmax;    //increase to get high precision 
		
		
		nprn=-1;		//debug switch
		init = -1;		//for different initial data mode,shouldn't affect result much
		avgi=sd=chi2a=0.0;



		vegas(dv->regn,dv->number_of_dimension,dv->afxn,init,ncall,itmax,nprn,&avgi,&sd,&chi2a);
		// printf("Number of iterations performed: %d\n",itmax);
		// printf("Integral, Standard Dev., Chi-sq. = %12.6f %12.6f% 12.6f\n",
			// avgi,sd,chi2a);


	//printf("Normal completion\n");
return  avgi;
}




class cluster_chain;
class cluster{
public:
cluster(){
next=NULL;
prev=NULL;
chain=NULL;
clear();
}

vector <int> x;
vector <int> z;
vector <double> q;

double xmean;
double zmean;
double qsum;
int size;
double ds;

cluster * next;
cluster * prev;
cluster_chain * chain;
void clear(){
x.clear();
z.clear();
next=NULL;
prev=NULL;
chain=NULL;
size=0;
xmean=0;
zmean=0;
};

void build(){
// assert(qsum!=0);

xmean=xmean/qsum;
zmean=zmean/qsum;
size=x.size();
}

void fill(int ix,int iz,int iq){
x.push_back(ix);
z.push_back(iz);
q.push_back(iq);
xmean+=double(ix*iq);
zmean+=double(iz*iq);
qsum+=double(iq);
}
int get_size(){
return size;
}
};

class cluster_chain{
public:
cluster_chain(){
clear();}
~cluster_chain(){}

void add(cluster *p){

    if(first==NULL){
	first=p;
	}else{
	p->prev=last;
	last->next=p;
	}
	last=p;
    last->chain=this;
	nhits++;

	if(nhits==2){
	fs=atan2(first->zmean-last->zmean,first->xmean-last->xmean);
	ls=atan2(last->zmean-first->zmean,last->xmean-first->xmean);
	}else if(nhits>2){
	ls=atan2(last->zmean-last->prev->zmean,last->xmean-last->prev->xmean);
	}
}

void add_back(cluster *p){

    if(last==NULL){
	last=p;
	}else{
	p->next=first;
	first->prev=p;
	}
	first=p;
    first->chain=this;
	nhits++;
	
	if(nhits==2){
	fs=atan2(first->zmean-last->zmean,first->xmean-last->xmean);
	ls=atan2(last->zmean-first->zmean,last->xmean-first->xmean);
	}else if(nhits>2){
	fs=atan2(first->zmean-first->next->zmean,first->xmean-first->next->xmean);
	}
}

void clear(){
first=NULL;last=NULL;nhits=0;}

void del(){
cluster * p=first;
for(int i=0;i<nhits;i++){
first=p->next;
p->next=NULL;
p->prev=NULL;
p=first;
}
clear();
}

cluster * first;
cluster * last;
double fs;
double ls;
int nhits;



   double    s11Xy  ;
   double    s12Xy  ;
   double    s22Xy  ;
   double    g1Xy   ;
   double    g2Xy   ;     	
   double    ddXy, a1Xy, a2Xy ;
   double    ddXx;
   int line_nHits;
    double    k;  
    double    b;  
void   reset_line(){
       s11Xy  =0;
       s12Xy  =0;
       s22Xy  =0;
       g1Xy   =0;
       g2Xy   =0;     	
       ddXy=0;
	   a1Xy=0;
	   a2Xy =0;
	   line_nHits=0;
}   
void   add_line(double wxy, double xp, double yp){
line_nHits++;
  s11Xy = s11Xy + wxy ;
  s12Xy = s12Xy + wxy * xp ;
  s22Xy = s22Xy + wxy * xp * xp ;
  g1Xy  = g1Xy  + wxy * yp ;
  g2Xy  = g2Xy  + wxy * xp * yp ;   
   }
void    build_line(){
     ddXy  = s11Xy * s22Xy - s12Xy * s12Xy ;
     if ( ddXy != 0 ) {
        a1Xy  = ( g1Xy * s22Xy - g2Xy * s12Xy ) / ddXy ;//b
		b=a1Xy;
        a2Xy  = ( g2Xy * s11Xy - g1Xy * s12Xy ) / ddXy ;//k
		k=a2Xy;
     }
     else {
		 //ddXx=s12Xy/nHits;
          //cout<<"ERR track:add: ddXy = 0"<<endl;
     }
   }  

   
   
double    check_line(double xp, double yp){
  double   temp;
  double   dxy;
     if ( ddXy != 0 ) {
   temp = (a2Xy * xp - yp + a1Xy) ;
   dxy  = temp * temp / ( a2Xy * a2Xy + 1.F ) ;
     }
     else {
   dxy=fabs(s12Xy/s11Xy-xp);
   dxy*=dxy;
     }
	 
	return dxy;   
	}
};




class digitization_tree_data{
public:
digitization_tree_data(int xd=126,int zd=626){
x_dim=xd;
z_dim=zd;
X_bin=new vector<int>;
Z_bin=new vector<int>;
Digit_adc=new vector<int>;
time1=new vector<int>;
data=new int[x_dim*z_dim];
set_max();
}
~digitization_tree_data(){
delete X_bin;
delete Z_bin;
delete Digit_adc;
delete time1;
delete []data;
}

void set_max(double ml=400,double ms=0.75*M_PI){
maxL=ml;
maxSita=ms;
}
void clear(){
    X_bin->clear();
    Z_bin->clear();
    Digit_adc->clear();
    time1->clear();
	clu.clear();
}

vector<int> *Z_bin;
vector<int> *X_bin;
vector<int> *Digit_adc;
vector<int> *time1;
int Eventid;

int * data;
vector <cluster> clu;
vector <cluster_chain> chain;
int x_dim;
int z_dim;
double maxL;
double maxSita;

double distance(cluster *p,cluster *q, double L){
double dx=fabs(p->xmean-q->xmean);
if(dx>L) return -1;
double dz=fabs(p->zmean-q->zmean);
if(dz>L) return -1;
return sqrt(dx*dx+dz*dz);
}

cluster * seed(){
	for(int i=0;i<clu.size();i++){
	if(clu[i].chain!=NULL) continue;
    clu[i].ds=-1;
    return 	&(clu[i]);
	}
    return 	NULL;
}

cluster * closest(cluster *p){
    cluster *mp=NULL;
	double min=1e10;

	double dis;
	for(int i=0;i<clu.size();i++){
	if(clu[i].chain!=NULL) continue;
	if((&(clu[i]))==p) continue;
    dis=distance(p,&(clu[i]),maxL);
	//cout<<"dis = "<<dis<<endl;
	if(dis<0) continue;
	if(min>dis) {
	min=dis;
	mp=&(clu[i]);
	mp->ds=-1;
	}
	}
return 	mp;
}

cluster * small_dis_angle(cluster *p,double s){
    cluster *mp=NULL;
	double min=1e10;
	double dis;
	double ls;
	for(int i=0;i<clu.size();i++){
	if(clu[i].chain!=NULL) continue;
	if((&(clu[i]))==p) continue;
	
    dis=distance(p,&(clu[i]),maxL);
	if(dis<0) continue;

	
	ls=atan2(clu[i].zmean-p->zmean,clu[i].xmean-p->xmean);
	ls-=s;
	if(ls>M_PI) ls-=2*M_PI;
	if(ls<=-M_PI) ls+=2*M_PI;
	ls=fabs(ls);
	if(ls>maxSita) continue;	
	
	if(min>dis*ls) {
	min=dis*ls;
	mp=&(clu[i]);
	mp->ds=ls;
	}
	}
return 	mp;
}

void build_chain(){
    chain.clear();

    cluster_chain *p=new cluster_chain();
	cluster *se;
	cluster *cl;
	cluster *sda;
	do{
	p->clear();
	se=seed();
	if(se==NULL) break;
	//cout<<"hehe  "<<p->nhits<<endl;	
	p->add(se);
	cl=closest(p->last);
	if(cl==NULL) break;
	//cout<<"hehe1  "<<p->nhits<<endl;	
	p->add(cl);

		do{
		sda=small_dis_angle(p->last,p->ls);
		if(sda==NULL) break;
		p->add(sda);
		}while(sda!=NULL);

		do{
		sda=small_dis_angle(p->first,p->fs);
		if(sda==NULL) break;
		p->add(sda);
		}while(sda!=NULL);
    
	chain.push_back(*p);
    	
	}while(se!=NULL);
	
}
void load_to_data(){
    for(int j=0;j<x_dim*z_dim;j++){
	data[j]=0;
    }    
	for(int i=0;i<Z_bin->size();i++){
	int px=(*X_bin)[i];
	int pz=(*Z_bin)[i];
	data[pz*x_dim+px]=(*Digit_adc)[i];
	}

}

int get_from_data(int* data,int x,int z){

//cout<<dat[x][z]<<"\t"<<x<<"\t"<<z<<endl;
return data[x+z*x_dim];
 
};

void  set_to_data(int* data,int x,int z,int val){
 data[x+z*x_dim]=val;
};

void seed_cluster(int * data, cluster *clu,int x,int z){
	if(x<0) return;
	if(x>=x_dim) return;
	if(z<0) return;
	if(z>=z_dim) return;
	int dat=get_from_data(data,x,z);
    if(dat==0) return;
	
	clu->fill(x,z,dat);
	//cout<<x<<" == "<<z<<" $$ "<<get_from_data(data,x,z)<<endl;
	
	set_to_data(data,x,z,0);
    
	seed_cluster(data,clu,x-1,z+0);
	seed_cluster(data,clu,x-1,z-1);
	seed_cluster(data,clu,x+0,z-1);
	seed_cluster(data,clu,x+1,z-1);
	seed_cluster(data,clu,x+1,z+0);
	seed_cluster(data,clu,x+1,z+1);
	seed_cluster(data,clu,x+0,z+1);
	seed_cluster(data,clu,x-1,z+1);

};
void finding(int minsize=10){
    clu.clear();

    cluster *p=new cluster;
	for(int i=0;i<Z_bin->size();i++){
	if(get_from_data(data,(*X_bin)[i],(*Z_bin)[i])!=0){
	p->clear();
	seed_cluster(data, p,(*X_bin)[i],(*Z_bin)[i]);
    if(p->x.size()>minsize){
	p->build();
	clu.push_back(*p);
	}
	}
	}


}
void set_digitization_tree(TTree *t1){
	t1->SetBranchAddress("Time", &time1);
    t1->SetBranchAddress("Digit_z", &Z_bin);
    t1->SetBranchAddress("Digit_x", &X_bin);
    t1->SetBranchAddress("Digit_adc", &Digit_adc); 
    t1->SetBranchAddress("Eventid", &Eventid);
}
};


class paras
{
public:
	paras(void);
	~paras(void);
int file_size(char* fn);
int para_load(char* fn,char* pd, int length);
int para_number(char* in);
int para_value(char* in,char**v);
int para(char *a);
int run(char* para_file);
char * v(int i);

int gargc;
int cp;
int help;
char** gargv;
char* para_data;
int para_size;

};

paras::paras(void)
{
	para_data=0;
	help=0;
}

paras::~paras(void)
{
	if(para_data!=0) {
		delete []para_data;
		delete []gargv;
	}
}

int paras::file_size(char* fn){FILE *p =fopen(fn,"rb");fseek(p,0L,SEEK_END);int fz=ftell(p);fclose(p);	return fz;}
int paras::para_load(char* fn,char* pd, int length){FILE *fp=fopen(fn,"rb");if(fp==NULL){printf("FILE %s Is Not Open",fn);return 0;}fread(pd,1,length,fp);fclose(fp);	return 1;}

int paras::para_number(char* in){int len=strlen(in);int i=0;int ret=0;for(i=0;i<len-1;i++){if(i==0&&in[i]!=' '&&in[i]!='\t'&&in[i]!='\r'&&in[i]!='\n'){ret++;
continue;}if(in[i  ]==' '||in[i  ]=='\t'||in[i  ]=='\r'||in[i  ]=='\n'){if(in[i+1]!=' '&&in[i+1]!='\t'&&in[i+1]!='\r'&&in[i+1]!='\n'){ret++;}}}	return ret;}

int paras::para_value(char* in,char**v){	int len=strlen(in);
	int i=0;	int ret=0;
	for(i=0;i<len-1;i++)	{
		if(i==0&&in[i]!=' '&&in[i]!='\t'&&in[i]!='\r'&&in[i]!='\n')
		{v[ret]=in+i;ret++;continue;}
		if(in[i  ]==' '||in[i  ]=='\t'||in[i  ]=='\r'||in[i  ]=='\n')
		{	in[i]=0;
			if(in[i+1]!=' '&&in[i+1]!='\t'&&in[i+1]!='\r'&&in[i+1]!='\n')
			{v[ret]=in+i+1;ret++;}
		}}return ret;}
int paras::para(char *a){if(help==1)cout<<a<<" ";int i;for(i=0;i<gargc;i++){if(strcmp(a,gargv[i])==0) { cp=i; return i;}}cp=-1;return -1;}

int paras::run(char* para_file){
	para_size=file_size(para_file);	
	para_data=new char[para_size+3];
	para_data[0]='e';
	para_data[1]=' ';
	para_data[para_size+2]=0;
	para_load(para_file,para_data+2,para_size);
	gargc=para_number(para_data);gargv=new char*[gargc];para_value(para_data,gargv);
	
	return 0;
}

char * paras::v(int i){
if(cp+i>=gargc) return 0;
return gargv[cp+i];
}



class weight
{

public:
	weight(void){};
	~weight(void){};
	weight * next;

    int z;
    int x;
    int adc;
    int time;

};


class ws_hash
{
public:
	ws_hash(void){};
	ws_hash(int NW,int NS,double WS,double WE,double SS,double SE){
	set( NW, NS, WS, WE, SS, SE);
	};

	~ws_hash(void){
	for(int i=0;i<nw;i++){
    delete [] pf[i];
    delete [] pl[i];
    delete [] num[i];
	}
    delete [] pf;
    delete [] pl;	
    delete [] num;	
	};
	void set(int NW,int NS,double WS,double WE,double SS,double SE);
	void set(ws_hash* a);
	int get_d1_bin(double w);
	int get_d2_bin(double s);
	void fill(weight *a);
	void clear();

	int nw;
	int ns;
	double ws;
	double we;
	double ss;
	double se;
	weight ***pf;
	weight ***pl;
	int **num;

};




	void ws_hash::set(ws_hash * a){
		nw=a->nw;
		ns=a->ns;
		ws=a->ws;
		we=a->we;
		ss=a->ss;
		se=a->se;

	pf=new weight**[nw];
	pl=new weight**[nw];
	num=new int*[nw];
	for(int i=0;i<nw;i++){
    pf[i]=new weight*[ns];
    pl[i]=new weight*[ns];
	num[i]=new int[ns];
	for(int j=0;j<ns;j++){
	pf[i][j]=NULL;
	pl[i][j]=NULL;
	num[i][j]=0;
	}//for(int j=0;
	}//for(int i=0;
	};//void set

	void ws_hash::set(int NW,int NS,double WS,double WE,double SS,double SE){
	nw=NW;
	ns=NS;
	ws=WS;
	we=WE;
	ss=SS;
	se=SE;
	pf=new weight**[nw];
	pl=new weight**[nw];
	num=new int*[nw];
	for(int i=0;i<nw;i++){
    pf[i]=new weight*[ns];
    pl[i]=new weight*[ns];
	num[i]=new int[ns];
	for(int j=0;j<ns;j++){
	pf[i][j]=NULL;
	pl[i][j]=NULL;
	num[i][j]=0;
	}//for(int j=0;
	}//for(int i=0;
	};//void set

	int ws_hash::get_d1_bin(double w){
	assert(we>ws);
	if(w>=we) return nw-1;
	if(w<ws) return 0;
	return (w-ws)*nw/(we-ws);
	}
	int ws_hash::get_d2_bin(double s){
	assert(se>ss);
	if(s>=se) return ns-1;
	if(s<ss) return 0;
	return (s-ss)*ns/(se-ss);
	}

	void ws_hash::fill(weight *a){
		int wi=get_d1_bin(a->x);
		int si=get_d2_bin(a->z);
	if(pf[wi][si]==NULL){
	pf[wi][si]=a;
	}else{
	pl[wi][si]->next=a;
	}//if(pf[wi][si]==NULL)
	pl[wi][si]=a;
	num[wi][si]++;
	}//void fill(weight *a)

	void ws_hash::clear(){
	for(int i=0;i<nw;i++){
	for(int j=0;j<ns;j++){
	pf[i][j]=NULL;
	pl[i][j]=NULL;
	num[i][j]=0;
	}//for(int j=0;
	}//for(int i=0;	
	}//void ws_hash::clear()








#endif	