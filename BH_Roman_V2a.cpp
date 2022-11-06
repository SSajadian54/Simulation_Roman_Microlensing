#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
//#include "VBBinaryLensingLibrary.h"
using namespace std;

#define Nx 7
#define Ny 3

const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double pi= M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double Hp=6.62607004*pow(10.0,-34);
const double KP=3.08568025*pow(10.,19); // in meter.
const double G=6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double M_sun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=365.2421875;
const double binary_fraction=double(2.0/3.0);
const double Avks=double(8.20922);
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;


///============================ Besancon constant ==============///
const double R_sun=8.5;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]={3.1,2.5,3.1,3.1};

///============================ WFIRST & OGLE &  KMTNet ===================
const int M=5;///number of filters  VIKH,W149
const double satu[M]={12.0, 12.0, 13.0, 13.0, 14.8}; //it should be changed
const double thre[M]={20.0, 21.0, 21.0, 21.0, 26.0};
const double FWHM[M]={0.33, 0.33, 0.33, 0.33 , 0.33};//3*pixel_size (0.11") of WFIRST, VIKH W149 Z087
const double AlAv[M]={1.009,0.600,0.118,0.184,0.225};///From besancon model[VIKH W149]
const double sigma[M]={0.022,0.022,0.02,0.025,0.025};//MOAاستفاده از مقاله کاردلی
const double Akv=0.118;
const double cade=double(15.16/(60.0*24.0));//W149_cadence in  day
const double dt=double(10.0/60.0/24.0);//days
const int coun= 41050;
////=======================================================
///**************  Number constant  ********************///
const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double od_disk,od_ThD,od_bulge,od_halo,opd;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart;
    double Nblend[M], blend[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M];
    double type, Tstar, logl, col, Rstar, mass,vs;
    double Romax,ros, deltao;
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double mus1, mus2, mul1, mul2;
    double pos1, pos2, Astar, xi;
    double fb, mbs;
    double ut; 
};
struct lens{
    int numl,struc;
    double Ml, Dl, vl , vs, Vt, xls;
    double rhomaxl,tE,RE;
    double mul, u0;
    double t0;
    double tetE;
    double piE; 
    double magw[Nw], errw[Nw];
};
struct astromet{
   double tetp;
   double omegae;
   double vearth;
   double fact;
   double Ve_n1, Ve_n2;
   double ue_n1, ue_n2;
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[5][N1],Rs_d[N1],mass_d[N1],type_d[N1]; int cl_d[N1];  ///thin disk
    double Teff_b[N2],logl_b[N2],Mab_b[5][N2],Rs_b[N2],mass_b[N2],type_b[N2]; int cl_b[N2];  /// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[5][N3],Rs_t[N3],mass_t[N3],type_t[N3]; int cl_t[N3];  ///thick disk
    double Teff_h[N4],logl_h[N4],Mab_h[5][N4],Rs_h[N4],mass_h[N4],type_h[N4]; int cl_h[N4];  /// halo
};
struct extinc{
   double dis[100];///distance
   double Extks[100];///ks-band extinction
   double Aks;
};
struct covarian{
   double FishM[Nx][Nx];
   double FishA[Ny][Ny];
};
///===================== FUNCTION ===============================================================
int Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm , extinc & ex);
void func_lens(lens & l, source & s);
void vrel(source & s , lens & l);
void Disk_model(source & s);
void optical_depth(source & s);
double Interpol(double ds, extinc & ex);
double RandN(double sigma, double Nn);
double RandR(double down, double up);
int ErrorRoman(lens & l, double mag);
void lightcurve(source & s, lens & l, astromet & as, double tim); 
//////////////////////////////////////////////////////
    time_t _timeNow;
      unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///////////////////////////////////////////////////////
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{

///****************************************************************************
//gettimeofday(&newTime,NULL);
//ftime(&tp);
	time(&_timeNow);
	_randStream = fopen("/dev/urandom", "r");
	_dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
    _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
	srand(_randVal);
///****************************************************************************
     //VBBinaryLensing vbb;
     //vbb.Tol=1.e-3;
     //vbb.LoadESPLTable("./files/ESPL.tbl");
     source s;
     lens l;
     astromet as;
     CMD cm;
     extinc ex;
     covarian co; 
     read_cmd(cm);
///===========================================================

    FILE * wfirstf;
    wfirstf=fopen("./files/sigma_WFIRST.txt","r");
    for(int i=0; i<Nw;++i){
    fscanf(wfirstf,"%lf  %lf\n",&l.magw[i],&l.errw[i]);}
    fclose(wfirstf);
    //cout<<"*****  sigma_wfirst was read *****  "<<endl;
///===========================================================


    FILE* fil1;
    fil1=fopen("./files/MONT/IBH_ROMAN1a.txt","a+");
    fclose(fil1);

    FILE* fil2;
    fil2=fopen("./files/MONT/IBH_ROMAN2a.txt","a+");
    fclose(fil2);
    
    FILE* fil3;
    fil3=fopen("./files/MONT/IBH_ROMAN3a.txt","a+");
    fclose(fil3);


    double obs[6][2]={0.0};
    obs[0][0]=0.0*year+1.0;    obs[0][1]=0.0*year +72.0;
    obs[1][0]=0.0*year+182.0;  obs[1][1]=0.0*year+254.0;
    obs[2][0]=1.0*year+1.0;    obs[2][1]=1.0*year+72.0;
    obs[3][0]=3.0*year+182.0;  obs[3][1]=3.0*year+254.0;
    obs[4][0]=4.0*year+1.0;    obs[4][1]=4.0*year+72.0;
    obs[5][0]=4.0*year+182.0;  obs[5][1]=4.0*year+254.0;


    as.tetp=double(M_PI/3.0);
    as.omegae=double(2.0*M_PI/year); ///radian per day
    as.vearth= as.omegae;
    as.fact=double(1000.0*3600.0*24.0/AU);
  
    
    int flag_det=0;
    int flagf,counter=0;
    int ndw,flago, w1;
    char filnam1[40], filnam2[40];
    double magnio, test;
    double magni[M];
    double timp1,lonn, chi, chi1, tim;
    float flag0, flag1, flag2;
    double dchi, mh;
    double lonn0, lat0, siga, snr, errx, erry;  
    double magw, derm1[2], derm2[2], dera1[2], derb1[2], dera2[2], derb2[2]; 
    double derm1f, derm2f, dera1f, derb1f, dera2f, derb2f; 
    double Delta1[Nx], Delta2[Ny];
    float sig[2]={+1.0, -1.0};  
   
   double *timn=new double[coun];
   double *magn=new double[coun];
   double *soux=new double[coun];
   double *souy=new double[coun];
   double *errm=new double[coun];
   double *erra=new double[coun];
   
   double nsim=0.0; 
   
 
     for(int icon=1; icon<100000; ++icon){
   //  cout<<"***** Step:  "<<icon<<endl;
     for (int li=1; li<=7; li++){
     if(li==1) {lonn0=1.3;  lat0=-0.875; }
     if(li==2) {lonn0=0.9;  lat0=-0.875; }
     if(li==3) {lonn0=1.3;  lat0=-1.625; }
     if(li==4) {lonn0=0.9;  lat0=-1.675; }
     if(li==5) {lonn0=0.5;  lat0=-1.675; }
     if(li==6) {lonn0=0.1;  lat0=-1.675; }
     if(li==7) {lonn0=-0.3; lat0=-1.675; }
     s.lat= lat0+ double(RandR(-0.375,0.375));
     lonn =lonn0+ double(RandR(-0.200,0.200));
  
     if(lonn<=0.0)   s.lon=360.0+lonn;
     else            s.lon=lonn;
     
     s.TET=(360.0-s.lon)/RA;///radian
     s.FI=s.lat/RA;///radian
     Disk_model(s);
     if(int(Extinction(ex,s))==1){

     do{
     func_source(s, cm, ex);
     func_lens(l, s);
     }while(l.tE<=0.1 or l.tE>30000.0);
     optical_depth(s);
    
    
     s.mus1= double(s.SV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
     s.mus2= double(s.SV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(s.Ds*AU);///mas/days
     s.mul1= double(s.LV_n1- s.VSun_n1)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
     s.mul2= double(s.LV_n2- s.VSun_n2)*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days 
     l.piE=double(1.0/l.Dl -1.0/s.Ds)/l.tetE; 
     
    // cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<endl; 
    // cout<<"mul1:  "<<s.mul1<<"\t mul2:  "<<s.mul2<<"\t l.piE:  "<<l.piE<<endl; 
    // cout<<"blending:  "<<s.fb<<"\t m_base:  "<<s.mbs<<endl;   
    
     test=(double)rand()/((double)(RAND_MAX)+(double)(1.0));
     if(test<=s.fb  and  s.mbs<=thre[4] ){
     counter+=1; 
     
     //flagf=1;
     //sprintf(filnam2,"./files/MONT/%c%c%c%c%d.dat",'d','a','t','_', counter);
     //data3=fopen(filnam2,"w");
     
   
///*********************************************************************************   
     for(int i=0; i<coun; ++i){
     timn[i]=0.0; magn[i]=0.0; soux[i]=0.0;  souy[i]=0.0;  errm[i]=0.0;  erra[i]=0.0; }
     flag0=flag1=flag2=0.0;
     ndw=flag_det=0;
     chi=0.0; chi1=0.0;
     timp1=tim=0.0;
     for(tim=0.0;  tim<double(5.0*year); tim+=dt){

     flago=0;
     for(int kl=0; kl<6; ++kl){
     if((obs[kl][0]-tim)*(obs[kl][1]-tim)<=0.0){flago=1; break;}}
     if(flago>0){
     
     timp1 += dt;    
     if(timp1>cade or timp1==cade){
     timp1=double(timp1-cade);
     
     test=(double)rand()*100.0/((double)(RAND_MAX)+(double)(1.0));
     if(test<90.0){
     
     lightcurve(s,l,as,tim); 
     s.Astar=double(s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
     
     for(int i=0;i<M; ++i){
     magni[i] =s.magb[i]-2.5*log10(s.Astar*s.blend[i] + 1.0 - s.blend[i]);}

     if(magni[4]>=satu[4] and magni[4]<=thre[4]){
   
     w1=ErrorRoman(l, magni[4]);
     magnio=magni[4] + RandN(l.errw[w1] ,3.0);
     
     chi  +=fabs( (magnio-     s.mbs)*(magnio-     s.mbs)/(l.errw[w1]*l.errw[w1]));
     chi1 +=fabs( (magnio-  magni[4])*(magnio-  magni[4])/(l.errw[w1]*l.errw[w1]));
     mh =  magni[3] + RandN(l.errw[w1]*3.0 , 3.0);///error in H-band
     snr= pow(10.0,0.4*(26.1-mh));
     siga= double(106.0/snr);///mili arcs
    
     //siga=1.0; ///milli arcsec 
     //errx= RandN(siga ,3.0);
     //erry= RandN(siga ,3.0);
     //fprintf(data3,"%.6lf  %.6lf  %.6lf    %.8lf    %.8lf    %.7lf\n",tim, magnio, l.errw[w1], s.pos1+errx, s.pos2+erry, siga);
     
     timn[ndw]=tim;
     magn[ndw]=magni[4];
     errm[ndw]=l.errw[w1];
     soux[ndw]=s.pos1;
     souy[ndw]=s.pos2;
     erra[ndw]=siga;
    
     
     flag2=0.0;
     if(fabs(magnio-s.mbs)>fabs(4.0*l.errw[w1]))  flag2=1.0;
     if(ndw>1 and float(flag0+flag1+flag2)>2.0)    flag_det=1;
     flag0=flag1;
     flag1=flag2;
     ndw+=1;
     if(ndw>=coun){cout<<"Error ndw:  "<<ndw<<"\t coun:  "<<coun<<endl;  int uue; cin>>uue; }
     }
     }}}}//loop time
    // fclose(data3);
    ///************ WFIRST *********************************************************************
    dchi=fabs(chi-chi1);
  //  cout<<"dchi:  "<<dchi<<"\t flag_det:  "<<flag_det<<"\t ndw:  "<<ndw<<endl;

    if(dchi>800.0 and flag_det>0 and ndw>5){
    fil1=fopen("./files/MONT/IBH_ROMAN1a.txt","a+");
    fprintf(fil1,
   "%d  %.5lf   %.5lf  "///3
   "%d  %.5lf   %.5lf   %.5lf  "///7
   "%d   %d  %.5lf   %.5lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf   %.4lf  %.4lf  "   //21
   "%.5lf  %.9lf  %.5lf  %.9lf   %.2lf  %.2lf   %.4lf  %.4lf   "  //29
   "%.8lf  %.6lf  %.6lf  %.8lf    %.7lf  %.8lf  %.4lf  %.8lf  %.9lf " ///38
   "%d  %d %.1lf  %d   %d   %.9lf  %.9lf   %.9lf  %.9lf  %.9lf   %.6lf\n", //49
   counter,s.lat, s.lon, //3
   l.struc, l.Ml, l.Dl, l.vl, //7
   s.struc, s.cl, s.mass, s.Ds, s.Tstar, s.Rstar, s.logl, s.type, s.col, s.vs, s.Mab[1], s.Mab[4], s.Map[1], s.Map[4],//21
   s.magb[1], s.mbs, s.blend[1], s.fb, s.Nblend[1], s.Nblend[4], s.Ai[1], s.Ai[4], //29
   l.tE, l.RE/AU, l.t0, l.mul, l.Vt, l.u0, s.opd*1.0e6, s.ros, l.tetE,//38
   flagf, flag_det, dchi, ndw,li, s.mus1, s.mus2, s.xi, s.mul1, s.mul2, l.piE);//49
   fclose(fil1);
  // cout<<"************* End of saving in the file ***********************"<<endl;

///####################################################################################
   Delta1[5]=0.05;///for mbs   
   Delta1[3]=0.1; //for xi in radian 
   if(l.t0>10.0) Delta1[0]=+10.0; 
   else          Delta1[0]=float(l.t0-2.0); 
   
   if(l.u0>0.1)  Delta1[1]=0.1; 
   else          Delta1[1]=float(l.u0*0.4); 
            
   if(l.tE>5.0)  Delta1[2]=5.0;          
   else          Delta1[2]=float(l.tE*0.4); 
   
   
   if(s.fb>0.1 and s.fb<0.9) Delta1[4]=0.1;
   else if (s.fb<=0.1)       Delta1[4]=float(s.fb*0.5);
   else if (s.fb>=0.9)       Delta1[4]=float(1.0-s.fb)*0.5;
   
   if(l.piE>0.01) Delta1[6]=0.01;
   else           Delta1[6]=float(l.piE*0.3);  
  
    
   if(l.tetE>0.2) Delta2[0]=0.2;
   else           Delta2[0]=l.tetE*0.4;   
   Delta2[1]= s.mus1*0.25;  
   Delta2[2]= s.mus2*0.25; 
   
   
   for(int j=0; j<Nx; ++j){
   for(int k=0; k<Nx; ++k){co.FishM[j][k]=0.0;}}
   for(int j=0; j<Ny; ++j){
   for(int k=0; k<Ny; ++k){co.FishA[j][k]=0.0;}}
  
  
  
   for(int i=0; i<ndw; ++i){
   
   for(int j=0; j<Nx; ++j){
   
   for(int h=0; h<2; ++h){
   if(j==0)  l.t0 += double(Delta1[0]*sig[h]);
   if(j==1)  l.u0 += double(Delta1[1]*sig[h]);
   if(j==2)  l.tE += double(Delta1[2]*sig[h]);
   if(j==3)  s.xi += double(Delta1[3]*sig[h]);
   if(j==4)  s.fb += double(Delta1[4]*sig[h]);
   if(j==5)  s.mbs+= double(Delta1[5]*sig[h]);
   if(j==6)  l.piE+= double(Delta1[6]*sig[h]);
   if(j==7)  s.ros+= double(Delta1[7]*sig[h]);
   
   lightcurve(s,l,as,timn[i]);
   
   //if(s.ut<float(20.0*s.ros)) s.Astar= vbb.ESPLMag2(s.ut , s.ros);
   //else                       
   s.Astar= (s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   
   magw= s.mbs - 2.5*log10(s.Astar*s.fb + 1.0 - s.fb); 
   derm1[h]= float(magw-magn[i])/(Delta1[j]*sig[h]);
   if(j==0) l.t0 -= double(Delta1[0]*sig[h]);  
   if(j==1) l.u0 -= double(Delta1[1]*sig[h]);  
   if(j==2) l.tE -= double(Delta1[2]*sig[h]);  
   if(j==3) s.xi -= double(Delta1[3]*sig[h]);  
   if(j==4) s.fb -= double(Delta1[4]*sig[h]);   
   if(j==5) s.mbs-= double(Delta1[5]*sig[h]); 
   if(j==6) l.piE-= double(Delta1[6]*sig[h]);
   if(j==7) s.ros-= double(Delta1[7]*sig[h]);}
   derm1f= double(derm1[0]+derm1[1])*0.5; 
   
   
   for(int k=0; k<=j; ++k){
   
   for(int h=0; h<2; ++h){
   if(k==0) l.t0+= double(Delta1[0]*sig[h]);
   if(k==1) l.u0+= double(Delta1[1]*sig[h]);
   if(k==2) l.tE+= double(Delta1[2]*sig[h]);
   if(k==3) s.xi+= double(Delta1[3]*sig[h]);
   if(k==4) s.fb+= double(Delta1[4]*sig[h]);
   if(k==5) s.mbs+=double(Delta1[5]*sig[h]);
   if(k==6) l.piE+=double(Delta1[6]*sig[h]);
   if(k==7) s.ros+=double(Delta1[7]*sig[h]);
   lightcurve(s,l,as,timn[i]);
   
   //if(s.ut<float(20.0*s.ros)) s.Astar= vbb.ESPLMag2(s.ut , s.ros);
   //else                       
   s.Astar= (s.ut*s.ut+2.0)/sqrt(s.ut*s.ut*(s.ut*s.ut+4.0));
   
   magw= s.mbs - 2.5*log10(s.Astar*s.fb + 1.0 - s.fb); 
   derm2[h]= float(magw-magn[i])/(Delta1[k]*sig[h]);
   if(k==0) l.t0 -= double(Delta1[0]*sig[h]);  
   if(k==1) l.u0 -= double(Delta1[1]*sig[h]);  
   if(k==2) l.tE -= double(Delta1[2]*sig[h]);  
   if(k==3) s.xi -= double(Delta1[3]*sig[h]);  
   if(k==4) s.fb -= double(Delta1[4]*sig[h]);   
   if(k==5) s.mbs-= double(Delta1[5]*sig[h]); 
   if(k==6) l.piE-= double(Delta1[6]*sig[h]);
   if(k==7) s.ros-= double(Delta1[7]*sig[h]);}
   derm2f= double(derm2[0]+derm2[1])*0.5;
   
   co.FishM[j][k] += derm1f*derm2f/(errm[i]*errm[i]);}
   
   }//end of for J
   }//end of data for
   
   for(int j=0; j<Nx; ++j){
   for(int k=(j+1);k<Nx; ++k) co.FishM[j][k]= co.FishM[k][j]; }
   //for(int j=0; j<Nx; ++j){
  // for(int k=0; k<=j; ++k)   cout << co.FishM[j][k]<<"    ";
  // cout << endl;}
   
  
    //inverse1(co);
    fil2=fopen("./files/MONT/IBH_ROMAN2a.txt","a+");
    fprintf(fil2,"%d  %.7lf  %.7lf %.7lf  %.7lf %.7lf   %.7lf  %.7lf  %.8e   %.8e %.8e  %.8e %.8e %.8e   %.8e %.8e %.8e %.8e   %.8e %.8e  %.8e %.8e %.8e   %.8e %.8e %.8e %.8e %.8e %.8e  %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
    counter,l.t0, l.u0, l.tE, s.xi, s.fb, s.mbs, l.piE, 
    co.FishM[0][0], 
    co.FishM[1][0], co.FishM[1][1], 
    co.FishM[2][0], co.FishM[2][1], co.FishM[2][2],
    co.FishM[3][0], co.FishM[3][1], co.FishM[3][2],co.FishM[3][3],
    co.FishM[4][0], co.FishM[4][1], co.FishM[4][2],co.FishM[4][3],co.FishM[4][4],
    co.FishM[5][0], co.FishM[5][1], co.FishM[5][2],co.FishM[5][3],co.FishM[5][4],co.FishM[5][5], 
    co.FishM[6][0], co.FishM[6][1], co.FishM[6][2],co.FishM[6][3],co.FishM[6][4],co.FishM[6][5], co.FishM[6][6]); 
    //co.FishM[7][0], co.FishM[7][1], co.FishM[7][2],co.FishM[7][3],co.FishM[7][4],co.FishM[7][5], co.FishM[7][6], co.FishM[7][7]);//44
    fclose(fil2); 
    //%.8e %.8e %.8e %.8e %.8e %.8e %.8e  %.8e\n",
   
   
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   

   
   for(int i=0; i<ndw; ++i){
   
   for (int j=0; j<Ny; ++j) {
   
   for(int h=0; h<2; ++h){
   if(j==0) l.tetE += double(Delta2[0]*sig[h]);
   if(j==1) s.mus1 += double(Delta2[1]*sig[h]);
   if(j==2) s.mus2 += double(Delta2[2]*sig[h]);
   lightcurve(s,l,as,timn[i]);
   dera1[h]= float(s.pos1-soux[i])/(Delta2[j]*sig[h]);  
   derb1[h]= float(s.pos2-souy[i])/(Delta2[j]*sig[h]);
   if(j==0) l.tetE -= double(Delta2[0]*sig[h]);
   if(j==1) s.mus1 -= double(Delta2[1]*sig[h]); 
   if(j==2) s.mus2 -= double(Delta2[2]*sig[h]);}
   dera1f=(dera1[0]+ dera1[1])*0.5;  
   derb1f=(derb1[0]+ derb1[1])*0.5;  
  
   
   for(int k=0; k<=j; ++k){   
   
   for(int h=0; h<2; ++h){
   if(k==0) l.tetE += double(Delta2[0]*sig[h]);
   if(k==1) s.mus1 += double(Delta2[1]*sig[h]);
   if(k==2) s.mus2 += double(Delta2[2]*sig[h]);
   lightcurve(s,l,as,timn[i]);
   dera2[h]= float(s.pos1-soux[i])/(Delta2[k]*sig[h]);  
   derb2[h]= float(s.pos2-souy[i])/(Delta2[k]*sig[h]);  
   if(k==0) l.tetE -= double(Delta2[0]*sig[h]);
   if(k==1) s.mus1 -= double(Delta2[1]*sig[h]); 
   if(k==2) s.mus2 -= double(Delta2[2]*sig[h]);}
   dera2f=(dera2[0]+ dera2[1])*0.5;  
   derb2f=(derb2[0]+ derb2[1])*0.5;  
   
   co.FishA[j][k] += (dera1f*dera2f + derb1f*derb2f)/(erra[i]*erra[i]);}
  
   }//end of loop J
   }//end of loop data
   
   for(int j=0; j<Ny; ++j){
   for(int k=(j+1);k<Ny; ++k) co.FishA[j][k]= co.FishA[k][j];}
 //  for(int j=0; j<Ny; ++j){
 //  for(int k=0; k<=j; ++k)  cout << co.FishA[j][k] << "     ";
 //  cout << endl;}
   
  
   //inverse2(co);
   fil3=fopen("./files/MONT/IBH_ROMAN3a.txt","a+");
   fprintf(fil3, "%d  %.7lf    %.7lf    %.7lf   %.9e    %.9e    %.9e   %.9e   %.9e   %.9e\n",
   counter,l.tetE, s.mus1, s.mus2,
   co.FishA[0][0], 
   co.FishA[1][0], co.FishA[1][1], 
   co.FishA[2][0], co.FishA[2][1], co.FishA[2][2]);
   fclose(fil3); 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/*
   cout<<"Inverse Matrix Magnification ****************************"<<endl;
   for(int i=0;i<Nx; ++i){
   for(int j=0;j<Nx; ++j)  cout << co.InvM[i][j] << "     ";
   cout << endl;}
   cout<<"Inverse Matrix Astrometry ***************** **************"<<endl;
   for(int i=0;i<Ny; ++i){
   for(int j=0;j<Ny; ++j)    cout << co.InvA[i][j] << "     ";
   cout << endl;}
   cout<<"************************************************"<<endl;*/
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH   
   
   
   if(int(counter)%100==0){
   cout<<"============================================================="<<endl;
   cout<<"counter:     "<<counter<<"\t icon:  "<<icon<<endl;
   cout<<"lat:  "<<s.lat<<"\t lon:  "<<s.lon<<endl;
   cout<<"********************** SOURCE **************************"<<endl;
   cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<"\t strucs:  "<<s.struc<<endl;
   cout<<"cl:  "<<s.cl<<"\t mass:  "<<s.mass<<"\t Tstar:  "<<s.Tstar<<endl;
   cout<<"Rstar:  "<<s.Rstar<<"\t logl:  "<<s.logl<<"\t type:  "<<s.type<<endl;
   cout<<"col:  "<<s.col<<"\t vs:  "<<s.vs<<"\t Mag_I:  "<<s.Mab[1]<<"\t Mag_W149:  "<<s.Mab[4]<<endl;
   cout<<"********************** LENS **************************"<<endl;
   cout<<"Dl:  "<<l.Dl<<"\t numl:  "<<l.numl<<"\t strucl:  "<<l.struc<<endl;
   cout<<"mass:  "<<l.Ml<<"vl:  "<<l.vl<<endl;
   cout<<"*********************** LENSING ***********************"<<endl;
   cout<<"tE:  "<<l.tE<<"\t RE(AU):  "<<l.RE/AU<<"\t t0:  "<<l.t0<<endl;
   cout<<"Vt:  "<<l.Vt<<"\t mul(mas/days):  "<<l.mul<<"\t u0:  "<<l.u0<<endl;
   cout<<"************ WFIRST & OGLE ***************************"<<endl;
   cout<<"flag_det:  "<<flag_det<<"\t dchi:  "<<dchi<<"\t ndw:  "<<ndw<<endl;
   cout<<"t0:  "<<l.t0<<"\t u0:  "<<l.u0<<"\t tE:  "<<l.tE<<endl;
   cout<<"ksi:  "<<s.xi<<"\t fb:  "<<s.fb<<"\t mbs:  "<<s.mbs<<endl;
   cout<<"piE:   "<<l.piE<<"\t l.tetE:  "<<l.tetE <<endl;
   cout<<"mus1:  "<<s.mus1<<"\t mus2:  "<<s.mus2<<endl;
   cout<<"==============================================================="<<endl; }

  
   }//if wfirst>0
   }//if magnitude of baseline 
   }//end of EXtinction
   }//loop icon
   } ///loop il
   
   delete [] magn, timn, soux, souy; 
   return(0);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
int ErrorRoman(lens & l, double mag){
     int w=-1;
     if(mag<l.magw[0] or mag==l.magw[0]) w=0;
     else if(mag>l.magw[Nw-1] or mag==l.magw[Nw-1]) w=Nw-1;
     else{
     for(int i=1; i<Nw; ++i){
     if(double((mag-l.magw[i])*(mag-l.magw[i-1]))<0.0 or  mag==l.magw[i-1]){w=i-1; break;}}}
     if(w==-1){cout<<"mapt[0]: "<<mag<<"\t w: "<<w<<endl;  int ww; cin>>ww;}
     return(w); 
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
//                                                                    //
//                         Light curves and Astrometry                //
//                                                                    //
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&//
void lightcurve(source & s, lens & l, astromet & as, double timh)
{ 
   double dvex, dvey, Ve_x, tt, ux, uy;
   double int1=0.0, int2=0.0; 
   for (int ig=0; ig<2; ++ig){
   if(ig==0) tt=timh; 
   if(ig==1) tt=0.0;    
   dvex= +as.vearth* sin(as.omegae* tt + M_PI/2.0) / as.omegae;//km/s
   dvey= -as.vearth* cos(as.omegae* tt + M_PI/2.0) / as.omegae;//km/s
   as.Ve_n1= cos(as.tetp)*dvex*sin(s.deltao) - dvey*cos(s.deltao);
   Ve_x=    -cos(as.tetp)*dvex*cos(s.deltao) - dvey*sin(s.deltao);
   as.Ve_n2=-sin(s.FI)*Ve_x + cos(s.FI)*sin(as.tetp)*dvex;
   if(ig==0) {int1 = as.Ve_n1; int2 = as.Ve_n2; }
   if(ig==1) {int1-= as.Ve_n1; int2-= as.Ve_n2; }}
   as.ue_n1= int1; 
   as.ue_n2= int2; 
      
   ux=-l.u0*sin(s.xi) + (timh-l.t0)*cos(s.xi)/l.tE  + l.piE * as.ue_n1;
   uy= l.u0*cos(s.xi) + (timh-l.t0)*sin(s.xi)/l.tE  + l.piE * as.ue_n2;
   
   s.ut=sqrt(ux*ux+ uy*uy);  
   if(s.ut==0.0)  s.ut=1.0e-50; 
   
   tt= (ux-l.piE*as.ue_n1)*(ux-l.piE*as.ue_n1) + (uy-l.piE*as.ue_n2)*(uy-l.piE*as.ue_n2);
   
   s.pos1= -l.u0*l.tetE*sin(s.xi) + s.mus1*(timh-l.t0) + (ux - l.piE*as.ue_n1)*l.tetE/(tt + 2.0); 
   s.pos2= +l.u0*l.tetE*cos(s.xi) + s.mus2*(timh-l.t0) + (uy - l.piE*as.ue_n2)*l.tetE/(tt + 2.0);   
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{ 
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD  &  cm, extinc & ex)
{
    int num,struc,nums,yye;
    double rho,rf;
    double Ds,Ai[M],Av;
    double Mab[M],Map[M];
    double maxnb=0.0;


    for(int i=0; i<M; ++i){
    s.Fluxb[i]=s.Nblend[i]=0.0;
    s.Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    s.Nblend[i]=s.Nblend[i]+RandN(sqrt(s.Nblend[i]),1.5);
    if(s.Nblend[i]<=1.0) s.Nblend[i]=1.0;
    if(s.Nblend[i]>maxnb) maxnb=s.Nblend[i];}
    for(int i=0; i<M; ++i){s.magb[i]=s.Ai[i]=s.Map[i]=s.Mab[i]=0.0;}
    //cout<<"Nblend[4]:  "<<s.Nblend[4]<<"\t maxnb: "<<maxnb<<endl;


    for(int k=1; k<=int(maxnb); ++k){

    do{
    num=int((Num-15.00)*(double)rand()/(double)(RAND_MAX+1.) +10.00);
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);///distance larger than 50.
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}
   // cout<<"Ds:  "<<Ds<<endl;


     rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
     if (rf<= s.rho_disk[nums])                    struc=0;///thin disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums]+ s.rho_halo[nums])) struc=3;///halo
    if(k==1) s.struc=struc;
    //cout<<"struc:  "<<struc<<endl;


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(struc==0){//thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<5; ++i){Mab[i]=cm.Mab_d[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_d[num];
    s.type= cm.type_d[num];
    s.mass= cm.mass_d[num];
    s.Tstar=cm.Teff_d[num];
    s.logl= cm.logl_d[num];
    s.col=Mab[0]-Mab[1];
    s.cl= cm.cl_d[num];}
    if(cm.mass_d[num]<0.0 or int(cm.cl_d[num])>5 or cm.Teff_d[num]<0.0 or cm.type_d[num]>8.0 or cm.mass_d[num]>1000.0){
    cout<<"Error(thin disk) temp: "<<cm.Teff_d[num]<<"\t mass: "<<cm.mass_d[num]<<"\t counter: "<<num<<endl; cin>>yye; }}


    if(struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_b[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_b[num];
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Tstar=cm.Teff_b[num];
    s.col= Mab[0]-Mab[1];
    s.logl=cm.logl_b[num];
    s.cl= cm.cl_b[num];}
    if(cm.mass_b[num]<0.0 or int(cm.cl_b[num])>5 or cm.Teff_b[num]<0.0 or cm.type_b[num]>8.0 or cm.mass_b[num]>10000.0){
    cout<<"Error(bulge) temp: "<<cm.Teff_b[num]<<"\t mass: "<<cm.mass_b[num]<<"\t counter: "<<num<<endl;   cin>>yye; }}




    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){
    s.Rstar= cm.Rs_t[num];
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Tstar=cm.Teff_t[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_t[num];
    s.cl= cm.cl_t[num];}
    if( cm.mass_t[num]<0.0 or int(cm.cl_t[num])>5 or cm.Teff_t[num]<0.0 or cm.type_t[num]>8.0 or cm.mass_t[num]>100000.0){
    cout<<"Error(thick disk) temp: "<<cm.Teff_t[num]<<"\t mass: "<<cm.mass_t[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    if(struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_h[i][num]; }
    if(k==1){
    s.Rstar=cm.Rs_h[num];
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Tstar=cm.Teff_h[num];
    s.col=Mab[0]-Mab[1];
    s.logl=cm.logl_h[num];
    s.cl= cm.cl_h[num];}
    if(cm.mass_h[num]<0.0 or int(cm.cl_h[num])>5 or cm.Teff_h[num]<0.0 or cm.type_h[num]>8.0 or cm.mass_h[num]>10000.0){
    cout<<"Error(Galactic halo) temp: "<<cm.Teff_h[num]<<"\t mass: "<<cm.mass_h[num]<<"\t counter: "<<num<<endl;  cin>>yye;}}



    Mab[4]=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149: (K+H+J)/3
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av<0.0)    Av=0.0;

    //if(Av>20.0 or Av<-0.00365 or Ds>MaxD or Ds<0.0){
    //cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl;
    //cout<<"Nblend[4]:  "<<s.Nblend[4]<<"\t Nblend[1]:  "<<s.Nblend[1]<<endl;
    //cout<<"Aks:  "<<ex.Aks<<"\t latitude:   "<<s.lat<<"\t longtide:  "<<s.lon<<endl; }


    for(int i=0;  i<M; ++i){
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i] + 5.0*log10(Ds*100.0) + Ai[i];
    if(s.Nblend[i]>=k){ s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));} }
    if(k==1){
    for(int i=0; i<M;  ++i){s.Ai[i]=Ai[i]; s.Map[i]=Map[i]; s.Mab[i]= Mab[i];}
    s.col=s.col+s.Ai[0]-s.Ai[1];}
    }///loop over the stars


    for(int i=0; i<M; ++i){
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.magb[i]= -2.5*log10(fabs(s.Fluxb[i]));    
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.009543 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; cin>>yye;} }
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0 or s.Rstar<0.0 or s.mass>10000.0 or
    s.nums>Num or s.nums<=0 or Av<0.0 or s.cl<0 or s.cl>=6){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<"\t cl:  "<<s.cl<<endl;   cin>>yye;}
   // cout<<"************** End of func_source  ****************"<<endl;
   s.fb= s.blend[4]; 
   s.mbs= s.magb[4]; 
   //cout<<"blending:  "<<s.fb<<"\t m_base:  "<<s.mbs<<endl;
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s)
{
    double f,test;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    int il, yye;
    double mmin=20.0;
    double mmax=50.0;

    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}


    do{
    l.numl=(int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test  =((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){
    cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;  int ue; cin>>ue;}
    }while(test>rholens[l.numl]);
    l.Dl=double(l.numl*step);



   double randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<=s.rho_disk[l.numl]) l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
  else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}


    //do{
    //l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(mmax - mmin)) + mmin;
    //test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
    //if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
    //if(l.Ml> 1.0) f=pow(l.Ml,-3.0);
    //}while(test>f);
    l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(mmax - mmin)) + mmin;
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ros=double(s.Rstar*Rsun*l.xls/l.RE);
    l.mul=l.Vt*1000.0*3600.0*24.0/(l.Dl*AU);///mas/days
    l.u0=RandR(0.0,1.0);
    if(l.u0==0.0) l.u0=1.0e-5; 
    l.tetE= double(l.RE/AU/l.Dl);///marcs
    l.t0=RandR(2.0,5.0*year-2.0);
    //l.pt1=0.0;
    //l.pt2=5.0*year; //+2.5*l.tE;
  
    //if(s.ros<=0.0 or l.tE<=0.0 or l.Dl>s.Ds or l.Vt<=0.0 or l.Ml<0.0){
  //  cout<<"ERROR ros:  "<<s.ros<<endl;
 //   cout<<"Vt:  "<<l.Vt<<endl;
 //   cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
 //   cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<endl;
 //   cout<<"numl:  "<<l.numl<<"\t Vt:  "<<l.Vt<<"\t mul:  "<<l.mul<<endl; // cin>>yye;}
 //   cout<<"ksi:   "<<s.xi<<"\t ros:  "<<s.ros<<"\t l.tetE:  "<<l.tetE<<endl;
 //   cout<<"************** End of func_Lens  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    //mass, Teff, Age, logL,  log(g),  Z,  Rs,  MB, MV, MI, MK, Cl, type (13)
    int yye, uui, k1, k2, h, g;
    double metal, age,gravity, MB;
    char filename[40];
    FILE *fp2;


    double Age1[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm1[YZ]={0.0};
    double Age2[YZ]={0.0}; double B2[YZ]={0.0};  double M2[YZ]={0.0};   double mm2[YZ]={0.0};
    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0};
    FILE *meta;
    meta=fopen("./files/CMD_WFIRST/metal.txt","r");
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf   %d  %d\n",&Metal[i],&count[i],&number[i]);
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta);
    FILE *hks;
    hks=fopen("./files/CMD_WFIRST/HKS.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(hks,"%lf  %lf  %lf  %lf\n",&Age1[i],&mm1[i],&B1[i],&M1[i]);
    if(Age1[i]<0.0 or mm1[i]<0.0 or fabs(B1[i])>0.3 or M1[i]<0.5 or Age1[i]>18.0){
    cout<<"ERROR Age(HKS): "<<Age1[i]<<"\t metal: "<<mm1[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl; cin>>uui; }}
    fclose(hks);
    FILE *ji;
    ji=fopen("./files/CMD_WFIRST/JI.txt", "r");
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age2[i],&mm2[i],&B2[i],&M2[i]);
    if(Age2[i]<0.0 or mm2[i]<0.0 or fabs(B2[i])>1.7 or M2[i]<0.5 or Age2[i]>18.0  or Age1[i]!=Age2[i] or mm1[i]!=mm2[i]){
    cout<<"ERROR Age(JI): "<<Age2[i]<<"\t metal: "<<mm2[i]<<"\t B[i]"<<B2[i]<<"\t M[i]: "<<M2[i]<<"\t i: "<<i<<endl;
    cout<<"Age1[i]:  "<<Age1[i]<<"\t mm1[i]:  "<<mm1[i]<<endl;  cin>>uui;}}
    fclose(ji);

////=================================== THIN DISK ==============================
    int j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','i','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTiW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&gravity,&metal,&cm.Rs_d[j],&MB,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.cl_d[j],&cm.type_d[j]);
    ///*******************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0])         h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*cm.Mab_d[2][j]*1.05263157894737); ///H-band  Mks/Mk=1.05263157894737
    cm.Mab_d[4][j]= double(B2[g]+M2[g]*cm.Mab_d[1][j]);   ///J-band
    if(fabs(cm.Mab_d[3][j]-cm.Mab_d[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_d[4][j]-cm.Mab_d[1][j])>1.5){
    cout<<"ERROR: Mab_d(y-band): "<<cm.Mab_d[4][j]<<"\t Mab_d(z-band): "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
    ///********************************************************
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.Teff_d[j]<0.0 or metal>0.1 or age>10 or
    int(cm.cl_d[j])==6 or float(cm.type_d[j])>8.0 or cm.type_d[j]<1.0){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;






////=================================== BULGE ==================================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','b','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDbW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&gravity,&metal,&cm.Rs_b[j],&MB,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.cl_b[j],&cm.type_b[j]);
    ///*****************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl; cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*cm.Mab_b[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_b[4][j]= double(B2[g]+M2[g]*cm.Mab_b[1][j]);   ///J-band
    if(fabs(cm.Mab_b[3][j]-cm.Mab_b[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_b[4][j]-cm.Mab_b[1][j])>1.5){
    cout<<"ERROR: Mab_b(y-band): "<<cm.Mab_b[4][j]<<"\t Mab_b(z-band): "<<cm.Mab_b[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///*****************************************************
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or age>10 or metal>0.9 or cm.cl_b[j]==6 or
    cm.type_b[j]>=8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) ){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;  cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;






////=================================== THICK DISK =============================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c%c.dat",'C','M','D','T','k','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTkW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&gravity,&metal,&cm.Rs_t[j],&MB,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.cl_t[j],&cm.type_t[j]);
    ///*********************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*cm.Mab_t[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_t[4][j]= double(B2[g]+M2[g]*cm.Mab_t[1][j]);   ///J-band
    if(fabs(cm.Mab_t[3][j]-cm.Mab_t[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_t[4][j]-cm.Mab_t[1][j])>1.5){
    cout<<"ERROR: Mab_t(y-band): "<<cm.Mab_t[4][j]<<"\t Mab_t(z-band): "<<cm.Mab_t[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///********************************************************
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or metal>0.2||cm.cl_t[j]==6|| cm.type_t[j]>=8.0 or
    (cm.cl_t[j]==5 and int(cm.type_t[j])>8) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;








////=================================== STELLAR HALO ===========================
    j=0;
    sprintf(filename,"./files/CMD_WFIRST/%c%c%c%c%c.dat",'C','M','D','h','W');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDhW.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&gravity,&metal,&cm.Rs_h[j],&MB,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.cl_h[j],&cm.type_h[j]);
    ///************************************************************
    h=-1;
    if(metal<Metal[0] || metal==Metal[0]) h=0;
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69) {cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;}
    g=-1;    k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm1[k1]!=mm1[k2-1] or number[h]==0){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl;
    cout<<"age: "<<age<<"\t Age1[k1]: "<<Age1[k1]<<"\t Age1[k2-1]: "<<Age1[k2-1]<<endl;
    cout<<"metal[k1]: "<<mm1[k1]<<"\t metal[k2-1]: "<<mm1[k2-1]<<endl;  cin>>uui;}
    if(age<Age1[k1]  or age==Age1[k1])  g=k1;
    else if(age>Age1[k2-1] or age==Age1[k2-1]) g=int(k2-1);
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age1[k-1]>Age1[k] or mm1[k-1]!=mm1[k]){
    cout<<"Bad error: Age1[k-1]: "<<Age1[k-1]<<"\t Age1[k]: "<<Age1[k]<<"\t mm[k-1]:"<<mm1[k-1]<<"\t mm[k]"<<mm1[k]<<endl;  cin>>uui;}
    if((age-Age1[k-1])*(age-Age1[k])<0.0  ||  age==Age1[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;}
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*cm.Mab_h[2][j]*1.05263157894737); ///H-band  ?=Mks/Mk
    cm.Mab_h[4][j]= double(B2[g]+M2[g]*cm.Mab_h[1][j]);   ///J-band
    if(fabs(cm.Mab_h[3][j]-cm.Mab_h[2][j])>1.5 or fabs(age-Age1[g])>3.0 or fabs(cm.Mab_h[4][j]-cm.Mab_h[1][j])>1.5){
    cout<<"ERROR: Mab_h(y-band): "<<cm.Mab_h[4][j]<<"\t Mab_h(z-band): "<<cm.Mab_h[3][j]<<"\t metal: "<<metal<<endl;
    cout<<"age: "<<age<<"\t Age1[g]: "<<Age1[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui; }
    ///***********************************************************
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<0 or cm.cl_h[j]<0  or cm.cl_h[j]==6  or  cm.Teff_h[j]<0.0 or
    metal>0.1 || cm.cl_h[j]>7|| cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>8) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
   cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
   cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;


}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars
   double fh=1.0;///No limitation
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
/*
I assume that the up-limit of the mass is indicated by the simulation. Because it depends on the structre, .... all details can be found in mass_averaged.cpp  code. */

   /*
    char filename[40];
    FILE *fill;
    sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
    fill=fopen(filename,"w");
    if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}
   */



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = R_sun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/R_sun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/R_sun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3]

s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 // fclose(fill);
}
///===========================================================================//
double RandN(double sigma, double nnd){
   double rr,f,frand;
   do{
   rr=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nnd; ///[-N sigma:N sigma]
   f= exp(-0.5*rr*rr/sigma/sigma);
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return rr;
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex){
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     //cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }}
     fclose(fpd);}
     //cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<flag<<endl;
     return(flag);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
  if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk,   sigma_T_Disk,  sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
 test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(R_sun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(R_sun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= double( pi-fabs(tetd) );
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    

    if(vls1>0.0 and vls2>=0.0)         s.xi=atan(fabs(vls2)/fabs(vls1));
    else if(vls1<=0.0 and vls2>0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + M_PI/2.0;
    else if(vls1<0.0 and vls2<=0.0)    s.xi=atan(fabs(vls2)/fabs(vls1)) + M_PI;
    else if(vls1>=0.0 and vls2<0.0)    s.xi=atan(fabs(vls1)/fabs(vls2)) + 3.0*M_PI/2.0;
    else if(s.xi>2.0*M_PI)             s.xi-= 2.0*M_PI; 
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}

    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
    cout<<"source type:  "<<s.type<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    //cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
