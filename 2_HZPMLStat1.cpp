#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <sys/timeb.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double pi= M_PI; 
const double step=MaxD/(double)Num/1.0;///step in kpc
const double Hp= 6.62607004*pow(10.0,-34); 
const double KP=3.08568025*pow(10.,19); // in meter.
const double G=6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double Msun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double Mearth= 5.972*pow(10.0,24.0);// kg
const double Mjupiter=1.898*pow(10,27.0);// kg
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=364.5; 
const double binary_fraction=double(2.0/3.0);
const double Avks=double(8.20922); 
///============================ Besancon constant ==============///
const double Dsun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]={3.1,2.5,3.1,3.1};

///============================ WFIRST & OGLE &  KMTNet ===================
const int M=5;///number of filters  VIKH,W149
const double satu[M]={12.0, 12.0, 13.0, 13.0, 14.8}; //it should be changed
const double thre[M]={20.0, 21.0, 21.0, 21.0, 26.0};
const double FWHM[M]={0.91, 0.91, 0.91, 0.91 , 0.33};//3*pixel_size (0.11") of WFIRST, VIKH W149 Z087
const double AlAv[M]={1.009,0.600,0.118,0.184,0.225};///From besancon model[VIKH W149]
const double sigma[M]={0.022,0.022,0.02,0.025,0.025};//MOAاستفاده از مقاله کاردلی
const double Akv=0.118;
const double cade=15.16/(60.0*24.0);///W149_cadence in  day
const double Tog=250.0;
const int Nog=988;
const int Nkm=1253; 
////=======================================================
///**************  Number constant  ********************///
const int nte=39;//// number of tE range 
const int Nw=int(123);///sigma_WFIRST.dat
const int YZ=3615;///number of rows in file yzma
const int N1=36224, N2=25000, N3=3818, N4=3500;///CMD_BESANCON, thinD, bulge, thickD, halo
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI, lat, lon;
    double od_disk,od_ThD,od_bulge,od_halo,opt;///optical  depth
    double od_dlmc, od_blmc, od_hlmc; 
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double rho_hlmc[Num], rho_dlmc[Num],rho_blmc[Num]; 
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart;
    double Nblend[M], blend[M], blendl[M], Fluxb[M], magb[M], Ai[M], Mab[M], Map[M]; 
    double type, Tstar, logl, col, Rstar, mass,vs;
    double Romax,ro_star; 
};
struct lens{
    int numl,struc, cl;
    double mass, Ml, Dl, vl , vs, Vt, xls;
    double rhomaxl,tE,RE;
    double mul, u0; 
    double pt1, pt2, dt, t0w, t0o; 
    double ksi, dis, q;
    double Tstar, logl, col,  Rstar, type, semi; 
    double Mab[M], Map[M], theta[M], Mp, Fl[M]; 
    int HZP[2];
    
};
struct dete { 
   double Eff[2]; ///OHZ,CHZ
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
///===================== FUNCTION ===============================================================
int Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void func_source(source & s, CMD & cm , extinc & ex);
void func_lens(lens & l, source & s, CMD & cm,  extinc & ex);
void vrel(source & s , lens & l);
void Disk_model(source & s);
void optical_depth(source & s);
double Interpol(double ds, extinc & ex);
double RandN(double sigma, double Nn);
double RandR(double down, double up);
double ErrorCal(double mag, int sur);
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

     source s;
     lens l;
     dete d; 
     CMD cm;  
     extinc ex;
     read_cmd(cm);  
     
     VBBinaryLensing vbb;
     vbb.Tol=1.e-3;
     vbb.LoadESPLTable("./files/ESPL.tbl");

///===========================================================   
    double magw[Nw], errw[Nw];
    FILE * wfirstf;  
    wfirstf=fopen("./files/sigma_WFIRST.txt","r");
    for(int i=0; i<Nw;++i){
    fscanf(wfirstf,"%lf  %lf\n",&magw[i],&errw[i]);}
    fclose(wfirstf);  
    cout<<"*****  sigma_wfirst was read *****  "<<endl;  
///===========================================================  
   FILE* oglef;  FILE* kmtnet; 
   double dto[Nog],  dtk[Nkm], dmo, mags; 
   oglef=fopen("./files/OGLE_delta.dat","r");
   for (int i=0; i<Nog; ++i){ fscanf(oglef,"%lf   %lf   %lf \n", &dto[i], &dmo, &mags); }
   kmtnet=fopen("./files/KMTNet_delta_2.dat","r");
   for (int i=0; i<Nkm; ++i){ fscanf(kmtnet,"%lf   %lf   %lf \n",&dtk[i], &dmo, &mags); }
   fclose(oglef);  fclose(kmtnet);
   cout<<"The OGLE & KMTNet  files were read *************************"<<endl; 
///===========================================================


    FILE* magg;  FILE* data2;    FILE* data3;   //FILE *caustic;  
    FILE* fil1;  //FILE* fil2;    
    fil1=fopen("./files/HZPparc1.txt","w");
   // fil2=fopen("./files/HZPEffi.txt","w");
    fclose(fil1);// fclose(fil2);
    
    double obs[6][2]={0.0};
    obs[0][0]=0.0*year+1.0;    obs[0][1]=0.0*year +72.0;
    obs[1][0]=0.0*year+182.0;  obs[1][1]=0.0*year+254.0;
    obs[2][0]=1.0*year+1.0;    obs[2][1]=1.0*year+72.0;
    obs[3][0]=3.0*year+182.0;  obs[3][1]=3.0*year+254.0;
    obs[4][0]=4.0*year+1.0;    obs[4][1]=4.0*year+72.0;
    obs[5][0]=4.0*year+182.0;  obs[5][1]=4.0*year+254.0;
    
    int status;
    int flad, flagf, nog, nkm, flag[2];
    int ogle, wfirst, ndw, ndo, flago, w1;
    char filnam1[40], filnam2[40], filnam3[40], filnam4[40];
    double xlens, ylens, xlens2, magnio, dism, erro, test,  uu2; 
    double Astar, Astar1, magni[M], magni1[M];
    double timp1, timpo, timpk, lonn, diss,  chi, chi1, chi2;  
    double shi, shi1, shi2, flag0, flag1, flag2, Flag0, Flag1, Flag2;
    double glag0, glag1, glag2, Glag0, Glag1, Glag2;  
    double t_0, t_1, T_0, T_1, deltaA, tt;
    double xcau[800], ycau[800], tcau; 
    double dchiP, dchiL, dshiP, dshiL; 
    double effic[2], Effic[2], Nmlw=0.0, Nmlo=0.0, Npw=0.0, Npo=0.0, NHZw=0.0, NHZo=0.0; 
    int flag_peak=0, Flag_peak=0, flag_det=0, Flag_det=0, flag_mlw=0, flag_mlo=0, glag_det=0, Glag_det=0;
    double close; 
    int typeb=0;
    int flagpw, flagpo;  
    double lml, dd1, dd2, stt, lfl, Npdet=0.0000000000085646465;
 
     
     d.Eff[0]=d.Eff[1]=0.0;  
 
 
 
     s.lat=-1.5;//RandR(,-4.0);
     lonn= 1.0 ;//RandR(+3.0,+6.0);
     //cout<<"lat:  "<<s.lat<<"\t lonn:  "<<lonn<<endl;
     //cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
     if(lonn<=0.0)   s.lon=360.0+lonn;
     else            s.lon=lonn;
     s.TET=(360.0-s.lon)/RA;///radian
     s.FI=s.lat/RA;///radian
     Disk_model(s);
     if(int(Extinction(ex,s))==1){

 
 
     for(int icon=1; icon<100000000; ++icon){ 
     cout<<"***** Step:  "<<icon<<endl;
     
    

     do{
     func_source(s, cm, ex);
     func_lens(l, s, cm, ex);     
     }while(l.tE<=0.1 or l.tE>300.0); 
     optical_depth(s);
     
     
     ogle=wfirst=0;
     //test= (double)rand()/((double)(RAND_MAX)+(double)(1.0));
     //if(test<=s.blendl[1] and s.magb[1]<=thre[1]) ogle=1; 
     test= (double)rand()/((double)(RAND_MAX)+(double)(1.0));
     if(test<=s.blendl[4] and s.magb[4]<=thre[4]) wfirst=1;

     
     
     if(wfirst>0){
     
     flagf=0;
     if(icon%100000000==0){  
     flagf=1;
     sprintf(filnam1,"./files/MONT/%c%c%c%c%d.dat",'m','a','g','_',icon);
     magg=fopen(filnam1,"w");
     //if(ogle>0){
     //sprintf(filnam2,"./files/MONT/%c%c%d%c%d.dat",'d','a',1,'_',icon);
     //data2=fopen(filnam2,"w");}
     if(wfirst>0){
     sprintf(filnam3,"./files/MONT/%c%c%d%c%d.dat",'d','a',0,'_',icon);
     data3=fopen(filnam3,"w");}}
     l.dt=double(10.0/60.0/24.0);
     
   /*vbb.PrintCau(l.dis,l.q,icon);
     cout<<"Caustic.txt file is made***************"<<endl; 
     sprintf(filnam4,"./files/MONT/%c%c%c%c%d.dat",'c','a','u','_',icon);
     caustic=fopen(filnam4,"r");
     if(!caustic){cout<<"File is not exit!!!! "<<"\t icoc:  "<<icon<<endl;}
     for(int i=0; i<800; ++i){
     fscanf(caustic, "%lf    %lf %d\n", &xcau[i],&ycau[i],&tcau);
     if(tcau==0)  typeb=1; //intermediate  1
     if(tcau==1)  typeb=2; //wide   2
     if(tcau==2)  typeb=3; //close   3
     }
     fclose(caustic); 
     status=1; 
     if(flagf<1) status = remove(filnam4);
     if(flagf>0 and status==0){
     cout<<"\nErorr file is removed mistakenly!"<<filnam4<<endl;;  int uue;  cin>>uue;} 
     cout<<"Here is end of caustic.txt file"<<endl; */
 
    
     ndw=0, ndo=0; 
     flag0=flag1=flag2=Flag0=Flag1=Flag2=0.0; 
     glag0=glag1=glag2=Glag0=Glag1=Glag2=0.0; 
     flag_peak=0, Flag_peak=0, flag_det=0, Flag_det=0, flag_mlw=0, flag_mlo=0, glag_det=0, Glag_det=0; 
     t_1=t_0=T_1=T_0=0.0; close=10.0; 
     chi=chi1=chi2=0.0; 
     shi=shi1=shi2=0.0; 
     timp1=timpo=timpk=0.0;
     nog=int((double)rand()*(Nog-5)/(double)(RAND_MAX+1.))+1;
     nkm=int((double)rand()*(Nkm-5)/(double)(RAND_MAX+1.))+1;
     
     for(double tim=l.pt1; tim<l.pt2; tim=tim+l.dt){
     xlens = tim/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
     ylens = tim/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
     xlens2= tim/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi) + fabs(l.dis*l.q/(1.0+l.q));
     uu2=sqrt(xlens2*xlens2+ylens*ylens); 
     
     
     /*dism=1000000.0; 
     for(int j=0; j<800; ++j){ 
     diss=sqrt((xlens-xcau[j])*(xlens-xcau[j])+(ylens-ycau[j])*(ylens-ycau[j])); 
     if(diss<dism)  dism=diss;}
     dism=double(dism/s.ro_star);     
     if(dism<close) close=dism;   */
     
     
     Astar1= vbb.ESPLMag2(uu2,s.ro_star);
     Astar=vbb.BinaryMag2(l.dis,l.q, xlens, ylens, s.ro_star);
     //cout<<"Astar1:   "<<Astar1<<"\t Astar:    "<<Astar<<endl;
     //cout<<"dis:  "<<l.dis<<"\t l.q:  "<<l.q<<"\t xlens:  "<<xlens<<"\t ylens:   "<<ylens<<endl;
     
     for(int i=0;i<M; ++i){
     magni1[i]=s.magb[i]-2.5*log10(Astar1*s.blendl[i] + 1.0 - s.blendl[i]);
     magni[i] =s.magb[i]-2.5*log10(Astar* s.blendl[i] + 1.0 - s.blendl[i]);}
     if(flagf>0) 
     fprintf(magg,"%.5lf  %.5lf  %.5lf  %.5lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf\n",
     tim,magni[4], magni1[4],magni[1], magni1[1],Astar,Astar1,xlens,ylens); 
     
     
     timp1 += l.dt;
     timpo += l.dt; 
     timpk += l.dt; 
///************************** WFIRST  ********************************************************************   
    if(wfirst>0){
    flago=0; 
    for(int ii=0; ii<6; ++ii){ 
    if((obs[ii][0]- tim-l.t0w)*(obs[ii][1]- tim-l.t0w)<=0.0) {flago=1; break;}}
    if(timp1>cade or timp1==cade){
    timp1=double(timp1-cade);
    if(magni[4]>=satu[4] and  magni[4]<=thre[4] and flago>0){
    w1=-1;
    if(magni[4]<magw[0] or magni[4]==magw[0]) w1=0;
    else if(magni[4]>magw[Nw-1] or magni[4]==magw[Nw-1]) w1=Nw-1;
    else{
    for(int i=1; i<Nw; ++i){
    if(double((magni[4]-magw[i])*(magni[4]-magw[i-1]))<0.0 or  magni[4]==magw[i-1]){w1=i-1; break;}}}
    if(w1==-1){cout<<"mapt[0]: "<<magni[4]<<"\t w: "<<w1<<endl;  int ww; cin>>ww;}
    
    magnio=magni[4] + RandN(errw[w1] ,3.0); 
    
    chi  +=fabs( (magnio- s.magb[4])*(magnio- s.magb[4])/(errw[w1]*errw[w1]));
    chi1 +=fabs( (magnio- magni1[4])*(magnio- magni1[4])/(errw[w1]*errw[w1]));
    chi2 +=fabs( (magnio-  magni[4])*(magnio-  magni[4])/(errw[w1]*errw[w1]));
    if(flagf>0) fprintf(data3,"%.5lf   %.7lf   %.7lf  %d\n",tim, magnio, errw[w1], 1);
    
    t_1=tim; 
    if(double(t_1 * t_0)<=0.0 and ndw>0 and fabs(t_0-t_1)<5.0*cade)  flag_peak=1;
    t_0=t_1; 
    
    
    flag2=0.0;
    if(fabs(magnio-s.magb[4])>fabs(4.0*errw[w1])) flag2=1.0;
    if(ndw>1 and float(flag0+flag1+flag2)>2.0)    flag_det=1;
    flag0=flag1;
    flag1=flag2;
    
    glag2=0.0; 
    if(fabs(magnio-magni1[4])>fabs(4.0*errw[w1])) glag2=1.0;
    if(ndw>1 and float(glag0+glag1+glag2)>2.0)    glag_det=1;
    glag0=glag1;
    glag1=glag2;
    
    ndw+=1;}}}
///************************** OGLE  ********************************************************************
    /*if(ogle>0){   
    flago=0; 
    if((tim+l.t0o)>0.0 and (tim+l.t0o)<=Tog) flago=1;
    flag[0]=flag[1]=-1;
    if(timpo>dto[nog]){//OGLE data
    timpo -= dto[nog];
    flag[0]=1; 
    nog=int((double)rand()*(Nog-3)/(double)(RAND_MAX+1.))+1;}
    if(timpk>dtk[nkm]){//KMTNet data
    timpk -= dtk[nkm];
    flag[1]=1; 
    nkm=int((double)rand()*(Nkm-3)/(double)(RAND_MAX+1.))+1;}
    for(int i=0; i<2; ++i){ 
    if(flag[i]>0 and magni[1]>=satu[1] and magni[1]<=thre[1] and flago>0){  
    erro= ErrorCal(magni[1],i); 
    magnio= magni[1] + RandN(erro,3.0);
    shi  +=fabs( (magnio- s.magb[1])*(magnio- s.magb[1])/(erro*erro));
    shi1 +=fabs( (magnio- magni1[1])*(magnio- magni1[1])/(erro*erro));
    shi2 +=fabs( (magnio-  magni[1])*(magnio-  magni[1])/(erro*erro));
    if(flagf>0) fprintf(data2,"%.5lf   %.7lf  %.7lf   %d\n",tim, magnio,erro,i);
    T_1=tim; 
    if(double(T_1*T_0)<=0.0 and ndo>0 and fabs(T_0-T_1)<1.5)  Flag_peak=1;
    T_0=T_1; 
    Flag2=0.0; 
    if(fabs(magnio-s.magb[1])>fabs(3.0*erro))     Flag2=1.0;
    if(ndo>1 and float(Flag0+Flag1+Flag2)>2.0)    Flag_det=1;
    Flag0=Flag1;
    Flag1=Flag2;
    Glag2=0.0; 
    if(fabs(magnio-magni1[1])>fabs(4.0*erro))     Glag2=1.0;
    if(ndo>1 and float(Glag0+Glag1+Glag2)>2.0)    Glag_det=1;
    Glag0=Glag1;
    Glag1=Glag2;
    ndo+=1;}}}*/
    }//loop time
    if(flagf>0){ fclose(magg);
    //if(ogle>0)   fclose(data2);  
    if(wfirst>0) fclose(data3); }
     cout<<"End of loop time ******************"<<endl;
    ///************ WFIRST *********************************************************************
    dchiP=fabs(chi2-chi1);   
    dchiL=fabs(chi2-chi );
    if(flag_peak>0 and dchiL>500.0 and flag_det>0 and ndw>5){Nmlw+=1.0;   flag_mlw=1;} 
    else flag_mlw=0;
    
    if(flag_mlw>0){
    if(dchiP>300.0 and glag_det>0)  flagpw=1;
    else flagpw=0;  }
   
  
   
    if(flag_mlw>0 and flagpw>0){
    Npdet+=1.0; 
    if(l.HZP[0]>0)  d.Eff[0]+=1.0; ///CHZ
    if(l.HZP[1]>0)  d.Eff[1]+=1.0; ///OHZ
    }
    if(flag_mlw>0){
    fil1=fopen("./files/HZPparc1.txt","a+");
    fprintf(fil1,
   "%d %.3lf  %.3lf  " ///3
   "%d  %d  %.5lf   %.4lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf   %.4lf  %.4lf  "  //17
   "%d  %d  %.5lf   %.4lf  %.4lf  %.4lf  %.4lf  %.2lf  %.3lf  %.4lf  %.4lf   %.4lf   %.4lf  %.4lf  "   //31
   "%.5lf  %.5lf  %.6lf  %.6lf  %.6lf  %.6lf  %.2lf  %.2lf  "  //39
   "%.7lf  %.6lf  %.5lf  %.5lf  %.3lf   %.5lf  %.5lf  %.6lf  %.4lf  %.5lf  %.6lf  %.7lf   " ///51
   "%d  %d  %d  "//54
   "%d  %d  %.1lf  %d  %.1lf  %d  %.7lf  %d  %.4lf   %d  %d  %d  %.8lf   %.8lf\n", //68
   icon,s.lat, s.lon, //3
   l.struc, l.cl, l.mass, l.Dl, l.Tstar, l.Rstar, l.logl, l.type, l.col, l.vl, l.Mab[1], l.Mab[4], l.Map[1], l.Map[4], //17
   s.struc, s.cl, s.mass, s.Ds, s.Tstar, s.Rstar, s.logl, s.type, s.col, s.vs, s.Mab[1], s.Mab[4], s.Map[1], s.Map[4], //31
   s.magb[1], s.magb[4], s.blend[1], s.blendl[1], s.blend[4], s.blendl[4], s.Nblend[1], s.Nblend[4], //39
   l.q*1.0e4, l.dis, l.tE, l.RE/AU, l.t0w, l.mul, l.Vt, l.u0, l.ksi, l.semi, s.opt*1.0e6, s.ro_star,//51
   l.theta[4], wfirst, flagf, //54
   flag_det, flag_peak, dchiL, flag_mlw, dchiP, ndw, l.Fl[4], glag_det, close, typeb,///64 
   l.HZP[0], l.HZP[1],  double(d.Eff[0]*100.0/Npdet), double(d.Eff[1]*100.0/Npdet));///68
   fclose(fil1); 
   cout<<"************* End of saving in the file ***********************"<<endl;}
   if(int(icon)%1==0){
   cout<<"============================================================="<<endl;
   cout<<"icon    "<<icon<<endl;
   cout<<"lat:  "<<s.lat<<"\t lon:  "<<s.lon<<endl;
   cout<<"********************** SOURCE **************************"<<endl; 
   cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<"\t strucs:  "<<s.struc<<endl;
   cout<<"cl:  "<<s.cl<<"\t mass:  "<<s.mass<<"\t Tstar:  "<<s.Tstar<<endl;
   cout<<"Rstar:  "<<s.Rstar<<"\t logl:  "<<s.logl<<"\t type:  "<<s.type<<endl;
   cout<<"col:  "<<s.col<<"\t vs:  "<<s.vs<<"\t Mag_I:  "<<s.Mab[1]<<"\t Mag_W149:  "<<s.Mab[4]<<endl;
   cout<<"********************** LENS **************************"<<endl; 
   cout<<"Dl:  "<<l.Dl<<"\t numl:  "<<l.numl<<"\t strucl:  "<<l.struc<<endl;
   cout<<"cl:  "<<l.cl<<"\t mass:  "<<l.mass<<"\t Tstar:  "<<l.Tstar<<endl;
   cout<<"Rstar:  "<<l.Rstar<<"\t logl:  "<<l.logl<<"\t type:  "<<l.type<<endl;
   cout<<"col:  "<<l.col<<"\t vs:  "<<l.vl<<"\t Mag_I:  "<<l.Mab[1]<<"\t Mag_W149:  "<<l.Mab[4]<<endl;
   cout<<"theta[1]:  "<<l.theta[1]<<"\t theta[4]:    "<<l.theta[4]<<endl;
   cout<<"*********************** LENSING ***********************"<<endl;
   cout<<"q:  "<<l.q<<"\t dis:  "<<l.dis<<"\t semi:  "<<l.semi<<endl; 
   cout<<"tE:  "<<l.tE<<"\t RE(AU):  "<<l.RE/AU<<"\t t0w:  "<<l.t0w<<"\t t0o:  "<<l.t0o<<endl;
   cout<<"Vt:  "<<l.Vt<<"\t mul:  "<<l.mul<<"\t u0:  "<<l.u0<<"\t ksi:  "<<l.ksi*180.0/pi<<endl;
   cout<<"************ WFIRST & OGLE ***************************"<<endl;
   cout<<"flag_mlw:  "<<flag_mlw<<"\t dchiL:  "<<dchiL<<"\t dchiP:  "<<dchiP<<endl;
   cout<<"flag_mlo:  "<<flag_mlo<<"\t dshiL:  "<<dshiL<<"\t dshiP:  "<<dshiP<<endl;  
   cout<<"effip_WFIRST:  "<<effic[0]<<"\t effip_OGLE/KMTNET:   "<<Effic[0]<<endl; 
   cout<<"effiHZ_WFIRST:  "<<effic[1]<<"\t effiHZ_OGLE/KMTNET:   "<<Effic[1]<<endl; 
   cout<<"Nmlw:    "<<Nmlw<<"\t   Nmlo:    "<<Nmlo<<endl;
   cout<<"HZP[0]:  "<<l.HZP[0]<<"\t HZP[1]:   "<<l.HZP[1]<<endl;
   cout<<"Fl[4]:   "<<l.Fl[4]<<"\t Npdet:  "<<Npdet<<endl;
   cout<<"Eff_OHZ:  "<<double(d.Eff[0]*100.0/Npdet)<<"\t Eff_CHZ:   "<<double(d.Eff[1]*100.0/Npdet)<<endl;
   cout<<"==============================================================="<<endl;} 
   
   }//if ogle>0 or wfirst>0
   }//end of EXtinction
   }//loop icon  
   return(0);
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
double ErrorCal(double mag, int sur){
   double error; 
   if(sur==0){ ///OGLE
   if(mag<15.0 ) error=0.003;   //8.62339038e-03 6.82867290e-03 2.27422710e-03 1.66494077e+01
   if(mag>=15.0) error=0.00862339038 + 0.00682867290*(mag-16.6494077) +  0.00227422710*(mag-16.6494077)*(mag-16.6494077);  
   if(mag<20.0 and  error >0.095) {cout<<"Error(OGLE) error:  "<<error<<"\t mag:  "<<mag<<endl;  int uue; cin>>uue;}}
   if(sur==1){//KMTNet
   if(mag<14.0 ) error= 0.0038;//7.36447617e-02 -1.61925360e-02  9.52245559e-04  5.61612331e+00
   if(mag>=14.0) error= 0.0736447617 -0.0161925360*(mag-5.61612331) +  0.000952245559*(mag-5.61612331)*(mag-5.61612331); 
   if(mag<20.0 and error>0.095) {cout<<"Error(KMTNet)  error:  "<<error<<"\t mag:  "<<mag<<endl;   int yye;  cin>>yye;}}
   return(error);
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///==============================================================//6
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*Msun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opt=0.0;
    s.od_dlmc=s.od_hlmc=s.od_blmc=0.0; 
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;
    s.od_dlmc +=  s.rho_dlmc[k]*x*(1.0-x)*dx*CC;
    s.od_blmc +=  s.rho_blmc[k]*x*(1.0-x)*dx*CC; 
    s.od_hlmc +=  s.rho_hlmc[k]*x*(1.0-x)*dx*CC; }
    s.opt= fabs(s.od_disk+s.od_ThD+s.od_bulge+s.od_halo+ s.od_dlmc + s.od_blmc + s.od_hlmc );///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
   cout<<"End of optical depth function **********************"<<endl;
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
   cout<<"Nblend[4]:  "<<s.Nblend[4]<<"\t maxnb: "<<maxnb<<endl;
   

    for(int k=1; k<=int(maxnb); ++k){
    
    do{
    num=int((Num-15.00)*(double)rand()/(double)(RAND_MAX+1.) +10.00);
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    Ds=(double)(num*step);
    }while(rho>s.Rostari[num] or Ds<0.1 or Ds>MaxD);///distance larger than 50.
    if(Ds>MaxD or Ds<0.1){cout<<"ERROR (1): Ds: "<<Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<num<<endl; cin>>yye;}
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}
    //cout<<"Ds:  "<<Ds<<endl;
    
    
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
    for(int i=0; i<5; ++i){ Mab[i]=cm.Mab_d[i][num]; }
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
    
    

    Mab[4]=(Mab[2]+Mab[3]+Mab[4])/3.0;//absolute magnitude in W149 (K+H+J)/3
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av<0.0)    Av=0.0;
    if(Av>6.0 or Av<-0.00365 or Ds>MaxD or Ds<0.0){cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl;   cin>>yye;}
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
    if(int(s.Nblend[i])<1.0 or (s.Nblend[i]==1.0 and s.blend[i]<1.0) or s.blend[i]>1.0 or s.blend[i]<0.0 or
    s.Fluxb[i]<-0.009543 or s.Ai[i]<0.0){
    cout<<"BIGG ERRROR nsbl: "<<s.Nblend[i]<<"\t Nblend[i]: "<<s.Nblend[i]<<"\t belnd: "<<s.blend[i]<<endl;
    cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<"\t No. filter:  "<<i<<"\t ext:  "<<s.Ai[i]<<endl; cin>>yye;} }
    if(s.type>8.0 or s.type<1.0 or s.mass<=0.0 or s.Tstar<0.0 or s.Rstar<0.0 or s.mass>10000.0 or 
    s.nums>Num or s.nums<=0 or Av<0.0 or s.cl<0 or s.cl>=6){
    cout<<"ERROR(source):  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<"\t cl:  "<<s.cl<<endl;   cin>>yye;}
    cout<<"************** End of func_source  ****************"<<endl;
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, CMD & cm,  extinc & ex)
{
    double f,test, Ai[M], dism, Av;
    double rholens[s.nums+2];
    l.rhomaxl=0.0;
    int il, yye;
    double x_rand, sign1, ksi, xrel, yrel,len1,len2, Gama, Beta, minm; 


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

  
  if(l.struc==0){///thin disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(5.0-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml> 1.0) f=pow(l.Ml,-3.0);
  }while(test>f);  
  minm=10000.0, il=-1;
  for(int i=0; i<N1; ++i){
  dism=fabs(cm.mass_d[i]-l.Ml);
  if(dism<minm and int(cm.cl_d[i])>3 and int(cm.cl_d[i])<7 ){il=i; minm=dism;}}
  for(int i=0; i<M; ++i) l.Mab[i]=cm.Mab_d[i][il];
  l.Mab[4]=(l.Mab[2]+l.Mab[3]+l.Mab[4])/3.0;
  l.cl=     cm.cl_d[il];  
  l.Rstar=  cm.Rs_d[il]; 
  l.type= cm.type_d[il];
  l.mass= cm.mass_d[il];
  l.Tstar=cm.Teff_d[il];
  l.logl= cm.logl_d[il];
  l.col=l.Mab[0]-l.Mab[1]; 
  if(l.mass<0.0 or l.Rstar<0.0 or l.type>8.0 or l.Tstar<0.0 or l.cl>6 or l.cl<4 or il<0 or fabs(l.mass-l.Ml)>5.0){
  cout<<"Error mass(thin disk): "<<l.mass<<"\t Tstar: "<<l.Tstar<<"\t l.cl: "<<l.cl<<"\t l.type:  "<<l.type<<endl;
  cout<<"il:  "<<il<<"\t diff:  "<<fabs(l.Ml-l.mass)<<"\t l.Ml:  "<<l.Ml<<endl; cin>>yye;} }
  l.Ml=l.mass; 


///**********************************************************************
  if(l.struc==1){///Galactic bulge
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(5.0-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*20.5);
  if(l.Ml<=0.7) f=pow(l.Ml,-1.0)*pow(0.7,-1.35);//Penney et al.  2019
  if(l.Ml>0.7 ) f=pow(l.Ml,-2.35);
  }while(test>f);  
  minm=10000.0, il=-1;
  for(int i=0; i<N2; ++i){
  dism=fabs(cm.mass_b[i]-l.Ml);
  if(dism<minm and int(cm.cl_b[i])>3 and int(cm.cl_b[i])<7){il=i; minm=dism;}}
  for(int i=0; i<M; ++i) l.Mab[i]=cm.Mab_b[i][il];
  l.Mab[4]=(l.Mab[2]+l.Mab[3]+l.Mab[4])/3.0;
  l.cl=     cm.cl_b[il];  
  l.Rstar=  cm.Rs_b[il]; 
  l.type= cm.type_b[il];
  l.mass= cm.mass_b[il];
  l.Tstar=cm.Teff_b[il];
  l.logl= cm.logl_b[il];
  l.col=l.Mab[0]-l.Mab[1]; 
  if(l.mass<0.0 or l.Rstar<0.0 or l.type>8.0 or l.Tstar<0.0 or l.cl>6 or l.cl<4 or il<0 or fabs(l.mass-l.Ml)>5.0){
  cout<<"Error mass(Galactic bulge): "<<l.mass<<"\t Tstar: "<<l.Tstar<<"\t l.cl: "<<l.cl<<"\t l.type:  "<<l.type<<endl;
  cout<<"il:  "<<il<<"\t diff:  "<<fabs(l.Ml-l.mass)<<"\t l.Ml:  "<<l.Ml<<endl;   cin>>yye;}}
  l.Ml=l.mass; 


///**********************************************************************
  if(l.struc==2){///thick disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(5.0-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*4.0);
  f=pow(l.Ml,-0.5);
  }while(test>f);  
  minm=10000.0, il=-1;
  for(int i=0; i<N3; ++i){
  dism=fabs(cm.mass_t[i]-l.Ml);
  if(dism<minm and int(cm.cl_t[i])>3 and int(cm.cl_t[i])<7){il=i; minm=dism;}}
  for(int i=0; i<M; ++i) l.Mab[i]=cm.Mab_t[i][il];
  l.Mab[4]=(l.Mab[2]+l.Mab[3]+l.Mab[4])/3.0;
  l.cl=     cm.cl_t[il];  
  l.Rstar=  cm.Rs_t[il]; 
  l.type= cm.type_t[il];
  l.mass= cm.mass_t[il];
  l.Tstar=cm.Teff_t[il];
  l.logl= cm.logl_t[il];
  l.col=l.Mab[0]-l.Mab[1]; 
  if(l.mass<0.0 or l.Rstar<0.0 or l.type>8.0 or l.Tstar<0.0 or l.cl>6 or l.cl<4 or il<0 or fabs(l.mass-l.Ml)>5.0){
  cout<<"Error mass(Galactic bulge): "<<l.mass<<"\t Tstar: "<<l.Tstar<<"\t l.cl: "<<l.cl<<"\t l.type:  "<<l.type<<endl;
  cout<<"il:  "<<il<<"\t diff:  "<<fabs(l.Ml-l.mass)<<"\t l.Ml:  "<<l.Ml<<endl;  cin>>yye;}}
  l.Ml= l.mass; 

///**********************************************************************
  if(l.struc==3){///stellar halo
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(5.0-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*4.0);
  f=pow(l.Ml,-0.5);
  }while(test>f);  
  minm=10000.0, il=-1;
  for(int i=0; i<N4; ++i){
  dism=fabs(cm.mass_h[i]-l.Ml);
  if(dism<minm and int(cm.cl_h[i])>3 and int(cm.cl_h[i])<7){il=i; minm=dism;}}
  for(int i=0; i<M; ++i) l.Mab[i]=cm.Mab_h[i][il];
  l.Mab[4]=(l.Mab[2]+l.Mab[3]+l.Mab[4])/3.0;
  l.cl=     cm.cl_h[il];  
  l.Rstar=  cm.Rs_h[il]; 
  l.type= cm.type_h[il];
  l.mass= cm.mass_h[il];
  l.Tstar=cm.Teff_h[il];
  l.logl= cm.logl_h[il];
  l.col=l.Mab[0]-l.Mab[1]; 
  if(l.mass<0.0 or l.Rstar<0.0 or l.type>8.0 or l.Tstar<0.0 or l.cl>6 or l.cl<4 or il<0 or fabs(l.mass-l.Ml)>5.0){
  cout<<"Error mass(Galactic bulge): "<<l.mass<<"\t Tstar: "<<l.Tstar<<"\t l.cl: "<<l.cl<<"\t l.type:  "<<l.type<<endl;
  cout<<"il:  "<<il<<"\t diff:  "<<fabs(l.Ml-l.mass)<<"\t l.Ml:  "<<l.Ml<<endl;  cin>>yye;}}
  l.Ml=l.mass;
 cout<<"l.Ml:  "<<l.Ml<<"\t l.mass:  "<<l.mass<<endl; 

///**********************************************************************
    ex.Aks=Interpol(l.Dl,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av<0.0)    Av=0.0;
    if(Av>6.0 or Av<-0.00365 or l.Dl>MaxD or l.Dl<0.0){cout<<"ERROR Dl:  "<<l.Dl<<" \t Av:  "<<Av<<endl; cin>>yye;}
    for(int i=0;  i<M;  ++i){    
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i], 1.5); //extinction in other bands
    if(Ai[i]<0.0) Ai[i]=0.0;
    l.Map[i]=l.Mab[i]+5.0*log10(l.Dl*100.0)+Ai[i];}
    l.col=l.col+Ai[0]-Ai[1]; 
    
    
     
    double qmin=double(0.00001);//in Earth mass
    double qmax=double(0.06);///0.08 Msun
    double q0=0.0002; 
    double hel, f1, f2, fac, dhz1, dhz2; 
    
do{  
    f1=pow(qmin/q0,-0.27); 
    f2=pow(qmax/q0,-1.93); 
    if(f1<f2){hel=f1;   f1=f2;  f2= hel;}
    do{ 
    l.q=(double)rand()/(double)(RAND_MAX+1.0) *(qmax - qmin) +qmin;
    test=fabs((double)rand()/(double)(RAND_MAX+1.)*(f1-f2) +f2);
    if(l.q<q0 )   f=pow(l.q/q0,-0.27);
    if(l.q>=q0)   f=pow(l.q/q0,-1.93);
    if(f>f1 or f<f2) {cout<<"Error F: "<<f<<"\t f1:  "<<f1<<"\t f2:  "<<f2<<endl;  int rre;  cin>>rre;}
    }while(test>f);  
    l.Mp= double(l.Ml*l.q*Msun/Mjupiter); /// in Jupiter mass
}while(l.Mp>15.0); 

    
    
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*(1.0+l.q)*Msun*s.Ds*KP)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///in day
    s.ro_star=s.Rstar*Rsun*l.xls/l.RE; 
    l.mul=l.Vt*1000.0*180.0*3600.0*1000.0*3600.0*24.0*364.5/(l.Dl*KP*M_PI);    
    l.t0w=RandR(1.0,5.0*year);
    l.t0o=(double)rand()/(double)(RAND_MAX+1.)*Tog; 
    l.pt1=-2.5*l.tE; 
    l.pt2=+2.5*l.tE; 
    l.ksi=RandR(0.0 , 359.0)*M_PI/180.0;
    l.u0=RandR(0.0,1.0);   
    

   
    Beta=fabs((double)rand()/(double)(RAND_MAX+1.)*89.0);
    sign1= fabs((double)rand()/(double)(RAND_MAX+1.));
    if(sign1>0.5) Beta= -1.0*Beta;
    Gama=fabs((double)rand()/(double)(RAND_MAX+1.)*89.0);
    sign1=fabs((double)rand()/(double)(RAND_MAX+1.));
    if(sign1<0.5) Gama=-1.0*Gama;
    Beta=Beta*M_PI/180.0;  
    Gama=Gama*M_PI/180.0;
    ksi=fabs((double)rand()/(double)(RAND_MAX+1.)*2.0*M_PI);
    do{
    x_rand= fabs((double)rand()/(double)(RAND_MAX+1.) * 99.99)+0.01;//AU
    f= pow(x_rand,-1.0);///dN/ds~1/s
    test= fabs((double)rand()/(double)(RAND_MAX+1.)*100.0);
    }while(test>f);
    l.semi=x_rand;
    xrel=l.semi*cos(ksi);
    yrel=l.semi*sin(ksi);
    len1=cos(Gama)*xrel+sin(Beta)*sin(Gama)*yrel;
    len2=cos(Beta)*yrel;
    l.dis=sqrt(len1*len1+len2*len2)*AU/l.RE; 
    
   
    fac= sqrt(pow(10.0,l.logl)); 
    dhz1= 0.95*fac; ///CHZ
    dhz2= 1.69*fac;
    if(l.semi>=dhz1 and l.semi<= dhz2) l.HZP[0]=1;
    else  l.HZP[0]=0;  
    
    dhz1= 0.5*fac; ///OHZ
    dhz2= 2.0*fac;
    if(l.semi>=dhz1 and l.semi<= dhz2) l.HZP[1]=1;
    else  l.HZP[1]=0; 
    ///if(l.HZP[1]<1) {cout<<"Big error HZP[1]: "<<l.HZP[1]<<endl;  int we; cin>>we; }
    
    

    for(int i=0; i<M; ++i){
    if(fabs((l.RE/l.Dl/KP)*180.0*3600.0/M_PI)<fabs(FWHM[i]/2.0)) l.theta[i]=1;
    else l.theta[i]=0;
    s.Fluxb[i]= double(s.Fluxb[i]+l.theta[i]*pow(10.0,-0.4*l.Map[i])); 
    s.blendl[i]=fabs(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    s.magb[i]=-2.5*log10(fabs(s.Fluxb[i]));
    l.Fl[i]=double(l.theta[i]*pow(10.0,-0.4*l.Map[i]))/s.Fluxb[i]; 
    if(s.blendl[i]<0.0  or s.blendl[i]>1.0 or s.Fluxb[i]<=0.0  or s.magb[i]<0.0 or l.theta[i]>1){
    cout<<"ERROR Func_lens:  blendl:  "<<s.blendl[i]<<"\t fluxb:  "<<s.Fluxb[i]<<"\t magb:  "<<s.magb[i]<<endl;
    cout<<"filter:  "<<i<<"\t theta:  "<<l.theta[i]<<endl; cin>>yye;}}   
  
 
    if(s.ro_star<=0.0 or l.tE<=0.0 or l.Dl>s.Ds or l.Vt<=0.0 or l.dis<0.0 or s.blendl[0]>1.0 or l.q<0.0  or l.q>1.0){
    cout<<"ERROR ro_star:  "<<s.ro_star<<endl;
    cout<<"l.semi:  "<<l.semi<<endl;
    cout<<"Vt:  "<<l.Vt<<"\t l.Rstar:  "<<l.Rstar<<endl; 
    cout<<"RE: "<<l.RE/AU<<"\t xls:  "<<l.xls<<"\t tE: "<<l.tE<<endl;
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t Ds:  "<<s.Ds<<endl;
    cout<<"q:  "<<l.q<<"\t blendl:  "<<s.blendl[0]<<endl;
    cout<<"numl:  "<<l.numl<<"\t Vt:  "<<l.Vt<<"\t mul:  "<<l.mul<<endl;  cin>>yye;}
    
    
     cout<<"l.q:  "<<l.q<<"\t l.semi:  "<<l.semi<<"\t u0:    "<<l.u0<<endl;
    cout<<"ksi:   "<<l.ksi<<"\t ro_star:  "<<s.ro_star<<"\t l.dis:  "<<l.dis<<endl;
    cout<<"************** End of func_Lens  ****************"<<endl;
    
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
   xb = Dsun-x*cos(s.FI)*cos(s.TET);
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
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-Dsun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/Dsun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/Dsun,-2.44);///M_sun/pc^3
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
 double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*l.Dl*cos(s.TET)*cos(s.FI));
 double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+Dsun*Dsun-2.*Dsun*s.Ds*cos(s.TET)*cos(s.FI));
 if(Rlc==0.0) Rlc=0.00034346123;
 if(Rsc==0.0) Rsc=0.0004762654134;  
 ///Source and Lens velocity components in Galactocentric cylindrical coordinates
 double SVT, SVR, SVZ, LVT, LVR, LVZ,SVt,LVt;
 ///Source and Lens velocity components in heliocenteric galactic coordinates
 double SVb, SVl, LVb, LVl;
 double fv, testfv;
 double VSunl,VSunt,VSunb,vls_b,vls_l;
 double betal,betas,deltal,deltas,deltao;

 double NN=3.0;
 double VSunR =-10.3;
 double VSunT =vro_sun*(1.00762+0.00712)+6.3;
 double VSunZ = 5.9;
 double sigma_R_Disk= 43.0,  sigma_T_Disk= 27.8, sigma_Z_Disk=17.5;
 double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
 double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
 double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
 double Rho[8]={00.0}; double maxr=0.0;
 for(int i=0; i<8; ++i){  Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}
 
  double v_R_lmc=-57.0;
  double v_T_lmc=-226.0; 
  double v_Z_lmc= 221.0;
  double sigma_LMC=20.2; 
  double err_rlmc= 13.0; ///error of global velocity
  double err_tlmc= 15.0; 
  double err_zlmc= 19.0; 
 

  double test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])     {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0;}
else if(test<=(Rho[0]+Rho[1])) {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]))  {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))  
                           {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))  
                           {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                           {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4;}
else if(test<=maxr)        {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5;}
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}  


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/// Generate Source velocity components in Glactocenteric cylindrical coordinates(x',y')
    SVR=SVT=SVZ=0.0;
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_Disk, NN); 
    SVT= RandN(sigma_T_Disk, NN);
    SVZ= RandN(sigma_Z_Disk, NN); 
    SVT =SVT+ vro_sun *(1.00762 * pow(Rsc/Dsun,0.0394) + 0.00712);}

   else if(s.struc==1){///Galactic bulge
   SVZ=RandN(sigma_Z_Bulge, NN); 
   SVR=RandN(sigma_R_Bulge, NN);
   SVT=RandN(sigma_T_Bulge, NN);}

   else if(s.struc==2){///thick disk
   SVR= RandN(sigma_R_TDisk, NN); 
   SVT= RandN(sigma_T_TDisk, NN); 
   SVZ= RandN(sigma_Z_TDisk, NN);
   SVT =SVT+ vro_sun *(1.00762*pow(Rsc/Dsun,0.0394) + 0.00712); }
   
   else if(s.struc==3){///stellar halo
   SVR= RandN(sigma_R_halo, NN); 
   SVT= RandN(sigma_T_halo, NN); 
   SVZ= RandN(sigma_Z_halo, NN);}
   
   else if(s.struc>3){
   SVR= RandN(sigma_LMC, NN); 
   SVT= RandN(sigma_LMC, NN); 
   SVZ= RandN(sigma_LMC, NN); 
   SVZ +=  v_Z_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_zlmc;
   SVR +=  v_R_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_rlmc;
   SVT +=  v_T_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_tlmc;}
   s.vs=sqrt(SVT*SVT+SVZ*SVZ+SVR*SVR);

///======================================================================================
/// Generate Lens velocity components in Glactocenteric cylindrical coordinates(x',y'
    LVR=LVT=LVZ=0.0;
    if(l.struc==0){///Galactic disk
    LVR= RandN(sigma_R_Disk, NN); 
    LVT= RandN(sigma_T_Disk, NN);
    LVZ= RandN(sigma_Z_Disk, NN); 
    LVT =LVT+ vro_sun *(1.00762 * pow(Rlc/Dsun,0.0394) + 0.00712);}

   else if(l.struc==1){///Galactic bulge
   LVZ=RandN(sigma_Z_Bulge, NN); 
   LVR=RandN(sigma_R_Bulge, NN);
   LVT=RandN(sigma_T_Bulge, NN);}

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN); 
   LVT= RandN(sigma_T_TDisk, NN); 
   LVZ= RandN(sigma_Z_TDisk, NN);
   LVT =LVT+ vro_sun *(1.00762*pow(Rlc/Dsun,0.0394) + 0.00712); }
   
   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN); 
   LVT= RandN(sigma_T_halo, NN); 
   LVZ= RandN(sigma_Z_halo, NN);}
   
   else if(l.struc>3){
   LVR= RandN(sigma_LMC, NN); 
   LVT= RandN(sigma_LMC, NN); 
   LVZ= RandN(sigma_LMC, NN); 
   LVZ +=  v_Z_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_zlmc;
   LVR +=  v_R_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_rlmc;
   LVT +=  v_T_lmc + ((double)rand()/(double)(RAND_MAX+1.)*2.0-1.0)*err_tlmc;}
   l.vl=sqrt(LVT*LVT+LVZ*LVZ+LVR*LVR);
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc-1.0)<0.05) betal=pi/2.0; 
   else if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc+1.0)<0.05) betal=-pi/2.0; 
   else  betal=asin(l.Dl*cos(s.FI)*sin(s.TET)/Rlc);///lens[-pi/2,pi/2]
   
   if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc-1.0)<0.05) betas=pi/2.0; 
   else if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc+1.0)<0.05) betas=-pi/2.0; 
   else  betas=asin(s.Ds*cos(s.FI)*sin(s.TET)/Rsc);///lens[-pi/2,pi/2]
   
   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc)>1.05 || fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc)>1.05 || Rlc==0.0 || Rsc==0.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(s.TET)/Rlc<<"\t sin(s): "<<s.Ds*cos(s.FI)*sin(s.TET)/Rsc<<endl;
   int ew; cin>>ew; }



   
   if(fabs(l.Dl*cos(s.FI))>sqrt(pow(Dsun,2.0)+pow(l.Dl*cos(s.FI)*sin(s.TET),2.0)) ) betal= pi-betal;
   if(fabs(s.Ds*cos(s.FI))>sqrt(pow(Dsun,2.0)+pow(s.Ds*cos(s.FI)*sin(s.TET),2.0)) ) betas= pi-betas;

    if(fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))-1.0)<0.05)   deltal=0.0; 
    else if (fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))+1.0)<0.05) deltal=pi; 
    else    deltal=acos((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI)));
    
    if(fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))-1.0)<0.05)   deltas=0.0; 
    else if (fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))+1.0)<0.05) deltas=pi; 
    else    deltas=acos((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI)));
    
    if(fabs((Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI)))>1.05 or fabs((Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI)))>1.05 or fabs(s.FI)==pi/2.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"betal: "<<betal<<"\t betas: "<<betas<<endl; 
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"cos(dl): "<<(Rlc-Dsun*cos(betal))/(l.Dl*cos(s.FI))<<endl;
    cout<<"cos(ds): "<<(Rsc-Dsun*cos(betas))/(s.Ds*cos(s.FI))<<endl;  int ew; cin>>ew; }


    deltao=pi/2.0;
    SVl=-SVR*    sin(deltas)+ SVT* cos(deltas);
    LVl=-LVR*    sin(deltal)+ LVT* cos(deltal);
    VSunl=-VSunR*sin(deltao)+VSunT*cos(deltao);

    SVt=  1.0*SVR*cos(deltas)+  SVT*sin(deltas);
    LVt=  1.0*LVR*cos(deltal)+  LVT*sin(deltal);
    VSunt=1.0*VSunR*cos(deltao)+VSunT*sin(deltao);

    SVb=-sin(s.FI)*(SVt) + cos(s.FI)* (SVZ);
    LVb=-sin(s.FI)*(LVt) + cos(s.FI)* (LVZ);
    VSunb=-sin(s.FI)*(VSunt)+cos(s.FI)*(VSunZ);

    vls_l= LVl-l.xls*SVl -(1.0-l.xls)*VSunl;
    vls_b= LVb-l.xls*SVb -(1.0-l.xls)*VSunb;
    l.Vt=sqrt(fabs(vls_l*vls_l+ vls_b*vls_b));
    if (l.Vt<0.0 || l.Vt>1.0e6 ){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
