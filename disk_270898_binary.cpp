/// In the name of GOD ###   5/8/98  Elahi Be Omide Tou
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;
///****************   constants  *********************///
const double pi= M_PI;
const double Pc=3.0857*pow(10.0,16);///[m]
const double G= 6.67408*pow(10.0,-11);/// [m^3/kg.s^2]
const double c =299792458.0;///[m/s]
const double Msun=1.989*pow(10.0,30.0);///[kg]
const double Rsun=6.9551*pow(10.0,8);///[m]
const double mp = 1.6726219*pow(10.0,-27);///[kg]
const double Wav[2]={1.469003,0.869048};///micro meter[W149, Z087]
const double sigmaR_vis=5.1*pow(10.0,-31.0);// Rhylie cross section in 532 nm
const double Au =1.495978707*pow(10.0,11);///[m]
const double m0=2.37*1.00797*mp;///[kg]
const double hplanck= 6.62607004*pow(10.0,-34.0);///m2 kg / s
const double kbol=1.38064852*pow(10.0,-23);/// boltezman constant [m2.kg./s2.K]
const double b_con= 2.897771955*pow(10.0,-3.0);///Wein displacement law
const double opacity=150.0;//Figure(8) of Inner region of protoplanetary disks(2010)
const double extinc=0.8; 
const double mzp[2]={27.615,26.387};//zero point magnitude W149, Z087
///**************  Number constant  ********************///
const int Nz=int(10);
const int Ndisk=int(430);
const int Nlens=int(400);
const int Ntem=int(8979);///number of rows in file Fstar_disk.dat
const int Nwa=int(10000);
const int Ntet=int(89);
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


struct disk{
    double Ri,Ri2,Ro,Rcon,epsi;
    double Mdisk,inc,sigmai;
    double It[2][Ndisk],Ir[2][Ndisk];
    double Rd[Ndisk],Td[Ndisk], lamb[Ndisk];
    double diskm[2][2];//,alf[Ndisk];
    double opt[Ndisk],Dsur[Ndisk]; 
    double Rs, h0, tau[2], powR, powZ;
    double epr, Itot[2];  
    double wav[Nwa], Lnu[3][Nwa];///disk+source radiation versus wave length
};
struct lensing{
    double q, dis; 
    double Ml, Dl, xrel;
    double tstar, pt1, pt2; 
    double Ms, Ds, Rstar, Tstar, ro_star, Fstar[2];
    double RE, tE, t0, Vt, proj, u0, ksi;
};
/// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void func_lens(lensing & l);
void func_disk(disk & d, lensing & l);
void Iref(disk & d, lensing & l);
void Ithre(disk & d, lensing & l);
double randN(int N, double sigma);
double randR(double down, double up);
double f(double r, double z,double Rstar);
double J(double r, double Rstar);
double K(double r, double Rstar);
double Planck(double lam,  double Tem); 
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///=========================================
// SEED Generation By signlesskid-Part 0
    time_t _timeNow;
    unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///=========================================
int main(){

    disk d;
    lensing l;
    time(&_timeNow);
    
  srand(847356347180);  
    
    printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
    VBBinaryLensing *vbb=new VBBinaryLensing;
    vbb->Tol=1.e-3; 
    vbb->LoadESPLTable("./files/ESPL.tbl");
    FILE * fii;
    fii=fopen("./files/magniB.txt","w");
    FILE * fii2;
    fii2=fopen("./files/diskB.txt","w");
    FILE * para;
    para=fopen("./files/paramB.txt","w");


    func_lens(l);
    func_disk(d,l);
    Iref(d,l);
    Ithre(d,l);
    
    int counter=0,flag, i;
    double xcent, ycent,t;
    double dt= fabs(l.pt2-l.pt1)/Nlens;
    double Mag_s, Mag_t[2], del_m[2], Flux[2], map[2];
    double mag1, mag2, mag3, mag4, Mdisk[2];
    double xx,yy,r,ym,x,y;
    double areas=double(pi*l.Rstar*l.Rstar);
    double dx=d.epr/sqrt(1.0+cos(d.inc)*cos(d.inc)); 
    double dy=fabs(dx*cos(d.inc) );
    double area=fabs(dx*dy);
    double rho=sqrt(area/pi)*l.proj; 
    
    //if(rho<1.0e-5) 
    
  //  rho=0.0;
    
cout<<"rho: "<<rho<<endl; 
//int yye;  cin>>yye;   
   
    
    for(int il=0; il<Nlens; ++il){
    t=double(l.pt1+ il*dt + l.t0);
    xcent = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
    ycent = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
    Mag_s = vbb->BinaryMag2(l.dis, l.q, xcent , ycent ,l.ro_star);
    
    ///***********************************************************************************
    d.diskm[0][0]=d.diskm[0][1]=d.diskm[1][0]=d.diskm[1][1]=0.0;
    for(x=0.0;  x<d.Ro;  x+= dx){
    ym= sqrt(d.Ro*d.Ro-x*x)*cos(d.inc); 
    for(y=0.0;  y<ym;    y+= dy){
    r= sqrt(x*x + y*y/(cos(d.inc)*cos(d.inc)));
    if((r>=d.Ri and r<=d.Ri2) or (r>=d.Rcon and r<d.Ro)){
    xx= x*l.proj;
    yy= y*l.proj; 
    if(sqrt(x*x+y*y)<=l.Rstar)   flag=0;
    else                         flag=1;
    mag1= vbb->BinaryMag2(l.dis , l.q , xcent + xx , ycent + yy ,rho);
    mag2= vbb->BinaryMag2(l.dis , l.q , xcent - xx , ycent + yy ,rho);
    mag3= vbb->BinaryMag2(l.dis , l.q , xcent - xx , ycent - yy ,rho);
    mag4= vbb->BinaryMag2(l.dis , l.q , xcent + xx , ycent - yy ,rho);
    Mdisk[1]= mag1*flag + mag2*flag + mag3 + mag4;
    Mdisk[0]=double( 2.0*flag + 2.0);
    i=int(double(r-d.Ri)/d.epr);
    if(i<0.0) i=0;
    else if(i>(Ndisk-1))  i=Ndisk-1;
    for(int k=0; k<2; ++k){
    d.diskm[k][1] +=  (d.Ir[k][i]+ d.It[k][i])*area*Mdisk[1]; 
    d.diskm[k][0] +=  (d.Ir[k][i]+ d.It[k][i])*area*Mdisk[0];}
    if(dx<0.0 or area<0.0 or rho<0.0 or Mdisk[1]<0.0 or Mdisk[0]<0.0){
    cout<<"ERROR dr: "<<dx<<"\t area: "<<area<<"\t rho: "<<rho<<"\t Mdisk1: "<<Mdisk[1]<<"\t Mdisk0: "<<Mdisk[0]<<endl; int ye;  cin>>ye; }
    }}}
    ///***********************************************************************************
    
    
    
    for(int k=0; k<2; ++k){
    Mag_t[k]=(d.diskm[k][1]+Mag_s*areas)/(d.diskm[k][0]+1.0*areas);
    del_m[k]=(Mag_t[k]-Mag_s);
    Flux[k]= (d.diskm[k][1]+Mag_s*areas)*l.Fstar[k]*pow(l.Ds*1000.0*Pc/Rsun,-2.0)/(4.0*pi); 
    map[k] = -mzp[k]-2.5*log10(Flux[k]) + extinc;}
    fprintf(fii,"%.5lf  %.5lf  %.5lf  %.8lf   %.5lf  %.5lf  %e  %e  %e   %e  %.5lf  %.5lf\n",
    (t-l.t0)/l.tE,xcent,ycent,Mag_s,Mag_t[0],Mag_t[1],del_m[0],del_m[1],Flux[0],Flux[1],map[0],map[1]); 
    
    if(il%500==0){
    cout<<"****************************************************************"<<endl;
    cout<<"step: "<<il<<"\t time(tE): "<<(t-l.t0)/l.tE<<endl;
    cout<<"time: "<<t<<"\t dt: "<<dt<<"\t pt1: "<<l.pt1<<"\t pt2:  "<<l.pt2<<endl;
    cout<<"**** Mag_total[0]: "<<Mag_t[0]<<"\t Mag_total[1]: "<<Mag_t[1]<<"\t Mag_s: "<<Mag_s<<endl;
    cout<<"del_mag[0]: "<<del_m[0]<<"\t del_mag[1]: "<<del_m[1]<<endl;
    cout<<"xcent: "<<xcent<<"\t ycent: "<<ycent<<endl;
    cout<<"Mag_disk[W149]: "<<d.diskm[0][1]<<"\t Mag_disk[Z087]: "<<d.diskm[1][1]<<endl;
    cout<<"Mag0_disk[W149]: "<<d.diskm[0][0]<<"\t Mag0_disk[Z087]: "<<d.diskm[1][0]<<endl;
    cout<<"Flux_tot[0]: "<<Flux[0]<<"\t Flux_tot[1]: "<<Flux[1]<<endl;
    cout<<"magnitude[0]: "<<map[0]<<"\t magnitude[1]: "<<map[1]<<endl;
    cout<<"****************************************************************"<<endl;} 
    
    }   
    for (double tet=0; tet<359.0; tet+=0.1){
    x=l.proj*cos(tet*pi/180.0);     
    y=l.proj*sin(tet*pi/180.0)*cos(d.inc);
    fprintf(fii2,"%e %e %e %e  %e %e %e %e %e %e\n",x*l.Rstar,y*l.Rstar,x*d.Ri,y*d.Ri,x*d.Ri2,y*d.Ri2,x*d.Rcon,y*d.Rcon,x*d.Ro,y*d.Ro);}
    fclose(fii2);
    fprintf(para,"%d  %.5lf  %.5lf %.5lf %.5lf  %.5lf  %.5lf  %e  %.5lf  %.5lf  %.5lf  %.5lf %e  %e  %e  %e  %e  %.5lf  %.5lf %.5lf  %.5lf\n",
    counter,l.ro_star,d.powR,d.powZ,d.Rs,d.h0,d.inc*180.0/pi,d.Mdisk/l.Ms/Msun,d.Ri,d.Ri2,d.Rcon,d.Ro,d.sigmai,d.tau[0],d.tau[1],l.Fstar[0],l.Fstar[1], l.q, l.dis, l.u0/l.ro_star, l.ksi);

    fclose(fii);
    fclose(para);  
    time(&_timeNow);
    printf("END_TIME >>>>>>>>  %s ",ctime(&_timeNow));
    delete vbb; 
    return(0);
    }
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_lensing                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void func_lens(lensing & l){
    l.Ml = 0.3;///[M_sun]
    l.Dl = 6.5;//* 1000.0* Pc;///[m]
    l.Ds = 8.0;//* 1000.0 *Pc;///[m]
    l.xrel=double(l.Dl/l.Ds);
    l.RE =sqrt(4.0*G*l.Ml*Msun*l.Ds*1000.0*Pc)*sqrt(l.xrel*(1.0-l.xrel))/c ;///[m]
    l.proj=double(l.xrel*Rsun/l.RE);
    l.Tstar=5000.69758;///[K]
    l.Rstar=1.0;/// [Rsun]
    l.Ms=1.0;///[Msun]
    l.t0=0.0;
    l.Vt= 200.0;///[km/s]
    l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0);///[days]
    l.ksi=0.0*pi/180.0;///[radian]
    l.ro_star=l.Rstar*l.proj;
    l.u0= l.ro_star*5.0;
    l.tstar= fabs(l.tE*l.ro_star);///days 
    l.pt1=-1.2*l.tE;///days 
    l.pt2=+1.2*l.tE;///days 
    
    l.q=1.0;
    l.dis=1.0; 
    
    cout<<"******************* FUNC_LENS *****************"<<endl;
    cout<<"Mass_lens[Msun]: "<<l.Ml<<"\t lens_dis[kpc]:  "<<l.Dl<<"\t RE[AU]: "<<l.RE/Au<<endl;
    cout<<"Ds[Kpc]: "<<l.Ds<<"\t xrel: "<<l.xrel<<"\t Rstar[Rsun]: "<<l.Rstar<<"\t Tstar[K]:  "<<l.Tstar<<endl;
    cout<<"ro_star: "<<l.ro_star<<"\t u0: "<<l.u0<<"\t tE[days]: "<<l.tE<<endl;
    cout<<"****************** END OF FUNCTION OF LENS *********"<<endl; 
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_disk                               ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void func_disk(disk & d, lensing & l)
{
    double H, ff; 
    d.powR=1.40+randR(-0.15,0.15);///[1.0:1.5] from Jami et al. 2017   
    d.powZ=1.40+randR(-0.15,0.15);//power of R in Z scale, protoplanetary_disk_review Hawaii
    d.Rs=double(100.0 + randR(-20.0,20.0) )*Au/Rsun;///[Rsun]
    d.h0= (10.0+randR(-2.0,2.0))*Au/Rsun;///[Rsun]
 d.inc=30.0*pi/180.0;// randR(1.0,89.0)*pi/180.0;///randian
    d.Mdisk=pow(10.0,-4.0-randN(1,1.5))*l.Ms*Msun;///[kg]
    d.Ri=  randR(0.02,0.025)*Au/Rsun;// Inner region of protoplanetary disks [Rsun]
    d.Ri2= randR(0.025,0.04)*Au/Rsun;// Inner region of protoplanetary disks [Rsun]
    d.epsi= randR(0.5,1.0);/// from review paper of protoplanetary disks
    d.Rcon=l.Rstar*pow(l.Tstar/1500.0,2.0)*0.5*pow(d.epsi,-0.5);//*sqrt(1.0+fac);/// consentration radius [Rsun] Jami 2017
    if(d.Ri2>d.Rcon)  d.Ri2=d.Rcon;
    d.Ro= pow(0.05,-1.0/d.powR)*d.Rcon;///outer radius [Rsun]
    d.epr=(d.Ro-d.Ri)/Ndisk/1.0;   ///log10(d.Ro/d.Ri)/Ndisk/1.0; 
    ff= (pow(d.Ro/d.Rcon,2.0-d.powR)-1.0)/(2.0-d.powR);
    d.sigmai= d.Mdisk/(2.0*pi*d.Rcon*d.Rcon*Rsun*Rsun)/ff;///[kg/m^2]
    d.tau[0]= d.sigmai*sigmaR_vis*pow(0.532/Wav[0],4.0)/m0;
    d.tau[1]= d.sigmai*sigmaR_vis*pow(0.532/Wav[1],4.0)/m0;


    
    for(int i=0; i<Ndisk; ++i){
    d.Rd[i]=d.Td[i]=d.lamb[i]=d.Dsur[i]=d.opt[i]=0.0; 
    d.Rd[i]=d.Ri + d.epr*i ;///pow(10.0,log10(d.Ri)+d.epr*i);
    H= d.h0*pow(d.Rd[i]/d.Rs,d.powZ); ///[Rsun]
    if(d.Rd[i]<d.Rcon) d.Td[i]=double(l.Tstar*pow(l.Rstar/d.Rd[i],0.5)); 
    else               d.Td[i]=double(l.Tstar*pow(l.Rstar*0.5/d.Rd[i],0.5)*pow(d.epsi,-0.25));      
    d.lamb[i]= double(b_con*1000000.0/d.Td[i]);///[micro-meter] 
    if(d.Rd[i]<=d.Ri2)  d.opt[i]=0.3;//from gasusian inner disk (Jami et al.  2017)
    else if(d.Rd[i]>d.Ri2 and d.Rd[i]<d.Rcon) d.opt[i]=0.0;/// nothing
    else {///dust
    d.Dsur[i]=d.sigmai*pow(d.Rcon/d.Rd[i],d.powR);
    d.opt[i]= 0.1*opacity*d.Dsur[i]/cos(d.inc);}
    if(H<0.0 or d.Td[i]<0.0 or d.opt[i]<0.0 or d.Rd[i]>d.Ro or d.Rd[i]<d.Ri or d.sigmai==0.0 or d.sigmai<0.0 or d.Dsur[i]<0.0 or ff<0.0){
    cout<<"ERROR step: "<<i<<endl;
    cout<<"Rd: "<<d.Rd[i]<<"\t H: "<<H<<"\t Sigma_sur: "<<d.Dsur[i]<<endl;
    cout<<"tem: "<<d.Td[i]<<"\t Wave: "<<d.lamb[i]<<"\t optical_depth:  "<<d.opt[i]<<endl;
    cout<<"*******************************************"<<endl;  int eew;  cin>>eew;}}
   
    for(int i=0; i<Nwa; ++i){
    d.wav[i]= pow(10.0,double(3.5*i/Nwa-0.5));///in micro meter and in interval [300nm,1mm]
    for(int j=0; j<3; ++j)  d.Lnu[j][i]=0.0;
    d.Lnu[0][i]= Planck(d.wav[i],l.Tstar)*2.0*pi*l.Rstar*l.Rstar;///[W/wave_length]/Rsun^2
    }

    cout<<"****************** DISK FUNCTION *******************"<<endl;
    cout<<"Rs: "<<d.Rs<<"\t h0: "<<d.h0<<"\t u0: "<<l.u0<<endl;
    cout<<"tau[0]: "<<d.tau[0]<<"\t tau[1]:  "<<d.tau[1]<<"\t sigmai: "<<d.sigmai<<endl;
    cout<<"incli(degree): "<<d.inc*180.0/pi<<"\t Mass_disk[M_star]: "<<d.Mdisk/(l.Ms*Msun)<<endl;
    cout<<"powRin density: "<<d.powR<<"power_Z_density: "<<d.powZ<<endl;
    cout<<"Inner_radius[Rsun]: "<<d.Ri<<"\t Inner radius_2[RSun]: "<<d.Ri2<<endl;
    cout<<"R_condensation: "<<d.Rcon<<"\t R_out: "<<d.Ro<<endl;
    cout<<"****************** END OF FUNCTION OF DISK *********"<<endl; 
    
 exit(0);   
    
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        INTENCITY_DISK                          ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Iref(disk & d, lensing & l)
{
   double H,up,down=0.0, dz,r;

   for(int i=0;  i<Ndisk; ++i){
   d.Ir[0][i]=d.Ir[1][i]=0.0;
   if(d.Rd[i]>=d.Rcon){
   H=d.h0*pow(d.Rd[i]/d.Rs,d.powZ);
   up=Nz*H;
   dz=(up-down)/1000.0;
   for(double z=down; z<=up; z+=dz){
   r= sqrt( z*z + d.Rd[i]*d.Rd[i] );
   d.Ir[0][i]+= f(r,z,l.Rstar)*exp(-0.5*z*z/H/H);
   d.Ir[1][i]+= f(r,z,l.Rstar)*exp(-0.5*z*z/H/H);}
   d.Ir[0][i] = d.Ir[0][i]*0.75*sigmaR_vis*pow(0.532/Wav[0],4.0)*d.Dsur[i]*dz*(l.Rstar*l.Rstar/d.Rd[i]/d.Rd[i])/(sqrt(2.0*pi)*H*m0); 
   d.Ir[1][i] = d.Ir[1][i]*0.75*sigmaR_vis*pow(0.532/Wav[1],4.0)*d.Dsur[i]*dz*(l.Rstar*l.Rstar/d.Rd[i]/d.Rd[i])/(sqrt(2.0*pi)*H*m0);}}

}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Ithre(disk & d, lensing & l)
{
   FILE * inten1;
   inten1= fopen("./files/Intensity_radiusB.txt","w");
   FILE * inten2;
   inten2= fopen("./files/Intensity_waveB.txt","w");
   FILE * tems;
   tems= fopen("./files/FstarT_disk.dat","r");
   
   double Mw[Ntem],Mz[Ntem],Tem[Ntem],frac,area;
   for(int i=0; i<Ntem; ++i){
   fscanf(tems,"%lf   %lf  %lf\n",&Tem[i],&Mw[i],&Mz[i]); 
   if(( (l.Tstar-Tem[i])*(l.Tstar-Tem[i-1])<0.0  and  i>0) or l.Tstar==Tem[i]){
   l.Fstar[0]= pow(10.0,-0.4*Mw[i]);
   l.Fstar[1]= pow(10.0,-0.4*Mz[i]);
   cout<<"Mstar[0]:  "<<Mw[i]<<"\t Mstar[1]:  "<<Mz[i]<<"\t l.Tstar: "<<l.Tstar<<endl;
   cout<<"Fstar[0]:  "<<l.Fstar[0]<<"\t Fstar[1]:  "<<l.Fstar[1]<<endl;}}
   fclose(tems);

   for(int i=1; i<Ndisk; ++i){
   for(int k=0; k<Nwa; ++k){
   frac=pow(Wav[0]/d.wav[k],4.0);
   area= 2.0*pi*d.Rd[i]*(d.Rd[i]-d.Rd[i-1])*cos(d.inc); 
   d.Lnu[1][k] +=area*(1.0-exp(-d.opt[i]))*Planck(d.wav[k],d.Td[i]);///thermal
   d.Lnu[2][k] +=area*d.Ir[0][i]*frac*Planck(d.wav[k],l.Tstar); }///reflected  

   int j=-1;
   if(d.Td[i]<Tem[0] or int(d.Td[i])==int(Tem[0]) ) j=0;   
   else if(d.Td[i]>Tem[Ntem-1]) j=int(Ntem-1); 
   else {
   for(int t=1; t<Ntem; ++t){
   if( (d.Td[i]-Tem[t])*(d.Td[i]-Tem[t-1])<0.0 or int(d.Td[i])==int(Tem[t])){j=t;   break;}}}
   if(j<0  or  j>(Ntem-1)){
   cout<<"Error Td[i]: "<<d.Td[i]<<"\t j:  "<<j<<endl;  int eew;  cin>>eew;}
      
   d.It[0][i]= pow(10.0,-0.4*Mw[j])*(1.0-exp(-d.opt[i]))/l.Fstar[0];
   d.It[1][i]= pow(10.0,-0.4*Mz[j])*(1.0-exp(-d.opt[i]))/l.Fstar[1];
   d.Itot[0]=   d.It[0][i] + d.Ir[0][i];   
   d.Itot[1]=   d.It[1][i] + d.Ir[1][i]; 
   fprintf(inten1,"%.5lf  %.2lf  %.5lf   %e   %e  %e   %e  %e %e  %.5lf  %.5lf  %e  %e\n",d.Rd[i],d.Td[i],d.lamb[i],d.It[0][i],d.It[1][i],d.Ir[0][i],d.Ir[1][i],d.Itot[0],d.Itot[1],d.opt[i],d.Dsur[i],1.0-exp(-d.opt[i]),pow(10.0,-0.4*Mw[j])/l.Fstar[0]);}
   fclose(inten1);
   
   for(int i=0; i<Nwa; ++i)
   fprintf(inten2,"%.5lf   %e   %e   %e\n",d.wav[i],d.Lnu[0][i],d.Lnu[1][i],d.Lnu[2][i]); 
   fclose(inten2);
   cout<<" ******************* End IF INTENCITY FUNCTION **********************"<<endl;
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        PLANCK FUNCTION                         ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double Planck(double lam,  double Tem){
    double b= double(8.0*pi*hplanck*c*c);
    double co=double(hplanck*c*1000000.0/kbol/Tem/lam); 
    double pp= b*pow(lam*0.000001,-5.0)/(exp(co)-1.0);//[W/m^2/m]
    if(pp<0.0  or co==0.0 or lam<0.0 or Tem<0.0 or b<0.0 or co<0.0){
    cout<<"ERROR planck dis:  "<<pp<<"\t wave: "<<lam<<"\t Tem: "<<Tem<<"\t power:  "<<co<<endl;
    int yye;  cin>>yye; }
    return(pp);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double f(double r, double z,double Rstar){
    return((3.0*J(r,Rstar)-K(r,Rstar))+(3.0*K(r,Rstar)-J(r,Rstar))*(z*z/r/r) );
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double J(double r,double rstar){
    if(r==rstar or r<rstar){
    cout<<"Error (J): r:  "<<r<<"\t R_star: "<<rstar<<endl;
    return(0.5) ; }
    return(0.5*(1.0-pow(1.0-rstar*rstar/r/r,0.5) ));
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double K(double r,double rstar){
    if(r==rstar or r<rstar){
    cout<<"Error (K): r:  "<<r<<"\t R_star: "<<rstar<<endl;
    return(1.0/6.0); }
    return((1.0-pow(1.0-rstar*rstar/r/r,1.5))/6.0);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double randN(int N, double sigma){
    double p, fp, frand;
    do{
    p=(double)rand()/((double)(RAND_MAX)+(double)(1))*N*sigma; ///[-N:N]
    fp=exp(-p*p/(2.*sigma*sigma));
    frand= (double)rand()/((double)(RAND_MAX)+(double)(1));
    }while(frand>fp);
    double sign= (double)rand()/((double)(RAND_MAX)+(double)(1));
    if(sign<0.5)     return(p);
    else             return(-p);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double randR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

