#include "tools.h"
#include "TGraph.h"
#include "TF1.h"
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "TMath.h"
#include <cmath>
tools::tools()
{

	//gROOT->ProcessLine(".L tools.cpp");
}
tools::~tools()
{
    
}
double tools::pol_n(double x,double *p,int i)
{
    double m = 0;
    int n;
    for(n=0;n<=i;n++){
        m=m+p[n]*pow(x,n);
    }
    return m;
}
double tools::gaus(double x,double *p){
    double m = 0;
    m = p[0]*TMath::Exp(-0.5*((x-p[1])/p[2])*((x-p[1])/p[2]));
    return m;
}
double tools::find_max(double *p,int i0)
{
    double x = 0;
    int k0;
    for(k0=0;k0<i0;k0++){
        if(p[k0]>x) x=p[k0];
    }
    return x;
}
double tools::find_min(double *p,int i0)
{
    double x = 100;
    int k0;
    for(k0=0;k0<i0;k0++){
        if(p[k0]<x) x=p[k0];
    }
    return x;
}
int tools::find_max_channel(double *p, int i0)
{
    double x = 0;
    int k0;
    for(k0=0;k0<i0;k0++){
        if(p[k0]>x) x=p[k0];
    }
    for(k0=0;k0<i0;k0++){
        if(p[k0]==x) break;
    }
    return k0;
}
int tools::find_min_channel(double *p, int i0)
{
    double x = 100;
    int k0;
    for(k0=0;k0<i0;k0++){
        if(p[k0]<x) x=p[k0];
    }
    for(k0=0;k0<i0;k0++){
        if(p[k0]==x) break;
    }
    return k0;
}
double tools::find_baseline(double *x,double *y,int i)
{
    /*TGraph *g = new TGraph(i,x,y);
    TF1 *f1 = new TF1("f1","[0]",-100,0);
    g->Fit("f1","Q");
    double par[3];
    f1->GetParameters(&par[0]);
    double m;
    m=par[0];
    delete g;
    g = NULL;
    delete f1;
    f1 = NULL;*/
    /*TGraph g(i,x,y);
    TF1 f1("f1","[0]",-100,0);
    g.Fit("f1","Q");
    double par[3];
    f1.GetParameters(&par[0]);
    double m;
    m=par[0];*/
    int n,mc;
    double total,m;
    total = 0;
    if(i > 0){
        for(n=0;n<i;n++){
            total = total + y[n];
        }
        m = total/i;
    }
    if(i == 0) m = 0;
    return m;
}
double tools::find_baseline(double *x,double *y,int i,int i0)//basenum,i
{
    int n,mc;
    double total,m;
    if(i > 0){
        mc = find_max_channel(y,i0);
        if(mc<i){
            //if(mc>40) i = mc-30;
            double min_0 = find_min(y,mc);
            double min_c = find_min_channel(y,mc);
            for(n=mc;n>0;n--){
                if(y[n]<0.01*(y[mc]-min_0)) break;
            }
            int t_num;
            t_num = (int)(5/(x[1]-x[0]));
            if((n-min_c)>t_num) i = n - (int)(3/(x[1]-x[0])); 
            else i = n - (int)((n-min_c)*0.6);
            //cout<<mc<<" "<<n<<" "<<min_0<<" "<<i<<endl;
        }
        total = 0;
        for(n=0;n<i;n++){
            total = total + y[n];
        }
        m = total/i;
    }
    if(i==0) m = 0;
    return m;
}
double tools::intergral(double *x,double *p,double bin,int i,int base_num){
    int k = 0;
    double a = 0;
    double base = 0;
    base = find_baseline(x,p,base_num,i);
    //cout<<base<<endl;
    for(k=0;k<i;k++){
        a = a + bin*(p[k]-base);
        //cout<<p[k]<<endl;
    }
    return a;
}
double tools::intergral_s(double *x,double *p,double bin,int i,int base_num,double intertime,double thr){
    double x0[7000],y0[7000],x1[7000],y1[7000],y2[7000],base,max,x_0,a,k01,k02;
    //tools *t1 = new tools();
    a = 0;
    int rise_m,n,k,m,xx,start_ch;
    base = find_baseline(x,p,base_num,i);
    start_ch = find_s_c(x,p,i,base_num);
    //width
    x_0 = 0;
    for(n=0;n<10;n++){
        x_0 = x_0 + (x[n+1]-x[n]);
    }
    x_0 = x_0*0.1;
    k01 = intertime/x_0;
    k02 = floor(k01);
    if(k02<k01) k = k02+1;
    else k = k01;
    for(m = start_ch;m<i;m++){
        if(p[m]>(base+thr)) break;
    }
    if((m+k)<i) xx=m+k;
    else xx=i;
    a = intergral(x,p,bin,xx,base_num);
    return a;
}
void tools::intergral_s_n(double *x,double *p,double bin,int i,int base_num,double *intertime,double thr,int in_num,double *e_n){
    double x0[7000],y0[7000],x1[7000],y1[7000],y2[7000],base,max,x_0,a,k01,k02;
    //tools *t1 = new tools();
    a = 0;
    int n_num = in_num;
    int rise_m,n,k[n_num],m,xx,start_ch;
    base = find_baseline(x,p,base_num,i);
    start_ch = find_s_c(x,p,i,base_num);
    //width
    x_0 = 0;
    for(n=0;n<10;n++){
        x_0 = x_0 + (x[n+1]-x[n]);
    }
    x_0 = x_0*0.1;
    for(n=0;n<in_num;n++){
        k01 = intertime[n]/x_0;
        k02 = floor(k01);
        if(k02<k01) k[n] = k02+1;
        else k[n] = k01;
    }
    for(m = start_ch;m<i;m++){
        if(p[m]>(base+thr)) break;
    }
    for(n=0;n<in_num;n++){
        if((m+k[n])<i) xx = m + k[n];
        else xx=i;
        e_n[n] = intergral(x,p,bin,xx,base_num);
    }
    //return a;
}
double tools::PSD(double *x,double *p,double bin,int i,int base_num,double intertime1,double intertime2,double thr){
    double a1,a2,psd_0;
    a1=intergral_s(x,p,bin,i,base_num,intertime1,thr);
    a2=intergral_s(x,p,bin,i,base_num,intertime2,thr);
    psd_0 = a1/a2;
    return psd_0;
}
int tools::find_s_c(double *x,double *p,int i,int base_num){
    int max_ch,k;
    double base;
    max_ch = find_max_channel(p,i);
    base = find_baseline(x,p,base_num,i);
    for(k=max_ch;k>0;k=k-1){
        if(p[k]<base) break;
    }
    return k;
}

double tools::getslope(double *x,double *p,int i,int base_num){
    int max_ch,k,start_ch;
    double base,slope,max;
    max_ch = find_max_channel(p,i);
    start_ch = find_s_c(x,p,i,base_num);
    base = find_baseline(x,p,base_num,i);
    max = find_max(p,i);
    slope = (max-base)/(x[max_ch]-x[start_ch]);
    return slope;
}

double tools::getslope_max(double *x,double *p,int i,int base_num){
    int max_ch,k,start_ch;
    double base,slope,max;
    max_ch = find_max_channel(p,i);
    start_ch = find_s_c(x,p,i,base_num);
    base = find_baseline(x,p,base_num,i);
    max = find_max(p,i);
    slope = 0;
    for(k=start_ch;k<max_ch;k++){
        double m;
        m = (p[k+1]-p[k])/(x[k+1]-x[k]);
        if(m>slope) slope = m;
    }
    return slope;
}

double tools::getslope_mean(double *x,double *p,int i,int base_num){
    int max_ch,k,start_ch;
    double base,slope,max,mean;
    max_ch = find_max_channel(p,i);
    start_ch = find_s_c(x,p,i,base_num);
    base = find_baseline(x,p,base_num,i);
    max = find_max(p,i);
    slope = 0;
    for(k=start_ch;k<max_ch;k++){
        double m;
        m = (p[k+1]-p[k])/(x[k+1]-x[k]);
        slope = slope+m;
    }
    mean = slope/(max_ch - start_ch);
    return mean;
}
double tools::getrisetime(double *x,double *p,int i,int base_num){
    int max_ch,k,start_ch,s;
    double base,slope,max,mean,t_0,t_1;
    max_ch = find_max_channel(p,i);
    start_ch = find_s_c(x,p,i,base_num);
    base = find_baseline(x,p,base_num,i);
    max = find_max(p,i);
    s = 0;
    for(k = start_ch;k<max_ch;k++){
        if(p[k]>(0.1*(max-base))&&s == 0){
            s++;
            t_0=(0.1*(max-base)-p[k-1])*(x[k]-x[k-1])/(p[k]-p[k-1]) + x[k-1];
        }
        if(p[k]>(0.9*(max-base))){
            t_1=(0.9*(max-base)-p[k-1])*(x[k]-x[k-1])/(p[k]-p[k-1]) + x[k-1];
            break;
        }
    }
    double risetime = t_1-t_0;
    return risetime;
}
