#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TF1.h"
#include "TVersionCheck.h"
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <typeinfo>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit.h>
#include "fitfunction.h"
#include "tools.h"

#pragma comment(lib,"libgsl_d.lib")
#pragma comment(lib,"libgslcblas_d.lib")
using namespace std;
Fitfunction::Fitfunction(){
	//gROOT->ProcessLine(".L tools.cpp");
}
Fitfunction::~Fitfunction(){
	
}
Double_t fitf(Double_t *x,Double_t *par) {
    Double_t fitval;
    //par[1]=-20;
    if(x[0]>par[1]){
        fitval = par[0]*TMath::Exp(-(x[0]-par[1])/par[2])*(1-TMath::Exp(-(x[0]-par[1])/par[3]))+par[4];
    }
    else{
        fitval = par[4];
    }
    return fitval;
}
Double_t fitf2(Double_t *x,Double_t *par){
    Double_t fitval;
    if(x[0]>par[1]){
        fitval = par[0]*(x[0]-par[1])+par[2];
    }
    else{
        fitval = par[2];
    }
    return fitval;
}
Double_t fitf3(Double_t *x,Double_t *par){
    Double_t fitval;
    if(x[0]>par[1]){
        fitval = par[0]*(x[0]-par[1])*(x[0]-par[2])+par[3];
    }
    else fitval = par[3];
    return fitval;
}
Double_t fitf4(Double_t *x,Double_t *par){
    Double_t fitval;
    if(x[0]>par[1]){
        fitval = (par[3]-par[2])*TMath::Exp(par[0]*(x[0]-par[1])) + par[2];
    }
    else fitval = par[3];
    return fitval;
}
Double_t fitf5(Double_t *x,Double_t *par){
    Double_t fitval;
    fitval = par[0]*x[0]+par[1];
    return fitval;
}
int find_rise(double *p,int i0,double a,double b){  //寻上升沿大致范围
    int k0,k1;
    for(k0=0;k0<i0;k0++){
        if(p[k0]>(a+b)) break;
    }
    for(k1=k0;k1<i0;k1++){
        if((p[k1+2]-p[k1+1])/(p[k1+1]-p[k1])<0.85) break;
    }
    return k1+1;
}
int find_rise2(double *p,int i0,double a,double b){    //寻上升沿道址
    int k0,k1;
    for(k0=0;k0<i0;k0++){
        if(p[k0]>(a+b)) break;
    }

    return k0+1;
}
double Fitfunction::fit_function_pol1(double *x,double *y,int i,int base_num,double thr){
    double x0[7000],y0[7000],base,max;
    int rise_m,n,k;

    base = t1.find_baseline(x,y,base_num,i);   //find baseline

    max = t1.find_max(y,i);   //find peak

    for(n=0;n<i;n++){
        if(y[n]>(0.1*(max-base)+base)) break;
    }

    for(k=0;k<i;k++){
        x0[k] = x[n];
        y0[k] = y[n];
        if(y0[k]>(0.2*(max-base)+base)) break;
        n++;
    }
    if(k==0){
        k++;
        n++;
        x0[k] = x[n];
        y0[k] = y[n];
    }
    double par0[2],par1[3];
    TGraph *gr0 = new TGraph(k+1,x0,y0);
    TF1 *f_1 = new TF1("f_1","pol1");
    gr0->Fit("f_1","Q");

    f_1->GetParameters(&par0[0]);
    double t_0,t_1;
    t_0 = (base-par0[0])/par0[1];

    int n0;
    n0 = find_rise(y,i,base,0.1*(max-base));
    TGraph *gr1 = new TGraph(n0,x,y);
    TF1 *f_2 = new TF1("f_2",fitf2,-100,100,3);
    f_2->SetParameters(par0[1],t_0,base);
    //f_2->FixParameter(2,base);
    //gr1->Draw();
    gr1->Fit("f_2","Q");

    f_2->GetParameters(&par1[0]);
    //cout<<par1[0];
    if(par1[0]==0) t_1=300;
    else t_1 = thr/par1[0] + par1[1];
    //cout<<t_1<<endl;
    delete f_1;
    delete f_2;
    delete gr0;
    delete gr1;
    return t_1;
}
//线性插值
double Fitfunction::fit_function_linear(double *x,double *y,int i,int base_num,double thr){
    double base,t_0;
    int m,n,start_c;
    //tools t1;
    start_c = t1.find_s_c(x,y,i,base_num);
    base = t1.find_baseline(x,y,base_num,i);
    for(m=start_c;m<i;m++){
        if(y[m]>(thr+base)) break;
    }
    t_0 = (thr+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    return t_0;
}
void Fitfunction::fit_function_linear_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx){
    double base,t_0;
    int m,n,start_c,d;
    //tools t1;
    start_c = t1.find_s_c(x,y,i,base_num);
    base = t1.find_baseline(x,y,base_num,i);
    d = 0;
    for(m=start_c;m<i;m++){
        if(y[m]>(thr[d]+base)){
            tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
            d++;
            m--;
            if(d==nx) break;
        }
    }
    //return t_0;
}
double Fitfunction::fit_function_pol2(double *x,double *y,int i,int base_num,double thr){
    double x0[7000],y0[7000],base,max;
    int rise_m,n,k;
    //tools t1;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    //cout<<max<<" "<<base<<endl;
    for(n=0;n<i;n++){
        if(y[n]>(0.03*(max-base)+base)) break;
    }

    for(k=0;k<i;k++){
        x0[k] = x[n];
        y0[k] = y[n];
        if(y0[k]>(0.06*(max-base)+base)) break;
        n++;
    }
    //cout<<k<<endl;
    if(k==0){
        k++;
        n++;
        x0[k] = x[n];
        y0[k] = y[n];
    }
    double par0[2],par1[3];
    TGraph *gr0 = new TGraph(k+1,x0,y0);
    TF1 *f_1 = new TF1("f_1","pol1");
    gr0->Fit("f_1","Q");
    f_1->GetParameters(&par0[0]);
    double t_0,t_1;
    t_0 = (base-par0[0])/par0[1];
    //cout<<t_0<<endl;
    int n0;
    n0 = find_rise2(y,i,base,0.2*(max-base));
    //cout<<n0<<endl;
    TGraph *gr1 = new TGraph(n0,x,y);
    TF1 *f_2 = new TF1("f_2",fitf3,-100,100,4);
    //f_2->SetParameters(par0[1],t_0,base);
    f_2->SetParameter(1,t_0);
    f_2->SetParameter(3,base);
    //f_2->FixParameter(2,base);
    //gr1->Draw();
    gr1->Fit("f_2","Q");

    f_2->GetParameters(&par1[0]);
    t_1 = thr/par1[0] + par1[1];
    //cout<<t_1<<endl;
    delete f_1;
    delete f_2;
    delete gr0;
    delete gr1;
    return t_1;
}
//gaus
double Fitfunction::fit_function_gaus(double *x,double *y,int i,int base_num,double thr,double range){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    int rise_m,n,k,p,zero_ch;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_3 = new TF1("f_3","gaus");
    gr1->Fit("f_3","Q");
    //gr1->Draw();
    f_3->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0&&par1[0]==0) t_1=300;
    else{
    for(x2=x[zero_ch]-0.15;x2<200;x2=x2+0.005){
        y2=t1.gaus(x2,par1);
        if(y2>base+thr){
            t_1=x2;
            break;
        }
    }
    }
    delete f_3;
    delete gr1;
    return t_1;
}
//gaus-linear
double Fitfunction::fit_function_gaus_l(double *x,double *y,int i,int base_num,double thr,double range){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    int rise_m,n,k,p,zero_ch,m;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_3 = new TF1("f_3","gaus");
    gr1->Fit("f_3","Q");
    //gr1->Draw();
    f_3->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0&&par1[0]==0){
        for(m=0;m<i;m++){
            if(y[m]>(thr+base)) break;
        }
        t_1 = (thr+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    }
    else{
    for(x2=x[zero_ch]-0.15;x2<200;x2=x2+0.005){
        y2=t1.gaus(x2,par1);
        if(y2>base+thr){
            t_1=x2;
            break;
        }
    }
    }
    delete f_3;
    delete gr1;
    return t_1;
}
//gaus thr
void Fitfunction::fit_function_gaus_l_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx,double range){
	int rise_m,n,k,p,m,d,zero_ch;
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
	TF1 *f_2 = new TF1("f_2","gaus");
	gr1->Fit("f_2","Q");
	f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0&&par1[0]==0){
        d=0;
        for(m=zero_ch;m<i;m++){
            if(y[m]>(thr[d]+base)){
                tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
                d++;
                m--;
                if(d==nx) break;
            }
        }
    }
    else{
        d=0;
        for(x2=t_0-0.15;x2<200;x2=x2+0.005){
            y2=t1.gaus(x2,par1);
            if(y2>base+thr[d]){
                double y2_0 = t1.gaus(x2-0.005,par1);
                tx[d]=(thr[d]+base-y2_0)*0.005/(y2-y2_0) + x2 - 0.005;
                d++;
                x2=x2-0.005;
                if(d==nx) break;
            }
        }
    }
    delete f_2;

    delete gr1;

}
//interpolation
double Fitfunction::fit_function_spline(double *x,double *y,int i,int base_num,double thr){
    double x0[50000],y0[50000],base,max,t_0;
    int nb,m,n;
    //n=3000;
    nb=i*10;
    base = t1.find_baseline(x,y,base_num,i);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    const gsl_interp_type *t = gsl_interp_cspline_periodic;
    gsl_spline *spline = gsl_spline_alloc(t,i);
    gsl_spline_init(spline,x,y,i);
    for ( Int_t k = 0; k < (nb-2); k++ )
    {
        x0[k]  = (Double_t) k*(x[i-1]-x[0])/(nb-1) + x[0];
        y0[k] = gsl_spline_eval(spline,x0[k],acc);
    }
    for(m=0;m<nb-2;m++){
        if(y0[m]>(thr+base)) break;
    }
    t_0 = (thr+base-y0[m-1])*(x0[m]-x0[m-1])/(y0[m]-y0[m-1]) + x0[m-1];
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    return t_0;
}

void Fitfunction::fit_function_spline_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx){
    //cout<<"asd"<<endl;
    double x0[50000],y0[50000],base,max,t_0;
    int nb,m,n,start_c;
    //n=3000;
    nb=i*10;
    base = t1.find_baseline(x,y,base_num,i);
    start_c = t1.find_s_c(x,y,i,base_num);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    const gsl_interp_type *t = gsl_interp_cspline_periodic;
    gsl_spline *spline = gsl_spline_alloc(t,i);
    gsl_spline_init(spline,x,y,i);
    for ( Int_t k = 0; k < (nb-2); k++ )
    {
        x0[k]  = (Double_t) k*(x[i-1]-x[0])/(nb-1) + x[0];
        y0[k] = gsl_spline_eval(spline,x0[k],acc);
    }
    int d = 0;
    for(m=(10*start_c-1);m<nb-2;m++){
        if(y0[m]>(thr[d]+base)){
            tx[d] = (thr[d]+base-y0[m-1])*(x0[m]-x0[m-1])/(y0[m]-y0[m-1]) + x0[m-1];
            d++;
            m--;
            if(d==nx) break;
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}
//polynomial fit
double Fitfunction::fit_function_poln(double *x,double *y,int i,int base_num,double thr,double range){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    int rise_m,n,k,p,zero_ch;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_2 = new TF1("f_2","pol4");
    //gr1->Draw();
    gr1->Fit("f_2","Q");
    f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0) t_1=300;
    else{
        for(x2=t_0;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr){
                t_1=x2;
                break;
            }
        }
        
    }
    
    delete f_2;
    delete gr1;
    return t_1;
}
//pol_n+linear
double Fitfunction::fit_function_poln_l(double *x,double *y,int i,int base_num,double thr,double range){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    int rise_m,n,k,p,zero_ch,m;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_2 = new TF1("f_2","pol4",x1[0],x1[p]);
    //gr1->Draw();
    gr1->Fit("f_2","Q");
    f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0){
        for(m=0;m<i;m++){
            if(y[m]>(thr+base)) break;
        }
        t_1 = (thr+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    }
    else{
        for(x2=t_0;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr){
                t_1=x2;
                break;
            }
        }
        
    }
    
    delete f_2;
    delete gr1;
    return t_1;
}
//pol4+linear thr
void Fitfunction::fit_function_poln_l_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx,double range){
	int rise_m,n,k,p,m,d,zero_ch;
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    if((p+1)>5){
        TGraph *gr1 = new TGraph(p+1,x1,y1);
        TF1 *f_2 = new TF1("f_2","pol4");
        //f_2->SetParameters(par0[1],t_0,base);
        //f_2->FixParameter(2,base);
        //gr1->Draw();
        gr1->Fit("f_2","Q");
        f_2->GetParameters(&par1[0]);
        delete f_2;
        delete gr1;
    }
    else{
        par1[1]=0;
        par1[2]=0;
    }
    double x2,y2;
    /*for(d=0;d<nx;d++){
    if(par1[1]==0&&par1[2]==0){
        for(m=0;m<i;m++){
            if(y[m]>(thr[d]+base)) break;
        }
        tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    }
    else{
        for(x2=t_0;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr[d]){
                tx[d]=x2;
                break;
            }
        }
    }
    }*/
    if(par1[1]==0&&par1[2]==0){
        d=0;
        for(m=zero_ch;m<i;m++){
            if(y[m]>(thr[d]+base)){
                tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
                d++;
                m--;
                if(d==nx) break;
            }
        }
    }
    else{
        d=0;
        for(x2=t_0;x2<200;x2=x2+0.005){
            y2=t1.pol_n(x2,par1,4);
            if(y2>base+thr[d]){
                double y2_0 = t1.pol_n(x2-0.005,par1,4);
                tx[d]=(thr[d]+base-y2_0)*0.005/(y2-y2_0) + x2 - 0.005;
                d++;
                x2=x2-0.005;
                if(d==nx) break;
            }
        }
    }
	//t_1 = thr/par1[0] + par1[1];
    //cout<<t_1<<endl;
    //return t_1;	
}
//poln gsl
void Fitfunction::fit_function_gsl_poln_l_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx,double range,int order){
	int rise_m,n,k,p,m,d,zero_ch;
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10],chisq;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)&&p>(order-1)) break;  //range
        n++;
    }
    //gsl_fit
    if((p+1)>order){
        gsl_matrix *x_n,*cov;
        gsl_vector *y_n,*c;
        x_n = gsl_matrix_alloc(p+1,order+1);
        y_n = gsl_vector_alloc(p+1);
        c = gsl_vector_alloc(order+1);
        cov = gsl_matrix_alloc(order+1,order+1);
        for(k=0;k < p+1;k++){
            for(m = 0;m < order+1;m++){
                gsl_matrix_set(x_n,k,m,pow(x1[k],m));
            }
            gsl_vector_set(y_n,k,y1[k]);
        }
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(p+1,order+1);
        gsl_multifit_linear(x_n,y_n,c,cov,&chisq,work);
        for(k=0;k<order+1;k++){
            par1[k] = gsl_vector_get(c,(k));
        }
        gsl_multifit_linear_free(work);
        gsl_matrix_free(x_n);
        gsl_matrix_free(cov);
        gsl_vector_free(y_n);
        gsl_vector_free(c);
    }
    else{
        par1[1]=0;
        par1[2]=0;
    }
    double x2,y2;
    if(par1[1]==0&&par1[2]==0){
        d=0;
        for(m=zero_ch;m<i;m++){
            if(y[m]>(thr[d]+base)){
                tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
                d++;
                m--;
                if(d==nx) break;
            }
        }
    }
    else{
        d=0;
        for(x2=t_0;x2<200;x2=x2+0.005){
            y2=t1.pol_n(x2,par1,order);
            if(y2>base+thr[d]){
                double y2_0 = t1.pol_n(x2-0.005,par1,order);
                tx[d]=(thr[d]+base-y2_0)*0.005/(y2-y2_0) + x2 - 0.005;
                d++;
                x2=x2-0.005;
                if(d==nx) break;
            }
        }
    }
}
//poln old
/*double Fitfunction::fit_function_poln(double *x,double *y,int i,int base_num,double thr){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max;
    int rise_m,n,k,p;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    for(n=0;n<i;n++){
        if(y[n]>(0.05*(max-base)+base)) break;
    }

    for(k=0;k<i;k++){
        x0[k] = x[n];
        y0[k] = y[n];
        if(y0[k]>(0.1*(max-base)+base)) break;
        n++;
    }
    if(k==0){
        k++;
        n++;
        x0[k] = x[n];
        y0[k] = y[n];
    }
    double par0[2],par1[10];
    TGraph *gr0 = new TGraph(k+1,x0,y0);
    TF1 *f_1 = new TF1("f_1","pol1");
    gr0->Fit("f_1","Q");
    f_1->GetParameters(&par0[0]);
    double t_0,t_1;
    t_0 = (base-par0[0])/par0[1];
    //cout<<t_0<<endl;
    int n0;
    n0 = find_rise2(y,i,base,0.3*(max-base));
    //cout<<n0<<endl;
    if(isnormal(t_0)){
    for(n=0;n<i;n++){
        if(x[n]>(t_0-0.2)) break;
    }
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(0.4*(max-base)+base)) break;  //0.4
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_2 = new TF1("f_2","pol4");
    //gr1->Draw();
    gr1->Fit("f_2","Q");
    f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0) t_1=300;
    else{
        for(x2=t_0-0.15;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr){
                t_1=x2;
                break;
            }
        }
        
    }
    
    delete f_2;
    delete gr1;
    }
    delete f_1;    
    delete gr0;
    return t_1;
}
//pol_n+linear
double Fitfunction::fit_function_poln_l(double *x,double *y,int i,int base_num,double thr){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max;
    int rise_m,n,k,p,m;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    for(n=0;n<i;n++){
        if(y[n]>(0.05*(max-base)+base)) break;
    }

    for(k=0;k<i;k++){
        x0[k] = x[n];
        y0[k] = y[n];
        if(y0[k]>(0.1*(max-base)+base)) break;
        n++;
    }
    if(k==0){
        k++;
        n++;
        x0[k] = x[n];
        y0[k] = y[n];
    }
    double par0[2],par1[10];
    TGraph *gr0 = new TGraph(k+1,x0,y0);
    TF1 *f_1 = new TF1("f_1","pol1");
    gr0->Fit("f_1","Q");
    f_1->GetParameters(&par0[0]);
    double t_0,t_1;
    t_0 = (base-par0[0])/par0[1];
    //cout<<t_0<<endl;
    int n0;
    n0 = find_rise2(y,i,base,0.3*(max-base));
    //cout<<n0<<endl;
    if(isnormal(t_0)){
    for(n=0;n<i;n++){
        if(x[n]>(t_0-0.2)) break;
    }
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(0.4*(max-base)+base)) break;  //0.4
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
    TF1 *f_2 = new TF1("f_2","pol4",x1[0],x1[p]);
    //gr1->Draw();
    gr1->Fit("f_2","Q");
    f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0){
        for(m=0;m<i;m++){
            if(y[m]>(thr+base)) break;
        }
        t_1 = (thr+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    }
    else{
        for(x2=t_0-0.15;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr){
                t_1=x2;
                break;
            }
        }
        
    }
    
    delete f_2;
    delete gr1;
    }
    delete f_1;    
    delete gr0;
    return t_1;
}
//pol4+linear thr
void Fitfunction::fit_function_poln_l_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx){
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max;
	int rise_m,n,k,p,m,d;
	base = t1.find_baseline(x,y,base_num,i);  //find baseline
	max = t1.find_max(y,i);   	//find peak
	//cout<<max<<" "<<base<<endl;
	for(n=0;n<i;n++){
		if(y[n]>(0.05*(max-base)+base)) break;
	}
	
	for(k=0;k<i;k++){
		x0[k] = x[n];
		y0[k] = y[n];
		if(y0[k]>(0.1*(max-base)+base)) break;
		n++;
	}
	//cout<<k<<endl;
	if(k==0){ 
        k++;
        n++;
        x0[k] = x[n];
		y0[k] = y[n];
    }
	double par0[2],par1[10];
	TGraph *gr0 = new TGraph(k+1,x0,y0);
	TF1 *f_1 = new TF1("f_1","pol1");
	gr0->Fit("f_1","Q");
	f_1->GetParameters(&par0[0]);
	double t_0,t_1;
	t_0 = (base-par0[0])/par0[1];
	//cout<<t_0<<endl;
	int n0;
	n0 = find_rise2(y,i,base,0.3*(max-base));
    //cout<<n0<<endl;
    if(isnormal(t_0)){
    for(n=0;n<i;n++){
        if(x[n]>(t_0-0.2)) break;
    }
    for(p=0;p<i;p++){
		x1[p] = x[n];
		y1[p] = y[n];
		if(y[n]>(0.4*(max-base)+base)) break;  //0.4
		n++;
	}
	TGraph *gr1 = new TGraph(p+1,x1,y1);
	TF1 *f_2 = new TF1("f_2","pol4");
	//f_2->SetParameters(par0[1],t_0,base);
	//f_2->FixParameter(2,base);
    //gr1->Draw();
	gr1->Fit("f_2","Q");
	f_2->GetParameters(&par1[0]);
    double x2,y2;
    for(d=0;d<nx;d++){
    if(par1[1]==0&&par1[2]==0){
        for(m=0;m<i;m++){
            if(y[m]>(thr[d]+base)) break;
        }
        tx[d] = (thr[d]+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    }
    else{
        for(x2=t_0-0.15;x2<200;x2=x2+0.005){
            tools *t2 = new tools();
            y2=t2->pol_n(x2,par1,4);
            delete t2;
            if(y2>base+thr[d]){
                tx[d]=x2;
                break;
            }
        }
    }
    }
	//t_1 = thr/par1[0] + par1[1];
    //cout<<t_1<<endl;
    delete f_1;
    delete f_2;
    delete gr0;
    delete gr1;
    }
    //return t_1;	
}*/
//LET
double Fitfunction::LET_1(double *x,double *y,int i,int base_num,double thr){
    double x0[7000],y0[7000],x1[7000],y1[7000],base,max;
    int rise_m,n,k,p;
    base = t1.find_baseline(x,y,base_num,i);
    for(n=0;n<i;n++){
        if(y[n]>(base+thr)) break;
    }
    if(n==0) return x[n];
    else return x[n-1];
}
void Fitfunction::LET_cthr(double *x,double *y,int i,int base_num,double *thr,double *tx,int nx){
    double base,t_0;
    int m,n,start_c,d;
    //tools t1;
    start_c = t1.find_s_c(x,y,i,base_num);
    base = t1.find_baseline(x,y,base_num,i);
    d = 0;
    for(m=start_c;m<i;m++){
        if(y[m]>(thr[d]+base)){
            if(m==0) tx[d] = x[m];
            else tx[d] = x[m-1];
            d++;
            if(m>0) m--;
            if(d==nx) break;
        }
    }
}
//CFD
double Fitfunction::CFD_1(double *x,double *y,int i,int base_num,double attenuation,double delay){
    double x0[7000],y0[7000],x1[7000],y1[7000],y2[7000],base,max,x_0,t_0;
    int rise_m,n,k,p,m;
    base = t1.find_baseline(x,y,base_num,i);
    for(n=0;n<i;n++){
        y[n]=y[n]-base;
    }
    for(n=0;n<i;n++){
        y0[n] = attenuation * y[n];
    }
    //width
    x_0 = 0;
    for(n=0;n<10;n++){
        x_0 = x_0 + (x[n+1]-x[n]);
    }
    x_0 = x_0*0.1;
    k = (int)(delay/x_0);
    //cout<<x_0<<" "<<k<<endl;
    for(n=0;n<i;n++){
        if(n<k) y1[n] = 0;
        else y1[n] = -y[n-k];
    }
    for(n=0;n<i;n++){
        y2[n] = y0[n] + y1[n];
    }
    /*TGraph *gr1 = new TGraph(i,x,y2);
    gr1->Draw();*/
    k = t1.find_max_channel(y2,i);
    for(m=k;m<i;m++){
        if(y2[m]<0) break;
    }
    t_0 = (0-y2[m-1])*(x[m]-x[m-1])/(y2[m]-y2[m-1]) + x[m-1];
    /*for(n=0;n<i;n++){
        y1[n] = -y1[n];
    }
    double y_0;
    y_0 = (t_0 - x[m-1])*(y1[m]-y1[m-1])/(x[m]-x[m-1]) + y1[m-1];
    max = t1.find_max(y1,i);
    y_0 = y_0/max;*/
    return t_0;
}
double Fitfunction::CFD_2(double *x,double *y,int i,int base_num,double fraction){
    double base,t_0,max;
    int m,n,start_c;
    //tools t1;
    start_c = t1.find_s_c(x,y,i,base_num);
    base = t1.find_baseline(x,y,base_num,i);
    max = t1.find_max(y,i);
    for(m=start_c;m<i;m++){
        if(y[m]>(fraction*(max-base)+base)) break;
    }
    t_0 = (fraction*(max-base)+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
    return t_0;
}
void Fitfunction::CFD_2_gsl_poln_l_cthr(double *x,double *y,int i,int base_num,double *fraction,double *tx,int nx,double range,int order){
	int rise_m,n,k,p,m,d,zero_ch;
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10],chisq;
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)&&p>(order-1)) break;  //range
        n++;
    }
    //gsl_fit
    if((p+1)>order){
        gsl_matrix *x_n,*cov;
        gsl_vector *y_n,*c;
        x_n = gsl_matrix_alloc(p+1,order+1);
        y_n = gsl_vector_alloc(p+1);
        c = gsl_vector_alloc(order+1);
        cov = gsl_matrix_alloc(order+1,order+1);
        for(k=0;k < p+1;k++){
            for(m = 0;m < order+1;m++){
                gsl_matrix_set(x_n,k,m,pow(x1[k],m));
            }
            gsl_vector_set(y_n,k,y1[k]);
        }
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(p+1,order+1);
        gsl_multifit_linear(x_n,y_n,c,cov,&chisq,work);
        for(k=0;k<order+1;k++){
            par1[k] = gsl_vector_get(c,(k));
        }
        gsl_multifit_linear_free(work);
        gsl_matrix_free(x_n);
        gsl_matrix_free(cov);
        gsl_vector_free(y_n);
        gsl_vector_free(c);
    }
    else{
        par1[1]=0;
        par1[2]=0;
    }
    double x2,y2;
    if(par1[1]==0&&par1[2]==0){
        d=0;
        for(m=zero_ch;m<i;m++){
            if(y[m]>(fraction[d]*(max-base)+base)){
                tx[d] = (fraction[d]*(max-base)+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
                d++;
                m--;
                if(d==nx) break;
            }
        }
    }
    else{
        d=0;
        for(x2=t_0;x2<200;x2=x2+0.005){
            y2=t1.pol_n(x2,par1,order);
            if(y2>(fraction[d]*(max-base)+base)){
                double y2_0 = t1.pol_n(x2-0.005,par1,order);
                tx[d]=(fraction[d]*(max-base)+base-y2_0)*0.005/(y2-y2_0) + x2 - 0.005;
                d++;
                x2=x2-0.005;
                if(d==nx) break;
            }
        }
    }
}
void Fitfunction::CFD_2_gaus_l_cthr(double *x,double *y,int i,int base_num,double *fraction,double *tx,int nx,double range){
	int rise_m,n,k,p,m,d,zero_ch;
	double x0[7000],y0[7000],x1[7000],y1[7000],base,max,par1[10];
    base = t1.find_baseline(x,y,base_num,i);  //find baseline
    max = t1.find_max(y,i);   	//find peak
    zero_ch = t1.find_s_c(x,y,i,base_num);
    double t_0,t_1;
    t_0 = x[zero_ch];
    n = zero_ch;
    for(p=0;p<i;p++){
        x1[p] = x[n];
        y1[p] = y[n];
        if(y[n]>(range*(max-base)+base)) break;  //range
        n++;
    }
    TGraph *gr1 = new TGraph(p+1,x1,y1);
	TF1 *f_2 = new TF1("f_2","gaus");
	gr1->Fit("f_2","Q");
	f_2->GetParameters(&par1[0]);
    double x2,y2;
    if(par1[1]==0&&par1[2]==0){
        d=0;
        for(m=zero_ch;m<i;m++){
            if(y[m]>(fraction[d]*(max-base)+base)){
                tx[d] = (fraction[d]*(max-base)+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
                d++;
                m--;
                if(d==nx) break;
            }
        }
    }
    else{
        d=0;
        for(x2=t_0;x2<200;x2=x2+0.005){
            y2=t1.gaus(x2,par1);
            if(y2>(fraction[d]*(max-base)+base)){
                double y2_0 = t1.gaus(x2-0.005,par1);
                tx[d]=(fraction[d]*(max-base)+base-y2_0)*0.005/(y2-y2_0) + x2 - 0.005;
                d++;
                x2=x2-0.005;
                if(d==nx) break;
            }
        }
    }
    delete f_2;

    delete gr1;

}
void Fitfunction::CFD_2_spline_cthr(double *x,double *y,int i,int base_num,double *fraction,double *tx,int nx){
    //cout<<"asd"<<endl;
    double x0[50000],y0[50000],base,max,t_0;
    int nb,m,n,start_c;
    //n=3000;
    nb=i*10;
    base = t1.find_baseline(x,y,base_num,i);
    start_c = t1.find_s_c(x,y,i,base_num);
    max = t1.find_max(y,i);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    const gsl_interp_type *t = gsl_interp_cspline_periodic;
    gsl_spline *spline = gsl_spline_alloc(t,i);
    gsl_spline_init(spline,x,y,i);
    for ( Int_t k = 0; k < (nb-2); k++ )
    {
        x0[k]  = (Double_t) k*(x[i-1]-x[0])/(nb-1) + x[0];
        y0[k] = gsl_spline_eval(spline,x0[k],acc);
    }
    int d = 0;
    for(m=(10*start_c-1);m<nb-2;m++){
        if(y0[m]>(fraction[d]*(max-base)+base)){
            tx[d] = (fraction[d]*(max-base)+base-y0[m-1])*(x0[m]-x0[m-1])/(y0[m]-y0[m-1]) + x0[m-1];
            d++;
            m--;
            if(d==nx) break;
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}
void Fitfunction::CFD_2_linear_cthr(double *x,double *y,int i,int base_num,double *fraction,double *tx,int nx){
    double base,t_0,max;
    int m,n,start_c,d;
    //tools t1;
    start_c = t1.find_s_c(x,y,i,base_num);
    base = t1.find_baseline(x,y,base_num,i);
    max = t1.find_max(y,i);
    d = 0;
    for(m=start_c;m<i;m++){
        if(y[m]>(fraction[d]*(max-base)+base)){
            tx[d] = (fraction[d]*(max-base)+base-y[m-1])*(x[m]-x[m-1])/(y[m]-y[m-1]) + x[m-1];
            d++;
            m--;
            if(d==nx) break;
        }
    }
    //return t_0;
}
