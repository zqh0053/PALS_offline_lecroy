#ifndef TOOLS_H
#define TOOLS_H


class tools
{
public:
    tools();
    ~tools();
    double pol_n(double ,double *,int );  //多项式计算
    double gaus(double ,double *);
    double find_max(double *,int );       //寻找峰值幅度
    double find_min(double *,int );       //寻找峰值幅度
    int find_max_channel(double *,int);   //寻找峰值道址
    int find_min_channel(double *,int);   //寻找峰值道址
    double find_baseline(double *,double *,int);   //
    double find_baseline(double *,double *,int,int);
    double intergral(double *,double *,double ,int ,int ); //x,y,bin,i,basenum
    double intergral_s(double *,double *,double ,int ,int ,double ,double );
    void intergral_s_n(double *,double *,double ,int ,int ,double *,double ,int ,double *);
    double PSD(double *,double *,double ,int ,int ,double ,double ,double );
    int find_s_c(double *,double *,int,int);
    double getslope(double *,double *,int,int);
    double getslope_max(double *,double *,int,int);
    double getslope_mean(double *,double *,int,int);
    double getrisetime(double *,double *,int,int);
};

#endif // TOOLS_H
