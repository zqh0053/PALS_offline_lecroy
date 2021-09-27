#ifndef FITFUNCTION_H
#define FITFUNCTION_H 1
#include "tools.h"
class Fitfunction
{
	public:
    //double fit_function_1(double*,double*,int ,int );
	Fitfunction();
	~Fitfunction();
    double fit_function_pol1(double*,double*,int ,int ,double);    //波形线性拟合
    double fit_function_linear(double*,double*,int ,int ,double);  //波形线性插值
    void fit_function_linear_cthr(double*,double*,int ,int ,double* ,double*,int); 
    double fit_function_spline(double*,double*,int ,int ,double);  //波形样条插值
    void fit_function_spline_cthr(double*,double*,int ,int ,double* ,double*,int); 
    double fit_function_pol2(double*,double*,int ,int ,double);    //二次拟合
    double fit_function_poln(double*,double*,int ,int ,double,double);    //高阶多项式拟合
    double fit_function_poln_l(double*,double*,int ,int ,double,double);    //高阶多项式拟合
    void fit_function_poln_l_cthr(double*,double*,int ,int ,double* ,double*,int,double); 
    void fit_function_gsl_poln_l_cthr(double*,double*,int ,int ,double* ,double*,int,double,int); 
    double fit_function_gaus(double*,double*,int ,int ,double,double);    //gaus
    void fit_function_gaus_l_cthr(double*,double*,int ,int ,double* ,double*,int,double); 
    double fit_function_gaus_l(double*,double*,int ,int ,double,double);
    double CFD_1(double *,double *,int ,int ,double ,double );  //CFD
    double CFD_2(double *,double *,int ,int ,double );  //CFD
    void CFD_2_gsl_poln_l_cthr(double*,double*,int ,int ,double* ,double*,int,double,int);
    void CFD_2_linear_cthr(double*,double*,int ,int ,double* ,double*,int);
    void CFD_2_spline_cthr(double*,double*,int ,int ,double* ,double*,int); 
    void CFD_2_gaus_l_cthr(double*,double*,int ,int ,double* ,double*,int,double); 
    double LET_1(double *,double *,int ,int ,double );
    void LET_cthr(double *,double *,int ,int ,double *,double *,int);
private:
    tools t1;
};
#endif
