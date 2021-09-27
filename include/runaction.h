#ifndef RUNACTION_H
#define RUNACTION_H
#include <string>
using namespace std;
const int Max_num = 100;
class RunAction
{
public:
    RunAction();
    ~RunAction();
    void run_csv(string ,string ,int ,int );
    void run_binary(char*,int);
    void run_1();
    string test_st;
    
    //parameters
    string FolderName;
    string OutFileName;
    string StopFileName;
    string FileName[10];
    int Signal_num;
    int SignalStartNum;
    int Signal_length;
    int SignalFileNum;
    int ThrNum;
    int FracNum;
    int QsNum;
    int QlNum;
    int ChNum;
    int BaseNum;
    int PolNum;
    int Polarity;
    int Let;
    int Cft;
    int Circ;
    double FitRange;
    double Frac_v[1000];
};
typedef struct {
    int thr_num;
    double thr[Max_num];
    double t_0[Max_num];  //poln
    double t_1[Max_num];  //gaus
    double t_2[Max_num];  //linear
    double t_3[Max_num];  //spline
    int frac_num;
    double frac[Max_num];
    double t_cfd0[Max_num]; //cfd
    double t_cfd1[Max_num];
    double t_cfd2[Max_num];
    double t_cfd3[Max_num];
    int L_num;
    double Q_l[Max_num];
    double E_l[Max_num];
    int S_num;
    double Q_s[Max_num];      
    double E_s[Max_num];
    double max_0;
    double energy;
    double slope;
    double slope_max;
    double slope_mean;
    double risetime;
} Channel;
#endif // RUNACTION_H
