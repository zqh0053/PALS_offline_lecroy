#ifndef READFILE_H
#define READFILE_H
#include <string>
using namespace std;

class readfile
{
public:
    readfile();
    ~readfile();
    void readfile_csv(string ,double *,double *,int *);
    void readfile_csv_1(string ,double *,double *,int *);
    void readfile_csv_ng(string ,double *,double *,int *);
    void readfile_csv_1_ng(string ,double *,double *,int *);
    void readfile_csv_n(string ,double *,double *,int *,int ,int );
    void readfile_csv_n_test(string ,double *,double *,int *,int ,int );
    void readfile_csv_n_ng(string ,double *,double *,int *,int ,int );
    void readfile_csv_n1(string name,double (*x0)[1200],double (*y0)[1200],int *i0,int sig_num,int length);
};

#endif // READFILE_H
