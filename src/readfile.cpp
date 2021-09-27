#include "readfile.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include <typeinfo>
#include <iostream>
#include <vector>
#include <sstream>
using namespace std;
readfile::readfile()
{

}

readfile::~readfile()
{

}

void readfile::readfile_csv(string name,double *x0,double *y0,int *i0)//读取单波形CSV文件。
{
    fstream file1(name,ios::in);
        //cout<<name<<endl;
        int mm=0,a,b;
        double x1[7000],y1[7000];
        a = *i0;

        string line;
        double z0[70000];
        while(getline(file1,line)){
            mm++;
            if(mm==14) break;  //5
        }
        int k=0;
        while(getline(file1,line)){
            istringstream sin(line);
            vector<string> fields;
            string field;
            int j = 0;
            string as[2];
            while (getline(sin, field, ',')){  //除去分割逗号
                //fields.push_back(field);
                //cout<<field<<endl;
                as[j] = field;
                //cout<<as[j]<<endl;
                j++;
            }
            stringstream sstr1(as[0]);
            stringstream sstr2(as[1]);
            sstr1>>z0[a];
            x0[a]=z0[a]*1000000000;
            sstr2>>y0[a];
            a++;

        }
        *i0=a;
}

void readfile::readfile_csv_1(string name,double *x0,double *y0,int *i0){//读取多波形csv文件，二维数组须转为一维数组
    fstream file1(name,ios::in);
    int mm=0,a[40]={0},b=0;

    string line;
    while(getline(file1,line)){
        mm++;
        if(mm==44) break;
    }
    int k=0,n=0,p=0;
    double z[60000];
    while(getline(file1,line)){
        istringstream sin(line);
        vector<string> fields;
        string field;
        int j = 0;
        string as[2];
        while (getline(sin, field, ',')){
            as[j] = field;
            j++;
        }
        stringstream sstr1(as[0]);
        stringstream sstr2(as[1]);
        sstr1>>z[p];
        if(p>0){
            if(z[p]<z[p-1]){
                n++;
            }
        }
        *(x0+n*3000+a[n])=z[p]*1000000000;
        sstr2>>*(y0+n*3000+a[n]);
        a[n]++;
        p++;
        //cout<<x[i]<<" "<<y[i]<<endl;

    }
    int as;
    for(as=0;as<40;as++){
        i0[as]=a[as];
    }

}
void readfile::readfile_csv_ng(string name,double *x0,double *y0,int *i0)//读取单波形CSV文件。
{
    fstream file1(name,ios::in);
        //cout<<name<<endl;
        int mm=0,a,b;
        double x1[7000],y1[7000];
        a = *i0;

        string line;
        double z0[70000];
        while(getline(file1,line)){
            mm++;
            if(mm==14) break;  //5
        }
        int k=0;
        while(getline(file1,line)){
            istringstream sin(line);
            vector<string> fields;
            string field;
            int j = 0;
            string as[2];
            while (getline(sin, field, ',')){  //除去分割逗号
                //fields.push_back(field);
                //cout<<field<<endl;
                as[j] = field;
                //cout<<as[j]<<endl;
                j++;
            }
            stringstream sstr1(as[0]);
            stringstream sstr2(as[1]);
            sstr1>>z0[a];
            x0[a]=z0[a]*1000000000;
            double yy;
            sstr2>>yy;
            y0[a]=-yy;
            a++;

        }
        *i0=a;
}
void readfile::readfile_csv_1_ng(string name,double *x0,double *y0,int *i0){//读取多波形csv文件，二维数组须转为一维数组
    fstream file1(name,ios::in);
    int mm=0,a[40]={0},b=0;

    string line;
    while(getline(file1,line)){
        mm++;
        if(mm==44) break;
    }
    int k=0,n=0,p=0;
    double z[60000];
    while(getline(file1,line)){
        istringstream sin(line);
        vector<string> fields;
        string field;
        int j = 0;
        string as[2];
        while (getline(sin, field, ',')){
            as[j] = field;
            j++;
        }
        stringstream sstr1(as[0]);
        stringstream sstr2(as[1]);
        sstr1>>z[p];
        if(p>0){
            if(z[p]<z[p-1]){
                n++;
            }
        }
        *(x0+n*3000+a[n])=z[p]*1000000000;
        double yy;
        sstr2>>yy;
        *(y0+n*3000+a[n]) = -yy;
        a[n]++;
        p++;
        //cout<<x[i]<<" "<<y[i]<<endl;

    }
    int as;
    for(as=0;as<40;as++){
        i0[as]=a[as];
    }

}
void readfile::readfile_csv_n(string name,double *x0,double *y0,int *i0,int sig_num,int length){//读取n波形csv文件，二维数组须转为一维数组,波形数
    
    fstream file1(name,ios::in);
    int mm=0,a[200]={0},b=0;
    string line;
    while(getline(file1,line)){
        mm++;
        if(mm==sig_num+4) break;
    }
    int k=0,n=0,p=0,nump;
    nump = sig_num*length;
    double z[nump];
    while(getline(file1,line)){
        istringstream sin(line);
        vector<string> fields;
        string field;
        int j = 0;
        string as[2];
        while (getline(sin, field, ',')){
            as[j] = field;
            j++;
        }
        stringstream sstr1(as[0]);
        stringstream sstr2(as[1]);
        sstr1>>z[p];
        if(p>0){
            if(z[p]<z[p-1]){
                n++;
            }
        }
        *(x0+n*length+a[n])=z[p]*1000000000;
        sstr2>>*(y0+n*length+a[n]);
        a[n]++;
        p++;
        //cout<<x[i]<<" "<<y[i]<<endl;

    }
    int as;
    for(as=0;as<sig_num;as++){
        i0[as]=a[as];
    }

}
void readfile::readfile_csv_n_ng(string name,double *x0,double *y0,int *i0,int sig_num,int length){//读取n波形csv文件，二维数组须转为一维数组,波形数
    
    fstream file1(name,ios::in);
    int mm=0,a[200]={0},b=0;
    string line;
    while(getline(file1,line)){
        mm++;
        if(mm==sig_num+4) break;
    }
    int k=0,n=0,p=0,nump;
    nump = sig_num*length;
    double z[nump];
    while(getline(file1,line)){
        istringstream sin(line);
        vector<string> fields;
        string field;
        int j = 0;
        string as[2];
        while (getline(sin, field, ',')){
            as[j] = field;
            j++;
        }
        stringstream sstr1(as[0]);
        stringstream sstr2(as[1]);
        sstr1>>z[p];
        if(p>0){
            if(z[p]<z[p-1]){
                n++;
            }
        }
        *(x0+n*length+a[n])=z[p]*1000000000;
        double yy;
        sstr2>>yy;
        *(y0+n*length+a[n]) = -yy;
        a[n]++;
        p++;
        //cout<<x[i]<<" "<<y[i]<<endl;

    }
    int as;
    for(as=0;as<sig_num;as++){
        i0[as]=a[as];
    }

}
void readfile::readfile_csv_n_test(string name,double *x0,double *y0,int *i0,int sig_num,int length){//读取n波形csv文件，二维数组须转为一维数组,波形数
    
    fstream file1(name,ios::in);
    int mm=0,a[200]={0},b=0;
    string line;
    while(getline(file1,line)){
        mm++;
        if(mm==sig_num+4) break;
    }
    int k=0,n=0,p=0,nump;
    nump = sig_num*length;
    double z[nump];
    int m_test = 0;
    while(getline(file1,line)){
        if(m_test%4 == 0){
        istringstream sin(line);
        vector<string> fields;
        string field;
        int j = 0;
        string as[2];
        while (getline(sin, field, ',')){
            as[j] = field;
            j++;
        }
        stringstream sstr1(as[0]);
        stringstream sstr2(as[1]);
        sstr1>>z[p];
        if(p>0){
            if(z[p]<z[p-1]){
                n++;
            }
        }
        *(x0+n*length+a[n])=z[p]*1000000000;
        sstr2>>*(y0+n*length+a[n]);
        a[n]++;
        p++;
        }
        m_test++;
        //cout<<x[i]<<" "<<y[i]<<endl;

    }
    int as;
    for(as=0;as<sig_num;as++){
        i0[as]=a[as];
    }

}
