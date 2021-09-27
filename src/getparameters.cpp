#include "getparameters.h"
#include "runaction.h"
#include "readfile.h"
#include "fitfunction.h"
#include "tools.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unistd.h>
using namespace std;
getparameters::getparameters(RunAction *rac)
:frunaction(rac)
{
    
}
getparameters::~getparameters(){
    
}
void getparameters::readfile(string name){
    const char *p = name.data();
    fstream file1(name,ios::in);
    string line;
    while(getline(file1,line)){
        istringstream sin(line);
        string n1;
        sin>>n1;
        if(n1 == "FolderName"){
            string fdname[10];
            int fn_num = 0;
            while(!sin.eof()){
                sin>>fdname[fn_num];
                fn_num++;
            }
            frunaction->FolderName = fdname[0];
            for(int i = 1;i<fn_num;i++) frunaction->FolderName = frunaction->FolderName + " " + fdname[i];
            cout<<frunaction->FolderName<<endl;
        }
        else if(n1 == "OutFileName"){
            sin>>frunaction->OutFileName;
        }
        else if(n1 == "StopFileName"){
            sin>>frunaction->StopFileName;
        }
        else if(n1 == "ChNum"){
            sin>>frunaction->ChNum;
        }
        else if(n1 == "FileName"){
            int i = 0;
            while(!sin.eof()){
                sin>>frunaction->FileName[i];
                i++;
            }
        }
        else if(n1 == "Signal"){
            sin>>frunaction->SignalFileNum>>frunaction->Signal_num>>frunaction->Signal_length;
        }
        else if(n1 == "SignalStartNum"){
            sin>>frunaction->SignalStartNum;
        }
        else if(n1 == "ThrNum"){
            sin>>frunaction->ThrNum;
        }
        else if(n1 == "FracNum"){
            sin>>frunaction->FracNum;
        }
        else if(n1 == "QsNum"){
            sin>>frunaction->QsNum;
        }
        else if(n1 == "QlNum"){
            sin>>frunaction->QlNum;
        }
        else if(n1 == "FitRange"){
            sin>>frunaction->FitRange;
        }
        else if(n1 == "BaseNum"){
            sin>>frunaction->BaseNum;
        }
        else if(n1 == "PolNum"){
            sin>>frunaction->PolNum;
        }
        else if(n1 == "Polarity"){
            sin>>frunaction->Polarity;
        }
        else if(n1 == "Let"){
            sin>>frunaction->Let;
        }
        else if(n1 == "Cft"){
            sin>>frunaction->Cft;
        }
        else if(n1 == "Circ"){
            sin>>frunaction->Circ;
        }
        else if(n1 == "FracValue"){
            int i = 0;
            while(!sin.eof()){
                sin>>frunaction->Frac_v[i];
                i++;
            }
        }
    }
}
