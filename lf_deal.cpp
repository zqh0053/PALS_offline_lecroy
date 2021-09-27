#include "fitfunction.h"
#include "readfile.h"
#include "getparameters.h"
#include "runaction.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <unistd.h>
#include <iostream>
#include <ctime>
using namespace std;
int main(int argc,char** argv)
{   
    clock_t startT,stopT;
    startT = clock();
    string steeringFileName=argv[1];
    RunAction *r1 = new RunAction; 
    getparameters *gp1 = new getparameters(r1);
    if(access(argv[1],0) == -1){
        cout<<"InputFile is not exist"<<endl;
    }
    else{
        
        gp1->readfile(steeringFileName);
        r1->run_1();
        
    }
    delete r1;
    delete gp1;
    stopT = clock();
    cout<<(double)(stopT-startT)/CLOCKS_PER_SEC<<endl;
    return 0;
}
