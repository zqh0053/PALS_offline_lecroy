#include "runaction.h"
#include "readfile.h"
#include "fitfunction.h"
#include "tools.h"
#include "TGraph.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "unistd.h"
#include "stdio.h"
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

RunAction::RunAction()
{
	FileName[0] = "C4--C1C2_1200M--";
    FileName[1] = "C2--C1C2_1200M--";
    StopFileName = "stop.txt";
    ThrNum = 15;
    FracNum = 60;
    QsNum = 10;
    QlNum = 10;
    OutFileName = "a1.root";
    SignalStartNum = 0;
    Signal_num = 100;
    Signal_length = 1500;
    FitRange = 0.4;
    ChNum = 2;
    BaseNum = 150;
    PolNum = 4;
    Polarity = 1;
    Let = 1;
    Cft = 1;
    Circ = 0;
    for(int i = 0 ;i < FracNum;i++) Frac_v[i] = 0.01*i + 0.01;
}
RunAction::~RunAction()
{

}
string num_str(int i)
{
    char ss[20];
    sprintf(ss,"%05d",i);
    return ss;
}

void RunAction::run_csv(string name0,string name1, int n0, int n1)
{
    
}

void RunAction::run_1(){
    const char* outf_name = OutFileName.data();
    TFile *tf1 =new TFile(outf_name,"recreate");
    //TFile *tf1 =new TFile("data/test3.root","recreate");
	TTree *t1 = new TTree("t1","time tree");
    TTree::SetMaxTreeSize(100000000000LL);      
    Channel ch0,ch1,ch[10];
    //ch0
    /*t1->Branch("ch0_thr_num",&ch0.thr_num,"thr_num/I");
    t1->Branch("ch0_thr",ch0.thr,"thr[thr_num]/D");
    t1->Branch("ch0_t_0",ch0.t_0,"t_0[thr_num]/D");
    t1->Branch("ch0_t_1",ch0.t_1,"t_1[thr_num]/D");
    t1->Branch("ch0_t_2",ch0.t_2,"t_2[thr_num]/D");
    t1->Branch("ch0_L_num",&ch0.L_num,"L_num/I");
    t1->Branch("ch0_S_num",&ch0.S_num,"S_num/I");
    t1->Branch("ch0_Q_l",ch0.Q_l,"Q_l[L_num]/D");
    t1->Branch("ch0_Q_s",ch0.Q_s,"Q_s[S_num]/D");
    t1->Branch("ch0_E_l",ch0.E_l,"E_l[L_num]/D");
    t1->Branch("ch0_E_s",ch0.E_s,"E_s[S_num]/D");
    t1->Branch("ch0_max_0",&ch0.max_0,"max_0/D");
    t1->Branch("ch0_energy",&ch0.energy,"energy/D");
    t1->Branch("ch0_slope",&ch0.slope,"slope/D");
    t1->Branch("ch0_slope_max",&ch0.slope_max,"slope_max/D");
    t1->Branch("ch0_slope_mean",&ch0.slope_mean,"slope_mean/D");
    //ch1
    t1->Branch("ch1_thr_num",&ch1.thr_num,"thr_num/I");
    t1->Branch("ch1_thr",ch1.thr,"thr[thr_num]/D");
    t1->Branch("ch1_t_0",ch1.t_0,"t_0[thr_num]/D");
    t1->Branch("ch1_t_1",ch1.t_1,"t_1[thr_num]/D");
    t1->Branch("ch1_t_2",ch1.t_2,"t_2[thr_num]/D");
    t1->Branch("ch1_L_num",&ch1.L_num,"L_num/I");
    t1->Branch("ch1_S_num",&ch1.S_num,"S_num/I");
    t1->Branch("ch1_Q_l",ch1.Q_l,"Q_l[L_num]/D");
    t1->Branch("ch1_Q_s",ch1.Q_s,"Q_s[S_num]/D");
    t1->Branch("ch1_E_l",ch1.E_l,"E_l[L_num]/D");
    t1->Branch("ch1_E_s",ch1.E_s,"E_s[S_num]/D");
    t1->Branch("ch1_max_0",&ch1.max_0,"max_0/D");
    t1->Branch("ch1_energy",&ch1.energy,"energy/D");
    t1->Branch("ch1_slope",&ch1.slope,"slope/D");
    t1->Branch("ch1_slope_max",&ch1.slope_max,"slope_max/D");
    t1->Branch("ch1_slope_mean",&ch1.slope_mean,"slope_mean/D");*/
	//
    for(int i = 0;i < ChNum;i++){
        char c[80];
        snprintf(c,sizeof(c),"%d",i);
        string c_str = c;
        string Branch_name;
        Branch_name = "ch" + c_str + "_thr_num";
        t1->Branch(Branch_name.c_str(),&ch[i].thr_num,"thr_num/I");
        Branch_name = "ch" + c_str + "_thr";
        t1->Branch(Branch_name.c_str(),ch[i].thr,"thr[thr_num]/D");
        Branch_name = "ch" + c_str + "_t_0";
        t1->Branch(Branch_name.c_str(),ch[i].t_0,"t_0[thr_num]/D");
        Branch_name = "ch" + c_str + "_t_1";
        t1->Branch(Branch_name.c_str(),ch[i].t_1,"t_1[thr_num]/D");
        Branch_name = "ch" + c_str + "_t_2";
        t1->Branch(Branch_name.c_str(),ch[i].t_2,"t_2[thr_num]/D");
        Branch_name = "ch" + c_str + "_t_3";
        t1->Branch(Branch_name.c_str(),ch[i].t_3,"t_3[thr_num]/D");
        Branch_name = "ch" + c_str + "_frac_num";
        t1->Branch(Branch_name.c_str(),&ch[i].frac_num,"frac_num/I");
        Branch_name = "ch" + c_str + "_frac";
        t1->Branch(Branch_name.c_str(),ch[i].frac,"frac[frac_num]/D");
        Branch_name = "ch" + c_str + "_t_cfd0";
        t1->Branch(Branch_name.c_str(),ch[i].t_cfd0,"t_cfd0[frac_num]/D");
        Branch_name = "ch" + c_str + "_t_cfd1";
        t1->Branch(Branch_name.c_str(),ch[i].t_cfd1,"t_cfd1[frac_num]/D");
        Branch_name = "ch" + c_str + "_t_cfd2";
        t1->Branch(Branch_name.c_str(),ch[i].t_cfd2,"t_cfd2[frac_num]/D");
        Branch_name = "ch" + c_str + "_t_cfd3";
        t1->Branch(Branch_name.c_str(),ch[i].t_cfd3,"t_cfd3[frac_num]/D");
        Branch_name = "ch" + c_str + "_L_num";
        t1->Branch(Branch_name.c_str(),&ch[i].L_num,"L_num/I");
        Branch_name = "ch" + c_str + "_S_num";
        t1->Branch(Branch_name.c_str(),&ch[i].S_num,"S_num/I");
        Branch_name = "ch" + c_str + "_Q_l";
        t1->Branch(Branch_name.c_str(),ch[i].Q_l,"Q_l[L_num]/D");
        Branch_name = "ch" + c_str + "_Q_s";
        t1->Branch(Branch_name.c_str(),ch[i].Q_s,"Q_s[S_num]/D");
        Branch_name = "ch" + c_str + "_E_l";
        t1->Branch(Branch_name.c_str(),ch[i].E_l,"E_l[L_num]/D");
        Branch_name = "ch" + c_str + "_E_s";
        t1->Branch(Branch_name.c_str(),ch[i].E_s,"E_s[S_num]/D");
        Branch_name = "ch" + c_str + "_max_0";
        t1->Branch(Branch_name.c_str(),&ch[i].max_0,"max_0/D");
        Branch_name = "ch" + c_str + "_energy";
        t1->Branch(Branch_name.c_str(),&ch[i].energy,"energy/D");
        Branch_name = "ch" + c_str + "_slope";
        t1->Branch(Branch_name.c_str(),&ch[i].slope,"slope/D");
        Branch_name = "ch" + c_str + "_slope_max";
        t1->Branch(Branch_name.c_str(),&ch[i].slope_max,"slope_max/D");
        Branch_name = "ch" + c_str + "_slope_mean";
        t1->Branch(Branch_name.c_str(),&ch[i].slope_mean,"slope_mean/D");
        Branch_name = "ch" + c_str + "_risetime";
        t1->Branch(Branch_name.c_str(),&ch[i].risetime,"risetime/D");
    }
    //cout<<steeringFileName<<endl;
    //FolderName = "/media/zqh/My Passport/coin_lfsf5mm_30035_511511_33_5v/";
    string name_1 = FolderName;
    string name_2 = ".csv";
	int num_0,num_1,an;
    num_1=0;
    an = 1;
	string data_name_0,data_name_1,data_name_2;
   /* for(num_0 = 0;num_0 < SignalFileNum;num_0++){
        readfile read1;
        char c[80];
        snprintf(c,sizeof(c),"%05d",num_0);
        string name_3 = c;
        int sig_num,sig_length;
        sig_num = Signal_num;
        sig_length = Signal_length;
        //channel0
        double x0[sig_num][sig_length],y0[sig_num][sig_length];
        int i0[sig_num];
        string name_channel0 = name_1+FileName[0]+name_3+name_2;
        int file_e;
        file_e = access(name_channel0.c_str(),0);
        if(file_e == -1) break;
        read1.readfile_csv_n(name_channel0,(double*)x0,(double*)y0,i0,sig_num,sig_length);
        //channel1
        double x1[sig_num][sig_length],y1[sig_num][sig_length];
        int i1[sig_num];
        string name_channel1 = name_1+FileName[1]+name_3+name_2;
        read1.readfile_csv_n(name_channel1,(double*)x1,(double*)y1,i1,sig_num,sig_length);
        for(num_1=0;num_1<sig_num;num_1++){
            tools t2;
            Fitfunction f1;
            ch0.thr_num=ch1.thr_num=ThrNum;
            ch0.S_num=ch1.S_num=15;
            ch0.L_num=ch1.L_num=10;
            ch0.energy = t2.intergral(x0[num_1],y0[num_1],0.2,i0[num_1],150);
            ch1.energy = t2.intergral(x1[num_1],y1[num_1],0.2,i1[num_1],150);
            ch0.max_0 = t2.find_max(y0[num_1],i0[num_1]);
            ch1.max_0 = t2.find_max(y1[num_1],i1[num_1]);
            ch0.Q_l[0]=ch1.Q_l[0]=10;
            ch0.Q_s[0]=ch1.Q_s[0]=0.5;
            ch0.thr[0]=ch1.thr[0]=0.000;
            ch0.slope = t2.getslope(x0[num_1],y0[num_1],i0[num_1],150);
            ch1.slope = t2.getslope(x1[num_1],y1[num_1],i1[num_1],150);
            ch0.slope_max = t2.getslope_max(x0[num_1],y0[num_1],i0[num_1],150);
            ch1.slope_max = t2.getslope_max(x1[num_1],y1[num_1],i1[num_1],150);
            ch0.slope_mean = t2.getslope_mean(x0[num_1],y0[num_1],i0[num_1],150);
            ch1.slope_mean = t2.getslope_mean(x1[num_1],y1[num_1],i1[num_1],150);
            for(int ix = 0;ix<ch0.thr_num;ix++){
                ch0.thr[ix]=ch1.thr[ix] = ch0.thr[0] + 0.001*ix;
            }
            for(int ix = 0;ix<ch0.S_num;ix++){
                ch0.Q_s[ix]=ch1.Q_s[ix] = ch0.Q_s[0] + 0.3*ix;                
            }
            for(int ix = 0;ix<ch0.L_num;ix++){
                ch0.Q_l[ix]=ch1.Q_l[ix] = ch0.Q_l[0] + 1.0*ix;                
            }
            t2.intergral_s_n(x0[num_1],y0[num_1],0.2,i0[num_1],150,ch0.Q_s,0.005,ch0.S_num,ch0.E_s);
            t2.intergral_s_n(x0[num_1],y0[num_1],0.2,i0[num_1],150,ch0.Q_l,0.005,ch0.L_num,ch0.E_l);
            t2.intergral_s_n(x1[num_1],y1[num_1],0.2,i1[num_1],150,ch1.Q_s,0.005,ch1.S_num,ch1.E_s);
            t2.intergral_s_n(x1[num_1],y1[num_1],0.2,i1[num_1],150,ch1.Q_l,0.005,ch1.L_num,ch1.E_l);
            if(ch0.thr_num>0){
                f1.fit_function_poln_l_cthr(x0[num_1],y0[num_1],i0[num_1],150,ch0.thr,ch0.t_0,ch0.thr_num,FitRange);
                f1.fit_function_poln_l_cthr(x1[num_1],y1[num_1],i1[num_1],150,ch1.thr,ch1.t_0,ch1.thr_num,FitRange);
                f1.fit_function_gsl_poln_l_cthr(x0[num_1],y0[num_1],i0[num_1],150,ch0.thr,ch0.t_1,ch0.thr_num,FitRange,6);
                f1.fit_function_gsl_poln_l_cthr(x1[num_1],y1[num_1],i1[num_1],150,ch1.thr,ch1.t_1,ch1.thr_num,FitRange,6);
                f1.fit_function_linear_cthr(x0[num_1],y0[num_1],i0[num_1],150,ch0.thr,ch0.t_2,ch0.thr_num);
                f1.fit_function_linear_cthr(x1[num_1],y1[num_1],i1[num_1],150,ch1.thr,ch1.t_2,ch1.thr_num);
            }
            t1->Fill();
        }
        if(num_0%1000 == 0) cout<<num_0<<endl;
    }*/
    cout<<FitRange<<endl;
    cout<<FracNum<<endl;
    for(num_0 = SignalStartNum;num_0 < SignalFileNum;num_0++){
        fstream stopfile(StopFileName,ios::app|ios::in);
        string stopcommand;
        stopfile>>stopcommand;
        readfile read1;
        char c[80];
        snprintf(c,sizeof(c),"%05d",num_0);
        string name_3 = c;
        int sig_num,sig_length;
        sig_num = Signal_num;
        sig_length = Signal_length;
        double x0[10][sig_num][sig_length],y0[10][sig_num][sig_length];
        int i0[10][sig_num];
        int file_e;
        
        for(int i = 0;i < ChNum;i++){
        //channel0
            string name_channel0 = name_1+FileName[i]+name_3+name_2;
			//cout<<name_channel0<<endl;
            file_e = access(name_channel0.c_str(),0);
            if(Circ == 1){
                while(file_e == -1){
                    file_e = access(name_channel0.c_str(),0);
                    sleep(0.1);
                }
            }
            if(file_e == -1) break;
            if(Polarity == -1){
                read1.readfile_csv_n_ng(name_channel0,(double*)x0[i],(double*)y0[i],i0[i],sig_num,sig_length);
            }
            else read1.readfile_csv_n(name_channel0,(double*)x0[i],(double*)y0[i],i0[i],sig_num,sig_length);
            //read1.readfile_csv_n_test(name_channel0,(double*)x0[i],(double*)y0[i],i0[i],sig_num,sig_length);
        }
        if(file_e == -1) break;
        if(stopcommand == "stop") break;
        for(num_1=0;num_1<sig_num;num_1++){
            for(int i = 0;i < ChNum;i++){
                tools t2;
                Fitfunction f1;
                ch[i].thr_num = ThrNum;
                ch[i].frac_num = FracNum;
                ch[i].S_num=QsNum;
                ch[i].L_num=QlNum;
                ch[i].energy = t2.intergral(x0[i][num_1],y0[i][num_1],0.2,i0[i][num_1],BaseNum);
                ch[i].max_0 = t2.find_max(y0[i][num_1],i0[i][num_1]);
                ch[i].Q_l[0]=10;
                ch[i].Q_s[0]=0.5;
                ch[i].thr[0]=0.000;
                ch[i].slope = t2.getslope(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum);
                ch[i].slope_max = t2.getslope_max(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum);
                ch[i].slope_mean = t2.getslope_mean(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum);
                ch[i].risetime = t2.getrisetime(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum);
                for(int ix = 0;ix<ch[i].thr_num;ix++){
                    ch[i].thr[ix] = ch[i].thr[0] + 0.001*ix;
                }
                for(int ix = 0;ix<ch[i].frac_num;ix++){
                    ch[i].frac[ix] = Frac_v[ix];
                }
                for(int ix = 0;ix<ch[i].S_num;ix++){
                    ch[i].Q_s[ix] = ch[i].Q_s[0] + 0.3*ix;                
                }
                for(int ix = 0;ix<ch[i].L_num;ix++){
                    ch[i].Q_l[ix] = ch[i].Q_l[0] + 1.0*ix;                
                }
                t2.intergral_s_n(x0[i][num_1],y0[i][num_1],0.2,i0[i][num_1],BaseNum,ch[i].Q_s,0.003,ch[i].S_num,ch[i].E_s);
                t2.intergral_s_n(x0[i][num_1],y0[i][num_1],0.2,i0[i][num_1],BaseNum,ch[i].Q_l,0.003,ch[i].L_num,ch[i].E_l);
                if(Let == 1){
                    f1.fit_function_gsl_poln_l_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].thr,ch[i].t_0,ch[i].thr_num,FitRange,PolNum);
                    f1.fit_function_gaus_l_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].thr,ch[i].t_1,ch[i].thr_num,FitRange);
                    f1.fit_function_linear_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].thr,ch[i].t_2,ch[i].thr_num);
                    //f1.LET_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].thr,ch[i].t_2,ch[i].thr_num);
                    f1.fit_function_spline_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].thr,ch[i].t_3,ch[i].thr_num);
                }
                //cout<<num_1<<endl;
                if(Cft == 1){
					//cout<<"asda"<<endl;
                    f1.CFD_2_gsl_poln_l_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].frac,ch[i].t_cfd0,ch[i].frac_num,FitRange,PolNum);
					//cout<<"asda"<<endl;
                    f1.CFD_2_gaus_l_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].frac,ch[i].t_cfd1,ch[i].frac_num,FitRange);
                    f1.CFD_2_linear_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].frac,ch[i].t_cfd2,ch[i].frac_num);
                    //f1.CFD_2_spline_cthr(x0[i][num_1],y0[i][num_1],i0[i][num_1],BaseNum,ch[i].frac,ch[i].t_cfd3,ch[i].frac_num);
                }
            }
            t1->Fill();
        }
        stopfile.close();
        if(num_0%100 == 0) cout<<num_0<<endl;
    }
    tf1->Write();
    sleep(5);
    remove(StopFileName.c_str());
}
