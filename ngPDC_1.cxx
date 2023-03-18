/*
    三次関数でng弁別, CFD較正
    my_analyze ディレクトリ から解析

*/

#include <stdlib.h>
#define N 256

void ngPDC_1(int RUNNUMBER){
    
    gROOT->Reset();
    gStyle->SetOptStat(1001110);

    TFile *f = TFile::Open(Form("./../root/%04d.root", RUNNUMBER), "R");
    TFile *g = TFile::Open(Form("./../root/%04d_ng1.root", RUNNUMBER), "recreate");

    //読み込みハコ
    int qdcraw[16]={}; //oldtree qdcraw
    int tdcraw[16]={}; //oldtree tdcraw
    double gain[7]={};
    double pedestal[16]={};
    double p0[7]={};
    double p1[7]={};
    double p2[7]={};
    double p3[7]={};  
    double a[7]={};
    double b[7]={};
    double c[7]={};
    double slope[9]={};
    double intercept[9]={};
    //書き込みハコ  tree の branch になる
    int allch[7]={0,1,3,5,9,11,13};   //qdc ch - 液シン の対応
    int fastch[7]={0,2,4,8,10,12,14};  //qdc ch - 液シン の対応
    double qdc[16]={};      //tree
    double all[7]={};       //tree
    double fast[7]={};      //tree
    //double qdc_ee[16]={};   //tree
    double all_ee[7]={};     //tree
    double fast_ee[7]={};    //tree
    int tdc[16]={};          //tree
    double tdc_calib[16]={}; //tree
    double function[7]={};   //tree
    double ngpsd[7]={};      //tree
    Double_t tdcsec[9]={};      //tree
    Double_t tdcvalue[9]={};    //tree


    // file reading ---------------------------------------------------------------
    FILE *fp;
    char fname1[] = "./../input/qdc_info_20230228.dat";
    fp = fopen(fname1, "r");
    if (fp == NULL) cout << fname1 << " file not open!" << endl;
    char header1[N];
    int qdc_ch, ch_s;
    int ch[16];
    double pedestal_s;
    double gainn_s;
    fgets(header1, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &qdc_ch, &ch_s, &pedestal_s, &gainn_s) !=
           EOF) {
        // 詰め込み
        ch[qdc_ch] = ch_s;
        pedestal[qdc_ch] = pedestal_s;
        gain[qdc_ch] = gainn_s;
    }
    fclose(fp);

    char fname2[] = "./../input/nfunc_info_20230228.dat";
    fp = fopen(fname2, "r");
    if (fp == NULL) cout << fname2 << " file not open!" << endl;
    char header2[N];
    double p0_s, p1_s, p2_s, p3_s;
    fgets(header2, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%lf%lf%lf%lf", &ch_s, &p0_s, &p1_s, &p2_s, &p3_s) !=
           EOF) {
        // 詰め込み
        p0[ch_s] = p0_s;
        p1[ch_s] = p1_s;
        p2[ch_s] = p2_s;
        p3[ch_s] = p3_s;
    }
    fclose(fp);

    char fname3[] = "./../input/cfd_calib_para_20230228.dat";
    fp = fopen(fname3, "r");
    if (fp == NULL) cout << fname3 << " file not open!" << endl;
    char header3[N];
    double a_s, b_s, c_s;
    fgets(header3, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%lf%lf%lf", &ch_s, &a_s, &b_s, &c_s) !=
           EOF) {
        // 詰め込み
        a[ch_s] = a_s;
        b[ch_s] = b_s;
        c[ch_s] = c_s;
    }
    fclose(fp);

    char fname4[] = "./../input/tdc_fit_para.dat";
    fp = fopen(fname4, "r");
    if (fp == NULL) printf("%s file not open!\n", fname4);
    char header4[N];
    int tdc_ch, sCh_t;
    double slope_s, inter_s;
    fgets(header4, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &tdc_ch, &sCh_t, &slope_s, &inter_s) !=
           EOF) {
        // 詰め込み
        //  ch1-6 is  sCh 1-6, RF is sCh0
        slope[sCh_t] = slope_s;
        intercept[sCh_t] = inter_s;
    }
    fclose(fp);

    //確認用　//ちゃんと読み込めていることを確認
    //for(int i=0; i<7; i++) { cout << i << " : " << p2[i] << " " << p3[i] << endl; }

    // tree setting ------------------------------------------------------------------
    TTree *oldtree = (TTree *)f->Get("tree");
    oldtree->SetBranchAddress("qdcraw", qdcraw);
    oldtree->SetBranchAddress("tdcraw", tdcraw);
    //new tree
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("qdc", qdc, "qdc[16]/D");
    tree->Branch("all", all, "all[7]/D");
    tree->Branch("fast", fast, "fast[7]/D");
    //tree->Branch("qdc_ee", qdc_ee, "qdc_ee[16]/D");
    tree->Branch("all_ee", all_ee, "all_ee[7]/D");
    tree->Branch("fast_ee", fast_ee, "fast_ee[7]/D");
    tree->Branch("tdc", tdc, "tdc[9]/I");
    tree->Branch("tdc_calib", tdc_calib, "tdc_calib[9]/D");
    tree->Branch("function", function, "function[7]/D");
    tree->Branch("ngpsd", ngpsd, "ngpsd[7]/D");
    tree->Branch("tdcsec", tdcsec, "tdcsec[9]/D");
    tree->Branch("tdcvalue", tdcvalue, "tdcvalue[9]/D");

    
    // data filing ---------------------------------------------
    int LOOP = oldtree->GetEntries();
    for(int i=0; i<LOOP; i++){
        oldtree->GetEntry(i);

        for(int j=0; j<16; j++){
            qdc[j] = qdcraw[j] - pedestal[j];
            tdc[j] = tdcraw[j];
            tdc_calib[j] = tdcraw[j];
        }
        for(int j=0; j<7; j++){
            // all & fast
            all[j] = qdc[allch[j]];
            fast[j] = qdc[fastch[j]];
            all_ee[j] = qdc[allch[j]]/gain[j];
            fast_ee[j] = qdc[fastch[j]]/gain[j];

            // function
            function[j] = p0[j] + p1[j]*all_ee[j] + p2[j]*pow(all_ee[j],2) + p3[j]*pow(all_ee[j],3) ;
            ngpsd[j] = fast_ee[j] - function[j];
            //ngpsd_f[j] = ngpsd[j]/function[j];

            //CFD calibration
            double diff = c[j] - (a[j]*all[j] + b[j]);
            if(diff >0 && tdc[j]!=0){
                tdc_calib[j] = tdc[j] + diff;
                tdc_calib[7] = tdc[7] + diff;
                tdc_calib[8] = tdc[8] + diff;
            }

        }
        for(int i=0; i<9; i++){
            tdcsec[i] = ((double)tdc_calib[i] + intercept[i]) / slope[i];

            if ((tdcsec[i] < 59 && tdcsec[i] > 54)|| i==7 || i== 8)
            {
                // RF(ch7)-tdcsec [ns]
                tdcvalue[i] = (double)(tdcsec[7] - tdcsec[i]);

                //if (tdcvalue[i] < 0) {tdcvalue[i] = 0.;}
            }
            else {
                tdcvalue[i] = -1000.;
            }
        }

        tree->Fill();
    }
    tree->Write();

    //確認用
    //for(int i=0; i<16; i++) { cout << i << " : " << pedestal[i] << " " << gain[i] << endl; }
}