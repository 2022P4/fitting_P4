/*
        tdcを用いたToFの計算
*/
#include <stdlib.h>

#include <iostream>
#define N 256
int push(void);

int ToF(int RUNNUMBER) {
    // ここは共通
    gROOT->Reset();
    gStyle->SetOptStat(1001110);
    //　本番用
    // int RUNNUMBER;
    TFile *f = TFile::Open(Form("/mnt/p4root_2022/%04d.root", RUNNUMBER), "R");
    TFile *g = TFile::Open(Form("./%04d_new.root", RUNNUMBER), "update");
    /*----------------------------------------------------------------------------*/

    int tdcraw[16];
    int sintiCh[7], sintiCh_s;
    int tdc[7];
    Double_t tdcsec[7];
    Double_t tdcsa[7];

    int evenum;
    int tdc_ch;
    double slope[7], slope_s, intercept[7], inter_s;
    //読み込むプログラム-------------------------------------------------
    FILE *fp;
    char fname[] = "./input/tdc_fit_para.dat";
    fp = fopen(fname, "r");
    if (fp == NULL) printf("%s file not open!\n", fname);
    char header[N];
    fgets(header, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &tdc_ch, &sintiCh_s, &slope_s, &inter_s) !=
           EOF) {
        //詰め込み
        // ch1-6 is  sintiCh 1-6, RF is sintiCh0
        slope[sintiCh_s] = slope_s;
        intercept[sintiCh_s] = inter_s;
        //確認
        printf("%d %f %f\n", tdc_ch, slope[sintiCh_s], intercept[sintiCh_s]);
    }
    fclose(fp);
    /*----------------------------------------------------------------*/
    TTree *oldtree = (TTree *)f->Get("tree");
    // deifine oldtree
    oldtree->SetBranchAddress("tdcraw", tdcraw);
    // set Br old info
    TTree *tree = new TTree("tree", "tree");
    /* Branchの設定 */
    tree->Branch("sintiCh", sintiCh, "sintiCh[7]/I");
    tree->Branch("tdc", tdc, "tdc[7]/I");
    tree->Branch("tdcsec", tdcsec, "tdcsec[7]/D");
    tree->Branch("tdcsa", tdcsa, "tdcsa[7]/D");

    int LOOP = oldtree->GetEntries();
    /*	データの詰め込み　*/
    for (int i = 0; i < LOOP; i++) {
        oldtree->GetEntry(i);
        evenum = i;
        for (int j = 0; j < 7; j++) {
            sintiCh_s = (j + 1) % 7;
            // if (i == 1) cout << sintiCh_s << endl;
            sintiCh[sintiCh_s] = sintiCh_s;
            tdc[sintiCh_s] = tdcraw[j];
            // sintiの番号を対応させる
            tdcsec[sintiCh_s] =
                slope[sintiCh_s] * (double)tdcraw[j] + intercept[sintiCh_s];
            // ch->time[ns]
            tdcsa[sintiCh_s] = tdcsec[0] - tdcsec[sintiCh_s];
        }
        tree->Fill();
    }
    tree->Write();
    return 0124;
}
