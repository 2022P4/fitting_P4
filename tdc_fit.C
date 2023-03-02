/*
        qdcを用いたn/gammaの識別
        rootfileを立ち上げる
*/
#include <math.h>
#include <stdlib.h>

#include <iostream>
#define N 256
#define Mn 939.56
#define L 39.912       //[m]
#define c 0.299792458  //[m/ns]

double Ene(double time);
double Beta(double time);

void tdc_fit(int RUNNUMBER) {
    // int RUNNUMBER = 115;
    gROOT->Reset();
    gStyle->SetOptStat(1001110);
    // 　本番用
    //  int RUNNUMBER = 150;

    TFile *f = TFile::Open(Form("./root/%04d_RF.root", RUNNUMBER), "R");
    TFile *g = TFile::Open(Form("./root/%04d_fit.root", RUNNUMBER),"recreate");
    //  sample test
    //TFile *f = TFile::Open("./sample_new.root", "R");
    //TFile *g = TFile::Open("./sample_new2.root", "recreate");
    /*	wada 用
    TFile *f = TFile::Open("./../../root/0070.root","R");
    TFile *g = TFile::Open("./0070_new.root","recreate");
    */

    int IsNeutron[16];
    int IsHit[9];
    double tdcsa_FA[9];
    double tdcsa_SA[9];
    // new area
    double tof_g = L/c;  //[ns]
    // tdc_gamma - tdc_n
    double tdc_gn_FA[7];
    double trf = 60;  //[ns]
    double ToFnraw_3_FA[7];
    double Tn_3_FA[7];

    double tdc_gn_SA[7];
    double ToFnraw_3_SA[7];
    double Tn_3_SA[7];
    double beta_FA[7];
    double beta_SA[7];
    // double ToFnraw_4[7];
    // double Tn_4[7];
    int Num = 0;
    // reset the data
    for (int i = 0; i < 7; i++) {
        tdcsa_FA[i] = 0.;
        tdc_gn_FA[i] = 0.;
        tdcsa_SA[i] = 0.;
        tdc_gn_SA[i] = 0.;
        ToFnraw_3_FA[i] = -1.;
        ToFnraw_3_SA[i] = -1.;
        beta_FA[i] = -1.;
        beta_SA[i] = -1.;
        Tn_3_SA[i] = -1;
        Tn_3_FA[i] = -1;
    }
    TTree *oldtree = (TTree *)f->Get("tree");

    oldtree->SetBranchAddress("IsNeutron", IsNeutron);
    oldtree->SetBranchAddress("IsHit", IsHit);
    oldtree->SetBranchAddress("tdcsa_FA", tdcsa_FA);
    oldtree->SetBranchAddress("tdcsa_SA", tdcsa_SA);

    TTree *tree = new TTree("tree", "tree");
    TTree *Etree = new TTree("Etree", "Etree");
    // old area =================================================//
    tree->Branch("IsNeutron", IsNeutron, "IsNeutron[7]/I");
    tree->Branch("IsHit", IsHit, "IsHit[7]/I");
    tree->Branch("tdcsa_FA", tdcsa_FA, "tdcsa_FA[7]/D");
    tree->Branch("tdcsa_SA", tdcsa_SA, "tdcsa_SA[7]/D");
    // new area =================================================//
    tree->Branch("tdc_gn_FA", tdc_gn_FA, "tdc_gn_FA[7]/D");
    tree->Branch("ToFnraw_3_FA", ToFnraw_3_FA, "ToFnraw_3_FA[7]/D");
    tree->Branch("Tn_3_FA", Tn_3_FA, "Tn_3_FA[7]/D");

    tree->Branch("tdc_gn_SA", tdc_gn_SA, "tdc_gn_SA[7]/D");
    tree->Branch("ToFnraw_3_SA", ToFnraw_3_SA, "ToFnraw_3_SA[7]/D");
    tree->Branch("Tn_3_SA", Tn_3_SA, "Tn_3_SA[7]/D");
    // Etree
    Etree->Branch("IsNeutron", IsNeutron, "IsNeutron[7]/I");
    Etree->Branch("IsHit", IsHit, "IsHit[7]/I");
    Etree->Branch("beta_FA", beta_FA, "beta_FA[7]/D");
    Etree->Branch("Tn_3_FA", Tn_3_FA, "Tn_3_FA[7]/D");
    Etree->Branch("beta_SA", beta_SA, "beta_SA[7]/D");
    Etree->Branch("Tn_3_SA", Tn_3_SA, "Tn_3_SA[7]/D");

    /* fitting area FA------------------------------------------------------*/
    TH1 *oldhist_FA[7];  // [0] is not used
    double MINCH = -2000.;
    double MAXCH = 2000.;
    // make the histgram
    for (int i = 1; i < 7; i++) {
        oldhist_FA[i] = new TH1F(
            Form("tdc_sa_FA %3d_sinti ch%d", RUNNUMBER, i),
            Form("#%3d double  gamma gaus fitting qdc ch%d;tdcsa[ns];Count",
                 RUNNUMBER, i),
            1000, MINCH, MAXCH);
    }
    // fill histgram
    for (int i = 1; i < 7; i++) {
        oldtree->Draw(
            Form("tdcsa_FA[%d]>>tdc_sa_FA %3d_sinti ch%d", i, RUNNUMBER, i),
            Form("IsNeutron[%d]==0", i));
    }
    // set fitting parameter
    // 変更があるならこの場所
    // ---------------------------------------------------------------------------------------//
    double init_mu_FA[7];
    init_mu_FA[0] = 0;
    for (int ch = 1; ch < 7; ch++) {
        init_mu_FA[ch] = 30.;
    }
    double range = 10.;
    // ---------------------------------------------------------------------------------------//
    // make gaus+constant
    TF1 *gaus111[7];
    for (int i = 1; i < 7; i++) {
        gaus111[i] = new TF1(Form("gaus%3d_%d", RUNNUMBER, i),
                             "[0]*exp(-0.5*pow((x-[1])/[2],2.0))+[4]",
                             init_mu_FA[i] - range, init_mu_FA[i] + range);

        gaus111[i]->SetParameters(1000., init_mu_FA[i], 5., 100.);
        gaus111[i]->SetNpx(1000);
    }

    // fit
    for (int i = 1; i < 7; i++) {
        oldhist_FA[i]->Fit(gaus111[i], "NRQ");
    }
    // get the parameter
    double height[7];
    double mean_FA[7];
    double mean_SA[7];
    double sigma[7];
    double constant[7];

    for (int i = 1; i < 7; i++) {
        height[i] = gaus111[i]->GetParameter(0);
        mean_FA[i] = gaus111[i]->GetParameter(1);
        sigma[i] = gaus111[i]->GetParameter(2);
        constant[i] = gaus111[i]->GetParameter(3);

        printf("mean_FA%d = %f\n", i, mean_FA[i]);
    }
    /* fitting area SA ------------------------------------------------------*/
    TH1 *oldhist_SA[7];  // [0] is not used
    MINCH = -2000.;
    MAXCH = 2000.;
    // make the histgram
    for (int i = 1; i < 7; i++) {
        oldhist_FA[i] = new TH1F(
            Form("tdc_sa_SA %3d_sinti ch%d", RUNNUMBER, i),
            Form("#%3d double  gamma gaus fitting qdc ch%d;tdcsa[ns];Count",
                 RUNNUMBER, i),
            1000, MINCH, MAXCH);
    }
    // fill histgram
    for (int i = 1; i < 7; i++) {
        oldtree->Draw(
            Form("tdcsa_SA[%d]>>tdc_sa_SA %3d_sinti ch%d", i, RUNNUMBER, i),
            Form("IsNeutron[%d]==0 && tdcsa_SA[%d] > 0", i, i));
    }
    // set fitting parameter
    // 変更があるならこの場所
    // ---------------------------------------------------------------------------------------//
    double init_mu_SA[7];
    init_mu_SA[0] = 0;
    for (int ch = 1; ch < 7; ch++) {
        init_mu_SA[ch] = 90.;
    }
    // init_mu_SA[5] = 60;
    //  ---------------------------------------------------------------------------------------//
    //  make gaus+constant
    //  TF1 *gaus111[7];
    for (int i = 1; i < 7; i++) {
        gaus111[i] = new TF1(Form("gaus%3d_%d", RUNNUMBER, i),
                             "[0]*exp(-0.5*pow((x-[1])/[2],2.0))+[4]",
                             init_mu_SA[i] - range, init_mu_SA[i] + range);

        gaus111[i]->SetParameters(1000., init_mu_SA[i], 5., 100.);
        gaus111[i]->SetNpx(1000);
    }

    // fit
    for (int i = 1; i < 7; i++) {
        oldhist_FA[i]->Fit(gaus111[i], "NRQ");
    }

    for (int i = 1; i < 7; i++) {
        height[i] = gaus111[i]->GetParameter(0);
        mean_SA[i] = gaus111[i]->GetParameter(1);
        sigma[i] = gaus111[i]->GetParameter(2);
        constant[i] = gaus111[i]->GetParameter(3);

        printf("mean_SA%d = %f\n", i, mean_SA[i]);
    }

    /*--------------------------------------------- */
    int LOOP = oldtree->GetEntries();
    int evenum;
    for (int i = 0; i < LOOP; i++) {
        oldtree->GetEntry(i);
        evenum = i;

        for (int j = 0; j < 7; j++) {
            //
            IsNeutron[j] = IsNeutron[j];
            IsHit[j] = IsHit[j];
            if (j == 0 || IsHit[j] == 0) {
                continue;
            }
            // tdc_gamma-neutron [ns]
            tdc_gn_FA[j] = mean_FA[j] - tdcsa_FA[j];
            tdc_gn_SA[j] = mean_SA[j] - tdcsa_SA[j];
            Num = 4;
// original is 3
            ToFnraw_3_FA[j] = tof_g + tdc_gn_FA[j] + Num * trf;
            ToFnraw_3_SA[j] = tof_g + tdc_gn_SA[j] + Num * trf;
            // Tn3 N= 3 Energy[MeV]
            beta_FA[j] = Beta(ToFnraw_3_FA[j]);
            beta_SA[j] = Beta(ToFnraw_3_SA[j]);
            Tn_3_FA[j] = Ene(ToFnraw_3_FA[j]);
            Tn_3_SA[j] = Ene(ToFnraw_3_SA[j]);
            if (Tn_3_FA[j] < 0) {
                Tn_3_FA[j] = 0.;
            }
            if (Tn_3_SA[j] < 0) {
                Tn_3_SA[j] = 0.;
            }
            // Tn4
            /*
            Num = 4;
            ToFnraw_4[j] = tof_g + tdc_gn[j] + Num * trf;
            Tn_4[j] = Ene(ToFnraw_4[j]);
            if (Tn_4[j] < 0 && Tn_4[j] > 70)
            {
                Tn_4[j] = 0.;
            }
            */
        }
        tree->Fill();
        Etree->Fill();
    }
    tree->Write();
    Etree->Write();
}
double Ene(double time) {
    double beta = L / (time * c);
    double beta2 = pow(beta, 2.0);
    double Tn = Mn * (-1.0 + pow(1 - beta2, -0.5));
    return Tn;
}
double Beta(double time) {
    double beta = L / (time * c);
    return beta;
};
