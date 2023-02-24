/*
        qdcを用いたn/gammaの識別
        rootfileを立ち上げる
*/
#include <math.h>
#include <stdlib.h>

#include <iostream>
#define N 256
#define Mn 939.56
#define L 33.0        //[m]
#define c 0.299792458 //[m/ns]

double Ene(double Ene);

void tdc_fit(void)
{
    int RUNNUMBER = 115;
    gROOT->Reset();
    gStyle->SetOptStat(1001110);
    // 　本番用
    //  int RUNNUMBER = 150;

    // TFile *f =
    //     TFile::Open(Form("/mnt/p4root_2022/%04d_new.root", RUNNUMBER), "R");
    // TFile *g = TFile::Open(Form("./%04d_new2.root", RUNNUMBER), "recreate");
    // test
    TFile *f = TFile::Open("./sample_new.root", "R");
    TFile *g = TFile::Open("./sample_new2.root", "recreate");
    /*	wada 用
    TFile *f = TFile::Open("./../../root/0070.root","R");
    TFile *g = TFile::Open("./0070_new.root","recreate");
    */

    int IsNeutron[16];
    double tdcsa[7];
    // new area
    double tof_g = 110.076; //[ns]
    // tdc_gamma - tdc_n
    double tdc_gn[7];
    double trf = 60; //[ns]
    double ToFnraw_3[7];
    double Tn_3[7];
    // double ToFnraw_4[7];
    // double Tn_4[7];
    int Num = 0;
    // reset the data
    for (int i = 1; i < 7; i++)
    {
        tdcsa[i] = 0.;
        tdc_gn[i] = 0.;
    }
    TTree *oldtree = (TTree *)f->Get("tree");

    oldtree->SetBranchAddress("IsNeutron", IsNeutron);
    oldtree->SetBranchAddress("tdcsa", tdcsa);

    TTree *tree = new TTree("tree", "tree");
    // old area =================================================//
    tree->Branch("IsNeutron", IsNeutron, "IsNeutron[13]/I");
    tree->Branch("tdcsa", tdcsa, "tdcsa[7]/D");
    // new area =================================================//
    tree->Branch("tdc_gn", tdc_gn, "tdc_gn[7]/D");
    tree->Branch("ToFnraw_3", ToFnraw_3, "ToFnraw_3[7]/D");
    tree->Branch("Tn_3", Tn_3, "Tn_3[7]/D");

    /* fitting area ---------------------------------------------------------*/
    TH1 *oldhist[7]; // [0] is not used
    double MINCH = -2000.;
    double MAXCH = 2000.;
    // make the histgram
    for (int i = 1; i < 7; i++)
    {
        oldhist[i] = new TH1F(
            Form("tdc_sa %3d_sinti ch%d", RUNNUMBER, i),
            Form("#%3d double  gamma gaus fitting qdc ch%d;tdcsa[ns];Count", RUNNUMBER, i),
            1000, MINCH, MAXCH);
    }
    // fill histgram
    for (int i = 1; i < 7; i++)
    {
        oldtree->Draw(Form("tdcsa[%d]>>tdc_sa %3d_sinti ch%d", i, RUNNUMBER, i),
                      Form("IsNeutron[%d]==0", i));
    }

    // set fitting parameter
    double init_mu[7];
    init_mu[0] = 0;
    for (int ch = 1; ch < 7; ch++)
    {
        init_mu[ch] = 110.;
    }
    double range = 100.;

    // make gaus+constant
    TF1 *gaus111[7];
    for (int i = 1; i < 7; i++)
    {
        gaus111[i] = new TF1(Form("gaus%3d_%d", RUNNUMBER, i),
                             "[0]*exp(-0.5*pow((x-[1])/[2],2.0))+[4]",
                             init_mu[i] - range, init_mu[i] + range);

        gaus111[i]->SetParameters(1000., init_mu[i], 50., 10.);
        gaus111[i]->SetNpx(1000);
    }

    // fit
    for (int i = 1; i < 7; i++)
    {
        oldhist[i]->Fit(gaus111[i], "NRQ");
    }
    // get the parameter
    double height[7];
    double mean[7];
    double sigma[7];
    double constant[7];

    for (int i = 1; i < 7; i++)
    {
        height[i] = gaus111[i]->GetParameter(0);
        mean[i] = gaus111[i]->GetParameter(1);
        sigma[i] = gaus111[i]->GetParameter(2);
        constant[i] = gaus111[i]->GetParameter(3);

        printf("mean%d = %f\n", i, mean[i]);
    }

    /* test area -------------------------------------
    int i = 1;
    oldhist[1] = new TH1F(
        Form("tdc_sa %3d_sinti ch%d", RUNNUMBER, i),
        Form("#%3d double  gamma gaus fitting qdc ch%d;tdcsa[ns];Count", RUNNUMBER, i),
        1000, MINCH, MAXCH);

    oldtree->Draw(Form("tdcsa[1]>>tdc_sa %3d_sinti ch%d", RUNNUMBER, i),
                  Form("IsNeutron[%d]==0", i));

    double init_mu[7];
    init_mu[0] = 0;
    for (int ch = 1; ch < 7; ch++)
    {
        init_mu[ch] = 110.;
    }

    double range = 100.;

    TF1 *gaus111[7];
    gaus111[1] = new TF1(Form("gaus%3d_%d", RUNNUMBER, i),
                         "[0]*exp(-0.5*pow((x-[1])/[2],2.0))+[4]",
                         init_mu[i] - range, init_mu[i] + range);

    gaus111[1]->SetParameters(1000., init_mu[i], 50., 10.);
    gaus111[1]->SetNpx(1000);

    oldhist[1]->Fit(gaus111[i], "NRQ");

    double height[7];
    double mean[7];
    double sigma[7];
    double constant[7];

    height[i] = gaus111[i]->GetParameter(0);
    mean[i] = gaus111[i]->GetParameter(1);
    sigma[i] = gaus111[i]->GetParameter(2);
    constant[i] = gaus111[i]->GetParameter(3);
    */
    /*
    int test = 1;
    printf("#%04d_ch%d\n", RUNNUMBER, test);
    cout << "height = " << height[test] << endl;
    cout << "mean = " << mean[test] << endl;
    cout << "sigma = " << sigma[test] << endl;
    cout << "const = " << constant[test] << endl;
    */
    /*--------------------------------------------- */
    int LOOP = oldtree->GetEntries();
    int evenum;
    for (int i = 0; i < LOOP; i++)
    {
        oldtree->GetEntry(i);
        evenum = i;

        for (int j = 0; j < 7; j++)
        {
            //
            IsNeutron[j] = IsNeutron[j];
            // tdc_gamma-neutron [ns]
            tdc_gn[j] = mean[j] - tdcsa[j];
            Num = 3;
            ToFnraw_3[j] = tof_g + tdc_gn[j] + Num * trf;
            // Tn3 N= 3 Energy[MeV]
            Tn_3[j] = Ene(ToFnraw_3[j]);
            if (Tn_3[j] < 0)
            {
                Tn_3[j] = 0.;
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
    }
    tree->Write();
}
double Ene(double time)
{
    double beta = L / (time * c);
    double beta2 = pow(beta, 2.0);
    double Tn = Mn * (-1.0 + pow(1 - beta2, -0.5));
    return Tn;
}
