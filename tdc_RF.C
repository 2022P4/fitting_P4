/*
        qdcを用いたn/gammaの識別
        rootfileを立ち上げる
*/
#include <stdlib.h>

#include <iostream>
#define N 256

int odd_even(int number);

int tdc_RF(void)
{
    gROOT->Reset();
    gStyle->SetOptStat(1001110);
    // 　本番用
    // TFile *f =
    //     TFile::Open(Form("/mnt/p4root_2022/%04d_new.root", RUNNUMBER), "R");
    // TFile *g = TFile::Open(Form("./%04d_new2.root", RUNNUMBER), "recreate");
    // test
    int RUNNUMBER = 150;
    TFile *f = TFile::Open("./sample.root", "R");
    TFile *g = TFile::Open("./sample_new.root", "recreate");

    /*	wada 用
    TFile *f = TFile::Open("./../../root/0070.root","R");
    TFile *g = TFile::Open("./0070_new.root","recreate");
    */
    double qdc[16];
    double hiritsu[16];
    // 0 + 2*6 ch
    int qdcraw[16]; // qdcrawのデータの箱
    // add 6ch
    int evenum;        // イベントNo.の管理
    int IsNeutron[16]; // n/gammaのタグ付け、n=1,gamma=0
    int ch[16], ch_s;
    int sCh_qdc[16], sCh_q;
    int sCh_tdc[7], sCh_t;

    double factor[16];               //[7]で十分
    double pedestal[16], pedestal_s; // ぺデスタル
    double base[16], base_s;         // n/gammaの識別の基準

    int tdcraw[16];
    // int sCh[7], sCh_s;
    int tdc[16];
    Double_t tdcsec[7];
    Double_t tdcsa[7];

    int tdc_ch;
    double slope[7], slope_s, intercept[7], inter_s;

    // 読み込むプログラム-------------------------------------------------
    FILE *fp;
    char fname1[] = "./input/qdc_pedestal.dat";
    fp = fopen(fname1, "r");
    if (fp == NULL)
        printf("%s file not open!\n", fname1);
    char header1[N];
    int qdc_ch;
    fgets(header1, N, fp); // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &qdc_ch, &ch_s, &pedestal_s, &base_s) !=
           EOF)
    {
        // 詰め込み
        ch[qdc_ch] = ch_s;
        pedestal[qdc_ch] = pedestal_s;
        base[qdc_ch] = base_s;
        // 確認
        //  printf("%d %d %f
        //  %f\n",qdc_ch,ch[qdc_ch],pedestal[qdc_ch],base[qdc_ch]);
    }
    fclose(fp);

    char fname2[] = "./input/tdc_fit_para.dat";
    fp = fopen(fname2, "r");
    if (fp == NULL)
        printf("%s file not open!\n", fname2);
    char header2[N];
    fgets(header2, N, fp); // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &tdc_ch, &sCh_t, &slope_s, &inter_s) !=
           EOF)
    {
        // 詰め込み
        //  ch1-6 is  sCh 1-6, RF is sCh0
        slope[sCh_t] = slope_s;
        intercept[sCh_t] = inter_s;
        // 確認
        printf("%d %f %f\n", tdc_ch, slope[sCh_t], intercept[sCh_t]);
    }
    fclose(fp);

    /*----------------------------------------------------------------*/

    TTree *oldtree = (TTree *)f->Get("tree");
    // deifine oldtree
    oldtree->SetBranchAddress("qdcraw", qdcraw);
    oldtree->SetBranchAddress("tdcraw", tdcraw);
    // set Br old info
    TTree *tree = new TTree("tree", "tree");
    /* Branchの設定 */
    // qdc area
    tree->Branch("qdc", qdc, "qdc[16]/D");
    tree->Branch("sCh_qdc", sCh_qdc, "sCh_qdc[16]/I");
    tree->Branch("hiritsu", hiritsu, "hiritsu[16]/D");
    tree->Branch("IsNeutron", IsNeutron, "IsNeutron[16]/I");
    // tdc area
    tree->Branch("sCh_tdc", sCh_tdc, "sCh_tdc[7]/I");
    tree->Branch("tdc", tdc, "tdc[7]/I");
    tree->Branch("tdcsec", tdcsec, "tdcsec[7]/D");
    tree->Branch("tdcsa", tdcsa, "tdcsa[7]/D");

    int LOOP = oldtree->GetEntries();
    int e_ch, o_ch;
    /*	データの詰め込み　*/
    for (int i = 0; i < LOOP; i++)
    {
        oldtree->GetEntry(i);
        evenum = i;
        for (ch_s = 1; ch_s < 7; ch_s++)
        {
            e_ch = 2 * ch_s;
            o_ch = 2 * ch_s - 1;
            // sinti_ch to qdc_Ch
            switch (ch_s)
            {
            case 1:
                factor[ch_s] = ((double)qdcraw[2] - pedestal[2]) /
                               ((double)qdcraw[1] - pedestal[1]);
                break;
            case 2:
                factor[ch_s] = ((double)qdcraw[4] - pedestal[4]) /
                               ((double)qdcraw[3] - pedestal[3]);
                break;
            case 3:
                factor[ch_s] = ((double)qdcraw[8] - pedestal[8]) /
                               ((double)qdcraw[5] - pedestal[5]);
                break;
            case 4:
                factor[ch_s] = ((double)qdcraw[10] - pedestal[10]) /
                               ((double)qdcraw[9] - pedestal[9]);
                break;
            case 5:
                factor[ch_s] = ((double)qdcraw[12] - pedestal[12]) /
                               ((double)qdcraw[11] - pedestal[11]);
                break;
            case 6:
                factor[ch_s] = ((double)qdcraw[14] - pedestal[14]) /
                               ((double)qdcraw[13] - pedestal[13]);
                break;
            case 0:
                factor[0] = 0;
                break;
            }
            // rewrite factor
            //
            // factor[ch_s] = ((double)qdcraw[e_ch] - pedestal[e_ch]) /
            //               ((double)qdcraw[o_ch] - pedestal[o_ch]);
            // fast/allのパラメータ
        }

        for (int j = 0; j < 16; j++)
        {
            qdc[j] = qdcraw[j] - pedestal[j];
            // pedestalの分を差し引いたqdcの値
            if (j < 7)
            {
                int j_odd = 2 * j - 1;
                switch (j)
                {
                case 4:
                    j_odd = 9;
                    break;
                case 5:
                    j_odd = 11;
                    break;
                case 6:
                    j_odd = 13;
                    break;
                default:
                    j_odd = 2 * j - 1;
                    break;
                }
                sCh_qdc[j] = j;
                hiritsu[j] = factor[j];
                if (hiritsu[j] > base[j] &&
                    qdcraw[j_odd] - pedestal[j_odd] >
                        500)
                { // n/gamma && minimum energy
                    IsNeutron[j] = j;
                }
                else
                {
                    IsNeutron[j] = 0;
                }
            }
            else
            {
                // for j>7
                sCh_qdc[j] = 100;
                hiritsu[j] = 0;
                IsNeutron[j] = 7;
            }
            /*
                // sintiの番号をtdcに対応させる
                if (factor[odd_even(j)] > base[odd_even(j)] &&
                    qdc[j] > 500) {  // factor is variable
                    // n/gammaのタグ付け
                    IsNeutron[j] = odd_even(j);
                } else {
                    IsNeutron[j] = 0;
                }
            */
        }
        // for tdc
        for (int j = 0; j < 7; j++)
        {
            // sCh_t = (j + 1) % 7;
            sCh_t = j;
            // if (i == 1) cout << sCh_t << endl;
            sCh_tdc[sCh_t] = sCh_t;
            tdc[sCh_t] = tdcraw[j];
            // sintiの番号を対応させる
            // tdcsec[sCh_t] = slope[sCh_t] * (double)tdcraw[j] + intercept[sCh_t];
            tdcsec[sCh_t] = ((double)tdcraw[j] + intercept[sCh_t]) / slope[sCh_t];
            // y = ax+b -> y = (x+b)
            // ch->time[ns]
            tdcsa[sCh_t] = tdcsec[0] - tdcsec[sCh_t];
            // RF(ch0)-tdcsec [ns]
        }
        tree->Fill();
    }
    tree->Write();

    /*	test */
    // tree->Draw("qdc[1]:qdc[0]>>hadc0(1000,0,4000,1000,0,2000)","IsNeutron==1");
    printf("neutron count is %lld\n", tree->GetEntries("IsNeutron==1"));
    printf("all count is %lld\n", tree->GetEntries(""));

    return 0124;
}
//	define the function
int odd_even(int number)
{
    if (number % 2 == 0)
    { // 偶数なら/2したもの
        return number / 2;
    }
    else
    { // 偶数なら(+1)/2したもの
        return (1 + number) / 2;
    }
}
