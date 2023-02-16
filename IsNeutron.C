/*
        qdcを用いたn/gammaの識別
        rootfileを立ち上げる
*/
#include <stdlib.h>

#include <iostream>
#define N 256

int odd_even(int number);

int IsNeutron(int RUNNUMBER) {
    gROOT->Reset();
    gStyle->SetOptStat(1001110);
    //　本番用
    // int RUNNUMBER = 150;
    TFile *f = TFile::Open(Form("/mnt/p4root_2022/%04d.root", RUNNUMBER), "R");
    TFile *g = TFile::Open(Form("./%04d_new.root", RUNNUMBER), "recreate");

    /*	wada 用
    TFile *f = TFile::Open("./../../root/0070.root","R");
    TFile *g = TFile::Open("./0070_new.root","recreate");
    */
    double qdc[13];
    // 0 + 2*6 ch
    int qdcraw[16];
    // add 6ch
    int evenum;         //イベントNo.の管理
    int IsNeutron[13];  // n/gammaのタグ付け、n=1,gamma=0
    int ch[13], ch_s;
    int sintiCh[13];
    double factor[13];                //[7]で十分
    double pedestal[13], pedestal_s;  //ぺデスタル
    double base[13], base_s;          // n/gammaの識別の基準
    //読み込むプログラム-------------------------------------------------
    FILE *fp;
    char fname[] = "./input/qdc_pedestal.dat";
    fp = fopen(fname, "r");
    if (fp == NULL) printf("%s file not open!\n", fname);
    char header[N];
    int qdc_ch;
    fgets(header, N, fp);  // headerの読み込み
    while (fscanf(fp, "%d%d%lf%lf", &qdc_ch, &ch_s, &pedestal_s, &base_s) !=
           EOF) {
        //詰め込み
        ch[qdc_ch] = ch_s;
        pedestal[qdc_ch] = pedestal_s;
        base[qdc_ch] = base_s;
        //確認
        // printf("%d %d %f
        // %f\n",qdc_ch,ch[qdc_ch],pedestal[qdc_ch],base[qdc_ch]);
    }
    fclose(fp);
    /*----------------------------------------------------------------*/

    TTree *oldtree = (TTree *)f->Get("tree");
    // deifine oldtree
    oldtree->SetBranchAddress("qdcraw", qdcraw);
    // set Br old info
    TTree *tree = new TTree("tree", "tree");
    /* Branchの設定 */
    tree->Branch("sintiCh", sintiCh, "sintiCh[13]/I");
    tree->Branch("qdc", qdc, "qdc[13]/D");
    tree->Branch("IsNeutron", IsNeutron, "IsNeutron[13]/I");

    int LOOP = oldtree->GetEntries();
    int e_ch, o_ch;
    /*	データの詰め込み　*/
    for (int i = 0; i < LOOP; i++) {
        oldtree->GetEntry(i);
        evenum = i;
        for (ch_s = 1; ch_s < 7; ch_s++) {
            e_ch = 2 * ch_s;
            o_ch = 2 * ch_s - 1;
            factor[ch_s] = ((double)qdcraw[e_ch] - pedestal[e_ch]) /
                           ((double)qdcraw[o_ch] - pedestal[o_ch]);
            // fast/allのパラメータ
        }
        for (int j = 0; j < 13; j++) {
            qdc[j] = qdcraw[j] - pedestal[j];
            // pedestalの分を差し引いたqdcの値
            sintiCh[j] = odd_even(j);
            // sintiの番号を対応させる
            if (factor[odd_even(j)] > base[odd_even(j)] &&
                qdc[j] > 500) {  // factor is variable
                // n/gammaのタグ付け
                IsNeutron[j] = odd_even(j);
            } else {
                IsNeutron[j] = 0;
            }
        }
        tree->Fill();
    }
    tree->Write();

    /*	test	*/
    // tree->Draw("qdc[1]:qdc[0]>>hadc0(1000,0,4000,1000,0,2000)","IsNeutron==1");
    printf("neutron count is %lld\n", tree->GetEntries("IsNeutron==1"));
    printf("all count is %lld\n", tree->GetEntries(""));

    return 0124;
}
//	define the function
int odd_even(int number) {
    if (number % 2 == 0) {  //偶数なら/2したもの
        return number / 2;
    } else {  //偶数なら(+1)/2したもの
        return (1 + number) / 2;
    }
}
