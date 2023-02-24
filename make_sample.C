#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define PI 3.14159265358979323846
int randn(double mean, double sd);

int make_sample()
{
    double r;
    double u1, u2, z;
    int qdcraw[16];
    int tdcraw[16];
    int IsNeutron[16];
    for (int i = 0; i < 16; i++)
    {
        IsNeutron[i] = 0;
        qdcraw[i] = 0;
        tdcraw[i] = 0;
    }
    TFile *f = TFile::Open("./sample.root", "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("qdcraw", qdcraw, "qdcraw[16]/I");
    tree->Branch("tdcraw", tdcraw, "tdcraw[16]/I");
    tree->Branch("IsNeutron", IsNeutron, "IsNeutron[16]/I");

    // 乱数の初期化
    srand((unsigned int)time(NULL));

    for (int i = 0; i < pow(10, 6); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            int even = 2 * (j + 1);
            int odd = 2 * (j + 1) - 1;
            int num = (rand() % 100);
            if (num > 1)
            {
                // gamma ray
                IsNeutron[j] = 0;
                qdcraw[odd] = randn(3500., 300.) % 4095;
                qdcraw[even] = (int)round((double)qdcraw[odd] * 0.15);
                if (j == 0) // RF
                {
                    tdcraw[j] = randn(3584., 100.);
                }
                else
                {
                    tdcraw[j] = randn(1344, 300);
                }
                /* 2022.02.24 rewrite
                            if (j > 0 && (rand() % 10 + 1) > 4)
                            {
                                tdcraw[j] = rand() % 4096;
                            }
                            else if (j > 0 && (rand() % 10 + 1) <= 4)
                            {
                                tdcraw[j] = randn(2500.0, 150.0);
                            }
                            else if (j == 0)
                            {
                                tdcraw[j] = randn(3000., 100.0);
                            }
                */
            }
            else if (num <= 1)
            {
                // neutron
                IsNeutron[j + 1] = j + 1;
                qdcraw[odd] = randn(3000., 300.) % 4095;
                qdcraw[even] = (int)round((double)qdcraw[odd] * 0.24);
                if (j == 0) // RF
                {
                    tdcraw[j] = randn(2716., 400.);
                }
                else
                {
                    tdcraw[j] = randn(1344, 300);
                }
                /* 2023.02.24
                if (j < 6 && (rand() % 10 + 1) > 3)
                {
                    tdcraw[j] = rand() % 4096;
                }
                else if (j < 6 && (rand() % 10 + 1) <= 3)
                {
                    tdcraw[j] = randn(1000.0, 150.0);
                }
                else if (j == 6)
                {
                    tdcraw[j] = randn(2000.0, 70.0);
                }
                else
                {
                    tdcraw[j] = 0;
                }
                */
            }
            // qdcraw

            // qdcraw[j] = randn(1000, 150);
            // tdcraw[j] = rand() % 4096;
            //  = (double)rand() / ((double)RAND_MAX + 1) * 0.5;
        }
        // 0からRAND_MAXまでの乱数を生成し、0から0.5にスケーリングする
        tree->Fill();
    }
    tree->Write();
    return 0;
}

int randn(double mean, double sd)
{
    double u1, u2, z;
    int r;

    // 乱数の初期化
    // srand((unsigned int)time(NULL));

    // 0から1までの一様分布に従う乱数を2つ生成
    u1 = (double)rand() / RAND_MAX;
    u2 = (double)rand() / RAND_MAX;

    // Box-Muller法により、標準正規分布に従う乱数を生成
    z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

    // 平均と標準偏差を考慮して、整数に変換
    r = (int)round(z * sd + mean);

    return r;
}