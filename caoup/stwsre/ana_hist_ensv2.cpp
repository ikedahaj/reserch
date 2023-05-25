#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define v0 1.0
#define R 80.
#define lo 0.5
#define dim 2
#define Np 4 * R *R *lo
// #define taketimes 5 // とるファイルの数;
#define takebit 1 // 隣り合うファイルの間のファイル数を入力;
#define takefst 0  // とるファイルの最初;
#define bithist 1.
#define moji "tyouwaenn"  // 粒子位置のファイル名;
#define folder "stwr80" // 粒子位置のフォルダ名;
#define folder2 "stwens" // ダスフォルダ名;
bool input_test(int k0,double tau,double mgn) {
    char filename[128];
    std::ifstream file;
    sprintf(filename,
            "./%s_coorlo%.2fv0%.1ftau%.3fm%.3f/%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder, lo, v0, tau, mgn, moji, lo, tau, mgn, k0);
    file.open(filename);
    if(!file) {
        std::cout << "Till here " << k0 << std::endl;
        return false;
    }
    

    file.close();
    return true;
}
bool input(double *r, int k0, double *lz, int times,double tau,double mgn) {
    char filename[128];
    std::ifstream file;
    sprintf(filename,
            "./%s_coorlo%.2fv0%.1ftau%.3fm%.3f/%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder, lo, v0, tau, mgn, moji, lo, tau, mgn, k0);
    file.open(filename);
    if(!file) {
        std::cout << "You Failed! " << times * takebit + takefst << std::endl;
        return true;
    }
    double t, x[dim], v[dim], lzi;
    for(int i = 0; i < Np; i++) {
        file >> t >> x[0] >> x[1] >> v[0] >> v[1];
        lzi = x[0] * v[1] - x[1] * v[0];
        r[(int)(i + times * Np)] = sqrt(x[0] * x[0] + x[1] * x[1]);
        lz[(int)(i + times * Np)] = lzi / r[(int)(i + times * Np)];
    }

    file.close();
    return false;
}
void make_hist(double *r, double *lz, double *lohist, double *lzhist,
               double *omhist,int taketimes) {
    double bithistinv = 1. / bithist;
    int ens = Np * taketimes, Nphist = R / bithist + 1;
    // double bun=bithistinv/ens;
    for(int i = 0; i < ens; i++) {
        if(r[i] < R) {
            lohist[(int)floor(r[i] * bithistinv)] += 1.;
            lzhist[(int)floor(r[i] * bithistinv)] += lz[i];
            omhist[(int)floor(r[i] * bithistinv)] += lz[i] / r[i];
        }
    }
    for(int j = 0; j < Nphist; j++) {
        if(lohist[j] != 0) {
            lzhist[j] = lzhist[j] / lohist[j];
            omhist[j] /= lohist[j];
        } else {
            lzhist[j] = 0.;
            omhist[j] = 0.;
        }
    }
}
void make_difhist(double *r, double *lz, double *lohist, double *lzhist,
                  double *difflz, double *omhist, double *diffom,int taketimes) {
    int ens = Np * taketimes, loc;
    double bithistinv = 1. / bithist, dif;
    for(int i = 0; i < ens; i++) {
        if(r[i] < R) {
            loc = (int)floor(r[i] * bithistinv);
            dif = lz[i] - lzhist[loc];
            difflz[loc] += dif * dif / lohist[loc];
            dif = lz[i] / r[i] - omhist[loc];
            diffom[loc] += dif * dif / lohist[loc];
        }
    }
}
void output(double *lohist, double *lzhist, double *difflz, double *omhist,
            double *diffom,double tau,double mgn,int taketimes) {
    int Nphist = R / bithist + 1;
    char filename[128];
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/vtdiff_v0%.1f_tau%.3f_m%.3f.dat",
            folder2, lo, v0, tau, mgn, v0, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for(int i = 0; i < Nphist; i++)
        file << (i + 0.5) * bithist << "\t" << lzhist[i] << "\t"
             << sqrt(difflz[i] / lohist[i]) << std::endl;
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/omdiff_v0%.1f_tau%.3f_m%.3f.dat",
            folder2, lo, v0, tau, mgn, v0, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for(int i = 0; i < Nphist; i++)
        file << (i + 0.5) * bithist << "\t" << omhist[i] << "\t"
             << sqrt(diffom[i] / lohist[i]) << std::endl;
    file.close();

    sprintf(
        filename,
        "./%slo%.2fv0%.1ftau%.3fm%.3f/lohist_lo%.3f_v0%.1f_tau%.3f_m%.3f.dat",
        folder2, lo, v0, tau, mgn, lo, v0, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for(int i = 0; i < Nphist; i++)
        file << (i + 0.5) * bithist << "\t"
             << lohist[i] / (8. * bithist * bithist * (i + 0.5)*Np*taketimes) << std::endl;
    file.close();
    sprintf(
        filename,
        "./%slo%.2fv0%.1ftau%.3fm%.3f/v_thetahist_lo%.3f_tau%.3f_m%.3f.dat",
        folder2, lo, v0, tau, mgn, lo,  tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for(int i = 0; i < Nphist; i++)
        file << (i + 0.5) * bithist << "\t" << lzhist[i] << "\t" << std::endl;
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/"
            "omegahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder2, lo, v0, tau, mgn, lo,  tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for(int i = 0; i < Nphist; i++)
        file << (i + 0.5) * bithist << "\t" << omhist[i] << "\t" << std::endl;
    file.close();
}
void ini_hist(double *hist) {
    int Nphist = R / bithist + 1;
    for(int i = 0; i < Nphist; i++)
        hist[i] = 0.;
}
int main() {
    double tau,mgn;
    int taketimes=0,k0=takefst;
    std::cout<<"tau and mgn is "<<std::endl;
    std::cin>>tau>>mgn;
    int Nphist = R / bithist + 1;
    while(input_test(k0,tau,mgn)){
        taketimes++;
        k0+=takebit;
    }
    std::cout<<taketimes<<std::endl;
    double r[(int)(Np * taketimes)], lz[(int)(Np * taketimes)], lohist[Nphist],
        lzhist[Nphist], omhist[Nphist], diffom[Nphist], difflz[Nphist];
     k0 = takefst;
    ini_hist(lohist);
    ini_hist(lzhist);
    ini_hist(difflz);
    ini_hist(omhist);
    ini_hist(diffom);
    char foldername[128];
    sprintf(foldername, "%slo%.2fv0%.1ftau%.3fm%.3f", folder2, lo, v0, tau,
            mgn);
    const char *fname = foldername;
    mkdir(fname, 0777);
    for(int i = 0; i < taketimes; i++) {
        if(input(r, k0, lz, i,tau,mgn)) {
            return -1;
        }

        k0 += takebit;
    }
    make_hist(r, lz, lohist, lzhist, omhist,taketimes);
    make_difhist(r, lz, lohist, lzhist, difflz, omhist, diffom,taketimes);
    output(lohist, lzhist, difflz, omhist, diffom,tau,mgn,taketimes);
    // std::cout<<"sucseed!"<<std::endl;
    return 0;
}
