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
bool input(double *r,  double *lz) {
    char filename[128];
    std::ifstream file;
    sprintf(filename,
            "lohis1t_lo0.500_tau50.000_m0.001.dat");
    file.open(filename);
    if(!file) {
        std::cout << "You Failed! " <<  takebit + takefst << std::endl;
        return true;
    }
    for(int i=0;i<80;i++){
        file >>r[i] >>lz[i];
    }

    file.close();
    std::ofstream file2;
    sprintf(filename,
            "lo1hist_lo0.500_tau50.000_m0.001.dat");
    file2.open(filename);
    for(int i=0;i<80;i++){
        file2<<r[i]<<" "<<lz[i]/817.<<std::endl;
    }
    file2.close();
    return false;
}

int main(){
    double r[80],ro[80];
    input(r,ro);
    return 0;
}
