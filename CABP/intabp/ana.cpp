#include <fstream>
// #include <iomanip>
#include <chrono>
#include <iostream>
#include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <sys/stat.h>

#define v0 1.0
#define R  20.
#define lo 0.5
#define Rbit 0.5 //-D
#define tau 40.
#define mgn 0.
#define dim 2
// #define Np      (4 * R *R *lo)
#define takebit 1 // 隣り合うファイルの間のファイル数を入力;
#define takefst 0 // とるファイルの最初;
#define bithist 1.
#define moji    "tyouwaenn" // 粒子位置のファイル名;
#define folder  "stw"       // 粒子位置のフォルダ名;
#define folder2 "stw"       // ダスフォルダ名;
constexpr double
usr_arccos(double theta) { // 1付近の誤差0.1程度.シミュレーションでは使うな;
    constexpr double waru[5] = {63. / 2816, 35. / 1152, 5. / 112, 3. / 40,
                                1. / 6};
    double           x = theta * theta;
    double asi = (((((waru[0] * x + waru[1]) * x + waru[2]) * x + waru[3]) * x +
                   waru[4]) *
                      x +
                  1.) *
                 theta;
    return M_PI_2 - asi;
}
constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 1000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}

static constexpr int    Nphist = R / bithist + 1;
static constexpr double center_left = -Rbit * 0.5 * R;
static constexpr double center_rignt = Rbit * 0.5 * R;
static constexpr double rbit_2 = Rbit * 0.5;
static constexpr double Npd =
    (lo * M_2_PI * 2. * R * R *
     (M_PI - usr_arccos(rbit_2) + rbit_2 * usr_sqrt(1 - rbit_2 * rbit_2))) *
    2;
static constexpr int Np = Npd;
inline double        dist2right(double *x) {
    double xb = x[0] - center_rignt;
    return xb * xb + x[1] * x[1];
}
inline double dist2left(double *x) {
    double xb = x[0] - center_left;
    return xb * xb + x[1] * x[1];
}
bool input_test(int k0) {
    char          filename[128];
    std::ifstream file;
    sprintf(filename,
            "./%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
            "%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder, R, lo, tau, mgn, Rbit, v0, moji, lo, tau, mgn, k0);
    file.open(filename);
    if (!file) {
        std::cout << "Till here " << k0 << std::endl;
        return false;
    }

    file.close();
    return true;
}
bool input(double *r, int k0, double *lz, int times, double *t, double *v,
           int (*rorl)[2]) {
    char          filename[128];
    std::ifstream file;
    sprintf(filename,
            "./%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
            "%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder, R, lo, tau, mgn, Rbit, v0, moji, lo, tau, mgn, k0);
    file.open(filename);
    if (!file) {
        std::cout << "You Failed! " << times * takebit + takefst << std::endl;
        return true;
    }
    double ti, x[Np][dim], vl[Np][dim], lzi;
    int    mid;
    for (int i = 0; i < Np; i++)
        file >> ti >> x[i][0] >> x[i][1] >> vl[i][0] >> vl[i][1];
    for (int i = 0; i < Np; i++) {
        if (x[i][0] > 0.) {
            x[i][0] -= center_rignt;
            lzi = x[i][0] * vl[i][1] - x[i][1] * vl[i][0];
            mid = (int) (i + times * Np);
            rorl[mid][0] = 1;
            rorl[mid][1] = (x[i][1] >= 0.) - (x[i][1] < 0.);
            r[mid] = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
            lz[mid] = lzi / (r[mid] * r[mid]);
            v[mid] = sqrt(vl[i][0] * vl[i][0] + vl[i][1] * vl[i][1]);
        } else if (x[i][0] < 0.) {
            x[i][0] -= center_left;
            lzi = x[i][0] * vl[i][1] - x[i][1] * vl[i][0];
            mid = (int) (i + times * Np);
            rorl[mid][0] = -1;
            rorl[mid][1] = (x[i][1] >= 0.) - (x[i][1] < 0.);
            r[mid] = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
            lz[mid] = lzi / (r[mid] * r[mid]);
            v[mid] = sqrt(vl[i][0] * vl[i][0] + vl[i][1] * vl[i][1]);
        } else {
            mid = (int) (i + times * Np);
            rorl[mid][0] = 0;
            rorl[mid][1] = (x[i][1] >= 0.) - (x[i][1] < 0.);
            r[mid] = abs(x[i][1]);
            lz[mid] = vl[i][0] / (r[mid] * r[mid]);
            v[mid] = sqrt(vl[i][0] * vl[i][0] + vl[i][1] * vl[i][1]);
        }
    }
    t[times] = ti;

    file.close();
    return false;
}