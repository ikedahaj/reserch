
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <chrono>
#include <deque>
#include <fstream>
#include <iostream>
#include <vector>

#include "BM.h"
// #define Np          12800 // 4の倍数であること;NP=4*r^2*lo
#define Nn                    50
#define tmax                  100  // 100000  // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg                0  // 5000  // 緩和時間は10たうとする;
#define tmaxani               100  // 900  //>tmaxの時プログラムを変更すること;
#define tbitani               1
#define dim                   2  // 変えるときはEomを変えること;
#define tmax_anill            10
#define skin                  1.5
#define dtlg                  1e-4
#define dt                    1e-4
#define gay_kappa             3.
#define folder_name           "ani_test"  //"h0" // 40文字程度で大きすぎ;
#define UPDATE_MAX(x, x_past) (x_past = (x > x_past) ? x : x_past)
#define msdBit                3
#define msdtill               2
#define msdini                0.01
#define temp                  5.0  // langevinで用いÏる温度;
#define ratf                  1.0
#define ratf_w                1.
#define w_list                ver_list  // ver_list or cell_list
#define w_list_anill          ver_list_anill
#define FLAG_MASS             0  // 1なら慣性あり0なら慣性なし;
#define cut                   1.122462048  // 3.

using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.
#ifndef lo
#define lo 0.1  // コンパイル時に代入する定数;
#endif
#define Rbit 1.8  // delta/R,Rにすると穴がなくなる;//
// コンパイル時に-D{変数名}={値}　例:-Dbit=80　とすること;
#define v0 1.
#ifndef TAU
#define TAU 100
#endif
#ifndef MS
#if FLAG_MASS == 1
#define MS 0.1
#else
#define MS 0.0000001
#endif
#endif
#ifndef Rs
#define Rs 10
#endif
static constexpr double tau = TAU;
static constexpr double mass = MS;
static constexpr double mgn = 0.;
static constexpr double R = Rs;  // 13.07252733;  // 固定;// ,0.1より大きいこと;
////parameters
constexpr double usr_arccos(double theta) {
    // 1付近の誤差0.1程度.シミュレーションでは使うな;
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
static constexpr double ini_rho = (0.6 / (gay_kappa * gay_kappa));
static constexpr double rat_anil =
    (gay_kappa - 0.1) / gay_kappa;       // アニールするときの比;
static constexpr double para3_tbit = 1;  // double;
static constexpr double rbit_2 = Rbit * 0.5;
static constexpr double Npd =
    (lo * 2. * M_2_PI * R * R / gay_kappa *
     (M_PI - usr_arccos(rbit_2) + rbit_2 * usr_sqrt(1 - rbit_2 * rbit_2))) *
    2.;
static constexpr int Np = (int) (Npd);
static constexpr int Np_act = Np;
static constexpr int Np_lange = Np - Np_act;
// static constexpr double R=usr_sqrt(Np/((lo * 2. * M_2_PI  *
//					(M_PI - usr_arccos(rbit_2) + rbit_2 *
////usr_sqrt(1 - rbit_2 * rbit_2))) *
//					2.));
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double M_1_PI2 = 1. / M_PI2;
static constexpr double Mg = mgn * dt;
static constexpr double Np_1 = 1. / Np;
static constexpr double center_left = -Rbit * 0.5 * R;
static constexpr double center_rignt = Rbit * 0.5 * R;
static constexpr double x0limit = R * usr_sqrt(1 - Rbit * Rbit * 0.25);
static constexpr int    tcoorch = (tau > 10) ? tau / dt : 10 / dt;
static constexpr double const_f = -48. * ratf;

static constexpr double msdbit = msdBit / dt;
static constexpr double msdtillch = msdtill / dt;
// χとかの出す時間について::１万アンサンブルぐらい取るようにパラメータを設定;
// 理想の出す数;１で固定するようにした;
static constexpr int para3_bitn = 1;
// 実際の出す間隔 intで;
static constexpr double para3_bitlonch =
    (para3_tbit / para3_bitn < dt) ? 1 : (int) (para3_tbit / para3_bitn / dt);
// till here;
char foldername1[128];
char foldername_ani[128];
char foldername_coor[128];
#include "ana_ft.h"
#include "analyzepart_single_circle.h"
#include "gay_bane_force_two_circle.hpp"

void usr_sincos(double kaku, double *x) {  // x[0]がcos,x[1]issin;
                                           // 制度は10^-13程度;
    constexpr double waru[8] = {1.0 / (3 * 4 * 5 * 6 * 7 * 8 * 9 * 10),
                                -1.0 / (3 * 4 * 5 * 6 * 7 * 8),
                                1.0 / (3 * 4 * 5 * 6),
                                -1.0 / (3 * 4),
                                1.0 / (2 * 3 * 4 * 5 * 6 * 7 * 8 * 9),
                                -1.0 / (2 * 3 * 4 * 5 * 6 * 7),
                                1.0 / (2 * 3 * 4 * 5),
                                -1.0 / (2 * 3)};
    kaku *= 0.0625;  // 0.03125;//kaku/=1/2^m;
    double c, s, z = kaku * kaku;
    c = ((((waru[0] * z + waru[1]) * z + waru[2]) * z + waru[3]) * z + 1.) * z;
    s = ((((waru[4] * z + waru[5]) * z + waru[6]) * z + waru[7]) * z + 1.) *
        kaku;
    for (int i = 0; i < 4; i++) {  // mmade;
        s = s * (2.0 - c);
        c = c * (4.0 - c);
    }
    x[0] = 1.0 - c * 0.5;
    x[1] = s;
}
inline double dist2right(double *x) {
    double xb = x[0] - center_rignt;
    return xb * xb + x[1] * x[1];
}
inline double dist2left(double *x) {
    double xb = x[0] - center_left;
    return xb * xb + x[1] * x[1];
}
inline double dist2right_r(double *x, double R_ini) {
    double xb = x[0] - R_ini * ( rbit_2);
    return xb * xb + x[1] * x[1];
}
inline double dist2left_r(double *x, double R_ini) {
    double xb = x[0] + R_ini * (rbit_2);
    return xb * xb + x[1] * x[1];
}

/// @brief 粒子座標の初期値;
/// @param x
/// @param flag 最初に通る時はfalse
/// @return
bool ini_coord_twocircles(double (*x)[dim], bool flag) {
    double R_ini = sqrt(
               Np * gay_kappa /
               ((ini_rho * 2. * M_2_PI *
                 (M_PI - acos(rbit_2) + rbit_2 * sqrt(1 - rbit_2 * rbit_2))) *
                2.)),
           rat = 1 + rbit_2;
    int Np_img = 4 * ini_rho * 4 * R_ini * R_ini * rat / (M_PI * gay_kappa);
    int nx, ny, cnt = 0;
    if (Np_img > 4 * lo * R * R * rat / M_PI / gay_kappa)
        nx = rat * sqrt(Np_img / rat) , ny = sqrt(Np_img / rat) ;
    else {
        Np_img = 4 * lo * R * R * rat / M_PI / gay_kappa,
        nx = rat * sqrt(Np / rat) , ny = sqrt(Np_img / rat) ;
        R_ini = R;
    }
    if(nx*ny<Np_img){
        ++nx;
        ++ny;
    }
    if (flag) {
        nx++;
        ny++;
    }
    double poj[2], r2 = (R_ini - gay_kappa * 0.5) * (R_ini - gay_kappa * 0.5);
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            poj[0] = 0.5 + i * 2 * R_ini * rat / nx - R_ini * (1 + rbit_2);
            poj[1] = 0.5 + j * 2 * R_ini / ny - R_ini;
            if (dist2left_r(poj,R_ini) < r2 || dist2right_r(poj,R_ini) < r2) {
                x[cnt][0] = poj[0];
                x[cnt][1] = poj[1];
                cnt++;
            }
            if (cnt == Np) break;
        }
        if (cnt == Np) break;
    }
    if (cnt != Np) {
        std::cout << cnt << "\t" << Np << endl;
        if (!flag && ini_coord_twocircles(x, true))
            return true;
        else
            return false;
    } else {
        std::cout << cnt << "\t"
                  << "Yes" << endl;
        return true;
    }
}

void set_diameter(double *a) {
    for (int i = 0; i < Np; ++i) a[i] = 0.5;
}
void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < dim; ++j) x[i][j] = 0.0;
}

void eom_abp9(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i, double *f_theta) {
    static constexpr double ddt = 1e-7;
    static constexpr double D = usr_sqrt(2. * ddt / tau), M_inv = ddt / mass;
    double                  ri, riw, w2, w6, dUr;
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        calc_force_wall_gay_bane(x[i], f[i], theta_i[i], &f_theta[i]);
        theta_i[i] += D * gaussian_rand() + Mg + 3 * f_theta[i] * ddt;
        theta_i[i] -= floor(theta_i[i] * M_1_PI2) * M_PI2;
// usr_sincos(theta_i[i], sico);
#if FLAG_MASS == 1
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1]) * M_inv;
#else
        v[i][0] = (v0 * cos(theta_i[i]) + f[i][0]);
        v[i][1] = (v0 * sin(theta_i[i]) + f[i][1]);
#endif
        x[i][0] += v[i][0] * ddt;
        x[i][1] += v[i][1] * ddt;
    }
}

void eom_langevin(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                  int (*list)[Nn], double *theta_i, double *f_theta,
                  double R_ini) {
    double zeta = 1.0, ddt = 1e-9;
    double fluc = sqrt(2. * zeta * 5. * ddt);
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        // /*force bitween wall;
        calc_wall_gay_bane_anill(x[i], f[i], theta_i[i], &f_theta[i], R_ini);
        f_theta[i] += f_theta[i] * ddt * 3;
        theta_i[i] -= (floor(theta_i[i] * M_1_PI2)) * M_PI2;
        for (int j = 0; j < dim; j++) {
            v[i][j] +=
                -zeta * v[i][j] * ddt + f[i][j] * ddt + fluc * gaussian_rand();
            x[i][j] += v[i][j] * ddt;
        }
    }
}
void eom_langevin_h(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                    int (*list)[Nn], double *theta, double *f_theta,
                    double temp0, double R_ini) {
    double        zeta = 1;
    static double fluc = sqrt(2. * zeta * temp0 * dtlg);
    double        ri, riw, w2, w6, dUr;
    calc_force_gay_bane(x, f, list, theta, f_theta);
    for (int i = 0; i < Np; i++) {
        // /*force bitween wall;
        calc_wall_gay_bane_anill(x[i], f[i], theta[i], &f_theta[i], R_ini);
        theta[i] += f_theta[i] * dtlg * 3;
        theta[i] -= (floor(theta[i] * M_1_PI2)) * M_PI2;
        for (int j = 0; j < dim; j++) {
            v[i][j] +=
                -v[i][j] * dtlg + f[i][j] * dtlg + fluc * gaussian_rand();
            x[i][j] += v[i][j] * dtlg;
        }
    }
}
void eom_8(double (*v)[dim], double (*x)[dim], double (*f)[dim],
           int (*list)[Nn], double *theta_i, double *f_theta) {
    // double                  sico[2];
    static constexpr double D = usr_sqrt(2. * dt / tau), M_inv = dt / mass;
    double                  ri, riw, w2, w6, dUr;
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        // /*force bitween wall;
        calc_force_wall_gay_bane(x[i], f[i], theta_i[i], &f_theta[i]);
        theta_i[i] += D * gaussian_rand() + Mg + 3. * f_theta[i] * dt;
        theta_i[i] -= floor(theta_i[i] * M_1_PI2) * M_PI2;
// usr_sincos(theta_i[i], sico);
#if FLAG_MASS == 1
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1]) * M_inv;
#else
        v[i][0] = (v0 * cos(theta_i[i]) + f[i][0]);
        v[i][1] = (v0 * sin(theta_i[i]) + f[i][1]);
#endif
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}
void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i, double *f_theta) {
    double                  ri, riw, w2, w6, dUr;
    static constexpr double D = usr_sqrt(2. * dt / tau), M_inv = dt / mass;
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        // /*force bitween wall;
        calc_force_wall_gay_bane(x[i], f[i], theta_i[i], &f_theta[i]);
        theta_i[i] += D * gaussian_rand() + Mg + f_theta[i] * dt;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
// usr_sincos(theta_i[i], sico);
#if FLAG_MASS == 1
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1]) * M_inv;
#else
        v[i][0] = (v0 * cos(theta_i[i]) + f[i][0]);
        v[i][1] = (v0 * sin(theta_i[i]) + f[i][1]);
#endif
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}
inline double usr_abs(double x) { return x * ((x > 0) - (x < 0)); }

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = unif_rand(0, 2) * M_PI;
    }
}
void output_ini(double (*v)[dim], double (*x)[dim]) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "./%s/tyokkei.dat", foldername_ani);
    file.open(filename /* std::ios::app*/);  // append
    file << tbitani << endl;
    for (int i = 0; i < Np; i++) file << 1 << endl;
    file.close();
    snprintf(filename, 128, "./%s/coor_R%.3f_m%.3f_t0.dat", foldername1, R,
             mgn);
    file.open(filename);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
}
void output(double (*v)[dim], double (*x)[dim]) {
    static int l = 0;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128, "./%s/coor_R%.3f_m%.3f_t%d.dat", foldername1, R,
             mgn, l);
    file.open(filename);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
void output_ani(double (*v)[dim], double (*x)[dim], double *theta_i) {
    static int l = 0;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128, "./%s/coor_m%.3f_t%d.dat", foldername_ani, mgn, l);
    file.open(filename /* std::ios::app*/);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << "\t" << theta_i[i] << endl;
    }
    file.close();
    l++;
}

bool out_setup() {  // filenameが１２８文字を超えていたらfalseを返す;
    char     filename[128];
    ofstream file;
    int      test =
        snprintf(filename, 128, "./%s/setupofst_Rs%.3f_tau%.3f_m%.3f_t%d.dat",
                 foldername1, lo, tau, mgn, tmax);
    std::cout << test << endl;
    file.open(filename, std::ios::app);  // append

    file << "dt=" << dt << endl;
    file << "cut" << cut << endl;
    file << "skin" << skin << endl;
    file << "Nn" << Nn << endl;
    file << "Np=" << Np << endl;
    file << "tamxlg=" << tmaxlg << endl;
    file << "tmax=" << tmax << endl;
    file << "v0=" << v0 << endl;
    file << "2DkaraD" << endl;
    file << "壁はWCA" << endl;
    file << "cell list" << endl;
    file << "usr_sincos" << endl;
    file << "自動Np" << endl;
    file << "x=0での壁を追加" << endl;
    file << "modNp" << endl;
    file << "χの平均を取る間隔" << para3_bitlonch << endl;
    file << "kaiの出す感覚" << para3_tbit << endl;
    file << "ratf=" << ratf << endl;
    file << "lo2="
         << Np / ((2. * M_2_PI * R * R *
                   (M_PI - usr_arccos(rbit_2) +
                    rbit_2 * usr_sqrt(1 - rbit_2 * rbit_2))) *
                  2.)
         << endl;
    file << "eta=1,t=1" << endl;
    file << "temp" << temp << endl;
    file.close();
    if (test == -1)
        return false;
    else
        return true;
}

void update(double (*x_update)[dim], double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++) x_update[i][j] = x[i][j];
}

void calc_disp_max(double *disp_max, double (*x)[dim],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp;
    for (int i = 0; i < Np; i++) {
        dx = x[i][0] - x_update[i][0];
        dy = x[i][1] - x_update[i][1];
        disp = dx * dx + dy * dy;
        if (disp > *disp_max) *disp_max = disp;
    }
}

void auto_list_update(double (*x)[dim], double (*x_update)[dim],
                      int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    static double           disp_max = skin2 + 100;
    calc_disp_max(&(disp_max), x, x_update);
    if (disp_max >= skin2) {
        w_list(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        disp_max = skinini;
        // count = 0;
    }
}

#include "two_wall_anill_elli.hpp"

int main() {
    std::chrono::system_clock::time_point start, end;  // 型は auto で可
    start = std::chrono::system_clock::now();          // 計測開始時間
    double x[Np][dim], v[Np][dim], theta[Np], f[Np][dim], f_theta[Np],
        x0[Np][dim], v1[Np][dim], x_update[Np][dim];
    // int(*list)[Nn] = new int[Np][Nn];
    int    list[Np][Nn];
    int    countout = 0;
    double tout = msdini;
    // set_diameter(a);
    if (!ini_coord_twocircles(x, false)) return -1;
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    ini_hist(theta, Np);
    snprintf(foldername1, 128, "%slo%.3fMs%.3ftau%.3fkap%.3fbit%.3fv0%.1f",
             folder_name, lo, mass, tau, gay_kappa, Rbit, v0);
    const char *fname = foldername1;
    mkdir(fname, 0777);
    snprintf(foldername_ani, 128,
             "%sR%.1f_animelo%.3fMs%.3ftau%.3fkap%.3fbit%.3fv0%.1f",
             folder_name, R, lo, mass, tau, gay_kappa, Rbit, v0);
    const char *fname3 = foldername_ani;
    mkdir(fname3, 0777);
    if (!out_setup()) {
        std::cout << "file name is too long" << endl;
        return -1;
    }
    std::cout << foldername_ani << endl;
    double R_ini =
        sqrt(Np * gay_kappa /
             ((ini_rho * 2. * M_2_PI *
               (M_PI - acos(rbit_2) + rbit_2 * sqrt(1 - rbit_2 * rbit_2))) *
              2.));
    if (R_ini > R) {
        // output_ani(v,x,theta);
        for (int i = 0; i < Np; i++) {
            if ((dist2left_r(x[i], R_ini) > R_ini * R_ini && x[i][0] < 0) ||
                dist2right_r(x[i], R_ini) > R_ini * R_ini && x[i][0] > 0)
                cout << "out of range"
                     << Np * gay_kappa * 0.25 / (R_ini * R_ini) << endl;
        }
        std::cout << R_ini << endl;
        for (int j = 0; j < 1; ++j) {
            auto_list_update_anill(x, x_update, list, R_ini);
            eom_langevin(v, x, f, list, theta, f_theta, R_ini);
        }

        for (int ch = 0; ch < Np; ch++) {
            if (x[ch][0] != x[ch][0] || x[ch][1] != x[ch][1]) {
                output_ini(v, x);
                std::cout << "hazure in kasanari" << ch << endl;
                return -1;
            }
        }
        for (int i = 0; i < Np; i++) {
            if ((dist2left_r(x[i], R_ini) > R_ini * R_ini && x[i][0] < 0) ||
                dist2right_r(x[i], R_ini) > R_ini * R_ini && x[i][0] > 0)
                cout << "out of range"
                     << Np * gay_kappa * 0.25 / (R_ini * R_ini) << endl;
        }
        std::cout << "passed kasanari!" << endl;
        // output_ini(v, x);
        unsigned long long int tmaxbefch = 0;  // 10 * R_ini / (dtlg);
        for (int j = 0; j < tmaxbefch; j++) {
            auto_list_update_anill(x, x_update, list, R_ini);
            eom_langevin_h(v, x, f, list, theta, f_theta, temp, R_ini);
        }
        constexpr double moku_R = R / rat_anil;

        std::cout << "passed kakimaze!" << endl;
        constexpr double tmax_anillch = tmax_anill / dtlg;
        for (int j = 0; j < tmax_anillch; ++j) {
            auto_list_update_anill(x, x_update, list, R_ini);
            eom_langevin_h(v, x, f, list, theta, f_theta, 0.1, R_ini);
        }
        for (int i = 0; i < Np; i++) {
            if ((dist2left_r(x[i], R_ini) > R_ini * R_ini && x[i][0] < 0) ||
                dist2right_r(x[i], R_ini) > R_ini * R_ini && x[i][0] > 0)
                cout << "out of range"
                     << Np * gay_kappa * 0.25 / (R_ini * R_ini) << endl;
        }
        while (R_ini > moku_R) {
            output_ani(v, x, theta);
            R_ini *= rat_anil;
            update_coord_aniil(x);
            std::cout << R_ini << endl;
            for (int j = 0; j < tmax_anillch; ++j) {
                auto_list_update_anill(x, x_update, list, R_ini);
                eom_langevin_h(v, x, f, list, theta, f_theta, 0.1, R_ini);
            }
            for (int i = 0; i < Np; i++) {
                if ((dist2left_r(x[i], R_ini) > R_ini * R_ini && x[i][0] < 0) ||
                    dist2right_r(x[i], R_ini) > R_ini * R_ini && x[i][0] > 0)
                    cout << "out of range"
                         << Np * gay_kappa * 0.25 / (R_ini * R_ini) << endl;
            }
        }
    }
    double rat_r = R / R_ini;
    for (int i = 0; i < Np; i++) {
        x[i][0] *= rat_r;
        x[i][1] *= rat_r;
    }

    for (int i = 0; i < 1e-8; i++) {
        auto_list_update(x, x_update, list);
        eom_abp9(v, x, f, list, theta, f_theta);
    }
    output_ani(v, x, theta);
    double tmaxbefch = tmaxlg / dt;
    // if (tau * 10 > tmaxlg) tmaxbefch = tau * 10 / dt;
    for (int j = 0; j < tmaxbefch; ++j) {
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta, f_theta);
    }
    for (int ch = 0; ch < Np; ch++) {
        if (x[ch][0] != x[ch][0] || x[ch][1] != x[ch][1]) {
            output_ini(v, x);
            std::cout << "hazure in owari" << ch << endl;
        }
    }
    output_ani(v, x, theta);
    std::cout << "passed owari!" << endl;
    constexpr double tmaxch = (tau / 10 < para3_tbit)
                                  ? (tmax * (tau / (para3_tbit * 10)) / dt)
                                  : (tmax / dt),
                     tanimaxch = tmaxani / dt;
    constexpr int tanibitch = tbitani / dt,
                  para3_tbitch = (tau / 10 < para3_tbit) ? (tau / 10 / dt)
                                                         : (para3_tbit / dt);

    tout = msdini / dt;
    unsigned long long int toutcoord = tcoorch;
    unsigned long long int para3_tbitco = para3_tbitch,
                           para3_bitlonco = para3_bitlonch;
    long long int          kanit = tanibitch;
    unsigned long long int j = 0;
    double                 faimax = -5., lzmax = -5., pibarmax = -5.;
    int                    k, count_fai = 0, ituibi = 0;
    calc_corrini(x, x0, v1, v);
    calc_fai_ini();
    output_ini(v, x);
    while (j < tanimaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta, f_theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        // eom_langevin_h(v, x, f, list);
        if (j >= para3_tbitco) {
            // calc_fai(x, v, theta, j);

            para3_tbitco += para3_tbitch;
        }

        if (j >= kanit) {
            output_ani(v, x, theta);
            kanit += tanibitch;
        }
        /*
                    if (j >= toutcoord) {
                        output(v, x);
                        toutcoord += tcoorch;
                    }
                } //*/
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j);
            if (j <= msdtillch)
                tout *= 1.1;
            else
                tout += msdbit;
        }
    }
    std::cout << "ani_end" << endl;
    while (j < tmaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta, f_theta);
        // eom_langevin_h(v, x, f, list);
        if (j >= para3_tbitco) {
            // calc_fai(x, v, theta, j);

            para3_tbitco += para3_tbitch;

            /*
            if (j >= toutcoord) {
                output(v, x);
                toutcoord += tcoorch;
            } //*/
        }
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j);

            tout += msdbit;
        }
    }
    int    counthazure = 0, maxnum = 0;
    double ave;
    for (int i = 0; i < Np; ++i) {
        ave += list[i][0] / (double) Np;
        if (list[i][0] > maxnum) maxnum = list[i][0];
        if (dist2left(x[i]) > R * R || dist2right(x[i]) > R * R) counthazure++;
    }
    end = std::chrono::system_clock::now();  // 計測終了時間
    char     filename[128];
    ofstream file2;
    snprintf(filename, 128,
             "./%slo%.3fMs%.3ftau%.3fbit%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
             folder_name, lo, mass, tau, Rbit, v0, lo, mgn);
    file2.open(filename, std::ios::app);  // append
    file2 << R << " " << counthazure << " " << ave << " " << maxnum << " "
          << endl;
    file2 << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()
          << endl;  // 処理に要した時間をミリ秒に変換
    file2.close();
    if (!do_fts(foldername1)) return -1;
    std::cout << "done" << endl;
    return 0;
}
