
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <deque>
#include <vector>

#include "BM.h"


// #define Np          12800 // 4の倍数であること;NP=4*r^2*lo
#define lo          0.4 // コンパイル時に代入する定数;
#define Nn          50
#define R           20. // 固定;// ,0.1より大きいこと;
#define tmax        2000 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg      1000 // 緩和時間は10たうとする;
#define v0          1.
#define tau         100. // コンパイル時に-D{変数名}={値}　例:-Dtau=80　とすること;
#define mgn         0. // Omega=omega/tau,mass,ここではomegaを入れること;
#define tmaxani     500 //>tmaxの時プログラムを変更すること;
#define tbitani     2
#define dim         2 // 変えるときはEomを変えること;
#define ratf        1.
#define cut         1.122462048 // 3.
#define skin        1.5
#define dtlg        0.00005
#define dt          0.00005
#define folder_name "stwmss"
#define msdbit      1.2
#define msdini      0.01
#define out_para3   1.        // double
#define wh_list     cell_list // ver_list // cell_list
#define FLAG_MASS   1         // 1だと慣性あり、０ならなし;
// #define polydispersity 0.2 コードも変える;
using std::endl;
using std::ofstream;
// #define radios 1.

#ifndef Ms
#if FLAG_MASS == 1
#define Ms 0.1
#else
#define Ms 0.000000000001
#endif
#endif
static constexpr double R_in = 2.5;
static constexpr double mass = Ms;
static constexpr int    Np = 4 * (R * R - R_in * R_in) * lo;
static constexpr double cut2 = cut * cut;
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double Mg = mgn * dt;
static constexpr double Np_1 = 1. / Np;
static constexpr double R_mid = (R + R_in) / 2;
static constexpr double const_f = -48 * ratf;
void usr_sincos(double kaku, double *x) { // x[0]がcos,x[1]issin;
    constexpr static double waru[6] = {
        -1.0 / (3 * 4 * 5 * 6 * 7 * 8), 1.0 / (3 * 4 * 5 * 6), -1.0 / (3 * 4),
        -1.0 / (2 * 3 * 4 * 5 * 6 * 7), 1.0 / (2 * 3 * 4 * 5), -1.0 / (2 * 3)};
    kaku *= 0.0625; // 0.03125;//kaku/=1/2^m;
    double c, s, z = kaku * kaku;
    c = (((waru[0] * z + waru[1]) * z + waru[2]) * z + 1.) * z;
    s = (((waru[3] * z + waru[4]) * z + waru[5]) * z + 1.) * kaku;
    for (int i = 0; i < 4; i++) { // mmade;
        s = s * (2.0 - c);
        c = c * (4.0 - c);
    }
    x[0] = 1.0 - c * 0.5;
    x[1] = s;
}
constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 5000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}
bool ini_coord_circle(double (*x)[dim]) {
    int    Np_m = lo * R * R * 16 * M_1_PI;
    double num_max = sqrt(Np_m) + 1;
    double bitween = 2 * (R) / num_max, R2 = R, R1 = R_in, poj[2], r;
    int    k = 0;
    for (int j = 0; j < num_max; ++j) {
        for (int i = 0; i < num_max; ++i) {
            poj[0] = i * bitween + 0.5 - R;
            poj[1] = j * bitween + 0.5 - R;
            r = poj[0] * poj[0] + poj[1] * poj[1];
            if (r <= R2 * R2 && R1 * R1 <= r) {
                x[k][0] = poj[0];
                x[k][1] = poj[1];
                ++k;
            } else {
                continue;
            }
            if (k >= Np)
                break;
        }
        if (k >= Np)
            break;
    }
    if (k != Np) {
        std::cerr << k << " " << Np << endl;
        return false;
    }

    else {
        std::cout << k << endl;
        return true;
    }
}

void set_diameter(double *a) {
    for (int i = 0; i < Np; ++i)
        a[i] = 0.5;
}

void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < dim; ++j)
            x[i][j] = 0.0;
}

void calc_force(double (*x)[dim], double (*f)[dim], int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, /*w12,*/ aij;
    ini_array(f);

    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];
            dr2 = dx * dx + dy * dy;
            if (dr2 < cut2) {
                w2 = 1. / dr2;
                w6 = w2 * w2 * w2;
                // w12 = w6 * w6;
                dUr = const_f * (w6 - 0.5) * w6 * w2 /* -12. * w12 / dr2*/;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}
inline void calc_force_2w(double *x, double *f) {
    double r2, dr, w2, w6, dUr = 0.;
    r2 = sqrt(x[0] * x[0] + x[1] * x[1]);
    dr = R + 0.5 - r2;
    if (dr < cut) {
        w2 = 1. / (dr * dr);
        w6 = w2 * w2 * w2;
        dUr = const_f * (w6 - 0.5) * w6 / (dr * r2);
    }
    f[0] = dUr * x[0];
    f[1] = dUr * x[1];
    dr = r2 - R_in + 0.5;
    if (dr < cut) {
        w2 = 1. / (dr * dr);
        w6 = w2 * w2 * w2;
        dUr = const_f * (w6 - 0.5) * w6 / (dr * r2);
        f[0] -= dUr * x[0];
        f[1] -= dUr * x[1];
    }
}

void eom_langevin(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                  int (*list)[Nn]) {
    double ddt = 0.0000001, fiw[dim], fluc = sqrt(6 * ddt);
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        for (int j = 0; j < dim; j++) {
            v[i][j] +=
                (-v[i][j] + f[i][j] + fiw[j]) * ddt + fluc * gaussian_rand();
            x[i][j] += v[i][j] * ddt;
        }
    }
}
void eom_langevin_t(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                    int (*list)[Nn]) {
    double fiw[dim], fluc = sqrt(2 * dtlg);
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        for (int j = 0; j < dim; j++) {
            v[i][j] +=
                (-v[i][j] + f[i][j] + fiw[j]) * dtlg + fluc * gaussian_rand();
            x[i][j] += v[i][j] * dtlg;
        }
    }
}
void eom_abp9(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i) {
    double           fiw[dim], sico[2];
    constexpr double ddt = 1e-7, D = usr_sqrt(2. * ddt / tau),
                     M_inv = ddt / mass;
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
#if FLAG_MASS
        v[i][0] += (-v[i][0] + v0 * sico[0] + f[i][0] + fiw[0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sico[1] + f[i][1] + fiw[1]) * M_inv;
#else
        v[i][0] = (v0 * sico[0] + f[i][0] + fiw[0]);
        v[i][1] = (v0 * sico[1] + f[i][1] + fiw[1]);
#endif
        x[i][0] += v[i][0] * ddt;
        x[i][1] += v[i][1] * ddt;
    }
}
void eom_abp8(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i) {
    double           fiw[dim], sico[2];
    constexpr double D = usr_sqrt(2. * dtlg / tau), M_inv = dtlg / mass;
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
#if FLAG_MASS
        v[i][0] += (-v[i][0] + v0 * sico[0] + f[i][0] + fiw[0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sico[1] + f[i][1] + fiw[1]) * M_inv;
#else
        v[i][0] = (v0 * sico[0] + f[i][0] + fiw[0]);
        v[i][1] = (v0 * sico[1] + f[i][1] + fiw[1]);
#endif
        x[i][0] += v[i][0] * dtlg;
        x[i][1] += v[i][1] * dtlg;
    }
}

void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i) {
    double              ri, riw, aij, w2, w6, dUr, fiw[dim];
    static const double D = sqrt(2. * dt / tau), M_inv = dt / mass;
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
#if FLAG_MASS
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0] + fiw[0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1] + fiw[1]) * M_inv;
#else
        v[i][0] = (v0 * cos(theta_i[i]) + f[i][0] + fiw[0]);
        v[i][1] = (v0 * sin(theta_i[i]) + f[i][1] + fiw[1]);
#endif
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = 0.;
    }
}
void calc_fai_ini() {

    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "fais_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename); // append
    file << "t fai pibar vtheta" << endl;
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "om_lim_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename);
    file << "t om_out om_in" << endl;
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "haikou_theta_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename);
    file << "t haikou_out haikou_in" << endl;
    file.close();
}
void calc_fai(double (*x)[dim], double (*v)[dim], double *theta_i,
              long long j) { // para[0]:fai para[1]:vt* para[2];om*;
    double r2, r, vt, sum_v = 0., sum_vt = 0., fai = 0., pibar = 0.,
                      om_out = 0., om_in = 0., haikou_out = 0., haikou_in = 0.;
    int              cnt_out = 0, cnt_in = 0;
    constexpr double bun_fai = 1. / (1 - M_2_PI);
    for (int i = 0; i < Np; i++) {
        r2 = x[i][0] * x[i][0] + x[i][1] * x[i][1];
        r = sqrt(r);
        vt = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / r2;
        sum_vt += vt * r * Np_1;
        sum_v += sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1]) * Np_1;
        pibar += ((vt > 0) - (vt < 0)) * Np_1;
        if (r >= R - 1) {
            om_out += vt;
            haikou_out +=
                M_PI2 * (int) ((theta_i[i] - atan2(x[i][1], x[i][0])) * M_1_PI);
            cnt_out++;
        }
        if (r <= R_in + 1) {
            om_in += vt;
            haikou_in +=
                M_PI2 * (int) ((theta_i[i] - atan2(x[i][1], x[i][0])) * M_1_PI);
            cnt_in++;
        }
    }
    fai = (sum_vt / sum_v - M_2_PI) * bun_fai;
    double in_cnt = 1. / cnt_in, out_cnt = 1. / cnt_out;
    om_out *= out_cnt;
    haikou_out *= out_cnt;
    om_in *= in_cnt;
    haikou_in *= in_cnt;
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "fais_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << fai << "\t" << pibar << "\t" << sum_vt << "\t"
         << endl;
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "om_lim_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename, std::ios::app);
    file << j * dt << "\t" << om_out << "\t" << om_in << endl;
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "haikou_theta_R%.3f.dat",
             folder_name, lo, mass, tau, 0, v0, R);
    file.open(filename, std::ios::app);
    file << j * dt << "\t" << haikou_out << "\t" << haikou_in << endl;
    file.close();
}
void output(double (*v)[dim], double (*x)[dim]) {
    static int l = 0;
    char       filename[128];
    ofstream   file;
    sprintf(filename,
            "./%sR%.1frs%.1f_coorlo%.2ftau%.3fMs%.2fv0%.1f/"
            "tyouwaenn_m%.3f_t%d.dat",
            folder_name, R, R_in, lo, tau, mass, v0, mgn, l);
    file.open(filename); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
void output_ani_ini(double (*v)[dim], double (*x)[dim]) {
    char     filename[128];
    ofstream file;
    sprintf(filename,
            "./%sR%.1frs%.1f_animelo%.2ftau%.3fMs%.3fv0%.1f/"
            "tyouwaenn_m%.3f_t%d.dat",
            folder_name, R, R_in, lo, tau, mass, v0, mgn, 0);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    sprintf(filename,
            "./%sR%.1frs%.1f_animelo%.2ftau%.3fMs%.3fv0%.1f/"
            "tyokkei_m%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, mgn);
    file.open(filename /* std::ios::app*/); // append
    file << tbitani << endl;
    for (int i = 0; i < Np; ++i) {
        file << 1 << endl;
    }
}
void output_ani(double (*v)[dim], double (*x)[dim]) {
    static int l = 1;
    char       filename[128];
    ofstream   file;
    sprintf(filename,
            "./%sR%.1frs%.1f_animelo%.2ftau%.3fMs%.3fv0%.1f/"
            "tyouwaenn_m%.3f_t%d.dat",
            folder_name, R, R_in, lo, tau, mass, v0, mgn, l);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
void out_setup() {
    char     filename[128];
    ofstream file;
    std::cout << snprintf(
        filename, 128,
        "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/setupr%fm%f.dat", folder_name,
        R, R_in, lo, tau, mass, v0, R, mgn);
    file.open(filename, std::ios::app); // append

    file << "dt=" << dt << endl;
    file << "cut" << cut << endl;
    file << "skin" << skin << endl;
    file << "Nn" << Nn << endl;
    file << "Np=" << Np << endl;
    file << "tmax=" << tmax << endl;
    file << "tmaxlg=" << tmaxlg << endl;
    file << "v0=" << v0 << endl;
    // file << "ens=" << ensemble << endl;
    file << "type=" << 2 << endl;
    file << "2DkaraD" << endl;
    file << "壁はWCA" << endl;
    file << "cell list" << endl;
    file << "usr_sincos" << endl;
    file.close();
}

void ini_corr(double (*x)[dim], double (*v)[dim], double (*x1)[dim],
              double (*v1)[dim]) {
    double xcor = 0., vcor = 0.;
    for (int i = 0; i < Np; i++) {
        x1[i][0] = x[i][0];
        x1[i][1] = x[i][1];
        xcor += (x[i][0] * x[i][0] + x[i][1] * x[i][1]) * Np_1;
        v1[i][0] = v[i][0];
        v1[i][1] = v[i][1];
        vcor += (v[i][0] * v[i][0] + v[i][1] * v[i][1]) * Np_1;
    }
    char     filename[128];
    ofstream file;
    sprintf(filename,
            "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/"
            "xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << 0 << "\t" << xcor << endl;
    file.close();
    sprintf(filename,
            "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/"
            "vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << 0 << "\t" << vcor << endl;
    file.close();
    sprintf(
        filename,
        "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/msd_lo%.3f_tau%.3f_m%.3f.dat",
        folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << 0 << "\t" << 0 << endl;
    file.close();
}
void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v1)[dim],
               double (*v)[dim], unsigned long long int j) {
    double dr, xcor = 0., vcor = 0., msd = 0.;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < dim; ++j) {
            xcor += x0[i][j] * x[i][j] * Np_1;
            vcor += v1[i][j] * v[i][j] * Np_1;
            dr = x[i][j] - x0[i][j];
            msd += dr * dr * Np_1;
        }
    }
    char     filename[128];
    ofstream file;
    sprintf(filename,
            "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/"
            "xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << xcor << endl;
    file.close();
    sprintf(filename,
            "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/"
            "vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << vcor << endl;
    file.close();
    sprintf(
        filename,
        "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/msd_lo%.3f_tau%.3f_m%.3f.dat",
        folder_name, R, R_in, lo, tau, mass, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << msd << endl;
    file.close();
}
inline int usr_max(int a, int b) { return ((a > b) ? a : b); }
inline int usr_min(int a, int b) { return ((a > b) ? b : a); }
void       cell_list(int (*list)[Nn], double (*x)[dim]) {
    int                     map_index, nx[Np][dim];
    static constexpr double Rbit = 0., xlen_2 = (2. * R + Rbit * R) / 2.,
                            threash2 = (cut + skin) * (cut + skin);
    static constexpr int Mx = (int) (xlen_2 * 2. / (cut + skin));
    static constexpr int My =
        (int) (2. * R / (cut + skin)); // M<=2R/(cutmax+skin)
    static constexpr int    m2 = Mx * My;
    static constexpr double R2 = 2. * R, bitx = Mx / (xlen_2 * 2.),
                            bity = My / (R2); // ひとつのせるの幅の逆数;
    // int(*map)[Np + 1] = new int[m2][Np + 1];
    std::vector<std::deque<int>> map(m2);

    for (int i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0] + xlen_2) * bitx);
        nx[i][1] = (int) ((x[i][1] + R) * bity);
        for (int m = usr_max(nx[i][1] - 1, 0),
                 mm = usr_min(nx[i][1] + 1, My - 1);
             m <= mm; ++m) {
            for (int l = usr_max(nx[i][0] - 1, 0),
                     lm = usr_min(nx[i][0] + 1, Mx - 1);
                 l <= lm; ++l) {
                map_index = l + Mx * m;
                map[map_index].emplace_front(i);
            }
        }
    }
    double dx, dy;
    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + Mx * nx[i][1];
        for (auto &j : map[map_index]) {
            // j = map[map_index][k];
            if (j > i) {
                dx = (x[i][0] - x[j][0]);
                dy = (x[i][1] - x[j][1]);
                if ((dx * dx + dy * dy) < threash2) {
                    list[i][0]++;
                    list[i][list[i][0]] = j;
                }
            }
        }
    }
    // delete[] map;
}
void ver_list(int (*list)[Nn], double (*x)[dim]) {
    double           dx, dy, dr2;
    constexpr double thresh2 = (cut + skin) * (cut + skin);
    for (int i = 0; i < Np; i++)
        list[i][0] = 0;

    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Np; j++) {
            if (j > i) {
                dx = x[i][0] - x[j][0];
                dy = x[i][1] - x[j][1];
                dr2 = dx * dx + dy * dy;
                if (dr2 < thresh2) {
                    list[i][0]++;
                    list[i][(int) list[i][0]] = j;
                }
            }
        }
}
void update(double (*x_update)[dim], double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x_update[i][j] = x[i][j];
}

void calc_disp_max(double *disp_max, double (*x)[dim],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp;
    for (int i = 0; i < Np; i++) {
        dx = x[i][0] - x_update[i][0];
        dy = x[i][1] - x_update[i][1];
        disp = dx * dx + dy * dy;
        if (disp > *disp_max)
            *disp_max = disp;
    }
}

void auto_list_update(double (*x)[dim], double (*x_update)[dim],
                      int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    static double    disp_max = skin2 + 100;
    calc_disp_max(&(disp_max), x, x_update);
    if (disp_max > skin2) {
        wh_list(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        disp_max = skinini;
        // count = 0;
    }
}

int main() {
    std::chrono::system_clock::time_point start, end; // 型は auto で可
    start = std::chrono::system_clock::now();         // 計測開始時間
    double x[Np][dim], v[Np][dim], theta[Np], f[Np][dim], x0[Np][dim],
        v1[Np][dim], x_update[Np][dim];
    // int(*list)[Nn] = new int[Np][Nn];
    int list[Np][Nn];
    int counthistv_theta = 0, countout = 0;
    if (!ini_coord_circle(x))
        return -1;
    ini_array(v);
    ini_array(f);
    ini_hist(theta, Np);
    char foldername[128];
    sprintf(foldername, "%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f", folder_name,
            R, R_in, lo, tau, mass, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%sR%.1frs%.1f_coorlo%.2ftau%.3fMs%.2fv0%.1f",
            folder_name, R, R_in, lo, tau, mass, v0);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    sprintf(foldername, "%sR%.1frs%.1f_animelo%.2ftau%.3fMs%.3fv0%.1f",
            folder_name, R, R_in, lo, tau, mass, v0);
    const char *fname3 = foldername;
    mkdir(fname3, 0777);

    out_setup();
    std::cout << foldername << endl;

    for (int j = 0; j < 1e7; j++) {
        auto_list_update(x, x_update, list);
        eom_langevin(v, x, f, list);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] * x[ch][0] + x[ch][1] * x[ch][1]) > (R * R)) {
            output(v, x);
            std::cout << "hazure in kasanari" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed kasanari!" << endl;
    for (int j = 0, jmax = R / dtlg; j < jmax; j++) {
        auto_list_update(x, x_update, list);
        eom_langevin_t(v, x, f, list);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] * x[ch][0] + x[ch][1] * x[ch][1]) > (R * R)) {
            output(v, x);
            std::cout << "hazure in kakimaze" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed kakimaze!" << endl;
    for (int j = 0; j < 1e7; j++) {
        auto_list_update(x, x_update, list);
        eom_abp9(v, x, f, list, theta);
    }

    int tmaxbefch = R / (dtlg * 5);
    for (int j = 0; j < tmaxbefch; j++) {
        auto_list_update(x, x_update, list);
        eom_abp8(v, x, f, list, theta);
    }

    tmaxbefch = tmaxlg / dt;
    for (int j = 0; j < tmaxbefch; j++) {
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] * x[ch][0] + x[ch][1] * x[ch][1]) > (R * R)) {
            output(v, x);
            std::cout << "hazure in owari" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed owari!" << endl;
    int ituibi = 0, tauch = tau / dt, tmaxch = tmax / dt,
        tanimaxch = tmaxani / dt, tanibitch = tbitani / dt,
        out_para3ch = out_para3 / dt;
    unsigned long long j = 0;
    double             tout = msdini / dt, toutcoord = 0, out_para3_seki = 0;
    long long int      kanit = 0;
    ini_corr(x, v, x0, v1);
    output_ani_ini(v, x);
    calc_fai_ini();
    while (j < tanimaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta);
        if (j >= out_para3_seki) {
            out_para3_seki += out_para3ch;
            calc_fai(x, v, theta, j);
            if (j >= kanit) {
                output_ani(v, x);
                kanit += tanibitch;
                if (j >= toutcoord) {
                    output(v, x);
                    toutcoord += tauch;
                }
            }
        }
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j);
            tout *= msdbit;
        }
    }
    while (j < tmaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta);
        if (j >= out_para3_seki) {
            out_para3_seki += out_para3ch;
            calc_fai(x, v, theta, j);
            if (j >= toutcoord) {
                output(v, x);
                toutcoord += tauch;
            }
        }
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j);
            tout *= msdbit;
        }
    }
    int    counthazure = 0, maxnum = 0;
    double ave;
    for (int i = 0; i < Np; ++i) {
        ave += list[i][0] / (double) Np;
        if (list[i][0] > maxnum)
            maxnum = list[i][0];
        if (x[i][0] * x[i][0] + x[i][1] * x[i][1] > R * R)
            counthazure++;
    }
    end = std::chrono::system_clock::now(); // 計測終了時間
    char filename[128];

    ofstream file;

    sprintf(filename,
            "./%sR%.1frs%.1flo%.2ftau%.3fMs%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
            folder_name, R, R_in, lo, tau, mass, v0, lo, mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << " " << endl;
    file << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()
         << endl; // 処理に要した時間をミリ秒に変換
    file.close();

    std::cout << "done" << endl;
    return 0;
}
