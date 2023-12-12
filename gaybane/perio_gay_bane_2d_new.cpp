
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

#define Np 30  // 累乗であること;
#define Nn 31

#define tmaxlg      0
#define tmaxkau     0     // 緩和時間;
#define tmaxani     50 //>tmaxの時プログラムを変更すること;
#define tmax        tmaxani
#define tbitani     0.1
#define ratf        1.
#define dim         2  // 変えるときはEomを変えること;
#define skin        1.5
#define gay_kappa   3.25
#define dtlg        1e-8
#define dt          1e-6
#define folder_name "gay_t"  // 40文字程度で大きすぎ;
#define rat         3      // Lx=Ly*rat;
#define msdbit      1.2
#define msdini      0.01
#define FLAG_MASS   0  // 1なら慣性あり;
#define FLAG_FORCE  0  // 1ならharmonic,0ならWCA;
#if FLAG_FORCE == 1
#define cut 1
#else
#define cut 1.122462048  // 3.; 
#endif
// #define polydispersity 0. // コードも変える;
using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.
#define lo 0.55  // コンパイル時に代入する定数;
// コンパイル時に-D{変数名}={値}　例:-Dbit=80　とすること;
#define v0 1.
#if !defined(MS)
#if FLAG_MASS == 1
#define MS 10000
#else
#define MS 0.000000000001
#endif
#endif  // MS
#if !defined(TAU)
#define TAU 0.0000000001
#endif  // TAU

static constexpr double tau = TAU;
static constexpr double mass = MS;
static constexpr double mgn = 0.;  // 有限にするときはコードを変える;
constexpr double cone_abs(double x){
    return (x>0)?x:-x;
}
constexpr double usr_sqrt(double x,bool flag) {
    if(cone_abs(2*dt/x-tau)<1e-6&&tau<=1e-6&&flag)return 0;
    double b = x;
    for (int i = 0; i < 5000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double Mg = mgn * dt;
static constexpr double Ly = usr_sqrt(Np * M_PI * 0.25 * gay_kappa / lo / rat,false);
static constexpr double Lx = Ly * rat;

static constexpr double Np_1 = 1. / Np;
static constexpr double fconst = ratf * 48.;
static constexpr int    dim2 = 2 * dim;
static constexpr double to = 10;     // (tau > 10) ? tau * 2. : 20;
static constexpr double tmaxka = 0;  // (tmaxkau > 5 * tau) ? tmaxkau : 5 * tau;

char foldername_out[128];
char foldername_corr[128];
char foldername_ani[128];
#include "gay_bane_perio_new.hpp"
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

void ini_hex(double (*x)[dim2]) {
    int num_x =Lx/(gay_kappa+cut-1);
    int num_y = Np/num_x;
    if (num_x * num_y != Np) {
        num_y++;
    }
    int    k = 0;
    double shift=0;
    for (int j = 0; j < num_y; j++) {
        for (int i = 0; i < num_x; i++) {
            // shift = (double) j * 0.5 - j / 2;
            x[i + num_x * j][0] = (shift + i) * Lx / (double) num_x;
            x[i + num_x * j][1] = j * Ly / (double) num_y;
            k++;
            if (k == Np) {
                break;
            }
        }
        if (k == Np) {
            break;
        }
    }
}

void set_diameter(double *a) {
    for (int i = 0; i < Np; ++i) a[i] = 0.5;
}

void ini_hist(double *his, int koo) {
    for (int i = 0; i < koo; i++) {
        his[i] = 0.;
    }
}
void ini_theta(double *theta) {
    for (int i = 0; i < Np; i++) {
        theta[i] = M_PI +M_PI*i+0.1;
        theta[i] -= (floor(theta[i] * M_1_PI * 0.5)) * M_PI2;
    }
}
inline double perio_x(double x) { return Lx * floor(x * Lx_inv); }
inline double perio_y(double x) { return Ly * floor(x * Ly_inv); }

void eom_abp9(double (*v)[dim], double (*x)[dim2], double (*f)[dim],
              double *f_theta, int (*list)[Nn], double *theta_i) {
    constexpr double ddt = 1e-7;
    double           sico[2];
    constexpr double D = usr_sqrt(2. * ddt / tau,true), M_inv = ddt / mass;
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        theta_i[i] += D * gaussian_rand() + f_theta[i] * 3 * ddt;
        theta_i[i] -= (floor(theta_i[i] * M_1_PI * 0.5)) * M_PI2;
        usr_sincos(theta_i[i], sico);
#if FLAG_MASS == 1
        v[i][0] += (-v[i][0] + v0 * sico[0] + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sico[1] + f[i][1]) * M_inv;
#else
        v[i][0] = v0 * sico[0] + f[i][0];
        v[i][1] = v0 * sico[1] + f[i][1];
#endif
        x[i][0] += v[i][0] * ddt;
        x[i][1] += v[i][1] * ddt;
        x[i][0] -= perio_x(x[i][0]);
        x[i][1] -= perio_y(x[i][1]);
    }
}

void ini_count(double (*x)[dim2]) {
    for (int i = 0; i < Np; i++) {
        x[i][2] = 0.;
        x[i][3] = 0.;
    }
}
void eom_abp1(double (*v)[dim], double (*x)[dim2], double (*f)[dim],
              double *f_theta, int (*list)[Nn], double *theta_i) {
    double           ov[2];
    constexpr double D  = usr_sqrt(2. * dt / tau,true), M_inv = dt / mass;
    calc_force_gay_bane(x, f, list, theta_i, f_theta);
    for (int i = 0; i < Np; i++) {
        theta_i[i] += D * gaussian_rand() + f_theta[i] * 3 * dt;
        theta_i[i] -= (floor(theta_i[i] * M_1_PI * 0.5)) * M_PI2;
// usr_sincos(theta_i[i], sico);
#if FLAG_MASS == 1
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1]) * M_inv;
#else
        v[i][0] = v0 * cos(theta_i[i]) + f[i][0];
        v[i][1] = v0 * sin(theta_i[i]) + f[i][1];
#endif
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
        ov[0] = -perio_x(x[i][0]);
        ov[1] = -perio_y(x[i][1]);
        x[i][0] += ov[0];
        x[i][1] += ov[1];
        x[i][2] += ov[0];
        x[i][3] += ov[1];
    }
}
inline double usr_abs(double x) { return x * ((x > 0) - (x < 0)); }

void output_ini(double (*v)[dim], double (*x)[dim2], double *theta) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s/"
             "tyokkei_m%.3f.dat",
             foldername_ani, mgn);
    file.open(filename /* std::ios::app*/);  // append
    file << tbitani << endl;
    file << gay_kappa << endl;
    for (int i = 0; i < Np; i++) file << 1 << endl;
    file.close();
    snprintf(filename, 128,
             "./%s/"
             "elli_m%.3f_t0.dat",
             foldername_corr, mgn);
    file.open(filename);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << "\t" << theta[i] << endl;
    }
    file.close();
}
void output(double (*v)[dim], double (*x)[dim2], double *theta) {
    static int l = 0;
    char       filename[128];
    ofstream   file;
    snprintf(filename, 128, "./%s/elli_m%.3f_t%d.dat", foldername_corr, mgn, l);
    file.open(filename);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << "\t" << theta[i] << endl;
    }
    file.close();
    l++;
}
void output_iniani(double (*v)[dim], double (*x)[dim2], double *theta,double *f_theta) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s/"
             "tyokkei_%.3f.dat",
             foldername_ani, mgn);
    file.open(filename /* std::ios::app*/);  // append
    file << tbitani << endl;
    file << gay_kappa << endl;
    for (int i = 0; i < Np; i++) file << 1 << endl;
    file.close();
    snprintf(filename, 128,
             "./%s/"
             "elli_m%.3f_t0.dat",
             foldername_ani, mgn);
    file.open(filename /* std::ios::app*/);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << "\t" << theta[i]<<"\t"<<f_theta[i] << endl;
    }
    file.close();
}
void output_ani(double (*v)[dim], double (*x)[dim2], double *theta,double *f_theta) {
    static int l = 1;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128, "./%s/elli_m%.3f_t%d.dat", foldername_ani, mgn, l);
    file.open(filename /* std::ios::app*/);  // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << "\t" << theta[i]<<"\t"<<f_theta[i] << endl;
    }
    file.close();
    l++;
}

bool out_setup() {  // filenameが１２８文字を超えていたらfalseを返す;
    char     filename[128];
    ofstream file;
    int      test = snprintf(filename, 128,
                             "./%s/"
                                  "setupofst_lo%.3f_tau%.3f_m%.3f_t%d.dat",
                             foldername_out, lo, tau, mgn, tmax);
    std::cout << test << endl;
    file.open(filename, std::ios::app);  // append

    file << "dt=" << dt << endl;
    file << "cut" << cut << endl;
    file << "skin" << skin << endl;
    file << "Nn" << Nn << endl;
    file << "Np=" << Np << endl;
    file << "tmax=" << tmax << endl;
    file << "tmaxlg=" << tmaxlg << endl;
    file << "v0=" << v0 << endl;
    file << "2DkaraD" << endl;
    file << "壁はWCA" << endl;
    file << "cell list" << endl;
    file << "usr_sincos" << endl;
    file << "自動Np" << endl;
    file << "x=0での壁を追加" << endl;
    file << "modNp" << endl;
    file << "L=" << Lx << endl;
    // file << "pi? " << L * L * lo / Np << endl;
    file << "ratio" << ratf << endl;
    file << "to" << to << endl;
    file.close();
    if (test == -1)
        return false;
    else
        return true;
}
void ini_corr(double (*x)[dim2], double (*x0)[dim], double (*v1)[dim],
              double (*v)[dim]) {
    double xcor = 0., vcor = 0.;
    for (int i = 0; i < Np; i++) {
        x0[i][0] = x[i][0];
        x0[i][1] = x[i][1];
        xcor += (x[i][0] * x[i][0] + x[i][1] * x[i][1]) * Np_1;
        x[i][2] = 0.;
        x[i][3] = 0.;
        v1[i][0] = v[i][0];
        v1[i][1] = v[i][1];
        vcor += (v[i][0] * v[i][0] + v[i][1] * v[i][1]) * Np_1;
    }
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "./%s/xcor_m%.3f.dat", foldername_out, mgn);
    file.open(filename);  // append
    file << 0 << "\t" << xcor << endl;
    file.close();
    snprintf(filename, 128, "./%s/vcor_m%.3f.dat", foldername_out, mgn);
    file.open(filename);  // append
    file << 0 << "\t" << vcor << endl;
    file.close();
    snprintf(filename, 128, "./%s/msd_m%.3f.dat", foldername_out, mgn);
    file.open(filename);  // append
    file << "#0 <<  << msd" << endl;
    file.close();
    snprintf(filename, 128, "./%s/msd2_m%.3f.dat", foldername_out, mgn);
    file.open(filename);  // append
    file << "#0 <<  << msd2" << endl;
    file.close();
}
void calc_corr(double (*x)[dim2], double (*x0)[dim], double (*v1)[dim],
               double (*v)[dim], unsigned long long int j) {
    double gc[2] = {0., 0.};
    for (int i = 0; i < Np; ++i) {
        gc[0] += (x[i][0] - x[i][2]) * Np_1;
        gc[1] += (x[i][1] - x[i][3]) * Np_1;
    }
    double xcor = 0., vcor = 0., msd = 0., msd2 = 0.;
    for (int i = 0; i < Np; i++) {
        double rsin[2] = {0., 0.};
        rsin[0] = x[i][0] - x[i][2] - gc[0];
        rsin[1] = x[i][1] - x[i][3] - gc[1];
        xcor += (x0[i][0] * rsin[0] + x0[i][1] * rsin[1]) * Np_1;
        vcor += (v1[i][0] * v[i][0] + v1[i][1] * v[i][1]) * Np_1;
        msd += ((rsin[0] - x0[i][0]) * (rsin[0] - x0[i][0]) +
                (rsin[1] - x0[i][1]) * (rsin[1] - x0[i][1])) *
               Np_1;
        msd2 = ((rsin[0] - x0[i][0] + gc[0]) * (rsin[0] - x0[i][0] + gc[0]) +
                (rsin[1] - x0[i][1] + gc[1]) * (rsin[1] - x0[i][1] + gc[1])) *
               Np_1;
    }
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "./%s/xcor_m%.3f.dat", foldername_out, mgn);
    file.open(filename, std::ios::app);  // append
    file << j * dt << "\t" << xcor << endl;
    file.close();
    snprintf(filename, 128, "./%s/vcor_m%.3f.dat", foldername_out, mgn);
    file.open(filename, std::ios::app);  // append
    file << j * dt << "\t" << vcor << endl;
    file.close();
    snprintf(filename, 128, "./%s/msd_m%.3f.dat", foldername_out, mgn);
    file.open(filename, std::ios::app);  // append
    file << j * dt << "\t" << msd << endl;
    file.close();
    snprintf(filename, 128, "./%s/msd2_m%.3f.dat", foldername_out, mgn);
    file.open(filename, std::ios::app);  // append
    file << j * dt << "\t" << msd2 << endl;
    file.close();
}

int main() {
    std::chrono::system_clock::time_point start, end;  // 型は auto で可
    start = std::chrono::system_clock::now();          // 計測開始時間
    double x[Np][dim2], v[Np][dim], theta[Np], f[Np][dim], f_theta[Np],
        x0[Np][dim], v1[Np][dim], x_update[Np][dim];
    // int(*list)[Nn] = new int[Np][Nn];
    int list[Np][Nn];
    ini_hex(x);
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    ini_theta(theta);
    snprintf(foldername_out, 128, "%skap%.2flo%.2fMs%.3ftau%.3fv0%.1f",
             folder_name, gay_kappa, lo, mass, tau, v0);
    const char *fname = foldername_out;
    mkdir(fname, 0777);
    snprintf(foldername_corr, 128, "%s_coorlo%.2fMs%.3ftau%.3fv0%.1f",
             folder_name, lo, mass, tau, v0);
    const char *fname2 = foldername_corr;
    mkdir(fname2, 0777);
    snprintf(foldername_ani, 128, "%s_animekap%.2flo%.2fMs%.3ftau%.3fv0%.1f",
             folder_name, gay_kappa, lo, mass, tau, v0);
    const char *fname3 = foldername_ani;
    mkdir(fname3, 0777);

    if (!out_setup()) {
        std::cout << "file name is too long" << endl;
        return -1;
    }
    std::cout << foldername_out << endl;
    // for (int j = 0; j < 1e7; ++j) {
    //     auto_list_update(x, x_update, list);
    //     eom_abp9(v, x, f, f_theta, list, theta);
    // }

    std::cout << "passed kakimaze!" << endl;

    int tmaxbefch = tmaxka / dt;
    for (long long int j = 0; j < tmaxbefch; ++j) {
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, f_theta, list, theta);
    }
    std::cout << "passed owari!" << endl;
    double ituibi = 0, toch = to / dt, tanibitch = tbitani / dt;
    constexpr unsigned long long int tmaxch = tmax / dt,
                                     tanimaxch = tmaxani / dt;
    double tout = msdini / dt;
    double toutcoord = 0.;

    int kanit = 0;
    ini_count(x);

    ini_corr(x, x0, v1, v);
    output_ini(v, x, theta);
    output_iniani(v, x, theta,f_theta);

    for (unsigned long long int j = 0; j < tanimaxch; ++j) {
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, f_theta, list, theta);
        double oouutt=j*dt;
        // make_v_thetahist(x, v, hist, hist2, lohist);
        if (j >= kanit) {
            output_ani(v, x, theta,f_theta);
            kanit += tanibitch;
            std::cout<<j*dt<<endl;
            if (j >= toutcoord) {
                output(v, x, theta);
                toutcoord += toch;
            }
        }  //*/
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j);
            tout *= msdbit;
        }
    }
    std::cout << "ani_end" << endl;
    for (unsigned long long int j = 0, tmaxch2 = tmaxch - tanimaxch;
         j < tmaxch2; ++j) {
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, f_theta, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);

        if (j >= toutcoord) {
            output(v, x, theta);
            toutcoord += toch;
            //*/
        }
        if (j >= tout) {
            calc_corr(x, x0, v1, v, j + tmaxch);
            tout *= msdbit;
        }
    }
    int    counthazure = 0, maxnum = 0;
    double ave;
    for (int i = 0; i < Np; ++i) {
        ave += list[i][0] / (double) Np;
        if (list[i][0] > maxnum) maxnum = list[i][0];
    }
    end = std::chrono::system_clock::now();  // 計測終了時間
    char     filename[128];
    ofstream file2;
    snprintf(filename, 128, "./%s/kekkalo%.3fm%.3f.dat", foldername_out, lo,
             mgn);
    file2.open(filename, std::ios::app);  // append
    file2 << ave << " " << maxnum << " " << endl;
    file2 << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()
          << endl;  // 処理に要した時間をミリ秒に変換
    file2.close();
    std::cout << "done" << endl;
    return 0;
}
