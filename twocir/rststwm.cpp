
#include <chrono>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define lo 0.7 // コンパイル時に代入する定数;
// #define Rbit 1.4 // delta/R,Rにすると穴がなくなる;//
// コンパイル時に-D{変数名}={値}　例:-Dbit=80　とすること;
#define v0 1.

static constexpr double tau = 40.;
static constexpr double mgn = 0.;
static constexpr double R = 10.; // 固定;// ,0.1より大きいこと;
// #define Np          12800 // 4の倍数であること;NP=4*r^2*lo
#define Nn             50
#define tmax           500 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg         2000 // 緩和時間は10たうとする;
#define tmaxani        500  //>tmaxの時プログラムを変更すること;
#define tbitani        1
#define dim            2           // 変えるときはEomを変えること;
#define cut            1.122462048 // 3.
#define skin           1.5
#define dtlg           0.0001
#define dt             0.000001
#define folder_name    "stw" // 40文字程度で大きすぎ;
#define msdbit         1.1
#define msdini         0.01
#define polydispersity 0. // コードも変える;
using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.

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

static constexpr double rbit_2 = Rbit * 0.5;
static constexpr double Npd =
    (lo * M_2_PI * 2. * R * R *
     (M_PI - usr_arccos(rbit_2) + rbit_2 * usr_sqrt(1 - rbit_2 * rbit_2))) *
    2.;
static constexpr int    Np = Npd;
static constexpr double cut2 = cut * cut;
static constexpr double cut2_1 = 1. / cut2;
static constexpr double cut_1 = 1. / cut;
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double Mg = mgn * dt;
static constexpr double Np_1 = 1. / Npd;
static constexpr double center_left = -Rbit * 0.5 * R;
static constexpr double center_rignt = Rbit * 0.5 * R;
static constexpr double x0limit = R * usr_sqrt(1 - Rbit * Rbit * 0.25);

void usr_sincos(double kaku, double *x) { // x[0]がcos,x[1]issin;
                                          // 制度は10^-13程度;
    constexpr static double waru[8] = {1.0 / (3 * 4 * 5 * 6 * 7 * 8 * 9 * 10),
                                       -1.0 / (3 * 4 * 5 * 6 * 7 * 8),
                                       1.0 / (3 * 4 * 5 * 6),
                                       -1.0 / (3 * 4),
                                       1.0 / (2 * 3 * 4 * 5 * 6 * 7 * 8 * 9),
                                       -1.0 / (2 * 3 * 4 * 5 * 6 * 7),
                                       1.0 / (2 * 3 * 4 * 5),
                                       -1.0 / (2 * 3)};
    kaku *= 0.0625; // 0.03125;//kaku/=1/2^m;
    double c, s, z = kaku * kaku;
    c = ((((waru[0] * z + waru[1]) * z + waru[2]) * z + waru[3]) * z + 1.) * z;
    s = ((((waru[4] * z + waru[5]) * z + waru[6]) * z + waru[7]) * z + 1.) *
        kaku;
    for (int i = 0; i < 4; i++) { // mmade;
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

bool ini_coord_twocircles(double (*x)[dim]) {
    double rbbit = rbit_2,
           bit = sqrt((lo * R * R *
                       (M_PI - usr_arccos(rbbit) +
                        rbbit * usr_sqrt(1 - rbbit * rbbit))) *
                      2. / Np),
           R2 = R - (0.5 +
                     polydispersity * 0.5); // radiousを変える時はここを変える;
    int namari = Np % 4;
    int nmax = Np / 4, k = 0;
    for (int i = 0; i < nmax; i++) {
        for (int j = 0; j < nmax; j++) {
            double r2[2] = {0.5 + i * bit, (0.5 + j * bit)};
            if (dist2right(r2) > R2 * R2)
                break;
            x[k][0] = 0.5 + i * bit;
            x[k][1] = 0.5 + j * bit;
            x[k + nmax][0] = 0.5 + i * bit;
            x[k + nmax][1] = -0.5 - j * bit;
            x[k + 2 * nmax][0] = -0.5 - i * bit;
            x[k + 2 * nmax][1] = 0.5 + j * bit;
            x[k + 3 * nmax][0] = -0.5 - i * bit;
            x[k + 3 * nmax][1] = -0.5 - j * bit;
            k++;
            if (k >= nmax)
                break;
        }
        if (k >= nmax)
            break;
    }
    for (int i = 0; i < namari; i++) {
        x[i + 4 * nmax][0] = (double) i;
        x[i + 4 * nmax][1] = 0.;
    }

    if (4 * k + namari == Np) {
        std::cout << k * 4 << " " << Np << endl;
        return true;
    } else {
        std::cout << "passed Np is" << k * 4 + namari << " " << Np << endl;
        return false;
    }
}

void set_diameter(double *a) {
    for (int i = 0; i < Np; ++i)
        a[i] = 0.5 + polydispersity * gaussian_rand() * 0.5;
}

void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < dim; ++j)
            x[i][j] = 0.0;
}

void calc_force(double (*x)[dim], double (*f)[dim], double *a,
                int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, aij /*w12*/;
    ini_array(f);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];
            dr2 = dx * dx + dy * dy;
            aij = (a[i] + a[list[i][j]]);
            w2 = aij * aij / dr2;
            if (cut2_1 < w2) {

                w6 = w2 * w2 * w2;
                // w12 = w6 * w6;
                dUr = (-48. * w6 + 24.) * w6 / dr2 /* -12. * w12 / dr2*/;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}

void eom_abp9(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i, int timei) {
    static constexpr double zeta = 1.0, ddt = 1e-7;
    static constexpr double fluc = usr_sqrt(8. * zeta * ddt);
    double                  vi[2], ri, riw, aij, w1, w2, w6, dUr, fiw[dim];
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[1] = dUr * x[i][1];
            }
        }
        for (int j = 0; j < dim; j++) {
            v[i][j] += -zeta * v[i][j] * ddt + f[i][j] * ddt + fiw[j] * ddt +
                       fluc * gaussian_rand();
            x[i][j] += v[i][j] * ddt;
        }
    }
}

void eom_langevin(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                  double *a, int (*list)[Nn], double temmp, double *theta_i) {
    double zeta = 1.0;
    double fluc = usr_sqrt(temmp * zeta * dtlg),D = usr_sqrt(2. * dtlg / tau);
    double vi[2], ri, riw, aij, w1, w2, w6, dUr, fiw[dim];
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[1] = dUr * x[i][1];
            }
        }
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        for (int j = 0; j < dim; j++) {
            v[i][j] += -zeta * v[i][j] * dtlg + f[i][j] * dtlg + fiw[j] * dtlg +
                       fluc * gaussian_rand();
            x[i][j] += v[i][j] * dtlg;
        }
    }
}

void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double                  ri, riw, aij, w1, w2, w6, dUr, fiw[dim], sico[2];
    constexpr static double D = usr_sqrt(2. * dt / tau);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            aij = 0.5 + a[i];
            w1 = aij / riw;
            if (cut_1 < w1) {
                w2 = w1 * w1;
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[1] = dUr * x[i][1];
            }
        }
        // till here*/
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
        v[i][0] = (-v[i][0] + v0 * sico[0] + f[i][0] + fiw[0]);
        v[i][1] = (-v[i][1] + v0 * sico[1] + f[i][1] + fiw[1]);
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = 0.;
    }
}
void output_ini(double (*v)[dim], double (*x)[dim], double *a) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1f_coorlo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyokkei.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0);
    file.open(filename /* std::ios::app*/); // append
    file << tmaxani << endl;
    for (int i = 0; i < Np; i++)
        file << a[i] * 2. << endl;
    file.close();
    snprintf(filename, 128,
             "./%sR%.1f_coorlo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t0.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
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

    snprintf(filename, 128,
             "./%sR%.1f_coorlo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
void output_iniani(double (*v)[dim], double (*x)[dim], double *a) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyokkei.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0);
    file.open(filename /* std::ios::app*/); // append
    file << tbitani << endl;
    for (int i = 0; i < Np; i++)
        file << a[i] * 2. << endl;
    file.close();
    snprintf(filename, 128,
             "./%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t0.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
}
void output_ani(double (*v)[dim], double (*x)[dim]) {
    static int l = 1;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128,
             "./%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn, l);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
bool out_setup() { // filenameが１２８文字を超えていたらfalseを返す;
    char     filename[128];
    ofstream file;
    int      test =
        snprintf(filename, 128,
                 "./%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
                 "setupofst_lo%.3f_tau%.3f_m%.3f_t%d.dat",
                 folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn, tmax);
    std::cout << test << endl;
    file.open(filename, std::ios::app); // append

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
    file << "polydispersity" << polydispersity << endl;
    file.close();
    if (test == -1)
        return false;
    else
        return true;
}

double calc_fai(double (*x)[dim], double (*v)[dim]) {
    double sum_vt = 0., sum_v = 0.;
    for (int i = 0; i < Np; i++) {
        if (x > 0) {
            sum_vt += ((x[i][0] - center_rignt) * v[i][1] - x[i][1] * v[i][0]) /
                      (dist2right(x[i]));
            sum_v += sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1]);
        } else if (x < 0) {
            sum_vt += ((x[i][0] - center_left) * v[i][1] - x[i][1] * v[i][0]) /
                      (dist2left(x[i]));
            sum_v += sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1]);
        } else if (x == 0) {
            sum_v += sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1]);
        }
    }
    return (sum_vt / sum_v - M_2_PI) / (1 - M_2_PI);
}
void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v1)[dim],
               double (*v)[dim], double *xcor, double *vcor, int k,
               double *msd) {
    double dr;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < dim; ++j) {
            xcor[k] += x0[i][j] * x[i][j] * Np_1;
            vcor[k] += v1[i][j] * v[i][j] * Np_1;
            dr = x[i][j] - x0[i][j];
            msd[k] += dr * dr * Np_1;
        }
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) {
    char     filename[128];
    double   v_theta;
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "xcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "vcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f/"
             "msd_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}
inline double usr_max(double a, double b) { return (a > b) ? a : b; }
inline double usr_min(double a, double b) { return (a > b) ? b : a; }
void          cell_list(int (*list)[Nn], double (*x)[dim], double *a) {
    int                     i, j, k, l, m, lm, mm, map_index, km, nx[Np][2];
    static constexpr double cutmax = cut * (1. + polydispersity);
    double                  dx, dy;
    static constexpr double xlen_2 = (2. * R + Rbit * R) / 2.;
    static constexpr int    Mx = (int) (xlen_2 * 2. / (cutmax + skin));
    static constexpr int    My =
        (int) (2. * R / (cutmax + skin)); // M<=2R/(cutmax+skin)
    static constexpr int    m2 = Mx * My;
    static constexpr double R2 = 2. * R, bitx = Mx / (xlen_2 * 2.),
                            bity = My / (R2); // ひとつのせるの幅の逆数;
    int(*map)[Np] = new int[m2][Np];

    for (i = 0; i < m2; ++i)
        map[i][0] = 0;

    for (i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0] + xlen_2) * bitx);
        nx[i][1] = (int) ((x[i][1] + R) * bity);
        for (m = usr_max(nx[i][1] - 1, 0), mm = usr_min(nx[i][1] + 1, My - 1);
             m <= mm; ++m) {
            for (l = usr_max(nx[i][0] - 1, 0),
                lm = usr_min(nx[i][0] + 1, Mx - 1);
                 l <= lm; ++l) {
                map_index = l + Mx * m;
                map[map_index][0]++;
                map[map_index][map[map_index][0]] = i;
            }
        }
    }
    double aij;
    for (i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + Mx * nx[i][1];
        for (k = 1, km = (map[map_index][0]); k <= km; ++k) {
            j = map[map_index][k];
            if (j > i) {
                dx = x[i][0] - x[j][0];
                dy = x[i][1] - x[j][1];
                // dx-=L*floor((dx+0.5*L)/L);
                // dy-=L*floor((dy+0.5*L)/L);
                aij = a[i] + a[j];
                if ((dx * dx + dy * dy) <
                    ((cut * aij + skin) * (cut * aij + skin))) {
                    list[i][0]++;
                    list[i][list[i][0]] = j;
                }
            }
        }
    }
    delete[] map;
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

void auto_list_update(double *disp_max, double (*x)[dim],
                      double (*x_update)[dim], int (*list)[Nn], double *a) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max >= skin2) {
        cell_list(list, x, a);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        *disp_max = skinini;
        // count = 0;
    }
}

int main() {
    std::chrono::system_clock::time_point start, end; // 型は auto で可
    start = std::chrono::system_clock::now();         // 計測開始時間
    double x[Np][dim], v[Np][dim], theta[Np], a[Np], f[Np][dim], x0[Np][dim],
        v1[Np][dim], x_update[Np][dim], disp_max = 0.;
    // int(*list)[Nn] = new int[Np][Nn];
    int                    list[Np][Nn];
    int                    counthistv_theta = 0, countout = 0;
    double                 tout = msdini, toutcoord = 0;
    unsigned long long int j = 0;
    int                    k = 0, kcoord = 0;
    set_diameter(a);
    if (!ini_coord_twocircles(x))
        return -1;
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    ini_hist(theta, Np);
    char foldername[128];
    snprintf(foldername, 128, "%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f",
             folder_name, R, lo, tau, mgn, Rbit, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    snprintf(foldername, 128, "%sR%.1f_coorlo%.2ftau%.3fm%.3fbit%.3fv0%.1f",
             folder_name, R, lo, tau, mgn, Rbit, v0);
    const char *fname2 = foldername;
    mkdir(fname2, 0777);
    snprintf(foldername, 128, "%sR%.1f_animelo%.2ftau%.3fm%.3fbit%.3fv0%.1f",
             folder_name, R, lo, tau, mgn, Rbit, v0);
    const char *fname3 = foldername;
    mkdir(fname3, 0777);

    if (!out_setup()) {
        std::cout << "file name is too long" << endl;
        return -1;
    }
    std::cout << foldername << endl;

    while (tout < tmax) {

        tout *= msdbit;

        countout++;
    }
    countout += 5;
    double msd[countout], t[countout], vcor[countout], msd2[countout];
    ini_hist(msd, countout);
    ini_hist(t, countout);
    ini_hist(msd2, countout);
    ini_hist(vcor, countout);
    j = 0;

    while (j < 1e7) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_abp9(v, x, f, a, list, theta, j);
    }
    j = 0;
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output_ini(v, x, a);
            std::cout << "hazure in kasanari" << ch << endl;
            // return -1;
        }
    }
    std::cout << "passed kasanari!" << endl;
    int tmaxbefch = R / (dtlg * 5);
    while (j < tmaxbefch * 0.5) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_langevin(v, x, f, a, list, 10);
    }
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_langevin(v, x, f, a, list, 1);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output_ini(v, x, a);
            std::cout << "hazure in kakimaze" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed kakimaze!" << endl;
    j = 0;
    tmaxbefch = tmaxlg / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_abp1(v, x, f, a, list, theta);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output_ini(v, x, a);
            std::cout << "hazure in owari" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed owari!" << endl;
    int ituibi = 0, tauch = tau / dt, tanibitch = tbitani / dt;
    unsigned long long int tmaxch = tmax / dt, tanimaxch = tmaxani / dt;
    for (int xnp = 0; xnp < Np; xnp++) {

        for (int xdim = 0; xdim < dim; xdim++) {
            x0[xnp][xdim] = x[xnp][xdim];
            v1[xnp][xdim] = v[xnp][xdim];
        }
    }
    j = 0;
    tout = msdini / dt;
    toutcoord = 0.;
    k = 0;
    kcoord = 0;
    int kani = 0;
    int kanit = 0;

    calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
    t[0] = 0.;
    ++kcoord;
    output_ini(v, x, a);
    output_iniani(v, x, a);
    ++k;
    while (j < tanimaxch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        if (j >= kanit) {
            output_ani(v, x);
            kanit += tanibitch;
            /*
            if (j >= toutcoord) {
                output( v, x, a);
                toutcoord += tauch;
            }
        }//*/
            if (j >= tout) {
                calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
                t[kcoord] = j * dt;
                ++kcoord;
                tout *= msdbit;
            }
        }
    }
    while (j < tmaxch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list, a);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        /*
                if (j >= toutcoord) {
                    output( v, x, a);
                    // outtuibi(x, toutcoord, v, ituibi);
                    toutcoord += tauch;
                }//*/
        if (j >= tout) {
            calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
            t[kcoord] = j * dt;
            ++kcoord;
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
    char     filename[128];
    ofstream file2;
    snprintf(filename, 128,
             "./%sR%.1flo%.2ftau%.3fm%.3fbit%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
             folder_name, R, lo, tau, mgn, Rbit, v0, lo, mgn);
    file2.open(filename, std::ios::app); // append
    file2 << counthistv_theta << " " << counthazure << " " << ave << " "
          << maxnum << " " << endl;
    file2 << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()
          << endl; // 処理に要した時間をミリ秒に変換
    file2.close();
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
