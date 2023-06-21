
#include <chrono>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

// #define Np          12800 // 4の倍数であること;NP=4*r^2*lo
#define lo          0.1 // コンパイル時に代入する定数;
#define Nn          100
#define R           10. // 固定;// ,0.1より大きいこと;
#define tmax        500 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg      800 // 緩和時間は10たうとする;
// #define Rbit        1.8 // delta/R,Rにすると穴がなくなる;
#define v0          1.
#define tau         40. // コンパイル時に-D{変数名}={値}　例:-Dtau=80　とすること;
#define mgn         0.  // Omega=omega/tau,ここではomegaを入れること;
#define tmaxani     500 //>tmaxの時プログラムを変更すること;
#define tbitani     1
#define dim         2           // 変えるときはEomを変えること;
#define cut         1.122462048 // 3.
#define skin        1.5
#define dtlg        0.0001
#define dt          0.000001
#define folder_name "stwmss" // 40文字程度で大きすぎ;
#define msdbit      1.1
#define msdini      0.01
#define M_ss        20
// #define polydispersity 0.2 コードも変える;
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
    double rbbit = rbit_2 - 0.15,
           bit = sqrt((R * R *
                       (M_PI - usr_arccos(rbbit) +
                        rbbit * usr_sqrt(1 - rbbit * rbbit))) *
                      2. / Np),
           R2 = R - 0.5; // radiousを変える時はここを変える;
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
        a[i] = 0.5;
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
            if (dr2 < cut2) {
                aij = (a[i] + a[list[i][j]]);
                w2 = aij * aij / dr2;
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
    double           vi[2], ri, riw, aij, w2, w6, dUr, fiw[dim], sinco[2];
    constexpr double ddt = 0.0000001, D = usr_sqrt(2. * ddt / 0.01);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[1] = dUr * x[i][1];
            }
        }
        // till here*/
        theta_i[i] += D * gaussian_rand();
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sinco);
        vi[0] = sinco[0] + f[i][0] + fiw[0];
        vi[1] = sinco[1] + f[i][1] + fiw[1];
        x[i][0] += vi[0] * ddt;
        x[i][1] += vi[1] * ddt;
    }
}

void eom_langevin(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                  double *a, int (*list)[Nn]) {
    double zeta = 1.0;
    double fluc = sqrt(8. * zeta * dt);
    double vi[2], ri, riw, aij, w2, w6, dUr, fiw[dim];
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[1] = dUr * x[i][1];
            }
        }
        for (int j = 0; j < dim; j++) {
            v[i][j] += -zeta * v[i][j] * dt + f[i][j] * dt + fiw[j] * dt +
                       fluc * gaussian_rand();
            x[i][j] += v[i][j] * dt;
        }
    }
}

void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double                  ri, riw, aij, w2, w6, dUr, fiw[dim], sico[2];
    constexpr static double D = usr_sqrt(2. * dt / tau), M_inv = dt / M_ss;
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        if (x[i][0] > 0.) {
            ri = sqrt(dist2right(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_rignt);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] < 0.) {
            ri = sqrt(dist2left(x[i]));
            riw = R + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
                w6 = w2 * w2 * w2;
                // w12=w6*w6;
                dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
                fiw[0] = dUr * (x[i][0] - center_left);
                fiw[1] = dUr * x[i][1];
            }
        } else if (x[i][0] == 0.) {
            ri = abs(x[i][1]);
            riw = x0limit + 0.5 - ri;
            if (riw < cut) {
                w2 = 1. / (riw * riw);
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
        v[i][0] += (-v[i][0] + v0 * sico[0] + f[i][0] + fiw[0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sico[1] + f[i][1] + fiw[1]) * M_inv;
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = 0.;
    }
}

void make_v_thetahist(double (*x)[dim], double (*v)[dim], double(*hist),
                      double *hist2, double *lohist) {
    // lohist  と一緒に運用し、outputでv_theta[i]/lo[i];
    // v_thetaとomegaを算出、histがｖhist2がΩ;
    double              v_t, dr;
    const static double bunbo = 1. / (floor(tmax / dt));
    int                 histint;
    for (int i = 0; i < Np; ++i) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / (dr * dr);
        if (dr < R) {
            histint = (int) dr;
            hist[histint] += v_t * bunbo * dr;
            hist2[histint] += v_t * bunbo;
            lohist[histint] += bunbo;
        }
    }
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1f_coorlo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << endl;
    }
    file.close();
}

void output_ani(int k, double (*v)[dim], double (*x)[dim], int l) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1f_animelo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn, l);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << endl;
    }
    file.close();
}
bool out_setup() { // filenameが１２８文字を超えていたらfalseを返す;
    char     filename[128];
    ofstream file;
    int      test =
        snprintf(filename, 128,
                 "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
                 "setupofst_lo%.3f_tau%.3f_m%.3f_t%d.dat",
                 folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn, tmax);
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
    file.close();
    if (test == -1)
        return false;
    else
        return true;
}

void outputhist(double *hist, int counthistv_theta, double *lohist,
                double *hist2) {
    char     filename[128];
    double   v_theta = 0.;
    double   bitthist = 1.;
    int      Nphist = (int) (R + 1.);
    double   rsyou = R - (int) R;
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "v_thetahist_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append

    if (lohist[0] != 0.) {
        file << (rsyou + 1.) * 0.5 << "\t" << (hist[0] / lohist[0]) << endl;

        v_theta += hist[0] * Np_1;
    } else {
        file << (rsyou + 1.) * 0.5 << "\t" << 0 << endl;
    }
    for (int i = 1; i < Nphist; ++i) {
        if (lohist[i] != 0.) {
            file << i + rsyou + 0.5 << "\t" << (hist[i] / lohist[i]) << endl;

            v_theta += hist[i] * Np_1;
        } else {
            file << i + rsyou + 0.5 << "\t" << 0 << endl;
        }
    }
    file.close();
    snprintf(
        filename, 128,
        "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/v_theta_lo%.3f_tau%.3f.dat",
        folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;

    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "omegahist_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append

    if (lohist[0] != 0.) {
        file << (rsyou + 1.) * 0.5 << "\t" << (hist2[0] / lohist[0]) << endl;
    } else {
        file << (rsyou + 1.) * 0.5 << "\t" << 0 << endl;
    }
    for (int i = 1; i < Nphist; ++i) {
        if (lohist[i] != 0.) {
            file << i + rsyou + 0.5 << "\t" << (hist2[i] / lohist[i]) << endl;
        } else {
            file << i + rsyou + 0.5 << "\t" << 0 << endl;
        }
    }
    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "lohist_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    file << (rsyou + 1.) * 0.5 << "\t"
         << (lohist[0] / (4. * (rsyou + 1.) * (rsyou + 1.))) << endl;
    for (int i = 1; i < Nphist; ++i) {
        file << i + 0.5 + rsyou << "\t"
             << (lohist[i] / (8. * (i + 0.5 + rsyou))) << endl;
    }
    file.close();
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
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "xcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "vcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/"
             "msd_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}
inline double usr_max(double a, double b) { return (a > b) ? a : b; }
inline double usr_min(double a, double b) { return (a > b) ? b : a; }
void          cell_list(int (*list)[Nn], double (*x)[dim]) {
    int                     i, j, k, l, m, lm, mm, map_index, km, nx[Np][2];
    static constexpr double thresh2 = (cut + skin) * (cut + skin);
    double                  dx, dy;
    static constexpr double xlen_2 = (2. * R + Rbit * R) / 2.;
    static constexpr int    Mx = (int) (xlen_2 * 2. / (cut + skin));
    static constexpr int My = (int) (2. * R / (cut + skin)); // M<=2R/(cut+skin)
    static constexpr int m2 = Mx * My;
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
                if ((dx * dx + dy * dy) < thresh2) {
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
                      double (*x_update)[dim], int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.8;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max >= skin2) {
        cell_list(list, x);
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
    int           list[Np][Nn];
    int           counthistv_theta = 0, countout = 0;
    int           Nphist = (int) (R + 1.);
    double        hist[Nphist], lohist[Nphist], hist2[Nphist];
    double        tout = msdini, toutcoord = 0;
    long long int j = 0;
    int           k = 0, kcoord = 0;
    set_diameter(a);
    if (!ini_coord_twocircles(x))
        return -1;
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    ini_hist(theta, Np);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    char foldername[128];
    snprintf(foldername, 128, "%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f",
             folder_name, R, lo, M_ss, tau, Rbit, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    snprintf(foldername, 128, "%sR%.1f_coorlo%.2fMs%.3ftau%.3fbit%.3fv0%.1f",
             folder_name, R, lo, M_ss, tau, Rbit, v0);
    const char *fname2 = foldername;
    mkdir(fname2, 0777);
    snprintf(foldername, 128, "%sR%.1f_animelo%.2fMs%.3ftau%.3fbit%.3fv0%.1f",
             folder_name, R, lo, M_ss, tau, Rbit, v0);
    const char *fname3 = foldername;
    mkdir(fname3, 0777);
    if (!out_setup()) {
        std::cout << "file name is too long" << endl;
        return -1;
    }
    output(-1, v, x, -1);
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
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp9(v, x, f, a, list, theta, j);
    }
    j = 0;
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output(-1, v, x, -1);
            std::cout << "hazure in kasanari" << ch << endl;
            // return -1;
        }
    }
    std::cout << "passed kasanari!" << endl;
    int tmaxbefch = R / (dtlg * 5);
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_langevin(v, x, f, a, list);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output(-1, v, x, -1);
            std::cout << "hazure in kakimaze" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed kakimaze!" << endl;
    j = 0;
    tmaxbefch = tmaxlg / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] > 0 && dist2right(x[ch]) > R * R) ||
            (x[ch][0] < 0 && dist2left(x[ch]) > R * R)) {
            output(-1, v, x, -1);
            std::cout << "hazure in owari" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed owari!" << endl;
    int           ituibi = 0, tauch = tau / dt, tanibitch = tbitani / dt;
    long long int tmaxch = tmax / dt, tanimaxch = tmaxani / dt;
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
    output(j, v, x, k);
    ++k;
    while (j < tanimaxch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        if (j >= kanit) {
            output_ani(j, v, x, kani);
            ++kani;
            kanit += tanibitch;
            if (j >= toutcoord) {
                output(j, v, x, k);
                toutcoord += tauch;
                ++k;
            }
        }
        if (j >= tout) {
            calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
            t[kcoord] = j * dt;
            ++kcoord;
            tout *= msdbit;
        }
    }
    while (j < tmaxch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);

        if (j >= toutcoord) {
            output(j, v, x, k);
            // outtuibi(x, toutcoord, v, ituibi);
            toutcoord += tauch;
            ++k;
        }
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
    ofstream file;
    snprintf(filename, 128,
             "./%sR%.1flo%.2fMs%.3ftau%.3fbit%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
             folder_name, R, lo, M_ss, tau, Rbit, v0, lo, mgn);
    file.open(filename, std::ios::app); // append
    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << " " << endl;
    file << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()
         << endl; // 処理に要した時間をミリ秒に変換
    file.close();

    outputhist(hist, counthistv_theta, lohist, hist2);
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
