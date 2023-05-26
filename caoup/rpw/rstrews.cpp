#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define Np     2048 // 4の倍数であること;NP=4*r^2*lo
#define R      80. // 固定;// ,0.1より大きいこと;
#define tmax   32000 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg 1600 // 緩和時間は10たうとする;
#define temp   160. // v0^2=2D/tau,ここではDを入れること;
#define tau    160.
#define mgn    16. // Omega=omega/tau,ここではomegaを入れること;
#define dim    2    // 変えるときはEomを変えること;
#define cut    1.122462048 // 3.
#define skin   1.5
#define dtlg   0.0001
#define dt     0.0001
// #define polydispersity 0.2 コードも変える;
#define folder_name "rewr80"
#define msdbit      1.1
#define msdini      0.01
#define histbit     0.5

using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.
// //粒径の平均値を変えるときはヒストグラムの変え方も変えること:現在は1;
//  v4:lohistをNpで割らなくした;
//  v五:ディレクトリ名でRを先にした;
// stwv2:ティレク取りメイをloに、ヒストグラムの取り方を外からに;
// rststw:壁との相互作用をWCAに、countbound廃止、corrのカウント方法変更;
// void tout_update(double *tout) { *tout *= 1.1; }

void ini_coord_circle(double (*x)[dim]) {
    double R2 = R - 0.5;
    double num_max = sqrt(Np / M_PI);
    double bitween = (R2 - 0.1) / num_max;
    int    n025 = Np * 0.25;
    int    i, j, k = 0;
    x[0][0] = 0.;
    x[0][1] = 0.;
    for (j = 1; j < num_max; ++j) {
        x[k][0] = j * bitween;
        x[k][1] = 0.;
        x[k + n025][0] = -j * bitween;
        x[k + n025][1] = 0.;
        x[k + n025 * 2][0] = 0.;
        x[k + n025 * 2][1] = j * bitween;
        x[k + n025 * 3][0] = 0.;
        x[k + n025 * 3][1] = -j * bitween;
        ++k;
    }
    for (j = 1; j < num_max; ++j) {
        for (i = 1; i < num_max; ++i) {
            if (i * bitween * i * bitween + j * j * bitween * bitween <
                R2 * R2) {
                x[k][0] = i * bitween;
                x[k][1] = j * bitween;
                x[k + n025][0] = -i * bitween;
                x[k + n025][1] = j * bitween;
                x[k + 2 * n025][0] = -i * bitween;
                x[k + 2 * n025][1] = -j * bitween;
                x[k + 3 * n025][0] = i * bitween;
                x[k + 3 * n025][1] = -j * bitween;
                ++k;
                if (k >= Np * 0.25)
                    break;
            } else {
                continue;
            }
        }
        if (k >= Np * 0.25)
            break;
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

void eom_aoup(double (*v)[dim], double (*x)[dim], double temp0,
              double (*F)[dim]) {
    double tauinv = dt / tau;
    double F0[dim], v0[2], fiw[2], fiwf[dim];
    double fluc = sqrt(temp0 / dt), ri, riw, aij, w2, w6, dUr, duf;
    for (int i = 0; i < Np; ++i) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        fiwf[0] = 0.;
        fiwf[1] = 0.;
        // /*force bitween wall;
        ri = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        riw = R + 1. - ri;
        if (riw < cut) {

            w2 = 1. / (riw * riw);
            w6 = w2 * w2 * w2;
            // w12=w6*w6;
            dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
            fiw[0] = dUr * x[i][0];
            fiw[1] = dUr * x[i][1];
            duf = -6. * w6 *tau/ (ri * ri);
            fiwf[0] = duf * F[i][0] * x[i][0] * x[i][0];
            fiwf[1] = duf * x[i][1] * F[i][1] * x[i][1];
        }
        /// till here;*/
        F0[0] = F[i][0];
        F0[1] = F[i][1];
        F[i][0] += (-F0[0] - F0[1] * mgn + fluc * gaussian_rand()+fiwf[0]) * tauinv;
        F[i][1] += (-F0[1] + F0[0] * mgn + fluc * gaussian_rand()+fiwf[1]) * tauinv;

        v[i][0] = F[i][0] + fiw[0];
        v[i][1] = F[i][1] + fiw[1];
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
    double v_t, dr, invhistbit = 1. / histbit, bunbo = 1. / floor(tmax / dt);
    int    histint;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / (dr * dr);
        if (dr <= R && dr != 0) {
            histint = (int) ceil(dr * invhistbit) - 1;
            hist[histint] += v_t * bunbo * dr;
            hist2[histint] += v_t * bunbo;
            lohist[histint] += bunbo;
        }
    }
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char     filename[128];
    double   Mgn = mgn / tau;
    ofstream file;
    sprintf(filename,
            "./%s_coorlo%.2fv0%.1ftau%.3fm%.3f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << endl;
    }
    file.close();
}

void out_setup() {
    char     filename[128];
    double   Mgn = mgn / tau, v0 = temp / tau;
    ofstream file;
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/setupr%fm%f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << "dt=" << dt << endl;
    file << "cut" << cut << endl;
    file << "skin" << skin << endl;
    // file << "Nn" << Nn << endl;
    file << "Np=" << Np << endl;
    file << "tmax=" << tmax << endl;
    file << "tmaxlg=" << tmaxlg << endl;
    file << "temp=" << temp << endl;
    // file << "ens=" << ensemble << endl;
    file << "type=" << 2 << endl;
    file << "2DkaraD" << endl;
    file << "壁はWCA" << endl;
    file << "cell list" << endl;
    file.close();
}

void outputhist(double *hist, int counthistv_theta, double *lohist,
                double *hist2) {
    char          filename[128];
    double        v_theta = 0.;
    double        bitthist = histbit;
    int           Nphist = (int) (R / bitthist + 1.);
    double        Mgn = mgn / tau, v0 = temp / tau;
    double        rsyou = R - (int) R;
    std::ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/v_thetahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append

    if (lohist[0] != 0.) {
        file << (rsyou + bitthist) * 0.5 << "\t" << (hist[0] / lohist[0])
             << std::endl;

        v_theta += hist[0] / Np;
    } else {
        file << (rsyou + bitthist) * 0.5 << "\t" << 0 << std::endl;
    }
    for (int i = 1; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << rsyou + (i + 0.5) * bitthist << "\t"
                 << (hist[i] / lohist[i]) << std::endl;

            v_theta += hist[i] / Np;
        } else {
            file << rsyou + (i + 0.5) * bitthist << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/v_theta_lo%.3f_tau%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << Mgn << "\t" << R << "\t" << v_theta << std::endl;
    file << tau << "\t" << Mgn << "\t" << R << "\t" << v_theta << std::endl;

    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/omegahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append

    if (lohist[0] != 0.) {
        file << (rsyou + bitthist) * 0.5 << "\t" << (hist2[0] / lohist[0])
             << std::endl;

    } else {
        file << (rsyou + bitthist) * 0.5 << "\t" << 0 << std::endl;
    }
    for (int i = 1; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << rsyou + (i + 0.5) * bitthist << "\t"
                 << (hist2[i] / lohist[i]) << std::endl;

        } else {
            file << bitthist * (i + 0.5) << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/lohist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    file << (rsyou + bitthist) * 0.5 << "\t"
         << (lohist[0] / (M_PI * (rsyou + bitthist) * (rsyou + bitthist)))
         << std::endl;
    for (int i = 1; i < Nphist; i++) {
        file << bitthist * (i + 0.5) + rsyou << "\t"
             << (lohist[i] /
                 (4 * bitthist * (2. * (i * bitthist + rsyou) + bitthist)))
             << std::endl;
    }
    file.close();
}

void outtuibi(double (*x)[dim], double t, double (*v)[dim], int cu) {
    char   filename[128];
    double v_theta;
    double Mgn = mgn / tau;
    v_theta = x[cu][0] * v[cu][1] - x[cu][1] * v[cu][0];

    ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/coordtui_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn);
    file.open(filename, std::ios::app); // append

    file << t << "\t" << x[cu][0] << "\t" << x[cu][1] << "\t" << v_theta
         << endl;

    file.close();
}

void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v0)[dim],
               double (*v)[dim], double *xcor, double *vcor, int k,
               double *msd) {
    double bunbo = 1. / Np, dr;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < dim; ++j) {
            xcor[k] += x0[i][j] * x[i][j] * bunbo;
            vcor[k] += v0[i][j] * v[i][j] * bunbo;
            dr = x[i][j] - x0[i][j];
            msd[k] += dr * dr * bunbo;
        }
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) {
    char     filename[128];
    double   v_theta;
    double   Mgn = mgn / tau, v0 = temp / tau;
    ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/msd_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, Np * 0.25 / (R * R),
            tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}

void outv_thetat(double *hist, int counthistv_theta, double tout) {
    char     filename[128];
    double   v_theta = 0.;
    double   Mgn = mgn / tau;
    int      Nphist = (int) (R + 1.);
    ofstream file;
    sprintf(filename,
            "./%slo%.2fv0%.1ftau%.3fm%.3f/v_thetatime_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn);
    file.open(filename, std::ios::app); // append
    double bunbo = floor(tmax / dt) / (counthistv_theta * Np);
    for (int i = 0; i < Nphist; ++i) {

        v_theta += hist[i] * bunbo;
    }
    file << tout << "\t" << v_theta << endl;
    file.close();
}

int main() {
    double x[Np][dim], v[Np][dim], F[Np][dim], x0[Np][dim],
        v0[Np][dim], disp_max = 0.;
    // int(*list)[Nn] = new int[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    int    Nphist = (int) (R/histbit + 1.);
    double hist[Nphist], lohist[Nphist], hist2[Nphist];
    double tout = msdini, toutcoord = 0;

    int j = 0, k = 0, kcoord = 0;

    ini_coord_circle(x);
    ini_array(v);
    ini_array(F);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    double Mgn = mgn / tau;
    char   foldername[128];
    sprintf(foldername, "%slo%.2fv0%.1ftau%.3fm%.3f", folder_name,
            Np * 0.25 / (R * R), temp / tau, tau, Mgn);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%s_coorlo%.2fv0%.1ftau%.3fm%.3f", folder_name,
            Np * 0.25 / (R * R), temp / tau, tau, Mgn);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);

    out_setup();
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
    double ttemp = 5. * temp;
    if (ttemp / tau < 5)
        ttemp = 5. * tau;
    j = 0;
    int tmaxbefch = 10 / dt;
    while (j < tmaxbefch) {
        ++j;
        eom_aoup(v, x, ttemp, F);
    }

    j = 0;
    tmaxbefch = tmaxlg / dt;
    while (j < tmaxbefch) {
        ++j;
        eom_aoup(v, x, temp, F);
    }

    int ituibi = 0, tauch = tau / dt, tmaxch = tmax / dt;

    for (int xnp = 0; xnp < Np; xnp++) {

        for (int xdim = 0; xdim < dim; xdim++) {
            x0[xnp][xdim] = x[xnp][xdim];
            v0[xnp][xdim] = v[xnp][xdim];
        }
    }
    j = 0;
    tout = msdini / dt;
    toutcoord = 0.;
    k = 0;
    kcoord = 0;
    int kani = 0;
    int kanit = 0;

    calc_corr(x, x0, v0, v, msd, vcor, kcoord, msd2);
    t[0] = 0.;
    ++kcoord;
    output(j, v, x, k);
    ++k;
    while (j < tmaxch) {
        ++j;
        eom_aoup(v, x, temp, F);
        make_v_thetahist(x, v, hist, hist2, lohist);

        if (j >= toutcoord) {
            output(j, v, x, k);
            // outtuibi(x, toutcoord, v, ituibi);
            toutcoord += tauch;
            ++k;
        }
        if (j >= tout) {
            calc_corr(x, x0, v0, v, msd, vcor, kcoord, msd2);
            t[kcoord] = j * dt;
            ++kcoord;
            tout *= msdbit;
        }
    }
    int    counthazure = 0, maxnum = 0;
    double ave;

    char filename[128];

    ofstream file;

    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/kekkalo%.3fm%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), Mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << endl;
    file.close();

    outputhist(hist, counthistv_theta, lohist, hist2);
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
