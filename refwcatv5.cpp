#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define Np     1024 // 8192// 4096// 4の倍数であること;
#define Nn     400
#define lowall 1
#define R                                                                      \
    22.62741699796952 // 32.//25.29822128//64.//128.// lo=0.5,N1000//3.
                      // //R=sqrt(np/4/lo);
// ,0.1より大きいこと;
#define tmax     500 // 973.686//2*100たうとするただし2000まで;
#define tmaxlg   100 // 緩和時間は10たうとする;
#define dtlg     0.0001
#define dt       0.0001
#define temp     1.          // v0^2=2D/tau,ここではDを入れること;
#define dim      2           // 変えるときはEomを変えること;
#define cut      1.122462048 // 3.
#define skin     1.5
#define tau      1.
#define ensemble 1
// #define polydispersity 0.2 コードも変える;
#define folder_name "reawt2"
#define mgn         1. // Omega=omega/tau,ここではomegaを入れること;
#define K           0.
#define Nphist      22 // 値は結果と相談すること;
// #define kiru        0.   // 0ならwall_boundaryが止まる;現在機能なし;
// v4:lohistをNpで割らなくした;
// v五:ディレクトリ名でRを先にした;
void tout_update(double *tout) { *tout *= 1.1; }

void ini_coord_circle(double (*x)[dim]) {
    double R2 = R - 0.5;
    double num_max = sqrt(Np / M_PI);
    double bitween = (R2 - 0.1) / num_max;
    int    n025 = Np * 0.25;
    int    i, j, k = 0;
    x[0][0] = 0.;
    x[0][1] = 0.;
    for (j = 1; j < num_max; j++) {
        x[k][0] = j * bitween;
        x[k][1] = 0.;
        x[k + n025][0] = -j * bitween;
        x[k + n025][1] = 0.;
        x[k + n025 * 2][0] = 0.;
        x[k + n025 * 2][1] = j * bitween;
        x[k + n025 * 3][0] = 0.;
        x[k + n025 * 3][1] = -j * bitween;
        k++;
    }
    for (j = 1; j < num_max; j++) {
        for (i = 1; i < num_max; i++) {
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
                k++;
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
void setcoord_cirwall(double (*x)[dim]) {
    double theta = 1. / (R * lowall);
    double Nwall = 2. * M_PI * R * lowall;
    for (int i = Np; i < Nwall; i++) {
        x[i][0] = R * cos(theta * i);
        x[i][1] = R * sin(theta * i);
    }
}
void set_diameter(double *a) {
    for (int i = 0; i < Np; i++)
        a[i] = 0.5;
}

void wall_boundary(double (*x)[dim], double (*v)[dim], double *a,
                   double *count) {
    double rpar2;
    double v0[dim];
    double bunbo = 1. / (tmax * ensemble), rwa, rr;
    for (int i = 0; i < Np; i++) {
        rpar2 = x[i][0] * x[i][0] + x[i][1] * x[i][1];
        rwa = R - a[i];
        if (rpar2 >= rwa * rwa) {
            count[i] += bunbo;
            rr = rwa / sqrt(rpar2);
            v0[0] = v[i][0];
            v0[1] = v[i][1];
            v[i][0] = ((-x[i][0] * x[i][0] + x[i][1] * x[i][1]) * v0[0] -
                       2. * x[i][0] * x[i][1] * v0[1]) /
                      rpar2;
            v[i][1] = (-2. * x[i][0] * x[i][1] * v0[0] +
                       (x[i][0] * x[i][0] - x[i][1] * x[i][1]) * v0[1]) /
                      rpar2;
            x[i][0] *= rr;
            x[i][1] *= rr;
        }
    }
}

void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x[i][j] = 0.0;
}

void list_verlet(int (*list)[Nn], double (*x)[dim]) {
    double dx, dy, dr2;
    double thresh = cut + skin;
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Nn; j++)
            list[i][j] = 0;

    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Np; j++) {
            if (j > i) {
                dx = x[i][0] - x[j][0];
                dy = x[i][1] - x[j][1];
                //   dx-=L*floor((dx+0.5*L)/L);
                //   dy-=L*floor((dy+0.5*L)/L);
                dr2 = dx * dx + dy * dy;
                if (dr2 < thresh * thresh) {
                    list[i][0]++;
                    list[i][(int) list[i][0]] = j;
                }
            }
        }
}

void calc_force(double (*x)[dim], double (*f)[dim], double *a,
                int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, w12, aij;
    ini_array(f);

    for (int i = 0; i < Np; i++)
        for (int j = 1; j <= list[i][0]; j++) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];

            dr2 = dx * dx + dy * dy;
            if (dr2 < cut * cut) {
                aij = (a[i] + a[list[i][j]]);
                w2 = aij * aij / dr2;
                w6 = w2 * w2 * w2;
                w12 = w6 * w6;
                dUr = (-48. * w12 + 24. * w6) / dr2 /* -12. * w12 / dr2*/;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}

void calc_force_harmonic(double (*x)[dim], double (*f)[dim]) {
    for (int i = 0; i < Np; i++) {
        f[i][0] = -K * x[i][0];
        f[i][1] = -K * x[i][1];
    }
}

void eom_aoup(
    double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
    double temp0, int (*list)[Nn],
    double (
        *F)[dim]) { // wallはここに入っていない。使うときに合わせて用いること;
    double tauinv = dt / tau;
    double F0[dim];
    double fluc = sqrt(2. * temp0 / dt);
    for (int i = 0; i < Np; i++) {

        F0[0] = F[i][0];
        F0[1] = F[i][1];
        F[i][0] += (-F0[0] - F0[1] * mgn + fluc * gaussian_rand()) * tauinv;
        F[i][1] += (-F0[1] + F0[0] * mgn + fluc * gaussian_rand()) * tauinv;
    }

    calc_force(x, f, a, list);
    // calc_force_harmonic(x, f);
    for (int i = 0; i < Np; i++) {
        v[i][0] = F[i][0] + f[i][0];
        v[i][1] = F[i][1] + f[i][1];
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; i++) {
        hist[i] = 0.;
    }
}

void make_v_thetahist(double (*x)[dim], double (*v)[dim], double(*hist),
                      double *hist2) {
    // lohist  と一緒に運用し、outputでv_theta[i]/lo[i];
    // v_thetaとomegaを算出、histがｖhist2がΩ;
    double v_t, dr, Rhistinv = Nphist / R,
                    bunbo = Nphist / (R * ensemble * floor(tmax / dt));
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / dr;
        if (dr < R) {
            hist[(int) floor(dr * Rhistinv)] += v_t * bunbo;
            hist2[(int) floor(dr * Rhistinv)] += v_t * bunbo / dr;
        }
    }
}

void make_lohist(double (*x)[dim],
                 double(*hist)) { // 関数を使うときは使用回数で割ること
    double dr, Rhistinv = Nphist / R,
               bunbo = Nphist / (R * ensemble * floor(tmax / dt));
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dr < R)
            hist[(int) floor(dr * Rhistinv)] += bunbo;
    }
}
void make_vt1hist(double (*x)[dim], double (*v)[dim], double(*hist)) {
    // lohist  と一緒に運用し、outputでv_theta[i]/lo[i];
    double v_t, dr, Rhistinv = Nphist / R, bunbo = Nphist / R;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        v_t = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / dr;
        if (dr < R)
            hist[(int) floor(dr * Rhistinv)] += v_t * bunbo;
    }
}

void make_lo1hist(double (*x)[dim], double(*hist)) {
    double dr, Rhistinv = Nphist / R, bunbo = Nphist / (R);
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dr < R)
            hist[(int) floor(dr * Rhistinv)] += bunbo;
    }
}
void make_lohistx0(double (*x)[dim], double(*hist),
                   double (*x0)[dim]) { // 関数を使うときは使用回数で割ること
    double dr, Rhistinv = Nphist / R;
    for (int i = 0; i < Np; i++) {
        dr = sqrt((x[i][0] - x0[i][0]) * (x[i][0] - x0[i][0]) +
                  (x[i][1] - x0[i][1]) * (x[i][1] - x0[i][1]));
        if (dr <= R)
            hist[(int) floor(dr * Rhistinv)] += 1.;
    }
}
void calc_diff(double (*x)[dim], double (*v)[dim], double *v_thetahist,
               double *lohist, double *diffv_theta, double *difflo) {
    double dr, v_theta, Rhistinv = Nphist / R, bub = 1. / Np, dif;
    for (int i = 0; i < Np; i++) {
        dr = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (dr < R) {
            v_theta = (x[i][0] * v[i][1] - x[i][1] * v[i][0]) / dr;
            dif = (dr - lohist[(int) floor(dr * Rhistinv)]);
            difflo[(int) floor(dr * Rhistinv)] += dif * dif / Np;
            dif = v_theta - v_thetahist[(int) floor(dr * Rhistinv)] /
                                lohist[(int) floor(dr * Rhistinv)];
            diffv_theta[(int) floor(dr * Rhistinv)] +=
                dif * dif / lohist[(int) floor(dr * Rhistinv)];
        }
    }
}
void out_gosahist(double (*x)[dim], double (*v)[dim]) {
    double v_thetahist[Nphist], lohist[Nphist], diffv_theta[Nphist],
        difflo[Nphist];
    double Mgn = mgn / tau, v0 = temp / tau;

    ini_hist(lohist, Nphist);
    ini_hist(v_thetahist, Nphist);
    ini_hist(difflo, Nphist);
    ini_hist(diffv_theta, Nphist);
    make_lo1hist(x, lohist);
    make_vt1hist(x, v, v_thetahist);
    calc_diff(x, v, v_thetahist, lohist, diffv_theta, difflo);
    char          filename[128];
    double        bitthist = R / Nphist;
    std::ofstream file;
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/lzdif_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << (i + 0.5) * bitthist << "\t" << (v_thetahist[i] / lohist[i])
                 << "\t" << diffv_theta[i] << std::endl;

        } else {
            file << i * bitthist << "\t" << 0 << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/lodif_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename, std::ios::app);
    for (int i = 0; i < Nphist; i++) {
        file << (i + 0.5) * bitthist << "\t"
             << lohist[i] / (2. * M_PI * (i + 0.5) * bitthist) << "\t"
             << difflo[i] << std::endl;
    }
    file.close();
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char   filename[128];
    double dr;
    double v_theta[Np];
    for (int i = 0; i < Np; i++) {
        // dr=sqrt(x[i][0]*x[i][0]+x[i][1]*x[i][1]);
        v_theta[i] = (x[i][0] * v[i][1] - x[i][1] * v[i][0]);
    }
    double        Mgn = mgn / tau;
    std::ofstream file;
    sprintf(filename,
            "./%s_coorR%.1fv0%.1ftau%.3fm%.3f/"
            "tyouwaenn_R%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, R, temp / tau, tau, Mgn, R, tau, Mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; i++) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << "\t" << v_theta[i] << std::endl;
    }
    file.close();
}
void out_setup() {
    char          filename[128];
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/setupr%fm%f.dat",
            folder_name, R, v0, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << "dt=" << dt << std::endl;
    file << "cut" << cut << std::endl;
    file << "skin" << skin << std::endl;
    file << "Nn" << Nn << std::endl;
    file << "Np=" << Np << std::endl;
    file << "tmax=" << tmax << std::endl;
    file << "tmaxlg=" << tmaxlg << std::endl;
    file << "temp=" << temp << std::endl;
    file << "ens=" << ensemble << std::endl;
    file << "nphist=" << Nphist << std::endl;
    file << "type=" << 2 << std::endl;
    file.close();
}

void outputhist(double *hist, int counthistv_theta, double *lohist,
                double *hist2) {
    char          filename[128];
    double        v_theta = 0.;
    double        bitthist = R / Nphist;
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/v_thetahist_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << i * bitthist << "\t" << (hist[i] / lohist[i]) << std::endl;

            v_theta += hist[i] * bitthist / Np;
        } else {
            file << i * bitthist << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/v_theta_R%.3f_tau%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << Mgn << "\t" << K << "\t" << R << "\t" << v_theta
         << std::endl;
    file << tau << "\t" << Mgn << "\t" << K << "\t" << R << "\t" << v_theta
         << std::endl;

    file.close();
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/omegahist_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < Nphist; i++) {
        if (lohist[i] != 0.) {
            file << i * bitthist << "\t" << (hist2[i] / lohist[i]) << std::endl;

        } else {
            file << i * bitthist << "\t" << 0 << std::endl;
        }
    }
    file.close();
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/lohist_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    double bunbo = bitthist * 8.;
    for (int i = 0; i < Nphist; i++) {
        file << (i + 0.5) * bitthist << "\t"
             << (lohist[i] / (bunbo * (i + 0.5))) << std::endl;
    }
    file.close();
}

void outtuibi(double (*x)[dim], double t, double (*v)[dim], int cu) {
    char   filename[128];
    double v_theta;
    double Mgn = mgn / tau;
    v_theta = x[cu][0] * v[cu][1] - x[cu][1] * v[cu][0];

    std::ofstream file;
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/coordtui_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, temp / tau, tau, Mgn, R, tau, Mgn);
    file.open(filename, std::ios::app); // append

    file << t << "\t" << x[cu][0] << "\t" << x[cu][1] << "\t" << v_theta
         << std::endl;

    file.close();
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
    static int count = 0;
    count++;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max > skin * skin * 0.25) {
        list_verlet(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        *disp_max = 0.0;
        count = 0;
    }
}

void calc_msd(double (*x)[dim], double (*x0)[dim], double *msd, int k) {
    double dr, bunbo = 1. / (Np * ensemble);
    for (int i = 0; i < Np; i++) {
        for (int j = 0; j < dim; j++) {
            dr = x[i][j] - x0[i][j];
            msd[k] += dr * dr * bunbo;
        }
    }
}
void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v0)[dim],
               double (*v)[dim], double *xcor, double *vcor, int k) {
    double bunbo = 1. / (Np * ensemble);
    for (int i = 0; i < Np; i++) {
        for (int j = 0; j < dim; j++) {
            xcor[k] += x0[i][j] * x[i][j] * bunbo;
            vcor[k] += v0[i][j] * v[i][j] * bunbo;
        }
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) { // takusannkaeru
    char          filename[128];
    double        v_theta;
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/xcor_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << msd[i] << std::endl;
    }
    file.close();
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/vcor_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << vcor[i] << std::endl;
    }
    file.close();
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/msd_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, v0, tau, Mgn, R, tau, Mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; i++) {
        file << t[i] << "\t" << msd2[i] << std::endl;
    }
    file.close();
}

void outv_thetat(double *hist, int counthistv_theta, double tout) {
    char          filename[128];
    double        v_theta = 0.;
    double        Mgn = mgn / tau;
    std::ofstream file;
    sprintf(filename,
            "./%sR%.1fv0%.1ftau%.3fm%.3f/v_thetatime_R%.3f_tau%.3f_m%.3f.dat",
            folder_name, R, temp / tau, tau, Mgn, R, tau, Mgn);
    file.open(filename, std::ios::app); // append
    double bunbo = (ensemble * floor(tmax / dt)) / (counthistv_theta * Np);
    for (int i = 0; i < Nphist; i++) {

        v_theta += hist[i] * bunbo;
    }
    file << tout << "\t" << v_theta << std::endl;
    file.close();
}
void ini_test(double (*x)[dim], double (*v)[dim]) {
    for (int i = 0; i < Np; i++) {
        x[i][0] = -1.;
        x[i][1] = 0.;
        v[i][0] = 0.;
        v[i][1] = 2.;
    }
}

void calc_bounddiff(double *bound, double *answer) {
    double dif;
    for (int i = 0; i < Np; i++)
        answer[0] += bound[i] / Np;
    for (int i = 0; i < Np; i++) {
        dif = answer[0] - bound[i];
        answer[1] += dif * dif / Np;
    }
}

int main() {
    double x[Np][dim], x_update[Np][dim], v[Np][dim], f[Np][dim], a[Np],
        F[Np][dim], x0[Np][dim];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    double hist[Nphist], lohist[Nphist], hist2[Nphist];
    double tout = 0.01, toutcoord = 0, U, disp_max = 0.0, countbound[Np];

    int j = 0, k = 0, kcoord = 0;
    set_diameter(a);
    ini_coord_circle(x);
    ini_array(v);
    ini_array(F);
    ini_array(f);
    ini_hist(countbound, Np);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    // ini_test(x,v);
    double Mgn = mgn / tau;
    char   foldername[128];
    sprintf(foldername, "%sR%.1fv0%.1ftau%.3fm%.3f", folder_name, R, temp / tau,
            tau, Mgn);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%s_coorR%.1fv0%.1ftau%.3fm%.3f", folder_name, ,
            R temp / tau, tau, Mgn);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    out_setup();
    while (j * dt < tmax) {
        j++;

        if (j * dt >= tout) {

            tout_update(&tout);

            countout++;
            // std::cout<<k;
        }
    }
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
    while (j * dtlg < 10) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, ttemp, list, F);
        wall_boundary(x, v, a, countbound);
    }
    double tcf[dim];

    // std::cout<<x[0][0];
    j = 0;
    double v0[Np][dim];
    for (int i2 = 0; i2 < Np; i2++)
        for (int se = 0; se < dim; se++)
            v0[i2][se] = v[i2][se];

    // double tcf[dim];
    tcf[0] = 0.;
    tcf[1] = 0.;
    for (int i2 = 0; i2 < Np; i2++) {
        tcf[0] += v0[i2][0] * v[i2][0] / Np;
        tcf[1] += v0[i2][0] * v[i2][0] / Np;
    }
    // std::cout << "tcfbef " << tcf[0] << " " << tcf[1] << std::endl;

    j = 0;
    while (j * dt < tmaxlg) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, temp, list, F);
        wall_boundary(x, v, a, countbound);
    }

    tcf[0] = 0.;
    tcf[1] = 0.;
    for (int i2 = 0; i2 < Np; i2++) {
        tcf[0] += v0[i2][0] * v[i2][0] / Np;
        tcf[1] += v0[i2][0] * v[i2][0] / Np;
    }
    // std::cout << "tcfaft " << tcf[0] << " " << tcf[1] << std::endl;

    // double tcf[dim];
    tcf[0] = 0.;
    tcf[1] = 0.;
    for (int i2 = 0; i2 < Np; i2++) {
        tcf[0] += v[i2][0] * v[i2][0] / Np;
        tcf[1] += floor((x[i2][0] * x[i2][0] + x[i2][1] * x[i2][1]) / (R * R));
    }
    // std::cout << "v^2 " << tcf[0] << " " << tcf[1] << std::endl;

    ini_hist(countbound, Np);
    double r2test, r2testmax = 0.;
    int    ituibi = 0;
    for (int ko = 0; ko < Np; ko++) {
        r2test = x[ko][0] * x[ko][0] + x[ko][1] * x[ko][1];
        if (r2test > r2testmax)
            r2testmax = r2test;
        ituibi = ko;
    }
    for (int i = 0; i < ensemble; i++) {
        for (int xnp = 0; xnp < Np; xnp++) {

            for (int xdim = 0; xdim < dim; xdim++) {
                x0[xnp][xdim] = x[xnp][xdim];
                v0[xnp][xdim] = v[xnp][xdim];
            }
        }
        j = 0;
        tout = 0.01;
        toutcoord = 0.;
        k = 0;
        kcoord = 0;

        calc_corr(x, x0, v0, v, msd, vcor, kcoord);
        t[0] = 0.;
        kcoord++;

        out_gosahist(x, v);
        output(j, v, x, k);
        // outlot(x,k,Nphist,v);
        k++;
        while (j * dt < tmax) {
            j++;
            auto_list_update(&disp_max, x, x_update, list);
            eom_aoup(v, x, f, a, temp, list, F);
            wall_boundary(x, v, a, countbound);
            make_v_thetahist(x, v, hist, hist2);
            make_lohist(x, lohist);
            counthistv_theta++;

            if (j * dt >= toutcoord) {
                output(j, v, x, k);
                outtuibi(x, toutcoord, v, ituibi);
                toutcoord += 0.1;

                // outv_thetat(hist, counthistv_theta, Nphist, toutcoord);

                //	outlot(x,k,Nphist,v);

                k++;
            }
            if (j * dt >= tout) {
                calc_corr(x, x0, v0, v, msd, vcor, kcoord);
                calc_msd(x, x0, msd2, kcoord);
                t[kcoord] = j * dt;
                kcoord++;
                tout_update(&tout);
                // std::cout<<k;
            }
        }
    }

    int counthazure = 0;
    for (int i = 0; i < Np; i++) {
        if (x[i][0] * x[i][0] + x[i][1] * x[i][1] > R * R)
            counthazure++;
    }
    char filename[128];

    std::ofstream file;
    double        sum = 0.;
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/boundR%.3fm%.3f.dat",
            folder_name, R, temp / tau, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; i++) {
        file << countbound[i] << std::endl;
        sum += countbound[i] / Np;
    }
    file.close();
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/kekkaR%.3fm%.3f.dat",
            folder_name, R, temp / tau, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << sum << " " << counthistv_theta << " " << counthazure << std::endl;
    file.close();
    double answer[2];
    answer[0] = 0.;
    answer[1] = 0.;
    calc_bounddiff(countbound, answer);
    sprintf(filename, "./%sR%.1fv0%.1ftau%.3fm%.3f/boundcuntR%.3fm%.3f.dat",
            folder_name, R, temp / tau, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << answer[0] << " " << answer[1] << std::endl;

    file.close();
    outputhist(hist, counthistv_theta, lohist, hist2);
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << sum << std::endl;
    return 0;
}
