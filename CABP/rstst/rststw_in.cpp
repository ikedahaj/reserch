
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

// #define Np          12800 // 4の倍数であること;NP=4*r^2*lo
#define lo          0.5 // コンパイル時に代入する定数;
#define Nn          50
#define R           10. // 固定;// ,0.1より大きいこと;
#define R_in        2.  // 内側の円の半径;
#define tmax        1600 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg      80// 緩和時間は10たうとする;
#define v0          1.
#define tau         10. // コンパイル時に-D{変数名}={値}　例:-Dtau=80　とすること;
#define mgn         0. // Omega=omega/tau,ここではomegaを入れること;
#define tmaxani     500 //>tmaxの時プログラムを変更すること;
#define tbitani     1
#define dim         2 // 変えるときはEomを変えること;
#define ratf        1.
#define cut         1.122462048 // 3.
#define skin        1.5
#define dtlg        0.0001
#define dt          0.0001
#define folder_name "stwr80"
#define msdbit      1.1
#define msdini      0.01
// #define polydispersity 0.2 コードも変える;
using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.

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
void ini_coord_circle(double (*x)[dim]) {
    int    Np_m = lo * R * R * 16 * M_1_PI;
    double num_max = sqrt(Np_m) + 1;
    double bitween = 2 * (R-0.5) / num_max, R2 = R - 0.5, R1 = R_in + 0.5, poj[2], r;
    int    k = 0;
    for (int j = 0; j < num_max; ++j) {
        for (int i = 0; i < num_max; ++i) {
            poj[0] = i * bitween + 0.5-R;
            poj[1] = j * bitween + 0.5-R;
            r = poj[0] * poj[0] + poj[1] * poj[1];
            if (r <= R2 * R2 && R1 * R1 <= r) {
                x[k][0] = i * bitween;
                x[k][1] = j * bitween;
                ++k;
            }else{
                continue;
            }
            if (k >= Np)
                break;
        }
        if (k >= Np)
            break;
    }
    if(k!=Np)std::cout<<k;
    else std::cout<<k;
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
    dr = r2 - R_in - 0.5;
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
    double ddt = 0.0000001, fiw[dim], sinco[2], fluc = sqrt(6 * ddt);
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        for (int j = 0; j < dim; j++) {
            v[i][j] += (-v[i][j] + f[i][j]) * ddt + fluc * gaussian_rand();
            x[i][j] += v[i][j] * ddt;
        }
    }
}
void eom_abp8(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i) {
    double        fiw[dim], sico[2];
    static double D = sqrt(2. * dtlg / tau);
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
        v[i][0] = sico[0] + f[i][0] + fiw[0];
        v[i][1] = sico[1] + f[i][1] + fiw[1];
        x[i][0] += v[i][0] * dtlg;
        x[i][1] += v[i][1] * dtlg;
    }
}

void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim],
              int (*list)[Nn], double *theta_i) {
    double              ri, riw, aij, w2, w6, dUr, fiw[dim], sico[2];
    static const double D = sqrt(2. * dt / tau);
    calc_force(x, f, list);
    for (int i = 0; i < Np; i++) {
        calc_force_2w(x[i], fiw);
        theta_i[i] += D * gaussian_rand() + Mg;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        v[i][0] = v0 * sin(theta_i[i]) + f[i][0] + fiw[0];
        v[i][1] = v0 * cos(theta_i[i]) + f[i][1] + fiw[1];
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
    const static double bunbo = 1. / tau;
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

void output(double (*v)[dim], double (*x)[dim]) {
    static int l = 1;
    char       filename[128];
    ofstream   file;
    sprintf(filename,
            "./%s_coorlo%.2ftau%.3fm%.3fv0%.1f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}

void output_ani(double (*v)[dim], double (*x)[dim]) {
    static int l = 0;
    char       filename[128];
    ofstream   file;
    sprintf(filename,
            "./%s_animelo%.2ftau%.3fm%.3fv0%.1f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn, l);
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
    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/setupr%fm%f.dat",
            folder_name, lo, tau, mgn, v0, R, mgn);
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

void outputhist(double *hist, int counthistv_theta, double *lohist,
                double *hist2) {
    char     filename[128];
    double   v_theta = 0.;
    double   bitthist = 1.;
    int      Nphist = (int) (R + 1.);
    double   rsyou = R - (int) R;
    ofstream file;
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/v_thetahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
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
    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/v_theta_lo%.3f_tau%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;

    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/omegahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
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
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/lohist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
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
            "./%slo%.2ftau%.3fm%.3fv0%.1f/xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << xcor << endl;
    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << vcor << endl;
    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/msd_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    file << j * dt << "\t" << msd << endl;
    file.close();
}

void cell_list(int (*list)[Nn], double (*x)[dim]) {
    int              map_index, km, nx[Np][2];
    constexpr int M=2*R/(cut+skin);
    constexpr double thresh2 = (cut + skin) * (cut + skin), bit = M / (2. * R);
    double           dx, dy;
    constexpr int    m2 = M * M;
    int(*map)[Np] = new int[m2][Np];

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            map[i + M * j][0] = 0;

    for (int i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0] + R) * bit);
        nx[i][1] = (int) ((x[i][1] + R) * bit);
        for (int m = max(nx[i][1] - 1, 0), mm = min(nx[i][1] + 1, M - 1);
             m <= mm; ++m) {
            for (int l = max(nx[i][0] - 1, 0), lm = min(nx[i][0] + 1, M - 1);
                 l <= lm; ++l) {
                map_index = l + M * m;
                map[map_index][map[map_index][0] + 1] = i;
                map[map_index][0]++;
            }
        }
    }

    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + M * nx[i][1];
        for (int k = 1, km = (map[map_index][0]); k <= km; ++k) {
            if (map[map_index][k] > i) {
                dx = x[i][0] - x[map[map_index][k]][0];
                dy = x[i][1] - x[map[map_index][k]][1];
                if ((dx * dx + dy * dy) < thresh2) {
                    list[i][0]++;
                    list[i][list[i][0]] = map[map_index][k];
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

void auto_list_update(double (*x)[dim], double (*x_update)[dim],
                      int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    static double    disp_max = skin2 + 100;
    calc_disp_max(&(disp_max), x, x_update);
    if (disp_max > skin2) {
        cell_list(list, x);
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
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    int    Nphist = (int) (R + 1.);
    double hist[Nphist], lohist[Nphist], hist2[Nphist];
    ini_coord_circle(x);
    ini_array(v);
    ini_array(f);
    ini_hist(theta, Np);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    char foldername[128];
    sprintf(foldername, "%slo%.2ftau%.3fm%.3fv0%.1f", folder_name, lo, tau, mgn,
            v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%s_coorlo%.2ftau%.3fm%.3fv0%.1f", folder_name, lo,
            tau, mgn, v0);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    sprintf(foldername, "%s_animelo%.2ftau%.3fm%.3fv0%.1f", folder_name, lo,
            tau, mgn, v0);
    const char *fname3 = foldername;
    mkdir(fname3, 0777);

    out_setup();
    std::cout << foldername << endl;

    for (int j = 0; j < 1e9; j++) {
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
    int tmaxbefch = R / (dtlg * 5);
    for (int j = 0; j < tmaxbefch; j++) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp8(v, x, f, list, theta);
    }
    for (int ch = 0; ch < Np; ch++) {
        if ((x[ch][0] * x[ch][0] + x[ch][1] * x[ch][1]) > (R * R)) {
            output(v, x);
            std::cout << "hazure in kakimaze" << ch << endl;
            return -1;
        }
    }
    std::cout << "passed kakimaze!" << endl;

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
        tanimaxch = tmaxani / dt, tanibitch = tbitani / dt;

    for (int xnp = 0; xnp < Np; xnp++) {
        for (int xdim = 0; xdim < dim; xdim++) {
            x0[xnp][xdim] = x[xnp][xdim];
            v1[xnp][xdim] = v[xnp][xdim];
        }
    }
    unsigned long long j = 0;
    double             tout = msdini / dt, toutcoord = 0;
    long long int      kanit = 0;
    output(v, x);
    while (j < tanimaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, list, theta);
        if (j >= kanit) {
            output_ani(v, x);
            kanit += tanibitch;
            if (j >= toutcoord) {
                output(v, x);
                toutcoord += tauch;
                make_v_thetahist(x, v, hist, hist2, lohist);
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
        if (j >= toutcoord) {
            output(v, x);
            make_v_thetahist(x, v, hist, hist2, lohist);
            toutcoord += tauch;
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

    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
            folder_name, lo, tau, mgn, v0, lo, mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << " " << endl;
    file << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()
         << endl; // 処理に要した時間をミリ秒に変換
    file.close();
    snprintf(filename, 128, "%s_animelo%.2ftau%.3fm%.3fv0%.1f/tyokkei.dat",
             folder_name, lo, tau, mgn, v0);
    file.open(filename);
    file << tbitani << endl;
    for (int i = 0; i < Np; i++) {
        file << 1 << endl;
    }

    outputhist(hist, counthistv_theta, lohist, hist2);
    std::cout << "done" << endl;
    return 0;
}
