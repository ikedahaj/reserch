#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define Np          16 // 4の倍数であること;NP=4*r^2*lo
#define Nn          10
#define R           80. // 固定;// ,0.1より大きいこと;
#define M           61  // M<=2R/(cut+skin)
#define tmax        16000 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg      800 // 緩和時間は10たうとする;
#define v0          1.
#define tau         80.
#define mgn         0.08 // Omega=omega/tau,ここではomegaを入れること;
#define tmaxani     500  //>tmaxの時プログラムを変更すること;
#define tbitani     1
#define dim         2           // 変えるときはEomを変えること;
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

void calc_force(double (*x)[dim], double (*f)[dim], double *a,
                int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, /*w12,*/ aij;
    ini_array(f);

    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];

            dr2 = dx * dx + dy * dy;
            if (dr2 < cut * cut) {
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
void eom_abp2(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              double vv0, int (*list)[Nn], double *theta_i) {
    double D = sqrt(2. / (tau * dt)), M_PI2 = 2. * M_PI, Co;
    for (int i = 0; i < Np; i++) {
        theta_i[i] += (D * gaussian_rand() + mgn) * dt;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        Co = cos(theta_i[i]);
        v[i][0] = v0 * Co + f[i][0];
        v[i][1] =
            v0 * ((theta_i[i] > 0) - (theta_i[i] < 0)) * sqrt(1. - Co * Co) +
            f[i][1];
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}
void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a, int (*list)[Nn], double *theta_i) {
    double D = sqrt(2. * dt / tau), Mg = mgn * dt, M_PI2 = 2. * M_PI, Co, ri,
           riw, aij, w2, w6, dUr, fiw[dim];
    for (int i = 0; i < Np; i++) {
        fiw[0] = 0.;
        fiw[1] = 0.;
        // /*force bitween wall;
        ri = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        riw = R + 0.5 - ri;
        if (riw < cut) {
            w2 = 1. / (riw * riw);
            w6 = w2 * w2 * w2;
            // w12=w6*w6;
            dUr = (-48. * w6 + 24.) * w6 / (riw * ri);
            fiw[0] = dUr * x[i][0];
            fiw[1] = dUr * x[i][1];
        }
        // till here*/
        theta_i[i] += (D * gaussian_rand() + mgn) * dt;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        v[i][0] = v0 * cos(theta_i[i]) + f[i][0];
        v[i][1] =
            v0 * sin(theta_i[i]) +
            f[i][1];
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
    double v_t, dr, bunbo = 1. / (floor(tmax / dt));
    int    histint;
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
    sprintf(filename,
            "./%s_coorlo%.2ftau%.3fm%.3fv0%.1f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn, l);
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
    sprintf(filename,
            "./%s_animelo%.2ftau%.3fm%.3fv0%.1f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn, l);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << endl;
    }
    file.close();
}
void out_setup() {
    char     filename[128];
    ofstream file;
    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/setupr%fm%f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, R, mgn);
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
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
    file.open(filename /*,std::ios::app*/); // append

    if (lohist[0] != 0.) {
        file << (rsyou + 1.) * 0.5 << "\t" << (hist[0] / lohist[0]) << endl;

        v_theta += hist[0] / Np;
    } else {
        file << (rsyou + 1.) * 0.5 << "\t" << 0 << endl;
    }
    for (int i = 1; i < Nphist; ++i) {
        if (lohist[i] != 0.) {
            file << i + rsyou + 0.5 << "\t" << (hist[i] / lohist[i]) << endl;

            v_theta += hist[i] / Np;
        } else {
            file << i + rsyou + 0.5 << "\t" << 0 << endl;
        }
    }
    file.close();
    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/v_theta_lo%.3f_tau%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau);
    file.open(filename, std::ios::app); // append
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;
    file << tau << "\t" << mgn << "\t" << R << "\t" << v_theta << endl;

    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/omegahist_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
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
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
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
    double bunbo = 1. / Np, dr;
    for (int i = 0; i < Np; ++i) {
        for (int j = 0; j < dim; ++j) {
            xcor[k] += x0[i][j] * x[i][j] * bunbo;
            vcor[k] += v1[i][j] * v[i][j] * bunbo;
            dr = x[i][j] - x0[i][j];
            msd[k] += dr * dr * bunbo;
        }
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) {
    char     filename[128];
    double   v_theta;
    ofstream file;
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/xcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/vcor_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    sprintf(filename,
            "./%slo%.2ftau%.3fm%.3fv0%.1f/msd_lo%.3f_tau%.3f_m%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}

void cell_list(int (*list)[Nn], double (*x)[dim]) {
    int    i, j, k, l, m, lm, mm, map_index, km, nx[Np][2];
    double dx, dy, thresh2 = (cut + skin) * (cut + skin), bit = M / (2. * R);

    int(*map)[Np] = new int[M * M][Np];

    for (i = 0; i < M; ++i)
        for (j = 0; j < M; ++j)
            map[i + M * j][0] = 0;

    for (i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0] + R) * bit);
        nx[i][1] = (int) ((x[i][1] + R) * bit);
        for (m = max(nx[i][1] - 1, 0), mm = min(nx[i][1] + 1, M - 1); m <= mm;
             ++m) {
            for (l = max(nx[i][0] - 1, 0), lm = min(nx[i][0] + 1, M - 1);
                 l <= lm; ++l) {
                map_index = l + M * m;
                map[map_index][map[map_index][0] + 1] = i;
                map[map_index][0]++;
            }
        }
    }

    for (i = 0; i < Np; ++i) {
        list[i][0] = 0;
        // nx = (int)((x[i][0]+R) * bit);
        // ny = (int)((x[i][1]+R) * bit);
        map_index = nx[i][0] + M * nx[i][1];
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
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max > skin * skin * 0.25) {
        cell_list(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        *disp_max = 0.0;
        // count = 0;
    }
}

int main() {
    double x[Np][dim], v[Np][dim], theta[Np], a[Np], f[Np][dim], x0[Np][dim],
        v1[Np][dim], x_update[Np][dim], disp_max = 0.;
    // int(*list)[Nn] = new int[Np][Nn];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    int    Nphist = (int) (R + 1.);
    double hist[Nphist], lohist[Nphist], hist2[Nphist];
    double tout = msdini, toutcoord = 0;

    int j = 0, k = 0, kcoord = 0;
    set_diameter(a);
    ini_coord_circle(x);
    ini_array(v);
    ini_array(f);
    ini_hist(theta, Np);
    ini_hist(hist, Nphist);
    ini_hist(lohist, Nphist);
    ini_hist(hist2, Nphist);
    char foldername[128];
    sprintf(foldername, "%slo%.2ftau%.3fm%.3fv0%.1f", folder_name,
            Np * 0.25 / (R * R), tau, mgn, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    sprintf(foldername2, "%s_coorlo%.2ftau%.3fm%.3fv0%.1f", folder_name,
            Np * 0.25 / (R * R), tau, mgn, v0);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    sprintf(foldername, "%s_animelo%.2ftau%.3fm%.3fv0%.1f", folder_name,
            Np * 0.25 / (R * R), tau, mgn, v0);
    const char *fname3 = foldername;
    mkdir(fname3, 0777);
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
    
    int tmaxbefch = 10 / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp1(v,x,f,a,list,theta);
    }

    j = 0;
    tmaxbefch = tmaxlg / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp2(v, x, f, a, v0, list, theta);
    }

    int ituibi = 0, tauch = tau / dt, tmaxch = tmax / dt,
        tanimaxch = tmaxani / dt, tanibitch = tbitani / dt;

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
        eom_abp1(v, x, f, a,  list, theta);
        make_v_thetahist(x, v, hist, hist2, lohist);
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
        make_v_thetahist(x, v, hist, hist2, lohist);

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
    char filename[128];

    ofstream file;

    sprintf(filename, "./%slo%.2ftau%.3fm%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
            folder_name, Np * 0.25 / (R * R), tau, mgn, v0, Np * 0.25 / (R * R),
            mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << endl;
    file.close();

    outputhist(hist, counthistv_theta, lohist, hist2);
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
