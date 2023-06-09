#include <algorithm>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ofstream;

#define Np       5120 // 8192// 4096// 4の倍数であること;
#define Nn       2000
#define R        80. //  //R=sqrt(np/4/lo);0.1より大きいこと;
#define M        61  // <=2R/(cut+skin);
#define tmax     500 // 適当;
#define tmaxlg   500 // 緩和時間は10たうとする;
#define tbit     1.  // ファイルを出すときの間隔;
#define dtlg     0.0001
#define dt       0.0001
#define temp     50.         // v0^2=2D/tau,ここではDを入れること;
#define dim      2           // 変えるときはEomを変えること;
#define cut      1.122462048 // 3.
#define skin     1.5
#define tau      50.
#define ensemble 1
// #define polydispersity 0.2 コードも変える;
#define folder_name "stwr80"
#define mgn         0.5 // Omega=omega/tau,ここではomegaを入れること;
// #define radios 1.
// //粒径の平均値を変えるときはヒストグラムの変え方も変えること:現在は1;
//  v4:lohistをNpで割らなくした;
//  v五:ディレクトリ名でRを先にした;
// stwv2:ティレク取りメイをloに、ヒストグラムの取り方を外からに;
void tout_update(double *tout) { *tout *= 1.1; }

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
                dUr = (-48. * w6 + 24.) * w6 / dr2;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}

void eom_aoup(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              double temp0, int (*list)[Nn], double (*F)[dim]) {
    double tauinv = dt / tau;
    double F0[dim], v0[2], fiw[Np][dim];
    double fluc = sqrt(temp0 / dt);
    calc_force(x, f, a, list);
    double ri, riw, w2, w6, aij, dUr;
    for (int i = 0; i < Np; ++i) {
        fiw[i][0] = 0.;
        fiw[i][1] = 0.;
        ////*force bitween wall;

        ri = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        riw = R + 0.5 - ri;
        if (riw < cut) {
            aij = 0.5 + a[i];
            w2 = aij * aij / (riw * riw);
            w6 = w2 * w2 * w2;
            dUr = (-48. * w6 + 24.) * w6 / (ri * riw);
            fiw[i][0] = dUr * x[i][0];
            fiw[i][1] = dUr * x[i][1];
        }
        /// till here;*/
        F0[0] = F[i][0];
        F0[1] = F[i][1];
        F[i][0] += (-F0[0] - F0[1] * mgn + fluc * gaussian_rand()) * tauinv;
        F[i][1] += (-F0[1] + F0[0] * mgn + fluc * gaussian_rand()) * tauinv;

        v[i][0] = F[i][0] + f[i][0] + fiw[i][0];
        v[i][1] = F[i][1] + f[i][1] + fiw[i][1];
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = 0.;
    }
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char     filename[128];
    double   Mgn = mgn / tau;
    ofstream file;
    sprintf(filename,
            "./%s_animelo%.2fv0%.1ftau%.3fm%.3f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn, l);
    file.open(filename /* std::ios::app*/); // append
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
    sprintf(filename, "./%s_animelo%.2fv0%.1ftau%.3fm%.3f/setupr%fm%f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, R, Mgn);
    file.open(filename, std::ios::app); // append

    file << "dt=" << dt << endl;
    file << "cut" << cut << endl;
    file << "skin" << skin << endl;
    file << "Nn" << Nn << endl;
    file << "Np=" << Np << endl;
    file << "tmax=" << tmax << endl;
    file << "tmaxlg=" << tmaxlg << endl;
    file << "temp=" << temp << endl;
    file << "ens=" << ensemble << endl;
    file << "type=" << 2 << endl;
    file << "2DkaraD" << endl;
    file << "M=" << M << endl;
    file << "cellbit=" << 2. * R / M << endl;
    file.close();
}
void cell_list(int (*list)[Nn], double (*x)[dim]) {
    int    i, j, k;
    int    nx[Np][2];
    int    l, m, lm, mm, map_index, km;
    double dx, dy;
    double thresh2 = (cut + skin) * (cut + skin), bit = M / (2. * R);

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

double calc_disp_max( double (*x)[dim],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp,disp_max=0.;
    for (int i = 0; i < Np; i++) {
        dx = x[i][0] - x_update[i][0];
        dy = x[i][1] - x_update[i][1];

        disp = dx * dx + dy * dy;
        if (disp > disp_max)
            disp_max = disp;
    }
    return disp_max;
}

void auto_list_update( double (*x)[dim],
                      double (*x_update)[dim], int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    
    if (calc_disp_max( x, x_update) > skin * skin * 0.25) {
        cell_list(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        // count = 0;
    }
}

int main() {
    double x[Np][dim], v[Np][dim], f[Np][dim], a[Np], F[Np][dim], x0[Np][dim],
        v0[Np][dim], x_update[Np][dim];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    double tout = 0.01, toutcoord = 0, U;

    int j = 0, k = 0, kcoord = 0;
    set_diameter(a);
    ini_coord_circle(x);
    ini_array(v);
    ini_array(F);
    ini_array(f);

    double Mgn = mgn / tau;
    char   foldername[128];
    sprintf(foldername, "%s_animelo%.2fv0%.1ftau%.3fm%.3f", folder_name,
            Np * 0.25 / (R * R), temp / tau, tau, Mgn);
    const char *fname = foldername;
    mkdir(fname, 0777);
    out_setup();
    cout << foldername << endl;
    j = 0;
    double ttemp = 5. * temp;
    if (ttemp / tau < 5)
        ttemp = 5. * tau;
    j = 0;
    double tmaxlgch = 10 / dt;
    while (j < tmaxlgch) {
        ++j;
        auto_list_update( x, x_update, list);
        eom_aoup(v, x, f, a, ttemp, list, F);
    }

    j = 0;
    tmaxlgch = tmaxlg / dt;
    while (j < tmaxlgch) {
        ++j;
        auto_list_update( x, x_update, list);
        eom_aoup(v, x, f, a, temp, list, F);
    }

    double tmaxch = tmax / dt, tbitch = tbit / dt;

    for (int i = 0; i < ensemble; ++i) {

        j = 0;
        toutcoord = 0.;
        k = 0;
        output(j, v, x, k);
        ++k;
        while (j < tmaxch) {
            ++j;
            auto_list_update( x, x_update, list);
            eom_aoup(v, x, f, a, temp, list, F);

            if (j >= toutcoord) {
                output(j, v, x, k);
                toutcoord += tbitch;
                ++k;
            }
        }
    }

    int counthazure = 0;
    for (int i = 0; i < Np; ++i) {
        if (x[i][0] * x[i][0] + x[i][1] * x[i][1] > R * R)
            counthazure++;
    }
    char filename[128];

    ofstream file;
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/kekkalo%.3fm%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), Mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << endl;
    file.close();

    cout << "done" << endl;
    return 0;
}
