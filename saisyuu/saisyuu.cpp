#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "BM.h"
#define Np    1000
#define Nn    50
#define L     50
#define cut   2.
#define skin  1.5
#define out1  10
#define out2  1500
#define out3  100
#define dt    0.01
#define dtmd  0.001
#define temp1 5.
#define temp2 0.2
#define dim   2

using namespace std;

constexpr double pow_int(double x, int q) {
    double ans = x;
    for (int i = 1; i < q; i++) {
        ans *= x;
    }
    return ans;
}
constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 5000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}

static constexpr double cut2 = cut * cut;
static constexpr double L_inv = 1. / L;
static constexpr double L_2 = L / 2.;
static constexpr double Np_1 = 1. / Np;
static constexpr double heikatu =
    48. / pow_int(cut, 13) - 24. / pow_int(cut, 7);



void ini_hex(double (*x)[dim]) {
    int    num_x = (int) sqrt(Np) + 1;
    int    num_y = (int) sqrt(Np) + 1;
    int    i, j, k = 0;
    double shift;
    for (j = 0; j < num_y; j++) {
        for (i = 0; i < num_x; i++) {
            shift = (double) j * 0.5 - j / 2;
            x[i + num_x * j][0] = (shift + i) * L / (double) num_x;
            x[i + num_x * j][1] = j * L / (double) num_y;
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
    for (int i = 0; i < Np; ++i)
        a[i] = 0.5;
}
void ini_array(double (*f)[dim]) {
    for (int i = 0; i < Np; i++) {
        f[i][0] = 0.;
        f[i][1] = 0.;
    }
}

inline double perio(double x) { return L * floor(x * L_inv); }
inline double pri_fce(double x) {
    x -= L * floor((x + L_2) * L_inv);
    return x;
}

void calc_force(double (*x)[dim], double (*f)[dim], double *a, double *U,
                int (*list)[Nn]) {
    double              dx, dy, dr2, dUr, w2, w6, aij /*w12*/;
    static const double Ucut = 4. / pow(cut, 12) - 4. / pow(cut, 6),
                        U0cut = -48. / pow(cut, 13) + 24. / pow(cut, 7);
    ini_array(f);
    *U = 0;
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = pri_fce(x[i][0] - x[list[i][j]][0]);
            dy = pri_fce(x[i][1] - x[list[i][j]][1]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < cut2) {
                w2 = 1. / dr2;
                w6 = w2 * w2 * w2;
                // w12 = w6 * w6;
                dUr = -48. * (w6 - 0.5) * w6 * w2 +
                      heikatu * sqrt(w2); // polydispersity;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
                *U += 4. * (w6 - 1.) * w6 - Ucut - U0cut * (sqrt(dr2) - cut);
            }
        }
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
                dx -= L * floor((dx + 0.5 * L) / L);
                dy -= L * floor((dy + 0.5 * L) / L);
                dr2 = dx * dx + dy * dy;
                if (dr2 < thresh * thresh) {
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

void auto_list_update(double *disp_max, double (*x)[dim],
                      double (*x_update)[dim], int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max >= skin2) {
        list_verlet(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        *disp_max = skinini;
        // count = 0;
    }
}
void p_boundary(double (*x)[dim]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x[i][j] -= L * floor(x[i][j] / L);
}
void eom_langevin(double (*v)[dim], double (*x)[dim], double (*f)[dim],
                  double *a, double *U, double temp0, int (*list)[Nn]) {
    double zeta = 1.0;
    double fluc = sqrt(2. * zeta * temp0 * dt);

    calc_force(x, f, a, &(*U), list);
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++) {
            v[i][j] +=
                (-zeta * v[i][j] + f[i][j]) * dt + fluc * gaussian_rand();
            x[i][j] += v[i][j] * dt;
        }
    p_boundary(x);
}

void eom_md(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
            double *U, int (*list)[Nn]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++) {
            x[i][j] += (v[i][j] + 0.5 * f[i][j] * dtmd) * dtmd;
            v[i][j] += 0.5 * f[i][j] * dtmd;
        }
    calc_force(x, f, a, &(*U), list);
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++) {
            v[i][j] += 0.5 * f[i][j] * dtmd;
        }
    p_boundary(x);
}
void output_1(double (*x)[dim]) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "coord_1.dat");
    file.open(filename /*std::ios::app*/);
    for (int i = 0; i < Np; i++) {
        file << x[i][0] << " " << x[i][1] << endl;
    }

    file.close();
}
void output_2(double (*x)[dim]) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "coord_2.dat");
    file.open(filename /*std::ios::app*/);
    for (int i = 0; i < Np; i++) {
        file << x[i][0] << " " << x[i][1] << endl;
    }

    file.close();
}
void output(int k, double (*v)[dim], double U) {
    char   filename[128];
    double K = 0.0;

    ofstream file;
    snprintf(filename, 128, "energy.dat");
    file.open(filename, ios::app); // append
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            K += 0.5 * v[i][j] * v[i][j];

    cout << setprecision(6) << k * dtmd << "\t" << K / Np << "\t" << U / Np
         << "\t" << (K + U) / Np << endl;
    file << setprecision(6) << k * dtmd << "\t" << K / Np << "\t" << U / Np
         << "\t" << (K + U) / Np << endl;
    file.close();
}
int main() {
    double a[Np], disp_max = 0., U = 0.;
    int(*list)[Nn] = new int[Np][Nn];
    double(*x)[dim] = new double[Np][dim];
    double(*v)[dim] = new double[Np][dim];
    double(*f)[dim] = new double[Np][dim];
    double(*x_update)[dim] = new double[Np][dim];
    // int list[Np][Nn];

    int j = 0;
    set_diameter(a);
    ini_hex(x);
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    // output_1(x);
    while (j * dt < out1) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_langevin(v, x, f, a, &U, temp1, list);
    }
    output_1(x);
    cout << 1 << endl;
    j = 0;
    while (j * dt < out2) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_langevin(v, x, f, a, &U, temp2, list);
    }
    output_2(x);
    cout << 2 << endl;
    int out_count = 0;
    j = 0;
    while (j * dtmd < out3) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_md(v, x, f, a, &U, list);
        if (j * dtmd > out_count) {
            output(j, v, U);
            out_count += 2;
        }
    }
    cout << 3 << endl;
    delete[] x;
    delete[] v;
    delete[] x_update;
    delete[] f;
    delete[] list;
    return 0;
}