
#include <chrono>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#define Np             65536 // 累乗であること;
#define Nn             100
#define tmtimes        400  // ファイルを出す回数;
#define tmaxlg         200  // 緩和時間は10たうとする;
#define tmaxka         1000 // 緩和時間;
#define tmaxani        500  //>tmaxの時プログラムを変更すること;
#define tbitani        2
#define ratf           1.
#define dim            2           // 変えるときはEomを変えること;
#define cut            1.122462048 // 3.
#define skin           1.5
#define dtlg           1e-4
#define dt             1e-2
#define folder_name    "Iapr2n2e16pir1" // 40文字程度で大きすぎ;
#define msdbit         1.1
#define msdini         0.01
#define polydispersity 0.3 // コードも変える;
using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.
#define lo 0.7 // コンパイル時に代入する定数;
// コンパイル時に-D{変数名}={値}　例:-Dbit=80　とすること;
#define v0 1.

static constexpr double tau = 40.;
static constexpr double mass = 80.;
static constexpr double mgn = 0.; // 有限にするときはコードを変える;
constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 5000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}

static constexpr double cut2 = cut * cut;
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double Mg = mgn * dt;
static constexpr double L = usr_sqrt(Np * M_PI * 0.25*(1.+polydispersity*polydispersity) / lo);
static constexpr double L_inv = 1. / L;
static constexpr double L_2 = L / 2.;
static constexpr double Np_1 = 1. / Np;
static constexpr double fconst = ratf * 48.;
static constexpr int    dim2 = 2 * dim;
static constexpr double to = (tau > 10) ? tau * 2. : 20;
static constexpr double tmax = to * tmtimes;

void usr_sincos(double kaku, double *x) { // x[0]がcos,x[1]issin;
                                          // 制度は10^-13程度;
    constexpr double waru[8] = {1.0 / (3 * 4 * 5 * 6 * 7 * 8 * 9 * 10),
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
double unif_rand(double left, double right) {
    return left + (right - left) * rand() / RAND_MAX;
}

double gaussian_rand(void) {
    static bool   iset = true;
    static double gset;
    double        fac, rsq, v1, v2;

    if (iset) {
        do {
            v1 = unif_rand(-1, 1);
            v2 = unif_rand(-1, 1);
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);

        gset = v1 * fac;
        iset = false;
        return v2 * fac;
    } else {
        iset = true;
        return gset;
    }
}
void ini_hex(double (*x)[dim2]) {
    int num_x = (int) sqrt(Np);
    int num_y = (int) sqrt(Np);
    if (num_x * num_y != Np) {
        num_x++;
        num_y++;
    }
    int    k = 0;
    double shift;
    for (int j = 0; j < num_y; j++) {
        for (int i = 0; i < num_x; i++) {
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
void ini_hist(double *his, int koo) {
    for (int i = 0; i < koo; i++) {
        his[i] = 0.;
    }
}

inline double perio(double x) { return L * floor(x * L_inv); }
inline double pri_fce(double x) {
    x -= L * floor((x + L_2) * L_inv);
    return x;
}
void calc_force(double (*x)[dim2], double (*f)[dim], double *a,
                int (*list)[Nn]) {
    double dx, dy, dr2, dUr, w2, w6, aij /*w12*/;
    ini_array(f);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = pri_fce(x[i][0] - x[list[i][j]][0]);
            dy = pri_fce(x[i][1] - x[list[i][j]][1]);
            dr2 = dx * dx + dy * dy;
            aij = (a[i] + a[list[i][j]]) * (a[i] + a[list[i][j]]);
            if (dr2 < cut2 * aij) {
                w2 = aij / dr2;
                w6 = w2 * w2 * w2;
                // w12 = w6 * w6;
                dUr = fconst * (-w6 + 0.5) * w6 / dr2; // polydispersity;
                f[i][0] -= dUr * dx;
                f[list[i][j]][0] += dUr * dx;
                f[i][1] -= dUr * dy;
                f[list[i][j]][1] += dUr * dy;
            }
        }
}

void eom_abp9(double (*v)[dim], double (*x)[dim2], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    constexpr double zeta = 1.0, ddt = dtlg;
    double           sico[2];
    constexpr double D = usr_sqrt(2. * ddt / tau), M_inv = ddt / mass;
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {

        theta_i[i] += D * gaussian_rand();
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
        v[i][0] += (-v[i][0] + v0 * sico[0] + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sico[1] + f[i][1]) * M_inv;
        x[i][0] += v[i][0] * ddt;
        x[i][1] += v[i][1] * ddt;
        x[i][0] -= perio(x[i][0]);
        x[i][1] -= perio(x[i][1]);
    }
}

void eom_langevin(double (*v)[dim], double (*x)[dim2], double (*f)[dim],
                  double *a, int (*list)[Nn], double *theta_i) {

    constexpr double zeta = 1.0, ddt = 1e-9;
    constexpr double fluc = usr_sqrt(4. * zeta * ddt);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {

        for (int j = 0; j < dim; j++) {
            v[i][j] += (-v[i][j] + f[i][j]) * ddt + fluc * gaussian_rand();
            x[i][j] += v[i][j] * ddt;
            x[i][j] -= perio(x[i][j]);
        }
    }
}
void eom_langevin_t(double (*v)[dim], double (*x)[dim2], double (*f)[dim],
                    double *a, int (*list)[Nn], double temp) {

    double zeta = 1.0;
    double fluc = sqrt(2. * zeta * temp * dtlg);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {

        for (int j = 0; j < dim; j++) {
            v[i][j] += (-v[i][j] + f[i][j]) * dtlg + fluc * gaussian_rand();
            x[i][j] += v[i][j] * dtlg;
            x[i][j] -= perio(x[i][j]);
        }
    }
}
void ini_count(double (*x)[dim2]) {
    for (int i = 0; i < Np; i++) {
        x[i][2] = 0.;
        x[i][3] = 0.;
    }
}
void eom_abp1(double (*v)[dim], double (*x)[dim2], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double           ov[2];
    constexpr double D = usr_sqrt(2. * dt / tau), M_inv = dt / mass;
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        theta_i[i] += D * gaussian_rand();
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        // usr_sincos(theta_i[i], sico);
        v[i][0] += (-v[i][0] + v0 * cos(theta_i[i]) + f[i][0]) * M_inv;
        v[i][1] += (-v[i][1] + v0 * sin(theta_i[i]) + f[i][1]) * M_inv;
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
        ov[0] = -perio(x[i][0]);
        ov[1] = -perio(x[i][1]);
        x[i][0] += ov[0];
        x[i][1] += ov[1];
        x[i][2] += ov[0];
        x[i][3] += ov[1];
    }
}
inline double usr_abs(double x) { return x * ((x > 0) - (x < 0)); }

void output_ini(double (*v)[dim], double (*x)[dim2], double *a) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s_animelo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyokkei.dat",
             folder_name, lo, mass, tau, v0);
    file.open(filename /* std::ios::app*/); // append
    file << tbitani << endl;
    for (int i = 0; i < Np; i++)
        file << a[i] * 2. << endl;
    file.close();
    snprintf(filename, 128,
             "./%s_coorlo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t0.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
}
void output(double (*v)[dim], double (*x)[dim2]) {
    static int l = 0;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128,
             "./%s_coorlo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
    l++;
}
void output_iniani(double (*v)[dim], double (*x)[dim2], double *a) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s_animelo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyokkei.dat",
             folder_name, lo, mass, tau, v0);
    file.open(filename /* std::ios::app*/); // append
    file << tbitani << endl;
    for (int i = 0; i < Np; i++)
        file << a[i] * 2. << endl;
    file.close();
    snprintf(filename, 128,
             "./%s_animelo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t0.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0] << "\t" << v[i][1]
             << endl;
    }
    file.close();
}
void output_ani(double (*v)[dim], double (*x)[dim2]) {
    static int l = 1;
    char       filename[128];
    ofstream   file;

    snprintf(filename, 128,
             "./%s_animelo%.2fMs%.3ftau%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn, l);
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
    int      test = snprintf(filename, 128,
                             "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
                                  "setupofst_lo%.3f_tau%.3f_m%.3f_t%d.dat",
                             folder_name, lo, mass, tau, v0, lo, tau, mgn, tmax);
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
    file << "cell list" << endl;
    file << "自動Np" << endl;
    file << "modNp" << endl;
    file << "L=" << L << endl;
    file << "pi? " << L * L * lo / Np << endl;
    file << "ratio" << ratf << endl;
    file << "to" << to << endl;
    file.close();
    if (test == -1)
        return false;
    else
        return true;
}
void ini_corr(double (*x)[dim2], double (*x0)[dim], double (*v1)[dim],
              double (*v)[dim], double *xcor, double *vcor) {
    for (int i = 0; i < Np; i++) {
        x0[i][0] = x[i][0];
        x0[i][1] = x[i][1];
        x[i][2] = 0.;
        x[i][3] = 0.;
        v1[i][0] = v[i][0];
        v1[i][1] = v[i][1];
    }
}
void calc_corr(double (*x)[dim2], double (*x0)[dim], double (*v1)[dim],
               double (*v)[dim], double *xcor, double *vcor, int k,
               double *msd) {
    double rsin[2], gc[2] = {0., 0.};
    for (int i = 0; i < Np; ++i) {
        gc[0] += (x[i][0] - x[i][2]) * Np_1;
        gc[1] += (x[i][1] - x[i][3]) * Np_1;
    }
    for (int i = 0; i < Np; i++) {
        rsin[0] = x[i][0] - x[i][2] - gc[0];
        rsin[1] = x[i][1] - x[i][3] - gc[1];
        xcor[k] += (x0[i][0] * rsin[0] + x0[i][1] * rsin[1]) * Np_1;
        vcor[k] += (v1[i][0] * v[i][0] + v1[i][1] * v[i][1]) * Np_1;
        msd[k] += ((rsin[0] - x0[i][0]) * (rsin[0] - x0[i][0]) +
                   (rsin[1] - x0[i][1]) * (rsin[1] - x0[i][1])) *
                  Np_1;
    }
}

void outputcorr(double *msd, double *vcor, double *t, int countout,
                double *msd2) {
    char     filename[128];
    double   v_theta;
    ofstream file;
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
             "xcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
             "vcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
             "msd_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, mass, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}
static constexpr int M = L / (cut + skin);
inline int           peri_cell(int m) {
    if (m < 0)
        return m + M;
    else if (m >= M)
        return m - M;
    else
        return m;
}
void cell_list(int (*list)[Nn], double (*x)[dim2]) {
    int              map_index, nx[Np][dim];
    constexpr int    m2 = M * M;
    constexpr double thresh2 = (cut * (1 + polydispersity) + skin) *
                               (cut * (1 + polydispersity) + skin),
                     bit = M / L;
    double dx, dy;
    int(*map)[Np] = new int[m2][Np];
    for (int i = 0; i < m2; ++i)
        map[i][0] = 0;

    for (int i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0]) * bit);
        nx[i][1] = (int) ((x[i][1]) * bit);
        for (int m = nx[i][1] - 1; m <= nx[i][1] + 1; ++m) {
            for (int l = nx[i][0] - 1; l <= nx[i][0] + 1; ++l) {
                map_index = peri_cell(l) + M * peri_cell(m);
                map[map_index][0]++;
                map[map_index][map[map_index][0]] = i;
            }
        }
    }
    // int km, j;
    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + M * nx[i][1];
        for (int k = 1; k <= map[map_index][0]; ++k) {
            // j = map[map_index][k];
            if (map[map_index][k] > i) {
                dx = pri_fce(x[i][0] - x[map[map_index][k]][0]);
                dy = pri_fce(x[i][1] - x[map[map_index][k]][1]);
                if ((dx * dx + dy * dy) < thresh2) {
                    list[i][0]++;
                    list[i][list[i][0]] = map[map_index][k];
                }
            }
        }
    }
    delete[] map;
}

void update(double (*x_update)[dim], double (*x)[dim2]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++)
            x_update[i][j] = x[i][j];
}

void calc_disp_max(double *disp_max, double (*x)[dim2],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp;
    for (int i = 0; i < Np; i++) {
        dx = pri_fce(x[i][0] - x_update[i][0]);
        dy = pri_fce(x[i][1] - x_update[i][1]);
        disp = dx * dx + dy * dy;
        if (disp > *disp_max)
            *disp_max = disp;
    }
}

void auto_list_update(double (*x)[dim2], double (*x_update)[dim],
                      int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static double    disp_max = skin * skin;
    constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    calc_disp_max(&disp_max, x, x_update);
    if (disp_max >= skin2) {
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
    double x[Np][dim2], v[Np][dim], theta[Np], a[Np], f[Np][dim], x0[Np][dim],
        v1[Np][dim], x_update[Np][dim], disp_max = 0;
    // int(*list)[Nn] = new int[Np][Nn];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    double tout = msdini, toutcoord = 0;
    int    k = 0, kcoord = 0;
    set_diameter(a);
    ini_hex(x);
    ini_array(v);
    ini_array(x_update);
    ini_array(f);
    ini_hist(theta, Np);
    char foldername[128];
    snprintf(foldername, 128, "%slo%.2fMs%.3ftau%.3fv0%.1f", folder_name, lo,
             mass, tau, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    snprintf(foldername, 128, "%s_coorlo%.2fMs%.3ftau%.3fv0%.1f", folder_name,
             lo, mass, tau, v0);
    const char *fname2 = foldername;
    mkdir(fname2, 0777);
    snprintf(foldername, 128, "%s_animelo%.2fMs%.3ftau%.3fv0%.1f", folder_name,
             lo, mass, tau, v0);
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
    long long int j = 0;

    while (j < 1e7) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_langevin(v, x, f, a, list, theta);
    }
    j = 0;

    std::cout << "passed kasanari!" << endl;
    unsigned long long int tmaxbefch = 5 / (dtlg);
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_langevin_t(v, x, f, a, list, 5.);
    }
    j = 0;
    while (j < 1e7) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp9(v, x, f, a, list, theta);
    }

    std::cout << "passed kakimaze!" << endl;
    j = 0;
    tmaxbefch = tmaxka / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
    }
    std::cout << "passed owari!" << endl;
    int                    ituibi = 0, toch = to / dt, tanibitch = tbitani / dt;
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
    ini_count(x);
    calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
    t[0] = 0.;
    ++kcoord;
    output_ini(v, x, a);
    output_iniani(v, x, a);
    ++k;
    while (j < tanimaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        if (j >= kanit) {
            output_ani(v, x);
            kanit += tanibitch;

            if (j >= toutcoord) {
                output(v, x);
                toutcoord += toch;
            }
        } //*/
        if (j >= tout) {
            calc_corr(x, x0, v1, v, msd, vcor, kcoord, msd2);
            t[kcoord] = j * dt;
            ++kcoord;
            tout *= msdbit;
        }
    }
    while (j < tmaxch) {
        ++j;
        auto_list_update(x, x_update, list);
        eom_abp1(v, x, f, a, list, theta);
        // make_v_thetahist(x, v, hist, hist2, lohist);

        if (j >= toutcoord) {
            output(v, x);
            toutcoord += toch;
            //*/
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
    }
    end = std::chrono::system_clock::now(); // 計測終了時間
    char     filename[128];
    ofstream file2;
    snprintf(filename, 128,
             "./%slo%.2fMs%.3ftau%.3fv0%.1f/kekkalo%.3fm%.3f.dat", folder_name,
             lo, mass, tau, v0, lo, mgn);
    file2.open(filename, std::ios::app); // append
    file2 << ave << " " << maxnum << " " << endl;
    file2 << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count()
          << endl; // 処理に要した時間をミリ秒に変換
    file2.close();
    outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
