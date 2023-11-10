
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


#define Np      32
#define lo      0.1
#define Nn      100
#define tmax    0 // 973.686//2*100たうとする;<tmaxaniの時気をつける;
#define tmaxlg  400 // 緩和時間は10たうとする;
#define v0      1.
#define tau     20. // コンパイル時に-D{変数名}={値}　例:-Dtau=80　とすること;
#define mgn     0.  // Omega=omega/tau,ここではomegaを入れること;
#define M_ss 80.
#define ratio_f 1./24
#define tmaxani 100 //>tmaxの時プログラムを変更すること;
#define tbitani 0.1
#define dim     2           // 変えるときはEomを変えること;
#define cut     1.122462048 // 3.
// #define cutal       2.5
#define skin 1.5
#define dtlg 0.0001
#define dt   0.0001
// #define ramda       0.333
// #define c_j         -1.
#define const_u     0.
#define folder_name "ttttr0"
// #define delta_om    0
#define msdbit 1.1
#define msdini 0.01
// #define polydispersity 0.2 コードも変える;
using std::endl;
using std::max;
using std::min;
using std::ofstream;
// #define radios 1.
constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 100; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}
static constexpr double L = usr_sqrt(M_PI_4 * Np / lo);
static constexpr double L_inv = 1. / L;
static constexpr double L_2 = L / 2.;
static constexpr double cut2 = cut * cut;
static constexpr double M_PI2 = 2. * M_PI;
static constexpr double Mg = mgn * dt;
static constexpr double Np_1 = 1. / Np;
// static constexpr double cutl = (cut > cutal) ? cut : cutal;
// static constexpr double cuts = (cut > cutal) ? cutal : cut;
static constexpr double cutl2 = cut * cut;
// static constexpr double cuts2 = cuts * cuts;
// static constexpr double cutal2 = cutal * cutal;
static constexpr int M = L / (cut + skin);
static constexpr double M_1=dt/M_ss;
static constexpr double fcons=ratio_f*48.;

// x[0]がcos,x[1]issin;
void usr_sincos(double kaku, double *x) {
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

void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < dim; ++j)
            x[i][j] = 0.0;
}
void ini_fce(double (*k)[3]) {
    for (int i = 0; i < Np; i++) {
        k[i][0] = 0.;
        k[i][1] = 0.;
        k[i][2] = 0.;
    }
}
inline double perio(double x) {
    return L * floor(x * L_inv);
}
inline double pri_fce(double x) {
    x -= L * floor((x + L_2) * L_inv);
    return x;
}
void calc_force(double (*x)[dim], double (*f)[dim], double *a, int (*list)[Nn]) { // dimdp
    double dx, dy, dr2, dUr,  w2, w6, /*w12,*/  aij;
    ini_array(f);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = pri_fce(x[i][0] - x[list[i][j]][0]);
            dy = pri_fce(x[i][1] - x[list[i][j]][1]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < cutl2) {
                aij = (a[i] + a[list[i][j]]);
                w2 = aij * aij / dr2;
                w6 = w2 * w2 * w2;
                // w12 = w6 * w6;
                if (dr2 < cut2) {
                    dUr = fcons*(- w6 + 0.5) * w6 / dr2 /* -12. * w12 / dr2*/;
                    f[i][0] -= dUr * dx;
                    f[list[i][j]][0] += dUr * dx;
                    f[i][1] -= dUr * dy;
                    f[list[i][j]][1] += dUr * dy;
                }
            }
        }
}

void eom_abp9(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double ddt = 0.0000001, D = sqrt(2. * ddt / 0.01), ri, riw, aij, w2, w6,
           dUr, fiw[dim], sinco[2];
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {

        // till here*/
        theta_i[i] += D * gaussian_rand();
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sinco);
        v[i][0] = sinco[0] + f[i][0];
        v[i][1] = sinco[1] + f[i][1];
        x[i][0] += v[i][0] * ddt;
        x[i][1] += v[i][1] * ddt;
        x[i][0] -= perio(x[i][0]);
        x[i][1] -= perio(x[i][1]);
    }
}
void eom_abp8(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double        ri, riw, aij, w2, w6, dUr, fiw[dim], sico[2];
    static double D = sqrt(2. * dt / 0.01);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {

        // till here*/
        theta_i[i] += D * gaussian_rand() + Mg ;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i], sico);
        v[i][0] += (-v[i][0]+5. * sico[0] + f[i][0] )*M_1;
        v[i][1] += (-v[i][0]+5. * sico[1] + f[i][1] )*M_1;
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
        x[i][0] = perio(x[i][0]);
        x[i][1] = perio(x[i][1]);
    }
}

void eom_abp1(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              int (*list)[Nn], double *theta_i) {
    double ri, riw, aij, w2, w6, dUr, fiw[dim], sico[2];
    constexpr static double D = usr_sqrt(2. * dt / tau);
    calc_force(x, f, a, list);
    for (int i = 0; i < Np; i++) {
        // till here*/
        theta_i[i] += D * gaussian_rand() + Mg ;
        theta_i[i] -= (int) (theta_i[i] * M_1_PI) * M_PI2;
        usr_sincos(theta_i[i] , sico);
        v[i][0] += (-v[i][0]+v0 * sico[0] + f[i][0])*M_1 ;
        v[i][1] += (-v[i][0]+v0 * sico[1] + f[i][1])*M_1 ;
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
        x[i][0] -= perio(x[i][0]);
        x[i][1] -= perio(x[i][1]);
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; ++i) {
        hist[i] = 0.;
    }
}

void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s_coorlo%.2ftau%.3fm%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, lo, tau, mgn, v0, lo, tau, mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << endl;
    }
    file.close();
}

void output_ani(int k, double (*v)[dim], double (*x)[dim], int l,
                double *theta) {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128,
             "./%s_animelo%.2ftau%.3fm%.3fv0%.1f/"
             "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
             folder_name, lo, tau, mgn, v0, lo, tau, mgn, l);
    file.open(filename /* std::ios::app*/); // append
    for (int i = 0; i < Np; ++i) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1] << "\t" << theta[i] << endl;
    }
    file.close();
}
void out_setup() {
    char     filename[128];
    ofstream file;
    snprintf(filename, 128, "./%slo%.2ftau%.3fm%.3fv0%.1f/setupl%fm%f.dat",
             folder_name, lo, tau, mgn, v0, L, mgn);
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
    file<<"L="<<L<<endl;
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
             "./%slo%.2ftau%.3fm%.3fv0%.1f/xcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%slo%.2ftau%.3fm%.3fv0%.1f/vcor_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << vcor[i] << endl;
    }
    file.close();
    snprintf(filename, 128,
             "./%slo%.2ftau%.3fm%.3fv0%.1f/msd_lo%.3f_tau%.3f_m%.3f.dat",
             folder_name, lo, tau, mgn, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/); // append
    for (int i = 0; i < countout; ++i) {
        file << t[i] << "\t" << msd2[i] << endl;
    }
    file.close();
}
inline int peri_cell(int m) {
    if (m < 0)
        return M - 1;
    else if (m >= M)
        return 0;
    else
        return m;
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
void cell_list(int (*list)[Nn], double (*x)[dim]) {
    int                     map_index, nx[Np][dim];
    static constexpr int    m2 = M * M;
    static constexpr double thresh2 = (cutl + skin) * (cutl + skin),
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
    int km, j;
    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + M * nx[i][1];
        for (int k = 1, km = (map[map_index][0]); k <= km; ++k) {
            j = map[map_index][k];
            if (j > i) {
                dx = pri_fce(x[i][0] - x[j][0]);
                dy = pri_fce(x[i][1] - x[j][1]);
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
        dx = pri_fce(x[i][0] - x_update[i][0]);
        dy = pri_fce(x[i][1] - x_update[i][1]);

        disp = dx * dx + dy * dy;
        if (disp > *disp_max)
            *disp_max = disp;
    }
}

void auto_list_update(double *disp_max, double (*x)[dim],
                      double (*x_update)[dim], int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = 0.8 * skin2;
    calc_disp_max(&(*disp_max), x, x_update);
    if (*disp_max > skin2) {
        list_verlet(list, x);
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
        v1[Np][dim], x_update[Np][dim], disp_max = 0., Fq[Np][3];
    // int(*list)[Nn] = new int[Np][Nn];
    int           list[Np][Nn];
    int           counthistv_theta = 0, countout = 0;
    double        tout = msdini, toutcoord = 0;
    int           k = 0, kcoord = 0;
    long long int j = 0;
    set_diameter(a);
    ini_hex(x);
    ini_array(v);
    ini_array(f);
    ini_array(x_update);
    ini_hist(theta, Np);
    char foldername[128];
    snprintf(foldername, 128, "%slo%.2ftau%.3fm%.3fv0%.1f", folder_name, lo,
             tau, mgn, v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    char foldername2[128];
    snprintf(foldername2, 128, "%s_coorlo%.2ftau%.3fm%.3fv0%.1f", folder_name,
             lo, tau, mgn, v0);
    const char *fname2 = foldername2;
    mkdir(fname2, 0777);
    snprintf(foldername, 128, "%s_animelo%.2ftau%.3fm%.3fv0%.1f", folder_name,
             lo, tau, mgn, v0);
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
    while (j < 1e7) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp9(v, x, f, a, list, theta, Fq);
    }
    j = 0;

    std::cout << "passed kasanari!" << endl;
    int tmaxbefch = L / (dt * 5);
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp8(v, x, f, a, list, theta, Fq);
    }

    std::cout << "passed kakimaze!" << endl;
    j = 0;
    tmaxbefch = tmaxlg / dt;
    while (j < tmaxbefch) {
        ++j;
        auto_list_update(&disp_max, x, x_update, list);
        eom_abp1(v, x, f, a, list, theta, Fq);
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
        eom_abp1(v, x, f, a, list, theta, Fq);
        // make_v_thetahist(x, v, hist, hist2, lohist);
        if (j >= kanit) {
            output_ani(j, v, x, kani, theta);
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
        eom_abp1(v, x, f, a, list, theta, Fq);
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
    end = std::chrono::system_clock::now(); // 計測終了時間
    char filename[128];

    ofstream file;

    snprintf(filename, 128, "./%slo%.2ftau%.3fm%.3fv0%.1f/kekkalo%.3fm%.3f.dat",
             folder_name, lo, tau, mgn, v0, lo, mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << " " << ave << " "
         << maxnum << " " << endl;
    file << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()
         << endl; // 処理に要した時間をミリ秒に変換
    file.close();

    // outputhist(hist, counthistv_theta, lohist, hist2);
    // outputcorr(msd, vcor, t, countout, msd2);
    std::cout << "done" << endl;
    return 0;
}
