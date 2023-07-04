#include <fstream>
// #include <iomanip>
#include <chrono>
#include <iostream>
#include <math.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <sys/stat.h>

#define v0      1.0
#define lo      0.5
#define tau     50.
#define mgn     0.
#define M_ss    80.
#define dim     2
#define Np      1000
#define takebit 1 // 隣り合うファイルの間のファイル数を入力;
#define takefst 0 // とるファイルの最初;
#define bithist 1.
#define moji    "tyouwaenn" // 粒子位置のファイル名;
#define folder  "ppn4"      // 粒子位置のフォルダ名;
#define folder2 "pe"        // ダスフォルダ名;
#define qmax    50

// using std::cos;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
// using std::sin;

constexpr double usr_sqrt(double x) {
    double b = x;
    for (int i = 0; i < 1000; i++) {
        b = (b * b + x) / (2. * b);
    }
    return b;
}

static constexpr double L = usr_sqrt(Np / lo);
static constexpr double qbit = 2 * M_PI / L;
static constexpr int    qnmax = qmax / qbit;
static const double     qbb = 2 * M_PI / sqrt(Np / lo);
bool                    input_test(int k0) {
    char     filename[128];
    ifstream file;
    sprintf(filename,
                               "./%s_coorlo%.2fMs%.3ftau%.3fv0%.1f/"
                                                  "%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
                               folder, lo, M_ss, tau, v0, moji, lo, tau, mgn, k0);
    file.open(filename);
    if (!file) {
        std::cout << "Till here " << k0 << endl;
        return false;
    }

    file.close();
    return true;
}
void input(double (*x)[dim], double (*vl)[dim]) {
    static int times = 0;
    static int k0 = takefst;
    char       filename[128];
    ifstream   file;
    sprintf(filename,
            "./%s_coorlo%.2fMs%.3ftau%.3fv0%.1f/"
            "%s_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder, lo, M_ss, tau, v0, moji, lo, tau, mgn, k0);
    file.open(filename);
    if (!file) {
        cout << "You Failed! " << times * takebit + takefst << " " << k0
             << filename << endl;
        //        return false;
    }

    for (int i = times * Np; i < (times + 1) * Np; i++)
        file >> x[i][0] >> x[i][1] >> vl[i][0] >> vl[i][1];

    file.close();
    times++;
    k0 += takebit;
}

inline double usr_abs(double x) { return x * ((x > 0) - (x < 0)); }

void output(double *s2, double *c2, double(*omp), double *omo) {
    // int           Nphist = R / bithist + 1;
    char     filename[128];
    ofstream file;
    sprintf(filename,
            "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
            "sq_lo%.3f_tau%.3f_m%.3f.dat",
            folder2, lo, M_ss, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for (int i = 0; i <= qnmax; i++)
        file << i * qbit << "\t" << (c2[i] + s2[i]) << endl;
    file.close();
    sprintf(filename,
            "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
            "omo_lo%.3f_tau%.3f_m%.3f.dat",
            folder2, lo, M_ss, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for (int i = 0; i <= qnmax; i++)
        file << i * qbit << " " << omo[i] << endl;
    file.close();
    sprintf(filename,
            "./%slo%.2fMs%.3ftau%.3fv0%.1f/"
            "omp_lo%.3f_tau%.3f_m%.3f.dat",
            folder2, lo, M_ss, tau, v0, lo, tau, mgn);
    file.open(filename /*,std::ios::app*/);
    for (int i = 0; i < qnmax; i++)
        file << i * qbit << "\t" << omp[i] << endl;
    file.close();
}
void calc_q_enn(double (*x)[2], double (*v)[2], double *c2, double *s2, int ens,
                double *omp, double *omoth) {
    int    zenn = ens * Np, nal, nal2, Nm;
    double cq, sq, jx[2] = {0., 0.}, jy[2] = {0., 0.}, co = 0., si = 0., qi,
                   jp[2] = {0., 0.}, jo[2] = {0., 0.}, men, sqe[2] = {0., 0.};
    for (int nx = 0; nx <= qnmax; nx++)
        for (int ny = 0; ny <= qnmax; ny++) {
            nal2 = (int) (nx * nx + ny * ny);
            if (nal2 < qnmax * qnmax + 1 && nal2 != 0) {
                nal = (int) sqrt(nal2);
                men = 1. / (M_PI_2 * (nal + 0.5) * zenn*nal2);
                for (int j = 0; j < ens; ++j) {
                    co = 0.;
                    si = 0.;
                    jx[0] = 0.;
                    jx[1] = 0.;
                    jy[0] = 0.;
                    jy[1] = 0.;
                    for (int i = Np * j, Nm = (j + 1) * Np; i < Nm; i++) {
                        ///(0.5*pi*int(sqrt(nx*nx+ny*ny)))/Np/ens
                        qi = qbit * (nx * x[i][0] + ny * x[i][1]);
                        cq = cos(qi);
                        sq = sin(qi);
                        jx[0] += v[i][0] * cq;
                        jx[1] += v[i][0] * sq;
                        jy[0] += v[i][1] * cq;
                        jy[1] += v[i][1] * sq;
                        co += cq;
                        si += sq;
                    }
                    c2[nal] += co * co * men*nal2;
                    s2[nal] += si * si * men*nal2;
                    jp[0] = jx[0] * nx + jy[0] * ny;
                    jp[1] = jx[1] * nx + jy[1] * ny;
                    jo[0] = jx[0] * ny - jy[0] * nx;
                    jo[1] = jx[1] * ny - jy[1] * nx;
                    omp[nal] += (jp[0] * jp[0] + jp[1] * jp[1]) * men ;
                    omoth[nal] += (jo[0] * jo[0] + jo[1] * jo[1]) * men;
                }
                // 円環平均の場合円環の幅の面積で割る;
            }
        }
}
void calc_q_x(double (*x)[2], double (*v)[2], double *c2, double *s2, int ens,
              double *omp, double *omoth) {
    int    zenn = ens * Np, Nm;
    double cq, sq, jx[2] = {0., 0.}, jy[2] = {0., 0.}, co = 0., si = 0., qi,
                   sqe[2] = {0., 0.};
    double bun = 1. / (zenn * Np);
    for (int ny = 0; ny <= qnmax; ny++) {
        sqe[0] = 0.;
        sqe[1] = 0.;

        for (int j = 0; j < ens; ++j) {
            jx[0] = 0.;
            jx[1] = 0.;
            jy[0] = 0.;
            jy[1] = 0.;
            co = 0.0;
            si = 0.0;
            for (int i = Np * j, Nm = (j + 1) * Np; i < Nm; i++) {
                ///(0.5*pi*int(sqrt(nx*nx+ny*ny)))/Np/ens
                qi = qbit * (ny * x[i][0]);
                cq = cos(qi);
                sq = sin(qi);
                jx[0] += v[i][0] * cq;
                jx[1] += v[i][0] * sq;
                jy[0] += v[i][1] * cq;
                jy[1] += v[i][1] * sq;
                co += cq;
                si += sq;
            }
            c2[ny] += co * co * bun;
            s2[ny] += si * si * bun;

            omp[ny] += (jx[0] * jx[0] + jx[1] * jx[1]) * bun;
            omoth[ny] += (jy[0] * jy[0] + jy[1] * jy[1]) * bun;
        }
        // 円環平均の場合円環の幅の面積で割る;
    }
}
int main(int argc, char *argv[]) {
    std::chrono::system_clock::time_point start, end; // 型は auto で可
    start = std::chrono::system_clock::now();         // 計測開始時間
    // double tau, Rbit;
    // tau=atof(argv[1]);
    // Rbit=atof(argv[2]);
    int taketimes = 0, k0 = takefst;
    cout << L << endl;
    // int Nphist = R / bithist + 1;
    while (input_test(k0)) {
        taketimes++;
        k0 += takebit;
    }
    if (taketimes == 0)
        return -1;
    cout << taketimes << endl;
    cout << "abs(2)" << usr_abs(2.) << endl;
    cout << "abs(-2)" << usr_abs(-2.) << endl;
    cout << "abs(0)" << usr_abs(0.) << endl;
    cout << (1 >= 0.) << endl;

    double(*x)[dim] = new double[(int) (Np * taketimes)][dim];
    double(*v)[dim] = new double[(int) (Np * taketimes)][dim];
    double c2[qnmax + 5], s2[qnmax + 5], omp[qnmax + 5], omoth[qnmax + 5];
    for (int i = 0; i < qnmax + 5; i++) {
        c2[i] = 0.;
        s2[i] = 0.;
        omp[i] = 0.;
        omoth[i] = 0.;
    }
    char foldername[128];
    sprintf(foldername, "%slo%.2fMs%.3ftau%.3fv0%.1f", folder2, lo, M_ss, tau,
            v0);
    const char *fname = foldername;
    mkdir(fname, 0777);
    for (int i = 0; i < taketimes; i++)
        input(x, v);
    cout << "in done" << endl;
    calc_q_enn(x, v, c2, s2, taketimes, omp, omoth);
    output(s2, c2, omp, omoth);
    end = std::chrono::system_clock::now(); // 計測終了時間
    std::cout << "succeed!"
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換 endl;
    delete[] x;
    delete[] v;

    return 0;
}
