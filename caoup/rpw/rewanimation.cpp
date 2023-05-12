#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "BM.h"

#define Np     4 // 8192// 4096// 4の倍数であること;
#define Nn     4
#define R   10.  //50.//28.86751346//  32.//64.//128.// lo=0.5,N1000//3.
                //  //R=sqrt(np/4/lo);
// ,0.1より大きいこと;
#define tmax     500 // 適当;
#define tmaxlg   0 // 緩和時間は10たうとする;
#define tbit 1.//ファイルを出すときの間隔;
#define dtlg     0.00001
#define dt       0.00001
#define temp     100.         // v0^2=2D/tau,ここではDを入れること;
#define dim      2           // 変えるときはEomを変えること;
#define cut      1.122462048 // 3.
#define skin     1.5
#define tau      100.
#define ensemble 1
// #define polydispersity 0.2 コードも変える;
#define folder_name "rew02"
#define mgn         0.// Omega=omega/tau,ここではomegaを入れること;
//#define radios 1. //粒径の平均値を変えるときはヒストグラムの変え方も変えること:現在は1;
// v4:lohistをNpで割らなくした;
// v五:ディレクトリ名でRを先にした;
//stwv2:ティレク取りメイをloに、ヒストグラムの取り方を外からに;
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
void set_diameter(double *a) {
    for (int i = 0; i < Np; i++)
        a[i] = 0.5;
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
double step_funk(double x){
    if(x>0)return x;
    else return 0;
}
void eom_aoup(double (*v)[dim], double (*x)[dim], double (*f)[dim], double *a,
              double temp0, int (*list)[Nn], double (*F)[dim]) {
    double tauinv = dt / tau;
    double F0[dim],v0[2],fiw[dim],ffw[2],du2;
    double fluc = sqrt(temp0 / dt),eri;
    //calc_force(x, f, a, list);
    double ri,riw,w2,w6,w12,aij,dUr;
    for (int i = 0; i < Np; i++) {
        fiw[0]=0.;
        fiw[1]=0.;
        ffw[0]=0.;
        ffw[1]=0.;
// /*force bitween wall;
        
        ri=sqrt(x[i][0]*x[i][0]+x[i][1]*x[i][1]);
        riw=R+0.5-ri;
        if(riw<cut){
            aij=0.5+a[i];

            w2=aij*aij/(riw*riw);
            w6=w2*w2*w2;
            w12=w6*w6;
            dUr=(-48.*w12+24.*w6)/(riw*ri);
            fiw[0]=dUr*x[i][0];
            fiw[1]=dUr*x[i][1];
            du2=-6.*w6/(ri*ri);
            ffw[0]=du2*F[i][0]*x[i][0]*x[i][0];
            ffw[1]=du2*x[i][1]*F[i][1]*x[i][1];
        }
///till here;*/
        F0[0] = F[i][0];
        F0[1] = F[i][1];
        F[i][0] += (-F0[0] - F0[1] * mgn + fluc * gaussian_rand()+ffw[0]*tau) * tauinv;
        F[i][1] += (-F0[1] + F0[0] * mgn + fluc * gaussian_rand()+ffw[1]*tau) * tauinv;

        v[i][0] = F[i][0] +fiw[0]+f[i][0];
        v[i][1] = F[i][1] +fiw[1]+f[i][1];
        x[i][0] += v[i][0] * dt;
        x[i][1] += v[i][1] * dt;
    }
}

void ini_hist(double *hist, int Nhist) {
    for (int i = 0; i < Nhist; i++) {
        hist[i] = 0.;
    }
}





void output(int k, double (*v)[dim], double (*x)[dim], int l) {
    char   filename[128];
    double        Mgn = mgn / tau;
    std::ofstream file;
    sprintf(filename,
            "./%s_animelo%.2fv0%.1ftau%.3fm%.3f/"
            "tyouwaenn_lo%.3f_tau%.3f_m%.3f_t%d.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), tau, Mgn, l);
    file.open(filename, std::ios::app); // append
    for (int i = 0; i < Np; i++) {
        file << k * dt << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << v[i][0]
             << "\t" << v[i][1]  << std::endl;
    }
    file.close();
}
void out_setup() {
    char          filename[128];
    double        Mgn = mgn / tau, v0 = temp / tau;
    std::ofstream file;
    sprintf(filename, "./%s_animelo%.2fv0%.1ftau%.3fm%.3f/setupr%fm%f.dat",
            folder_name, Np * 0.25 / (R * R), v0, tau, Mgn, R, Mgn);
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
    file << "type=" << 2 << std::endl;
    file<<"2DkaraD"<<std::endl;
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


int main() {
    double x[Np][dim], x_update[Np][dim], v[Np][dim], f[Np][dim], a[Np],
        F[Np][dim], x0[Np][dim],v0[Np][dim];
    int    list[Np][Nn];
    int    counthistv_theta = 0, countout = 0;
    double tout = 0.01, toutcoord = 0, U, disp_max = 0.0;

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

    j = 0;
    double ttemp = 5. * temp;
    if (ttemp / tau < 5)
        ttemp = 5. * tau;
    j = 0;
    while (j * dtlg < 0) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, ttemp, list, F);
        
    }
 
    j = 0;
    while (j * dt < tmaxlg) {
        j++;
        auto_list_update(&disp_max, x, x_update, list);
        eom_aoup(v, x, f, a, temp, list, F);
    }

double tmaxch=tmax/dt,tbitch=tbit/dt;
    
    for (int i = 0; i < ensemble; i++) {
        
        j = 0;
        toutcoord = 0.;
        k = 0;
        output(j, F, x, k);
        k++;
        while (j < tmaxch) {
            j++;
            auto_list_update(&disp_max, x, x_update, list);
            eom_aoup(v, x, f, a, temp, list, F);

            if (j  >= toutcoord) {
                output(j, F, x, k);
                toutcoord += tbitch;
                k++;
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
    sprintf(filename, "./%slo%.2fv0%.1ftau%.3fm%.3f/kekkalo%.3fm%.3f.dat",
            folder_name, Np * 0.25 / (R * R), temp / tau, tau, Mgn,
            Np * 0.25 / (R * R), Mgn);
    file.open(filename, std::ios::app); // append

    file << counthistv_theta << " " << counthazure << std::endl;
    file.close();
    
    std::cout << f[0][0]  << std::endl;
    return 0;
}
