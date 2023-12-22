#ifndef ANALYZEPART_SINGLE_PARTICLE
#define ANALYZEPART_SINGLE_PARTICLE
char file_xcor[128], file_vcor[128], file_msd[128], file_lzcor[128];
// 要求する変数はNp,R,mgn,foldername1;

static constexpr double Np_act_1 = 1. / Np_act;
static double           lz0[Np];

void calc_corrini(double (*x)[dim], double (*x0)[dim], double (*v1)[dim],
                  double (*v)[dim]) {
    double xcor = 0., vcor = 0., lzcor = 0.;
    for (int i = 0; i < Np_act; i++) {
        x0[i][0] = x[i][0];
        x0[i][1] = x[i][1];
        xcor += (x[i][0] * x[i][0] + x[i][1] * x[i][1]) * Np_act_1;
        v1[i][0] = v[i][0];
        v1[i][1] = v[i][1];
        vcor += (v[i][0] * v[i][0] + v[i][1] * v[i][1]) * Np_act_1;
        lz0[i] = (x[i][0] * v[i][1] - x[i][1] * v[i][0]);
        lzcor += lz0[i] * lz0[i] * Np_act_1;
    }
    ofstream file;
    snprintf(file_xcor, 128, "./%s/xcor_R%.3f_m%.3f.dat", foldername1, R, mgn);
    file.open(file_xcor); // append
    file << 0 << "\t" << xcor << endl;
    file.close();
    snprintf(file_vcor, 128, "./%s/vcor_R%.3f_m%.3f.dat", foldername1, R, mgn);
    file.open(file_vcor); // append
    file << 0 << "\t" << vcor << endl;
    file.close();
    snprintf(file_msd, 128, "./%s/msd_R%.3f_m%.3f.dat", foldername1, R, mgn);
    file.open(file_msd); // append
    file << "#t msd" << endl;
    file.close();
    snprintf(file_lzcor, 128, "./%s/lzcor_R%.3f_m%.3f.dat", foldername1, R,
             mgn);
    file.open(file_lzcor); // append
    file << 0 << "\t" << lzcor << endl;
    file.close();
}

// iniに加えてdt;
void calc_corr(double (*x)[dim], double (*x0)[dim], double (*v1)[dim],
               double (*v)[dim], unsigned long long j) {
    double dr, xcor = 0., vcor = 0., msd = 0., lzcor = 0.;

    for (int i = 0; i < Np_act; ++i) {
        for (int j = 0; j < dim; ++j) {
            xcor += x0[i][j] * x[i][j] * Np_act_1;
            vcor += v1[i][j] * v[i][j] * Np_act_1;
            dr = x[i][j] - x0[i][j];
            msd += dr * dr * Np_act_1;
            lzcor +=
                lz0[i] * (x[i][0] * v[i][1] - x[i][1] * v[i][0]) * Np_act_1;
        }
    }
    ofstream file;
    file.open(file_xcor, std::ios::app); // append
    file << j * dt << "\t" << xcor << endl;
    file.close();
    file.open(file_vcor, std::ios::app); // append
    file << j * dt << "\t" << vcor << endl;
    file.close();
    file.open(file_msd, std::ios::app); // append
    file << j * dt << "\t" << msd << endl;
    file.close();
    file.open(file_lzcor, std::ios::app); // append
    file << j * dt << "\t" << lzcor << endl;
    file.close();
}
char file_fais[128], file_del_theta[128], file_haikou_theta[128],
    file_nhat[128], file_nhalf[128];
// 要求する変数はfoldername1,R;
void calc_fai_ini() {
    ofstream file;
    snprintf(file_fais, 128, "./%s/fais_R%.3f.dat", foldername1, R);
    file.open(file_fais);
    file << "# t fai lzb pib" << endl;
    file.close();
    snprintf(file_del_theta, 128, "./%s/del_theta_R%.3f.dat", foldername1, R);
    file.open(file_del_theta);
    file << "#t del_theta" << endl;
    file.close();
    snprintf(file_haikou_theta, 128, "./%s/haikou_theta_R%.3f.dat", foldername1,
             R);
    file.open(file_haikou_theta);
    file << "#t del_theta_t" << endl;
    file.close();
    snprintf(file_nhat, 128, "./%s/nhat_R%.3f.dat", foldername1, R);
    file.open(file_nhat);
    file << "#t nhat/np" << endl;
    file.close();
#if FLAG_2PARA_EQ == 1
    snprintf(file_nhalf, 128, "./%s/nhalf_R%.3f.dat", foldername1, R);
    file.open(file_nhalf);
    file << "#t nfst nscd" << endl;
    file.close();
#endif
}
inline double usrr_abs(double x) { return (x * ((x > 0.) - (x < 0.))); }
// 要求する変数はiniに加えてNp,R,Np_act_1,M_PI2,v0,dt;
void calc_fai(double (*x)[dim], double (*v)[dim], double *theta_i,
              long long j) { // para[0]:fai para[1]:vt* para[2];om*;
    double sum_vt = 0., sum_v = 0., vt, r, r2, sum_vrl[2] = {0., 0.},
           sum_lzrl[2] = {0., 0.}, para3[3] = {0, 0, 0}, del_theta_td = 0.,
           nhat = 0., sum_om = 0., haikou = 0.;
    int           cnt = 0, cnt_half = 0;
    static double theta_past[Np_act];
    static double delta_theta = 0.;
    // int                     count[2] = {0, 0};
    // constexpr double bun = 1 / (para3_tbit / para3_bitlonch / dt);
    constexpr double bun_kai = 1 / (1 - M_2_PI), R_1cor2 = (R - 1) * (R - 1),
                     R_2cor2 = (R - 2) * (R - 2), Rhalf = R / M_SQRT2;
    for (int i = 0; i < Np_act; ++i) {
        r2 = x[i][0] * x[i][0] + x[i][1] * x[i][1];
        r = sqrt(r2);
        vt = ((x[i][0]) * v[i][1] - x[i][1] * v[i][0]) / r2;
        sum_vt += usrr_abs(vt * r);
        sum_v += sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1]);
        sum_vrl[0] += ((vt > 0) - (vt < 0)) * Np_act_1;
        sum_lzrl[0] += vt * r2 * Np_act_1;
        double atan_i = atan2(x[i][1], x[i][0]);
        nhat += v0 * r * sin(theta_i[i] - atan_i) * Np_act_1;
        if (r >= Rhalf) {
            cnt_half++;
        }
        if (r2 > R_1cor2) {
            cnt++;
            del_theta_td += atan_i - theta_past[i] -
                            M_PI2 * (int) ((atan_i - theta_past[i]) * M_1_PI);
            theta_past[i] = atan_i;
            haikou +=theta_i[i] - atan_i- M_PI2 * (int) ((theta_i[i] - atan_i) * M_1_PI);
            sum_om += vt;
        } else if (r2 > R_2cor2) {
            theta_past[i] = atan2(x[i][1], x[i][0]);
        }
    }
#if FLAG_2PARA_EQ == 1
    int cnt_falf2 = 0;
    for (int i = Np_act; i < Np; i++) {
        r = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        if (r >= Rhalf)
            cnt_falf2++;
    }
#endif
    para3[0] += (sum_vt / sum_v - M_2_PI) * bun_kai;
    para3[1] += sum_vrl[0];
    para3[2] += sum_lzrl[0];
    double bun = 1. / cnt;
    del_theta_td *= bun;
    delta_theta += del_theta_td;
    sum_om *= bun;
    haikou *= bun;
    ofstream file;
    file.open(file_fais, std::ios::app); // append
    file << j * dt << "\t" << para3[0] << "\t" << para3[1] << "\t" << para3[2]
         << "\t" << sum_om << endl;
    file.close();
    file.open(file_del_theta, std::ios::app);
    file << j * dt << "\t" << delta_theta << endl;
    file.close();
    file.open(file_haikou_theta, std::ios::app);
    file << j * dt << "\t" << haikou << endl;
    file.close();
    file.open(file_nhat, std::ios::app);
    file << j * dt << "\t" << nhat << endl;
    file.close();
#if FLAG_2PARA_EQ == 1
    file.open(file_nhalf, std::ios::app);
    file << j * dt << "\t" << cnt_half * Np_act_1 << "\t"
         << cnt_falf2 / Np_lange << endl;
    file.close();
#endif
}
#endif