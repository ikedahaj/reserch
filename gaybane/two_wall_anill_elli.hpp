#ifndef TWO_WALL_ANILL
#define TWO_WALL_ANILL 1
void cell_list_anill(int (*list)[Nn], double (*x)[dim], double R_ini) {
    int              map_index, nx[dim];
    constexpr double threash2 =
        (cut + gay_kappa - 1 + skin) * (cut + gay_kappa - 1 + skin);
    double xlen_2 = (2. * R_ini + Rbit * R_ini) / 2.;
    int    Mx = (int) (xlen_2 * 2. / (cut + gay_kappa - 1 + skin));
    int    My = (int) (2. * R_ini /
                    (cut + gay_kappa - 1 + skin));  // M<=2R/(cutmax+skin)
    int    m2 = Mx * My;
    double R2 = 2. * R_ini, bitx = Mx / (xlen_2 * 2.),
           bity = My / (R2);  // ひとつのせるの幅の逆数;
    double dx, dy;
    // int(*map)[Np + 1] = new int[m2][Np + 1];
    std::vector<std::vector<int>> map(m2);
    for (int i = 0; i < m2; ++i) {
        map[i].reserve(20);
    }

    for (int i = 0; i < Np; ++i) {
        nx[0] = (int) ((x[i][0] + xlen_2) * bitx);
        nx[1] = (int) ((x[i][1] + R_ini) * bity);
        for (int m = max(nx[1] - 1, 0), mm = min(nx[1] + 1, My - 1);
             m <= mm; ++m) {
            for (int l = max(nx[0] - 1, 0), lm = min(nx[0] + 1, Mx - 1);
                 l <= lm; ++l) {
                map_index = l + Mx * m;
                map[map_index].emplace_back(i);
            }
        }
    }
    // int km, j;

    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = (int) ((x[i][0] + xlen_2) * bitx) +
                    Mx * (int) ((x[i][1] + R_ini) * bity);
        for (int j = 0, jm = map[map_index].size(); j < jm; ++j) {
            if (map[map_index][j] > i) {
                dx = (x[i][0] - x[map[map_index][j]][0]);
                dy = (x[i][1] - x[map[map_index][j]][1]);
                if ((dx * dx + dy * dy) < threash2) {
                    list[i][0]++;
                    list[i][list[i][0]] = map[map_index][j];
                }
            }
        }
    }
    // delete[] map;
}
void ver_list(int (*list)[Nn], double (*x)[dim], double R_ini) {
    double           dx, dy, dr2;
    constexpr double thresh2 =
        (cut + gay_kappa - 1 + skin) * (cut + gay_kappa - 1 + skin);
    for (int i = 0; i < Np; i++) list[i][0] = 0;

    for (int i = 0; i < Np; i++)
        for (int j = 0; j < Np; j++) {
            if (j > i) {
                dx = x[i][0] - x[j][0];
                dy = x[i][1] - x[j][1];
                dr2 = dx * dx + dy * dy;
                if (dr2 < thresh2) {
                    list[i][0]++;
                    list[i][(int) list[i][0]] = j;
                }
            }
        }
}

void auto_list_update_anill(double (*x)[dim], double (*x_update)[dim],
                            int (*list)[Nn], double R_ini) {
    // static int count = 0;
    // count++;
    static constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    static double           disp_max = skin2 + 100;
    calc_disp_max(&(disp_max), x, x_update);
    if (disp_max >= skin2) {
        w_list_anill(list, x, R_ini);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        disp_max = skinini;
        // count = 0;
    }
}
void update_coord_aniil(double (*x)[dim]) {
    for (int i = 0; i < Np; i++) {
        x[i][0] *= rat_anil;
        x[i][1] *= rat_anil;
    }
}
/// @brief 壁とgay-bane粒子の相互作用.2次元用;
/// @param x 粒子座標。粒子iについての二次元ザ行を送ること;
/// @param f 粒子の力。同上;
/// @param theta 粒子の角度.値私で良い;
/// @param f_theta 確度にかかる力.(&f[i])などと渡すこと.
inline void calc_force_wall_gay_bane_anill(double *x, double *f, double theta,
                                     double *f_theta,double R_ini) {
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) /
                                 (gay_kappa * gay_kappa + 1.),
                     gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) /
                                   (nth_root(gay_kappa_p, gay_M) + 1),
                     
                     gay_alpha2 = gay_alpha * gay_alpha,
                     gay_alpha_p2 = gay_alpha_p * gay_alpha_p;
    double r2, sa,cutwf2 =
                         (R_ini + 1.5 - cut - 1. / sqrt(1 - gay_alpha)) *
                         (R_ini + 1.5 - cut - 1. / sqrt(1 - gay_alpha));
    if (x[0] > 0)
        sa = center_rignt;
    else
        sa = center_left;
    x[0] -= sa;
    r2 = x[0] * x[0] + x[1] * x[1];
    if (r2 < cutwf2) {
        double r, sini, cosi, dr_dot_et, dt_dr_dot_dt, bun_alpha, sigma, xi;
        r = sqrt(r2);
        sini = sin(theta);
        cosi = cos(theta);
        dr_dot_et = x[0] * cosi + x[1] * sini;
        dt_dr_dot_dt = x[0] * sini - x[1] * cosi;
        bun_alpha = 1. / (r2 - gay_alpha2 * dt_dr_dot_dt * dt_dr_dot_dt);
        sigma = 1. / sqrt(1 - gay_alpha * (dr_dot_et * dr_dot_et * bun_alpha));
        xi = r - sigma + 1.;
        if (xi < cut) {
            double sigma3, bun_alpha_p, r_1, r2_1, epsilon13, xi_1, xi2, xi6,
                UWCA, dUWCA, dx_xi, dy_xi, dt_xi, epsilon1, epsilon1_N_1,
                dx_epsilon1, dy_epsilon1, dt_epsilon1, epsilon2, epsilon2_M_1,
                dr_epsilon2, dt_epsilon2;
            sigma3 = sigma * sigma * sigma;
            bun_alpha_p =
                1. / (r2 - gay_alpha_p2 * dt_dr_dot_dt * dt_dr_dot_dt);
            r_1 = 1. / r;
            r2_1 = r_1 * r_1;
            xi_1 = 1. / xi;
            xi2 = xi_1 * xi_1;
            xi6 = xi2 * xi2 * xi2;
            UWCA = 4. * xi6 * (xi6 - 1);
            dUWCA = 48. * xi6 * (0.5 - xi6) * xi_1;
            dx_xi = r_1 * x[0] + x[1] * sigma3 * gay_alpha * (1 - gay_alpha2) *
                                     dr_dot_et * dt_dr_dot_dt * bun_alpha *
                                     bun_alpha;
            dy_xi = r_1 * x[1] - x[0] * sigma3 * gay_alpha * (1 - gay_alpha2) *
                                     dr_dot_et * dt_dr_dot_dt * bun_alpha *
                                     bun_alpha;
            dt_xi = gay_alpha * sigma3 * (1 - gay_alpha2) * r2 * dr_dot_et *
                    dt_dr_dot_dt * bun_alpha * bun_alpha;
            epsilon1 = 1 - gay_alpha2 * r2_1 * dt_dr_dot_dt * dt_dr_dot_dt;
            epsilon1_N_1 = fast_power(epsilon1, gay_N - 1);
            epsilon13 = epsilon1 * epsilon1 * epsilon1;
            dx_epsilon1 = epsilon13 * gay_alpha2 * r2_1 *
                          (sini - x[0] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dy_epsilon1 = -epsilon13 * gay_alpha2 * r2_1 *
                          (cosi + x[1] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dt_epsilon1 = -epsilon13 * gay_alpha2 * dt_dr_dot_dt * dr_dot_et;
            epsilon2 = 1 - gay_alpha_p * dr_dot_et * dr_dot_et * bun_alpha_p;
            epsilon2_M_1 = fast_power(epsilon2, gay_M - 1);
            dr_epsilon2 = gay_alpha_p * (1 - gay_alpha_p2) * dr_dot_et *
                          dt_dr_dot_dt * bun_alpha_p *
                          bun_alpha_p;  // xの時y,yの時-xをかける;
            dt_epsilon2 = 2 * gay_alpha_p * r2 * dr_dot_et * dt_dr_dot_dt *
                          (1 - gay_alpha_p2);
            f[0] -= (dUWCA * dx_xi * epsilon1 * epsilon2 +
                     UWCA * (dx_epsilon1 * gay_N * epsilon2 +
                             x[1] * dr_epsilon2 * epsilon1 * gay_M)) *
                    epsilon1_N_1 * epsilon2_M_1;
            f[1] -= (dUWCA * dy_xi * epsilon1 * epsilon2 +
                     UWCA * (dy_epsilon1 * gay_N * epsilon2 -
                             x[0] * dr_epsilon2 * epsilon1 * gay_M)) *
                    epsilon1_N_1 * epsilon2_M_1;
            *f_theta -= (dUWCA * dt_xi * epsilon1 * epsilon2 +
                         UWCA * (dt_epsilon1 * gay_N * epsilon2 +
                                 dt_epsilon2 * epsilon1 * gay_M)) *
                        epsilon1_N_1 * epsilon2_M_1;
        }
    }
    x[0] += sa;
}
#endif
