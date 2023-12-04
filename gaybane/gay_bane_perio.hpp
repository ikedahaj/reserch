// Nn,dim=2d,gay_kapaが必要:2次元要のポテンシャル計算です。3次元に使わないでください;
#ifndef GAY_BANE_PERIO
#define GAY_BANE_PERIO

#define gay_kappa_p 5.
#define gay_M       2  // intでなくなるときはこーどを変更すること;
#define gay_N       1
#ifndef gay_kappa
#define gay_kappa 2
#endif
// baseのpower乗をかえす関数,power>=0;
inline double fast_power(double base, int power) {
    // if(power<0)base=1./base;
    double result = 1;
    while (power > 0) {
        if (power & 1) {
            // Can also use (power & 1) to make code even faster
            result = (result * base);
        }
        base = (base * base);
        power >>= 1;
        // Can also use power >>= 1; to make code even faster
    }
    return result;
}
constexpr double fast_power_cone(double base, int power) {
    double result = 1;
    while (power > 0) {
        if (power % 2 == 1) {
            // Can also use (power & 1) to make code even faster
            result = (result * base);
        }
        base = (base * base);
        power = power / 2;
        // Can also use power >>= 1; to make code even faster
    }
    return result;
}
// xのt乗根を返す;
constexpr double nth_root(double x, int t) {
    double x0 = x, ans = (x > 1) ? x : 1;
    if (t <= 0 || x < 0) return 1;
    for (int i = 0; i < 1000; i++) {
        x0 = ans;
        ans = x0 -
              (fast_power_cone(x0, t) - x) / (t * fast_power_cone(x0, t - 1));
        if (ans > x0) break;
    }
    return ans;
}
void ini_array(double (*x)[dim]) {
    for (int i = 0; i < Np; ++i)
        for (int j = 0; j < dim; ++j) x[i][j] = 0.0;
}
void ini_force(double (*f)[dim], double *f_theta) {
    for (int i = 0; i < Np; i++) {
        f[i][0] = 0;
        f[i][1] = 0.;
        f_theta[i] = 0.;
    }
}
static constexpr double Ly_inv = 1. / Ly;
static constexpr double Lx_inv = 1. / Lx;

inline double pri_fce_x(double x) {
    constexpr double Lx_2 = Lx / 2.;
    x -= Lx * floor((x + Lx_2) * Lx_inv);
    return x;
}
inline double pri_fce_y(double x) {
    constexpr double Ly_2 = Ly / 2.;
    x -= Ly * floor((x + Ly_2) * Ly_inv);
    return x;
}
void calc_force_gay_bane(double (*x)[dim2], double (*f)[dim], int (*list)[Nn],
                         double *theta, double *f_theta) {
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) /
                                 (gay_kappa * gay_kappa + 1.),
                     gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) /
                                   (nth_root(gay_kappa_p, gay_M) + 1),
                     cutf2 = cut * cut * gay_kappa * gay_kappa;
    // kがi,lがj;
    double dx, dy, dr2, cosk, sink, cosl, sinl, cos_kl_plus, cos_kl_minus,
        sin_kl_plus, sin_kl_minus, f_alpha, f_alpha_p, dr_dot_ek, dr_dot_el,
        sigma_hat, dr, xi, xi2, xi6, UWCA, dUWCA, dr_1, dr2_1, dx_xi, dy_xi,
        d_theta_xi_k, d_theta_xi_l, d_theta_epsilon1, epsilon1,
        epsilon1_N_1 /*epsilon1のN-1乗*/, dx_epsilon2, dy_epsilon2,
        dtheta_epsilon2_k, dtheta_epsilon2_l, epsilon2,
        epsilon2_M_1 /*epsilon2のM-1乗*/, Fx, Fy;
    ini_force(f, f_theta);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = pri_fce_x(x[i][0] - x[list[i][j]][0]);
            dy = pri_fce_y(x[i][1] - x[list[i][j]][1]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < cutf2) {  // 長径に合わせたカットおふ;
                cosk = cos(theta[i]);
                sink = sin(theta[i]);
                cosl = cos(theta[list[i][j]]);
                sinl = sin(theta[list[i][j]]);
                cos_kl_minus = cosk * cosl + sink * sinl;
                dr_dot_ek = dx * cosk + dy * sink;
                dr_dot_el = dx * cosl + dy * sinl;
                dr = sqrt(dr2);
                dr_1 = 1. / dr;
                dr2_1 = dr_1 * dr_1;
                f_alpha = (dr_dot_ek + dr_dot_el) * (dr_dot_ek + dr_dot_el) /
                              (1 + gay_alpha * cos_kl_minus) +
                          (dr_dot_ek - dr_dot_el) * (dr_dot_ek - dr_dot_el) /
                              (1 - gay_alpha * cos_kl_minus);
                sigma_hat = 1. / sqrt(1 - gay_alpha * 0.5 * dr2_1 * f_alpha);
                xi = dr - sigma_hat + 1;
                if (xi < cut) {
                    xi2 = 1. / (xi * xi);
                    xi6 = xi2 * xi2 * xi2;
                    UWCA = 4. * xi6 * (xi6 - 1.);
                    dUWCA = 48. * xi6 * (0.5 - xi6) / xi;
                    sin_kl_minus = sink * cosl - cosk * sinl;
                    sin_kl_plus = sink * cosl + cosk * sinl;
                    dx_xi = dx * dr_1 +
                            sigma_hat * sigma_hat * sigma_hat * gay_alpha *
                                dr2_1 *
                                (-dx * dr2_1 * f_alpha +
                                 (dx * (cosk + cosl) * (cosk + cosl) +
                                  dy * sin_kl_plus * (1 + cos_kl_minus)) /
                                     (1 + gay_alpha * cos_kl_minus) +
                                 (dx * (cosk - cosl) * (cosk - cosl) +
                                  dy * sin_kl_plus * (cos_kl_minus - 1)) /
                                     (1 - gay_alpha * cos_kl_minus));
                    dy_xi = dy * dr_1 +
                            sigma_hat * sigma_hat * sigma_hat * gay_alpha *
                                dr2_1 *
                                (-dy * dr2_1 * f_alpha +
                                 (dy * (sink + sinl) * (sink + sinl) +
                                  dx * sin_kl_plus * (1 + cos_kl_minus)) /
                                     (1 + gay_alpha * cos_kl_minus) +
                                 (dy * (sink - sinl) * (sink - sinl) +
                                  dx * sin_kl_plus * (cos_kl_minus - 1)) /
                                     (1 - gay_alpha * cos_kl_minus));

                    d_theta_xi_k =
                        0.5 * gay_alpha * dr2_1 *
                        ((2 * (dr_dot_ek + dr_dot_el) *
                              (-dx * sink + dy * cosk) *
                              (1 + gay_alpha * cos_kl_minus) +
                          gay_alpha * sin_kl_minus * (dr_dot_ek + dr_dot_el) *
                              (dr_dot_ek + dr_dot_el)) /
                             ((1 + gay_alpha * cos_kl_minus) *
                              (1 + gay_alpha * cos_kl_minus)) +
                         (2 * (dr_dot_ek - dr_dot_el) *
                              (dx * sink - dy * cosk) *
                              (1 - gay_alpha * cos_kl_minus) -
                          gay_alpha * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha * cos_kl_minus) *
                              (1 - gay_alpha * cos_kl_minus)));
                    d_theta_xi_l =
                        0.5 * gay_alpha * dr2_1 *
                        ((2 * (dr_dot_ek + dr_dot_el) *
                              (-dx * sinl + dy * cosl) *
                              (1 + gay_alpha * cos_kl_minus) -
                          gay_alpha * sin_kl_minus * (dr_dot_ek + dr_dot_el) *
                              (dr_dot_ek + dr_dot_el)) /
                             ((1 + gay_alpha * cos_kl_minus) *
                              (1 + gay_alpha * cos_kl_minus)) +
                         (2 * (dr_dot_ek - dr_dot_el) *
                              (dx * sinl - dy * cosl) *
                              (1 - gay_alpha * cos_kl_minus) +
                          gay_alpha * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha * cos_kl_minus) *
                              (1 - gay_alpha * cos_kl_minus)));
                    epsilon1 = 1. / sqrt(1 - gay_alpha * gay_alpha *
                                                 cos_kl_minus * cos_kl_minus);
                    epsilon1_N_1 = fast_power(epsilon1, gay_N - 1);
                    d_theta_epsilon1 = -gay_alpha * gay_alpha * epsilon1 *
                                       epsilon1 * epsilon1 * cos_kl_minus *
                                       sin_kl_minus;
                    f_alpha_p =
                        (dr_dot_ek + dr_dot_el) * (dr_dot_ek + dr_dot_el) /
                            (1 + gay_alpha_p * cos_kl_minus) +
                        (dr_dot_ek - dr_dot_el) * (dr_dot_ek - dr_dot_el) /
                            (1 - gay_alpha_p * cos_kl_minus);
                    epsilon2 = 1 - gay_alpha_p * 0.5 * dr2_1 * f_alpha_p;
                    epsilon2_M_1 = fast_power(epsilon2, gay_M - 1);
                    dx_epsilon2 = -gay_alpha_p * dr2_1 *
                                  (-dx * dr2_1 * f_alpha_p +
                                   (dx * (cosk + cosl) * (cosk + cosl) +
                                    dy * sin_kl_plus * (1 + cos_kl_minus)) /
                                       (1 + gay_alpha_p * cos_kl_minus) +
                                   (dx * (cosk - cosl) * (cosk - cosl) +
                                    dy * sin_kl_plus * (cos_kl_minus - 1)) /
                                       (1 - gay_alpha_p * cos_kl_minus));
                    dy_epsilon2 = -gay_alpha_p * dr2_1 *
                                  (-dy * dr2_1 * f_alpha_p +
                                   (dy * (sink + sinl) * (sink + sinl) +
                                    dx * sin_kl_plus * (1 + cos_kl_minus)) /
                                       (1 + gay_alpha_p * cos_kl_minus) +
                                   (dy * (sink - sinl) * (sink - sinl) +
                                    dx * sin_kl_plus * (cos_kl_minus - 1)) /
                                       (1 - gay_alpha_p * cos_kl_minus));
                    dtheta_epsilon2_k =
                        -0.5 * gay_alpha_p * dr2_1 *
                        ((2 * (dr_dot_ek + dr_dot_el) *
                              (-dx * sink + dy * cosk) *
                              (1 + gay_alpha_p * cos_kl_minus) +
                          gay_alpha_p * sin_kl_minus * (dr_dot_ek + dr_dot_el) *
                              (dr_dot_ek + dr_dot_el)) /
                             ((1 + gay_alpha_p * cos_kl_minus) *
                              (1 + gay_alpha_p * cos_kl_minus)) +
                         (2 * (dr_dot_ek - dr_dot_el) *
                              (dx * sink - dy * cosk) *
                              (1 - gay_alpha_p * cos_kl_minus) -
                          gay_alpha_p * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha_p * cos_kl_minus) *
                              (1 - gay_alpha_p * cos_kl_minus)));
                    dtheta_epsilon2_l =
                        -0.5 * gay_alpha_p * dr2_1 *
                        ((2 * (dr_dot_ek + dr_dot_el) *
                              (-dx * sinl + dy * cosl) *
                              (1 + gay_alpha_p * cos_kl_minus) -
                          gay_alpha_p * sin_kl_minus * (dr_dot_ek + dr_dot_el) *
                              (dr_dot_ek + dr_dot_el)) /
                             ((1 + gay_alpha_p * cos_kl_minus) *
                              (1 + gay_alpha_p * cos_kl_minus)) +
                         (2 * (dr_dot_ek - dr_dot_el) *
                              (dx * sinl - dy * cosl) *
                              (1 - gay_alpha_p * cos_kl_minus) +
                          gay_alpha_p * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha_p * cos_kl_minus) *
                              (1 - gay_alpha_p * cos_kl_minus)));
                    Fx = (dUWCA * dx_xi * epsilon2 +
                          gay_M * dx_epsilon2 * UWCA) *
                         epsilon1_N_1 * epsilon1 * epsilon2_M_1;
                    f[i][0] -= Fx;
                    f[list[i][j]][0] += Fx;
                    Fy = (dUWCA * dy_xi * epsilon2 +
                          gay_M * dy_epsilon2 * UWCA) *
                         epsilon1_N_1 * epsilon1 * epsilon2_M_1;
                    f[i][1] -= Fy;
                    f[list[i][j]][1] += Fy;
                    // double f_epsi_x =
                    //            (dUWCA * d_theta_xi_k * epsilon1 * epsilon2 +
                    //             (gay_M * dtheta_epsilon2_k * epsilon1 +
                    //              gay_N * d_theta_epsilon1 * epsilon2) *
                    //                 UWCA) *
                    //            epsilon1_N_1 * epsilon2_M_1,
                    //        f_epsi_l =
                    //            (dUWCA * d_theta_xi_l * epsilon1 * epsilon2 +
                    //             (gay_M * dtheta_epsilon2_l * epsilon1 -
                    //              gay_N * d_theta_epsilon1 * epsilon2) *
                    //                 UWCA) *
                    //            epsilon1_N_1 * epsilon2_M_1;

                    f_theta[i] -= (dUWCA * d_theta_xi_k * epsilon1 * epsilon2 +
                                   (gay_M * dtheta_epsilon2_k * epsilon1 +
                                    gay_N * d_theta_epsilon1 * epsilon2) *
                                       UWCA) *
                                  epsilon1_N_1 * epsilon2_M_1;
                    f_theta[list[i][j]] -=
                        (dUWCA * d_theta_xi_l * epsilon1 * epsilon2 +
                         (gay_M * dtheta_epsilon2_l * epsilon1 +
                          gay_N * d_theta_epsilon1 * epsilon2) *
                             UWCA) *
                        epsilon1_N_1 * epsilon2_M_1;
                }
            }
        }
}

static constexpr int Mx = Lx / (cut * gay_kappa + skin);
static constexpr int My = Ly / (cut * gay_kappa + skin);
inline int           peri_cell_x(int m) {
    if (m < 0)
        return m + Mx;
    else if (m >= Mx)
        return m - Mx;
    else
        return m;
}
inline int peri_cell_y(int m) {
    if (m < 0)
        return m + My;
    else if (m >= My)
        return m - My;
    else
        return m;
}
void cell_list(int (*list)[Nn], double (*x)[dim2]) {
    int              map_index, nx[Np][dim];
    constexpr int    m2 = Mx * My;
    constexpr double thresh2 =
                         (cut * gay_kappa + skin) * (cut * gay_kappa + skin),
                     bitx = Mx / Lx, bity = My / Ly;
    double dx, dy;
    // int(*map)[Np + 1] = new int[m2][Np + 1];
    std::vector<std::deque<int>> map(m2);

    for (int i = 0; i < Np; ++i) {
        nx[i][0] = (int) ((x[i][0]) * bitx);
        nx[i][1] = (int) ((x[i][1]) * bity);
        for (int m = nx[i][1] - 1; m <= nx[i][1] + 1; ++m) {
            for (int l = nx[i][0] - 1; l <= nx[i][0] + 1; ++l) {
                map_index = peri_cell_x(l) + My * peri_cell_y(m);
                map[map_index].emplace_front(i);
            }
        }
    }
    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index = nx[i][0] + My * nx[i][1];
        for (auto &j : map[map_index]) {
            if (j > i) {
                dx = pri_fce_x(x[i][0] - x[j][0]);
                dy = pri_fce_y(x[i][1] - x[j][1]);
                if ((dx * dx + dy * dy) < thresh2) {
                    list[i][0]++;
                    list[i][list[i][0]] = j;
                }
            }
        }
    }
    // delete[] map;
}
void ver_list(int (*list)[Nn], double (*x)[dim2]) {
    double           dx, dy, dr2;
    constexpr double thresh2 =
        (cut * gay_kappa + skin) * (cut * gay_kappa + skin);
    for (int i = 0; i < Np; i++) {
        list[i][0] = 0;
    }

    for (int i = 0; i < Np; i++) {
        for (int j = i + 1; j < Np; j++) {
            dx = pri_fce_x(x[i][0] - x[j][0]);
            dy = pri_fce_y(x[i][1] - x[j][1]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < thresh2) {
                list[i][0]++;
                list[i][list[i][0]] = j;
            }
        }
    }
}
void update(double (*x_update)[dim], double (*x)[dim2]) {
    for (int i = 0; i < Np; i++)
        for (int j = 0; j < dim; j++) x_update[i][j] = x[i][j];
}

void calc_disp_max(double *disp_max, double (*x)[dim2],
                   double (*x_update)[dim]) {
    double dx, dy;
    double disp;
    for (int i = 0; i < Np; i++) {
        dx = pri_fce_x(x[i][0] - x_update[i][0]);
        dy = pri_fce_y(x[i][1] - x_update[i][1]);
        disp = dx * dx + dy * dy;
        if (disp > *disp_max) *disp_max = disp;
    }
}

void auto_list_update(double (*x)[dim2], double (*x_update)[dim],
                      int (*list)[Nn]) {
    // static int count = 0;
    // count++;
    static double    disp_max = skin * skin+100;
    constexpr double skin2 = skin * skin * 0.25, skinini = skin2 * 0.9;
    calc_disp_max(&disp_max, x, x_update);
    if (disp_max >= skin2) {
        ver_list(list, x);
        update(x_update, x);
        //    std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
        disp_max = skinini;
        // count = 0;
    }
}

#endif