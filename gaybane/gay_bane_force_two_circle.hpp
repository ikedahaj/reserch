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
        base  = (base * base);
        power = power / 2;
        // Can also use power >>= 1; to make code even faster
    }
    return result;
}
// xのt乗根を返す;
constexpr double nth_root(double x, int t) {
    double x0 = x, ans = (x > 1) ? x : 1;
    if (t == 0 || x < 0)
        return 1;
    else if (t < 0)
        return 1. / nth_root(x, -t);
    for (int i = 0; i < 1000; i++) {
        x0  = ans;
        ans = x0 - (fast_power_cone(x0, t) - x) / (t * fast_power_cone(x0, t - 1));
        if (ans > x0)
            break;
    }
    return ans;
}

void ini_force(double (*f)[dim], double *f_theta) {
    for (int i = 0; i < Np; i++) {
        f[i][0]    = 0;
        f[i][1]    = 0.;
        f_theta[i] = 0.;
    }
}

/// @brief
/// gay-baneポテンシャルの力の2次元計算.周期境界条件を考慮していない.グローバル変数としてdim,Nn,cutが必要;
/// @param x 粒子座標 Np行dim列配列;
/// @param f 力      同上;
/// @param list 帳簿法のリスト　Np行Nn列配列;
/// @param theta 粒子角度　Np列配列;
/// @param f_theta 角度に掛かる力　同上;
void calc_force_gay_bane(double (*x)[dim], double (*f)[dim], int (*list)[Nn], double *theta, double *f_theta) {
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) / (gay_kappa * gay_kappa + 1.), gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) / (nth_root(gay_kappa_p, gay_M) + 1),
                     cutf2 = (cut + gay_kappa - 1) * (cut + gay_kappa - 1), gay_alpha2 = gay_alpha * gay_alpha, gay_alpha_p2 = gay_alpha_p * gay_alpha_p;
    // kがi,lがj;
    double dx, dy, dr2, cosk, sink, cosj, sinj, cos_kl_plus, cos_kl_minus, cos_kl_minus2, sin_kl_plus, sin_kl_minus, dx_f_bunsi, dy_f_bunsi, f_alpha, f_bun_alpha, f_bun_alpha_p, f_syou_alpha_p,
        dr_dot_ek, dr_dot_el, dk_dr_dot_ek, dl_dr_dot_el, sigma_hat, sigma_hat3, dr, xi, xi_1, xi2, xi6, UWCA, dUWCA, dr_1, dr2_1, dx_xi, dy_xi, d_theta_xi_k, d_theta_xi_l, d_theta_epsilon1, epsilon1,
        epsilon1_N_1 /*epsilon1のN-1乗*/, dx_epsilon2, dy_epsilon2, d_theta_epsilon2_k, d_theta_epsilon2_l, epsilon2, epsilon2_M_1 /*epsilon2のM-1乗*/, Fx, Fy;
    ini_force(f, f_theta);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx  = x[i][0] - x[list[i][j]][0];
            dy  = x[i][1] - x[list[i][j]][1];
            dr2 = dx * dx + dy * dy;
            if (dr2 < cutf2) {
                cosk          = cos(theta[i]);
                sink          = sin(theta[i]);
                cosj          = cos(theta[list[i][j]]);
                sinj          = sin(theta[list[i][j]]);
                dr_dot_ek     = dx * cosk + dy * sink;
                dr_dot_el     = dx * cosj + dy * sinj;
                cos_kl_minus  = cosk * cosj + sink * sinj;
                cos_kl_minus2 = cos_kl_minus * cos_kl_minus;
                sin_kl_plus   = cosk * sinj + sink * cosj;
                f_bun_alpha   = 1. / (1 - gay_alpha2 * cos_kl_minus2);
                f_alpha       = 2 * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el - 2. * gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_el)) * f_bun_alpha;
                dr            = sqrt(dr2);
                dr_1          = 1. / dr;
                dr2_1         = dr_1 * dr_1;
                sigma_hat     = 1. / sqrt(1. - gay_alpha * 0.5 * dr2_1 * f_alpha);
                xi            = dr - sigma_hat + 1.;
                if (xi < cut) {
                    xi_1         = 1. / xi;
                    xi2          = xi_1 * xi_1;
                    xi6          = xi2 * xi2 * xi2;
                    UWCA         = 4. * xi6 * (xi6 - 1.);
                    dUWCA        = 48. * xi6 * (0.5 - xi6) * xi_1;
                    dk_dr_dot_ek = -dx * sink + dy * cosk;
                    dl_dr_dot_el = -dx * sinj + dy * cosj;
                    sin_kl_minus = sink * cosj - cosk * sinj;
                    sigma_hat3   = sigma_hat * sigma_hat * sigma_hat;
                    dx_xi =
                        dx * dr_1 +
                        sigma_hat3 * gay_alpha * dr2_1 *
                            (dx * dr2_1 * f_alpha * 0.5 -
                             gay_alpha * ((dx * (cosk * cosk + cosj * cosj) + dy * (cosk * sink + cosj * sinj)) - (2. * dx * cosk * cosj + dy * sin_kl_plus) * cos_kl_minus * gay_alpha) * f_bun_alpha);
                    dy_xi =
                        dy * dr_1 +
                        sigma_hat3 * gay_alpha * dr2_1 *
                            (dy * dr2_1 * f_alpha * 0.5 -
                             gay_alpha * ((dy * (sink * sink + sinj * sinj) + dx * (cosk * sink + cosj * sinj)) - (2. * dy * sink * sinj + dx * sin_kl_plus) * cos_kl_minus * gay_alpha) * f_bun_alpha);
                    d_theta_xi_k =
                        -gay_alpha * sigma_hat3 * dr2_1 *
                        ((dr_dot_ek - dr_dot_el * gay_alpha * cos_kl_minus) * dk_dr_dot_ek * (1. - gay_alpha2 * cos_kl_minus2) +
                         gay_alpha * sin_kl_minus * ((1 + gay_alpha2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                        f_bun_alpha * f_bun_alpha;
                    d_theta_xi_l =
                        -gay_alpha * sigma_hat3 * dr2_1 *
                        ((dr_dot_el - dr_dot_ek * gay_alpha * cos_kl_minus) * dl_dr_dot_el * (1. - gay_alpha2 * cos_kl_minus2) -
                         gay_alpha * sin_kl_minus * ((1 + gay_alpha2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                        f_bun_alpha * f_bun_alpha;
                    epsilon1         = sqrt(f_bun_alpha);  // 内容を変えるとき注意;
                    epsilon1_N_1     = fast_power(epsilon1, gay_N - 1);
                    d_theta_epsilon1 = -epsilon1 * epsilon1 * epsilon1 * gay_alpha2 * cos_kl_minus * sin_kl_minus;
                    f_bun_alpha_p    = 1. / (1 - gay_alpha_p2 * cos_kl_minus2);
                    f_syou_alpha_p   = 2. * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el - 2. * gay_alpha_p * cos_kl_minus * dr_dot_ek * dr_dot_el);
                    epsilon2         = 1 - gay_alpha_p * dr2_1 * 0.5 * f_syou_alpha_p * f_bun_alpha_p;
                    epsilon2_M_1     = fast_power(epsilon2, gay_M - 1);
                    dx_epsilon2      = dr2_1 * gay_alpha_p *
                                  (dr2_1 * dx * f_syou_alpha_p -
                                   2. * (dx * (cosk * cosk + cosj * cosj) + dy * (cosk * sink + cosj * sinj) - (dx * cosk * cosj * 2 + dy * sin_kl_plus) * gay_alpha_p * cos_kl_minus)) *
                                  f_bun_alpha_p;
                    dy_epsilon2 = dr2_1 * gay_alpha_p *
                                  (dr2_1 * dy * f_syou_alpha_p -
                                   2 * (dy * (sink * sink + sinj * sinj) + dx * (cosk * sink + cosj * sinj) - (dy * sink * sinj * 2 + dx * sin_kl_plus) * gay_alpha_p * cos_kl_minus)) *
                                  f_bun_alpha_p;
                    d_theta_epsilon2_k =
                        -2. * gay_alpha_p * dr2_1 *
                        ((dr_dot_ek - dr_dot_el * gay_alpha_p * cos_kl_minus) * dk_dr_dot_ek * (1 - gay_alpha_p2 * cos_kl_minus2) +
                         gay_alpha_p * sin_kl_minus * ((1 + gay_alpha_p2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha_p * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                        f_bun_alpha_p * f_bun_alpha_p;
                    d_theta_epsilon2_l =
                        -2. * gay_alpha_p * dr2_1 *
                        ((dr_dot_el - dr_dot_ek * gay_alpha_p * cos_kl_minus) * dl_dr_dot_el * (1 - gay_alpha_p2 * cos_kl_minus2) -
                         gay_alpha_p * sin_kl_minus * ((1 + gay_alpha_p2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha_p * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                        f_bun_alpha_p * f_bun_alpha_p;
                    Fx = (dUWCA * dx_xi * epsilon2 + gay_M * dx_epsilon2 * UWCA) * epsilon1_N_1 * epsilon1 * epsilon2_M_1;
                    f[i][0] -= Fx;
                    f[list[i][j]][0] += Fx;
                    Fy = (dUWCA * dy_xi * epsilon2 + gay_M * dy_epsilon2 * UWCA) * epsilon1_N_1 * epsilon1 * epsilon2_M_1;
                    f[i][1] -= Fy;
                    f[list[i][j]][1] += Fy;
                    f_theta[i] -= (dUWCA * d_theta_xi_k * epsilon1 * epsilon2 + (gay_M * d_theta_epsilon2_k * epsilon1 + gay_N * d_theta_epsilon1 * epsilon2) * UWCA) * epsilon1_N_1 * epsilon2_M_1;
                    f_theta[list[i][j]] -=
                        (dUWCA * d_theta_xi_l * epsilon1 * epsilon2 + (gay_M * d_theta_epsilon2_l * epsilon1 - gay_N * d_theta_epsilon1 * epsilon2) * UWCA) * epsilon1_N_1 * epsilon2_M_1;
                }
            }
        }
}
inline void calc_force_wall_between_wall_force_part(double *x, double *f, double theta, double *f_theta) {
    constexpr double R_P = R + 0.5, y0 = (-gay_kappa * R * Rbit - 2 * R_P * fast_power_cone(gay_kappa * gay_kappa - R * R * Rbit * Rbit + 4 * R_P * R_P, 2)) / (gay_kappa * gay_kappa + 4 * R_P * R_P),
                     x0        = -gay_kappa * y0 / (2 * R_P);
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) / (gay_kappa * gay_kappa + 1.), gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) / (nth_root(gay_kappa_p, gay_M) + 1),
                     cutf2 = (cut + gay_kappa - 1) * (cut + gay_kappa - 1), gay_alpha2 = gay_alpha * gay_alpha, gay_alpha_p2 = gay_alpha_p * gay_alpha_p;
    constexpr double cosj = -y0 / R_P, sinj = x0 / R_P;
    // kがi,lがj;
    double dx, dy, dr2, cosk, sink, cos_kl_plus, cos_kl_minus, cos_kl_minus2, sin_kl_plus, sin_kl_minus, dx_f_bunsi, dy_f_bunsi, f_alpha, f_bun_alpha, f_bun_alpha_p, f_syou_alpha_p, dr_dot_ek,
        dr_dot_el, dk_dr_dot_ek, dl_dr_dot_el, sigma_hat, sigma_hat3, dr, xi, xi_1, xi2, xi6, UWCA, dUWCA, dr_1, dr2_1, dx_xi, dy_xi, d_theta_xi_k, d_theta_xi_l, d_theta_epsilon1, epsilon1,
        epsilon1_N_1 /*epsilon1のN-1乗*/, dx_epsilon2, dy_epsilon2, d_theta_epsilon2_k, d_theta_epsilon2_l, epsilon2, epsilon2_M_1 /*epsilon2のM-1乗*/, Fx, Fy;

    dx  = x[0] - x0;
    dy  = x[1] - y0;
    dr2 = dx * dx + dy * dy;
    if (dr2 < cutf2) {
        cosk          = cos(theta);
        sink          = sin(theta);
        dr_dot_ek     = dx * cosk + dy * sink;
        dr_dot_el     = dx * cosj + dy * sinj;
        cos_kl_minus  = cosk * cosj + sink * sinj;
        cos_kl_minus2 = cos_kl_minus * cos_kl_minus;
        sin_kl_plus   = cosk * sinj + sink * cosj;
        f_bun_alpha   = 1. / (1 - gay_alpha2 * cos_kl_minus2);
        f_alpha       = 2 * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el - 2. * gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_el)) * f_bun_alpha;
        dr            = sqrt(dr2);
        dr_1          = 1. / dr;
        dr2_1         = dr_1 * dr_1;
        sigma_hat     = 1. / sqrt(1. - gay_alpha * 0.5 * dr2_1 * f_alpha);
        xi            = dr - sigma_hat + 1.;
        if (xi < cut) {
            xi_1         = 1. / xi;
            xi2          = xi_1 * xi_1;
            xi6          = xi2 * xi2 * xi2;
            UWCA         = 4. * xi6 * (xi6 - 1.);
            dUWCA        = 48. * xi6 * (0.5 - xi6) * xi_1;
            dk_dr_dot_ek = -dx * sink + dy * cosk;
            dl_dr_dot_el = -dx * sinj + dy * cosj;
            sin_kl_minus = sink * cosj - cosk * sinj;
            sigma_hat3   = sigma_hat * sigma_hat * sigma_hat;
            dx_xi        = dx * dr_1 +
                    sigma_hat3 * gay_alpha * dr2_1 *
                        (dx * dr2_1 * f_alpha * 0.5 -
                         gay_alpha * ((dx * (cosk * cosk + cosj * cosj) + dy * (cosk * sink + cosj * sinj)) - (2. * dx * cosk * cosj + dy * sin_kl_plus) * cos_kl_minus * gay_alpha) * f_bun_alpha);
            dy_xi = dy * dr_1 +
                    sigma_hat3 * gay_alpha * dr2_1 *
                        (dy * dr2_1 * f_alpha * 0.5 -
                         gay_alpha * ((dy * (sink * sink + sinj * sinj) + dx * (cosk * sink + cosj * sinj)) - (2. * dy * sink * sinj + dx * sin_kl_plus) * cos_kl_minus * gay_alpha) * f_bun_alpha);
            d_theta_xi_k = -gay_alpha * sigma_hat3 * dr2_1 *
                           ((dr_dot_ek - dr_dot_el * gay_alpha * cos_kl_minus) * dk_dr_dot_ek * (1. - gay_alpha2 * cos_kl_minus2) +
                            gay_alpha * sin_kl_minus * ((1 + gay_alpha2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                           f_bun_alpha * f_bun_alpha;
            d_theta_xi_l = -gay_alpha * sigma_hat3 * dr2_1 *
                           ((dr_dot_el - dr_dot_ek * gay_alpha * cos_kl_minus) * dl_dr_dot_el * (1. - gay_alpha2 * cos_kl_minus2) -
                            gay_alpha * sin_kl_minus * ((1 + gay_alpha2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                           f_bun_alpha * f_bun_alpha;
            epsilon1         = sqrt(f_bun_alpha);  // 内容を変えるとき注意;
            epsilon1_N_1     = fast_power(epsilon1, gay_N - 1);
            d_theta_epsilon1 = -epsilon1 * epsilon1 * epsilon1 * gay_alpha2 * cos_kl_minus * sin_kl_minus;
            f_bun_alpha_p    = 1. / (1 - gay_alpha_p2 * cos_kl_minus2);
            f_syou_alpha_p   = 2. * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el - 2. * gay_alpha_p * cos_kl_minus * dr_dot_ek * dr_dot_el);
            epsilon2         = 1 - gay_alpha_p * dr2_1 * 0.5 * f_syou_alpha_p * f_bun_alpha_p;
            epsilon2_M_1     = fast_power(epsilon2, gay_M - 1);
            dx_epsilon2 =
                dr2_1 * gay_alpha_p *
                (dr2_1 * dx * f_syou_alpha_p - 2. * (dx * (cosk * cosk + cosj * cosj) + dy * (cosk * sink + cosj * sinj) - (dx * cosk * cosj * 2 + dy * sin_kl_plus) * gay_alpha_p * cos_kl_minus)) *
                f_bun_alpha_p;
            dy_epsilon2 =
                dr2_1 * gay_alpha_p *
                (dr2_1 * dy * f_syou_alpha_p - 2 * (dy * (sink * sink + sinj * sinj) + dx * (cosk * sink + cosj * sinj) - (dy * sink * sinj * 2 + dx * sin_kl_plus) * gay_alpha_p * cos_kl_minus)) *
                f_bun_alpha_p;
            d_theta_epsilon2_k =
                -2. * gay_alpha_p * dr2_1 *
                ((dr_dot_ek - dr_dot_el * gay_alpha_p * cos_kl_minus) * dk_dr_dot_ek * (1 - gay_alpha_p2 * cos_kl_minus2) +
                 gay_alpha_p * sin_kl_minus * ((1 + gay_alpha_p2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha_p * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                f_bun_alpha_p * f_bun_alpha_p;
            d_theta_epsilon2_l =
                -2. * gay_alpha_p * dr2_1 *
                ((dr_dot_el - dr_dot_ek * gay_alpha_p * cos_kl_minus) * dl_dr_dot_el * (1 - gay_alpha_p2 * cos_kl_minus2) -
                 gay_alpha_p * sin_kl_minus * ((1 + gay_alpha_p2 * cos_kl_minus2) * dr_dot_ek * dr_dot_el - gay_alpha_p * cos_kl_minus * (dr_dot_ek * dr_dot_ek + dr_dot_el * dr_dot_el))) *
                f_bun_alpha_p * f_bun_alpha_p;
            Fx = (dUWCA * dx_xi * epsilon2 + gay_M * dx_epsilon2 * UWCA) * epsilon1_N_1 * epsilon1 * epsilon2_M_1;
            f[0] -= Fx;
            Fy = (dUWCA * dy_xi * epsilon2 + gay_M * dy_epsilon2 * UWCA) * epsilon1_N_1 * epsilon1 * epsilon2_M_1;
            f[1] -= Fy;
            *f_theta -= (dUWCA * d_theta_xi_k * epsilon1 * epsilon2 + (gay_M * d_theta_epsilon2_k * epsilon1 + gay_N * d_theta_epsilon1 * epsilon2) * UWCA) * epsilon1_N_1 * epsilon2_M_1;
        }
    }
}
inline void calc_force_wall_bitween_wall(double *x, double *f, double theta, double *f_theta) {
    double flag_minus = 1;
    if (x[1] > 0)
        flag_minus = -1;
    x[1] *= flag_minus;
    // 右側のダミー粒子との力;
    calc_force_wall_between_wall_force_part(x, f, theta, f_theta);
    x[0] *= -1;
    calc_force_wall_between_wall_force_part(x, f, theta, f_theta);
    x[0] *= -1;
    x[1] *= flag_minus;
}
/// @brief 壁とgay-bane粒子の相互作用.2次元用;
/// @param x 粒子座標。粒子iについての二次元ザ行を送ること;
/// @param f 粒子の力。同上;
/// @param theta 粒子の角度.値私で良い;
/// @param f_theta 確度にかかる力.(&f[i])などと渡すこと.
inline void calc_force_wall_gay_bane(double *x, double *f, double theta, double *f_theta) {
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) / (gay_kappa * gay_kappa + 1.), gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) / (nth_root(gay_kappa_p, gay_M) + 1),
                     cutwf2 = (R + 1.5 - cut - 1. / nth_root(1 - gay_alpha, 2)) * (R + 1.5 - cut - 1. / nth_root(1 - gay_alpha, 2)), gay_alpha2 = gay_alpha * gay_alpha,
                     gay_alpha_p2 = gay_alpha_p * gay_alpha_p;
    double r2, sa;
    if (x[0] > 0)
        sa = center_rignt;
    else
        sa = center_left;
    x[0] -= sa;
    r2 = x[0] * x[0] + x[1] * x[1];
    if (r2 > cutwf2) {
        double r, sini, cosi, dr_dot_et, dt_dr_dot_dt, bun_alpha, sigma, xi;
        r            = sqrt(r2);
        sini         = sin(theta);
        cosi         = cos(theta);
        dr_dot_et    = x[0] * cosi + x[1] * sini;
        dt_dr_dot_dt = x[0] * sini - x[1] * cosi;
        bun_alpha    = 1. / (r2 - gay_alpha2 * dt_dr_dot_dt * dt_dr_dot_dt);
        sigma        = 1. / sqrt(1 - gay_alpha * dr_dot_et * dr_dot_et * bun_alpha);
        xi           = R + 0.5 - r - sigma + 1.;
        if (xi < cut) {
            double sigma3, bun_alpha_p, r_1, r2_1, epsilon13, xi_1, xi2, xi6, UWCA, dUWCA, dx_xi, dy_xi, dt_xi, epsilon1, epsilon1_N_1, dx_epsilon1, dy_epsilon1, dt_epsilon1, epsilon2, epsilon2_M_1,
                dr_epsilon2, dt_epsilon2;
            sigma3       = sigma * sigma * sigma;
            bun_alpha_p  = 1. / (r2 - gay_alpha_p2 * dt_dr_dot_dt * dt_dr_dot_dt);
            r_1          = -1. / r;
            r2_1         = r_1 * r_1;
            xi_1         = 1. / xi;
            xi2          = xi_1 * xi_1;
            xi6          = xi2 * xi2 * xi2;
            UWCA         = 4. * xi6 * (xi6 - 1);
            dUWCA        = 48. * xi6 * (0.5 - xi6) * xi_1;
            dx_xi        = r_1 * x[0] + x[1] * sigma3 * gay_alpha * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            dy_xi        = r_1 * x[1] - x[0] * sigma3 * gay_alpha * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            dt_xi        = gay_alpha * sigma3 * r2 * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            epsilon1     = sqrt(bun_alpha) * r;
            epsilon1_N_1 = fast_power(epsilon1, gay_N - 1);
            epsilon13    = epsilon1 * epsilon1 * epsilon1;
            dx_epsilon1  = epsilon13 * gay_alpha2 * r2_1 * (sini - x[0] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dy_epsilon1  = -epsilon13 * gay_alpha2 * r2_1 * (cosi + x[1] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dt_epsilon1  = -epsilon13 * gay_alpha2 * dt_dr_dot_dt * dr_dot_et * r2_1;
            epsilon2     = 1 - gay_alpha_p * dr_dot_et * dr_dot_et * bun_alpha_p;
            epsilon2_M_1 = fast_power(epsilon2, gay_M - 1);
            dr_epsilon2  = gay_alpha_p * (-1 - gay_alpha_p2) * dr_dot_et * dt_dr_dot_dt * bun_alpha_p * bun_alpha_p;  // xの時y,yの時-xをかける;
            dt_epsilon2  = 2 * gay_alpha_p * dr_dot_et * dt_dr_dot_dt * r2 * (1 - gay_alpha_p2) * bun_alpha_p * bun_alpha_p;
            f[0] -= (dUWCA * dx_xi * epsilon1 * epsilon2 + UWCA * (dx_epsilon1 * gay_N * epsilon2 + x[1] * dr_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
            f[1] -= (dUWCA * dy_xi * epsilon1 * epsilon2 + UWCA * (dy_epsilon1 * gay_N * epsilon2 - x[0] * dr_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
            *f_theta -= (dUWCA * dt_xi * epsilon1 * epsilon2 + UWCA * (dt_epsilon1 * gay_N * epsilon2 + dt_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
        }
    }
    x[0] += sa;
}
/// @brief 壁とgay-bane粒子の相互作用.2次元用,R_INI!=Rの時用;
/// @param x 粒子座標。粒子iについての二次元ザ行を送ること;
/// @param f 粒子の力。同上;
/// @param theta 粒子の角度.値私で良い;
/// @param f_theta 確度にかかる力.(&f[i])などと渡すこと.
/// @param R_ini そのときの半径;
inline void calc_wall_gay_bane_anill(double *x, double *f, double theta, double *f_theta, double R_ini) {
    constexpr double gay_alpha = (gay_kappa * gay_kappa - 1.) / (gay_kappa * gay_kappa + 1.), gay_alpha_p = (nth_root(gay_kappa_p, gay_M) - 1) / (nth_root(gay_kappa_p, gay_M) + 1),

                     gay_alpha2 = gay_alpha * gay_alpha, gay_alpha_p2 = gay_alpha_p * gay_alpha_p;
    double r2, sa, cutwf2 = (R_ini + 1.5 - cut - 1. / nth_root(1 - gay_alpha, 2)) * (R_ini + 1.5 - cut - 1. / nth_root(1 - gay_alpha, 2));
    if (x[0] > 0)
        sa = R_ini * rbit_2;
    else
        sa = -R_ini * rbit_2;
    x[0] -= sa;
    r2 = x[0] * x[0] + x[1] * x[1];
    if (r2 > cutwf2) {
        double r, sini, cosi, dr_dot_et, dt_dr_dot_dt, bun_alpha, sigma, xi;
        r            = sqrt(r2);
        sini         = sin(theta);
        cosi         = cos(theta);
        dr_dot_et    = x[0] * cosi + x[1] * sini;
        dt_dr_dot_dt = x[0] * sini - x[1] * cosi;
        bun_alpha    = 1. / (r2 - gay_alpha2 * dt_dr_dot_dt * dt_dr_dot_dt);
        sigma        = 1. / sqrt(1 - gay_alpha * dr_dot_et * dr_dot_et * bun_alpha);
        xi           = R_ini + 0.5 - r - sigma + 1.;
        if (xi < cut) {
            double sigma3, bun_alpha_p, r_1, r2_1, epsilon13, xi_1, xi2, xi6, UWCA, dUWCA, dx_xi, dy_xi, dt_xi, epsilon1, epsilon1_N_1, dx_epsilon1, dy_epsilon1, dt_epsilon1, epsilon2, epsilon2_M_1,
                dr_epsilon2, dt_epsilon2;
            sigma3       = sigma * sigma * sigma;
            bun_alpha_p  = 1. / (r2 - gay_alpha_p2 * dt_dr_dot_dt * dt_dr_dot_dt);
            r_1          = -1. / r;
            r2_1         = r_1 * r_1;
            xi_1         = 1. / xi;
            xi2          = xi_1 * xi_1;
            xi6          = xi2 * xi2 * xi2;
            UWCA         = 4. * xi6 * (xi6 - 1);
            dUWCA        = 48. * xi6 * (0.5 - xi6) * xi_1;
            dx_xi        = r_1 * x[0] + x[1] * sigma3 * gay_alpha * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            dy_xi        = r_1 * x[1] - x[0] * sigma3 * gay_alpha * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            dt_xi        = gay_alpha * sigma3 * r2 * (1 - gay_alpha2) * dr_dot_et * dt_dr_dot_dt * bun_alpha * bun_alpha;
            epsilon1     = sqrt(bun_alpha) * r;
            epsilon1_N_1 = fast_power(epsilon1, gay_N - 1);
            epsilon13    = epsilon1 * epsilon1 * epsilon1;
            dx_epsilon1  = epsilon13 * gay_alpha2 * r2_1 * (sini - x[0] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dy_epsilon1  = -epsilon13 * gay_alpha2 * r2_1 * (cosi + x[1] * r2_1 * dt_dr_dot_dt) * dt_dr_dot_dt;
            dt_epsilon1  = -epsilon13 * gay_alpha2 * dt_dr_dot_dt * dr_dot_et * r2_1;
            epsilon2     = 1 - gay_alpha_p * dr_dot_et * dr_dot_et * bun_alpha_p;
            epsilon2_M_1 = fast_power(epsilon2, gay_M - 1);
            dr_epsilon2  = gay_alpha_p * (-1 - gay_alpha_p2) * dr_dot_et * dt_dr_dot_dt * bun_alpha_p * bun_alpha_p;  // xの時y,yの時-xをかける;
            dt_epsilon2  = 2 * gay_alpha_p * dr_dot_et * dt_dr_dot_dt * r2 * (1 - gay_alpha_p2) * bun_alpha_p * bun_alpha_p;
            f[0] -= (dUWCA * dx_xi * epsilon1 * epsilon2 + UWCA * (dx_epsilon1 * gay_N * epsilon2 + x[1] * dr_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
            f[1] -= (dUWCA * dy_xi * epsilon1 * epsilon2 + UWCA * (dy_epsilon1 * gay_N * epsilon2 - x[0] * dr_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
            *f_theta -= (dUWCA * dt_xi * epsilon1 * epsilon2 + UWCA * (dt_epsilon1 * gay_N * epsilon2 + dt_epsilon2 * epsilon1 * gay_M)) * epsilon1_N_1 * epsilon2_M_1;
        }
    }
    x[0] += sa;
}
void cell_list(int (*list)[Nn], double (*x)[dim]) {
    int              map_index, nx[dim];
    constexpr double threash2 = (cut + gay_kappa - 1 + skin) * (cut + gay_kappa - 1 + skin);
    constexpr double xlen_2   = (2. * R + Rbit * R) / 2.;
    constexpr int    Mx       = (int) (xlen_2 * 2. / (cut + gay_kappa - 1 + skin));
    constexpr int    My       = (int) (2. * R / (cut + gay_kappa - 1 + skin));  // M<=2R/(cutmax+skin)
    constexpr int    m2       = Mx * My;
    constexpr double R2 = 2. * R, bitx = Mx / (xlen_2 * 2.),
                     bity = My / (R2);  // ひとつのせるの幅の逆数;
    double dx, dy;
    // int(*map)[Np + 1] = new int[m2][Np + 1];
    std::vector<std::vector<int>> map(m2);
    for (int i = 0; i < m2; ++i) {
        map[i].reserve(20);
    }

    for (int i = 0; i < Np; ++i) {
        nx[0] = (int) ((x[i][0] + xlen_2) * bitx);
        nx[1] = (int) ((x[i][1] + R) * bity);
        for (int m = max(nx[1] - 1, 0), mm = min(nx[1] + 1, My - 1); m <= mm; ++m) {
            for (int l = max(nx[0] - 1, 0), lm = min(nx[0] + 1, Mx - 1); l <= lm; ++l) {
                map_index = l + Mx * m;
                map[map_index].emplace_back(i);
            }
        }
    }
    // int km, j;

    for (int i = 0; i < Np; ++i) {
        list[i][0] = 0;
        map_index  = (int) ((x[i][0] + xlen_2) * bitx) + Mx * (int) ((x[i][1] + R) * bity);
        for (int j = 0, jm = map[map_index].size(); j < jm; ++j) {
            if (map[map_index][j] > i) {
                dx = x[i][0] - x[map[map_index][j]][0];
                dy = x[i][1] - x[map[map_index][j]][1];
                if ((dx * dx + dy * dy) < threash2) {
                    list[i][0]++;
                    list[i][list[i][0]] = map[map_index][j];
                }
            }
        }
    }
    // delete[] map;
}
void ver_list(int (*list)[Nn], double (*x)[dim]) {
    double           dx, dy, dr2;
    constexpr double thresh2 = (cut + gay_kappa - 1 + skin) * (cut + gay_kappa - 1 + skin);
    for (int i = 0; i < Np; i++) {
        list[i][0] = 0;
    }

    for (int i = 0; i < Np; i++) {
        for (int j = i + 1; j < Np; j++) {
            dx  = (x[i][0] - x[j][0]);
            dy  = (x[i][1] - x[j][1]);
            dr2 = dx * dx + dy * dy;
            if (dr2 < thresh2) {
                list[i][0]++;
                list[i][list[i][0]] = j;
            }
        }
    }
}

#endif