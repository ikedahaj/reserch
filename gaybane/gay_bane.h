#ifndef GAY_BANE
#define GAY_BANE

#define gay_kapppa   2
#define gay_kapppa_p 5
#define gay_M        2
#define gay_N        1

// baseのpower乗を和えす関数;
inline double fast_power(double base, int power) {
    double result = 1;
    while (power > 0) {
        if (power & 1) {  // Can also use (power & 1) to make code even faster
            result = (result * base);
        }
        base = (base * base);
        power >>= 1;  // Can also use power >>= 1; to make code even faster
    }
    return result;
}
constexpr double fast_power_cone(double base, int power) {
    double result = 1;
    while (power > 0) {
        if (power % 2 ==
            1) {  // Can also use (power & 1) to make code even faster
            result = (result * base);
        }
        base = (base * base);
        power =
            power / 2;  // Can also use power >>= 1; to make code even faster
    }
    return result;
}
constexpr double nth_root(double x, int t) {
    double x0 = x, ans = (x > 1) ? x : 1;
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

void calc_force_gay_bane(double (*x)[dim], double (*f)[dim], int (*list)[Nn],
                         double *theta, double *f_theta) {
    constexpr double gay_alpha = (gay_kapppa * gay_kapppa - 1) /
                                 (gay_kapppa * gay_kapppa + 1),
                     gay_alpha_p = (nth_root(gay_kapppa_p, 1. / gay_M) - 1) /
                                   (nth_root(gay_kapppa_p, 1. / gay_M)),
                     cutf2 = cut * cut * gay_kapppa * gay_kapppa;
    // kがi,lがj;
    double dx, dy, dr2, cosk, sink, cosl, sinl, cos_kl_plus, cos_kl_minus,
        sin_kl_plus, sin_kl_minus, f_alpha, f_alpha_p, dr_dot_ek, dr_dot_el,
        sigma_hat, dr, xi, xi2, xi6, UWCA, dUWCA, dr_1, dr2_1, dx_xi, dy_xi,
        d_theta_xi_k, d_theta_xi_l, epsilon1, epsilon1_N_1 /*epsilon1のN-1乗*/,
        epsilon2, epsilon2_M_1 /*epsilon2のM-1乗*/;
    ini_array(f);
    for (int i = 0; i < Np; ++i)
        for (int j = 1; j <= list[i][0]; ++j) {
            dx = x[i][0] - x[list[i][j]][0];
            dy = x[i][1] - x[list[i][j]][1];
            dr2 = dx * dx + dy * dy;
            // aij = (a[i] + a[list[i][j]]) * (a[i] + a[list[i][j]]);
            if (dr2 < cutf2) {  // 長径に合わせたカットおふ;
                cosk = cos(theta[i]);
                sink = sin(theta[i]);
                cosl = cos(theta[list[i][j]]);
                sinl = sin(theta[list[i][j]]);
                cos_kl_minus = cosk * cosl + sink * sinl;
                dr_dot_ek = dx * cosk + dy * sink;
                dr_dot_el = dx * cosl + dy * sinl;
                f_alpha = (dr_dot_ek + dr_dot_el) * (dr_dot_ek + dr_dot_el) /
                              (1 + gay_alpha * cos_kl_minus) +
                          (dr_dot_ek - dr_dot_el) * (dr_dot_ek - dr_dot_el) /
                              (1 - gay_alpha * cos_kl_minus);
                sigma_hat = 1. / sqrt(1 - gay_alpha * 0.5 * (f_alpha));
                dr = sqrt(dr2);
                xi = dr - sigma_hat + 1;
                if (xi < cut) {
                    xi2 = 1. / (xi * xi);
                    xi6 = xi2 * xi2 * xi2;
                    UWCA = xi6 * (xi6 - 1.);
                    dUWCA = 12 * xi6 * (0.5 - xi6) / xi;
                    dr_1 = 1. / dr;
                    dr2_1 = dr_1 * dr_1;
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
                    sin_kl_minus = sink * cosl - cosk * sinl;
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
                              (-dx * sink + dy * cosk) *
                              (1 - gay_alpha * cos_kl_minus) -
                          gay_alpha * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha * cos_kl_minus) *
                              (1 - gay_alpha * cos_kl_minus)));
                    d_theta_xi_l =
                        0.5 * gay_alpha * dr2_1 *
                        ((2 * (dr_dot_ek + dr_dot_el) *
                              (-dx * sinl + dy * cosl) *
                              (1 + gay_alpha * cos_kl_minus) +
                          gay_alpha * sin_kl_minus * (dr_dot_ek + dr_dot_el) *
                              (dr_dot_ek + dr_dot_el)) /
                             ((1 + gay_alpha * cos_kl_minus) *
                              (1 + gay_alpha * cos_kl_minus)) +
                         (2 * (dr_dot_ek - dr_dot_el) *
                              (-dx * sinl + dy * cosl) *
                              (1 - gay_alpha * cos_kl_minus) -
                          gay_alpha * sin_kl_minus * (dr_dot_ek - dr_dot_el) *
                              (dr_dot_ek - dr_dot_el)) /
                             ((1 - gay_alpha * cos_kl_minus) *
                              (1 - gay_alpha * cos_kl_minus)));
                    // 力の計算をするときにマイナスをつけること;
                    
                }
            }
        }
}
#endif