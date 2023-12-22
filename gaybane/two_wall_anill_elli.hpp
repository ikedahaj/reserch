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
void ver_list_anill(int (*list)[Nn], double (*x)[dim], double R_ini) {
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
#endif
