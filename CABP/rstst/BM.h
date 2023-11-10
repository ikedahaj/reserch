double unif_rand(double left, double right) {
    static constexpr double RAND_MAX_1=1./RAND_MAX;
    return left + (right - left) * rand() *RAND_MAX_1;
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
