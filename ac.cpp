#include <bits/stdc++.h>

using namespace std;

int main() {
    int         n, m;
    vector<int> a(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
    }

    sort(a.begin(), a.end());
    unsigned long long int ans = 0;
    for (int i = 0; i < n - m; i++) {
        ans += (a[n - m + i] + a[n - m - 1 - i]) *
               (a[n - m + i] + a[n - m - 1 - i]);
    }
    for (int i = 2 * (n - m); i < n; i++) {
        ans += a[i] * a[i];
    }
    cout << ans << endl;

    return 0;
}