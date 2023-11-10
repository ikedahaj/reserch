#include <chrono>
#include <iostream>

#define pll 0.00001
using std::endl;
class params {
  public:
    static constexpr double some = 0.00001;
    inline constexpr double ret_son() { return 0.00001; }
};
double                  param = 0.00001;
static constexpr double pp = 0.00001;

int main() {
    double                                sum;
    std::chrono::system_clock::time_point start, end; // 型は auto で可;
    start = std::chrono::system_clock::now();         // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += param;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "global double time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換;
    // std::chrono::system_clock::time_point start, end; // 型は auto で可;
    start = std::chrono::system_clock::now(); // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += params::some;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "class static constexper time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換
    // std::chrono::system_clock::time_point start, end; // 型は auto で可;
    start = std::chrono::system_clock::now(); // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += pp;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "global static constexper time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換
    static constexpr double pl = 0.00001;
    start = std::chrono::system_clock::now(); // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += pl;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "main() static constexper time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換
    start = std::chrono::system_clock::now(); // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += pll;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "#define time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換
    start = std::chrono::system_clock::now(); // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += 0.00001;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "direct time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換
    start = std::chrono::system_clock::now();         // 計測開始時間;
    sum = 1;
    for (int i = 0; i < 1e9; i++) {
        // for (int j = 0; j <  1e9; j++)
            sum += param;
    }
    std::cout << sum << std::endl;
    end = std::chrono::system_clock::now(); // 計測終了時間;
    std::cout << "global double time"
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
                                                                       start)
                     .count()
              << endl; // 処理に要した時間をミリ秒に変換;
}