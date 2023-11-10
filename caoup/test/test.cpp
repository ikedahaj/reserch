#include <iostream>

#define R 80.
class testofclass{
    public:
        constexpr static int ohimesama=R*50.;
};
class TestClass
{
    int num;

    //静的メンバ変数
    static int sNum;
    
public:
    // const static double k;
    TestClass(int n = 0) { num = n; }
    int get() { return num; }

    int getStatic() { return sNum; }
    void setStatic(int n) { sNum = n; }
};
const static double k=5;
//静的メンバ変数の定義
int TestClass::sNum = 0;
double fumk_ret(){
    return compiler;
}

int main()
{
    TestClass tc1(10);
    TestClass tc2(20);
    constexpr double tett=compiler+55;
    tc1.setStatic(50);

    std::cout << tc1.get() << std::endl; //10
    std::cout << tc2.get() << std::endl; //20
    std::cout<<fumk_ret();
    // std::cout << tc2.getStatic() << std::endl; //50
    std::cout << tc2.getStatic() << std::endl; //50
    std::cout<<1./testofclass::ohimesama<<std::endl;
    // std::cin.get();
    return 0;
}