#include <iostream>

using namespace std;

static constexpr int some=5;

#include "test2.cpp"
int pop[5];
void output(){
    for (int i = 0; i < 5; i++)
    {
        cout<<pop[i]<<endl;
    }
    
}
void time(){
    pop[4]*=2;
}
int *ini(){
    int naka[5]={0,1,2,3,5};
    return naka;
}
int main(){

    
    cout << add_some(9) << endl;
    return 0;
}

