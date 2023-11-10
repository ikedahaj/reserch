#include <bits/stdc++.h>

using namespace std;

int set(int a,int b);
int set_some(int *x){
    for (int i = 0; i < 10; i++)
    {
        x[i]=i;
    }
    cout<<"pass!"<<endl;
    return 0;
}
int calc_sum(int x){
    int sum=0;
    int x0[10];
    static int vv=set_some(x0);
    for (int i = 9; i >=0; i--)
    {
        sum+=(x+i-x0[9-i])*(x+i-x0[9-i]);
    }
    cout<<sum<<vv<<endl;
    return sum;
}
int main(){
    int x[10];
    cout <<"test"<< '\t' <<"test"<< endl;
    return 0;
}
int set(int a,int b){

}