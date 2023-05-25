
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <BM.h>

#define Np 4
#define dim 2
#define tau 1.
#define om 0.
#define v0 1.
#define dt 0.1


void eom_abp(double (*x)[dim],double (*v)[dim],double *theta_i,double (*f)[dim]){
    double D=sqrt(2./(tau*dt)),M_PI2=2.*M_PI;
    for(int i=0;i<Np;i++){
        theta_i[i]+=(D*gaussian_rand()+om)*dt;
        theta_i[i]-=(int)(theta_i[i]*M_1_PI)*M_PI2;
        v[i][0]=v0*cos(theta_i[i])+f[i][0];
        v[i][1]=v0*sin(theta_i[i])+f[i][1];
        x[i][0]+=v[i][0]*dt;
        x[i][1]+=v[i][1]*dt;
    }

}

void eom_abp2(double (*x)[dim],double (*v)[dim],double *theta_i,double (*f)[dim]){
    double D=sqrt(2./(tau*dt)),M_PI2=2.*M_PI,Co;
    for(int i=0;i<Np;i++){
        theta_i[i]+=(D*gaussian_rand()+om)*dt;
        theta_i[i]-=(int)(theta_i[i]*M_1_PI)*M_PI2;
        Co=cos(theta_i[i]);
        v[i][0]=v0*Co+f[i][0];
        v[i][1]=v0*(theta_i[i]>0-theta_i[i]<0)*sqrt(1.-Co*Co)+f[i][1];
        x[i][0]+=v[i][0]*dt;
        x[i][1]+=v[i][1]*dt;
    }

}