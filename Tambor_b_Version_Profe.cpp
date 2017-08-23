#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

double f1(double x1,double x2,double r){return x2;}

double f2(double lambda,double x1,double x2,double r){
  return -( (1/r)*x2 + pow(lambda,2)*x1 );
}

void UnPasoDeRungeKuta4(double lambda,double & x1,double & x2, double & r,double dr){
  
  double dx11,dx21,dx31,dx41;                double dx12,dx22,dx32,dx42;
  dx11=dr*f1(x1,x2,r);                       dx12=dr*f2(lambda,x1,x2,r);
  dx21=dr*f1(x1+dx11/2,x2+dx12/2,r+dr/2);    dx22=dr*f2(lambda,x1+dx11/2,x2+dx12/2,r+dr/2);
  dx31=dr*f1(x1+dx21/2,x2+dx22/2,r+dr/2);    dx32=dr*f2(lambda,x1+dx21/2,x2+dx22/2,r+dr/2);
  dx41=dr*f1(x1+dx31,x2+dx32,r+dr);          dx42=dr*f2(lambda,x1+dx31,x2+dx32,r+dr);
       
  r+=dr; x1+=(dx11+2*dx21+2*dx31+dx41)/6;    x2+=(dx12+2*dx22+2*dx32+dx42)/6;
}

void f(double  lambda,double & R){  
  double x1=1,x2=0;
  double r;  const double dr=0.01;
  for (r=0.01;r<=5;){
    UnPasoDeRungeKuta4(lambda,x1,x2,r,dr);
    R=x1;
  }
}


int main(void){
  double lambda; const double dlambda=0.01;
  double R=0.0;
  
  for (lambda=0.0;lambda<80.0;lambda+=dlambda ){
    f(lambda,R);
    cout<<lambda<<"  "<<R<<endl;
  }
  return 0;
}
