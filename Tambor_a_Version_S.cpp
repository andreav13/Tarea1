#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;


// r²*R''+r*R'+lambda²r²*R=0  ; x1=R  x2=R'  f1=x2 f2= -(1/r)*x2 - lambda²x1
// C.I => x1(0)=1 ; x2(0)= 0

const double lambda=1;


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



int main(void){
  double x1=1,x2=0,r;
  double dr=0.01;
  for(r=0.01;r<10;)
    {
      cout<<r<<"  "<<x1<<endl ;
      UnPasoDeRungeKuta4(lambda,x1,x2,r,dr);

    }
  
   return 0;
}
