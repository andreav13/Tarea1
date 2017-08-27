#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

const double dr=0.01;

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

double f(double  lambda){  
  double x1=1,x2=0;
  double r;
  for (r=0.01;r<=5;){
    UnPasoDeRungeKuta4(lambda,x1,x2,r,dr);
  }
  return x1;

}





void ZerosPorBis(double a, double b,double & zero){
  const double Err=1e-5;
  double m;
  
  if (f(a)*f(b)<0){
    do{
      m=(a+b)/2;	
      if(f(a)*f(m)<0){b=m;}
      else{a=m;}
    }while(abs(f(m))>Err);
    zero=m;
  }
  else if(f(a)*f(b)==0){
    if(f(a)==0){zero=f(a);}
    else{zero=f(b);}
	
  }
  else{cout<<"No hay un cero en el intervalo"<<endl;}
}

void  CerosFinder(int n){
  double x,k=0.001,y;
  double zero=0.0;
  int ceros=0;
  double Matriz[n];
  for(x=0;ceros<n;x+=k){
    y=x+k;
    if (f(x)*f(y)<=0){
      ZerosPorBis(x,y,zero);
      cout<<"R(r=5)=0 en lambda = "<<zero<<endl;
      Matriz[ceros]=zero;
      ceros+=1;
      }
    else{ceros=ceros;}   
  }


  //-----------Código para la obtención de los modos de vibración para R(a=5) en diferentes archivos.dat------//
  
    for(int i=0;i<n;i++){
      double x1=1,x2=0;
      
      ofstream R_r;
      R_r.open("R_r_"+std::to_string(Matriz[i])+".dat");
    
      for(double r=0.01;r<=5.0;){
	UnPasoDeRungeKuta4(Matriz[i],x1,x2,r,dr);
	R_r<<r<<" "<<x1<< "\n";
      }
      R_r.close();

    }
  
  }




int main(void){
  int n=5;
  CerosFinder(n);
  return 0;
}
