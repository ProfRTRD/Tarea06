#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>


struct functor_f{ //friccion
  double  V = 45.0, RE , mu =1.79e-5, D, ep = 14e-6;
  double operator()(double x)
  {
    double a = (ep/3.7*D);
    double b = (2.51/(RE* std::sqrt(x)));
    double c = (1.0/std::sqrt(x));
    return (2.0*std::log10(a+b))+c;
  }
};

struct functor_g{//reynolds
  double ro = 1.23, V =45.0, D, mu= 1.79e-5;
  double operator()()
  {
    return (ro * V * D)/mu;
  }
};

struct functor_h{//delta p
  double Fr, L = 0.2, ro = 1.23, V = 45.0, D;
  double operator ()(){
    return (Fr*L*ro*std::pow(V,2))/(2.0*D);
  }
};

  
template <typename func_t>
double newton(double x0, double eps, func_t func, int nmax, int & nsteps);

int main (int argc, char *argv[])
{
  std::ofstream fout{"datos.txt"};
  double eps = std::atof(argv[1]);//1.0e-5
  int NMAX = 1000;
  fout.precision(15); 
  fout.setf(std::ios::scientific);
  int steps = 0;

  // root for function f
  functor_f fr;
  int n = 30;
  double D_min = 0.001, D_max = 0.030;
  double D_delta = (D_max - D_min)/n;
  for (int ii = 0; ii < (n-1); ++ii) {
       double Dt  = D_min + ii*D_delta;
    //reynolds
    functor_g reynolds;
    reynolds.D = Dt;
    double Re = reynolds();
    
    //friccion
    functor_f friccion;
    friccion.RE = Re;
    double X0 = 0.316/std::pow(Re,0.25);
    double root = newton(X0, eps, friccion, NMAX, steps);
 
    //deltap
    functor_h deltap;
    deltap.D = Dt;
    deltap.Fr = root;
    double result = deltap();
    
    fout << Dt <<"\t" << result <<"\n";
  }
  fout.close();
}




template <typename func_t>
double newton(double x0, double eps, func_t friccion, int nmax, int & nsteps)
{
  nsteps = 0;
  double xr = x0;
  while(nsteps <= nmax) {
    if (std::fabs(friccion(xr)) < eps) {
      break;
    } else {
      double h = 0.001;
      double deriv = (friccion(xr+h/2) - friccion(xr-h/2))/h;
      xr = xr - friccion(xr)/deriv;
    }
    nsteps++;
  }

  return xr;
}


