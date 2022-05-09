template <typename func_t>
double newton(double x0, double eps, func_t func, int nmax, int & nsteps)
{
  nsteps = 0;
  double xr = x0;
  while(nsteps <= nmax) {
    if (std::fabs(func(xr)) < eps) {
      break;
    } else {
      double h = 0.001;
      double deriv = (func(xr+h/2) - func(xr-h/2))/h;
      xr = xr - func(xr)/deriv;
    }
    nsteps++;
  }

  return xr;
}
