#include "RK.hpp"
#include <vector> 


TGraph RK1Solve(
  double (*dy_dx)(double x, double y),
	double x0, 
	double y0,
	const unsigned int nsteps, 
	double xmax
){
  
  std::vector<double> X{x0}, Y{y0}; 
  X.reserve(nsteps); 
  Y.reserve(nsteps);

  const double dx = (xmax - x0)/((double)nsteps); 

  double x = x0; 
  double y = y0; 

  for (unsigned int i=1; i<nsteps; i++) {
  
    y += dy_dx(x,y)*dx; 
    x += dx; 
    
    X.push_back(x);
    Y.push_back(y);
  }

  return TGraph(nsteps, X.data(), Y.data()); 
}

TGraph RK2Solve(
  double (*dy_dx)(double x, double y),
	double x0, 
	double y0,
	const unsigned int nsteps, 
	double xmax
){

  std::vector<double> X{x0}, Y{y0}; 
  X.reserve(nsteps); 
  Y.reserve(nsteps);

  const double dx = (xmax - x0)/((double)nsteps); 

  double x = x0; 
  double y = y0; 

  for (unsigned int i=1; i<nsteps; i++) {

    double k1 = dy_dx(x,      y); 
    double k2 = dy_dx(x + dx, y + dx*k1); 

    y += dx*(k1 + k2)/2.; 
    x += dx; 

    X.push_back(x);
    Y.push_back(y);
  }

  return TGraph(nsteps, X.data(), Y.data()); 
}

TGraph RK4Solve(
  double (*dy_dx)(double x, double y),
	double x0, 
	double y0,
	const unsigned int nsteps, 
	double xmax
){

  std::vector<double> X{x0}, Y{y0}; 
  X.reserve(nsteps); 
  Y.reserve(nsteps);

  const double dx = (xmax - x0)/((double)nsteps); 

  double x = x0; 
  double y = y0; 

  for (unsigned int i=1; i<nsteps; i++) {

    double k1 = dy_dx(x,         y); 
    double k2 = dy_dx(x + dx/2., y + k1*dx/2.); 
    double k3 = dy_dx(x + dx/2., y + k2*dx/2.);
    double k4 = dy_dx(x + dx,    y + k3*dx);

    y += dx*(k1 + k2*2. + k3*2. + k4)/6.; 
    x += dx; 

    X.push_back(x);
    Y.push_back(y);
  }

  return TGraph(nsteps, X.data(), Y.data()); 
}
