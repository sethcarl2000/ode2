#include "RK.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cmath> 
#include <TAxis.h> 
#include <TLine.h>

using namespace std;

// differential equation to be solved
double fun1(double x, double y){
  (void)x;              // prevent unused variable warning
  return -2*y;          // f = y'(x,y) = -2 * y(x)  
}                       // solution: y(x) = 3 * exp(-2*x) ; with initial condition y(0)=3

double fun2(double x, double y){
  return -y/x-2/(x*x);  // f = y'(x,y) = -y(x)/x - 2/x^2 
}                       // -2*log(|x|)/x+2/x  ; with initial condition y(0)=2

//bernouli equation
double f_bernouli(double x, double y){
  return y*(2./x) - (x*x*y*y); 
}
//

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  const double x_min = 0.1; 
  const double x_max = 4.0;

  const double y0    = 0.1;
  
  const unsigned int npts = 15; 

  // solve our DEQ using RK1 or RK2 methods!
  // Two examples are given.  Choose a function for testing
  TGraph tg1=RK1Solve(f_bernouli, x_min, y0, npts, x_max);                     // initial condition y(0)=3
  TGraph tg2=RK2Solve(f_bernouli, x_min, y0, npts, x_max);
  TGraph tg4=RK4Solve(f_bernouli, x_min, y0, npts, x_max);
  TF1 fun_sol=TF1("fun_sol","pow(x,2)/(pow(x,5)/5. + 0.1)",x_min,x_max);           // exact solution
  //TGraph tg1=RK1Solve(fun2,2,100,1,100);                // initial condition y(1)=2
  //TGraph tg2=RK2Solve(fun2,2,100,1,100);
  //TF1 fun_sol=TF1("fun_sol","-2*log(x)/x+2/x",1,100);   // exact solution

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  TCanvas *c1 = new TCanvas("c1","DEQ solutions",dw,dh);
  c1->Divide(1,2); 

  tg1.SetMarkerSize(0.015*dh/8);  // size scale: 1 = 8 pixels, so here we choose the size to be 1.5% of the window height
  tg2.SetMarkerSize(0.015*dh/8);
  tg4.SetMarkerSize(0.015*dh/8);

  tg4.GetXaxis()->SetRangeUser(0., 4.);
  

  tg1.SetMarkerStyle(kFullTriangleUp);
  tg2.SetMarkerStyle(kFullTriangleDown);
  tg4.SetMarkerStyle(kOpenSquare);
  
  tg1.SetMarkerColor(kRed);     tg1.SetLineColor(kRed);
  tg2.SetMarkerColor(kGreen-2); tg2.SetLineColor(kGreen-2);
  tg4.SetMarkerColor(kBlue);    tg4.SetLineColor(kBlue);
  
  fun_sol.SetLineColor(kBlack);
  fun_sol.SetLineStyle(2);
  
  // plot the results
  c1->cd(1);
  tg4.SetTitle("ODE demo: #partial y = 2y/x - (xy)^{2};x;y");
  tg4.Draw("AP SAME"); 
  tg1.Draw("P SAME");
  tg2.Draw("P SAME");
  fun_sol.Draw("same");
  
  TLegend *tl = new TLegend(0.6,0.7,0.9,0.9);
  tl->AddEntry(&tg1,"RK1 Solution","p");
  tl->AddEntry(&tg2,"RK2 Solution","p");
  tl->AddEntry(&tg4,"RK4 Solution","p"); 
  tl->AddEntry(&fun_sol,"Exact Solution","l");
  tl->Draw();
  c1->Draw();
  c1->Update();
  c1->Print("OED_cpp.png");

  // retreive the data from the graphs and write to a file
  FILE *fp=fopen("RKdemo.dat","w");
  double *x, *y1, *y2;
  x=tg1.GetX();
  y1=tg1.GetY();
  y2=tg2.GetY();
  fprintf(fp,"#%8s %9s %9s %9s\n","x","RK1","RK2","Exact");
  for (int i=0; i<tg1.GetN(); i++){
    fprintf(fp,"%9.4lf %9.4lf %9.4lf %9.4lf\n",x[i],y1[i],y2[i],fun_sol.Eval(x[i]));
  }
  fclose(fp);

  //______________________________________________________________________________________________
  auto Compute_residuals = [npts, &fun_sol](const TGraph &g)
  {
    //make copies of our data
    const double* Xptr = g.GetX(); 
    const double* Yptr = g.GetY(); 
    vector<double> X( Xptr, Xptr+npts-1 ); 
    vector<double> Y( Yptr, Yptr+npts-1 ); 

    TGraph g_out(g);

    for (int i=0; i<npts; i++) {
      g_out.SetPointY(i, Y[i] - fun_sol.Eval(X[i])); 
    }
    
    return g_out; 
  };
  //______________________________________________________________________________________________
  
  auto tg1_resid = Compute_residuals(tg1);
  auto tg2_resid = Compute_residuals(tg2); 
  auto tg4_resid = Compute_residuals(tg4); 
  
  c1->cd(2);
  tg1_resid.SetTitle(";x;Residuals - y");
  tg1_resid.Draw("ALP");
  tg4_resid.Draw("LP"); 
  tg2_resid.Draw("LP");


  cout << "Press ^c to exit" << endl;
  //theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

