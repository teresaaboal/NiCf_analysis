#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPaveStats.h"

void analysis(TString infile,
	      const int     start_lo = +99999,
	      const int     start_hi = -99999,
	      const char *  hname = "likelihood")
{
  TFile * f = TFile::Open(infile.Data());
  TH1D  * h;
  f->GetObject(hname, h);

  const double xmin_content = h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove());
  const double xmin = h->GetXaxis()->GetXmin();
  const double xmax = h->GetXaxis()->GetXmax();

  std::cout << h->GetNbinsX() << " bins from " << xmin << " to " << xmax << std::endl
	    << xmin_content << " is lower edge of content" << std::endl;

  //expo(0) = exp([0]+[1]*x)
  //gaus(0) = [0]*exp(-0.5*((x-[1])/[2])**2)
  const double first_range_min = (start_lo < +99998) ? start_lo : xmin_content;
  const double first_range_max = (start_hi > -99998) ? start_hi : xmax;
  TF1 *fit = new TF1("fit", "gaus(0)*(x<1.5)+(x>=1.5)*(expo(3)+expo(5)+expo(7))", first_range_min, first_range_max);
  fit->SetLineColor(kRed);

  TCanvas * can = new TCanvas();
  can->SetLogy(1);
  //can->SetTopMargin(0.2);

  //ensure the bin errors are sensible - make them sqrt content
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) {
    const double x = h->GetBinContent(ix);
    if(x > 1E-6)
      h->SetBinError(ix, TMath::Sqrt(x));
  }
  h->SetMarkerStyle(7);
  h->Draw();

  //First time in the loop, set initial values and attempt to fit a smaller range (if given)
  //Second time, free the range to the entire histogram
  for(int i = 0; i < 2; i++) {
 
    //set initial parameters
    if(i == 0) {
      fit->SetParameter(0,5.85655e+04);
      fit->SetParameter(1,4.72705e-01);
      fit->SetParameter(2,2.68512e+00);
      fit->SetParameter(3,7.53496e+00);
      fit->SetParameter(4,-9.49172e-03);
      fit->SetParameter(5,1.05580e+01);
      fit->SetParameter(6,-2.81987e-01);
    }
    if(i == 1) {
      fit->SetRange(xmin_content, xmax);
    }
    //and fit
    h->Fit("fit", "RME");

    //reset and fit again
    fit->SetParameter(7,0);
    fit->SetParameter(8,0);
    h->Fit("fit", "RME");

    //reset and fit again
    fit->SetParameter(7,0);
    fit->SetParameter(8,0);
    h->Fit("fit", "RME");
  }

  h->ls();
  h->Draw();
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h->FindObject("stats");
  st->SetX1NDC(0.3); //new x start position
  st->SetX2NDC(0.7); //new x end position
  st->SetY1NDC(0.2); //new y start position
  st->SetY2NDC(0.6); //new y end position
  st->SetTextFont(42); //less bold
  st->Paint();
  gPad->Update();

  //save the canvas
  can->SaveAs(infile.ReplaceAll(".root", "_fit.pdf"));

  //std::cout <<" dx=((i*TBIN+" << fit->GetParameter(1) << ")/" << fit->GetParameter(2) << " );" << std::endl;
  std::cout << std::endl << std::endl
       << "The following should be copied into create_like.cc" << std::endl
       << std::endl
       <<" dx=(i*TBIN/" << fit->GetParameter(2) << " );" << std::endl
       <<" if (x>-1.5)"<< std::endl
       <<" dhisto[i+nneg]+=1*exp(-0.5*dx*dx);"<< std::endl
       <<" if (x<=-1.5)"<< std::endl
       <<" {"<< std::endl
       <<"   dhisto[i+nneg]+=" << exp(fit->GetParameter(3))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(4)) << "*x);"<< std::endl
       <<"   dhisto[i+nneg]+=" << exp(fit->GetParameter(5))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(6)) << "*x);"<< std::endl
       <<"   dhisto[i+nneg]+=" << exp(fit->GetParameter(7))/fit->GetParameter(0) << "*exp(" << -1*(fit->GetParameter(8)) << "*x);"<< std::endl
       <<" }"<< std::endl;
}
