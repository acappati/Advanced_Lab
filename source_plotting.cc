#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>

#include <TSystem.h>
#include <TMatrixF.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TStopwatch.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>

using namespace std;

Double_t Func_Track_Gaus(Double_t *, Double_t *);
Double_t Func_Sum_Track_Gaus(Double_t *, Double_t *);

const int dimX1 = 9;
const int dimX2 = 6;
const int dimX1X2 = 15;
TF2 *funcTrackGaus[dimX1X2];

void source_plotting(){
  gStyle -> SetOptStat(0);
  TH2D *histGrid = new TH2D("histGrid","",100,-7.,7.,100,-7.,7.);

  double theta = 0;

  double trackNorm1[dimX1] = {18622.468179,20646.549073,20703.446134,21402.483777,19284.718767,17878.501612,19759.205335,19426.203887,21887.383};
  double trackY1[dimX1] = {9.790181,9.238331,7.666330,5.735672,4.496143,4.369718,5.439571,7.229848,8.852554};
  double trackErrY1[dimX1] = {0.273795,0.318383,0.329923,0.300995,0.339809,0.367801,0.343508,0.338200,0.272218};
  double trackTheta1[dimX1] = {0,40,80,120,160,200,240,280,320};

  TF1 *funcTrack1[dimX1];
  for(int i = 0;i < dimX1;i++){
    trackY1[i] = trackY1[i] - 7.;
    theta = (trackTheta1[i]/180.)*TMath::Pi();

    TMatrixF matrixRotation(2,2);
    matrixRotation(0,0) = TMath::Cos(theta);
    matrixRotation(0,1) = -TMath::Sin(theta);
    matrixRotation(1,0) = TMath::Sin(theta);
    matrixRotation(1,1) = TMath::Cos(theta);

    TMatrixF matrixXY(2,1);
    matrixXY(0,0) = 12.;
    matrixXY(1,0) = trackY1[i];

    TMatrixF matrixXYRotation(2,1);
    matrixXYRotation = matrixRotation*matrixXY;
    printf("Before rotation (%f,%f) - After rotation (%f,%f) \n",matrixXY(0,0),matrixXY(1,0),matrixXYRotation(0,0),matrixXYRotation(1,0));

    funcTrack1[i] = new TF1("funcTrack","[0] + [1]*(x - [2])",-12.,12.);
    funcTrack1[i] -> SetParameter(0,matrixXYRotation(1,0));
    funcTrack1[i] -> SetParameter(1,TMath::Tan(theta));
    funcTrack1[i] -> SetParameter(2,matrixXYRotation(0,0));
    funcTrack1[i] -> SetLineColor(kRed);

    funcTrackGaus[i] = new TF2("funcTrackGaus",Func_Track_Gaus,-5,5,-5,5,5);
    funcTrackGaus[i] -> SetNpx(500);
    funcTrackGaus[i] -> SetNpy(500);
    funcTrackGaus[i] -> SetParameter(0,trackNorm1[i]);
    funcTrackGaus[i] -> SetParameter(1,matrixXYRotation(1,0));
    funcTrackGaus[i] -> SetParameter(2,TMath::Tan(theta));
    funcTrackGaus[i] -> SetParameter(3,matrixXYRotation(0,0));
    funcTrackGaus[i] -> SetParameter(4,trackErrY1[i]);
  }

  printf("======================================================================== \n");
  double trackNorm2[dimX2] = {196.135459,197.381551,149.640702,182.017845,193.836493,175.018336};
  double trackY2[dimX2] = {5.317608,4.664446,5.196628,8.166169,9.238523,9.407446};
  double trackErrY2[dimX2] = {0.345687,0.286275,0.540735,0.373195,0.355325,0.396841};
  double trackTheta2[dimX2] = {0,40,80,160,200,240};

  TF1 *funcTrack2[dimX2];
  for(int i = 0;i < dimX2;i++){
    trackY2[i] = trackY2[i] - 7.;
    theta = (trackTheta2[i]/180.)*TMath::Pi();

    TMatrixF matrixRotation(2,2);
    matrixRotation(0,0) = TMath::Cos(theta);
    matrixRotation(0,1) = -TMath::Sin(theta);
    matrixRotation(1,0) = TMath::Sin(theta);
    matrixRotation(1,1) = TMath::Cos(theta);

    TMatrixF matrixXY(2,1);
    matrixXY(0,0) = 12.;
    matrixXY(1,0) = trackY2[i];

    TMatrixF matrixXYRotation(2,1);
    matrixXYRotation = matrixRotation*matrixXY;
    printf("Before rotation (%f,%f) - After rotation (%f,%f) \n",matrixXY(0,0),matrixXY(1,0),matrixXYRotation(0,0),matrixXYRotation(1,0));

    funcTrack2[i] = new TF1("funcTrack","[0] + [1]*(x - [2])",-12.,12.);
    funcTrack2[i] -> SetParameter(0,matrixXYRotation(1,0));
    funcTrack2[i] -> SetParameter(1,TMath::Tan(theta));
    funcTrack2[i] -> SetParameter(2,matrixXYRotation(0,0));
    funcTrack2[i] -> SetLineColor(kBlue);

    funcTrackGaus[i+9] = new TF2("funcTrackGaus",Func_Track_Gaus,-5,5,-5,5,5);
    funcTrackGaus[i+9] -> SetNpx(500);
    funcTrackGaus[i+9] -> SetNpy(500);
    funcTrackGaus[i+9] -> SetParameter(0,trackNorm2[i]);
    funcTrackGaus[i+9] -> SetParameter(1,matrixXYRotation(1,0));
    funcTrackGaus[i+9] -> SetParameter(2,TMath::Tan(theta));
    funcTrackGaus[i+9] -> SetParameter(3,matrixXYRotation(0,0));
    funcTrackGaus[i+9] -> SetParameter(4,trackErrY2[i]);
  }

  TCanvas *canvasTrack = new TCanvas("canvasTrack","canvasTrack",20,20,600,600);
  histGrid -> Draw();
  for(int i = 0;i < dimX1;i++){funcTrack1[i] -> Draw("same");}
  for(int i = 0;i < dimX2;i++){funcTrack2[i] -> Draw("same");}

  //============================================================================
  //Int_t colors[] = {0, 1, 2, 3, 4, 5, 6}; // #colors >= #levels - 1
  //gStyle -> SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  //Double_t levels[] = {-3.4e38, 1.17e-38, 0.90, 0.95, 1.00, 1.05, 1.10, 3.4e38};

  TF2 *FuncSumTrackGaus = new TF2("FuncSumTrackGaus",Func_Sum_Track_Gaus,-7,7,-7,7,0);
  //FuncSumTrackGaus -> SetContour((sizeof(levels)/sizeof(Double_t)), levels);
  FuncSumTrackGaus -> SetNpx(500);
  FuncSumTrackGaus -> SetNpy(500);

  TCanvas *canvasTrackGaus = new TCanvas("canvasTrackGaus","canvasTrackGaus",20,20,600,600);
  histGrid -> Draw();
  FuncSumTrackGaus -> Draw("sameCOLZ");


}
////////////////////////////////////////////////////////////////////////////////
Double_t Func_Track_Gaus(Double_t *x, Double_t *par){
  double PI = TMath::Pi();
  double nX = par[0];
  double meanX = par[1] + par[2]*(x[0] - par[3]);
  double sigmaX = par[4];

  return (1/(TMath::Sqrt(2*PI)*sigmaX))*TMath::Exp(-(TMath::Power(x[1] - meanX,2))/(2*sigmaX*sigmaX));
}
//==============================================================================
Double_t Func_Sum_Track_Gaus(Double_t *x, Double_t *par){
  double funcSumTrackGaus = 0;
  for(int i = 0;i < dimX1X2;i++){funcSumTrackGaus += funcTrackGaus[i] -> EvalPar(x,par);}
  return funcSumTrackGaus;
}
