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

using namespace std;

void source_plotting(){
  gStyle -> SetOptStat(0);
  TH2D *histGrid = new TH2D("histGrid","",100,-7.,7.,100,-7.,7.);

  double theta = 0;

  const int dimX1 = 9;
  double trackY1[dimX1] = {9.790181,9.238331,7.666330,5.735672,4.496143,4.369718,5.439571,7.229848,8.852554};
  double trackErrY1[dimX1] = {0.273795,0.318383,0.329923,0.300995,0.339809,0.367801,0.343508,0.338200,0.272218};
  double trackTheta1[dimX1] = {0,40,80,120,160,200,240,280,320};

  TF1 *funcTrack1[dimX1];
  for(int i = 0;i < dimX1;i++){
    trackY1[i] = trackY1[i] - 7.;
    theta = (trackTheta1[i]/180.)*TMath::Pi();
    //if(trackTheta1[i] > 90 && trackTheta1[i] < 270){
      //theta = ((trackTheta1[i] - 180)/180.)*TMath::Pi();
    //}

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
  }

  printf("======================================================================== \n");
  const int dimX2 = 6;
  double trackY2[dimX2] = {5.317608,4.664446,5.196628,8.166169,9.238523,9.407446};
  double trackErrY2[dimX2] = {0.345687,0.286275,0.540735,0.373195,0.355325,0.396841};
  double trackTheta2[dimX2] = {0,40,80,160,200,240};

  TF1 *funcTrack2[dimX2];
  for(int i = 0;i < dimX2;i++){
    trackY2[i] = trackY2[i] - 7.;
    theta = (trackTheta2[i]/180.)*TMath::Pi();
    //if(trackTheta2[i] > 90 && trackTheta2[i] < 270){
      //theta = ((trackTheta2[i] - 180)/180.)*TMath::Pi();
    //}

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
  }

  TCanvas *canvasTrack = new TCanvas("canvasTrack","canvasTrack",20,20,600,600);
  histGrid -> Draw();
  for(int i = 0;i < dimX1;i++){funcTrack1[i] -> Draw("same");}
  for(int i = 0;i < dimX2;i++){funcTrack2[i] -> Draw("same");}
}
