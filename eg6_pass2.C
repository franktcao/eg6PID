#define eg6_pass2_cxx
#include "eg6_pass2.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <iomanip>
#include <time.h>
#include <TVector3.h>
#include <TMath.h>
#include "TCutG.h"
#include <TLine.h>
#include <TLorentzVector.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <TGraphErrors.h>
#include "TLegend.h"
#include "SpecsFID.hh"


SpecsFID FID;
SpecsGEO GEO;

using namespace std;

int *Decimal_Binary_Converter(int);
int Find_HWP(int);
int Find_helicity(int, int);
float PT2_H_Find(TLorentzVector, TLorentzVector);
float T_H_Find(TLorentzVector, TLorentzVector);
float PHI_Find(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
bool DC_Mid_Cut(float, float);
bool IC_theta_En_Cut (float, float);
float CalcTwoPhotonInvarMass(TLorentzVector, TLorentzVector);
float Edge2(float);
float Edge1(float);
bool CC_Fid(float, float);
int basics_vcrpl(float*, float*, float*, float*, float*);
bool IC_DC_shadow_fiducial (float , float);
bool DCff_e(float, float, int);
bool IC_Not_Hot_Channel(float,float);
bool ICFiducialCut(float, float);
TVector3 EC_XYZ_UVW(TVector3);
double EG6ev2file(const int runno,const int envo);
double EG6ECsampcorr(const int sector,const double rf150);
double EG6ECsampcorr(const int sector,const int runno,const int evno);
double EG6VertexCorrCLAS(const int runno, const double vz, const double cx,const double cy,const double cz);
double EG6VertexCorrCLAS(const int runno, const double vz, const double theta, const double phi);
float dedx_cal(float p, int charge, float mass);


// new fundtion for helicity
int    GetHelicity(const int h10hel,const int runno);
int    GetHWP(const int runno);
float Find_Error(float x, float y);

// RTPC accptance
bool AcceptRTPC(float vz, float theta, float phi);

void eg6_pass2::Loop()
 {
  gROOT->Reset();
  gStyle->SetOptStat(kFALSE);

//  gStyle->SetLineWidth(1);
//  gStyle->SetTextSize(0.9);
 //gStyle->SetPalette(55);
 gStyle->SetLabelSize(0.03,"xyz"); // size of axis value font 
 gStyle->SetTitleSize(0.035,"xyz"); // size of axis title font 
 gStyle->SetTitleFont(22,"xyz"); // font option 
 gStyle->SetLabelFont(22,"xyz"); 
 gStyle->SetTitleOffset(1.2,"y");
 gStyle->SetCanvasBorderMode(0); 
 gStyle->SetCanvasBorderSize(0); 
 gStyle->SetPadBottomMargin(0.16); //margins... 
 gStyle->SetPadTopMargin(0.16); 
 gStyle->SetPadLeftMargin(0.16); 
 gStyle->SetPadRightMargin(0.16); 
 gStyle->SetFrameBorderMode(0); 
 gStyle->SetPaperSize(20,24); 
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.06,"xy");


  TFile *skim_file = new  TFile("eg6_pass2_all_p_el_700.root", "Recreate");
  TTree *tree = new TTree("TT","the root file");
 
  TLine *l = new TLine();
  l->SetLineStyle(7);
  l->SetLineWidth(4);
  l->SetLineColor(2);  

  // constants //
  float TODEG = 180./TMath::Pi();
  float TORAD = TMath::Pi()/180.;
  float Eb= 6.064;       //in Gev
  float M_4He= 3.727;    //in GeV   
  float M_p= 0.93827;      //in GeV
  float M_pi= 0.13957;    // mass of pi(+/-) in GeV
  float M_e= 0.000511;
  float LightSpeed = 29.9792458;	// cm per ns 

  int n_run = 450;
  float N_p_Run[n_run];
  float N_e_Run[n_run];
  for(int i=0; i<n_run; i++){
     N_p_Run[i] = 0.0;
     N_e_Run[i] = 0.0;
    }

/////////////////// find the helicity signals for all the electrons //////////////////

  float N_pos_e[n_run];
  float N_neg_e[n_run];
 for(int i=0; i<n_run; i++)
   {
     N_pos_e[i] = 0.0;
     N_neg_e[i] = 0.0;
    }



  TVector3 EC_XYZ;       //Hit position in Electromagnetic Calorimeter  (XYZ coordinates)
  TVector3 EC_UVW;        //Hit position in Electromagnetic Calorimeter (UVW coordinates)

  TLorentzVector  InEl4Vector;     //incident electron 4 vector
  TLorentzVector  Target4Vector;   // rest proton target
  TLorentzVector  GammaStar4Vector;

  int event_id;
  int event_helicity;
  int event_HWP;
  int RunNumber;

  // electron variables //
  int n_e_ec; 
  TLorentzVector Elec4Vector[40];
  float nu[40], Q2[40], xB[40], yy[40], W_4He[40], W_p[40];
  float El_Px[40], El_Py[40], El_Pz[40], El_E[40], El_Phi[40], El_Theta[40], El_P[40];
  float El_z[40];
  int El_s[40];
  int El_nphe[40];
 
  // FX variabels for IC fiducial cut for selecting the electron//
  int  e_s;
  float e_ic_x, e_ic_y;

  // photon variables
  int n_phot_EC;
  int n_phot_IC;
  int n_photon; 
  float Ph_EC_Px[40], Ph_EC_Py[40], Ph_EC_Pz[40], Ph_EC_P[40], Ph_EC_E[40], Ph_EC_Theta[40], Ph_EC_Phi[40];
  float Ph_IC_Px[40], Ph_IC_Py[40], Ph_IC_Pz[40], Ph_IC_P[40], Ph_IC_E[40], Ph_IC_Theta[40], Ph_IC_Phi[40];
  float Ph_Px[40], Ph_Py[40], Ph_Pz[40], Ph_P[40], Ph_E[40], Ph_Theta[40], Ph_Phi[40];
  TLorentzVector Phot4VectEC[40]; // EC photons 4 vector
  TLorentzVector Phot4VectIC[40];     // IC photon 4 vector
  TLorentzVector Phot4VectIC_before[40];
  TLorentzVector Phot4Vect[40];       // photon 4 vecotr
  TVector3 ICPhot3Vect[40];
    
  /// pi0 variables
  int n_pi0;
  TLorentzVector Pi04Vect[40];      
  float Pi0_Px[40], Pi0_Py[40], Pi0_Pz[40], Pi0_E[40], Pi0_Phi[40], Pi0_Theta[40], Pi0_P[40], Pi0_M[40];
  float Pi0_InvMass2photECEC=0.;
  float Pi0_InvMass2photECIC=0.;
  float Pi0_InvMass2photICIC=0.;

  /// proton variabels
  int n_proton;
  int dc_count;
  float DC_IC_stat[40];
  TLorentzVector Proton4Vector[40];
  float Proton_z[40], Proton_Px[40], Proton_Py[40], Proton_Pz[40], Proton_P[40], Proton_E[40], Proton_Theta[40], Proton_Phi[40], Proton_M[40];

  /// Helium variables
  int n_tpc, npd_track_tpc[40], npd_event_tpc[40], pid_tpc[40], index_tpc[40], npts_tpc[40], fiterr_tpc[40], tothits_tpc[40], bonus_bits_tpc[40];
  float sdist_tpc[40], edist_tpc[40], r0_tpc[40], x2_tpc[40];
  float x_tpc[40], y_tpc[40], z_tpc[40], x_start_tpc[40], y_start_tpc[40], z_start_tpc[40], x_end_tpc[40], y_end_tpc[40], z_end_tpc[40];
  float px_tpc[40], py_tpc[40], pz_tpc[40], ptot_tpc[40], phi_tpc[40], theta_tpc[40], dedx_tpc[40];
  float charge_tpc[40], qtot_tpc[40], vtl_tpc[40], dca_tpc[40];

  // correcter tpc bank
  int rtpc_id1_c[40], rtpc_id2_c[40],rtpc_id3_c[40],rtpc_id4_c[40],rtpc_id5_c[40];
  float rtpc_p1_c[40], rtpc_p2_c[40],rtpc_p3_c[40],rtpc_p4_c[40],rtpc_p5_c[40];
  float rtpc_poverq_c[40], rtpc_dedx_c[40], rtpc_dedx2_c[40], rtpc_theta_c[40], rtpc_phi_c[40], rtpc_vz_c[40], rtpc_bad_c[40], rtpc_gcpb_c[40];



  tree->Branch("event_id",&event_id, "event_id/I");
  tree->Branch("event_helicity",&event_helicity, "event_helicity/I");
  tree->Branch("event_HWP",&event_HWP, "event_HWP/I");
  tree->Branch("RunNumber",&RunNumber, "RunNumber/I");
  // electron branches
  tree->Branch("n_e_ec",&n_e_ec, "n_e_ec/I");
  tree->Branch("El_Px",El_Px,"El_Px[n_e_ec]/F");
  tree->Branch("El_Py",El_Py,"El_Py[n_e_ec]/F");
  tree->Branch("El_Pz",El_Pz,"El_Pz[n_e_ec]/F");
  tree->Branch("El_P",El_P,"El_P[n_e_ec]/F");
  tree->Branch("El_E",El_E,"El_E[n_e_ec]/F");
  tree->Branch("El_Phi",El_Phi,"El_Phi[n_e_ec]/F");
  tree->Branch("El_Theta",El_Theta,"El_Theta[n_e_ec]/F");
  tree->Branch("El_z",El_z,"El_z[n_e_ec]/F");
  tree->Branch("El_s",El_s,"El_s[n_e_ec]/I");
  tree->Branch("El_nphe",El_nphe,"El_nphe[n_e_ec]/I");

  tree->Branch("nu", nu, "nu[n_e_ec]/F");
  tree->Branch("Q2", Q2, "Q2[n_e_ec]/F");
  tree->Branch("xB", xB, "xB[n_e_ec]/F");
  tree->Branch("yy", yy, "yy[n_e_ec]/F");
  tree->Branch("W_4He", W_4He, "W_4He[n_e_ec]/F");
  tree->Branch("W_p", W_p, "W_p[n_e_ec]/F");

   // Photon EC
  tree->Branch("n_phot_EC",&n_phot_EC,"n_phot_EC/I");
  tree->Branch("Ph_EC_Px",Ph_EC_Px,"Ph_EC_Px[n_phot_EC]/F");
  tree->Branch("Ph_EC_Py",Ph_EC_Py,"Ph_EC_Py[n_phot_EC]/F");
  tree->Branch("Ph_EC_Pz",Ph_EC_Pz,"Ph_EC_Pz[n_phot_EC]/F");
  tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[n_phot_EC]/F");
  tree->Branch("Ph_EC_E",Ph_EC_E,"Ph_EC_E[n_phot_EC]/F");
  tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[n_phot_EC]/F");
  tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[n_phot_EC]/F");

   // Photon IC
  tree->Branch("n_phot_IC",&n_phot_IC,"n_phot_IC/I");
  tree->Branch("Ph_IC_Px",Ph_IC_Px,"Ph_IC_Px[n_phot_IC]/F");
  tree->Branch("Ph_IC_Py",Ph_IC_Py,"Ph_IC_Py[n_phot_IC]/F");
  tree->Branch("Ph_IC_Pz",Ph_IC_Pz,"Ph_IC_Pz[n_phot_IC]/F");
  tree->Branch("Ph_IC_P",Ph_IC_P,"Ph_IC_P[n_phot_IC]/F");
  tree->Branch("Ph_IC_E",Ph_IC_E,"Ph_IC_E[n_phot_IC]/F");
  tree->Branch("Ph_IC_Theta",Ph_IC_Theta,"Ph_IC_Theta[n_phot_IC]/F");
  tree->Branch("Ph_IC_Phi",Ph_IC_Phi,"Ph_IC_Phi[n_phot_IC]/F");

  // EC+IC variables
  tree->Branch("n_photon",&n_photon,"n_photon/I"); 
  tree->Branch("Ph_Px",Ph_Px,"Ph_Px[n_photon]/F");
  tree->Branch("Ph_Py",Ph_Py,"Ph_Py[n_photon]/F");
  tree->Branch("Ph_Pz",Ph_Pz,"Ph_Pz[n_photon]/F");
  tree->Branch("Ph_P",Ph_P,"Ph_P[n_photon]/F");
  tree->Branch("Ph_E",Ph_E,"Ph_E[n_photon]/F");
  tree->Branch("Ph_Theta",Ph_Theta,"Ph_Theta[n_photon]/F");
  tree->Branch("Ph_Phi",Ph_Phi,"Ph_Phi[n_photon]/F");

  // pi0  branches
  tree->Branch("n_pi0",&n_pi0,"n_pi0/I");
  tree->Branch("Pi0_Px",Pi0_Px,"Pi0_Px[n_pi0]/F");
  tree->Branch("Pi0_Py",Pi0_Py,"Pi0_Py[n_pi0]/F");
  tree->Branch("Pi0_Pz",Pi0_Pz,"Pi0_Pz[n_pi0]/F");
  tree->Branch("Pi0_P",Pi0_P,"Pi0_P[n_pi0]/F");
  tree->Branch("Pi0_E",Pi0_E,"Pi0_E[n_pi0]/F");
  tree->Branch("Pi0_Phi",Pi0_Phi,"Pi0_Phi[n_pi0]/F");
  tree->Branch("Pi0_Theta",Pi0_Theta,"Pi0_Theta[n_pi0]/F");
  tree->Branch("Pi0_M",Pi0_M,"Pi0_M[n_pi0]/F");

  tree->Branch("Pi0_InvMass2photECEC",&Pi0_InvMass2photECEC,"Pi0_InvMass2photECEC/F");
  tree->Branch("Pi0_InvMass2photICIC",&Pi0_InvMass2photICIC,"Pi0_InvMass2photICIC/F");
  tree->Branch("Pi0_InvMass2photECIC",&Pi0_InvMass2photECIC,"Pi0_InvMass2photECIC/F");

   // proton branches
  tree->Branch("n_proton",&n_proton, "n_proton/I");
  tree->Branch("Proton_z",Proton_z,"Proton_z[n_proton]/F");
  tree->Branch("Proton_Px",Proton_Px,"Proton_Px[n_proton]/F");
  tree->Branch("Proton_Py",Proton_Py,"Proton_Py[n_proton]/F");
  tree->Branch("Proton_Pz",Proton_Pz,"Proton_Pz[n_proton]/F");
  tree->Branch("Proton_P",Proton_P,"Proton_P[n_proton]/F");
  tree->Branch("Proton_E",Proton_E,"Proton_E[n_proton]/F");
  tree->Branch("Proton_Phi",Proton_Phi,"Proton_Phi[n_proton]/F");
  tree->Branch("Proton_Theta",Proton_Theta,"Proton_Theta[n_proton]/F");
  tree->Branch("Proton_M",Proton_M,"Proton_M[n_proton]/F");

  tree->Branch("dc_count",&dc_count, "dc_count/I");
  tree->Branch("DC_IC_stat",DC_IC_stat,"DC_IC_stat[dc_count]/F");

  /// rtpc branches
  tree->Branch("n_tpc",&n_tpc, "n_tpc/I");
  tree->Branch("npd_track_tpc",npd_track_tpc,"npd_track_tpc[n_tpc]/I");
  tree->Branch("npd_event_tpc",npd_event_tpc,"npd_event_tpc[n_tpc]/I");
  tree->Branch("pid_tpc",pid_tpc,"pid_tpc[n_tpc]/I");
  tree->Branch("index_tpc",index_tpc,"index_tpc[n_tpc]/I");
  tree->Branch("npts_tpc",npts_tpc,"npts_tpc[n_tpc]/I");
  tree->Branch("fiterr_tpc",fiterr_tpc,"fiterr_tpc[n_tpc]/I");
  tree->Branch("tothits_tpc",tothits_tpc,"tothits_tpc[n_tpc]/I");
  tree->Branch("bonus_bits_tpc",bonus_bits_tpc,"bonus_bits_tpc[n_tpc]/I");
  tree->Branch("sdist_tpc",sdist_tpc,"sdist_tpc[n_tpc]/F");
  tree->Branch("edist_tpc",edist_tpc,"edist_tpc[n_tpc]/F");
  tree->Branch("r0_tpc",r0_tpc,"r0_tpc[n_tpc]/F");
  tree->Branch("x2_tpc",x2_tpc,"x2_tpc[n_tpc]/F");
  tree->Branch("x_tpc",x_tpc,"x_tpc[n_tpc]/F");
  tree->Branch("y_tpc",y_tpc,"y_tpc[n_tpc]/F");
  tree->Branch("z_tpc",z_tpc,"z_tpc[n_tpc]/F");
  tree->Branch("x_start_tpc",x_start_tpc,"x_start_tpc[n_tpc]/F");
  tree->Branch("y_start_tpc",y_start_tpc,"y_start_tpc[n_tpc]/F");
  tree->Branch("z_start_tpc",z_start_tpc,"z_start_tpc[n_tpc]/F");
  tree->Branch("x_end_tpc",x_end_tpc,"x_end_tpc[n_tpc]/F");
  tree->Branch("y_end_tpc",y_end_tpc,"y_end_tpc[n_tpc]/F");
  tree->Branch("z_end_tpc",z_end_tpc,"z_end_tpc[n_tpc]/F");
  tree->Branch("px_tpc",px_tpc,"px_tpc[n_tpc]/F");
  tree->Branch("py_tpc",py_tpc,"py_tpc[n_tpc]/F");
  tree->Branch("pz_tpc",pz_tpc,"pz_tpc[n_tpc]/F");
  tree->Branch("ptot_tpc",ptot_tpc,"ptot_tpc[n_tpc]/F");
  tree->Branch("phi_tpc",phi_tpc,"phi_tpc[n_tpc]/F");
  tree->Branch("theta_tpc",theta_tpc,"theta_tpc[n_tpc]/F");
  tree->Branch("dedx_tpc",dedx_tpc,"dedx_tpc[n_tpc]/F");
  tree->Branch("charge_tpc",charge_tpc,"charge_tpc[n_tpc]/F");
  tree->Branch("dca_tpc",dca_tpc,"dca_tpc[n_tpc]/F");
  tree->Branch("qtot_tpc",qtot_tpc,"qtot_tpc[n_tpc]/F");
  tree->Branch("vtl_tpc",vtl_tpc,"vtl_tpc[n_tpc]/F");


  // corrected bank of tpc
  tree->Branch("rtpc_id1_c",rtpc_id1_c,"rtpc_id1_c[n_tpc]/I");
  tree->Branch("rtpc_id2_c",rtpc_id2_c,"rtpc_id2_c[n_tpc]/I");
  tree->Branch("rtpc_id3_c",rtpc_id3_c,"rtpc_id3_c[n_tpc]/I");
  tree->Branch("rtpc_id4_c",rtpc_id4_c,"rtpc_id4_c[n_tpc]/I");
  tree->Branch("rtpc_id5_c",rtpc_id5_c,"rtpc_id5_c[n_tpc]/I");
  tree->Branch("rtpc_p1_c",rtpc_p1_c,"rtpc_p1_c[n_tpc]/F");
  tree->Branch("rtpc_p2_c",rtpc_p2_c,"rtpc_p2_c[n_tpc]/F");
  tree->Branch("rtpc_p3_c",rtpc_p3_c,"rtpc_p3_c[n_tpc]/F");
  tree->Branch("rtpc_p4_c",rtpc_p4_c,"rtpc_p4_c[n_tpc]/F");
  tree->Branch("rtpc_p5_c",rtpc_p5_c,"rtpc_p5_c[n_tpc]/F");
  tree->Branch("rtpc_poverq_c",rtpc_poverq_c,"rtpc_poverq_c[n_tpc]/F");
  tree->Branch("rtpc_dedx_c",rtpc_dedx_c,"rtpc_dedx_c[n_tpc]/F");
  tree->Branch("rtpc_dedx2_c",rtpc_dedx2_c,"rtpc_dedx2_c[n_tpc]/F");
  tree->Branch("rtpc_theta_c",rtpc_theta_c,"rtpc_theta_c[n_tpc]/F");
  tree->Branch("rtpc_phi_c",rtpc_phi_c,"rtpc_phi_c[n_tpc]/F");
  tree->Branch("rtpc_vz_c",rtpc_vz_c,"rtpc_vz_c[n_tpc]/F");
  tree->Branch("rtpc_bad_c",rtpc_bad_c,"rtpc_bad_c[n_tpc]/F");
  tree->Branch("rtpc_gcpb_c",rtpc_gcpb_c,"rtpc_gcpb_c[n_tpc]/F");


  // Bethe-Bloche formula
/*
  // proton
  TF1 *fp = new TF1("fp","(0.307075*1.03*[0]*[0]*66/126.787)*((0.5 *(pow(x*[0],2)+[1]*[1])/pow(x*[0],2))*log(pow(2*[2]*(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))/(1-(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))),2)/(pow(99.794E-6,2)*(1+(2*[2]/[1])/sqrt(1- pow(x*[0],2)/(pow(x*[0],2)+[1]*[1])) +pow([2]/[1],2))))-1)",50, 250);
  fp->SetLineColor(28);
  fp->SetLineWidth(2);
  fp->SetLineStyle(1);
  fp->SetParameters(1, 938, 0.511);
  // deutrium
   TF1 *fd = new TF1("fd","(0.307075*1.03*[0]*[0]*66/126.787)*((0.5 *(pow(x*[0],2)+[1]*[1])/pow(x*[0],2))*log(pow(2*[2]*(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))/(1-(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))),2)/(pow(99.794E-6,2)*(1+(2*[2]/[1])/sqrt(1- pow(x*[0],2)/(pow(x*[0],2)+[1]*[1])) +pow([2]/[1],2))))-1)",50, 250);
  fd->SetLineColor(6);
  fd->SetLineWidth(2);
  fd->SetLineStyle(1);
  fd->SetParameters(1, 1875.6, 0.511);
  // 4He
  TF1 *f4He = new TF1("f4He","(0.307075*1.03*[0]*[0]*66/126.787)*((0.5 *(pow(x*[0],2)+[1]*[1])/pow(x*[0],2))*log(pow(2*[2]*(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))/(1-(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))),2)/(pow(99.794E-6,2)*(1+(2*[2]/[1])/sqrt(1- pow(x*[0],2)/(pow(x*[0],2)+[1]*[1])) +pow([2]/[1],2))))-1)",50, 250);
  f4He->SetLineColor(2);
  f4He->SetLineWidth(2);
  f4He->SetLineStyle(1);
  f4He->SetParameters(2, 3727.38, 0.511);
  // 3He
 TF1 *f3He = new TF1("f3He","(0.307075*1.03*[0]*[0]*66/126.787)*((0.5 *(pow(x*[0],2)+[1]*[1])/pow(x*[0],2))*log(pow(2*[2]*(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))/(1-(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))),2)/(pow(99.794E-6,2)*(1+(2*[2]/[1])/sqrt(1- pow(x*[0],2)/(pow(x*[0],2)+[1]*[1])) +pow([2]/[1],2))))-1)",50, 250);
  f3He->SetLineColor(kBlack);
  f3He->SetLineWidth(2);
  f3He->SetLineStyle(1);
  f3He->SetParameters(2, 2808.921, 0.511);

 //3H
  TF1 *f3H = new TF1("f3H","(0.307075*1.03*[0]*[0]*66/126.787)*((0.5 *(pow(x*[0],2)+[1]*[1])/pow(x*[0],2))*log(pow(2*[2]*(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))/(1-(pow(x*[0],2)/(pow(x*[0],2)+[1]*[1]))),2)/(po     w(99.794E-6,2)*(1+(2*[2]/[1])/sqrt(1- pow(x*[0],2)/(pow(x*[0],2)+[1]*[1])) +pow([2]/[1],2))))-1)",50, 250);
  f3H->SetLineColor(4);
  f3H->SetLineWidth(2);
  f3H->SetLineStyle(1);
  f3H->SetParameters(1, 2808.921, 0.511);
*/



 TH1D *h_A_RunNum = new TH1D("h_A_RunNum","A_{C} vs. Run number", 450, 61500, 61945);
       h_A_RunNum->SetLineWidth(2);
       h_A_RunNum->SetXTitle("Run Number");
       h_A_RunNum->SetYTitle("A_{C}");


  int N = 250;
  float mom[N] , mom_err[N];
  float DEDx[N];
  float DEDx_err[N];
  for(int i =0; i<N; i++){
     mom[i] = 50 + i;
     mom_err[i] = .0;
     DEDx_err[i] = .0;   }

  // proton
  for(int i =0; i<N; i++)  DEDx[i] = dedx_cal(mom[i], 1, 938.0);
  TGraphErrors *fp = new TGraphErrors(N, mom, DEDx, mom_err, DEDx_err);
  fp->SetLineColor(28);
  fp->SetLineWidth(3);
  fp->SetLineStyle(1);

  // deutrium
  for(int i =0; i<N; i++)  DEDx[i] = dedx_cal(mom[i], 1, 1875.6);
  TGraphErrors *fd = new TGraphErrors(N, mom, DEDx, mom_err, DEDx_err);
  fd->SetLineColor(6);
  fd->SetLineWidth(3);
  fd->SetLineStyle(1);

 //3H
  for(int i =0; i<N; i++)  DEDx[i] = dedx_cal(mom[i], 1, 2808.921);
  TGraphErrors *f3H = new TGraphErrors(N, mom, DEDx, mom_err, DEDx_err);
  f3H->SetLineColor(4);
  f3H->SetLineWidth(3);
  f3H->SetLineStyle(1);

  // 3He
  for(int i =0; i<N; i++)  DEDx[i] = dedx_cal(mom[i], 2, 2808.921);
  TGraphErrors *f3He = new TGraphErrors(N, mom, DEDx, mom_err, DEDx_err);
  f3He->SetLineColor(kBlack);
  f3He->SetLineWidth(3);
  f3He->SetLineStyle(1);

  // 4He
  for(int i =0; i<N; i++)  DEDx[i] = dedx_cal(mom[i], 2, 3727.38);
  TGraphErrors *f4He = new TGraphErrors(N, mom, DEDx, mom_err, DEDx_err);
  f4He->SetLineColor(2);
  f4He->SetLineWidth(3);
  f4He->SetLineStyle(1);


  //** electron distributions **//
  TH2D *h_EC_XY_el_1 = new TH2D("h_EC_XY_el_1","Negative charged particles BEFORE UVW cuts", 200, -400, 400, 200, -400, 400);
        h_EC_XY_el_1->SetXTitle("X [cm]"); 
        h_EC_XY_el_1->SetYTitle("Y [cm]");
  TH2D *h_EC_XY_el_2 = new TH2D("h_EC_XY_el_2","Negative charged particles AFTER UVW cuts", 200, -400, 400, 200, -400, 400);
        h_EC_XY_el_2->SetXTitle("X [cm]"); 
        h_EC_XY_el_2->SetYTitle("Y [cm]");

  TH2D *h_DC_el_1 = new TH2D("h_DC_el_1","Negative charged particles BEFORE DC fiducial cuts", 200, -80, 80, 200, -80, 80);
        h_DC_el_1->SetXTitle("X [cm]"); 
        h_DC_el_1->SetYTitle("Y [cm]");
  TH2D *h_DC_el_2 = new TH2D("h_DC_el_2","Negative charged particles AFTER DC fiducial cuts", 200, -80, 80, 200, -80, 80);
        h_DC_el_2->SetXTitle("X [cm]"); 
        h_DC_el_2->SetYTitle("Y [cm]");

  TH2D *h_e_phi_theta_1 = new TH2D("h_e_phi_theta_1", "#phi_{e^{-}} vs. #theta_{e^{-}} before fiducial cuts", 100, 0, 70, 200, -180, 180);
        h_e_phi_theta_1->SetXTitle("#theta_{e^{-}} [deg.]"); 
        h_e_phi_theta_1->SetYTitle("#phi_{e^{-}} [deg.]"); 
  TH2D *h_e_phi_theta_2 = new TH2D("h_e_phi_theta_2", "#phi_{e^{-}} vs. #theta_{e^{-}} after uvw fiducial cuts", 100, 0, 70, 200, -180, 180);
        h_e_phi_theta_2->SetXTitle("#theta_{e^{-}} [deg.]");  
        h_e_phi_theta_2->SetYTitle("#phi_{e^{-}} [deg.]"); 
  TH2D *h_e_phi_theta_3 = new TH2D("h_e_phi_theta_3", "#phi_{e^{-}} vs. #theta_{e^{-}} after uvw and IC and DC fiducial cuts", 100, 0, 70, 200, -180, 180);
        h_e_phi_theta_3->SetXTitle("#theta_{e^{-}} [deg.]");   
        h_e_phi_theta_3->SetYTitle("#phi_{e^{-}} [deg.]"); 
  TH2D *h_e_phi_theta_4 = new TH2D("h_e_phi_theta_4", "#phi_{e^{-}} vs. #theta_{e^{-}} for the selected e^{-} after all the cuts", 100, 0, 70, 200, -180, 180);
        h_e_phi_theta_4->SetXTitle("#theta_{e^{-}} [deg.]");    
        h_e_phi_theta_4->SetYTitle("#phi_{e^{-}} [deg.]");

 
  TH2D *h_el_DC_IC_XY_1 = new TH2D("h_el_DC_IC_XY_1","Negative charged particles BEFORE IC shadow cut", 300, -65, 65, 300, -65, 65);
        h_el_DC_IC_XY_1->SetXTitle("X [cm]"); 
        h_el_DC_IC_XY_1->SetYTitle("Y [cm]");
  TH2D *h_el_DC_IC_XY_2 = new TH2D("h_el_DC_IC_XY_2","Negative charged particles AFTER IC shadow cut", 300, -65, 65, 300, -65, 65);
        h_el_DC_IC_XY_2->SetXTitle("X [cm]");
        h_el_DC_IC_XY_2->SetYTitle("Y [cm]");


  TH2D *h_CC_el_phi_theta_1 = new TH2D("h_CC_el_phi_theta_1", "Negative charged particles BEFORE CC fiducial cuts", 100, 0, 50, 200, -40, 40);
        h_CC_el_phi_theta_1->SetXTitle("#theta_{e^{-}} [deg.]");   
        h_CC_el_phi_theta_1->SetYTitle("#phi_{e^{-}} [deg.]");
  TH2D *h_CC_el_phi_theta_2 = new TH2D("h_CC_el_phi_theta_2", "Negative charged particles AFTER CC fiducial cuts", 100, 0, 50, 200, -40, 40);
        h_CC_el_phi_theta_2->SetXTitle("#theta_{e^{-}} [deg.]");   
        h_CC_el_phi_theta_2->SetYTitle("#phi_{e^{-}} [deg.]");


  TH2D *h_etot_p = new TH2D("h_etot_p", "e_{tot}/p vs. p)", 200, 0.7, 5, 200, 0, 1);
        h_etot_p->SetXTitle("p [GeV]");
        h_etot_p->SetYTitle("e_{tot}/p ");
  TH2D *h_etot_p_3 = new TH2D("h_etot_p_3", "e_{tot}/p vs. p", 200, 0, 5, 200, 0, 1);
        h_etot_p_3->SetXTitle("p_{e} [GeV]");
        h_etot_p_3->SetYTitle("e_{tot}/p ");
  TH2D *h_etot_p_2 = new TH2D("h_etot_p_2", "Rejected e^{-} : (e_{tot}/p vs. p)", 200, 0, 5, 200, 0, 1);
        h_etot_p_2->SetXTitle("p_{e} [GeV]");
        h_etot_p_2->SetYTitle("e_{tot}/p ");

  TH2D *h_eco_eci = new TH2D("h_eco_eci","EC_{outer} vs. EC_{inner}", 200, 0, 0.6, 200, 0, 0.5);
        h_eco_eci->SetXTitle("EC_{in} [GeV]"); 
        h_eco_eci->SetYTitle("EC_{out} [GeV]");

  TH1D *h_W_4He = new TH1D("h_W_4He","W Distribution for 4He",300, 3.5,8);
        h_W_4He->SetXTitle("W_{4He} [GeV]");
        h_W_4He->SetLineWidth(2);
  TH1D *h_W_p = new TH1D("h_W_p","W Distribution for p",300, 0, 4.5);
        h_W_p->SetXTitle("W_{p} [GeV]");
        h_W_p->SetLineWidth(2);
  TH1D *h_Q2 = new TH1D("h_Q2","Q2 Distribution",300,0,7);
        h_Q2->SetXTitle("Q^{2} [GeV^{2}]");
        h_Q2->SetLineWidth(2);
  TH1D *h_xB = new TH1D("h_xB","xB Distribution", 300, 0, 2);
  TH1D *h_yy = new TH1D("h_yy","y Distribution",300,0,1);
  TH1D *h_nu = new TH1D("h_nu","nu Distribution",300, 0, 7);
  TH2D *h_Q2_xB= new TH2D("h_Q2_xB", "Q2 vs. xB Distribution", 300, 0, 1, 300, 0, 5);
  TH1D *h_z_e = new TH1D("h_z_e","z_{vertex} for negative charged particles", 300, -95, -25);
        h_z_e->SetXTitle("z_{vertex} [cm]");
        h_z_e->SetLineWidth(2);


  //** photoelectrons histograms **//
  TH1D *h_e_1_nphe = new TH1D("h_e_1_nphe","Nphe", 250, 0, 250); 
        h_e_1_nphe->SetLineColor(kBlack);
        h_e_1_nphe->SetLineWidth(2);
        h_e_1_nphe->SetXTitle("10*nphe");
  TH1D *h_e_2_nphe = new TH1D("h_e_2_nphe","Nphe", 250, 0, 250); 
        h_e_2_nphe->SetLineWidth(2);
        h_e_2_nphe->SetLineColor(kBlue);
  TH1D *h_e_3_nphe = new TH1D("h_e_3_nphe","Nphe", 250, 0, 250); 
        h_e_3_nphe->SetLineWidth(2);
        h_e_3_nphe->SetLineColor(kRed);

  TH1D *h_e_nphe_1[7];
  TH1D *h_e_nphe_2[7];
  TH1D *h_e_nphe_3[7];
  TH2D *h_ece_over_p_nphe[7];
  TH2D *h_ch2cc_nphe[7];
  for(int i = 0; i<7; i++)
  {
      h_e_nphe_1[i] = new TH1D(Form("h_e_nphe_1[%u]",i),Form("e^{-}: nphe : sector %u",i), 250, 0, 250);
       h_e_nphe_1[i]->SetLineColor(kBlack);
      h_e_nphe_2[i] = new TH1D(Form("h_e_nphe_2[%u]",i),Form("e^{-}: nphe : sector %u",i), 250, 0, 250);
           h_e_nphe_2[i]->SetLineColor(kBlue);
      h_e_nphe_3[i] = new TH1D(Form("h_e_nphe_3[%u]",i),Form("e^{-}: nphe : sector %u",i), 250, 0, 250);
           h_e_nphe_3[i]->SetLineColor(kRed);

      h_ece_over_p_nphe[i] = new TH2D(Form("h_ece_over_p_nphe[%u]",i),
                                Form("e^{-}: etot/p vs.  nphe: sector %u",i), 250, 0, 250, 200, 0., 0.6 );
           h_ece_over_p_nphe[i]->SetXTitle("nphe");
           h_ece_over_p_nphe[i]->SetYTitle("etot/p");
      h_ch2cc_nphe[i] = new TH2D(Form("h_ch2cc_nphe[%u]",i),
                                Form("e^{-}: CC_{chi^{2}} vs.  nphe: sector %u",i), 250, 0, 250, 200, -0.05, 0.3 );
           h_ch2cc_nphe[i]->SetXTitle("nphe");
         1  h_ch2cc_nphe[i]->SetYTitle("cc_chi2 [rad]");
      }


  //**** EC coordinates histograms ***//
  TH1D *h_EC_el_U  = new TH1D("h_EC_el_U", "Electrons U distribution in the ECs", 150, 0, 450);
        h_EC_el_U->SetLineWidth(2);
        h_EC_el_U->SetXTitle("U [cm]");
  TH1D *h_EC_el_V  = new TH1D("h_EC_el_V", "Electrons V distribution in the ECs", 150, 0, 450);
        h_EC_el_V->SetLineWidth(2);
        h_EC_el_V->SetXTitle("V [cm]");
  TH1D *h_EC_el_W  = new TH1D("h_EC_el_W", "Electrons W distribution in the ECs", 150, 0, 450);
        h_EC_el_W->SetLineWidth(2);
        h_EC_el_W->SetXTitle("W [cm]");

  TH1D *h_EC_photon_U  = new TH1D("h_EC_photon_U", "Photons U distribution in the ECs", 150, 0, 450);
        h_EC_photon_U->SetLineWidth(2);
        h_EC_photon_U->SetXTitle("U [cm]");
  TH1D *h_EC_photon_V  = new TH1D("h_EC_photon_V", "Photons V coordinate in the ECs", 150, 0, 450);
        h_EC_photon_V->SetLineWidth(2);
        h_EC_photon_V->SetXTitle("V [cm]");
  TH1D *h_EC_photon_W  = new TH1D("h_EC_photon_W", "Photons W coordinate in the ECs", 150, 0, 450);
        h_EC_photon_W->SetLineWidth(2);
        h_EC_photon_W->SetXTitle("W [cm]");


  // photon histograms ///
  TH1D *h_photon_beta = new TH1D("h_photon_beta","", 300, 0.6, 1.4);
        h_photon_beta->SetLineWidth(2);
        h_photon_beta->SetXTitle("#beta_{#gamma}");
  TH1D *h_eng_photon = new TH1D("h_eng_photon","Photons energy distribution", 300, 0, 2); 
        h_eng_photon->SetLineWidth(2);
        h_eng_photon->SetXTitle("E_{#gamma} [GeV]");
  TH1D *h_ic_time = new TH1D("h_ic_time","IC reconstructed cluster time", 200, -500, 500);
        h_ic_time->SetLineWidth(2);
        h_ic_time->SetXTitle("IC_time");
  TH1D *h_ic_theta = new TH1D("h_ic_theta","#theta_{IC} before cuts", 200, 0, 20);
        h_ic_theta->SetLineWidth(2);
        h_ic_theta->SetXTitle("#theta_{IC} [deg.]");

  TH2D *h_IC_theta_En_1 = new TH2D("h_IC_theta_En_1","IC photons: #theta_{#gamma} vs . E_{#gamma} after IC fiducial cut", 150, 0., 6. , 300, 2.5, 14.5 );
        h_IC_theta_En_1->SetYTitle("#theta_{#gamma} [deg.]");
        h_IC_theta_En_1->SetXTitle("E_{#gamma} [GeV]");
  TH2D *h_IC_theta_En_2 = new TH2D("h_IC_theta_En_2","IC photons: #theta_{#gamma} vs . E_{#gamma}", 150, 0., 6. , 300, 2.5, 14.5 );
        h_IC_theta_En_2->SetYTitle("#theta_{#gamma} [deg.]");  
        h_IC_theta_En_2->SetXTitle("E_{#gamma} [GeV]");

  TH2D *h_photon_phi_theta_1 = new TH2D("h_photon_phi_theta_1", "#phi_{#gamma} vs. #theta_{#gamma} before fiducial cuts", 100, 0, 70, 200, -180, 180);
        h_photon_phi_theta_1->SetXTitle("#theta_{#gamma} [deg.]"); 
        h_photon_phi_theta_1->SetYTitle("#phi_{#gamma} [deg.]");
  TH2D *h_photon_phi_theta_2 = new TH2D("h_photon_phi_theta_2", "#phi_{#gamma} vs. #theta_{#gamma} after all the cuts", 100, 0, 70, 200, -180, 180);
        h_photon_phi_theta_2->SetXTitle("#theta_{#gamma} [deg.]");
        h_photon_phi_theta_2->SetYTitle("#phi_{#gamma} [deg.]");

  TH2D *h_IC_XY_1 = new TH2D("h_IC_XY_1","XY projection of IC photons ", 100, -20, 20, 100, -20, 20);
        h_IC_XY_1->SetXTitle("X [cm]"); 
        h_IC_XY_1->SetYTitle("Y [cm]");
  TH2D *h_IC_XY_2 = new TH2D("h_IC_XY_2","XY projection of IC photons ", 100, -20, 20, 100, -20, 20);
        h_IC_XY_2->SetXTitle("X [cm]"); h_IC_XY_2->SetYTitle("Y [cm]");
  TH2D *h_IC_XY_center = new TH2D("h_IC_XY_center","XY projection of central hit coord. in IC after fid. cut", 100, -20, 20, 100, -20, 20);
        h_IC_XY_center->SetXTitle("X [cm]"); 
        h_IC_XY_center->SetYTitle("Y [cm]");

  TH2D *h_EC_XY_photon_1 = new TH2D("h_EC_XY_photon_1","Photons XY projection in the EC ", 200, -400, 400, 200, -400, 400);
        h_EC_XY_photon_1->SetXTitle("X [cm]"); 
        h_EC_XY_photon_1->SetYTitle("Y [cm]");
  TH2D *h_EC_XY_photon_2 = new TH2D("h_EC_XY_photon_2","Photons XY projection in the EC", 200, -400, 400, 200, -400, 400);
        h_EC_XY_photon_2->SetXTitle("X [cm]"); 
        h_EC_XY_photon_2->SetYTitle("Y [cm]");

   //// pi0 final state
  TH1D *h_pi0_ECEC = new TH1D("h_pi0_ECEC","PI0 Inv.Mass: ECEC", 200, 0, 1.0);
        h_pi0_ECEC->SetLineWidth(2);
        h_pi0_ECEC->SetXTitle("M_{#gamma #gamma} [GeV]");
  TH1D *h_pi0_ICIC = new TH1D("h_pi0_ICIC","PI0 Inv.Mass: ICIC", 200, 0, 0.5);
        h_pi0_ICIC->SetLineWidth(2);
        h_pi0_ECEC->SetXTitle("M_{#gamma #gamma} [GeV]");
  TH1D *h_pi0_ECIC = new TH1D("h_pi0_ECIC","PI0 Inv.Mass:ECIC", 200, 0, 0.7);
        h_pi0_ECIC->SetLineWidth(2);
        h_pi0_ECIC->SetXTitle("M_{#gamma #gamma} [GeV]");


   /// proton and positive partciles 
 TH1D *h_delta_z_e_prot = new TH1D("h_delta_z_e_prot"," ", 200, -7, 7);
       h_delta_z_e_prot->SetLineWidth(2);
       h_delta_z_e_prot->SetXTitle("#Delta z (= z_{e} - z_{p}) [cm]");

  TH1D *h_proton_z = new TH1D("h_proton_z","z-vertez for the positive charged particles", 150, -120, -20);
        h_proton_z->SetXTitle("z_{p} [cm]");
        h_proton_z->SetLineWidth(2);
  TH2D *h_prot_DeltaT_p = new TH2D("h_prot_DeltaT_p","#Delta T vs p for the positive particles",300,0,3.,100,-10.0,10.0);
        h_prot_DeltaT_p->SetXTitle("p (GeV)");
        h_prot_DeltaT_p->SetYTitle("#Delta T (ns)");
  TH2D *h_prot_DeltaBeta_p = new TH2D("h_prot_DeltaBeta_p","",300,0,3.,140,-0.7,0.7);
        h_prot_DeltaBeta_p->SetXTitle("p [GeV/c]"); 
        h_prot_DeltaBeta_p->SetYTitle("#Delta#beta");
 TH1D *h_prot_DeltaBeta =new TH1D("h_prot_DeltaBeta","",140,-0.1,0.1);
       h_prot_DeltaBeta->SetXTitle("#Delta#beta");
       h_prot_DeltaBeta->SetLineWidth(2);

  TH2D *h_prot_Beta_p = new TH2D("h_prot_Beta_p","#beta vs p for the positive particles",300,0,3.,100,0.,1.2);
        h_prot_Beta_p->SetXTitle("p [GeV/c]");
        h_prot_Beta_p->SetYTitle("#beta");

  TH2D *h_prot_Theta_Phi_1=new TH2D("h_prot_Theta_Phi_1","#theta vs #phi for all the positive charged particles before fiducial cuts",140,0,70,360,-180,180);
        h_prot_Theta_Phi_1->SetXTitle("#theta_{p} [Deg.]");
        h_prot_Theta_Phi_1->SetYTitle("#phi_{p} [Deg.]");
  TH2D *h_prot_Theta_Phi_2=new TH2D("h_prot_Theta_Phi_2","#theta vs #phi for all the collected protons",140,0,70,360,-180,180);
        h_prot_Theta_Phi_2->SetXTitle("#theta_{p} [Deg.]");
        h_prot_Theta_Phi_2->SetYTitle("#phi_{p} [Deg.]");

  TH2D *h_proton_DC_IC_XY_1 = new TH2D("h_proton_DC_IC_XY_1","All the positive charged particles before the fiducial cuts", 300, -90, 90, 300, -90, 90);
        h_proton_DC_IC_XY_1->SetXTitle("X [cm]");
        h_proton_DC_IC_XY_1->SetYTitle("Y [cm]");
  TH2D *h_proton_DC_IC_XY_2 = new TH2D("h_proton_DC_IC_XY_2","Positive charged particles rejected by the fiducial cuts", 300, -90, 90, 300, -90, 90);
        h_proton_DC_IC_XY_2->SetXTitle("X [cm]");
        h_proton_DC_IC_XY_2->SetYTitle("Y [cm]");
  TH2D *h_proton_DC_IC_XY_3 = new TH2D("h_proton_DC_IC_XY_3","Positive charged particles passed the fiducial cuts", 300, -90, 90, 300, -90, 90);
        h_proton_DC_IC_XY_3->SetXTitle("X [cm]");
        h_proton_DC_IC_XY_3->SetYTitle("Y [cm]");
  TH2D *h_proton_DC_IC_XY_4 = new TH2D("h_proton_DC_IC_XY_4","XY projection in DC1 for the collected protons", 300, -90, 90, 300, -90, 90);
        h_proton_DC_IC_XY_4->SetXTitle("X [cm]");
        h_proton_DC_IC_XY_4->SetYTitle("Y [cm]");

  /// 4He Histograms
  TH1D *h_rtpc_npd_l = new TH1D("h_rtpc_npd_l", "RTPC: Number of active readout pads in each track", 20, 0, 20);
        h_rtpc_npd_l->SetLineWidth(2);
  TH1D *h_rtpc_sdist_l = new TH1D("h_rtpc_sdist_l", "RTPC: sdist", 250, -6, 6);
        h_rtpc_sdist_l->SetLineWidth(2);
        h_rtpc_sdist_l->SetXTitle("sdist [mm]");
  TH1D *h_rtpc_edist_l = new TH1D("h_rtpc_edist_l", "RTPC: edist", 250, -5, 15);
        h_rtpc_edist_l->SetLineWidth(2);
        h_rtpc_edist_l->SetXTitle("edist [mm]");
  TH1D *h_rtpc_r0_l = new TH1D("h_rtpc_r0_l", "RTPC: Track's redius of  curvature", 250, -200, 200);
        h_rtpc_r0_l->SetLineWidth(2);
        h_rtpc_r0_l->SetXTitle("r_{0} [mm]");        
  TH1D *h_rtpc_X2_l = new TH1D("h_rtpc_X2_l", "RTPC: Track #chi^{2}", 200, 0, 10);
        h_rtpc_X2_l->SetLineWidth(2);
        h_rtpc_X2_l->SetXTitle("#chi^{2}");
  TH1D *h_rtpc_z_l = new TH1D("h_rtpc_z_l", "RTPC: z_{^{4}He}", 300, -170, 210); 
        h_rtpc_z_l->SetLineWidth(2);
        h_rtpc_z_l->SetXTitle("^{4}He_z [mm]");
  TH1D *h_delta_z_l = new TH1D("h_delta_z_l", "#Delta_z = z_{e} - z_{^{4}He}", 250, -50, 50);
        h_delta_z_l->SetLineWidth(2);
        h_delta_z_l->SetXTitle("#Delta_z [mm]");


  TH1D *h_rtpc_npd_r = new TH1D("h_rtpc_npd_r", "RTPC: Number of active readout pads in each track", 20, 0, 20);
        h_rtpc_npd_r->SetLineWidth(2);
        h_rtpc_npd_r->SetLineColor(kRed);
        h_rtpc_npd_r->SetMinimum(0.0);
        h_rtpc_npd_r->SetXTitle("npd_track");
  TH1D *h_rtpc_sdist_r = new TH1D("h_rtpc_sdist_r", "RTPC: sdist", 250, -6, 6);
        h_rtpc_sdist_r->SetLineWidth(2);
        h_rtpc_sdist_r->SetLineColor(kRed);
        h_rtpc_sdist_r->SetMinimum(0.0);
        h_rtpc_sdist_r->SetXTitle("sdist [mm]");
  TH1D *h_rtpc_edist_r = new TH1D("h_rtpc_edist_r", "RTPC: edist", 250, -5, 15);
        h_rtpc_edist_r->SetLineWidth(2);
        h_rtpc_edist_r->SetLineColor(kRed);
        h_rtpc_edist_r->SetMinimum(0.0);
        h_rtpc_edist_r->SetXTitle("edist [mm]");
  TH1D *h_rtpc_r0_r = new TH1D("h_rtpc_r0_r", "RTPC: Track's radius of curvature", 250, -200, 200);
        h_rtpc_r0_r->SetLineColor(kRed);
        h_rtpc_r0_r->SetLineWidth(2);
        h_rtpc_r0_r->SetMinimum(0.0);
        h_rtpc_r0_r->SetXTitle("r_{0} [mm]");
  TH1D *h_rtpc_X2_r = new TH1D("h_rtpc_X2_r", "RTPC: Track #chi^{2}", 200, 0, 10);
        h_rtpc_X2_r->SetLineColor(kRed);
        h_rtpc_X2_r->SetLineWidth(2);
        h_rtpc_X2_r->SetMinimum(0.0);
        h_rtpc_X2_r->SetXTitle("#chi^{2}");
  TH1D *h_rtpc_z_r = new TH1D("h_rtpc_z_r", "RTPC: z_{^{4}He}", 300, -170, 210);
        h_rtpc_z_r->SetLineColor(kRed);
        h_rtpc_z_r->SetLineWidth(2);
        h_rtpc_z_r->SetMinimum(0.0);
        h_rtpc_z_r->SetXTitle("^{4}He_z [mm]");
  TH1D *h_delta_z_r = new TH1D("h_delta_z_r", "#Delta z = z_{e} - z_{^{4}He}", 250, -50, 50);
        h_delta_z_r->SetLineColor(kRed);
        h_delta_z_r->SetLineWidth(2);
        h_delta_z_r->SetMinimum(0.0);
        h_delta_z_r->SetXTitle("#Delta_z [mm]");

  TH2D *h_rtpc_phi_theta_1 = new TH2D("h_rtpc_phi_theta_1", "RTPC: #phi vs. #theta BEFORE the quality cuts", 
                                       350, 10, 170, 360, 0, 360);
        h_rtpc_phi_theta_1->SetLineWidth(2);
        h_rtpc_phi_theta_1->SetXTitle("#theta [deg.]");
        h_rtpc_phi_theta_1->SetYTitle("#phi [deg.]");

  TH2D *h_rtpc_phi_theta_2 = new TH2D("h_rtpc_phi_theta_2", "RTPC: #phi vs. #theta AFTER the quality cuts",
                                       350, 10, 170, 360, 0, 360);
        h_rtpc_phi_theta_2->SetLineWidth(2);
        h_rtpc_phi_theta_2->SetXTitle("#theta [deg.]");
        h_rtpc_phi_theta_2->SetYTitle("#phi [deg.]");


  TH2D *h_dedx_p_Coh_l = new TH2D("h_dedx_p_Coh_l","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l->SetYTitle("dE/dx");
        h_dedx_p_Coh_l->SetContour(50);
  TH2D *h_dedx_p_Coh_r = new TH2D("h_dedx_p_Coh_r","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r->SetYTitle("dE/dx");
        h_dedx_p_Coh_r->SetContour(50);

  TH2D *h_dedx_p_Coh_l_0 = new TH2D("h_dedx_p_Coh_l_0","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_0->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_0->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_0->SetContour(50);
  TH2D *h_dedx_p_Coh_r_0 = new TH2D("h_dedx_p_Coh_r_0","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r_0->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_0->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_0->SetContour(50);


  TH2D *h_dedx_p_Coh_l_1 = new TH2D("h_dedx_p_Coh_l_1","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_1->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_1->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_1->SetContour(50);
  TH2D *h_dedx_p_Coh_r_1 = new TH2D("h_dedx_p_Coh_r_1","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r_1->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_1->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_1->SetContour(50);

  TH2D *h_dedx_p_Coh_l_2 = new TH2D("h_dedx_p_Coh_l_2","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_2->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_2->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_2->SetContour(50);
  TH2D *h_dedx_p_Coh_r_2 = new TH2D("h_dedx_p_Coh_r_2","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r_2->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_2->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_2->SetContour(50);

  TH2D *h_dedx_p_Coh_l_3 = new TH2D("h_dedx_p_Coh_l_3","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_3->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_3->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_3->SetContour(50);
  TH2D *h_dedx_p_Coh_r_3 = new TH2D("h_dedx_p_Coh_r_3","Right side: dE/dx vs. p/q", 250, 40, 300, 250, 0, 1200);
        h_dedx_p_Coh_r_3->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_3->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_3->SetContour(50);

  TH2D *h_dedx_p_Coh_l_4 = new TH2D("h_dedx_p_Coh_l_4","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_4->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_4->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_4->SetContour(50);
  TH2D *h_dedx_p_Coh_r_4 = new TH2D("h_dedx_p_Coh_r_4","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r_4->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_4->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_4->SetContour(50);

  TH2D *h_dedx_p_Coh_l_5 = new TH2D("h_dedx_p_Coh_l_5","Left side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_l_5->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_l_5->SetYTitle("dE/dx");
        h_dedx_p_Coh_l_5->SetContour(50);
  TH2D *h_dedx_p_Coh_r_5 = new TH2D("h_dedx_p_Coh_r_5","Right side: dE/dx vs. p/q", 250, 40, 250, 250, 0, 1200);
        h_dedx_p_Coh_r_5->SetXTitle("p/q [MeV/c]");
        h_dedx_p_Coh_r_5->SetYTitle("dE/dx");
        h_dedx_p_Coh_r_5->SetContour(50);


  TH2D *h_dedx_p_s[7];
    for(int ii=0; ii<7; ii++)
       h_dedx_p_s[ii] = new TH2D(Form("h_dedx_p_s[%d]",ii),Form("dEdx vs. p/q :: Sector %d",ii), 250, 40, 250, 250, 0, 1200);     


  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int initial_number_of_e = 0;
  int k_entry = 0;


//////////////////// loop over the events ///////////////////////////////////////////////


  TCanvas *c1 = new TCanvas("c1");
  TCanvas *c2 = new TCanvas("c2","",1200,800);  c2->Divide(3,2);

  TCanvas *c5 = new TCanvas("c5");
   TLegend* leg = new TLegend(0.65,0.65,0.84,0.84);
            leg-> SetNColumns(1);
            leg->AddEntry(h_rtpc_npd_l,"Left side","L");
            leg->AddEntry(h_rtpc_npd_r,"Right side","L");

  TLegend* dedx_leg = new TLegend(0.75,0.5,0.84,0.84);
           dedx_leg-> SetNColumns(1);
           dedx_leg->AddEntry(f4He," ^{4}He","L");
           dedx_leg->AddEntry(f3He," ^{3}He","L");
           dedx_leg->AddEntry(f3H," ^{3}H","L");
           dedx_leg->AddEntry(fd," d","L");
           dedx_leg->AddEntry(fp," p","L");

  for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
    k_entry = jentry;
    if (jentry>0 && jentry% 1000000 == 0) {
    printf("still running %d \n",(int)jentry);

 c1->cd();

    h_delta_z_e_prot->Draw();
     l->DrawLine(-3*7.89741e-01 +1.19861e-01,0,-3*7.89741e-01 +1.19861e-01,h_delta_z_e_prot->GetMaximum());
     l->DrawLine(3*7.89741e-01 +1.19861e-01,0,3*7.89741e-01 +1.19861e-01, h_delta_z_e_prot->GetMaximum());
     c1->Print("fig/proton_delta_z_e_prot.png");
     c1->Print("fig/proton_delta_z_e_prot.C");

   for(int i = 1; i<7; i++){
      c2->cd(i);
      h_ece_over_p_nphe[i]->Draw("colz");
     }
     c2->Print("fig/etot_over_p_vs_nphe.png");

   for(int i = 1; i<7; i++){
      c2->cd(i);
      h_ch2cc_nphe[i]->Draw("colz");
     }
     c2->Print("fig/ch2cc_vs_nphe.png");

   for(int i = 1; i<7; i++)
   {
    c2->cd(i);
    h_e_nphe_1[i]->Draw();  
    h_e_nphe_2[i]->Draw("same"); 
    h_e_nphe_3[i]->Draw("same");   
//    l->DrawLine(20,0,20,h_e_nphe_1[i]->GetMaximum());
   }
   c2->Print("fig/e_nphe.png");

    
    }
//      if(jentry>3000000) break;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    // to find the helicity
    RunNumber = runnb; 
    event_HWP = Find_HWP(RunNumber);



//    int* Binary_helicity = Decimal_Binary_Converter(helicity_cor);
//    int helicity_cor2 = Binary_helicity[0];
//    event_helicity = Find_helicity(helicity_cor2, event_HWP);

    int myhel = GetHelicity(helicity_cor, RunNumber);
    if(myhel == 0) myhel = -1;

     event_helicity = myhel;
     if(event_helicity ==1) N_pos_e[RunNumber-61500] += 1.0;
       else if(event_helicity == -1) N_neg_e[RunNumber-61500] += 1.0;


//    bool bool_helicity= true; // make sure that the bits 5-24 in the helicity are not zeros
//    if(Binary_helicity[5] == 0 &&  Binary_helicity[6] == 0 &&   Binary_helicity[7] == 0 &&
//       Binary_helicity[8] == 0 &&  Binary_helicity[9] == 0 &&   Binary_helicity[10] == 0 &&
//       Binary_helicity[10] == 0 &&  Binary_helicity[11] == 0 &&   Binary_helicity[12] == 0 &&
//       Binary_helicity[13] == 0 &&  Binary_helicity[14] == 0 &&   Binary_helicity[15] == 0 &&
//       Binary_helicity[16] == 0 &&  Binary_helicity[17] == 0 &&   Binary_helicity[18] == 0 &&
//       Binary_helicity[19] == 0 &&  Binary_helicity[20] == 0 &&   Binary_helicity[21] == 0 && 
//       Binary_helicity[22] == 0 &&  Binary_helicity[23] == 0 &&   Binary_helicity[24] == 0 )
//       bool_helicity= false;          

    // condition on the event from the helicity
  if(myhel == 1 || myhel == -1)
   {
    int id_e_ec[40];      
    int Iam_el= -1;
    int Iam_proton= -1;
    int Iam_photon= -1;
 
    float el_time = 0.; // electron time of flight
    event_id= evntid;

    n_e_ec= 0;
    n_phot_EC= 0;
    n_phot_IC= 0;
    n_photon= 0; // for both EC and IC
    n_pi0= 0; 
    n_proton= 0;
    dc_count= 0;
    n_tpc= 0; 
    float Theta= 0., Phi= 0., beta_cal= 0., beta_meas=0. ;

    // initialisation of Lorentz vectors 
    for(int pp=0; pp<40;pp++)
     {
      Elec4Vector[pp].SetPxPyPzE(0.,0.,0.,0.);
      Pi04Vect[pp].SetPxPyPzE(0.,0.,0.,0.);
      Proton4Vector[pp].SetPxPyPzE(0.,0.,0.,0.);
      Phot4VectEC[pp].SetPxPyPzE(0.,0.,0.,0.);
      Phot4VectIC[pp].SetPxPyPzE(0.,0.,0.,0.);
      Phot4Vect[pp].SetPxPyPzE(0.,0.,0.,0.);
      Phot4VectIC_before[pp].SetPxPyPzE(0.,0.,0.,0.);
      ICPhot3Vect[pp].SetXYZ(0.,0.,0.);
      }

    InEl4Vector.SetPxPyPzE(0, 0, Eb, Eb);
    Target4Vector.SetPxPyPzE(0, 0, 0, M_p);
    GammaStar4Vector.SetPxPyPzE(0., 0., 0., 0.);


 //////////////////////// electron selection in EC////////////////////////////////////////////
                      // gpart bank  for electrons //

   for(int i= 0; i< TMath::Min(25, gpart); i++)
    { 

    bool flag = false;
    if( p[i]>0.7 && p[i]<6.0 && id[i] ==11 && stat[i]>0 && 
        ec_ei[ec[i]-1] >0.0 && ec_eo[ec[i]-1] >0.0  && 
        dc_stat[dc[i]-1]>0 && nphe[cc[i]-1] >0 ){
        // Do the vertez correction of Nathan
        vz[i] = EG6VertexCorrCLAS(RunNumber, vz[i], cx[i], cy[i], cz[i]);
        flag = true;      
        h_z_e->Fill(vz[i]);
       }

  if( flag == true)// && -77.<vz[i] && vz[i]<-50.)
    {
     initial_number_of_e++;
     float rr = sqrt(tl1_x[dc[i]-1]*tl1_x[dc[i]-1] +
                     tl1_y[dc[i]-1]*tl1_y[dc[i]-1] + 
                     pow((tl1_z[dc[i]-1]-vz[i]),2));
     float dc_theta = (TMath::ACos((tl1_z[dc[i]-1]-vz[i])/rr))*TODEG;
     float dc_phi = (TMath::ATan2(tl1_y[dc[i]-1], tl1_x[dc[i]-1])) *TODEG;
     float ec_theta=(TMath::ACos(cz[i]))*TODEG;
     float ec_phi= (TMath::ATan2 (cy[i],cx[i])) *TODEG;
     float ece = TMath::Max(ec_ei[ec[i]-1]+ec_eo[ec[i]-1], etot[ec[i]-1]);
     float ece_p;
     if ( 0.1 <EG6ECsampcorr(dc_sect[dc[i]-1], RunNumber, evntid)) 
         ece_p = (ece/p[i])* (0.31/EG6ECsampcorr(dc_sect[dc[i]-1], RunNumber, evntid)); 
      else ece_p = (ece/p[i]);


     ////////Set the EC_XYZ to EC_UVW  ---------------------------------
     EC_XYZ.SetXYZ(ech_x[ec[i]-1],ech_y[ec[i]-1],ech_z[ec[i]-1]);
     EC_UVW=EC_XYZ_UVW(EC_XYZ);
     float ec_u= EC_UVW(0);
     float ec_v= EC_UVW(1); 
     float ec_w= EC_UVW(2); 
       h_EC_el_U->Fill(ec_u);
       h_EC_el_V->Fill(ec_v);
       h_EC_el_W->Fill(ec_w);

     bool UVW_Fiducial_cut= false;
     if(60.0<ec_u && ec_u<350 && ec_v <370. && ec_w < 400.0) UVW_Fiducial_cut=true;

     ////////// FX backword method to find the shadow of ic on dc --------------
     // track position at L1 DC
     float Xi= tl1_x[dc[i]-1];
     float Yi= tl1_y[dc[i]-1];
     float Zi= tl1_z[dc[i]-1];
     // track direction (cosines) at L1 DC
     float cXi = tl1_cx[dc[i]-1];
     float cYi = tl1_cy[dc[i]-1];
     float cZi = tl1_cz[dc[i]-1];
     // extrapolate track to Z=+16
     TVector3 VI(Xi,Yi,Zi);
     TVector3 PI(cXi,cYi,cZi);
     e_ic_x= 0;
     e_ic_y= 0;
     e_s= 0;
     // get the shift on the track between the IC and the first region of DC         
     if(cZi>0 && Zi!=0) {
       TVector3 Shift= ((Zi-16)/cZi)*PI;
       e_ic_x= (float)(VI-Shift).X();
       e_ic_y= (float)(VI-Shift).Y();
       e_s= dc_sect[dc[i]-1];
       }          

      /////// DC fiducial cut  --------------------------------------
      bool DC_Fiducial_cut= false;
      DC_Fiducial_cut= DCff_e(e_ic_x, e_ic_y , e_s );

       h_e_1_nphe->Fill(nphe[cc[i]-1]);
       h_e_nphe_1[dc_sect[dc[i]-1]]->Fill(nphe[cc[i]-1]);
       h_e_phi_theta_1->Fill(dc_theta, dc_phi);

       ////////// CC coordinates fiducial cuts
       int nt;
       float point[3], dir[3], xx[3], dist, r, s;
       float cc_pln[3] = { - 0.0007840784063, 0., - 0.001681461571 }; //was static in Vlassov's code
       point[0] = dc_xsc[dc[i]-1];
       point[1] = dc_ysc[dc[i]-1];
       point[2] = dc_zsc[dc[i]-1];
       dir[0] = dc_cxsc[dc[i]-1];
       dir[1] = dc_cysc[dc[i]-1];
       dir[2] = dc_czsc[dc[i]-1];
       nt = basics_vcrpl(point, dir, cc_pln, &dist, xx);
          r= sqrt(xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2]);
          s= sqrt(xx[0]*xx[0] + xx[1]*xx[1]);
          float CC_Th_e = TMath::RadToDeg()*TMath::ACos(xx[2]/r);
          float CC_Ph_e = TMath::RadToDeg()*TMath::ATan2(xx[1]/s,xx[0]/s);


       ////// CC fiducial cut --------------------------------
       bool CC_Fiducial_cut= false;
       if(nt){ if(CC_Fid(CC_Th_e,CC_Ph_e)) CC_Fiducial_cut= true;  }//found CC plane intersection

       /////// IC shadow on DC fiducial cut ------------------------
       bool IC_DC_Shadow_Fiducial_cut = false;
       IC_DC_Shadow_Fiducial_cut = IC_DC_shadow_fiducial (e_ic_x, e_ic_y);


       h_EC_XY_el_1->Fill(ech_x[ec[i]-1],ech_y[ec[i]-1]);
       //h_DC_el_1->Fill(e_ic_x, e_ic_y);
       h_CC_el_phi_theta_1->Fill(CC_Th_e, CC_Ph_e); 
       h_el_DC_IC_XY_1->Fill(e_ic_x, e_ic_y);

        if(UVW_Fiducial_cut==true)        h_EC_XY_el_2->Fill(ech_x[ec[i]-1],ech_y[ec[i]-1]);
        if( CC_Fiducial_cut == true)      h_CC_el_phi_theta_2->Fill(CC_Th_e, CC_Ph_e);
        if(IC_DC_Shadow_Fiducial_cut == true)          h_el_DC_IC_XY_2->Fill(e_ic_x, e_ic_y);
        if (DC_Fiducial_cut == true && IC_DC_Shadow_Fiducial_cut == true )  h_DC_el_1->Fill(e_ic_x, e_ic_y);
        if (DC_Fiducial_cut == false || IC_DC_Shadow_Fiducial_cut == false ) h_DC_el_2->Fill(e_ic_x, e_ic_y);

       /// uvw fiducial cut
       if(UVW_Fiducial_cut==true)
        {
         float mean_ece_p = 2.56084e-01 +4.32374e-02*p[i] -9.14180e-03*pow(p[i],2) +8.15895e-04*pow(p[i],3);
         float sigma_p = 0.0572976 -0.0272689*p[i] +0.00857596*pow(p[i],2) -0.000979978*pow(p[i],3);

         h_e_phi_theta_2->Fill(dc_theta, dc_phi);

         //// CC fiducial cut + IC shadow on DC fiducial cut + the DC sectors cut
         if(CC_Fiducial_cut= true && IC_DC_Shadow_Fiducial_cut == true && DC_Fiducial_cut ==true)
          {
           h_e_nphe_2[dc_sect[dc[i]-1]]->Fill(nphe[cc[i]-1]);
           h_e_2_nphe->Fill(nphe[cc[i]-1]);
           h_e_phi_theta_3->Fill(dc_theta, dc_phi);
           h_eco_eci->Fill(ec_ei[ec[i]-1], ec_eo[ec[i]-1]);
           h_etot_p->Fill(p[i], ece_p);
           

           if(abs(ece_p - mean_ece_p)> 2.5 * sigma_p) h_etot_p_2->Fill(p[i], ece_p);

           if(ec_ei[ec[i]-1]>0.06 && abs(ece_p - mean_ece_p)< 2.5 * sigma_p)
            {
             h_e_3_nphe->Fill(nphe[cc[i]-1]);
             h_etot_p_3->Fill(p[i], ece_p);
             h_e_nphe_3[dc_sect[dc[i]-1]]->Fill(nphe[cc[i]-1]);              
             h_ece_over_p_nphe[dc_sect[dc[i]-1]]->Fill(nphe[cc[i]-1], ece_p);              
             h_ch2cc_nphe[dc_sect[dc[i]-1]]->Fill(nphe[cc[i]-1], cc_c2[cc[i]-1]);              

            // if (nphe[cc[i]-1] >20 )
            //  {
               Q2[n_e_ec] = 4*p[i]*Eb*pow(sin(ec_theta*TORAD/2), 2);
               nu[n_e_ec] = (Eb - p[i]);
               yy[n_e_ec] = nu[n_e_ec]/Eb;
               W_4He[n_e_ec]= sqrt(M_4He*M_4He + 2*M_4He*nu[n_e_ec] - Q2[n_e_ec]);
               W_p[n_e_ec]= sqrt(M_p*M_p + 2*M_p*nu[n_e_ec] - Q2[n_e_ec]);
               xB[n_e_ec]= Q2[n_e_ec]/(2*M_p*nu[n_e_ec]);

               h_e_phi_theta_4->Fill(dc_theta, dc_phi);                             
               h_nu->Fill(nu[n_e_ec]);
               h_Q2->Fill(Q2[n_e_ec]);
               h_W_4He->Fill(W_4He[n_e_ec]);
               h_W_p->Fill(W_p[n_e_ec]);
               h_xB->Fill(xB[n_e_ec]);
               h_yy->Fill(yy[n_e_ec]);
               h_Q2_xB->Fill(xB[n_e_ec], Q2[n_e_ec]);
              
               Elec4Vector[n_e_ec].SetPxPyPzE(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i], p[i]);
/////// HERE!!!
	               El_Px[n_e_ec] = Elec4Vector[n_e_ec].Px();
        	       El_Py[n_e_ec] = Elec4Vector[n_e_ec].Py();
	               El_Pz[n_e_ec] = Elec4Vector[n_e_ec].Pz();
        	       El_P[n_e_ec] = Elec4Vector[n_e_ec].P();
	               El_E[n_e_ec] = Elec4Vector[n_e_ec].E();
        	       El_Phi[n_e_ec] = Elec4Vector[n_e_ec].Phi()*TODEG;
	               El_Theta[n_e_ec] = Elec4Vector[n_e_ec].Theta()*TODEG;
        	       El_z[n_e_ec] =vz[i];
	               El_s[n_e_ec] = dc_sect[dc[i]-1];
	               El_nphe[n_e_ec] = nphe[cc[i]-1];
               
               Iam_el= i;
               n_e_ec++;
               id_e_ec[n_e_ec - 1] = i;

               N_e_Run[runnb-61500] += 1.0;

             //  }
              }
             }
           }
          }
        } /////// end gpart bank for electrons from EC ------------------------------------------------------



///////////////////////////////////////////////////////////////////////////////////////////////////////
         //////////// gpart bank  for charged particels and photons from EC //////////// 

 if(n_e_ec ==1)
  {
   GammaStar4Vector= InEl4Vector - Elec4Vector[0];

   for(int i= 0; i< TMath::Min(25, gpart); i++)
    { 
     float ece = TMath::Max(ec_ei[ec[i]-1]+ec_eo[ec[i]-1], etot[ec[i]-1]);
  
     //Set the EC_XYZ to EC_UVW
     EC_XYZ.SetXYZ(ech_x[ec[i]-1],ech_y[ec[i]-1],ech_z[ec[i]-1]);
     EC_UVW=EC_XYZ_UVW(EC_XYZ);
     float ec_u= EC_UVW(0);
     float ec_v= EC_UVW(1); 
     float ec_w= EC_UVW(2); 

     ////////// FX backword method to find the shadow of ic on dc --------------
     // track position at L1 DC
     float Xi= tl1_x[dc[i]-1];
     float Yi= tl1_y[dc[i]-1];
     float Zi= tl1_z[dc[i]-1];
     // track direction (cosines) at L1 DC
     float cXi = tl1_cx[dc[i]-1];
     float cYi = tl1_cy[dc[i]-1];
     float cZi = tl1_cz[dc[i]-1];
     // extrapolate track to Z=+16
     TVector3 VI(Xi,Yi,Zi);
     TVector3 PI(cXi,cYi,cZi);
     e_ic_x= 0;
     e_ic_y= 0;
     e_s= 0;
     // get the shift on the track between the IC and the first region of DC         
     if(cZi>0 && Zi!=0)
      {
       TVector3 Shift= ((Zi-16)/cZi)*PI;
       e_ic_x= (float)(VI-Shift).X();
       e_ic_y= (float)(VI-Shift).Y();
       e_s= dc_sect[dc[i]-1];
       }   

 //////------------------------------- proton selection --------------------------------///////

    // ...Define Electron time to evaluate beta of the other particles ...
    el_time = sc_t[sc[id_e_ec[0]] - 1] - sc_r[sc[id_e_ec[0]] - 1]/LightSpeed;

    if( b[i]!= 0. && q[i]>0 && sc_stat[sc[i]-1]>0 && dc_stat[dc[i]-1]>0 && DC_Mid_Cut(e_ic_x, e_ic_y) )
     {
      
      if(TMath::Abs(beta_meas - beta_cal)<0.025 ) h_proton_z->Fill(vz[i]);
       float rr = sqrt(tl1_x[dc[i]-1]*tl1_x[dc[i]-1] + 
                       tl1_y[dc[i]-1]*tl1_y[dc[i]-1] + 
                       pow((tl1_z[dc[i]-1]-vz[id_e_ec[0]]),2));
       Theta = (TMath::ACos((tl1_z[dc[i]-1]-vz[id_e_ec[0]])/rr))*TODEG;
       Phi =   (TMath::ATan2 (tl1_y[dc[i]-1], tl1_x[dc[i]-1])) *TODEG;
      //Theta = (TMath::ACos(cz[i]))*TODEG;
      //Phi = (TMath::ATan2 (cy[i],cx[i])) *TODEG;
     if(-77.<vz[i] && vz[i]<-50. )
      {
	      // Define Delta Beta for the proton
	      beta_cal=p[i]/sqrt(p[i]*p[i] + M_p*M_p);        // Bcalc= p/E
	      beta_meas = b[i];                               // Bmeasured= l/c.t
         
              h_prot_Theta_Phi_1->Fill(Theta, Phi);          
	      h_proton_DC_IC_XY_1->Fill(e_ic_x, e_ic_y);
	      // ic shadow on DC fiducial cut
	      if( !IC_DC_shadow_fiducial (e_ic_x, e_ic_y) || !
			DCff_e(e_ic_x, e_ic_y , e_s ) ) 
                  h_proton_DC_IC_XY_2->Fill(e_ic_x, e_ic_y);
           
	      //// IC shadow on DC fiducial cut + the DC sectors cut
	      if(IC_DC_shadow_fiducial (e_ic_x, e_ic_y) && DCff_e(e_ic_x, e_ic_y , e_s ) )
	       {
		     h_proton_DC_IC_XY_3->Fill(e_ic_x, e_ic_y);
	   	     DC_IC_stat[dc_count]= +11.0; // protons which pass the DC+IC fiducial cuts
	  	     dc_count++;

	    	h_prot_DeltaT_p->Fill(p[i], (sc_t[i] -el_time) - sc_r[i]*sqrt(p[i]*p[i] + M_p*M_p)/(p[i]*LightSpeed));
		h_prot_DeltaBeta_p->Fill(p[i],beta_meas - beta_cal);
	        h_prot_DeltaBeta->Fill(beta_meas - beta_cal);

		h_prot_Beta_p->Fill(p[i],b[i]);
                 if(TMath::Abs(beta_meas - beta_cal)<0.025 ) h_delta_z_e_prot->Fill(vz[id_e_ec[0]] - vz[i]); 
	        // delta beta cut
	        if(TMath::Abs(beta_meas - beta_cal)<0.025  &&  id[i]==2212 && abs(vz[id_e_ec[0]] - vz[i])< 2.0 )
	         {
		          h_proton_DC_IC_XY_4->Fill(e_ic_x, e_ic_y);
	        	  h_prot_Theta_Phi_2->Fill(Theta, Phi);

		          Proton4Vector[n_proton].SetPxPyPzE(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i],
                                                             sqrt(p[i]*p[i]+M_p*M_p));


	        	  Proton_Px[n_proton] = Proton4Vector[n_proton].Px();
		          Proton_Py[n_proton] = Proton4Vector[n_proton].Py();
        		  Proton_Pz[n_proton] = Proton4Vector[n_proton].Pz();
		          Proton_P[n_proton] = Proton4Vector[n_proton].P();
        		  Proton_E[n_proton] = Proton4Vector[n_proton].E();
		          Proton_Theta[n_proton] = Proton4Vector[n_proton].Theta()*TODEG;
        		  Proton_Phi[n_proton] = Proton4Vector[n_proton].Phi()*TODEG;
		          Proton_M[n_proton] = Proton4Vector[n_proton].M();
       		          Proton_z[n_proton] = vz[i];         


		          Iam_proton= i;
		          n_proton++;
    	  	          N_p_Run[runnb-61500] += 1.0;

                   }
                 }

            else
              {
                DC_IC_stat[dc_count]= -11.0; // protonss not passing DC+IC fid cuts
               dc_count++;
              }
          }
        } // end of proton selection


 /////-------------------------- photon selection in EC ----------------------------------------//////
      
    if( q[i]== 0 && b[i]!= 0) 
     {
     Theta = TMath::ATan2(TMath::Sqrt( TMath::Power(cx[i],2) + TMath::Power(cy[i],2)), cz[i])*TODEG; 
     Phi = TMath::ATan2(cy[i],cx[i])*TODEG;

      h_photon_phi_theta_1->Fill(Theta, Phi);
      bool photons_uvw_cut = false;
      if(100.<ec_u && ec_u<390. && ec_v <360. && ec_w < 390.)
        {
          photons_uvw_cut = true;
          h_EC_XY_photon_2->Fill(ech_x[ec[i]-1],ech_y[ec[i]-1]);
         }
       else h_EC_XY_photon_1->Fill(ech_x[ec[i]-1],ech_y[ec[i]-1]);

      if( ece/0.273>0.3 && 100.<ec_u && ec_u<390. && ec_v <360. && ec_w < 390.) h_photon_beta->Fill(b[i]);
 
      h_EC_photon_U->Fill(ec_u);
      h_EC_photon_V->Fill(ec_v);
      h_EC_photon_W->Fill(ec_w);

      /// uvw fiducial cut 
      if( ece/0.273>0.3 && 0.93<b[i] && b[i]<1.07 && 100.<ec_u && ec_u<390. && ec_v <360. && ec_w < 390.)
       {
        
	        Phot4VectEC[n_phot_EC].SetPxPyPzE(ece*cx[i]/0.273,ece*cy[i]/0.273,ece*cz[i]/0.273,ece/0.273);
        	Phot4Vect[n_photon].SetPxPyPzE(ece*cx[i]/0.273,ece*cy[i]/0.273,ece*cz[i]/0.273,ece/0.273);
	        Ph_EC_Px[n_phot_EC]=Phot4VectEC[n_phot_EC].Px();
        	Ph_EC_Py[n_phot_EC]=Phot4VectEC[n_phot_EC].Py();
		Ph_EC_Pz[n_phot_EC]=Phot4VectEC[n_phot_EC].Pz();
        	Ph_EC_P[n_phot_EC]=Phot4VectEC[n_phot_EC].P();
	        Ph_EC_E[n_phot_EC]=Phot4VectEC[n_phot_EC].E();			       
        	Ph_EC_Theta[n_phot_EC]=Phot4VectEC[n_phot_EC].Theta()*TODEG;
	        Ph_EC_Phi[n_phot_EC]=Phot4VectEC[n_phot_EC].Phi()*TODEG;

        	Phot4Vect[n_photon].SetPxPyPzE(ece*cx[i]/0.273,ece*cy[i]/0.273,ece*cz[i]/0.273,ece/0.273);
     
        	Ph_Px[n_photon]=Phot4Vect[n_photon].Px();
	        Ph_Py[n_photon]=Phot4Vect[n_photon].Py();
        	Ph_Pz[n_photon]=Phot4Vect[n_photon].Pz();
	        Ph_P[n_photon]=Phot4Vect[n_photon].P();
        	Ph_E[n_photon]=Phot4Vect[n_photon].E();
	        Ph_Theta[n_photon]=Phot4Vect[n_photon].Theta()*TODEG;
        	Ph_Phi[n_photon]=Phot4Vect[n_photon].Phi()*TODEG;

       	        h_photon_phi_theta_2->Fill(Theta, Phi);
       
	        n_phot_EC++;
        	n_photon++;
	        Iam_photon= i;
        }
       } // end of photon selction in EC



     }// end gpart bank for other particles


////////////////////////////////// photon selection in IC///////////////////////////////////////////
                         ////////     ic_part bank      /////////

  
   // use ic_part bank for the fiducial cut and icpart for energy an. angels 
   int N=0; // relative index between icpart to ic_part bank 
   for(int j=0; j<TMath::Min(icpart,25); j++) 
     {
      N= (statc[j]-statc[j]%10000)/10000-1;
      //N=j;
      ICPhot3Vect[j].SetXYZ(ich_x[N], ich_y[N], -vz[id_e_ec[0]]);
      Phot4VectIC_before[j].SetPxPyPzE(etc[j]*(ICPhot3Vect[j].Unit()).X(), 
                                       etc[j]*(ICPhot3Vect[j].Unit()).Y(),
                                       etc[j]*(ICPhot3Vect[j].Unit()).Z(),
                                       etc[j]);

         h_photon_phi_theta_1->Fill(Phot4VectIC_before[j].Theta()*TODEG, Phot4VectIC_before[j].Phi()*TODEG);
         if( !ICFiducialCut(xc[j], yc[j])) h_IC_XY_1->Fill(xc[j], yc[j]);
          else if( ICFiducialCut(xc[j], yc[j]) && IC_Not_Hot_Channel(ich_xgl[N], ich_ygl[N]) == true) 
              {
                h_IC_XY_2->Fill(xc[j], yc[j]);
                h_IC_theta_En_1->Fill(etc[j], Phot4VectIC_before[j].Theta()*TODEG);
                }                                                                       
 
      h_eng_photon->Fill(etc[j]);
     // take out hot channels
     if(IC_Not_Hot_Channel(ich_xgl[N], ich_ygl[N]) == true)
       {    
        h_ic_time->Fill(time[N]);
        h_ic_theta->Fill(Phot4VectIC_before[j].Theta()*TODEG);

       if(IC_theta_En_Cut(etc[j], Phot4VectIC_before[j].Theta()*TODEG) && 
          ICFiducialCut(xc[j], yc[j])  && 
          (Phot4VectIC_before[j].Theta()*TODEG)<14.0 && 
           etc[j]>0.3 ) 
        {
     		    h_IC_theta_En_2->Fill(etc[j], Phot4VectIC_before[j].Theta()*TODEG);
                    h_IC_XY_center->Fill(ich_xgl[N],ich_ygl[N]);
 
                    float E_gamma = etc[j];

                    Phot4VectIC[n_phot_IC].SetPxPyPzE(E_gamma*(ICPhot3Vect[j].Unit()).X(), 
                                                      E_gamma*(ICPhot3Vect[j].Unit()).Y(),
                                                      E_gamma*(ICPhot3Vect[j].Unit()).Z(), 
                                                      E_gamma);
		         Ph_IC_Px[n_phot_IC]=Phot4VectIC[n_phot_IC].Px();
        		 Ph_IC_Py[n_phot_IC]=Phot4VectIC[n_phot_IC].Py();
	        	 Ph_IC_Pz[n_phot_IC]=Phot4VectIC[n_phot_IC].Pz();
	        	 Ph_IC_P[n_phot_IC]=Phot4VectIC[n_phot_IC].P();
		         Ph_IC_E[n_phot_IC]=Phot4VectIC[n_phot_IC].E();
        		 Ph_IC_Theta[n_phot_IC]=Phot4VectIC[n_phot_IC].Theta()*TODEG;
		         Ph_IC_Phi[n_phot_IC]=Phot4VectIC[n_phot_IC].Phi()*TODEG;

                     h_photon_phi_theta TLorentzVector( _2->Fill(Ph_IC_Theta[n_phot_IC], Ph_IC_Phi[n_phot_IC]);
                     Phot4Vect[n_photon].SetPxPyPzE(E_gamma*(ICPhot3Vect[j].Unit()).X(), 
                                                    E_gamma*(ICPhot3Vect[j].Unit()).Y(),
                                                    E_gamma*(ICPhot3Vect[j].Unit()).Z(),
                                                    E_gamma);

        		 Ph_Px[n_photon]=Phot4Vect[n_photon].Px();
	        	 Ph_Py[n_photon]=Phot4Vect[n_photon].Py();
	        	 Ph_Pz[n_photon]=Phot4Vect[n_photon].Pz();
		         Ph_P[n_photon]=Phot4Vect[n_photon].P();
        		 Ph_E[n_photon]=Phot4Vect[n_photon].E();
		         Ph_Theta[n_photon]=Phot4Vect[n_photon].Theta()*TODEG;
		         Ph_Phi[n_photon]=Phot4Vect[n_photon].Phi()*TODEG;

                      n_photon++; 
                      n_phot_IC++;
          }
        }
      } ////////endl ic_part bank ---------------------------------------------------       

//  if(n_phot_IC <0) cout<<n_phot_IC<<endl;

/////////////////////////////////////////////////////////////////////////////////////////////////////
///                              helium selection in RTPC                                       /////
//                                                                                             ///// 
/////////////////////////////////////////////////////////////////////////////////////////////////////

  int test_rtpc = 0;
  for(int i=0; i<TMath::Min(30, gcpart); i++)
   {
    float Delta_z= ((vz[id_e_ec[0]]+64)*10 - z[i] );  //! -FC 
    if (abs(Delta_z)<50 && r_0[i]> 0 ) h_rtpc_phi_theta_1->Fill(theta[i]*TODEG, phi[i]*TODEG);

  if ( AcceptRTPC(vz[id_e_ec[0]], theta[i], phi[i]))
  {
   if (abs(Delta_z)<50)
    {
  	  if(90.0<phi[i]*TODEG && phi[i]*TODEG <270.0)  h_rtpc_r0_r->Fill(r_0[i]);
	       else h_rtpc_r0_l->Fill(r_0[i]);

	    if( 90.0<phi[i]*TODEG && phi[i]*TODEG <270.0)
	     {
	        if(r_0[i]> 0)   h_rtpc_npd_r->Fill(npd_track[i]);
	        if(r_0[i]> 0 && npd_track[i]>3)
	         {
                   h_rtpc_z_r->Fill(z[i]);
                   h_rtpc_X2_r->Fill(x2[i]);
                   if( abs(z[i]) <100.0 && x2[i] <3.5 )
                     {
	          	 h_rtpc_sdist_r->Fill(sdist[i]);
	          	 h_rtpc_edist_r->Fill(edist[i]);
 	              	 if (abs(sdist[i]) < 2.0 && edist[i] < 5.0) h_delta_z_r->Fill(Delta_z);
                     }
	          }
	        }
	     else if( 90.0>phi[i]*TODEG || phi[i]*TODEG >270.0)
	      {
	        if(r_0[i]> 0) h_rtpc_npd_l->Fill(npd_track[i]);
	        if(r_0[i]> 0 && npd_track[i]>3)
	         {
	           h_rtpc_X2_l->Fill(x2[i]);
	           h_rtpc_z_l->Fill(z[i]);
                   if( abs(z[i]) <100.0 && x2[i] <3.5 )
                    {
		        h_rtpc_sdist_l->Fill(sdist[i]);
		        h_rtpc_edist_l->Fill(edist[i]);
		        if (abs(sdist[i])< 2.0 && edist[i] < 5.0) h_delta_z_l->Fill(Delta_z); 
                     }
	          }      
	       }
        } // end if (abs(Delta_z)<50)

      if(p_tot[i]<1000.0  && -100.0<z[i] && z[i]<100.0  && x2[i]<3.5  &&
         r_0[i]>0. && npd_track[i]>3  && abs(sdist[i]) < 2.0 &&
         -1.0<edist[i] && edist[i]<5.0 && abs(Delta_z)<20.0)
         {
          test_rtpc++;
          h_rtpc_phi_theta_2->Fill(theta[i]*TODEG, phi[i]*TODEG);
         }
       }

    }
 
   if(test_rtpc ==1)
    {
       // rtpc origional back
       for(int i=0; i<TMath::Min(30, gcpart); i++)
        {
          float Delta_z= (z[i] - (vz[id_e_ec[0]]+64)*10); // in mm
          if( AcceptRTPC(vz[id_e_ec[0]], theta[i], phi[i]) &&
             p_tot[i]<1000.0  && -100.0<z[i] && z[i]<100.0  && x2[i]<3.5  &&
             r_0[i]>0. && npd_track[i]>3  && abs(sdist[i]) < 2.0 &&
             -1.0<edist[i] && edist[i]<5.0 && abs(Delta_z)<20.0)
             {
               if(90.0<phi[i]*TODEG && phi[i]*TODEG <270.0) h_dedx_p_Coh_r->Fill(p_tot[i], dedx[i]);
                else if(90.0>phi[i]*TODEG || phi[i]*TODEG >270.0)  h_dedx_p_Coh_l->Fill(p_tot[i], dedx[i]);
              }
          }
     
        // rtpc corrcted bank
       for(int i = 0; i < TMath::Min(30, rtpc_npart); i++)
        {
          float Delta_z= (z[rtpc_gcpb[i]] - (vz[id_e_ec[0]]+64)*10); // in mm
          if(AcceptRTPC(vz[id_e_ec[0]], theta[rtpc_gcpb[i]], phi[rtpc_gcpb[i]])  &&
             p_tot[rtpc_gcpb[i]]<1000.0  &&  abs(z[rtpc_gcpb[i]])<100.0  &&
             x2[rtpc_gcpb[i]]<3.5   &&  r_0[rtpc_gcpb[i]]>0.0  &&
             npd_track[rtpc_gcpb[i]]>3 && abs(sdist[rtpc_gcpb[i]])<2.0 &&
             -1.0<edist[rtpc_gcpb[i]] && edist[rtpc_gcpb[i]]<5.0 && abs(Delta_z)<20.0)
            {
              
              h_dedx_p_s[El_s[0]]->Fill(rtpc_poverq[i], rtpc_dedx2[i]);

              if(90.0<rtpc_phi[i]*TODEG && rtpc_phi[i]*TODEG <270.0)
               {
                  h_dedx_p_Coh_r_0->Fill(rtpc_poverq[i], rtpc_dedx[i]);
                  h_dedx_p_Coh_r_1->Fill(rtpc_poverq[i], rtpc_dedx2[i]);
                  h_dedx_p_Coh_r_2->Fill(rtpc_poverq[i], rtpc_dedxa[i]);
                  h_dedx_p_Coh_r_3->Fill(rtpc_poverq[i], rtpc_dedxl[i]);
                  h_dedx_p_Coh_r_4->Fill(rtpc_poverq[i], rtpc_dedxal[i]);
                  h_dedx_p_Coh_r_5->Fill(rtpc_poverq[i], rtpc_dedxs[i]);
                }
              else if( 90.0>rtpc_phi[i]*TODEG || rtpc_phi[i]*TODEG >270.0)
                {
                  h_dedx_p_Coh_l_0->Fill(rtpc_poverq[i], rtpc_dedx[i]);
                  h_dedx_p_Coh_l_1->Fill(rtpc_poverq[i], rtpc_dedx2[i]);
                  h_dedx_p_Coh_l_2->Fill(rtpc_poverq[i], rtpc_dedxa[i]);
                  h_dedx_p_Coh_l_3->Fill(rtpc_poverq[i], rtpc_dedxl[i]);
                  h_dedx_p_Coh_l_4->Fill(rtpc_poverq[i], rtpc_dedxal[i]);
                  h_dedx_p_Coh_l_5->Fill(rtpc_poverq[i], rtpc_dedxs[i]);
                 }                                                          
            }
        } 
      } // end  if(test_rtpc ==1)


  n_tpc = 0;
  for(int i = 0; i < TMath::Min(30, rtpc_npart); i++)
   {
         // variables from gcpb bank where we need the index conversion :rtpc_gcpb[i]
	 pid_tpc[n_tpc]= pid[rtpc_gcpb[i]];
         x_tpc[n_tpc]= x[rtpc_gcpb[i]];
	 y_tpc[n_tpc]= y[rtpc_gcpb[i]];
	 z_tpc[n_tpc]= z[rtpc_gcpb[i]];
         dedx_tpc[n_tpc]= dedx[rtpc_gcpb[i]];
	 px_tpc[n_tpc]= px[rtpc_gcpb[i]];
         py_tpc[n_tpc]= py[rtpc_gcpb[i]];
	 pz_tpc[n_tpc]= pz[rtpc_gcpb[i]];
         ptot_tpc[n_tpc]= p_tot[rtpc_gcpb[i]];
	 x2_tpc[n_tpc]= x2[rtpc_gcpb[i]];
         theta_tpc[n_tpc]= theta[rtpc_gcpb[i]];
	 charge_tpc[n_tpc]= charge[rtpc_gcpb[i]];
         dca_tpc[n_tpc]= dca[rtpc_gcpb[i]];
	 index_tpc[n_tpc]= index[rtpc_gcpb[i]];
         phi_tpc[n_tpc]= phi[rtpc_gcpb[i]];
	 vtl_tpc[n_tpc]= vtl[rtpc_gcpb[i]];
         sdist_tpc[n_tpc]= sdist[rtpc_gcpb[i]];
	 edist_tpc[n_tpc]= edist[rtpc_gcpb[i]];
         npts_tpc[n_tpc]= npts[rtpc_gcpb[i]];
	 r0_tpc[n_tpc]= r_0[rtpc_gcpb[i]];
         fiterr_tpc[n_tpc]= fiterr[rtpc_gcpb[i]];
	 tothits_tpc[n_tpc]= tothits[rtpc_gcpb[i]];
         npd_track_tpc[n_tpc]= npd_track[rtpc_gcpb[i]];
	 npd_event_tpc[n_tpc]= npd_event[rtpc_gcpb[i]];
         bonus_bits_tpc[n_tpc]= bonus_bits[rtpc_gcpb[i]];
	 qtot_tpc[n_tpc]= q_tot[rtpc_gcpb[i]];
         x_start_tpc[n_tpc]= x_start[rtpc_gcpb[i]];
	 y_start_tpc[n_tpc]= y_start[rtpc_gcpb[i]];
         z_start_tpc[n_tpc]= z_start[rtpc_gcpb[i]];
	 x_end_tpc[n_tpc]= x_end[rtpc_gcpb[i]];
         y_end_tpc[n_tpc]= y_end[rtpc_gcpb[i]];
	 z_end_tpc[n_tpc]= z_end[rtpc_gcpb[i]];

        // varibles from the corrected rtpc bank
	 rtpc_id1_c[n_tpc]= rtpc_id1[i];
         rtpc_id2_c[n_tpc]= rtpc_id2[i];
	 rtpc_id3_c[n_tpc]= rtpc_id3[i];
         rtpc_id4_c[n_tpc]= rtpc_id4[i];
	 rtpc_id5_c[n_tpc]= rtpc_id5[i];
         rtpc_p1_c[n_tpc]= rtpc_p1[i];
	 rtpc_p2_c[n_tpc]= rtpc_p2[i];
         rtpc_p3_c[n_tpc]= rtpc_p3[i];
	 rtpc_p4_c[n_tpc]= rtpc_p4[i];
         rtpc_p5_c[n_tpc]= rtpc_p5[i];
	 rtpc_poverq_c[n_tpc]= rtpc_poverq[i];
         rtpc_dedx_c[n_tpc]= rtpc_dedx[i];
	 rtpc_dedx2_c[n_tpc]= rtpc_dedx2[i];
         rtpc_theta_c[n_tpc]= rtpc_theta[i];
	 rtpc_phi_c[n_tpc]= rtpc_phi[i];
         rtpc_vz_c[n_tpc]= rtpc_vz[i];
	 rtpc_bad_c[n_tpc]= rtpc_bad[i];
         rtpc_gcpb_c[n_tpc]= rtpc_gcpb[i];

         n_tpc++;
     } // end of:::> for(int i = 0; i < gcpart; i++)
 


/////////////////////////////////////////////////////////////////////////////////////////////////
///////////                  pi_zero selection                              /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////  
 
  if(n_photon >=2)
   {
    int i=0;
    float Mass_1= 0., Mass_2= 0., Mass_3= 0.;
    for(int kk=0; kk<n_photon; kk++){
      for (int mm=0; mm<n_photon; mm++){
        Pi04Vect[i].SetPxPyPzE(0.,0.,0.,0.);
       // both from EC
        if(n_phot_EC>=2 && mm<n_phot_EC && kk<n_phot_EC && mm>kk){
          Pi04Vect[i]= Phot4Vect[mm] + Phot4Vect[kk];
          Mass_1= CalcTwoPhotonInvarMass(Phot4Vect[mm], Phot4Vect[kk]);
          if(0.01<Mass_1 && Mass_1<2.){
            Pi0_InvMass2photECEC= Mass_1;
            h_pi0_ECEC->Fill(Mass_1);
            i++;
            }
          }
        // both from IC
        if( n_phot_IC>=2 && mm>=n_phot_EC && kk>=n_phot_EC &&  mm>kk){
          Pi04Vect[i]= Phot4Vect[mm] + Phot4Vect[kk];
          Mass_2= CalcTwoPhotonInvarMass(Phot4Vect[mm], Phot4Vect[kk]);
          if(0.01<Mass_2 && Mass_2<2.){
            Pi0_InvMass2photICIC= Mass_2;
            h_pi0_ICIC->Fill(Mass_2);
            i++;
            }
          }
        // one from IC and one from EC -ordered
        if(n_phot_IC !=0 && n_phot_EC !=0 && mm>n_phot_EC && kk<n_phot_EC && mm>kk){
          Pi04Vect[i]= Phot4Vect[mm] + Phot4Vect[kk];
          Mass_3= CalcTwoPhotonInvarMass(Phot4Vect[mm], Phot4Vect[kk]);
          if(0.01<Mass_3 && Mass_3<2.){
            Pi0_InvMass2photECIC= Mass_3;
            h_pi0_ECIC->Fill(Mass_3);
            i++;
            }
          }
 
        }
       } 
      
     n_pi0 = i;
     if(n_pi0 <20) {
       for(int j=0; j<n_pi0; j++) {
         Pi0_Px[j]= Pi04Vect[j].Px();
         Pi0_Py[j]= Pi04Vect[j].Py();
         Pi0_Pz[j]= Pi04Vect[j].Pz();
         Pi0_E[j]=  Pi04Vect[j].E();
         Pi0_Phi[j]= Pi04Vect[j].Phi()*TODEG;
         Pi0_Theta[j]= Pi04Vect[j].Theta()*TODEG;
         Pi0_P[j]=  Pi04Vect[j].P();
         Pi0_M[j]=  Pi04Vect[j].M();
         }
       }      
    } // end if(nphot >=2)

   tree->Fill();
   }  // end of  if(n_e_ec ==1) 
  }  // end of   if(helicity_cor>0 && 5_24_helicity == true)
  } // end of for (Long64_t jentry=0; jentry<nentries;jentry++)
 
//  tree->Write();
  tree->AutoSave();
  skim_file->Write();

cout<<"  num_entries = "<<k_entry<<endl;
cout<<"  initial_number_of_e = "<<initial_number_of_e<<endl;
l->SetLineColor(2);


//////////////////////// helicity vs. run for all the electrons //////////////////////
 TCanvas *c44 = new TCanvas("c44", "", 1300, 550); c44->cd();
 for(int i =0; i<n_run; i++)
    {
      int RUN = (i + 61500);
      float AA = 0.0;
      if(N_pos_e[i]> 0. || N_neg_e[i]> 0.)
        {
         AA = (N_pos_e[i] - N_neg_e[i]) / (N_pos_e[i] + N_neg_e[i]);
         h_A_RunNum->Fill(RUN, AA);
         }
      cout<<RUN<<"      "<<N_pos_e[i]<<"      "<<N_neg_e[i]<<"      "<<AA<<endl;
    }
   h_A_RunNum->Draw(); c44->Print("ALU_RunNum.png");    c44->Print("ALU_RunNum.C");

///////////////////////////////////////////////////////////////////////////////////////////////////
// find the number of the p/e vs. Run number ///////

  TCanvas *c4 = new TCanvas("c4", "", 1000, 550);
  c4->cd();

  float Run_Number[n_run];
  float Run_Number_err[n_run];
  float p_over_e[n_run];
  float p_over_e_err[n_run];

  for(int i =0; i<n_run; i++)
     {
      Run_Number[i] = (float) (i + 61500);
      Run_Number_err[i] = 0.;
      p_over_e[i] = 0.;
      p_over_e_err[i] = 0.;
      if( N_p_Run[i]> 0. && N_e_Run[i] > 0.)
       {
         p_over_e[i] =  N_p_Run[i]/N_e_Run[i];
         p_over_e_err[i] = Find_Error(N_p_Run[i], N_e_Run[i]);
        }
cout<<Run_Number[i]<<"      "<<p_over_e[i]<<"     "<<p_over_e_err[i]<<endl;
      }

    TGraphErrors *plot_p_over_e_Run = new TGraphErrors(n_run, Run_Number, p_over_e, Run_Number_err, p_over_e_err);
    plot_p_over_e_Run->SetTitle("<p/e> vs. RunNumber");
    plot_p_over_e_Run->GetYaxis()->SetTitle("<p/e>");
    plot_p_over_e_Run->GetXaxis()->SetTitle("RunNumber");
    plot_p_over_e_Run->SetMarkerStyle(21);
    plot_p_over_e_Run->Draw("AP");
    c4->Print("prot_over_e_Run.png");
    c4->Print("prot_over_e_Run.C");


  c1->cd();

  /// electron -------------------------------------------------------------------
  h_z_e->Draw();
    l->DrawLine(-77,0,-77,h_z_e->GetMaximum());
    l->DrawLine(-50,0,-50,h_z_e->GetMaximum());
    c1->Print("fig/Z_e.png");
  h_EC_el_U->Draw();
    l->DrawLine(350,0,350,h_EC_el_U->GetMaximum());
    l->DrawLine(60,0,60,h_EC_el_U->GetMaximum());
    c1->Print("fig/EC_el_U.png");
  h_EC_el_V->Draw();
    l->DrawLine(370,0,370,h_EC_el_V->GetMaximum());
    c1->Print("fig/EC_el_V.png");
  h_EC_el_W->Draw();
    l->DrawLine(400,0,400,h_EC_el_W->GetMaximum());
    c1->Print("fig/EC_el_W.png");
  h_e_1_nphe->Draw();
     h_e_2_nphe->Draw("same");  h_e_3_nphe->Draw("same");
     l->DrawLine(20,0,20,h_e_1_nphe->GetMaximum());
     c1->Print("fig/nphe_all.png");

  h_EC_XY_el_1->Draw("colz");   c1->Print("fig/EC_XY_el_1.png");
  h_EC_XY_el_2->Draw("colz");   c1->Print("fig/EC_XY_el_2.png");
  h_DC_el_1->Draw("colz");   c1->Print("fig/DC_el_1.png");
  h_DC_el_2->Draw("colz");   c1->Print("fig/DC_el_2.png");
  h_e_phi_theta_1->Draw("colz");  c1->Print("fig/e_phi_theta_1.png");
  h_e_phi_theta_2->Draw("colz");  c1->Print("fig/e_phi_theta_2.png");
  h_e_phi_theta_3->Draw("colz");  c1->Print("fig/e_phi_theta_3.png");
  h_e_phi_theta_4->Draw("colz");  c1->Print("fig/e_phi_theta_4.png");
  h_el_DC_IC_XY_1->Draw("colz"); c1->Print("fig/el_DC_IC_XY_1.png");
  h_el_DC_IC_XY_2->Draw("colz"); c1->Print("fig/el_DC_IC_XY_2.png");
  h_CC_el_phi_theta_1->Draw("colz"); c1->Print("fig/CC_el_phi_theta_1.png");
  h_CC_el_phi_theta_2->Draw("colz"); c1->Print("fig/CC_el_phi_theta_2.png");
  h_etot_p->Draw("colz");   c1->Print("fig/etot_p.png");
  h_etot_p_3->Draw("colz");   c1->Print("fig/etot_p_3.png");
  h_etot_p_2->Draw("colz");   c1->Print("fig/etot_p_2.png");
  h_eco_eci->Draw("colz");  l->DrawLine(0.06,0,0.06,0.45); c1->Print("fig/eco_eci.png");



  h_W_4He->Draw();  c1->Print("fig/W_4He.png");
  h_W_p->Draw();   c1->Print("fig/W_p.png");
  h_Q2->Draw();   c1->Print("fig/Q2.png");
  h_xB->Draw();   c1->Print("fig/xB.png");
  h_yy->Draw();   c1->Print("fig/yy.png");
  h_nu->Draw();   c1->Print("fig/nu.png");
  h_Q2_xB->Draw("colz");   c1->Print("fig/Q2_xB.png");


  for(int i = 1; i<7; i++)
   {
    c2->cd(i);
    h_e_nphe_1[i]->Draw();  
    h_e_nphe_2[i]->Draw("same"); 
    h_e_nphe_3[i]->Draw("same");   
    l->DrawLine(20,0,20,h_e_nphe_1[i]->GetMaximum());
   }
   c2->Print("fig/e_nphe.png");


  for(int i = 1; i<7; i++)
   {
     c2->cd(i);
     h_ece_over_p_nphe[i]->Draw("colz");
     }
   c2->Print("fig/etot_over_p_vs_nphe.png");




  c1->cd();


        

  /// pi(0)-----------------------------------------------------------------------
  h_pi0_ECEC->Draw(); c1->Print("fig/pi0_ECEC.png");
  h_pi0_ICIC->Draw(); c1->Print("fig/pi0_ICIC.png");
  h_pi0_ECIC->Draw(); c1->Print("fig/pi0_ECIC.png");

  /// proton --------------------------------------------------------------------
  h_prot_DeltaBeta_p->Draw("colz"); 
    l->DrawLine(0,3*7.89741e-01 +1.19861e-01, 2.95, 3*7.89741e-01 +1.19861e-01);
    l->DrawLine(0,-3*7.89741e-01 +1.19861e-01, 2.95, -3*7.89741e-01 +1.19861e-01);
    c1->Print("fig/prot_DeltaBeta_p.png");
   h_delta_z_e_prot->Draw();
    l->DrawLine(-3*7.89741e-01 +1.19861e-01,0,-3*7.89741e-01 +1.19861e-01,h_delta_z_e_prot->GetMaximum());
    l->DrawLine(3*7.89741e-01 +1.19861e-01,0,3*7.89741e-01 +1.19861e-01, h_delta_z_e_prot->GetMaximum());
    c1->Print("fig/proton_delta_z_e_prot.png");
   h_prot_DeltaBeta->Draw();
     l->DrawLine(-3*1.07340e-02,0, -3*1.07340e-02, h_prot_DeltaBeta->GetMaximum());
     l->DrawLine(3*1.07340e-02,0, 3*1.07340e-02, h_prot_DeltaBeta->GetMaximum());
     c1->Print("fig/prot_DeltaBeta.png");

  h_prot_Beta_p->Draw("colz"); c1->Print("fig/prot_Beta_p.png");
  h_prot_Theta_Phi_1->Draw("colz"); c1->Print("fig/prot_Theta_Phi_1.png");
  h_prot_Theta_Phi_2->Draw("colz"); c1->Print("fig/prot_Theta_Phi_2.png");
  h_proton_DC_IC_XY_1->Draw("colz"); c1->Print("fig/proton_DC_IC_XY_1.png");
  h_proton_DC_IC_XY_2->Draw("colz"); c1->Print("fig/proton_DC_IC_XY_2.png");
  h_proton_DC_IC_XY_3->Draw("colz"); c1->Print("fig/proton_DC_IC_XY_3.png");
  h_proton_DC_IC_XY_4->Draw("colz"); c1->Print("fig/proton_DC_IC_XY_4.png");
  h_proton_z->Draw();                              
    l->DrawLine(-77,0,-77,h_proton_z->GetMaximum());
    l->DrawLine(-50,0,-50,h_proton_z->GetMaximum());
    c1->Print("fig/proton_z.png");                                            


  /// rtpc --------------------------------------------------------------------
  float scale_factor = h_rtpc_npd_r->GetEntries()/h_rtpc_npd_l->GetEntries();

 l->SetLineColor(kRed);

   h_dedx_p_Coh_l->Draw("colz"); dedx_leg->Draw();                                                         
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh.png");
   h_dedx_p_Coh_r->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh.png");


   h_dedx_p_Coh_l_1->Draw("colz"); dedx_leg->Draw();                                                  
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_1.png");
         c1->Print("fig/rtpc_l_dedx_p_Coh_1.C");
   h_dedx_p_Coh_r_1->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_1.png");                                                     
        c1->Print("fig/rtpc_r_dedx_p_Coh_1.C");                                                      


   h_dedx_p_Coh_l_0->Draw("colz"); dedx_leg->Draw();
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_0.png");
   h_dedx_p_Coh_r_0->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_0.png");


   h_dedx_p_Coh_l_2->Draw("colz");  dedx_leg->Draw();                                                      
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_2.png");
   h_dedx_p_Coh_r_2->Draw("colz");  dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_2.png");


   h_dedx_p_Coh_l_3->Draw("colz"); dedx_leg->Draw();                                                       
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_3.png");
   h_dedx_p_Coh_r_3->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_3.png");


   h_dedx_p_Coh_l_4->Draw("colz");  dedx_leg->Draw();
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_4.png");         
   h_dedx_p_Coh_r_4->Draw("colz");  dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_4.png");


   h_dedx_p_Coh_l_5->Draw("colz"); dedx_leg->Draw();
        fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
         c1->Print("fig/rtpc_l_dedx_p_Coh_5.png");
   h_dedx_p_Coh_r_5->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
        c1->Print("fig/rtpc_r_dedx_p_Coh_5.png");


   for(int i=1; i<7; i++)
     {
      h_dedx_p_s[i]->Draw("colz"); dedx_leg->Draw();
      fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
      c1->Print(Form("fig/rtpc_dedx_p_s_%d.png",i));
      c1->Print(Form("fig/rtpc_dedx_p_s_%d.C",i));
      }

//  h_rtpc_npd_l->Scale(scale_factor); 
      h_rtpc_npd_r->Draw();    h_rtpc_npd_l->Draw("same");   leg->Draw(); 
      c1->Print("fig/rtpc_npd.png");

//  h_rtpc_sdist_l->Scale(scale_factor);
      h_rtpc_sdist_r->Draw();  h_rtpc_sdist_l->Draw("same"); leg->Draw(); 
      l->DrawLine(-2,0,-2, h_rtpc_sdist_r->GetMaximum());
      l->DrawLine(2,0,2, h_rtpc_sdist_r->GetMaximum());
      c1->Print("fig/rtpc_sdist.png");

//  h_rtpc_npd_l->Scale(scale_factor);
      h_rtpc_edist_r->Draw();  h_rtpc_edist_l->Draw("same"); leg->Draw();
      l->DrawLine(-1,0,-1, h_rtpc_edist_r->GetMaximum());
      l->DrawLine(5,0,5, h_rtpc_edist_r->GetMaximum());
      c1->Print("fig/rtpc_edist.png");

//  h_rtpc_r0_l->Scale(scale_factor);
      h_rtpc_r0_r->Draw();     h_rtpc_r0_l->Draw("same");    leg->Draw();
      c1->Print("fig/rtpc_r0.png");

//  h_rtpc_X2_l->Scale(scale_factor); 
      h_rtpc_X2_r->Draw();     h_rtpc_X2_l->Draw("same");    leg->Draw();
      l->DrawLine(3.5, 0,3.5, h_rtpc_X2_r->GetMaximum());
      c1->Print("fig/rtpc_X2.png");

//  h_rtpc_z_l->Scale(scale_factor);
      h_rtpc_z_r->Draw();      h_rtpc_z_l->Draw("same");     leg->Draw();
      l->DrawLine(-100, 0,-100, h_rtpc_z_r->GetMaximum());
      l->DrawLine( 100, 0, 100, h_rtpc_z_r->GetMaximum());
      c1->Print("fig/rtpc_z.png");

//  h_delta_z_l->Scale(scale_factor); 
      h_delta_z_r->Draw();     h_delta_z_l->Draw("same");
      l->DrawLine(-20, 0,-20, h_delta_z_r->GetMaximum());
      l->DrawLine( 20, 0, 20, h_delta_z_r->GetMaximum());
      leg->Draw();
      c1->Print("fig/rtpc_delta_z.png");

  h_rtpc_phi_theta_1->Draw("colz"); c1->Print("fig/rtpc_phi_theta_1.png");
  h_rtpc_phi_theta_2->Draw("colz"); c1->Print("fig/rtpc_phi_theta_2.png");

  /// photons --------------------------------------------------------------------
 l->SetLineColor(kRed);
  h_EC_XY_photon_1->Draw("colz");   c1->Print("fig/photons_EC_XY_photon_1.png");
  h_EC_XY_photon_2->Draw("colz");   h_EC_XY_photon_1->Draw("same");  c1->Print("fig/photons_EC_XY_photon_2.png");
  h_EC_photon_U->Draw();   
      l->DrawLine(100,0,100, h_EC_photon_U->GetMaximum());  
      l->DrawLine(390,0,390,h_EC_photon_U->GetMaximum()); 
      c1->Print("fig/photons_EC_photon_U.png");
  h_EC_photon_V->Draw();
      l->DrawLine(370,0,370,h_EC_photon_V->GetMaximum());   
      c1->Print("fig/photons_EC_photon_V.png");
  h_EC_photon_W->Draw(); 
      l->DrawLine(390,0,390,h_EC_photon_W->GetMaximum());
      c1->Print("fig/photons_EC_photon_W.png");
  h_eng_photon->Draw(); 
   c1->Print("fig/photons_eng_photon.png");
  h_photon_beta->Draw();   
   l->DrawLine(-3*2.04674e-02+1.00707,0,-3*2.04674e-02+1.00707, h_photon_beta->GetMaximum());
   l->DrawLine(3*2.04674e-02+1.00707,0,  3*2.04674e-02+1.00707, h_photon_beta->GetMaximum());
   c1->Print("fig/photons_beta_EC.png");
  h_IC_XY_center->Draw("colz");   c1->Print("fig/photons_IC_XY_center.png");
  h_ic_theta->Draw();  c1->Print("fig/photons_ic_theta.png");
  h_ic_time->Draw();  c1->Print("fig/photons_ic_time.png");
 
 l->SetLineColor(kBlack);
 c5->cd();
 c5->SetLogz();
  h_IC_theta_En_1->Draw("colz"); 
    l->DrawLine(0.3,6,0.3,14.5);
    l->DrawLine(0.8,2.2,0.3,6.0);   
    c5->Print("fig/photons_IC_theta_En_1.png");
  h_IC_theta_En_2->Draw("colz");       c5->Print("fig/photons_IC_theta_En_2.png");               
  h_IC_XY_2->Draw("colz");  h_IC_XY_1->Draw("same");   c5->Print("fig/photons_IC_XY_2.png");
  h_photon_phi_theta_1->Draw("colz");   c5->Print("fig/photon_phi_theta_1.png");
  h_photon_phi_theta_2->Draw("colz");   c5->Print("fig/photon_phi_theta_2.png");


    for(int i=1; i<7; i++)
      {
       c2->cd(i);
       h_dedx_p_s[i]->Draw("colz"); dedx_leg->Draw();
       fp->Draw("same"); fd->Draw("same"); f3H->Draw("same"); f4He->Draw("same"); f3He->Draw("same");
       }

      c2->Print("fig/rtpc_dedx_p_s.png");                                                            
      c2->Print("fig/rtpc_dedx_p_s.C");                                                            


}
// LOOP ENDS HERE - Frank
//////////////////////////////////////////////////////////////
///      bethe-block
// input momentum (p) in MeV, charge, and mass in MeV
float dedx_cal( float p, int charge, float mass)
{
 float mom=p*charge;
 float beta=mom/sqrt(mom*mom+mass*mass);
 float DENSITY=1.03e-3; // density (g/cm^3)
 float AEFF=126.787;    // effective atomic number (g/mol)
 float ZEFF=66;         // effective atomic charge
 float IEFF=99.794E-6;  // effective ionization energy (MeV)
 float KAPPA=0.307075;  // PDG constant (MeV/g*cm^2)

 float melec=0.511;
   if (beta<0.00001 || beta>0.99999) return 0;
   float gamma=1/sqrt(1-beta*beta);
   float Tmax=2*melec*pow(beta*gamma,2)/(1+2*gamma*melec/mass+pow(melec/mass,2)); // MeV
   float coeff=KAPPA*ZEFF/AEFF*pow(charge/beta,2);
   float bb=DENSITY*coeff*(log(2*melec*pow(beta*gamma,2)*Tmax/pow(IEFF,2))/2-beta*beta); // MeV/cm
   return bb>0 ? 900*bb : 0;

}


/////////////////////////////////////////////////////////////
/// Accp of the RTPC
bool AcceptRTPC( float vz, float theta, float phi)
{
     float rinner=3.0;
     float router=6.0;
     float R2D = 180./TMath::Pi();
     TVector3 vert(0,0,vz+64);

     TVector3 dir(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta));

    // make sure track intersects cathode and gem inside RTPC:
    TVector3 r1=vert+dir.Unit()*(rinner/fabs(sin(dir.Theta())));
    TVector3 r2=vert+dir.Unit()*(router/fabs(sin(dir.Theta())));
    if (fabs(r1.Z())>10) return 0;
    if (fabs(r2.Z())>10) return 0;

    // kill top/bottom support regions:
    float phi2 = GEO.GetPhi(phi,0);
    if (fabs(phi2- 90)<30) return 0;
    if (fabs(phi2-270)<30) return 0;

    // kill if track goes through target holder at upstream end:
    if (FID.RTPCnose(vz,cos(theta))) return 0;

    return 1;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
  /// function to find the error of n_x/n_y ////////////
float Find_Error(float x, float y)
{
 float statistical_error = 0.0;
       statistical_error = (sqrt(x)*y + x*sqrt(y))/pow(y,2);

  return statistical_error;
 }

/////////////////////////////////////////////////////////////////////////////////////////////
// new functions for the helicity
int GetHelicity(const int h10hel,const int runno)
{
    // returns 0/1 if helicity is certain, else -99999

    // if it's negative, the true helicity could not be determined:
    if (h10hel<0) return -99999;

    // if bits 5-24 are zero, this is an incomplete quartet:
    if ( !(h10hel & 0xFFFFF0) ) return -99999;

    // HWP position 0/1 for OUT/IN, else unknown:
    const int hwp=GetHWP(runno);
    if (hwp!=0 && hwp!=1) return -99999;

    return (hwp ^ h10hel) & 1;
}

int GetHWP(const int runno)
{
    // For EG6 Run Numbers, return position of Half-Wave Plate:
    //  1 is IN
    //  0 is OUT
    //  else UNKNOWN

    // before first HWP change, HWP is OUT(0):
    const int initialstate = 0;

    // these are the first runs AFTER a HWP change:
    static const int nchanges=20;
    static const int changes[nchanges]=
    {
        61626,61642,61670,61688,61691,
        61717,61734,61745,61759,61770,
        61791,61793,61802,61819,61836,
        61852,61873,61886,61906,61931
    };

    // FIRST check for bad run numbers where HWP is uncertain:
    if (runno < 61510 || runno > 61930)
    {
        cerr<<"HWP unknown for run "<<runno<<endl;
        return -99999;
    }
    else if (runno == 61733)
    {
        cerr<<"HWP changed in middle of run "<<runno<<endl;
        return -99999;
    }
    else if (runno > 61779 && runno < 61791)
    {
        cerr<<"HWP unknown for run "<<runno<<endl;
        return -99999;
    }

    // THEN find the HWP position:   (should replace with binary search)
    for (int ii=0; ii<nchanges; ii++)
        if (runno < changes[ii])
            return (initialstate ^ ii) & 1;

    return -99999;
}











///////////////////////////////////////////////////////////////////////////////////////////////
                  /// Nathan CLAS vertez correction //
double EG6VertexCorrCLAS(const int runno,
                         const double vz,
                         const double theta,const double phi)
{
    // INPUT:
    // runno = EG6 Run Number
    // vz    = CLAS z-vertex (cm)
    // theta = CLAS theta (rad)
    // phi   = CLAS phi (rad)
    //
    // OUTPUT: 
    // corrected z-vertex (cm)
    //
    static const int nr=4;
    static const int runranges[nr]={61483,61580,61850,99999};
    static const double xx[nr]={0.155, 0.237, 0.27,  0.30}; // x beam position (cm)
    static const double yy[nr]={0.029,-0.040, 0.04, -0.04}; // y beam position (cm)
    if (fabs(sin(theta)) < 1e-3) return -9999;
    int therange=0;
    for (int ii=0; ii<nr; ii++) {
        if (runno < runranges[ii]) {
            therange=ii;
            break;
        }
    }
    const double rr=sqrt(pow(xx[therange],2)+pow(yy[therange],2));
    const double phi0=atan2(yy[therange],xx[therange])+3.141593;
    return vz-rr*cos(phi-phi0)/tan(theta);
}
double EG6VertexCorrCLAS(const int runno,
                         const double vz,
                         const double cx,const double cy,const double cz)
{
    return EG6VertexCorrCLAS(runno,vz,acos(cz),atan2(cy,cx));
}

////////////////////////////////////////////////////////////////////////////////////////////////
                  //// Nathan sampling fraction corrections ///

double EG6ev2file(const int runno,const int evno)
{
//    runno = Run Number (61510-61930)
//    evno  = Event Number
//
//    Returns the interpolated file number.
//    Negative return value is an error.  (e.g. run not found)
//    Reads parameters from EG6ev2file.dat
//
    static bool initialized=0;
    static const int ncols=4;

    // table columns:
    static vector <int> tRUN,tFILE,tEVMIN,tEVMAX;

    // read the data table:
    if (!initialized)
    {
        const char* filename="EG6ev2file.dat";
        FILE *fin=fopen(filename,"r");
        if (!fin)
        {
            fprintf(stderr,"EG6ev2file:  Missing Input File:  %s\n",filename);
            return -9999;
        }

        char buf[256];
        double tmp[ncols];

        while ((fgets(buf,256,fin)) != NULL)
        {
            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ev2file:  Error Reading File:  %s\n",filename);
                    return -9999;
                }
            }

            tRUN .push_back(tmp[0]);
            tFILE.push_back(tmp[1]);
            tEVMIN.push_back(tmp[2]);
            tEVMAX.push_back(tmp[3]);
        }

        initialized=1;
        fclose(fin);
    }
 
    // check for invalid input:
    if ( runno<tRUN[0] || runno > tRUN[tRUN.size()-1])
    {
        fprintf(stderr,"EG6ev2file:  Run Out of Range:  %d\n",runno);
        return -9999;
    }
    if ( evno<0 )
    {
        fprintf(stderr,"EG6ev2file:  Invalid Event #:  %d\n",evno);
        return -9999;
    }

    // binary search:
    int imax=tRUN.size()-1;
    int imin=0;
    int jj=-1;

    while (imax >= imin)
    {
        jj = (imax+imin)/2;

        if      (runno > tRUN[jj]) imin=jj+1;
        else if (runno < tRUN[jj]) imax=jj-1;
        else
        {
            if      (evno > tEVMAX[jj]) imin=jj+1;
            else if (evno < tEVMIN[jj]) imax=jj-1;
            else break;
        }
    }

    if (runno != tRUN[jj])
    {
        fprintf(stderr,"EG6ev2file:  Run Not Found:  %d\n",runno);
        return -9999;
    }

    // interpolate:
    const double nevents=tEVMAX[jj]-tEVMIN[jj]+1;
    return tFILE[jj] + (evno-tEVMIN[jj])/nevents;
}


double EG6ECsampcorr(const int sector,const double rf150)
{
//    sector == Sector (1-6)
//    rf150  == Run# + File# / 150
//
//    Returns the sampling fraction parameterized by piecewise functions of the form:
//    E0 + A0 * ( exp(-ALPHA1*(x-X0)) + exp(-ALPHA2*(x-tX0)) )
//    where x = Run# + File#/150
//    Negative return value is an error.
//    Reads parameters from EG6ECsampcorr.dat
//
    static bool initialized=0;
    static const int npar=5;
    static const int ncols=2+6*npar;

    // table columns:
    static vector <float> tXLO,tXHI;  // <--- Run# + File#/150
    static vector <float> tE0[6],tX0[6],tA0[6],tALPH1[6],tALPH2[6];

    // read the data table:
    if (!initialized)
    {
        const char* filename="EG6ECsampcorr.dat";
        FILE *fin=fopen(filename,"r");
        if (!fin)
        {
            fprintf(stderr,"EG6ECsampcorr:  Missing Input File:  %s\n",filename);
            return -9999;
        }

        char buf[1024];
        double tmp[ncols];

        while ((fgets(buf,1024,fin)) != NULL)
        {
            for (int ii=0; ii<ncols; ii++)
            {
                if (ii==0) tmp[ii] = atof((char*)strtok(buf," "));
                else       tmp[ii] = atof((char*)strtok(NULL," "));
            
                if (isnan(tmp[ii]))
                {
                    fprintf(stderr,"EG6ECsampcorr:  Error Reading File:  %s\n",filename);
                    return -9999;
                }
            }

            tXLO.push_back(tmp[0]);
            tXHI.push_back(tmp[1]);

            for (int sec=0; sec<6; sec++)
            {
                tE0[sec]   .push_back(tmp[2+npar*sec+0]);
                tX0[sec]   .push_back(tmp[2+npar*sec+1]);
                tA0[sec]   .push_back(tmp[2+npar*sec+2]);
                tALPH1[sec].push_back(tmp[2+npar*sec+3]);
                tALPH2[sec].push_back(tmp[2+npar*sec+4]);
            }
        }

        initialized=1;
        fclose(fin);
    }

    // check for invalid input:
    if (sector<1 || sector>6)
    {
        fprintf(stderr,"EG6ECsampcorr:  Invalid sector:  %d\n",sector);
        return -9999;
    }
    if (rf150 < tXLO[0] || rf150 >= tXHI[tXHI.size()-1])
    {
        fprintf(stderr,"EG6ECsampcorr:  Invalid rf150:  %f\n",rf150);
        return -9999;
    }

    // binary search:
    int imax=tXLO.size()-1;
    int imin=0;
    int jj=-1;

    while (imax >= imin)
    {
        jj = (imax+imin)/2;
        if      (rf150 >= tXHI[jj]) imin=jj+1;
        else if (rf150 <  tXLO[jj]) imax=jj-1;
        else break;
    }

    // the parameterization:
    const int ss=sector-1;
    const double corr=
        tE0[ss][jj] +
        tA0[ss][jj] * exp(-tALPH1[ss][jj]*(rf150-tX0[ss][jj])) +
        tA0[ss][jj] * exp(-tALPH2[ss][jj]*(rf150-tX0[ss][jj])) ;

    return corr;
}

double EG6ECsampcorr(const int sector,const int runno,const int evno)
{
//    sector == Sector (1-6)
//    runno  == Run Number (61510-61930)
//    evno   == Event Number
//
    const double file=EG6ev2file(runno,evno);
    if (file<0) return -9999;
    const double rf150 = runno + file/150.0;
//    const double rf150 = runno;
    return EG6ECsampcorr(sector,rf150);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// function to convert the decimal to binary //
int *Decimal_Binary_Converter(int decimal)
 { 
  int remainder, digits = 0, dividend = decimal;
  while(dividend != 0)
   {
    dividend = dividend / 2;
    digits++;
    }
 
  static int *pointer;
  static int array[25]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  pointer = array;
  dividend = decimal;
  for(int i =0; i< digits; i++)
   {
    remainder = dividend % 2;
    array[i] = remainder;
    dividend = dividend / 2;
    }

/*   
    cout << setw(3) << decimal << " in binary is" << setw(3);
    for(int j = 0; j < 25; j++)
     {  cout << array[j];}
    cout << endl;  
*/        
//   cout<<array[0]<<endl;
  return pointer;
 }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
                       /// function to find the helicity ////////
int Find_helicity(int helicity_cor, int event_HWP)
 {
  int Helicity= 0;
  Helicity = (event_HWP^helicity_cor)&1;
  if (Helicity==0) Helicity= -1;
  if (Helicity==1) Helicity= 1;

   return Helicity;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
                     // function to find the half wave plate position //////
int Find_HWP(int runno)
{
    // For EG6 Run Numbers, return position of Half-Wave Plate:
    //  1 is IN
    //  0 is OUT
    //  else UNKNOWN

    // before first HWP change, HWP is OUT(0):
    const int initialstate = 0;

    // these are the first runs AFTER a HWP change:
    static const int nchanges=20;
    static const int changes[nchanges]=
    {
        61626,61642,61670,61688,61691,
        61717,61734,61745,61759,61770,
        61791,61793,61802,61819,61836,
        61852,61873,61886,61906,61931
    };
    
    // FIRST check for bad run numbers where HWP is uncertain:
    if (runno < 61510 || runno > 61930) 
    {
     cerr<<"HWP unknown for run "<<runno<<endl;
     return -99999;
    }
    else if (runno == 61733)
    {
     cerr<<"HWP changed in middle of run "<<runno<<endl;
     return -99999;
    }
    else if (runno > 61779 && runno < 61791) 
    {
     cerr<<"HWP unknown for run "<<runno<<endl;
     return -99999;
    }

    // THEN find the HWP position:
    for (int ii=0; ii<nchanges; ii++)
        if (runno < changes[ii])
            return (initialstate ^ ii) & 1;

    return -99999;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// function to find the transverse momentum of the hadron ///////////
float PT2_H_Find(TLorentzVector h_prim, TLorentzVector Gamma_star)
 {
  float PT2_H= 0.;
  float PT2_L= 0.;
  PT2_L = pow((Gamma_star).Dot(h_prim),2)/((Gamma_star).Dot(Gamma_star));  
  PT2_H = (h_prim).Dot(h_prim) - PT2_L;
     
  return PT2_H;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
           //// function to find t_h which is the transfere momentum from the virtual photon to the hadron //
float T_H_Find(TLorentzVector h, TLorentzVector h_prim)
 {
  float T_H= 0.;
  TLorentzVector Diff;
  Diff.SetPxPyPzE(0., 0., 0., 0.);
  Diff = h_prim - h;
  T_H = (Diff).Dot(Diff);

  return T_H;   
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
           /////// function to fine the phi between the leptonic and Hadronic planes //////////
float PHI_Find( TLorentzVector l, TLorentzVector l_prim, TLorentzVector Gamma_star, TLorentzVector h_prim)
 {
  float PHI_h= 0.;
  float TODEG= 180./TMath::Pi();
  TVector3 LeptonicPlane = (Gamma_star.Vect()).Cross(l_prim.Vect());
  TVector3 HadronicPlane = (h_prim.Vect()).Cross(Gamma_star.Vect());
  PHI_h = LeptonicPlane.Angle(HadronicPlane)*TODEG;
  if(LeptonicPlane.Dot(h_prim.Vect())>0.) PHI_h = 360.0 - PHI_h;

  return PHI_h;
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////
              /// function to remove the middel hot region in DC XY projection /////////////
bool DC_Mid_Cut(float X, float Y)
 {
  bool result= true;
  if( -0.2<X && X<0.6 && -0.2<Y && Y<0.6)
   result= false;
  
   return result;
  }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// IC cut: theta vs. energy ///////////
bool IC_theta_En_Cut (float X, float Y)
 {
  bool result= true;

  float points[6][2];
  points[0][0] = 0.00; points[0][1] = 0.00;      //1  2
  points[1][0] = 0.00; points[1][1] = 15.0;      //
  points[2][0] = 0.3;  points[2][1] = 15.0;     //   3
  points[3][0] = 0.3;  points[3][1] = 6.00;     //
  points[4][0] = 0.80;  points[4][1] =  0.00;    // 0   4     
  points[5][0] = 0.00;   points[5][1] = 0.00;    // the same initial point        

  TCutG *geocut = new TCutG("geocut", 6);
  float cutpoints[2][6] = {{points[0][0], points[1][0], points[2][0], points[3][0], points[4][0], points[5][0]},
                            {points[0][1], points[1][1], points[2][1], points[3][1], points[4][1], points[5][1]}};

  for(Int_t j = 0; j < 6; j++)
   {
    geocut->SetPoint(j, cutpoints[0][j], cutpoints[1][j]);
    }
   if(geocut->IsInside(X, Y))
    result = false;

   return result;
  }

//////////////////////////////////////////////////////////////////////////////////////////
             //////  calculate the invarient mass of two photonos///////////

float CalcTwoPhotonInvarMass(TLorentzVector Phot1, TLorentzVector Phot2)
 {
  float TwoPhotInvMass=(Phot1+Phot2).M(); 
  if(TwoPhotInvMass>0.01 && TwoPhotInvMass<1) return TwoPhotInvMass;
  else return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
              ////// CC Fiducial cut //////////////

int basics_vcrpl(float *r0, float *dir, float *plane_par, float *dist, float *cross_point)
 {
/*
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : crossing of the stright line(R0,d)
C-                         with a plane
C-
C-   Inputs  :   r0(3) - initial point of line
C-               dir(3) - vector direction: r = R0 + s*D
C-               plane_par(3) - array of plane parameters:
C-      plane_par(1)*x + plane_par(2)*y + plane_par(3)*z + 1 = 0
C-
C-   Outputs :   cc_vcrpl =  0 - no cross with the plane.
C-                           1 - cross in positive direction
C-                          -1 - cross in negative direction
C-               dist    =  Distance to the cross point
C-               cross_point(3) =  Cross point coordinates.
C-
C-   Created    23-NOV-1998   Alexander V. Vlassov
C-   Modified
C-
C----------------------------------------------------------------------

*/
  double a,b,t,c,d[3];
  const double un = 1.0000000000;
  const float vsmall = 0.000001;
  int i, ires;

  a = b = c = 0.;
  *dist = 0.;
  for(i=0;i<3;i++)
    c += un*dir[i]*dir[i];
  c = sqrt(c);
  if(c <= vsmall)
    {ires = 0;
      return(ires);}

  for(i=0;i<3;i++)
    d[i] = un*dir[i]/c;
  for(i=0;i<3;i++)
    { a += un*plane_par[i]*d[i]; b += un*plane_par[i]*r0[i]; }
  b += un;
  if(fabs(b) <= vsmall)
    {
    for(i=0;i<3;i++)
      cross_point[i] = r0[i];
    ires = 1;
    }
  else
    {
    if(fabs(a) <= vsmall)
      {
      for(i=0;i<3;i++)
        cross_point[i] = 0.;
      ires = 0;
      }
    else
      {
      t = -b/a;
      for(i=0;i<3;i++)
        cross_point[i] = t*dir[i] + r0[i];
      ires = 1;
      *dist = t;
      if(t < 0.)
        { *dist = -t; ires = -1;}
      }
        }
  return(ires);
  }

 // CC fiducial
bool CC_Fid(Float_t Theta, Float_t Phi)
 {
  bool IsInFid = false;
  while(Phi<-30.)Phi+=60.;
  while(Phi> 30.)Phi-=60.;
  if(TMath::Abs(Phi)<Edge1(Theta) &&TMath::Abs(Phi)>Edge2(Theta)) IsInFid = true;
  return IsInFid;
  }

float Edge1(float Theta)
 {
  return -6.332792e+01+Theta*1.105609e+01-Theta*Theta*6.344957e-01+Theta*Theta*Theta*1.873895e-02
         -2.762131e-04*Theta*Theta*Theta*Theta+1.604035e-06*Theta*Theta*Theta*Theta*Theta;
  }
float Edge2(float Theta)
 {
  if(Theta<43.)Theta=43.;
  return 20.*TMath::Sqrt((Theta-43.)/2.);
  }


/////////////////////////////////////////////////////////////////////////////////////////////
              ///// fiducial cut from the shadow of IC on DC ///////
bool IC_DC_shadow_fiducial (float X, float Y)
 {
  bool result= true; 

  float points[11][2];
  points[0][0] = -11.15; points[0][1] = -26.07;              //
  points[1][0] = -11.15; points[1][1] = -23.1;              //
  points[2][0] = -23.1;  points[2][1] = -12.85;            //      4        5
  points[3][0] = -23.1;  points[3][1] = 11.5;               //
  points[4][0] = -10.3;  points[4][1] = 22.95;              // 3                    6
  points[5][0] = 9.91;   points[5][1] = 22.95;              //
  points[6][0] = 23.73;   points[6][1] = 13.1;             // 2                    7
  points[7][0] = 23.73;   points[7][1] = -12.4 ;           //      1        8
  points[8][0] = 12.3;    points[8][1] = -22.36;           //      0        9
  points[9][0] = 12.3;     points[9][1] = -26.07;
  points[10][0] = -11.15; points[10][1] = -26.07;              // the same initial point 	

  TCutG *geocut = new TCutG("geocut", 10);
  float cutpoints[2][11] = {{points[0][0], points[1][0], points[2][0], points[3][0], points[4][0], points[5][0], points[6][0], points[7][0], points[8][0], points[9][0], points[10][0]},
		            {points[0][1], points[1][1], points[2][1], points[3][1], points[4][1], points[5][1], points[6][1], points[7][1], points[8][1], points[9][1], points[10][1]}};
	
  for(Int_t j = 0; j < 11; j++)
   {
    geocut->SetPoint(j, cutpoints[0][j], cutpoints[1][j]);
    }
   if(geocut->IsInside(X, Y))
    result = false;

   return result;
  }

/////////////////////////////////////////////////////////////////////////////////////////////
                ///// fiducial cut of DC for selecting the electrons in EC /////
bool DCff_e(float X, float Y, int S)
 { 
  bool result= false;
  if( (S==3 || S==4 || S==5 || (Y>X*TMath::Tan(TMath::Pi()*((S-1)/3.-1./9)) && Y<X*TMath::Tan(TMath::Pi()*((S-1)/3.+1./9))))
      && (S==1 || S==2 || S==6 || (Y<X*TMath::Tan(TMath::Pi()*((S-1)/3.-1./9)) && Y>X*TMath::Tan(TMath::Pi()*((S-1)/3.+1./9)))) ) result= true;
  
 return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
                       ///// take out the hot chnnels in IC ///////
bool IC_Not_Hot_Channel(float ic_x,float ic_y)
 {
  bool result= true;
  if (
       (-11.0<ic_x && ic_x<-10.3 && -3.0<ic_y && ic_y<-2.2)||
       (-5.8<ic_x  && ic_x<-5.1 && -8.5<ic_y  && ic_y<-7.9)||
       (-1.7<ic_x && ic_x<-1.1 && -11.3<ic_y && ic_y<-10.7)||
       (-3.0<ic_x && ic_x<-2.3 && -8.5<ic_y && ic_y<-7.9)  ||
       (-7.5<ic_x && ic_x<-6.0 && 10.5<ic_y && ic_y<11.5)   || 
       (-12.8<ic_x && ic_x<-11.5 && -8.5<ic_y && ic_y<-7.5)  ||
//       (3.9<ic_x && ic_x<4.5 && -3.0<ic_y && ic_y<-2.3)    ||
//       (3.9<ic_x && ic_x<4.5 && -0.1<ic_y && ic_y<0.5)     ||
       (3.9<ic_x && ic_x<4.5 && -14.1<ic_y && ic_y<-13.5)  
//       (-0.1<ic_x && ic_x<0.5 && 3.9<ic_y && ic_y<4.5)
      )  
   result= false;

   return result;
  } 

///////////////////////////////////////////////////////////////////////////////////////////////
                   ////// fiducial cuts for selecting the photons in IC /////
/*
 bool ICFiducialCut(float IC_X, float IC_Y)
  {
   bool result= true; // assume initially that the his is in the good region
   float stepX = 1.346;
   float stepY = 1.360;
   float InnerICEdge = 3.25;
   float OuterICEdge = 10.75;
   float absX= TMath::Abs(IC_X);
   float absY= TMath::Abs(IC_Y);
  
  if (
      // check if the hit is inside the inner edge 
      ( ( absX / stepX ) < InnerICEdge                                 &&
        ( absY / stepY ) < InnerICEdge                                 &&
        ( absX / stepX ) + (absY / stepY) < InnerICEdge*TMath::Sqrt(2) &&
        ( absX / stepX ) - (absY / stepY) < InnerICEdge*TMath::Sqrt(2) 
        ) 
      
      ||
      
      // check if the hit is outside the outer edge 
      ( ( absX / stepX ) > OuterICEdge                                 ||
        ( absY / stepY ) > OuterICEdge                                 ||
        ( absX / stepX ) + (absY / stepY) > OuterICEdge*TMath::Sqrt(2) ||
        ( absX / stepX ) - (absY / stepY) > OuterICEdge*TMath::Sqrt(2)  
        )
      ) result= false;
  
  return result;
*/
bool ICFiducialCut(float xx, float yy)
{
// from fx
    // inputs are xc,yc from ICPB
    // this is to reject gammas near the inner/outer edges of the IC
    static const float dx=1.346; // cm
    static const float dy=1.360; // cm
    static const float nin=3.25;
    static const float nout=10.75;
    static const float root2=sqrt(2);

    // INNER:
    if (fabs(xx)/dx <= nin  &&
        fabs(yy)/dy <= nin  &&
        fabs(xx/dx - yy/dy) <= nin*root2 &&
        fabs(xx/dx + yy/dy) <= nin*root2 ) return 0;

    // OUTER:
    if (fabs(xx)/dx >= nout  ||
        fabs(yy)/dy >= nout  ||
        fabs(xx/dx - yy/dy) >= nout*root2 ||
        fabs(xx/dx + yy/dy) >= nout*root2 ) return 0;

    return 1;

  }

///////////////////////////////////////////////////////////////////////////////////////////////
                  /// convert xyz coordinate to uvw coordinate for EC ////
 TVector3 EC_XYZ_UVW(TVector3 xyz)
  {
    // Converts x,y,z EC hit in CLAS coordinate system
    // into u,v,w distances of the EC hit.

   float ex=0.;
   float wy=0.;
   float zd=0.;
   float yu=0.;
   float ve=0.;
   float wu=0.;
   float xi=0., yi=0., zi=0.;
   float ec_phy = 0.;
   float phy = 0.;
   float rot[3][3];

   // Parameters
   float ec_the = 0.4363323;
   float ylow = -182.974;
   float yhi = 189.956;
   float tgrho = 1.95325;
   float sinrho = 0.8901256;
   float cosrho = 0.455715;

   // Variables
   ex = xyz[0];
   wy = xyz[1];
   zd = xyz[2];

   phy = atan2(wy,ex)*57.29578;
   if(phy<0.){phy = phy + 360;}
   phy = phy+30.;
   if(phy>360.){phy = phy-360.;}

   ec_phy = ((Int_t) (phy/60.))*1.0471975;

  rot[0][0] = TMath::Cos(ec_the)*TMath::Cos(ec_phy);
  rot[0][1] = -TMath::Sin(ec_phy);
  rot[0][2] = TMath::Sin(ec_the)*TMath::Cos(ec_phy);
  rot[1][0] = TMath::Cos(ec_the)*TMath::Sin(ec_phy);
  rot[1][1] = TMath::Cos(ec_phy);
  rot[1][2] = TMath::Sin(ec_the)*TMath::Sin(ec_phy);
  rot[2][0] = -TMath::Sin(ec_the);
  rot[2][1] = 0.;
  rot[2][2] = TMath::Cos(ec_the);

   yi = ex*rot[0][0]+wy*rot[1][0]+zd*rot[2][0];
   xi = ex*rot[0][1]+wy*rot[1][1]+zd*rot[2][1];
   zi = ex*rot[0][2]+wy*rot[1][2]+zd*rot[2][2];
   zi = zi-510.32 ;

   yu = (yi-ylow)/sinrho;
   ve = (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
   wu = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;

    TVector3 result(yu,ve,wu);

   return result;
} 

/////////////////////////////////////////////////////////////////////////////////////////////////


