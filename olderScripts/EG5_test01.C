#define EG5_test_cxx
#include "EG5_test.h"
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <string>
#include <math.h>
#include "TMath.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
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
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <cmath>
#include "TRandom3.h"

using namespace std;

void EG5_test::isElectron(Int_t particleIndex){
	if( q[particleIndex] < 0)
	if( ){
		return true;
	}	
	return false;
}
void EG5_test::Begin()
{
	Float_t *binformation2[] = {200,	0.25, 4, 	200,0.25,1.1,		200,0,70,	200,-180,55,		200,0,0.6,		200,0,0.5,		200,0,5,		200,0,1,		200,0,1,		200,0,1};
	char *hist2Names[] =  {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};	
	char *hist2Titles[] = {"Beta vs Momentum", "Phi vs Theta", "Eo vs Ei", "Etot/p vs p", "Etot vs Ei+Eo"};
	cout << hist2Titles[0] << endl;
	h2length = sizeof(hist2Names)/sizeof(hist2Names[0]); 
	cout << "h2length is " << (int) h2length << endl; 
	
	Float_t *binformation1[] = {200,	0.25,	1.1,		200, 	0.25, 4,		200,	-100, -20,		200, 	-180, 180,		200,	-180,	180};
	char *hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi"};
	char *hist1Titles[] = {"Beta", "Momentum", "Vertex Z", "Theta", "Phi"};
	cout << hist1Names[0] << endl;
	h1length = sizeof(hist1Names)/sizeof(hist1Names[0]); 
	cout << "h1length is " << (int) h1length << endl; 


	for(int ii = 0; ii < (int) h2length; ii++){ 
	hCuts2[ii] = new TH2D(hist2Names[ii], hist2Titles[ii], (int) binformation2[6*ii+0], (Float_t) binformation2[6*ii+1], (Float_t) binformation2[6*ii+2], (int)binformation2[6*ii+3],(Float_t)binformation2[6*ii+4], (Float_t)binformation2[6*ii+5]);
	  cout << "ok " << ii << " / " << h2length << endl;
	}; 

	for(int ii = 0; ii < (int) h1length; ii++){ 
	  hCuts1[ii] = new TH1D(hist1Names[ii], hist1Titles[ii], (int) binformation1[3*ii+0], (Float_t) binformation1[3*ii+1], (Float_t) binformation1[3*ii+2]);
	  cout << "ok " << ii << " / " << h2length << endl;
	}; 
}

void EG5_test::printNumEvents()
{
	Long64_t nentries = fChain->GetEntriesFast();
	cout << std::scientific;
	cout  << nentries << endl;
}

void EG5_test::Loop()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry < nentries && jentry < 100000; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		int ient = 760000;
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << ient/1000 << "\t k     \t=\t" << (Float_t) jentry/(ient+1)*100 << "\t %" << endl;}
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
			  //TString cat = "==";
			  //for( int idc_part = 0; idc_part < dc_part; idc_part++){
			  //		cout << cat << dc_stat[idc_part] << endl;
			  //		cat += "===========";
			  //	}
		//	if (Cut(ientry) < 0) continue;
// =====			ELECTRON ID 		====
		// Check if all the detectors are collecting good data
		if(stat[gpart] > 0)
		if(dc_stat[dc_part-1] > 0)
		if(ec_stat[ec_part-1] > 0)
		if(sc_stat[sc_part-1] > 0){
		using namespace TMath;
			for( int ipart = 0; ipart < gpart; ipart++){
			// Initial Cuts
				if(q[ipart] < 0)			// Electrons are negatively charged, ya dummy
			// CC Cut
				if(nphe[ipart] > 2)
			// Initial Cuts
				if(p[ipart] > 0.8) 		// Make sure we're not getting delta electrons	
				if(vz[ipart] < -54)
				if(vz[ipart] > -74){ 		// If all good, go and plot stuff, stupid
			// Fiducial Cuts
					TVector3 pos(cx[ipart],cy[ipart],cz[ipart]);
					
					Double_t theta = pos.Theta();
					Double_t phi = pos.Phi();

					Double_t phimax = 0;
					Double_t phicoefs[] = {-63.32792, 11.05609, -0.6344957, 1.873895*pow(10.0,-2), 2.762131*pow(10.0,-2), 1.604035*pow(10.0,-2)};
					for(int ii = 0; ii < sizeof(phicoefs)/sizeof(phicoefs[0]); ii++){phimax += phicoefs[ii]*pow(theta*TMath::Pi()/180,ii);}					
					if(phi > phimax){
							  hCuts2[1]->Fill(theta*180/Pi(),phi*180/Pi());
							  hCuts1[3]->Fill(theta*180/Pi());
							  hCuts1[4]->Fill(phi*180/Pi());
					}	
			// EC Energy Cuts
				// Minimum Energy Required
					Double_t etote = (etot[ipart] < ec_ei[ipart] + ec_eo[ipart]) ? ec_ei[ipart] + ec_eo[ipart] : etot[ipart];
				// Sampling Fraction
					Float_t looseness = 2.5;

					Double_t ucoefs[] = { 0.2560840,  0.0432374, -0.00914180,  0.00081589}; 
					Double_t ocoefs[] = { 0.0572976, -0.0272689,  0.00857600, -0.00097998}; 
					 
					Double_t mu = 0;
					Double_t sigma = 0;
					for(int ii = 0; ii < sizeof(ucoefs)/sizeof(ucoefs[0]); ii++){mu += ucoefs[ii]*pow(p[ipart],ii); sigma += ocoefs[ii]*pow(p[ipart],ii);}					
					if(etote/p[ipart] < mu + looseness * sigma)
					if(etote/p[ipart] > mu - looseness * sigma){ 
						hCuts2[3]->Fill(p[ipart],etote/p[ipart]);
					}	
					hCuts2[0]->Fill(p[ipart],b[ipart]);	
					hCuts2[2]->Fill(ec_ei[ipart],ec_eo[ipart]);
					hCuts2[4]->Fill(ec_ei[ipart]+ec_eo[ipart],etote);

					hCuts1[0]->Fill(b[ipart]);
					hCuts1[1]->Fill(p[ipart]);
					hCuts1[2]->Fill(vz[ipart]);
				}
		/////////// CUTS REFERENCE
		//string hist2Names[] = {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};
		//string hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi"};
			}
		}
	}
}
void EG5_test::writeFiles()
{
	// Let's open a TFile
	TFile outfile("hShort.root","recreate");//,"overwrite");
	//Write the histogram in the file
	Hlist1.Add(hCuts1);
	Hlist2.Add(hCuts2);

	Hlist1.Write();
	Hlist2.Write();
	//	outfile.mkdir("2D Histograms");
//	outfile.cd("2D Histograms");
//	hbvp->Write();
//	h2phivthe->Write();
//	h2etotpvp->Write();

//	outfile.cd("..");
//	outfile.mkdir("1D Histograms");
//	outfile.cd("1D Histograms");
//	h1p->Write();
//	h1b->Write();
//	h1vz->Write();
//	h1theta->Write();
//	h1phi->Write();
	// Close the file
	outfile.Close();
}

void EG5_test::drawHist()
{
	char name[10], title[20];
	TCanvas *c;
	for( int ii = 0; ii < (int) h2length; ii++ ){
		if( ii % 4 == 0){
			sprintf(name,"c%d",ii);
			sprintf(title,"Canvas %d",ii);
			c = new TCanvas(name,title,600,400);
			c->Divide(2,2);
		}
		c->cd(ii%4+1);
		hCuts2[ii]->Draw("colz");
	}
	char name1[10], title1[20];
	TCanvas *c1;
	for( int ii = 0; ii < (int) h1length; ii++ ){
		if( ii % 4 == 0){
			sprintf(name1,"c1%d",ii);
			sprintf(title1,"Canvas1 %d",ii);
			c1 = new TCanvas(name1,title1,600,400);
			c1->Divide(2,2);
		}
		c1->cd(ii%4+1);
		hCuts1[ii]->Draw("colz");
	}
}
