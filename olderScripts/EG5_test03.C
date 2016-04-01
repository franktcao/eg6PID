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

//Bool_t EG5_test::
//Bool_t isEVertexCut(Int_t partIndex);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  < Electron ID Cuts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////   Electron ID Cuts >
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  < Electron ID Cut Functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<vector<float>> readDataTable(string fileName0, const int bufSize,  const int nCols)
{	// Reads data from table, puts them into vectors, checks to see if vectors are not nonsensical, binary searches for certain value, returns the required double 
	//	NEEDS TO BE CLEANED UP!!
   // Declare table columns:
	vector<vector<float>>	outVec;
	fileName = fileName0.c_str();
	if( strcmp(fileName, "EG6ev2file.dat") )
	{
	}
	else if( strcmp(fileName, "EG6ECsampcorr.dat") )
	{
	}
	else {cout << "NOT A GOOD FILE" <<endl; return outVec;}
	
	// Read the data table from file and put them into appropriate vectors:
		FILE *fin=fopen(fileName,"r");
     	if (!fin)
     	{
     	    fprintf(stderr," Missing Input File:  %s\n",fileName);
     	    return outVec;
     	}

		char buf[bufSize];
     	double tmp[nCols];
		// Row by row
		while ((fgets(buf,bufSize,fin)) != NULL)
     	{
			vector<float> row;
			// Element by element
			for(int ii=0; ii < nCols; ii++)
         {
				if (ii==0) tmp[ii] = (double) atof((char*)strtok(buf," "));
            else       tmp[ii] = (double) atof((char*)strtok(NULL," "));
				//if( isnan( tmp[ii]) )
            //{
				//	fprintf(stderr,"Error Reading File:  %s\n",filename);
            //   return;
				//}
				row.push_back((float) tmp[ii]);
			}
			outVec.push_back(row);
		}
   	fclose(fin);
	return outVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::Begin()
{
	// Set up histograms
	initHists();
	// Set up data tables
	//		Initiate the file number interpolator 
	initFileIntrp();
	cout << "INITIALIZED FILE INTERPO " << endl;	
	//		Initiate the sampling fraction corrector
	initSampCorr();
	cout << "INITIALIZED SAMP FRAC CORR " << endl;	
}



void EG5_test::Loop()
{
	using namespace TMath;
	if (fChain == 0) return;
	//const char *currentFN = ((TChain*)(EG5_test::fChain))->GetFile()->GetName();
	//cout << currentFN<< endl;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry < nentries && jentry < 1e5; jentry++) {
//		return;
		int ient = 760000;
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << ient/1000 << "\t k     \t=\t" << (Float_t) jentry/(ient+1)*100 << "\t %" << endl;}

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
		//	if (Cut(ientry) < 0) continue;
		if( goodDetectors() ){
			for( Int_t ipart = 0; ipart < gpart; ipart++){
				TLorentzVector pDir;	//	Current particle's 4-momentum vector	
				pDir.SetXYZT(p[ipart]*cx[ipart],p[ipart]*cy[ipart],p[ipart]*cz[ipart], p[ipart]);	
				myTheta = RadToDeg()*pDir.Theta();
				myPhi = RadToDeg()*pDir.Phi();
				sector = dc_sect[dc[ipart]];
			// EC Info	
				etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];

			// CC Info
				if( isElectron(ipart) )
				{ 
				/////////// CUTS REFERENCE
				//string hist2Names[] = {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};
				//string hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi"};
						
					hCuts2[3]->Fill(p[ipart],etote/p[ipart]);
						
					hCuts2[0]->Fill(p[ipart],b[ipart]);	
					hCuts2[2]->Fill(ec_ei[ec[ipart]],ec_eo[ec[ipart]]);
					hCuts2[4]->Fill(ec_ei[ec[ipart]]+ec_eo[ec[ipart]],etote);

					hCuts1[0]->Fill(b[ipart]);
					hCuts1[1]->Fill(p[ipart]);
					hCuts1[2]->Fill(vz[ipart]);
					// CC Info
					trackECIntPos0 = new TVector3(dc_xsc[dc[ipart]], dc_ysc[dc[ipart]], dc_zsc[dc[ipart]]);		// Point of track and EC-plane intersection 

					Double_t ccTheta0 = RadToDeg()*trackECIntPos0->Theta(); 
					Double_t ccPhi0 = RadToDeg()*trackECIntPos0->Phi();
					hCuts2[1]->Fill(myTheta,myPhi);
					hCuts1[3]->Fill(myTheta);
					hCuts1[4]->Fill(myPhi);
				}
			}
		}
	}
}

void EG5_test::writeFiles()
{
	// Let's open a TFile
	//Write the histogram in the file
	TList * h2list = new TList();
	TList * h1list = new TList();

	using namespace TMath;
	Int_t maxHistLength = (Int_t) Max((Int_t) h2length, (Int_t) h1length);
	for( Int_t ii = 0; ii < maxHistLength; ii ++){
		if(ii < h2length){ h2list->Add(hCuts2[ii]);	}
		if(ii < h1length){ h1list->Add(hCuts1[ii]);	}
	}
	TFile outfile1("hist1.root", "RECREATE");//,"overwrite");
	h2list->Write("h1", TObject::kSingleKey);
	outfile1.Close();

	TFile outfile2("hist2.root", "RECREATE");//,"overwrite");
	h1list->Write("h2", TObject::kSingleKey);
	outfile2.Close();
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
}

void EG5_test::drawHist()
{
	char name2[10], title2[20];
	TCanvas *c2;
	char name1[10], title1[20];
	TCanvas *c1;
	Int_t maxHistLength = (Int_t) TMath::Max((Int_t) h2length, (Int_t) h1length);
	cout << maxHistLength << endl;
	for( Int_t ii = 0; ii < maxHistLength; ii++ ){
		if( ii < (Int_t) h2length ){
			if( ii % 4 == 0){
				sprintf(name2,"c2 %d",ii);
				sprintf(title2,"Canvas2 %d",ii);
				c2 = new TCanvas(name2,title2,600,400);
				c2->Divide(2,2);

			}
			c2->cd(ii%4+1);
			hCuts2[ii]->Draw("colz");
		}
		if( ii < (Int_t) h1length ){
			if( ii % 4 == 0){
				sprintf(name1,"c1 %d",ii);
				sprintf(title1,"Canvas1 %d",ii);
				c1 = new TCanvas(name1,title1,600,400);
				c1->Divide(2,2);
			}
			c1->cd(ii%4+1);
			hCuts1[ii]->Draw("colz");
		}
	}
}

void EG5_test::printNumEvents()
{
	Long64_t nentries = fChain->GetEntriesFast();
	cout << std::scientific;
	cout  << nentries << endl;
}
