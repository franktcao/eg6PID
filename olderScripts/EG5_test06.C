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
			// NOT NEEDED ANYMORE
			//		Initiate the file number interpolator 
			//initFileIntrp();
			//cout << "INITIALIZED FILE INTERPO " << endl;	
	//		Initiate the sampling fraction corrector
	initSampCorr();
	cout << "INITIALIZED SAMP FRAC CORR " << endl;	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::Loop()
{
	using namespace TMath;
	if (fChain == 0) return;
//	cout << currentFN<< endl;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry < nentries && jentry < 1e6; jentry++) {
		curFile = (char*) ((TChain*)(EG5_test::fChain))->GetFile()->GetName();

		int ient = 760000;
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << ient/1000 << "\t k     \t=\t" << (Float_t) jentry/(ient+1)*100 << "\t %" << endl;}

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
		//	if (Cut(ientry) < 0) continue;
		if( goodDetectors() ){
			eIndex = -999;
			vector<Int_t> eIndices;
			vector<TLorentzVector> electrons;					// In case I need all electrons found
			vector<TLorentzVector> eVerticesTrigTimes;

			Int_t				iiTheElectron;							// Index of the fastest electron with respect to the vector<TLorentzVector> electrons
			Double_t			eBetaFastest = 0;	
			// EVNT BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//		This loops over most of the particles 
			for( Int_t ipart = 0; ipart < gpart; ipart++){
//				TLorentzVector pDir = p[ipart] * TLorentzVector(cx[ipart], cy[ipart], cz[ipart],0);
				if( isElectron(ipart) )
				{ 
					eIndices.push_back(ipart);

					eTrigTime = sc_t[sc[ipart]] - sc_r[sc[ipart]]/30. ; // The electron time is the scintillator time minus the time it takes to traverses the sc (speed of light is 30 cm/ns)	
					TLorentzVector eVertexT = TLorentzVector(vx[ipart], vy[ipart], eVz, eTrigTime); 		// Gets the electron's vertex and trigger time;
					eVerticesTrigTimes.push_back(eVertexT);


					// EC Info	
					etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];
					TLorentzVector e4Vector = p[ipart] * TLorentzVector( cx[ipart], cy[ipart], cz[ipart], 0);
					e4Vector.SetE(etote);
					
					electrons.push_back(e4Vector);
					if( b[ipart] > eBetaFastest ){
						eBetaFastest = b[ipart];
						iiTheElectron = (Int_t) electrons.size() - 1;
					}

					Double_t dcTheta =  e4Vector.Theta();
					Double_t dcPhi =  e4Vector.Phi();
					///////////////////////////////////////////	Histogram filling stage	/////////////////////////////////////////////////////////////////
					/////////// CUTS REFERENCE
					//string hist2Names[] = {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};
					//string hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi", "fail"};
					hCuts2[0].Fill(p[ipart],b[ipart]);	
					hCuts2[2].Fill(ec_ei[ec[ipart]],ec_eo[ec[ipart]]);
					hCuts2[3].Fill(p[ipart],etote/p[ipart]);
					hCuts2[4].Fill(ec_ei[ec[ipart]]+ec_eo[ec[ipart]],etote);
					hCuts1[0].Fill(b[ipart]);
					hCuts1[1].Fill(p[ipart]);
					hCuts1[2].Fill(vz[ipart]);
					// CC Info
					// TVector3 trackECIntPos0(dc_xsc[dc[ipart]], dc_ysc[dc[ipart]], dc_zsc[dc[ipart]]);		// Point of track and EC-plane intersection 
					// Double_t ccTheta0 = TMath::RadToDeg()*trackECIntPos0.Theta(); 
					// Double_t ccPhi0 = TMath::RadToDeg()*trackECIntPos0.Phi();
					hCuts2[1].Fill(dcTheta,dcPhi);
					hCuts1[3].Fill(dcTheta);
					hCuts1[4].Fill(dcPhi);
					///////////////////////////////////////////	Histogram filling stage	end	/////////////////////////////////////////////////////////////////
				}
			}	
			// We got an electron, now to find the rest of the particles
			eIndex = eIndices[iiTheElectron];
			Int_t nElectrons = (Int_t) electrons.size();

			vector<Int_t> pIndices;
			vector<TLorentzVector> protons;					// In case I need all electrons found
			
			if( nElectrons > 0 ){
				for( Int_t ipart = eIndex; ipart < gpart; ipart++){				// Index should be greater since electron triggers all other particles
					if( isProton(ipart) ){
						pIndices.push_back(ipart);

						etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];
						TLorentzVector proton = p[ipart] * TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 0);
						proton.SetE(etote);
						protons.push_back(proton);
						// We should probably grab the index, the event number, the run number, etc for future use
						//	TVector3 dcHitPos = TVector3( tl1_x[dc[ipart]], tl1_y[dc[ipart]], tl1_z[dc[ipart]] - eVz);	// Position in layer 1 of the DC	
					}
					
					vector <TLorentzVector> photonsEC;
					if( isPhotonEC(ipart) ){
						// We should probably grab the index, the event number, 
						//	the run number, etc for future use
						Double_t gECNrg = Max(ec_ei[ec[ipart]] + ec_eo[ec[ipart]], etot[ec[ipart]]);					// From here Nrg means energy
						gECNrg = gECNrg / 0.273 ;																						// Rescale the total energy 
						TLorentzVector gEC4Vector = gECNrg * TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 1.);

						photonsEC.push_back( gEC4Vector );  
					}
					
					// ICHB BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// 	This loops over photons in the IC
				
					vector<TLorentzVector> photonsIC;
					
					for( Int_t ipartIC = 0; ipartIC < icpart; ipartIC++){
						if(isPhotonIC(ipartIC)){
							// We should probably grab the index, the event number, 
							//	the run number, etc for future use
							TVector3 gICVertex = TVector3(ich_x[ipartIC], ich_y[ipartIC], -eVz );		// IC Photon's vertex
							TVector3 gICDir = gICVertex.Unit();														// Normalize vertex
							Double_t gICNrg = etc[ipartIC];														// Reconstructed energy inside IC

							TLorentzVector gIC4Vector;
							gIC4Vector.SetVect(gICDir);															// Sets the first 3 components of the 4-vector to the normalized 3-vector
							gIC4Vector.SetE(1.);																			
							gIC4Vector = gICNrg * gIC4Vector;													// Scale the vector back up by the energy

							photonsIC.push_back( gIC4Vector );
						}
					}
					// ICHB BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
					// So we found the electron and used it to find all photons (in EC and IC). Now that we have the photons, we can use it to find the pions
					vector<vector<TLorentzVector>> allPi0s = getPi0s(photonsEC, photonsIC);
					vector<TLorentzVector> pi0sECEC = allPi0s[0];
					vector<TLorentzVector> pi0sECIC = allPi0s[1];
					vector<TLorentzVector> pi0sICIC = allPi0s[2];

					
					// GCPB BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
					vector<TLorentzVector> helia;				// Multiple helium, no?
					for( Int_t ipartGC = 0; ipartGC < gcpart; ipartGC++){
						if( isHelium4(ipartGC) ){
							TLorentzVector helium = TLorentzVector( px[ipartGC], py[ipartGC], pz[ipartGC], p_tot[ipartGC]);
							helia.push_back(helium);
						}
					}
					// GCPB BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}		// Ends the if detected electron block
			// EVNT BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}			// Ends loop over particles in EVNT bank
		}				// Ends the good detectors if block
	}					// Ends the loop over events
}						// Ends Loop function

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::drawHist()
{
	TH2D cat = hCuts2[0];
	cat.Draw();
	/*
	char name2[10], title2[20];
	TCanvas *c2;
	char name1[10], title1[20];
	TCanvas *c1;
	Int_t maxHistLength = (Int_t) TMath::Max((Int_t) h2length, (Int_t) h1length);
	cout << maxHistLength << endl;
	for( Int_t ii = 0; ii < maxHistLength; ii++ ){
		if( ii <  h2length ){
			if( ii % 4 == 0){
				sprintf(name2,"c2 %d",ii);
				sprintf(title2,"Canvas2 %d",ii);
				c2 = new TCanvas(name2,title2,600,400);
				c2->Divide(2,2);
			}
			c2->cd(ii%4+1);
			hCuts2[ii].Draw("colz");
		}
		if( ii < (Int_t) h1length ){
			if( ii % 4 == 0){
				sprintf(name1,"c1 %d",ii);
				sprintf(title1,"Canvas1 %d",ii);
				c1 = new TCanvas(name1,title1,600,400);
				c1->Divide(2,2);
			}
			hCuts1[ii].Draw("colz");
			c1->cd(ii%4+1);
		}
	}
	*/
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::writeFiles()
{
	// Let's open a TFile
	//Write the histogram in the file
	//TList * h2list = new TList();
	//TList * h1list = new TList();

	using namespace TMath;
	TFile outfile1("hists2.root", "RECREATE");//,"overwrite");
	Int_t maxHistLength = (Int_t) Max((Int_t) h2length, (Int_t) h1length);
	for( Int_t ii = 0; ii < maxHistLength; ii ++){
		if(ii < h2length){ hCuts2[ii].Write(); }
		if(ii < h1length){ hCuts1[ii].Write();	}
	}
	outfile1.Close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::printNumEvents()
{
	Long64_t nentries = fChain->GetEntriesFast();
	cout << std::scientific;
	cout  << nentries << endl;
}
