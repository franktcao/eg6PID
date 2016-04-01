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

void EG5_test::Begin()
{
	// Set up histograms
	initHists();
	//		Initiate the sampling fraction corrector
	initSampCorr();
	cout << "INITIALIZED SAMP FRAC CORR " << endl;	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::Loop()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	vector<Int_t> nullInt;
	vector<TLorentzVector> nullLV;					// In case I need all electrons found

	vector<Int_t> eIndices;
	vector<TLorentzVector> electrons;					// In case I need all electrons found
	vector<TLorentzVector> eVerticesTrigTimes;
	
	vector<Int_t> pIndices;
	vector<TLorentzVector> protons;					// In case I need all electrons found

	vector<Int_t> gECIndices;
	vector<TLorentzVector> photonsEC;
	vector<Int_t> gICIndices;
	vector<TLorentzVector> photonsIC;

	vector<Int_t> he4Indices;				
	vector<TLorentzVector> helia;				// Multiple helium, no?

	vector<TLorentzVector> pi0sECEC;
	vector<TLorentzVector> pi0sECIC;
	vector<TLorentzVector> pi0sICIC;
	
	Int_t				iiTheElectron = -999;							// Index of the fastest electron with respect to the vector<TLorentzVector> electrons
	Double_t			eBetaFastest = 0;	
	
	using namespace TMath;
	for (Long64_t jentry=0; jentry < nentries && jentry < 1e6; jentry++) {
		curFile = (char*) ((TChain*)(EG5_test::fChain))->GetFile()->GetName();
		//cout << "CURRENT FILE IS " << curFile << endl;
		fileNumber = getFileNum(curFile);
		//cout<< "FILE NUMBER " << fileNumber << endl;
		int ient = 760000;
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << ient/1000 << "\t k     \t=\t" << (Float_t) jentry/(ient+1)*100 << "\t %" << endl;}
	
		eIndices = nullInt;
		electrons = nullLV;
		eVerticesTrigTimes = nullLV;	
	
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
		//	if (Cut(ientry) < 0) continue;
//		if( goodDetectors() ){
		eIndex = -999;
		Int_t nElectrons = -999; 
		// EVNT BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//		This loops over most of the particles 
		for( Int_t ipart = 0; ipart < gpart; ipart++ ){
			// Begin Electron ID /////////////////////////////////////////////////////////////
			if( isElectron(ipart) )
			{ 
				eIndices.push_back(ipart);

				// SC Info /////////////////
				eTrigTime = sc_t[sc[ipart]] - sc_r[sc[ipart]]/30. ; // The electron time is the scintillator time minus the time it takes to traverses the sc (speed of light is 30 cm/ns)	
				TLorentzVector eVertexT = TLorentzVector(vx[ipart], vy[ipart], eVz, eTrigTime); 		// Gets the electron's vertex and trigger time;
				eVerticesTrigTimes.push_back(eVertexT);


				// EC Info /////////////////	
				etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];
				TLorentzVector e4Vector = p[ipart] * TLorentzVector( cx[ipart], cy[ipart], cz[ipart], 0);
				e4Vector.SetE(etote);
				Double_t dcTheta = RadToDeg() * e4Vector.Theta();
				Double_t dcPhi = RadToDeg() *  e4Vector.Phi();

				// CC Info /////////////////
				// TVector3 trackECIntPos0(dc_xsc[dc[ipart]], dc_ysc[dc[ipart]], dc_zsc[dc[ipart]]);		// Point of track and EC-plane intersection 
				// Double_t ccTheta0 = TMath::RadToDeg()*trackECIntPos0.Theta(); 
				// Double_t ccPhi0 = TMath::RadToDeg()*trackECIntPos0.Phi();

				electrons.push_back(e4Vector);

				if( b[ipart] > eBetaFastest ){
					eBetaFastest = b[ipart];
					iiTheElectron = (Int_t) electrons.size() - 1;
				}

				///////////////////////////////////////////	Histogram filling stage	/////////////////////////////////////////////////////////////////
				/////////// CUTS REFERENCE
				//string hist2Names[] = {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};
				//string hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi", "fail"};
				hCuts2[1]->Fill(dcTheta,dcPhi);
				hCuts2[2]->Fill(ec_ei[ec[ipart]],ec_eo[ec[ipart]]);
				hCuts2[3]->Fill(p[ipart],etote/p[ipart]);
				hCuts2[4]->Fill(ec_ei[ec[ipart]]+ec_eo[ec[ipart]],etote);
				hCuts1[0]->Fill(b[ipart]);
				hCuts1[1]->Fill(p[ipart]);
				hCuts1[2]->Fill(vz[ipart]);
				hCuts1[3]->Fill(dcTheta);
				hCuts1[4]->Fill(dcPhi);
				///////////////////////////////////////////	Histogram filling stage	end	/////////////////////////////////////////////////////////////////
			}
			// End Electron ID /////////////////////////////////////////////////////////////
		}	
		
		nElectrons = (Int_t) electrons.size();
		// We got an electron, now to find the rest of the particles
		if( nElectrons > 0 ) { 
			eIndex = eIndices[iiTheElectron];
			eVz = eVerticesTrigTimes[iiTheElectron].Z();

			pIndices = nullInt;
			protons = nullLV;					// In case I need all electrons found

			gECIndices = nullInt;
			photonsEC = nullLV;

			for( Int_t ipart = 0; ipart < gpart; ipart++){				// Index should be greater since electron triggers all other particles
				if( isProton(ipart) ){
					pIndices.push_back(ipart);

					etote = Max( etot[ec[ipart]] , ec_ei[ec[ipart]] + ec_eo[ec[ipart]] );
					TLorentzVector proton = p[ipart] * TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 0);
					proton.SetE(etote);
					protons.push_back(proton);
					// We should probably grab the index, the event number, the run number, etc for future use
					//	TVector3 dcHitPos = TVector3( tl1_x[dc[ipart]], tl1_y[dc[ipart]], tl1_z[dc[ipart]] - eVz);	// Position in layer 1 of the DC	
					hCuts2[0]->Fill(p[ipart],b[ipart]-getCalcBeta(p[ipart],0.93827));	
				}
				else if( isPhotonEC(ipart) ){
					// We should probably grab the index, the event number, 
					//	the run number, etc for future use
					gECIndices.push_back(ipart);
					Double_t gECNrg = Max( ec_ei[ec[ipart]] + ec_eo[ec[ipart]], etot[ec[ipart]] );					// From here Nrg means energy
					gECNrg = gECNrg / 0.273 ;																						// Rescale the total energy 
					TLorentzVector photonEC = gECNrg * TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 1.);

					photonsEC.push_back( photonEC );  
				}
				else{ 
					//cout<< "FAIL PHOTON, PROTON " << endl << endl; 
					continue;
					}
			}	
			// ICHB BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// 	This loops over photons in the IC

			photonsIC = nullLV;
			gICIndices = nullInt;

			for( Int_t ipartIC = 0; ipartIC < icpart; ipartIC++){
				if(isPhotonIC(ipartIC)){
					// We should probably grab the index, the event number, 
					//	the run number, etc for future use
					Int_t icHitID = (statc[ipartIC] -statc[ipartIC]%10000)/10000 - 1 ;
					TVector3 gICVertex = TVector3( ich_x[icHitID], ich_y[icHitID], -eVz );		// IC Photon's vertex
					TVector3 gICDir = gICVertex.Unit();														// Normalize vertex
					Double_t gICNrg = etc[ipartIC];														// Reconstructed energy inside IC

					TLorentzVector gIC4Vector;
					gIC4Vector.SetVect(gICDir);															// Sets the first 3 components of the 4-vector to the normalized 3-vector
					gIC4Vector.SetE(1.);																			
					gIC4Vector = gICNrg * gIC4Vector;													// Scale the vector back up by the energy

					gECIndices.push_back(ipartIC);
					photonsIC.push_back( gIC4Vector );
				}
			}
			// ICHB BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// GCPB BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			he4Indices = nullInt;				
			helia = nullLV;				// Multiple helium, no?
			for( Int_t ipartGC = 0; ipartGC < gcpart; ipartGC++){
				if( isHelium4(ipartGC) ){
					TLorentzVector helium = TLorentzVector( px[ipartGC], py[ipartGC], pz[ipartGC], p_tot[ipartGC]);
					helia.push_back(helium);
					he4Indices.push_back(ipartGC);
				}
			}

			// So we found the electron and used it to find all photons (in EC and IC). Now that we have the photons, we can use it to find the pions
			//vector<vector<TLorentzVector>> allPi0s = getPi0s(photonsEC, photonsIC);
			pi0sECEC = nullLV;
			pi0sECIC = nullLV;
			pi0sICIC = nullLV;

			if( (Int_t) photonsEC.size() > 0)
				pi0sECEC =  getPi0sFromPhotons(photonsEC, nullLV);
			else{ cout << "=================== no EC photons yo"  << endl; }

			if( (Int_t) photonsEC.size() > 0)
			if( (Int_t) photonsIC.size() > 0)
					pi0sECIC =  getPi0sFromPhotons(photonsEC, photonsIC);
			if( (Int_t) photonsEC.size() == 0)
			if( (Int_t) photonsIC.size() == 0)
				{ cout << "============================ no photons AT ALL  yo"  << endl; }

			if( (Int_t) photonsIC.size() > 0)
				pi0sICIC =  getPi0sFromPhotons(photonsIC, nullLV);
			else{ cout << "======================================no IC photons yo"  << endl; }
			cout << endl << endl;
			cout << "===================================" << endl;
			cout << "ELECTRONs :  " << nElectrons << " / " << gpart  << endl; 
			cout << "PROTONs :  " << (Int_t) protons.size() << endl;
			cout << "EC PHOTONS's :  " << (Int_t) photonsEC.size() << endl;
			cout << "IC PHOTONS's :  " << (Int_t) photonsIC.size() << " / " << icpart << endl;
			cout << "=======EC pi0's :  " << (Int_t) pi0sECEC.size() << endl;
			cout << "=======ECIC pi0's :  " << (Int_t) pi0sECIC.size() << endl;
			cout << "=======IC pi0's :  " << (Int_t) pi0sICIC.size() << endl;
			cout << "HELIA :  " << (Int_t) helia.size() << " / " << gcpart << endl;
			cout << "===================================" << endl;
			cout << endl;
			cout << endl;
			// GCPB BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// EVNT BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}			// Ends loop if electron detected block
		else { 
			// cout << "No electrons, we can skip this event " << endl; 
			continue;
			} // Where this is is Important!!
		//	}				// Ends the good detectors if block
	}					// Ends the loop over events
}						// Ends Loop function

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::drawHist()
{
	// TH2D *cat = hCuts2[1];
	// cat->Draw("colz");
	char name2[10], title2[20];
	TCanvas *c2;
	char name1[10], title1[20];
	TCanvas *c1;
	Int_t maxHistLength = (Int_t) TMath::Max((Int_t) h2length, (Int_t) h1length);
	// cout << maxHistLength << endl;
	for( Int_t ii = 0; ii < maxHistLength; ii++ ){
		if( ii <  h2length ){
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::writeFiles()
{
	// Let's open a TFile
	//Write the histogram in the file
	//TList * h2list = new TList();
	//TList * h1list = new TList();

	using namespace TMath;
	TFile outfile1("hists1.root", "RECREATE");//,"overwrite");
	Int_t maxHistLength = (Int_t) Max((Int_t) h2length, (Int_t) h1length);
	for( Int_t ii = 0; ii < maxHistLength; ii ++){
		if(ii < h2length){ hCuts2[ii]->Write(); }
		if(ii < h1length){ hCuts1[ii]->Write();	}
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
