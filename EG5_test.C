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
	// cout << " IN BEGIN " << endl;
	// Set up histograms
	initHists();
	cout << "HISTOGRAMS INITIALIZED" << endl;
	//		Initiate the sampling fraction corrector
	initSampCorr();
	cout << "INITIALIZED SAMP FRAC CORR " << endl;	
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::Loop()
{
	//cout << " IN LOOP " << endl;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	// vector<Int_t> nullInt, eIndices, pIndices, gECIndices, gICIndices, he4Indices;
	// vector<TLorentzVector> nullLV, electrons, protons, photonsEC, photonsIC, helia, pi0sECEC, pi0sECIC, pi0sICIC;					// In case I need all electrons found

	Int_t				iiTheElectron = -999;							// Index of the fastest electron with respect to the vector<TLorentzVector> electrons
	Double_t			eBetaFastest = 0;	
	
	using namespace TMath;
	
	
	stopEvent = 1e6;
	stopEvPow = (Int_t) Log10( (Double_t) stopEvent);


	for (Long64_t jentry=0; jentry < nentries && jentry < stopEvent; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
		
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << stopEvent/1000 << "\t k     \t=\t" << (Float_t) jentry/(stopEvent+1)*100 << "\t %" << endl;}
		
		curFile = (char*) ((TChain*)(EG5_test::fChain))->GetFile()->GetName();
		fileNumber = getFileNum(curFile);
	
		vector<Int_t> eIndices;
		vector<TLorentzVector> electrons, eVerticesTrigTimes;
		eIndex = -999;
		Int_t nElectrons = -999;

		// EVNT BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//		This loops over most of the particles 
		for( Int_t ipart = 0; ipart < gpart; ipart++ ){
			// Begin Electron ID /////////////////////////////////////////////////////////////
			if( stat[ipart] > 0 )							// Checks to see if it is a good run
			if( dc_stat[dc[ipart]-1] > 0 )				// Checks to see if DC is good for this run
			if( sc_stat[sc[ipart]-1] > 0 )				// Checks to see if SC is good for this run
			if( ec_stat[ec[ipart]-1] > 10100 )			// Checks to see if EC is good for this run ( that particular number is outer, inner, and total number of hits in the EC)
			if( q[ipart] < 0 )								// Negative track 
			{	
				// Things needed for electron cuts (might not need) ////////////////////////	
				// 	DC INFO //////////////////////////////////////////////////////////////
				
				TVector3 pDir = TVector3( cx[ipart], cy[ipart], cz[ipart]);
					double vzShift = getCorrectedVzShift(ipart);
					if(vzShift <= -999)
					{	continue;	}

				eVz = vz[ipart] - vzShift;

				dcHitPos = TVector3(tl1_x[dc[ipart]-1], tl1_y[dc[ipart]-1], tl1_z[dc[ipart]-1]);
				dcHitDir = TVector3(tl1_cx[dc[ipart]-1], tl1_cy[dc[ipart]-1], tl1_cz[dc[ipart]-1]);
					dcHitPosVzShifted = dcHitPos - TVector3( 0, 0, eVz);
						dcTheta = RadToDeg() * dcHitPosVzShifted.Theta();
						dcPhi = RadToDeg() * dcHitPosVzShifted.Phi();
					TVector3 shiftICtoDC = getICtoDCShift(dcHitPos, dcHitDir);
						TVector3 icPos = dcHitPos - shiftICtoDC;
							icX = icPos.X();
							icY = icPos.Y();
							dcX = icX;		// This is kinda pointless but I need it for now
							dcY = icY;		// 	Same with this
				
				Int_t sector = dc_sect[dc[ipart]-1];
				pp = p[ipart];
				
				/////////////////////////////////////////////////// DC INFO /////////////////

				// 	CC Info ///////////////////////////////////////////////////////////////
				
				nPhE = nphe[cc[ipart] - 1];
				ccHitECPos = TVector3(dc_xsc[dc[ipart]-1], dc_ysc[dc[ipart]-1], dc_zsc[dc[ipart]-1]);		// Point of track and EC-plane intersection 
					ccTheta = RadToDeg() * ccHitECPos.Theta(); 
					ccPhi = RadToDeg() * ccHitECPos.Phi();

				///////////////////////////////////////////////// CC Info /////////////////
			
				// 	EC Info /////////////////////////////////////////////////////////////
				
				ecHitPos = TVector3(ech_x[ec[ipart]-1],ech_y[ec[ipart]-1],ech_z[ec[ipart]-1]);			// EC Hit Position
					ecX = ecHitPos.X();
					ecY = ecHitPos.Y();

				eO = ec_eo[ec[ipart]-1];
				eI = ec_ei[ec[ipart]-1];
				Double_t etot0 = etot[ec[ipart]-1];
					eTot = Max( (Double_t) etot0, eO + eI);
					Double_t corrScale = getCorrectedSampFrac((Int_t) runnb, (Int_t) evntid, (Int_t) sector); 
						eSampFrac = eTot/pp;
						if( corrScale > 0.1) 
						{ eSampFrac = eSampFrac * (0.31/corrScale); }
			
				/////////////////////////////////////////////// EC Info /////////////////

				//////////////////// Things needed for electron cuts (might not need) ///	

				// Fill all negative tracks 
				fillHists( 0, 0, 0);
				
				if( isElectron(ipart) )
				{ 
					eIndices.push_back(ipart);

					// SC Info /////////////////
					eTrigTime = sc_t[sc[ipart]-1] - sc_r[sc[ipart]-1]/30. ; // The electron time is the scintillator time minus the time it takes to traverses the sc (speed of light is 30 cm/ns)	
					TLorentzVector eVertexT = TLorentzVector(vx[ipart], vy[ipart], eVz, eTrigTime); 		// Gets the electron's vertex and trigger time;
					eVerticesTrigTimes.push_back(eVertexT);

					TLorentzVector e4Vector = pp * TLorentzVector( cx[ipart], cy[ipart], cz[ipart], 0);
					e4Vector.SetE(eTot);

					electrons.push_back(e4Vector);

					if( b[ipart] > eBetaFastest ){
						eBetaFastest = b[ipart];
						iiTheElectron = (Int_t) electrons.size() - 1;
					}
					for( Int_t itest = 0; itest < e_ntests; itest++ ){
						if( eReportCard[itest] > 0 ){
							// if passes this test, fill all histograms with this test in it
							fillHists( 0, 1, itest );	
						}
						else if( eReportCard[itest] == 0 ){
							// if this fails, fill all other histograms except the ones with this test
							fillHists( 0, 2, itest );	
						}
						else{
							cout << "DIDN'T PASS THE GOOD DETECTORS " << endl;
							continue;
						}
					}
				}
				
				vector<Float_t> allPass(e_ntests,1);
				if( eReportCard == allPass )
				{
					fillHists( 0, 3, 0);
				}
			}
			////////////////////////////////////////////////////////////// Electron ID ///////
		}
		/////////////////////////////////////////	This loops over most of the particles 	////
		
		nElectrons = (Int_t) electrons.size();
		// We got an electron, now to find the rest of the particles
		if( nElectrons > 0 ) { 
			eIndex = eIndices[iiTheElectron];
			eVz = eVerticesTrigTimes[iiTheElectron].Z();
			Double_t eTime = eVerticesTrigTimes[iiTheElectron].T(); 
			
			for( Int_t ipart = 0; ipart < gpart; ipart++){				// Index should be greater since electron triggers all other particles
				if( 	b[ipart] != 0					)	  			// Nonzero velocity 
				if( 	stat[ipart] > 0				)				// Checks to see if it is a good run
				if( 	dc_stat[dc[ipart]-1] > 0	)				// Checks to see if DC is good for this run
				if( ec_stat[ec[ipart]-1] > 10100 )			// Checks to see if EC is good for this run ( that particular number is outer, inner, and total number of hits in the EC)
				if( 	sc_stat[sc[ipart]-1] > 0	){				// Checks to see if SC is good for this run
					TLorentzVector particleDir = TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 0);
					TLorentzVector particleP4V = pp * particleDir; 
					eTot = Max( etot[ec[ipart]-1] , ec_ei[ec[ipart]-1] + ec_eo[ec[ipart]-1] );
					particleP4V.SetE(eTot);

					dcHitPos = TVector3(tl1_x[dc[ipart]-1], tl1_y[dc[ipart]-1], tl1_z[dc[ipart]-1]);
					dcHitDir = TVector3(tl1_cx[dc[ipart]-1], tl1_cy[dc[ipart]-1], tl1_cz[dc[ipart]-1]);
					dcHitPosVzShifted = dcHitPos - TVector3( 0, 0, eVz);
					dcTheta = RadToDeg() * dcHitPosVzShifted.Theta();
					dcPhi = RadToDeg() * dcHitPosVzShifted.Phi();
					TVector3 shiftICtoDC = getICtoDCShift(dcHitPos, dcHitDir);
					TVector3 icPos = dcHitPos - shiftICtoDC;
					dcX = icPos.X();
					dcY = icPos.Y();

					pp = p[ipart];
					dB = b[ipart] - getCalcBeta(pp, 0.93827);

					// Fill all positive tracks
					if( q[ipart] > 0 )
					{	
						fillHists( 1, 0, 0 );	
					}


					if( isProton(ipart) ){
						pIndices.push_back(ipart);
						protons.push_back(particleP4V);
						// We should probably grab the index, the event number, the run number, etc for future use
						// Fill Histograms

						for( Int_t itest = 0; itest < p_ntests; itest++ ){
							if( pReportCard[itest] > 0 ){
								// if passes this test, fill all histograms with this test in it
								fillHists( 1, 1, itest );
							}
							else if( pReportCard[itest] == 0 ){
								// if this fails, fill all other histograms except the ones with this test
								fillHists( 1, 2, itest );	
							}
						}
						vector<Float_t> allPass(p_ntests,1);
						if( pReportCard == allPass )
						{
							fillHists( 1, 3, 0);
						}
					}
					/*
						else if( isPhotonEC(ipart) ){
					// We should probably grab the index, the event number, 
					//	the run number, etc for future use
					gECIndices.push_back(ipart);

					Double_t gECNrg = eTot / 0.273 ;																						// Rescale the total energy 
					TLorentzVector photonEC = particleDir; 
					photonEC.SetE(1.);
					photonEC = gECNrg * photonEC;

					photonsEC.push_back( particleP4V );  
					// Fill Histograms
					}
					else{
					// No proton or EC photon for this electron
					//cout<< "FAIL PHOTON, PROTON " << endl << endl; 
					continue;
					}
					 */
				}
			}
				// ICHB BANK PARTICLE ID BEGIN ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// 	This loops over photons in the IC
	/*
			//photonsIC = nullLV;
			//gICIndices = nullInt;
			vector<Int_t> gICIndices;
			vector<TLorentzVector> photonsIC;

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
			//he4Indices = nullInt;				
			//helia = nullLV;				
			vector<Int_t> he4Indices;
			vector<TLorentzVector> helia;	// Multiple helium, no?
			for( Int_t ipartGC = 0; ipartGC < gcpart; ipartGC++){
				if( isHelium4(ipartGC) ){
					TLorentzVector helium = TLorentzVector( px[ipartGC], py[ipartGC], pz[ipartGC], p_tot[ipartGC]);
					helia.push_back(helium);
					he4Indices.push_back(ipartGC);
				}
			}

			// So we found the electron and used it to find all photons (in EC and IC). Now that we have the photons, we can use it to find the pions
			//vector<vector<TLorentzVector>> allPi0s = getPi0s(photonsEC, photonsIC);
			//pi0sECEC = nullLV;
			//pi0sECIC = nullLV;
			//pi0sICIC = nullLV;
			vector<TLorentzVector> pi0sECEC, pi0sECIC, pi0sICIC;	// pi0's in each "topology"
			vector<TLorentzVector> nullLV;	// null LorentzVector 
			// Event by event report
			/// *
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
			// *  /
			// GCPB BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// EVNT BANK PARTICLE ID END ////////////////////////////////////////////////////////////////////////////////////////////////////////////
	*/	
		}			// Ends loop if electron detected block
		else
		{ 
			/// Where this is is Important!!
			/// cout << "No electrons, we can skip this event " << endl; 
			continue;
		} 
	}					// Ends the loop over events
}						// Ends Loop function

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::fillHists(Int_t ipid, Int_t iIO, Int_t itest)
{
	// The Electron histograms
	if( ipid == 0 ){
		if( iIO == 2 ){
			for( Int_t ii = 0; ii < e_ntests ; ii++){
				if( ii != itest){
					e_2Cuts[0][iIO][ii]->Fill(icX, icY);
					e_2Cuts[1][iIO][ii]->Fill(dcX, dcY);
					e_2Cuts[2][iIO][ii]->Fill(ccTheta, ccPhi);
					e_2Cuts[3][iIO][ii]->Fill(ecX, ecY);
					e_2Cuts[4][iIO][ii]->Fill(eI, eO);
					e_2Cuts[5][iIO][ii]->Fill(pp, eSampFrac);
					e_2Cuts[6][iIO][ii]->Fill(dcTheta,dcPhi);

					e_1Cuts[0][iIO][ii]->Fill(eVz);
					e_1Cuts[1][iIO][ii]->Fill(nPhE);
					e_1Cuts[2][iIO][ii]->Fill(dcPhi);
				}
			}
		}
		if( iIO == 1 ){
			e_2Cuts[0][iIO][itest]->Fill(icX, icY);
			e_2Cuts[1][iIO][itest]->Fill(dcX, dcY);
			e_2Cuts[2][iIO][itest]->Fill(ccTheta, ccPhi);
			e_2Cuts[3][iIO][itest]->Fill(ecX, ecY);
			e_2Cuts[4][iIO][itest]->Fill(eI, eO);
			e_2Cuts[5][iIO][itest]->Fill(pp, eSampFrac);
			e_2Cuts[6][iIO][itest]->Fill(dcTheta,dcPhi);

			e_1Cuts[0][iIO][itest]->Fill(eVz);
			e_1Cuts[1][iIO][itest]->Fill(nPhE);
			e_1Cuts[2][iIO][itest]->Fill(dcPhi);
		}
		else{
			e_2Cuts[0][iIO][0]->Fill(icX, icY);
			e_2Cuts[1][iIO][0]->Fill(dcX, dcY);
			e_2Cuts[2][iIO][0]->Fill(ccTheta, ccPhi);
			e_2Cuts[3][iIO][0]->Fill(ecX, ecY);
			e_2Cuts[4][iIO][0]->Fill(eI, eO);
			e_2Cuts[5][iIO][0]->Fill(pp, eSampFrac);
			e_2Cuts[6][iIO][0]->Fill(dcTheta,dcPhi);

			e_1Cuts[0][iIO][0]->Fill(eVz);
			e_1Cuts[1][iIO][0]->Fill(nPhE);
			e_1Cuts[2][iIO][0]->Fill(dcPhi);
		}
	}
	if( ipid == 1 ){
			if( iIO == 2 ){
				for( Int_t ii = 0; ii < p_ntests ; ii++){
					if( ii != itest){
						p_2Cuts[0][iIO][ii]->Fill(dcX, dcY);
						p_2Cuts[1][iIO][ii]->Fill(pp, dB);
						p_2Cuts[2][iIO][ii]->Fill(dcTheta,dcPhi);

						p_1Cuts[0][iIO][ii]->Fill(dB);
						p_1Cuts[1][iIO][ii]->Fill(dVz);
					}
				}
			}
			if( iIO == 1 ){
				p_2Cuts[0][iIO][itest]->Fill(dcX, dcY);
				p_2Cuts[1][iIO][itest]->Fill(pp, dB);
				p_2Cuts[2][iIO][itest]->Fill(dcTheta,dcPhi);

				p_1Cuts[0][iIO][itest]->Fill(dB);
				p_1Cuts[1][iIO][itest]->Fill(dVz);
			}
			else{
				p_2Cuts[0][iIO][0]->Fill(dcX, dcY);
				p_2Cuts[1][iIO][0]->Fill(pp, dB);
				p_2Cuts[2][iIO][0]->Fill(dcTheta,dcPhi);

				p_1Cuts[0][iIO][0]->Fill(dB);
				p_1Cuts[1][iIO][0]->Fill(dVz);
			}
		}
}

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
	char png2Name[20];
	char png1Name[20];
	// Draw Protons
	for( Int_t iplot = 0; iplot < (Int_t) p_2Cuts.size(); iplot++ ){
			for( Int_t iIO = 0; iIO < 4; iIO++ ){
				Int_t ntests = p_ntests;
				for( Int_t itest = 0; itest < ntests; itest++ ){
					if( iIO  != 2){
						Int_t itestName = itest+1;
						Int_t nplotsPcanv = 2;		// number of plots per canvas
						if( iIO == 0 || iIO == 3)
						{ ntests = 1; nplotsPcanv =1; itestName = iIO*33;}

						sprintf(name2,"c2 %d %d", iplot,itestName);
						sprintf(title2,"Canvas2 %d %d ", iplot, itestName);
						c2 = new TCanvas(name2,title2,675*nplotsPcanv,400);
						c2->Divide(nplotsPcanv,1);
						c2->cd( 1 );
						p_2Cuts[iplot][iIO][itest]->Draw("colz");
						if( iIO == 1 ){
							c2->cd(2);
							p_2Cuts[iplot][iIO+1][itest]->Draw("colz");
						}
						sprintf(png2Name, "images/1e%d/p_2Cuts_%d_%d.png", stopEvPow, iplot, itestName);
						c2->Print(png2Name,"png");
					}
				}
				}
			}
	for( Int_t iplot = 0; iplot < (Int_t) p_1Cuts.size(); iplot++ ){
			for( Int_t iIO = 0; iIO < 4; iIO++ ){
				Int_t ntests = p_ntests;
				for( Int_t itest = 0; itest < ntests; itest++ ){
					if( iIO  != 2){
						Int_t itestName = itest+1;
						Int_t nplotsPcanv = 2;		// number of plots per canvas
						if( iIO == 0 || iIO == 3)
						{ ntests = 1; nplotsPcanv =1; itestName = iIO*33;}

						sprintf(name2,"c1 %d %d", iplot,itestName);
						sprintf(title2,"Canvas1 %d %d ", iplot, itestName);
						c2 = new TCanvas(name2,title2,675*nplotsPcanv,400);
						c2->Divide(nplotsPcanv,1);
						c2->cd( 1 );
						p_1Cuts[iplot][iIO][itest]->Draw("colz");
						if( iIO == 1 ){
							c2->cd(2);
							p_1Cuts[iplot][iIO+1][itest]->Draw("colz");
						}
						sprintf(png2Name, "images/1e%d/p_1Cuts_%d_%d.png", stopEvPow, iplot, itestName);
						c2->Print(png2Name,"png");
					}
				}
			}
	}
	// Draw Electrons
	for( Int_t iplot = 0; iplot < (Int_t) e_2Cuts.size(); iplot++ ){
			for( Int_t iIO = 0; iIO < 4; iIO++ ){
				Int_t ntests = e_ntests;
				for( Int_t itest = 0; itest < ntests; itest++ ){
					if( iIO  != 2){
						Int_t itestName = itest+1;
						Int_t nplotsPcanv = 2;		// number of plots per canvas
						if( iIO == 0 || iIO == 3)
						{ ntests = 1; nplotsPcanv =1; itestName = iIO*33;}

						sprintf(name2,"c2 %d %d", iplot,itestName);
						sprintf(title2,"Canvas2 %d %d ", iplot, itestName);
						c2 = new TCanvas(name2,title2,675*nplotsPcanv,400);
						c2->Divide(nplotsPcanv,1);
						c2->cd( 1 );
						e_2Cuts[iplot][iIO][itest]->Draw("colz");
						if( iIO == 1 ){
							c2->cd(2);
							e_2Cuts[iplot][iIO+1][itest]->Draw("colz");
						}
						sprintf(png2Name, "images/1e%d/e_2Cuts_%d_%d.png", stopEvPow, iplot, itestName);
						c2->Print(png2Name,"png");
					}
				}
				}
			}
	for( Int_t iplot = 0; iplot < (Int_t) e_1Cuts.size(); iplot++ ){
			for( Int_t iIO = 0; iIO < 4; iIO++ ){
				Int_t ntests = e_ntests;
				for( Int_t itest = 0; itest < ntests; itest++ ){
					if( iIO  != 2){
						Int_t itestName = itest+1;
						Int_t nplotsPcanv = 2;		// number of plots per canvas
						if( iIO == 0 || iIO == 3)
						{ ntests = 1; nplotsPcanv =1; itestName = iIO*33;}

						sprintf(name2,"c1 %d %d", iplot,itestName);
						sprintf(title2,"Canvas1 %d %d ", iplot, itestName);
						c2 = new TCanvas(name2,title2,675*nplotsPcanv,400);
						c2->Divide(nplotsPcanv,1);
						c2->cd( 1 );
						e_1Cuts[iplot][iIO][itest]->Draw("colz");
						if( iIO == 1 ){
							c2->cd(2);
							e_1Cuts[iplot][iIO+1][itest]->Draw("colz");
						}
						sprintf(png2Name, "images/1e%d/e_1Cuts_%d_%d.png", stopEvPow, iplot, itestName);
						c2->Print(png2Name,"png");
					}
				}
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

	using namespace TMath;/*
	TFile outfile1("hists1.root", "RECREATE");//,"overwrite");
	Int_t maxHistLength = (Int_t) Max((Int_t) h2length, (Int_t) h1length);
	for( Int_t ii = 0; ii < maxHistLength; ii ++){
		if(ii < h2length){ hCuts2[ii]->Write(); }
		if(ii < h1length){ hCuts1[ii]->Write();	}
	}
	outfile1.Close();
*/
	TFile outfile2("eID.root", "RECREATE");//,"overwrite");
	for( Int_t ii = 0; ii < e_n2plots; ii ++){
		for( Int_t jj = 0; jj < 4; jj ++){
			Int_t ntests = e_ntests;
			if( jj == 0 || jj == 3)
			{	ntests = 1;	}
			for( Int_t kk = 0; kk < ntests; kk ++){
				e_2Cuts[ii][jj][kk]->Write(); 
			}
		}
	}
	for( Int_t ii = 0; ii < e_n1plots; ii ++){
		for( Int_t jj = 0; jj < 4; jj ++){
			Int_t ntests = e_ntests;
			if( jj == 0 || jj == 3)
			{	ntests = 1;	}
			for( Int_t kk = 0; kk < ntests; kk ++){
				e_1Cuts[ii][jj][kk]->Write(); 
			}
		}
	}
	outfile2.Close();
	TFile outfile3("pID.root", "RECREATE");//,"overwrite");
	for( Int_t ii = 0; ii < p_n2plots; ii ++){
		for( Int_t jj = 0; jj < 4; jj ++){
			Int_t ntests = p_ntests;
			if( jj == 0 || jj == 3)
			{	ntests = 1;	}
			for( Int_t kk = 0; kk < ntests; kk ++){
				p_2Cuts[ii][jj][kk]->Write(); 
			}
		}
	}
	for( Int_t ii = 0; ii < p_n1plots; ii ++){
		for( Int_t jj = 0; jj < 4; jj ++){
			Int_t ntests = p_ntests;
			if( jj == 0 || jj == 3)
			{	ntests = 1;	}
			for( Int_t kk = 0; kk < ntests; kk ++){
				p_1Cuts[ii][kk][jj]->Write(); 
			}
		}
	}
	outfile3.Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

