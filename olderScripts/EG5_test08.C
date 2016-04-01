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
	cout << " IN BEGIN " << endl;
	// Set up histograms
	initHists();
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

	vector<Int_t> nullInt, eIndices, pIndices, gECIndices, gICIndices, he4Indices;
	vector<TLorentzVector> nullLV, electrons, protons, photonsEC, photonsIC, helia, pi0sECEC, pi0sECIC, pi0sICIC;					// In case I need all electrons found

	Int_t				iiTheElectron = -999;							// Index of the fastest electron with respect to the vector<TLorentzVector> electrons
	Double_t			eBetaFastest = 0;	
	
	using namespace TMath;
	Int_t stopEvent = 1e6;
	for (Long64_t jentry=0; jentry < nentries && jentry < stopEvent; jentry++) {
		curFile = (char*) ((TChain*)(EG5_test::fChain))->GetFile()->GetName();
		fileNumber = getFileNum(curFile);
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << stopEvent/1000 << "\t k     \t=\t" << (Float_t) jentry/(stopEvent+1)*100 << "\t %" << endl;}
	
		eIndices = nullInt;
		electrons = nullLV;
		eVerticesTrigTimes = nullLV;	
	
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
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
				// DC INFO /////////////////////////////////////////////
				
				TVector3 pDir = TVector3( cx[ipart], cy[ipart], cz[ipart]);
					double pTheta = pDir.Theta();
					double pPhi = pDir.Phi();
					double vzShift = getCorrectedVzShift(pTheta, pPhi);
					if(vzShift <= -999)
					{	continue;	}

				eVz = vz[ipart] - vzShift;

				dcHitPos = TVector3(tl1_x[dc[ipart]-1], tl1_y[dc[ipart]-1], tl1_z[dc[ipart]-1]);
				dcHitDir = TVector3(tl1_cx[dc[ipart]-1], tl1_cy[dc[ipart]-1], tl1_cz[dc[ipart]-1]);
				//	dcX = dcHitPos.X();
				//	dcY = dcHitPos.Y();
					dcTheta = RadToDeg() * dcHitPos.Theta();
					dcPhi = RadToDeg() * dcHitPos.Phi();
				
					TVector3 shiftICtoDC = getICtoDCShift(dcHitPos, dcHitDir);
					TVector3 dcShiftedPos = dcHitPos - shiftICtoDC;

						dcICx = dcShiftedPos.X();
						dcICy = dcShiftedPos.Y();
						dcX = dcICx;		// This is kinda pointless but I need it for now
						dcY = dcICy;		// 	Same with this
				Int_t sector = dc_sect[dc[ipart]-1];
				pp = p[ipart];
				
				////////////////////////////// DC INFO /////////////////

				// CC Info /////////////////////////////////////////////
				
				nPhE = nphe[cc[ipart] - 1];
				ccHitECPos = TVector3(dc_xsc[dc[ipart]-1], dc_ysc[dc[ipart]-1], dc_zsc[dc[ipart]-1]);		// Point of track and EC-plane intersection 
					ccTheta = RadToDeg() * ccHitECPos.Theta(); 
					ccPhi = RadToDeg() * ccHitECPos.Phi();

				////////////////////////////// CC Info /////////////////
			
				// EC Info /////////////////////////////////////////////
				ecHitPos = TVector3(ech_x[ec[ipart]-1],ech_y[ec[ipart]-1],ech_z[ec[ipart]-1]);			// EC Hit Position
					ecX = ecHitPos.X();
					ecY = ecHitPos.Y();

				eO = ec_eo[ec[ipart]-1];
				eI = ec_ei[ec[ipart]-1];
				etot0 = etot[ec[ipart]-1];
					eTot = Max( eTot, eO + eI);
					Double_t corrScale = getCorrectedSampFrac((Int_t) runnb, (Int_t) evntid, (Int_t) sector); 
						eSampFrac = eTot/pp * 0.31/corrScale;
			
				////////////////////////////// EC Info /////////////////

				// Fill all negative tracks
				fillHists( 0, 1, 0 );
				
				if( isElectron(ipart) )
				{ 
					eIndices.push_back(ipart);

					// SC Info /////////////////
					eTrigTime = sc_t[sc[ipart]] - sc_r[sc[ipart]]/30. ; // The electron time is the scintillator time minus the time it takes to traverses the sc (speed of light is 30 cm/ns)	
					TLorentzVector eVertexT = TLorentzVector(vx[ipart], vy[ipart], eVz, eTrigTime); 		// Gets the electron's vertex and trigger time;
					eVerticesTrigTimes.push_back(eVertexT);

					TLorentzVector e4Vector = pp * TLorentzVector( cx[ipart], cy[ipart], cz[ipart], 0);
					e4Vector.SetE(eTot);

					electrons.push_back(e4Vector);

					if( b[ipart] > eBetaFastest ){
						eBetaFastest = b[ipart];
						iiTheElectron = (Int_t) electrons.size() - 1;
					}
/*
yy					///////////////////////////////////////////	Histogram filling stage	/////////////////////////////////////////////////////////////////
					/////////// CUTS REFERENCE
					//string hist2Names[] = {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};
					//string hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi", "fail"};
					hCuts2[1]->Fill(dcTheta,dcPhi);
					hCuts2[2]->Fill(ec_ei[ec[ipart]],ec_eo[ec[ipart]]);
					hCuts2[3]->Fill(p[ipart],eTot/p[ipart]);
					hCuts2[4]->Fill(ec_ei[ec[ipart]]+ec_eo[ec[ipart]],eTot);
					hCuts1[0]->Fill(b[ipart]);
					hCuts1[1]->Fill(p[ipart]);
					hCuts1[2]->Fill(vz[ipart]);
					hCuts1[3]->Fill(dcTheta);
					hCuts1[4]->Fill(dcPhi);
					///////////////////////////////////////////	Histogram filling stage	end	/////////////////////////////////////////////////////////////////
*/
//!
					for( Int_t ii = 0; ii < e_ntests; ii++ ){
						if( eReportCard[ii] > 0 ){
							fillHists( 0, 0 , 1 );	
						}
						else if( eReportCard[ii] == 0 ){
							fillHists( 0, 0 , 0 );	
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
					fillHists( 0, 1, 1);
				}
			}
		//	else{	// Fill all negative tracks (this else goes with the if( isElectron ) block )
		//	}
				// End Electron ID /////////////////////////////////////////////////////////////
			}
		}
	/*	
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

					eTot = Max( etot[ec[ipart]] , ec_ei[ec[ipart]] + ec_eo[ec[ipart]] );
					TLorentzVector proton = p[ipart] * TLorentzVector(cx[ipart], cy[ipart], cz[ipart], 0);
					proton.SetE(eTot);
					protons.push_back(proton);
					// We should probably grab the index, the event number, the run number, etc for future use
					//	TVector3 dcHitPos = TVector3( tl1_x[dc[ipart]], tl1_y[dc[ipart]], tl1_z[dc[ipart]] - eVz);	// Position in layer 1 of the DC	
					// hCuts2[0]->Fill(p[ipart],b[ipart]-getCalcBeta(p[ipart],0.93827));	
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
		}			// Ends loop if electron detected block
		else
		{ 
			/// Where this is is Important!!
			/// cout << "No electrons, we can skip this event " << endl; 
			continue;
		} 
			}				// Ends the good detectors if block
	}					// Ends the loop over events
		*/
}						// Ends Loop function

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::fillHists(Int_t ipid, Int_t iIoAn, Int_t iIO)
{
	// The Electron histograms
	if( ipid == 0 ){

		if( iIoAn == 0 ){
			for( Int_t ii = 0; ii < e_ntests; ii++){
				e_2CutsIO[0][ii][iIO]->Fill(dcICx, dcICy);
				e_2CutsIO[1][ii][iIO]->Fill(dcX, dcY);
				e_2CutsIO[2][ii][iIO]->Fill(ccTheta,ccPhi);
				e_2CutsIO[3][ii][iIO]->Fill(ecX, ecY);
				e_2CutsIO[4][ii][iIO]->Fill(eI, eO);
				e_2CutsIO[5][ii][iIO]->Fill(pp, eSampFrac);
				e_2CutsIO[6][ii][iIO]->Fill(dcTheta, dcPhi);

				e_1CutsIO[0][ii][iIO]->Fill(eVz);
				e_1CutsIO[1][ii][iIO]->Fill(nPhE);
			}
		}
		if( iIoAn == 1 ){
				e_2CutsAllNone[0][iIO]->Fill(dcICx, dcICy);
				e_2CutsAllNone[1][iIO]->Fill(dcX, dcY);
				e_2CutsAllNone[2][iIO]->Fill(ccTheta,ccPhi);
				e_2CutsAllNone[3][iIO]->Fill(ecX, ecY);
				e_2CutsAllNone[4][iIO]->Fill(eI, eO);
				e_2CutsAllNone[5][iIO]->Fill(pp, eSampFrac);
				e_2CutsAllNone[6][iIO]->Fill(dcTheta, dcPhi);

				e_1CutsAllNone[0][iIO]->Fill(eVz);
				e_1CutsAllNone[1][iIO]->Fill(nPhE);
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
	// cout << maxHistLength << endl;
	/*
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
			if( ii % 4 == 3 || ii == h2length - 1 ){
				sprintf(png2Name, "images/hists2_%d.png", ii/4+1);
				c2->Print(png2Name,"png");
			}
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
			if( ii % 4 == 3 || ii == h1length - 1 ){
				sprintf(png1Name, "images/hists1_%d.png", ii/4+1);
				c1->Print(png1Name,"png");
			}
		}
	}i
	*/
	for( Int_t iplot = 0; iplot < TMath::Max(e1length,e2length); iplot++ ){
		for( Int_t itest = 0; itest < (Int_t) eReportCard.size(); itest++ ){
			for( Int_t iIO = 0; iIO < 2; iIO++ ){
				if( iplot <  (Int_t) e2length && (iplot == 1 || iplot == 0) && false){
					if( iIO % 2 == 0 ){
						sprintf(name2,"c2 %d %d", iplot, itest+1);
						sprintf(title2,"Canvas2 %d %d", iplot, itest+1);
						c2 = new TCanvas(name2,title2,1350,400);
						c2->Divide(2,1);
					}	
					c2->cd( iIO % 2 + 1 );
					e_2CutsIO[iplot][itest][iIO]->Draw("colz");
					if( iIO % 2 == 1 || iplot == e2length - 1 ){
						sprintf(png2Name, "images/e_2Cuts_%d_%d.png", iplot, itest+1);
						c2->Print(png2Name,"png");
					}
				}

				if( iplot < (Int_t) e1length && false ){
					if( iIO % 2 == 0 ){
						sprintf(name1,"c1 %d %d", iplot, itest+1);
						sprintf(title1,"Canvas1 %d %d", iplot, itest+1);
						c1 = new TCanvas(name1,title1,1350,400);
						c1->Divide(2,1);
					}
					c1->cd( iIO % 2 + 1 );
					e_1CutsIO[iplot][itest][iIO]->Draw("colz");
					if( iIO % 2 == 1 || iplot == h1length - 1 ){
						sprintf(png1Name, "images/e_1Cuts_%d_%d.png", iplot, itest+1);
						c1->Print(png1Name,"png");
					}
				}
			}
		}
	}
	for( Int_t iplot = 0; iplot < TMath::Max(e1length,e2length); iplot++ ){
			for( Int_t iIO = 0; iIO < 2; iIO++ ){
				if( iplot <  (Int_t) e2length && (iplot == 1 || iplot == 0)){
					sprintf(name2,"c2 %d %d", iplot, iIO);
					sprintf(title2,"Canvas2 %d %d", iplot, iIO);
					c2 = new TCanvas(name2,title2, 675, 400);
					e_2CutsAllNone[iplot][iIO]->Draw("colz");
					sprintf(png2Name, "images/e_2Cuts_%d_%d.png", iplot, 99*iIO);
					c2->Print(png2Name,"png");
				}
				if( iplot < (Int_t) e1length && false){
					sprintf(name1,"c1 %d %d", iplot, iIO);
					sprintf(title1,"Canvas1 %d %d", iplot, iIO);
					c1 = new TCanvas(name1,title1,675,400);
					e_1CutsAllNone[iplot][iIO]->Draw("colz");
					sprintf(png1Name, "images/e_1Cuts_%d_%d.png", iplot, 99*iIO);
					c1->Print(png1Name,"png");
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
	TFile outfile2("eID2.root", "RECREATE");//,"overwrite");
	for( Int_t ii = 0; ii < e2length; ii ++){
		for( Int_t jj = 0; jj < (Int_t) eReportCard.size(); jj ++){
			for( Int_t kk = 0; kk < 2; kk ++){
				e_2CutsIO[ii][jj][kk]->Write(); 
				if( jj == 0 )
				{	e_2CutsAllNone[ii][kk]->Write();	}
			}
		}
	}
	for( Int_t ii = 0; ii < e1length; ii ++){
		for( Int_t jj = 0; jj < (Int_t) eReportCard.size(); jj ++){
			for( Int_t kk = 0; kk < 2; kk ++){
				e_1CutsIO[ii][jj][kk]->Write(); 
				if( jj == 0 )
				{	e_1CutsAllNone[ii][kk]->Write();	}
			}
		}
	}
	outfile2.Close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::printNumEvents()
{
	Long64_t nentries = fChain->GetEntriesFast();
	cout << std::scientific;
	cout  << nentries << endl;
}
