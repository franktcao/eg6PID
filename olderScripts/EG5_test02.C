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

Bool_t EG5_test::goodDetectors(){
	if(dc_stat[dc_part] > 0)
	if(ec_stat[ec_part] > 0)
	if(sc_stat[sc_part] > 0){
		return true;
	}
	return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  < Electron ID Cuts
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t EG5_test::isE_VertexCut(Int_t partIndex){
//	 return true;
	// Now to correct for z-vertex since beam center changes from run to run
	TVector3 pDir = TVector3( cx[partIndex], cy[partIndex], cz[partIndex]);
	double pTheta = pDir.Theta();
	double pPhi = pDir.Phi();
	double vzShift = getCorrectedVzShift(pTheta, pPhi);
	if(vzShift <= -999){
		return false;
	}
	vz[partIndex] = vz[partIndex] - vzShift;
	if(vz[partIndex] < -54.0)
	if(vz[partIndex] > -74.0){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t EG5_test::isE_ECFiducialCut(Int_t partIndex){
	return true;
   ////////Set the EC_XYZ to EC_UVW  ---------------------------------
	TVector3 ecHitXYZ;
	ecHitXYZ.SetXYZ(ech_x[ec[partIndex]],ech_y[ec[partIndex]],ech_z[ec[partIndex]]);			// EC Hit Position
   TVector3 ecHitUVW;
	ecHitUVW = getUVWfromXYZ(ecHitXYZ);																		// Convert XYZ coordinates to UVW

	float uu = ecHitUVW(0);
   float vv = ecHitUVW(1);
   float ww = ecHitUVW(2);
	if( uu > 60.	)
	if( uu < 350.	)
	if( vv < 370.	)
	if( ww < 400.	){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t EG5_test::isE_ECEnergyCut(Int_t partIndex){
	return true;
	// EC Energy Cuts
	// Minimum Energy Required
	// Sampling Fraction
	if( ec_eo[ec[partIndex]] > 0 )				//	Energy Positive
	if( ec_ei[ec[partIndex]] > 0 )				//		(in GeV)
	if( etot[ec[partIndex]]	 > 0 )
	{
		Float_t looseness = 2.5;

		Double_t ucoefs[] = { 0.2560840,  0.0432374, -0.00914180,  0.00081589}; 
		Double_t ocoefs[] = { 0.0572976, -0.0272689,  0.00857600, -0.00097998}; 
		
		Double_t mu = 0;
		Double_t sigma = 0;
		for(int ii = 0; ii < sizeof(ucoefs)/sizeof(ucoefs[0]); ii++){mu += ucoefs[ii]*pow(p[partIndex],ii); sigma += ocoefs[ii]*pow(p[partIndex],ii);}					
		Double_t sampFrac = etote/p[partIndex];
		Int_t sector = dc_sect[dc[partIndex]];
		sampFrac = getCorrectedSampFrac(runnb, evntid, sector);	
		if(sampFrac < mu + looseness * sigma)
		if(sampFrac > mu - looseness * sigma){
			return true;
		}
		return false;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t EG5_test::isE_DCFiducialCut(Int_t partIndex){
//	return true;
	// DC Fiducial Cut
	TVector3 dcTrackPos = TVector3(tl1_x[partIndex], tl1_y[partIndex], tl1_z[partIndex]);
	TVector3 dcTrackDir = TVector3(tl1_cx[partIndex], tl1_cy[partIndex], tl1_cz[partIndex]);
	TVector3 shift = getICtoDCShift(dcTrackPos, dcTrackDir);

	TVector3 shiftedPos = dcTrackPos - shift;
	
	Double_t X = shiftedPos.X();
	Double_t Y = shiftedPos.Y();
	Int_t		S = dc_sect[dc[partIndex]];

	if( !isInsideIConDCShadow(X, Y) ){											// Checks to see if the shifted position is in the IC shadow
		Double_t goodRelAngL  = ((S-1)-1./3)*Pi()/3;
		Double_t goodRelAngR = ((S-1)+1./3)*Pi()/3;
		if( S==3 || S==4 || S==5 || Y>X*Tan(goodRelAngL) )
		if( Y<X*Tan(goodRelAngR) ){
			return true;
		}
		if( S==1 || S==2 || S==6 || Y<X*Tan(goodRelAngL) )
		if( Y>X*Tan(goodRelAngR) ){
			return true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isInsideIConDCShadow(Double_t X, Double_t Y)
{
	///// NOTE: This cut requires the shifted DC track X and Y positions. 
	///// fiducial cut from the shadow of IC on DC ///////
	          //      4        5
	          //   
	          // 3                    6
	          //   
	          // 2                    7
	          //      1        8
	          //      0        9
	          // the same initial point  
	TCutG *geoCut = new TCutG("geocut", 10); 
	// First row of cutPoints is x-coordinates, second is y-coordinates
	float cutPoints[2][11] = {	{-11.15, -11.15, -23.1, -23.1, -10.3, 9.91, 23.73, 23.73, 12.3, 12.3, -11.15},
										{-26.07, -23.1, -12.85, 11.5, 22.95, 22.95, 13.1, -12.4, -22.36, -26.07, -26.07}	};
	
	for(Int_t ii = 0; ii < 11; ii++) 
	{
		geoCut->SetPoint(ii, cutPoints[0][ii], cutPoints[1][ii]);
	}
	if(geoCut->IsInside(X, Y)){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t EG5_test::isE_CCFiducialCut(Int_t partIndex){
//	return true;
	// CC Fiducial Cuts
	TVector3 trackECIntPos = TVector3(dc_xsc[dc[partIndex]], dc_ysc[dc[partIndex]], dc_zsc[dc[partIndex]]);		// Point of track and EC-plane intersection 
	TVector3 trackECIntDir = TVector3(dc_cxsc[dc[partIndex]], dc_cysc[dc[partIndex]], dc_czsc[dc[partIndex]]);		// Direction of track and EC-plane intersection 
	
	TVector3 ccThetaPhi = getCCThetaPhi(trackECIntPos, trackECIntDir);		//	Gets CC Theta and Phi in Degrees
	if(ccThetaPhi.X() <= -999){
		return false;
	}

	Double_t ccTheta = ccThetaPhi.X();		//	In Degrees
	Double_t ccPhi = ccThetaPhi.Y();			// In Degrees

	//	Go to relative phi
	while(ccPhi<-30.)ccPhi+=60.;			
   while(ccPhi> 30.)ccPhi-=60.;
	// Make sure cc's phi is in fiducial region. 
	//		phiCoefs are coefficients of 5th degree polynomial fit from Mohammad's thesis
	Double_t phiCoefs[] = {-63.32792, 11.05609, -0.6344957, 1.873895*pow(10.0,-2), 2.762131*pow(10.0,-2), 1.604035*pow(10.0,-2)};
	Double_t phiEdge1 = 0.;
	Double_t phiEdge2 = 0.;
	for(int ii = 0; ii < sizeof(phiCoefs)/sizeof(phiCoefs[0]); ii++){phiEdge1 += phiCoefs[ii]*pow(ccTheta,ii);}					
	if( ccTheta > 43. ){ 
		phiEdge2 = 20*Sqrt(0.5*(ccTheta - 43.)); }

	if( fabs(ccPhi) < phiEdge1 )
	if( fabs(ccPhi) > phiEdge2 ){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////   Electron ID Cuts >
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  < Electron ID Cut Functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<vector<float>> readDataTable(char* fileName, const int bufSize,  const int nCols)
{	// Reads data from table, puts them into vectors, checks to see if vectors are not nonsensical, binary searches for certain value, returns the required double 
	//	NEEDS TO BE CLEANED UP!!
   // Declare table columns:
	if( strcmp(fileName, "EG6ev2file.dat") )
	{
	}
	else if( strcmp(fileName, "EG6ECsampcorr.dat") )
	{
	}
	else {cout << "NOT A GOOD FILE" <<endl; return;}
	
	vector<vector<float>>	outVec;
	int initialized = 0;
	// Read the data table from file and put them into appropriate vectors:
   if (!initialized)
   {
		FILE *fin=fopen(fileName,"r");
     	if (!fin)
     	{
     	    fprintf(stderr," Missing Input File:  %s\n",filename);
     	    return;
     	}
		static const int bufSize = bufSize0;
		static const int nCols = nCols0;

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
				//cout << tmp[ii] << endl;
				row.push_back(tmp[ii]);
//				cout << "ROW INFOOROOOOOOO " << ii << "/t\t"<< row[ii] << endl;
			}
			outVec.push_back(row);
		}
   	initialized=1;
   	fclose(fin);
	}
//	cout << "OUT VEC STUFF: OUTVEC SIZE " << outVec.size() << endl;
//	cout << "OUT VEC : HOW BIG IS THE 0th VECTOR " << outVec[1].size() << endl;
	return outVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t EG5_test::getCorrectedSampFrac(Int_t irun, Int_t ievent, Int_t sector)
{
	// Gets corrected Sampling Fraction based on run and event number.
	//    sector == Sector (1-6)
	//    rf150  == Run# + File# / 150
	//
	//    Returns the sampling fraction parameterized by piecewise functions of the form:
	//    E0 + A0 * ( exp(-ALPHA1*(x-X0)) + exp(-ALPHA2*(x-tX0)) )
	//    where x = Run# + File#/150
	//    Negative return value is an error.
	//    Reads parameters from EG6ECsampcorr.dat

/*
	struct paramCols{
		paramCols(){
			static vector<float> param0;
			static vector<float> param1;
			static vector<float> param2;
			static vector<float> param3;
			static vector<float> param4;
			static vector<float> param5;
			static vector<float> param6;
			static vector<float> param7;
			static vector<float> param8;
		}

		const vector<float> getParam( vector<float> &param )
		{	return param;	}

	}
*/
   if (sector < 1 || sector > 6)
	{
		fprintf(stderr,"EG6ECsampcorr:  Invalid sector:  %d\n",sector);
      return -999;
   }
 	if ( ievent<0 )
 	{
   	fprintf(stderr,"EG6ev2file:  Invalid Event #:  %d\n",ievent);
    	return -999;
 	}


	static bool initialized=0;
	const double file = getInterpFileN(irun, ievent);
	if (file < 0) return -999;

	const double rf150 = irun + file/150.0;
	// check for invalid input:
	
	static bool initialized=0;
	if(!initSampCorr) 	// If everything is not initialized, get all files from .dat file
	{
		vector<vector<float>> fileParams = readDataTable("EG6ECsampcorr.dat", 1024, 2+5*6); 
		int vecSize0 = (int) fileParams.size();
		const int vecSize = vecSize0;

		//cout <<"== "<< vecSize0 << endl;
		//cout <<"==== "<<  vecSize << endl;

		Float_t tXLO[vecSize], tXHI[vecSize];  // <--- Run# + File#/150
   	Float_t tE0[6][vecSize],tX0[6][vecSize],tA0[6][vecSize],tALPH1[6][vecSize],tALPH2[6][vecSize];

		int nPars = 5;
		int jj = -999;
   	const int ss=sector-1;
		for( int ii = 0; ii < vecSize; ii++)
		{	
			// xLO and xHI are file run ranges: run \in range(xLO, xHI)
			tXLO[ii] = fileParams[ii][0];
			tXHI[ii] = fileParams[ii][1];

			// Get index of where rf150 is inbetween. No need for an obscure binary search
			if( rf150 <= tXHI[ii] ) 	// Since we are going from low to high, this will fail first
			if( rf150 >= tXLO[ii] )
			{	jj = ii;	}

			for(int sec = 0; sec < 6; sec++)
			{	
				tE0[sec][ii]	 = 	fileParams[ii][2+nPars*sec+0];
				tX0[sec][ii]	 = 	fileParams[ii][2+nPars*sec+1];
				tA0[sec][ii]	 = 	fileParams[ii][2+nPars*sec+2];
				tALPH1[sec][ii] = 	fileParams[ii][2+nPars*sec+3];
				tALPH2[sec][ii] = 	fileParams[ii][2+nPars*sec+4];
			}
		}
		initSampCorr = true;
	}
	//cout << rf150 << endl;
	//cout << tXLO[0] <<endl;
	//cout<< tXHI[vecSize-1] << endl;
	if (rf150 < tXLO[0] || rf150 >= tXHI[vecSize-1])
	{
  		fprintf(stderr,"EG6ECsampcorr:  Invalid rf150:  %f\n",rf150);
  	   return -999;
  	}
	
   // the parameterization:
   return tE0[ss][jj] + tA0[ss][jj] * (	exp(-tALPH1[ss][jj]*(rf150-tX0[ss][jj])) + exp(-tALPH2[ss][jj]*(rf150-tX0[ss][jj]))		) ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t EG5_test::getInterpFileN( Int_t irun, Int_t ievent)
{
	// Gets Interpolated File number
	if(!initFileInt)
	{
		vector<vector<float>> fileParams = readDataTable("EG6ev2file.dat", 256, 4); 
		int vecSize0 = (int) fileParams.size();
		const int vecSize = vecSize0;

		static int tRUN[vecSize], tFILE[vecSize], tEVMIN[vecSize], tEVMAX[vecSize];
		int jj = -999;

		for( int ii = 0 ; ii < vecSize; ii++)
		{
			tRUN[ii] 	= 	(int)	fileParams[ii][0];
			tFILE[ii] 	= 	(int) fileParams[ii][1];
   		tEVMIN[ii]	= 	(int) fileParams[ii][2];
   		tEVMAX[ii]	=	(int) fileParams[ii][3];
			
			// Grab index that corresponds to the right run and event number correction
			if( irun == tRUN[ii]		)
			if( ievent <= tEVMAX[ii]	)
			if( ievent >= tEVMIN[ii]	)
			{	jj = ii;	}
		}
		if ( irun < tRUN[0] || irun > tRUN[vecSize-1])
 		{
   		fprintf(stderr,"EG6ev2file:  Run Out of Range:  %d\n",irun);
   		 return -999;
 		}

		if (irun != tRUN[jj])
		{
			fprintf(stderr,"EG6ev2file:  Run Not Found:  %d\n",irun);
   		return -999;
		}
		initFileInt = true;
	}
   const double nEvents = tEVMAX[jj]-tEVMIN[jj]+1;
	return tFILE[jj] + (ievent-tEVMIN[jj])/nEvents;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t EG5_test::getCorrectedVzShift(Double_t theta, Double_t phi){
	// Corrects the z-vertex to account for the change in the beam center for different run numbers
	static const int nr=4;
	static const int runranges[nr]={61483,61580,61850,99999};
	static const double xx[nr]={0.155, 0.237, 0.27,  0.30}; // x beam position (cm)
	static const double yy[nr]={0.029,-0.040, 0.04, -0.04}; // y beam position (cm)
  	if (fabs(sin(theta)) < 1e-3) return -999;
	int therange=0;
	for (int ii=0; ii<nr; ii++) {
		if (runnb < runranges[ii]) {
			therange=ii;
			break;
		}
	}
	const double rr = sqrt(pow(xx[therange],2)+pow(yy[therange],2));
	const double phi0 = atan2(yy[therange],xx[therange])+TMath::Pi();
	return rr*cos(phi-phi0)/tan(theta);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 EG5_test::getICtoDCShift(TVector3 trackPos, TVector3 trackDir){
	if(trackPos(2) != 0)
	if(trackDir(2)  > 0){
		TVector3 shift = ((trackPos(2)-16)/trackDir(2))*trackDir;
		return shift;
	}
	return TVector3(0.,0.,0.);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 EG5_test::getUVWfromXYZ(TVector3 xyz){
 	// Converts x,y,z EC hit in CLAS coordinate system into u,v,w distances of the EC hit.
	using namespace TMath;
   // Parameters
   float ec_the = 0.4363323;
   float ylow = -182.974;
   float yhi = 189.956;
   float sinrho = 0.8901256;
   float cosrho = 0.455715;
	float tanrho = sinrho/cosrho;
   // Variables
	float ex = xyz[0];
   float wy = xyz[1];
   float zd = xyz[2];

   float phy = xyz.Phi()*180.0/TMath::Pi();
   if(phy <0.){ phy = phy + 360.;}
   phy = phy+30.;
   if(phy>360.){ phy = phy - 360.;}

   float ec_phy = ((Int_t) (phy/60.))*TMath::Pi()/3; //		2pi/6

   float rot[3][3] = {	{Cos(ec_the)*Cos(ec_phy), -Sin(ec_phy), 	Sin(ec_the)*Cos(ec_phy)	},
								{Cos(ec_the)*Sin(ec_phy),	Cos(ec_phy),	Sin(ec_the)*Sin(ec_phy)	},
								{				-Sin(ec_the),				0.,					Cos(ec_the)	}};
	
	float yi = ex*rot[0][0]+wy*rot[1][0]+zd*rot[2][0];
	float xi = ex*rot[0][1]+wy*rot[1][1]+zd*rot[2][1];
	float zi = ex*rot[0][2]+wy*rot[1][2]+zd*rot[2][2];
	zi = zi-510.32 ;
	
	float yu = (yi-ylow)/sinrho;
	float ve = (yhi-ylow)/tanrho - xi + (yhi-yi)/tanrho;
	float wu = ((yhi-ylow)/tanrho + xi + (yhi-yi)/tanrho)/2./cosrho;

	TVector3 uvw(yu,ve,wu);
	return uvw;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 EG5_test::getCCThetaPhi(TVector3 pos, TVector3 dir){
	TVector3 crossPoint = pos;
	TVector3 thetaPhi = TVector3(-999.,-999.,-999.);		// If not changed, CC Fid Cut fails.

	for(int ii = 0; ii<3; ii++)
		dir(ii) = dir(ii)/dir.Mag(); // Make sure the direction is normalized to unity

	if( dir.Mag() <= 0.000001) return thetaPhi;

	TVector3 ccPlane =  TVector3(- 0.0007840784063, 0., - 0.001681461571); //was static in Vlassov's code
	Double_t aa, bb;

	bb = 1.;
	for(int ii = 0; ii < 3; ii++)
		{ aa += ccPlane(ii)*dir(ii); bb += ccPlane(ii)*pos(ii); }
	if(fabs(bb) > 0.000001){
		if(fabs(aa) < 0.000001){ return thetaPhi; }
		Double_t tt = -bb/aa;
		crossPoint = tt*dir + pos; 			
	}


	thetaPhi = TVector3(	RadToDeg()*crossPoint.Theta(),	RadToDeg()*crossPoint.Phi(),	0.);		// 0th and 1st components are theta and phi respectively	
	
	return thetaPhi;																																														  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   Electron ID Cut Functions >
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t EG5_test::isElectron(Int_t partIndex){
	if( 	stat[partIndex] > 0				)				// Checks to see if it is a good run
	if( 	dc_stat[dc[partIndex]] > 0		)				// Checks to see if DC is good for this run
	if( 	q[partIndex] < 0					)	  			// Charge Cut
	if( 	nphe[cc[partIndex]] > 2.0 		) 				// CC Cut
	if( 	p[partIndex] > 0.7 				)				// Momentum Cut for delta electrons
	if( 	p[partIndex] < 6.0 				)				// 	(in GeV/c)
	if( 	isE_VertexCut(partIndex) 		)				// Make sure particle is not electrons coming from window
	if( 	isE_ECFiducialCut(partIndex) 	)				// Make sure particle is in the EC fiducial region 
	if( 	isE_ECEnergyCut(partIndex) 	)				// Make sure enough energy is deposited in the EC
	if( 	isE_CCFiducialCut(partIndex) 	)				// Make sure particle is in the CC fiducial region
	if( 	isE_DCFiducialCut(partIndex) 	)				// Make sure particle is in DC fiducial region
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EG5_test::Begin()
{
	Float_t *binformation2[] = {200,	0.25, 4, 	200,0.25,1.1,		200,0,70,	200,-180,55,		200,0,0.6,		200,0,0.5,		200,0,5,		200,0,1,		200,0,1,		200,0,1};
	char *hist2Names[] =  {"hbvp", "h2phivthe", "h2eovei", "h2etotpvp", "h2etotveieo"};	
	char *hist2Titles[] = {"Beta vs Momentum", "Phi vs Theta", "Eo vs Ei", "Etot/p vs p", "Etot vs Ei+Eo"};
	h2length = sizeof(hist2Names)/sizeof(hist2Names[0]); 
	
	Float_t *binformation1[] = {200,	0.25,	1.1,		200, 	0.25, 4,		200,	-100, -20,		200, 	-180, 180,		200,	-180,	180};
	char *hist1Names[] = {"h1b", "h1p", "h1vz", "h1theta", "h1phi"};
	char *hist1Titles[] = {"Beta", "Momentum", "Vertex Z", "Theta", "Phi"};
	h1length = sizeof(hist1Names)/sizeof(hist1Names[0]); 


	for(int ii = 0; ii < (int) h2length; ii++){ 
	hCuts2[ii] = new TH2D(hist2Names[ii], hist2Titles[ii], (int) binformation2[6*ii+0], (Float_t) binformation2[6*ii+1], (Float_t) binformation2[6*ii+2], (int)binformation2[6*ii+3],(Float_t)binformation2[6*ii+4], (Float_t)binformation2[6*ii+5]);
	}; 

	for(int ii = 0; ii < (int) h1length; ii++){ 
	  hCuts1[ii] = new TH1D(hist1Names[ii], hist1Titles[ii], (int) binformation1[3*ii+0], (Float_t) binformation1[3*ii+1], (Float_t) binformation1[3*ii+2]);
	};


}


void EG5_test::Loop()
{
	using namespace TMath;
	if (fChain == 0) return;
	//cout << ( fChain->GetCurrentFile()->GetName()).c_str() << endl;
	const char *currentFN = ((TChain*)(EG5_test::fChain))->GetFile()->GetName();
	cout << currentFN<< endl;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry < nentries && jentry < 10000000000000000000; jentry++) {
		return;
		int ient = 760000;
		if (jentry % 10000 == 0) {cout << "\t Event: \t" << jentry/1000 << " k \t  / \t " << ient/1000 << "\t k     \t=\t" << (Float_t) jentry/(ient+1)*100 << "\t %" << endl;}

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   
		nbytes += nb;
		//	if (Cut(ientry) < 0) continue;
		if( goodDetectors() ){
			for( int ipart = 0; ipart < gpart; ipart++){
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
	TFile outfile("hShort.root");//,"overwrite");
	//Write the histogram in the file
	Hlist1.Add(hCuts1);
	Hlist2.Add(hCuts2);
	Int_t maxHistLength = (Int_t) TMath::Max((Int_t) h2length, (Int_t) h1length);
	for( Int_t ii = 0; ii < maxHistLength; ii ++){
		if(ii < h1length){hCuts1[ii]->Write();};
		if(ii < h2length){hCuts2[ii]->Write();};
	}
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
