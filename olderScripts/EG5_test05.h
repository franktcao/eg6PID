//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 15 14:45:23 2016 by ROOT version 6.05/02
// from TChain ch/EG6_test
//////////////////////////////////////////////////////////

#ifndef EG5_test_h
#define EG5_test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
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
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
// Header file for the classes stored in the TTree if any.

class EG5_test {
   public :
   TTree         *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t          fCurrent; //!current Tree number in a TChain

   TH2D          *hbvp;
   TH2D          *h2phivthe;
   TH2D          *h2eovei;
   TH2D          *h2etotpvp;
	TH2D			  *h2etotveieo;

   TH1D			  *h1b;
   TH1D			  *h1p;
   TH1D			  *h1vz;
   TH1D			  *h1theta;
   TH1D			  *h1phi;

	int	h2length;
	int	h1length;
	
	int	eIndex;
	Double_t eTrigTime;

	TVector3 eV;
	TVector3 pV;
	TVector3 shiftICtoDC; 
	Double_t eVz;
	


	vector<TH2D> 	hCuts2;
	vector<TH1D> 	hCuts1;
	Double_t			myTheta;
	Double_t			myPhi;
	Double_t			etote;
	Int_t				sector;
	Int_t				nRowsSampFrac;
	Int_t				nRowsFileIntrp;
	//int				*tRun;
	//int				*tEvPerFile;
	Float_t			*tXLO;
	Float_t			*tXHI;
	Float_t			**tE0;
	Float_t			**tX0;
	Float_t			**tA0;
	Float_t			**tALPH1;
	Float_t			**tALPH2;
	char				*curFile;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UChar_t         npart;
   UInt_t          runnb;
   UInt_t          evntid;
   UChar_t         evstat;
   Char_t          evntype;
   Int_t           evntclas;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time;
   Int_t           l2bit;
   Int_t           l3bit;
   UChar_t         helicity;
   Int_t           hlsc;
   Int_t           intt;
   Int_t           helicity_cor;
   Int_t           gpart;
   Int_t           id[40];   //[gpart]
   Int_t           stat[40];   //[gpart]
   Int_t           dc[40];   //[gpart]
   Int_t           cc[40];   //[gpart]
   Int_t           sc[40];   //[gpart]
   Int_t           ec[40];   //[gpart]
   Int_t           lec[40];   //[gpart]
   Int_t           st[40];   //[gpart]
   Float_t         p[40];   //[gpart]
   Float_t         m[40];   //[gpart]
   Int_t           q[40];   //[gpart]
   Float_t         b[40];   //[gpart]
   Float_t         cx[40];   //[gpart]
   Float_t         cy[40];   //[gpart]
   Float_t         cz[40];   //[gpart]
   Float_t         vx[40];   //[gpart]
   Float_t         vy[40];   //[gpart]
   Float_t         vz[40];   //[gpart]
   Int_t           dc_part;
   Int_t           dc_sect[40];   //[dc_part]
   Int_t           dc_trk[40];   //[dc_part]
   Int_t           dc_stat[40];   //[dc_part]
   Int_t           tb_st[40];   //[dc_part]
   Float_t         dc_xsc[40];   //[dc_part]
   Float_t         dc_ysc[40];   //[dc_part]
   Float_t         dc_zsc[40];   //[dc_part]
   Float_t         dc_cxsc[40];   //[dc_part]
   Float_t         dc_cysc[40];   //[dc_part]
   Float_t         dc_czsc[40];   //[dc_part]
   Float_t         dc_vx[40];   //[dc_part]
   Float_t         dc_vy[40];   //[dc_part]
   Float_t         dc_vz[40];   //[dc_part]
   Float_t         dc_vr[40];   //[dc_part]
   Float_t         tl1_cx[40];   //[dc_part]
   Float_t         tl1_cy[40];   //[dc_part]
   Float_t         tl1_cz[40];   //[dc_part]
   Float_t         tl1_x[40];   //[dc_part]
   Float_t         tl1_y[40];   //[dc_part]
   Float_t         tl1_z[40];   //[dc_part]
   Float_t         tl1_r[40];   //[dc_part]
   Float_t         dc_c2[40];   //[dc_part]
   Int_t           ec_part;
   Int_t           ec_stat[40];   //[ec_part]
   Int_t           ec_sect[40];   //[ec_part]
   Int_t           ec_whol[40];   //[ec_part]
   Int_t           ec_inst[40];   //[ec_part]
   Int_t           ec_oust[40];   //[ec_part]
   Float_t         etot[40];   //[ec_part]
   Float_t         ec_ei[40];   //[ec_part]
   Float_t         ec_eo[40];   //[ec_part]
   Float_t         ec_t[40];   //[ec_part]
   Float_t         ec_r[40];   //[ec_part]
   Float_t         ech_x[40];   //[ec_part]
   Float_t         ech_y[40];   //[ec_part]
   Float_t         ech_z[40];   //[ec_part]
   Float_t         ec_m2[40];   //[ec_part]
   Float_t         ec_m3[40];   //[ec_part]
   Float_t         ec_m4[40];   //[ec_part]
   Float_t         ec_c2[40];   //[ec_part]
   Int_t           sc_part;
   Int_t           sc_sect[40];   //[sc_part]
   Int_t           sc_hit[40];   //[sc_part]
   Int_t           sc_pd[40];   //[sc_part]
   Int_t           sc_stat[40];   //[sc_part]
   Float_t         edep[40];   //[sc_part]
   Float_t         sc_t[40];   //[sc_part]
   Float_t         sc_r[40];   //[sc_part]
   Float_t         sc_c2[40];   //[sc_part]
   Int_t           cc_part;
   Int_t           cc_sect[40];   //[cc_part]
   Int_t           cc_hit[40];   //[cc_part]
   Int_t           cc_segm[40];   //[cc_part]
   Int_t           nphe[40];   //[cc_part]
   Float_t         cc_t[40];   //[cc_part]
   Float_t         cc_r[40];   //[cc_part]
   Float_t         cc_c2[40];   //[cc_part]
   Int_t           ic_part;
   Float_t         et[200];   //[ic_part]
   Float_t         egl[200];   //[ic_part]
   Float_t         time[200];   //[ic_part]
   Float_t         time_next[200];   //[ic_part]
   Float_t         ich_x[200];   //[ic_part]
   Float_t         ich_y[200];   //[ic_part]
   Float_t         ich_z[200];   //[ic_part]
   Float_t         ich_xgl[200];   //[ic_part]
   Float_t         ich_ygl[200];   //[ic_part]
   Float_t         ich_xwid[200];   //[ic_part]
   Float_t         ich_ywid[200];   //[ic_part]
   Float_t         ich_xm3[200];   //[ic_part]
   Float_t         ich_ym3[200];   //[ic_part]
   Int_t           ic_stat[200];   //[ic_part]
   Int_t           gcpart;
   Int_t           pid[800];   //[gcpart]
   Float_t         x[800];   //[gcpart]
   Float_t         y[800];   //[gcpart]
   Float_t         z[800];   //[gcpart]
   Float_t         dedx[800];   //[gcpart]
   Float_t         px[800];   //[gcpart]
   Float_t         py[800];   //[gcpart]
   Float_t         pz[800];   //[gcpart]
   Float_t         p_tot[800];   //[gcpart]
   Float_t         x2[800];   //[gcpart]
   Float_t         theta[800];   //[gcpart]
   Float_t         charge[800];   //[gcpart]
   Float_t         dca[800];   //[gcpart]
   Int_t           index[800];   //[gcpart]
   Float_t         phi[800];   //[gcpart]
   Float_t         vtl[800];   //[gcpart]
   Float_t         sdist[800];   //[gcpart]
   Float_t         edist[800];   //[gcpart]
   Int_t           npts[800];   //[gcpart]
   Float_t         r_0[800];   //[gcpart]
   Int_t           fiterr[800];   //[gcpart]
   Int_t           tothits[800];   //[gcpart]
   Int_t           npd_track[800];   //[gcpart]
   Int_t           npd_event[800];   //[gcpart]
   Int_t           bonus_bits[800];   //[gcpart]
   Float_t         q_tot[800];   //[gcpart]
   Float_t         x_start[800];   //[gcpart]
   Float_t         y_start[800];   //[gcpart]
   Float_t         z_start[800];   //[gcpart]
   Float_t         x_end[800];   //[gcpart]
   Float_t         y_end[800];   //[gcpart]
   Float_t         z_end[800];   //[gcpart]
   Int_t           rtpc_npart;
   Int_t           rtpc_id1[50];   //[rtpc_npart]
   Int_t           rtpc_id2[50];   //[rtpc_npart]
   Int_t           rtpc_id3[50];   //[rtpc_npart]
   Int_t           rtpc_id4[50];   //[rtpc_npart]
   Int_t           rtpc_id5[50];   //[rtpc_npart]
   Float_t         rtpc_p1[50];   //[rtpc_npart]
   Float_t         rtpc_p2[50];   //[rtpc_npart]
   Float_t         rtpc_p3[50];   //[rtpc_npart]
   Float_t         rtpc_p4[50];   //[rtpc_npart]
   Float_t         rtpc_p5[50];   //[rtpc_npart]
   Float_t         rtpc_poverq[50];   //[rtpc_npart]
   Float_t         rtpc_dedx[50];   //[rtpc_npart]
   Float_t         rtpc_dedx2[50];   //[rtpc_npart]
   Float_t         rtpc_dedxa[50];   //[rtpc_npart]
   Float_t         rtpc_dedxl[50];   //[rtpc_npart]
   Float_t         rtpc_dedxal[50];   //[rtpc_npart]
   Float_t         rtpc_rxy[50];   //[rtpc_npart]
   Float_t         rtpc_slope[50];   //[rtpc_npart]
   Float_t         rtpc_chisq[50];   //[rtpc_npart]
   Float_t         rtpc_dedxs[50];   //[rtpc_npart]
   Float_t         rtpc_theta[50];   //[rtpc_npart]
   Float_t         rtpc_phi[50];   //[rtpc_npart]
   Float_t         rtpc_vz[50];   //[rtpc_npart]
   Float_t         rtpc_bad[50];   //[rtpc_npart]
   Int_t           rtpc_gcpb[50];   //[rtpc_npart]
   Int_t           icpart;
   Float_t         etc[30];   //[icpart]
   Float_t         ecc[30];   //[icpart]
   Float_t         tc[30];   //[icpart]
   Float_t         tn[30];   //[icpart]
   Float_t         xc[30];   //[icpart]
   Float_t         yc[30];   //[icpart]
   Float_t         zc[30];   //[icpart]
   Float_t         m2c[30];   //[icpart]
   Float_t         m3c[30];   //[icpart]
   Int_t           statc[30];   //[icpart]
   Int_t           shh_part;
   Int_t           shh_id[72];   //[shh_part]
   Float_t         shh_x[72];   //[shh_part]
   Float_t         shh_y[72];   //[shh_part]
   Float_t         shh_z[72];   //[shh_part]
   Float_t         shh_nphe[72];   //[shh_part]
   Float_t         shh_time[72];   //[shh_part]
   Int_t           shh_stat[72];   //[shh_part]
   Int_t           shpart;
   Int_t           shid[54];   //[shpart]
   Float_t         shx[54];   //[shpart]
   Float_t         shy[54];   //[shpart]
   Float_t         shz[54];   //[shpart]
   Float_t         shnphe[54];   //[shpart]
   Float_t         shtime[54];   //[shpart]
   Int_t           shstat[54];   //[shpart]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_runnb;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time;   //!
   TBranch        *b_l2bit;   //!
   TBranch        *b_l3bit;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_hlsc;   //!
   TBranch        *b_intt;   //!
   TBranch        *b_helicity_cor;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_lec;   //!
   TBranch        *b_st;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_tb_st;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_vx;   //!
   TBranch        *b_dc_vy;   //!
   TBranch        *b_dc_vz;   //!
   TBranch        *b_dc_vr;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_tl1_cz;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl1_r;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_ic_part;   //!
   TBranch        *b_et;   //!
   TBranch        *b_egl;   //!
   TBranch        *b_time;   //!
   TBranch        *b_time_next;   //!
   TBranch        *b_ich_x;   //!
   TBranch        *b_ich_y;   //!
   TBranch        *b_ich_z;   //!
   TBranch        *b_ich_xgl;   //!
   TBranch        *b_ich_ygl;   //!
   TBranch        *b_ich_xwid;   //!
   TBranch        *b_ich_ywid;   //!
   TBranch        *b_ich_xm3;   //!
   TBranch        *b_ich_ym3;   //!
   TBranch        *b_ic_stat;   //!
   TBranch        *b_gcpart;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_p_tot;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_dca;   //!
   TBranch        *b_index;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_vtl;   //!
   TBranch        *b_sdist;   //!
   TBranch        *b_edist;   //!
   TBranch        *b_npts;   //!
   TBranch        *b_r_0;   //!
   TBranch        *b_fiterr;   //!
   TBranch        *b_tothits;   //!
   TBranch        *b_npd_track;   //!
   TBranch        *b_npd_event;   //!
   TBranch        *b_bonus_bits;   //!
   TBranch        *b_q_tot;   //!
   TBranch        *b_x_start;   //!
   TBranch        *b_y_start;   //!
   TBranch        *b_z_start;   //!
   TBranch        *b_x_end;   //!
   TBranch        *b_y_end;   //!
   TBranch        *b_z_end;   //!
   TBranch        *b_rtpc_npart;   //!
   TBranch        *b_rtpc_id1;   //!
   TBranch        *b_rtpc_id2;   //!
   TBranch        *b_rtpc_id3;   //!
   TBranch        *b_rtpc_id4;   //!
   TBranch        *b_rtpc_id5;   //!
   TBranch        *b_rtpc_p1;   //!
   TBranch        *b_rtpc_p2;   //!
   TBranch        *b_rtpc_p3;   //!
   TBranch        *b_rtpc_p4;   //!
   TBranch        *b_rtpc_p5;   //!
   TBranch        *b_rtpc_poverq;   //!
   TBranch        *b_rtpc_dedx;   //!
   TBranch        *b_rtpc_dedx2;   //!
   TBranch        *b_rtpc_dedxa;   //!
   TBranch        *b_rtpc_dedxl;   //!
   TBranch        *b_rtpc_dedxal;   //!
   TBranch        *b_rtpc_rxy;   //!
   TBranch        *b_rtpc_slope;   //!
   TBranch        *b_rtpc_chisq;   //!
   TBranch        *b_rtpc_dedxs;   //!
   TBranch        *b_rtpc_theta;   //!
   TBranch        *b_rtpc_phi;   //!
   TBranch        *b_rtpc_vz;   //!
   TBranch        *b_rtpc_bad;   //!
   TBranch        *b_rtpc_gcpb;   //!
   TBranch        *b_icpart;   //!
   TBranch        *b_etc;   //!
   TBranch        *b_ecc;   //!
   TBranch        *b_tc;   //!
   TBranch        *b_tn;   //!
   TBranch        *b_xc;   //!
   TBranch        *b_yc;   //!
   TBranch        *b_zc;   //!
   TBranch        *b_m2c;   //!
   TBranch        *b_m3c;   //!
   TBranch        *b_statc;   //!
   TBranch        *b_shh_part;   //!
   TBranch        *b_shh_id;   //!
   TBranch        *b_shh_x;   //!
   TBranch        *b_shh_y;   //!
   TBranch        *b_shh_z;   //!
   TBranch        *b_shh_nphe;   //!
   TBranch        *b_shh_time;   //!
   TBranch        *b_shh_stat;   //!
   TBranch        *b_shpart;   //!
   TBranch        *b_shid;   //!
   TBranch        *b_shx;   //!
   TBranch        *b_shy;   //!
   TBranch        *b_shz;   //!
   TBranch        *b_shnphe;   //!
   TBranch        *b_shtime;   //!
   TBranch        *b_shstat;   //!

   EG5_test(TTree *tree=0);
   virtual ~EG5_test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Begin();
   virtual void     Loop();
   virtual void     printNumEvents();
   virtual void     drawHist();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
	//virtual TVector3 getCCThetaPhi(TVector3 pos, TVector3 dir);
	virtual void	  writeFiles();
	//virtual Bool_t   isEVertexCut(Int_t ipart);
	//virtual Bool_t	  isElectron(Int_t ipart);
	Bool_t goodDetectors();
	//vector<vector<float>> readDataTable(TString fileName, const int bufSize,  const int nCols);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<vector<float>> readDataTable(TString fileName0, const int bufSize,  const int nCols)
{	// Reads data from table, puts them into vectors, checks to see if vectors are not nonsensical, binary searches for certain value, returns the required double 
   // Declare table columns:
	vector<vector<float>>	outVec;
	
	const char* fileName;
	fileName = fileName0.Data();

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
						 row.push_back(tmp[ii]);
		}
		outVec.push_back(row);
	}
	fclose(fin);
	return outVec;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initHists(){
	// Initiate Histograms that will be needed
	// 	Binfo2: Bin Information looks like binning in x, x min, x max, binning in y, y min, y max.
	vector<vector<Float_t>> binfo2; 
		binfo2.push_back({200.,	0.25, 4. , 		200., 0.25, 1.1});
		binfo2.push_back({200.,	0.	 , 70.,		200.,-180., 55.});		
		binfo2.push_back({200.,	0.	 ,	0.6,		200., 0.,	0.5});		
		binfo2.push_back({200.,	0.	 ,	5. ,		200., 0.,	1.});		
		binfo2.push_back({200.,	0.	 ,	1. ,		200., 0.,	1.});
	vector<string> hist2Names; 
		hist2Names.push_back("hbvp");
		hist2Names.push_back("h2phivthe");
		hist2Names.push_back("h2eovei");
		hist2Names.push_back("h2etotpvp");
		hist2Names.push_back("h2etotveieo");	
	vector<string> hist2Titles;
		hist2Titles.push_back("Beta vs Momentum"); 
		hist2Titles.push_back("Phi vs Theta"); 
		hist2Titles.push_back("Eo vs Ei"); 
		hist2Titles.push_back("Etot/p vs p"); 
		hist2Titles.push_back("Etot vs Ei+Eo");
	h2length = hist2Names.size(); 

	// Binfo2: Bin Information looks like binning in x, x min, x max, binning in y, y min, y max.
	vector<vector<Float_t>> binfo1; 
		binfo1.push_back({200.,	0.25,	1.1});	
		binfo1.push_back({200., 	0.25, 4.});		
		binfo1.push_back({200.,	-100., -20.});		
		binfo1.push_back({200., 	-180., 180.});		
		binfo1.push_back({200.,	-180.,	180.});
		binfo1.push_back({40.,	0.,	40.});

	vector<string> hist1Names; 
		hist1Names.push_back("h1b"); 
		hist1Names.push_back("h1p"); 
		hist1Names.push_back("h1vz"); 
		hist1Names.push_back("h1theta"); 
		hist1Names.push_back("h1phi");
		hist1Names.push_back("fail");

	vector<string> hist1Titles;
		hist1Titles.push_back("Beta"); 
		hist1Titles.push_back("Momentum"); 
		hist1Titles.push_back("Vertex Z"); 
		hist1Titles.push_back("Theta"); 
		hist1Titles.push_back("Phi");
		hist1Titles.push_back("Failed Tests");

	h1length = hist1Names.size(); 

	static const int h2length0 = h2length;
	static const int h1length0 = h1length;

	for(int ii = 0; ii < h2length; ii++){ 
		TH2D curHist = TH2D(hist2Names[ii].c_str(), hist2Titles[ii].c_str(), (int) binfo2[ii][0], binfo2[ii][1], binfo2[ii][2], (int) binfo2[ii][3], binfo2[ii][4], binfo2[ii][5]);
		hCuts2.push_back(curHist);	
		}


	for(int ii = 0; ii < (int) h1length; ii++){ 
		TH1D curHist = TH1D(hist1Names[ii].c_str(), hist1Titles[ii].c_str(), (int) binfo1[ii][0], binfo1[ii][1], binfo1[ii][2]);
		hCuts1.push_back(curHist);	
		} 
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void initFileIntrp(){
	// INITIALIZE EG6 File Number Interpolation
	vector<vector<float>> fileParams = readDataTable("EG6ev2file.dat", 256, 4); 
	vector<vector<int>> redFileParams;				// Reduced File Parameters. I think this is all we need.
	int tempRunNum = -999;

	int nRows = (int) fileParams.size();
	for( int ii = 0 ; ii < nRows; ii++)
	{
		vector<int> row;
		if(tempRunNum != (int) fileParams[ii][0]) 
		{
			row.push_back((int) fileParams[ii][0]);
			row.push_back((int) fileParams[ii][3]);
			redFileParams.push_back(row); 
			tempRunNum = (int) fileParams[ii][0];
		}
	}
	
	// Redfine nRows, because I can.
	nRows = (int) redFileParams.size();
	nRowsFileIntrp = nRows;
	//cout << "number of rows from reduced parameters" <<  nRows << endl;
	//cout << redFileParams.capacity() << endl;
	const int nRows0  = nRows;

	//const int tRUN[vecSize], tEvPerFile[vecSize];
	tRun = new int[nRows0];
	tEvPerFile = new int[nRows0];
	
	for( int ii = 0; ii < nRows0; ii++)
	{
		tRun[ii] = (int) redFileParams[ii][0];
		tEvPerFile[ii] = (int) redFileParams[ii][1];
	}
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initSampCorr(){
	// INITIALIZE EG6EC sampling fraction corrections
	vector<vector<float>> fileParams00 = readDataTable("EG6ECsampcorr.dat", 1024, 2+5*6); 
	int nRows = (int) fileParams00.size();
	const int nRows00 = nRows;
	nRowsSampFrac = nRows;

	tXLO = new Float_t[nRows00];
	tXHI = new Float_t[nRows00];
	tE0 = new Float_t*[nRows00];
	for( int ii = 0 ; ii < nRows00; ii++)
	{	tE0[ii] = new Float_t[6];	}
	
	tX0 = new Float_t*[nRows00];
	for( int ii = 0 ; ii < nRows00; ii++)
	{	tX0[ii] = new Float_t[6];	}
	
	tA0 = new Float_t*[nRows00];
	for( int ii = 0 ; ii < nRows00; ii++)
	{	tA0[ii] = new Float_t[6];	}
	
	tALPH1 = new Float_t*[nRows00];
	for( int ii = 0 ; ii < nRows00; ii++)
	{	tALPH1[ii] = new Float_t[6];	}
	
	tALPH2 = new Float_t*[nRows00];
	for( int ii = 0 ; ii < nRows00; ii++)
	{	tALPH2[ii] = new Float_t[6];	}
	
	int nPars = 5;
	for( int ii = 0; ii < nRows00; ii++)
	{	
		tXLO[ii] = (Float_t) fileParams00[ii][0];
		tXHI[ii] = (Float_t) fileParams00[ii][1];
		for(int sec = 0; sec < 6; sec++)
		{	
			tE0[ii][sec]	 = (Float_t) 	fileParams00[ii][2+nPars*sec+0];
			tX0[ii][sec]	 = (Float_t)	fileParams00[ii][2+nPars*sec+1];
			tA0[ii][sec] 	 = (Float_t)	fileParams00[ii][2+nPars*sec+2];
			tALPH1[ii][sec] = (Float_t)	fileParams00[ii][2+nPars*sec+3];
			tALPH2[ii][sec] = (Float_t)	fileParams00[ii][2+nPars*sec+4];
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t getFileNum(char* fullFileName0){
	// Splits full file name into tokens. 
	// 	Splits string up by "/" so each token is a directory or filename
	TString fullFileName = (TString) fullFileName0;
	vector<char*> fullFileNameV;
	char *tokens = strtok( (char*) fullFileName.Data(),  " / ");
	while (tokens) {
		fullFileNameV.push_back(tokens);
		tokens = strtok(NULL, "/");
	}
	//		Splits resulting string by "." and then "_" to get at the file number	
	char* file = fullFileNameV[(int) fullFileNameV.size() - 1];
	char* fileName = strtok( (char*) file ,  ".");
	tokens = strtok( (char*) fileName,  "_");
	vector<char*> fileNameV;
	while (tokens) {
		fileNameV.push_back(tokens);
		tokens = strtok(NULL, "_");
	}
	return (Int_t) *fileNameV[3];
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t getCorrectedSampFrac(Int_t irun, Int_t ievent, Int_t sector, char* curFile)
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
	double file = getFileNum(curFile);
	//InterpFileN(irun, ievent);

	double rf150 = irun + file/150.0;
	//cout << "FILE NUMBER IS     " << file << endl;
	//cout << "rf150 IS     " << rf150 << endl;
	// check for invalid input:
	
	int ilast = nRowsSampFrac-1;
	// cout << ilast << endl;
	if (rf150 < tXLO[0] || rf150 >= tXHI[ilast])
	{
  		fprintf(stderr,"EG6ECsampcorr:  Invalid rf150:  %f\n",rf150);
  	   return -999;
  	}
	 // binary search:
	int imax = ilast;
   int imin = 0;
   int jj = -1;
   
	while (imax >= imin)
   {    
   	jj = (imax+imin)/2;
  	   if      (rf150 >= tXHI[jj]) imin=jj+1;
  		else if (rf150 <  tXLO[jj]) imax=jj-1;
		else break;			//{ cout << "Break with " << jj << endl; break;}
   }    
	if( jj > ilast || jj < 0 )
	{	return -999;	}

	int ss = sector - 1;
   // the parameterization:
   return tE0[jj][ss] + tA0[jj][ss] * (	exp(-tALPH1[jj][ss]*(rf150-tX0[jj][ss])) + exp(-tALPH2[jj][ss]*(rf150-tX0[jj][ss]))		) ;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* OLD, NOT USED
Double_t getInterpFileN( Int_t irun, Int_t ievent)
{
	// Gets Interpolated File number
    // binary search:
	Int_t iMaxx = nRowsFileIntrp;//sizeof(tRun) / sizeof(*tRun) ;
	//cout << "trun (2) is " << tRun[3] << endl;
	//cout << "iMAXXXXX " << iMaxx << endl;
	Int_t iMinn = 0;
	Int_t jjj = -1;

	while (iMaxx >= iMinn)
	{    
		jjj = (iMaxx+iMinn)/2;
		if      (irun >= tRun[jjj]) iMinn=jjj+1;
		else if (irun <  tRun[jjj]) iMaxx=jjj-1;
		else break;
	}
	if( jjj < 0 || jjj > nRowsFileIntrp)
	{	cout << "JJJ NO GOOD" << endl; return	-999;}
	//cout << "jjj" <<  jjj << endl;
	//cout << "==irun" << irun << endl;
	//cout << "===tRun " << tRun[jjj] << endl;
	//cout << "====ievent" << ievent << endl;
	//cout << "===== Events per file" << tEvPerFile[jjj] << endl;
	//cout << " JJJ IS : ====== " << jjj << endl;
	double  nEvents = 1 + tEvPerFile[jjj];//1 + tRun[jjj]-tRun[jjj-1];
	return ievent/tEvPerFile[jjj] + (ievent % tEvPerFile[jjj])/nEvents;    		// Integer division of ievent/tEvPerFile is my estimate for file number
	//return irun + ievent/((float) tEvPerFile[jjj]);//floor(ievent/((float) tEvPerFile[jjj]));
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t getCorrectedVzShift(Double_t theta, Double_t phi){
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

TVector3 getICtoDCShift(TVector3 trackPos, TVector3 trackDir){
	if(trackPos(2) != 0)
	if(trackDir(2)  > 0){
		TVector3 shift = ((trackPos(2)-16)/trackDir(2))*trackDir;
		return shift;
	}
	return TVector3(0.,0.,0.);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 getUVWfromXYZ(TVector3 xyz){
 	// Converts x,y,z EC hit in CLAS coordinate system into u,v,w distances of the EC hit.
	using namespace TMath;
   // Parameters
   Double_t ec_the = 0.4363323;
   Double_t ylow = -182.974;
   Double_t yhi = 189.956;
   Double_t sinrho = 0.8901256;
   Double_t cosrho = 0.455715;
	Double_t tanrho = sinrho/cosrho;
   // Variables
	Double_t ex = xyz[0];
   Double_t wy = xyz[1];
   Double_t zd = xyz[2];

   Double_t phy = xyz.Phi()*180.0/TMath::Pi();
   if(phy <0.){ phy = phy + 360.;}
   phy = phy+30.;
   if(phy>360.){ phy = phy - 360.;}

   Double_t ec_phy = ((Int_t) (phy/60.))*TMath::Pi()/3; //		2pi/6

   Double_t rot[3][3] = {	{Cos(ec_the)*Cos(ec_phy), -Sin(ec_phy), 	Sin(ec_the)*Cos(ec_phy)	},
								{Cos(ec_the)*Sin(ec_phy),	Cos(ec_phy),	Sin(ec_the)*Sin(ec_phy)	},
								{				-Sin(ec_the),				0.,					Cos(ec_the)	}};
	
	Double_t yi = ex*rot[0][0]+wy*rot[1][0]+zd*rot[2][0];
	Double_t xi = ex*rot[0][1]+wy*rot[1][1]+zd*rot[2][1];
	Double_t zi = ex*rot[0][2]+wy*rot[1][2]+zd*rot[2][2];
	zi = zi-510.32 ;
	
	Double_t yu = (yi-ylow)/sinrho;
	Double_t ve = (yhi-ylow)/tanrho - xi + (yhi-yi)/tanrho;
	Double_t wu = ((yhi-ylow)/tanrho + xi + (yhi-yi)/tanrho)/2./cosrho;

	TVector3 uvw(yu,ve,wu);
	return uvw;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 getCCThetaPhi(TVector3 pos, TVector3 dir){
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
	if( fabs(bb) > 0.000001 ){
		if(fabs(aa) < 0.000001){ return thetaPhi; }
		Double_t tt = -bb/aa;
		crossPoint = tt*dir + pos; 			
	}
	
	using namespace TMath;
	thetaPhi = TVector3(	RadToDeg()*crossPoint.Theta(),	RadToDeg()*crossPoint.Phi(),	0.);		// 0th and 1st components are theta and phi respectively	
	
	return thetaPhi;																																														  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isVertexMatch(Float_t partVz){
//	 return true;
	// Makes sure the particle comes from the same vertex as the electron so we know they are from the same event 
	Double_t dVz = eVz - partVz;
	if( Abs(dVz) < 2.0)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t getDeltaT(Int_t ipart){
//	 return true;
	// Calculates the time difference between given particle and electron time 
	Double_t measTime = sc_t[sc[ipart]];
	Double_t relTime = measTime - eTrigTime;

	Double_t pp = p[ipart];
	Double_t mm = m[ipart];
	Double_t calcTime = sc_r[ipart] /  getCalcBeta(pp, mm) / 30.;		//	Calculated time is (sc pathlength)/(sc velocity)
	
	Double_t deltaT = relTime - calcTime;
	return deltaT;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t getCalcBeta(Double_t pp, Double_t mass){
	return pp/Sqrt(pow(pp,2) + pow(mass,2));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   Electron ID Cut Functions >
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OPTIMIZING
/*
Bool_t isElectron(Int_t ipart){
	if( 	!(stat[ipart] > 0)					)		hCuts1[5].Fill(0);		// Checks to see if it is a good run
	if( 	!(dc_stat[dc[ipart]] > 0)		)		hCuts1[5].Fill(1);		// Checks to see if DC is good for this run
	if( 	!(q[ipart] < 0)						)	  	hCuts1[5].Fill(2);		// Charge Cut
	if( 	!(nphe[cc[ipart]] > 2.0) 		) 		hCuts1[5].Fill(3);		// CC Cut
	if( 	!(p[ipart] > 0.7) 					)		hCuts1[5].Fill(4);		// Momentum Cut for delta electrons
	if( 	!(p[ipart] < 6.0) 					)		hCuts1[5].Fill(5);		// 	(in GeV/c)
	if( 	!isE_VertexCut(ipart) 			)		hCuts1[5].Fill(6);		// Make sure particle is not electrons coming from window
	if( 	!isE_ECFiducialCut(ipart) 	)		hCuts1[5].Fill(7);		// Make sure particle is in the EC fiducial region 
	if( 	!isE_ECEnergyCut(ipart) 		)		hCuts1[5].Fill(8);		// Make sure enough energy is deposited in the EC
	if( 	!isE_CCFiducialCut(ipart) 	)		hCuts1[5].Fill(9);		// Make sure particle is in the CC fiducial region
	if( 	!isE_DCFiducialCut(ipart) 	)		hCuts1[5].Fill(10);		// Make sure particle is in DC fiducial region
	if(0){
		return true;
	}
	return false;
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isElectron(Int_t ipart){
	if( 	q[ipart] < 0					)	  			// Charge Cut
	if( 	stat[ipart] > 0				)				// Checks to see if it is a good run
	if( 	dc_stat[dc[ipart]] > 0		)				// Checks to see if DC is good for this run
	if( 	isE_ECEnergyCut(ipart) 		)				// Make sure enough energy is deposited in the EC
	if( 	isE_VertexCut(ipart) 		)				// Make sure particle is not electrons coming from window
	if( 	isDCFiducialCut(ipart)	 	)				// Make sure particle is in DC fiducial region
	if( 	p[ipart] > 0.7 				)				// Momentum Cut for delta electrons
	if( 	isE_CCFiducialCut(ipart) 	)				// Make sure particle is in the CC fiducial region
	if( 	p[ipart] < 6.0 				)				// 	(in GeV/c)
	if( 	nphe[cc[ipart]] > 2.0 		) 				// CC Cut
	if( 	isE_ECFiducialCut(ipart) 	)				// Make sure particle is in the EC fiducial region 
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isProton(Int_t ipart){
	if( 	q[ipart] > 0					)	  			// Charge Cut
	if( 	b[ipart] != 0					)	  			// Nonzero velocity 
	if( 	stat[ipart] > 0				)				// Checks to see if it is a good run
	if( 	dc_stat[dc[ipart]] > 0		)				// Checks to see if DC is good for this run
	if( 	sc_stat[sc[ipart]] > 0		)				// Checks to see if SC is good for this run
	if(	hasElectronWith(ipart)		)
	if(	isP_vertexCut(ipart)			)
	if(	isDCFiducialCut(ipart)		)
	if(	isP_betaCut(ipart)			)
	if( 	pid[ipart] == 2212			)	  			// Nonzero velocity 
	if(	isP_thetaPhiCut(ipart)		)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isPhoton(Int_t ipart){
	if( 	stat[ipart] > 0				)				// Checks to see if it is a good run
	if( 	!(dc_stat[dc[ipart]] > 0)	)				// Checks to see if DC is good for this run
	if( 	sc_stat[sc[ipart]] > 0		)				// Checks to see if DC is good for this run
	if( 	q[ipart] == 0					)	  			// Charge Cut
	if( 	b[ipart] != 0					)	  			// Nonzero Beta Cut
	if(	hasElectronWith(ipart)		)
	if(	isG_ECEnergyCut(ipart)		)				// EC Energy Cut
	if(	isG_ECFiducialCut(ipart)	)				// EC Fiducial Cut
	if(	isG_vertexCut(ipart)			)		
	if(	isG_betaCut(ipart)			)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4(Int_t ipart){
	if( 	q[ipart] > 0								)	  			// Charge Cut
	if( 	dc_stat[dc[ipart]] > 0					)				// Checks to see if DC is good for this run
	if( 	sc_stat[sc[ipart]] > 0					)				// Checks to see if SC is good for this run
	if(	hasElectronWith(ipart)					)
	if(	isHe4_enoughPads(ipart)							)
	if(	isHe4_vertexCut(ipart)					)
	if(	isHe4_trackHelix(ipart)					)
	if(	isHe4_sDist(ipart)						)
	if(	isHe4_eDist(ipart)						)
	if(	isHe4_RTPCFiducialCut(ipart)			)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isPi0(Int_t ipart){
	if( 	q[ipart] == 0								)	  			// Charge Cut
	if( 	!(dc_stat[dc[ipart]] > 0)				)				// Checks to see if DC is good for this run
	if(	hasElectronWith(ipart)					)
	if(	isPi0_ICICCheck(ipart)					)
	if(	isPi0_ICECCheck(ipart)					)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t hasElectronWith(Int_t ipart){
//	return true;
	// Checks to see if there is an electron in the event	
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isPi0_ICECCheck(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isPi0_ICICCheck(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_RTPCFiducialCut(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_eDist(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_sDist(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_trackHelix(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_vertexCut(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(	isVertexMatch(ipart) 		)
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isHe4_enoughPads(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isG_ECEnergyCut(Int_t ipart){
//	return true;
	// EC Energy Cuts
	// Checks to see if the particle has enough energy deposited into the EC 
	// Sampling Fraction
	if( ec_eo[ec[ipart]] > 0 )				//	Energy Positive
	if( ec_ei[ec[ipart]] > 0 )				//		(in GeV)
	if( etot[ec[ipart]]	 > 0 )
	{
		Double_t etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];
		Double_t corrScale = etote/0.273; 
		if( corrScale > 0.3)
		{
			return true;
		}
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isG_ICFiducialCut(Int_t ipart){
//	return true;
	//	Cut photons produced by Møller electrons
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isG_betaCut(Int_t ipart){
//	return true;
	// Discriminates the photon from very fast electrons that both deposit energy into EC
//	TVector3 ecHitPos = TVector3(ech_x[ec[ipart]],ech_y[ec[ipart]],ech_z[ec[ipart]]);			// EC Hit Position
//	TVector3 eVertex = TVector3(vx[ipart], vy[ipart], eVz);
//	Double_t pathLength = (ecHitPos - eVertex).Mag();

//	Double_t ecTime = ec_t[ec_part];
//	Double_t trigTime = tr_time;
//	Double_t travTime = ecTime - trigTime;
//	Double_t beta = pathLength / 30 / travTime;		//  b = l / ct
	
	//!
	if(b[ipart] > 0.93 )
	if(b[ipart] < 1.07 )
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isG_ECFiducialCut(Int_t ipart){
//	return true;
	// Cuts all photons at edges of EC 
	TVector3 ecHitXYZ;
	ecHitXYZ.SetXYZ(ech_x[ec[ipart]],ech_y[ec[ipart]],ech_z[ec[ipart]]);			// EC Hit Position
   TVector3 ecHitUVW;
	ecHitUVW = getUVWfromXYZ(ecHitXYZ);																		// Convert XYZ coordinates to UVW

	float uu = ecHitUVW(0);
   float vv = ecHitUVW(1);
   float ww = ecHitUVW(2);
	if( uu > 100.	)
	if( uu < 390.	)
	if( vv < 360.	)
	if( ww < 390.	){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isG_vertexCut(Int_t ipart){
//	return true;
	// Correct for particle's z-vertex (if needed) and check to see if it coincides with an electron's z-vertex
	Double_t partVz = vz[ipart];	// Calculate particle's z vertex 
	if(	isVertexMatch(partVz) 	){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isP_betaCut(Int_t ipart){
//	 return true;
	// Discriminate proton from much faster pions 

	Double_t pp = p[ipart]; 
	// Double_t tof = sc_t[ipart];
	// Double_t pathLength = sc_r[ipart];
	Double_t betaSC = b[ipart];												// beta measured from sc (pathLength/30/tof;)
	Double_t betaDC = pp/sqrt(pow(pp,2) + pow(0.93827,2));		// beta calculated from DC
	Double_t dBeta = betaSC - betaDC;	
	//!
	if( Abs(dBeta) < 0.025 )
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isP_vertexCut(Int_t ipart){
//	 return true;
	// Make sure the protons are coming from where we expect them to be 
	TVector3 pDir = TVector3( cx[ipart], cy[ipart], cz[ipart]);
	double pTheta = pDir.Theta();
	double pPhi = pDir.Phi();
	double vzShift = getCorrectedVzShift(pTheta, pPhi);

	if(vzShift <= -999){
		return false;
	}
	
	pVz = vz[ipart] - vzShift;

	if( pVz < -50.0 )
	if( pVz > -77.0 )
	if( isVertexMatch(pVz) )		// Checks to see if the proton's vertex coincides with the electron vertex ( to within 2 cm )
	{	
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isP_thetaPhiCut(Int_t ipart){
//	return true;
	// Make sure the proton wasn't reconstructed from a region that gives bad data
	//!
	if(true)
	{
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isE_VertexCut(Int_t ipart){
//	 return true;
	// Corrects the electron vertex based on run number then makes sure the electrons are within the target walls 
	TVector3 pDir = TVector3( cx[ipart], cy[ipart], cz[ipart]);
	double pTheta = pDir.Theta();
	double pPhi = pDir.Phi();
	double vzShift = getCorrectedVzShift(pTheta, pPhi);
	if(vzShift <= -999){
		return false;
	}
	eVz = vz[ipart] - vzShift;

	if( eVz < -54.0)
	if( eVz > -74.0){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isE_ECFiducialCut(Int_t ipart){
//	return true;
	// Makes sure reconstructed electrons are not coming from edges of the EC
   ////////Set the EC_XYZ to EC_UVW  ---------------------------------
	TVector3 ecHitXYZ;
	ecHitXYZ.SetXYZ(ech_x[ec[ipart]],ech_y[ec[ipart]],ech_z[ec[ipart]]);			// EC Hit Position
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

Bool_t isE_ECEnergyCut(Int_t ipart){
//	return true;
	// EC Energy Cuts
	// Checks to see if the particle has enough energy deposited into the EC 
	// Sampling Fraction
	if( ec_eo[ec[ipart]] > 0 )				//	Energy Positive
	if( ec_ei[ec[ipart]] > 0 )				//		(in GeV)
	if( etot[ec[ipart]]	 > 0 )
	{
		Float_t looseness = 2.5;

		Double_t ucoefs[] = { 0.2560840,  0.0432374, -0.00914180,  0.00081589}; 
		Double_t ocoefs[] = { 0.0572976, -0.0272689,  0.00857600, -0.00097998}; 
		
		Double_t mu = 0;
		Double_t sigma = 0;
		for(int ii = 0; ii < sizeof(ucoefs)/sizeof(ucoefs[0]); ii++){mu += ucoefs[ii]*pow(p[ipart],ii); sigma += ocoefs[ii]*pow(p[ipart],ii);}					
		Double_t etote = (etot[ec[ipart]] < ec_ei[ec[ipart]] + ec_eo[ec[ipart]]) ? ec_ei[ec[ipart]] + ec_eo[ec[ipart]] : etot[ec[ipart]];
		Double_t sampFrac = etote/p[ipart];

		Int_t sector = dc_sect[dc[ipart]];
		Double_t corrScale = getCorrectedSampFrac((Int_t) runnb, (Int_t) evntid, (Int_t) sector, curFile); 
		if( corrScale > 0.1)
		{ 
         sampFrac = (sampFrac)* (0.31/corrScale);
		}
		if( sampFrac < mu + looseness * sigma)
		if( sampFrac > mu - looseness * sigma){
			return true;
		}
		return false;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isDCFiducialCut(Int_t ipart){
//	return true;
	// DC Fiducial Cut: Checks to see if the particle is in the fiducial region for the DC
	TVector3 dcTrackPos = TVector3(tl1_x[ipart], tl1_y[ipart], tl1_z[ipart]);
	TVector3 dcTrackDir = TVector3(tl1_cx[ipart], tl1_cy[ipart], tl1_cz[ipart]);
	shiftICtoDC = getICtoDCShift(dcTrackPos, dcTrackDir);

	TVector3 shiftedPos = dcTrackPos - shiftICtoDC;
	
	Double_t X = shiftedPos.X();
	Double_t Y = shiftedPos.Y();
	Int_t		S = dc_sect[dc[ipart]];

	if(q[ipart] > 0)			// This is the DC Mid Cut used on protons only
	{
  		if( X > -0.2 )
  		if( X <  0.6 )
  		if( Y > -0.2 ) 
		if( Y <  0.6 )
		{	return false;	}
	}

	using namespace TMath;
	if( !isInsideIConDCShadow(X, Y) ){											// Checks to see if the shifted position is in the IC shadow
		Double_t goodRelAngL  = ((S-1)-1./3)*Pi()/3;							//		Makes sure angle is within the left and right good relative angles depending on section
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
//	return true;
	// Checks to see if the particle is inside the IC shadow
	///// NOTE: This cut requires the shifted DC track X and Y positions. 
	          //      4        5
	          //   
	          // 3                    6
	          //   
	          // 2                    7
	          //      1        8
	          //      0        9
	          // the same initial point  
//	TCutG *geoCut = new TCutG("geocut", 10); 
	TCutG geoCut = TCutG("geocut", 10); 
	// First row of cutPoints is x-coordinates, second is y-coordinates
	float cutPoints[2][11] = {	{-11.15, -11.15, -23.1, -23.1, -10.3, 9.91, 23.73, 23.73, 12.3, 12.3, -11.15},
										{-26.07, -23.1, -12.85, 11.5, 22.95, 22.95, 13.1, -12.4, -22.36, -26.07, -26.07}	};
	
	for(Int_t ii = 0; ii < 11; ii++) 
	{
		//geoCut->SetPoint(ii, cutPoints[0][ii], cutPoints[1][ii]);
		geoCut.SetPoint(ii, cutPoints[0][ii], cutPoints[1][ii]);
	}
	if(geoCut.IsInside(X, Y)){
		return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Bool_t isE_CCFiducialCut(Int_t ipart){
//	return true;
	// CC Fiducial Cuts
	// Makes sure the particles are not hitting the edges of the CC

	TVector3 trackECIntPos = TVector3(dc_xsc[dc[ipart]], dc_ysc[dc[ipart]], dc_zsc[dc[ipart]]);		// Point of track and EC-plane intersection 
	TVector3 trackECIntDir = TVector3(dc_cxsc[dc[ipart]], dc_cysc[dc[ipart]], dc_czsc[dc[ipart]]);		// Direction of track and EC-plane intersection 
	
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
		phiEdge2 = 20*TMath::Sqrt(0.5*(ccTheta - 43.)); }

	if( fabs(ccPhi) < phiEdge1 )
	if( fabs(ccPhi) > phiEdge2 ){
		return true;
	}
	return false;
}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Bool_t EG5_test::goodDetectors(){
// return true;
	// Makes sure all of the current event's status are good
	if(dc_stat[dc_part] > 0)
	if(ec_stat[ec_part] > 0)
	if(sc_stat[sc_part] > 0){
		return true;
	}
	return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif

#ifdef EG5_test_cxx
EG5_test::EG5_test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      //f->GetObject("ch",tree);
      f->GetObject("ch",tree);
#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("h10_1e","EG6_test");
     // chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_61550*.root/h10_1e"); //small
     chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_6155*.root/h10_1e"); //medium
	  //	chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_61*.root/h10_1e"); //large
     //chain->Add("data/*"); //small
	  tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

EG5_test::~EG5_test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EG5_test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EG5_test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EG5_test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("runnb", &runnb, &b_runnb);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
   fChain->SetBranchAddress("evntype", &evntype, &b_evntype);
   fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("rf_time", &rf_time, &b_rf_time);
   fChain->SetBranchAddress("l2bit", &l2bit, &b_l2bit);
   fChain->SetBranchAddress("l3bit", &l3bit, &b_l3bit);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("hlsc", &hlsc, &b_hlsc);
   fChain->SetBranchAddress("intt", &intt, &b_intt);
   fChain->SetBranchAddress("helicity_cor", &helicity_cor, &b_helicity_cor);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
   fChain->SetBranchAddress("lec", lec, &b_lec);
   fChain->SetBranchAddress("st", st, &b_st);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("b", b, &b_b);
   fChain->SetBranchAddress("cx", cx, &b_cx);
   fChain->SetBranchAddress("cy", cy, &b_cy);
   fChain->SetBranchAddress("cz", cz, &b_cz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("tb_st", tb_st, &b_tb_st);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
   fChain->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
   fChain->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
   fChain->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
   fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
   fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
   fChain->SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
   fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
   fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
   fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
   fChain->SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
   fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
   fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
   fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
   fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
   fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
   fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
   fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("edep", edep, &b_edep);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);
   fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
   fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
   fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
   fChain->SetBranchAddress("ic_part", &ic_part, &b_ic_part);
   fChain->SetBranchAddress("et", et, &b_et);
   fChain->SetBranchAddress("egl", egl, &b_egl);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("time_next", time_next, &b_time_next);
   fChain->SetBranchAddress("ich_x", ich_x, &b_ich_x);
   fChain->SetBranchAddress("ich_y", ich_y, &b_ich_y);
   fChain->SetBranchAddress("ich_z", ich_z, &b_ich_z);
   fChain->SetBranchAddress("ich_xgl", ich_xgl, &b_ich_xgl);
   fChain->SetBranchAddress("ich_ygl", ich_ygl, &b_ich_ygl);
   fChain->SetBranchAddress("ich_xwid", ich_xwid, &b_ich_xwid);
   fChain->SetBranchAddress("ich_ywid", ich_ywid, &b_ich_ywid);
   fChain->SetBranchAddress("ich_xm3", ich_xm3, &b_ich_xm3);
   fChain->SetBranchAddress("ich_ym3", ich_ym3, &b_ich_ym3);
   fChain->SetBranchAddress("ic_stat", ic_stat, &b_ic_stat);
   fChain->SetBranchAddress("gcpart", &gcpart, &b_gcpart);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("dedx", dedx, &b_dedx);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("p_tot", p_tot, &b_p_tot);
   fChain->SetBranchAddress("x2", x2, &b_x2);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("dca", dca, &b_dca);
   fChain->SetBranchAddress("index", index, &b_index);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("vtl", vtl, &b_vtl);
   fChain->SetBranchAddress("sdist", sdist, &b_sdist);
   fChain->SetBranchAddress("edist", edist, &b_edist);
   fChain->SetBranchAddress("npts", npts, &b_npts);
   fChain->SetBranchAddress("r_0", r_0, &b_r_0);
   fChain->SetBranchAddress("fiterr", fiterr, &b_fiterr);
   fChain->SetBranchAddress("tothits", tothits, &b_tothits);
   fChain->SetBranchAddress("npd_track", npd_track, &b_npd_track);
   fChain->SetBranchAddress("npd_event", npd_event, &b_npd_event);
   fChain->SetBranchAddress("bonus_bits", bonus_bits, &b_bonus_bits);
   fChain->SetBranchAddress("q_tot", q_tot, &b_q_tot);
   fChain->SetBranchAddress("x_start", x_start, &b_x_start);
   fChain->SetBranchAddress("y_start", y_start, &b_y_start);
   fChain->SetBranchAddress("z_start", z_start, &b_z_start);
   fChain->SetBranchAddress("x_end", x_end, &b_x_end);
   fChain->SetBranchAddress("y_end", y_end, &b_y_end);
   fChain->SetBranchAddress("z_end", z_end, &b_z_end);
   fChain->SetBranchAddress("rtpc_npart", &rtpc_npart, &b_rtpc_npart);
   fChain->SetBranchAddress("rtpc_id1", rtpc_id1, &b_rtpc_id1);
   fChain->SetBranchAddress("rtpc_id2", rtpc_id2, &b_rtpc_id2);
   fChain->SetBranchAddress("rtpc_id3", rtpc_id3, &b_rtpc_id3);
   fChain->SetBranchAddress("rtpc_id4", rtpc_id4, &b_rtpc_id4);
   fChain->SetBranchAddress("rtpc_id5", rtpc_id5, &b_rtpc_id5);
   fChain->SetBranchAddress("rtpc_p1", rtpc_p1, &b_rtpc_p1);
   fChain->SetBranchAddress("rtpc_p2", rtpc_p2, &b_rtpc_p2);
   fChain->SetBranchAddress("rtpc_p3", rtpc_p3, &b_rtpc_p3);
   fChain->SetBranchAddress("rtpc_p4", rtpc_p4, &b_rtpc_p4);
   fChain->SetBranchAddress("rtpc_p5", rtpc_p5, &b_rtpc_p5);
   fChain->SetBranchAddress("rtpc_poverq", rtpc_poverq, &b_rtpc_poverq);
   fChain->SetBranchAddress("rtpc_dedx", rtpc_dedx, &b_rtpc_dedx);
   fChain->SetBranchAddress("rtpc_dedx2", rtpc_dedx2, &b_rtpc_dedx2);
   fChain->SetBranchAddress("rtpc_dedxa", rtpc_dedxa, &b_rtpc_dedxa);
   fChain->SetBranchAddress("rtpc_dedxl", rtpc_dedxl, &b_rtpc_dedxl);
   fChain->SetBranchAddress("rtpc_dedxal", rtpc_dedxal, &b_rtpc_dedxal);
   fChain->SetBranchAddress("rtpc_rxy", rtpc_rxy, &b_rtpc_rxy);
   fChain->SetBranchAddress("rtpc_slope", rtpc_slope, &b_rtpc_slope);
   fChain->SetBranchAddress("rtpc_chisq", rtpc_chisq, &b_rtpc_chisq);
   fChain->SetBranchAddress("rtpc_dedxs", rtpc_dedxs, &b_rtpc_dedxs);
   fChain->SetBranchAddress("rtpc_theta", rtpc_theta, &b_rtpc_theta);
   fChain->SetBranchAddress("rtpc_phi", rtpc_phi, &b_rtpc_phi);
   fChain->SetBranchAddress("rtpc_vz", rtpc_vz, &b_rtpc_vz);
   fChain->SetBranchAddress("rtpc_bad", rtpc_bad, &b_rtpc_bad);
   fChain->SetBranchAddress("rtpc_gcpb", rtpc_gcpb, &b_rtpc_gcpb);
   fChain->SetBranchAddress("icpart", &icpart, &b_icpart);
   fChain->SetBranchAddress("etc", etc, &b_etc);
   fChain->SetBranchAddress("ecc", ecc, &b_ecc);
   fChain->SetBranchAddress("tc", tc, &b_tc);
   fChain->SetBranchAddress("tn", tn, &b_tn);
   fChain->SetBranchAddress("xc", xc, &b_xc);
   fChain->SetBranchAddress("yc", yc, &b_yc);
   fChain->SetBranchAddress("zc", zc, &b_zc);
   fChain->SetBranchAddress("m2c", m2c, &b_m2c);
   fChain->SetBranchAddress("m3c", m3c, &b_m3c);
   fChain->SetBranchAddress("statc", statc, &b_statc);
   fChain->SetBranchAddress("shh_part", &shh_part, &b_shh_part);
   fChain->SetBranchAddress("shh_id", shh_id, &b_shh_id);
   fChain->SetBranchAddress("shh_x", shh_x, &b_shh_x);
   fChain->SetBranchAddress("shh_y", shh_y, &b_shh_y);
   fChain->SetBranchAddress("shh_z", shh_z, &b_shh_z);
   fChain->SetBranchAddress("shh_nphe", shh_nphe, &b_shh_nphe);
   fChain->SetBranchAddress("shh_time", shh_time, &b_shh_time);
   fChain->SetBranchAddress("shh_stat", shh_stat, &b_shh_stat);
   fChain->SetBranchAddress("shpart", &shpart, &b_shpart);
   fChain->SetBranchAddress("shid", shid, &b_shid);
   fChain->SetBranchAddress("shx", shx, &b_shx);
   fChain->SetBranchAddress("shy", shy, &b_shy);
   fChain->SetBranchAddress("shz", shz, &b_shz);
   fChain->SetBranchAddress("shnphe", shnphe, &b_shnphe);
   fChain->SetBranchAddress("shtime", shtime, &b_shtime);
   fChain->SetBranchAddress("shstat", shstat, &b_shstat);
   Notify();
}

Bool_t EG5_test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EG5_test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EG5_test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EG5_test_cxx
