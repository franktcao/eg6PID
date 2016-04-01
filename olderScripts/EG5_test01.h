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
#include <deque> 

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

	TH2D 			*hCuts2[5];
	TH1D 			*hCuts1[5];
	const int	*h2length;
	const int	*h1length;
	TObjArray	Hlist2;
	TObjArray	Hlist1;


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
	virtual void	  writeFiles();
	virtual Bool_t	  isElectron(Int_t particleIndex);
};

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
      f->GetObject("ch",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("ch","EG6_test");
      chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_61550*.root/h10_1e"); //small
     // chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_6155*.root/h10_1e"); //medium
	//	chain->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_61*.root/h10_1e"); //large
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
