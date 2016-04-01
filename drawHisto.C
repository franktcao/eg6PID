#define drawHist_cxx
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TCanvas.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

void drawHist(TFile inputfile){
	TH1F* h = nullptr;
	inputFile.GetObject("hbvp",h);
	TCanvas can = new TCanvas;
	h->Draw();
}
