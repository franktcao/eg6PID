void printHist(String filename){
// Let's open the TFile
TFile in_file(filename);
// Get the Histogram out
TH2F* h;
in_file.GetObject("hbvp",h);
// Draw it
h->Draw();
}
