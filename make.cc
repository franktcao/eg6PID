#include <TChain.h>


void make()
{
  TChain *ch = new TChain("ch", "EG6_test");
  ch->Add("/music/clas1/clas/eg6-data/pass2/6GeV/pass2v0/hroot_1e_61551*.root/h10_1e"); 
  ch->MakeClass("EG5_test");
}
