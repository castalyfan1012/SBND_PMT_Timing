//////////////////////////////////////////////////////////
// PMTWaveformTree.h
//
// MODIFIED from auto-generated version:
//   - Constructor auto-loads all EventBuilder*.root via
//     wildcard TChain when no tree is supplied.
//   - Loop() takes diffuser_filter argument.
//////////////////////////////////////////////////////////

#ifndef PMTWaveformTree_h
#define PMTWaveformTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <iostream>
#include <numeric>
#include <algorithm>
#include "vector"

class PMTWaveformTree {
public :
   TTree          *fChain;
   Int_t           fCurrent;

   Int_t           Run;
   Int_t           Event;
   Int_t           Channel;
   Int_t           StartTick_Waveform;
   vector<int>     *Waveform;
   Int_t           StartTick_Derivative;
   vector<float>   *Derivative;

   TBranch        *b_Run;
   TBranch        *b_Event;
   TBranch        *b_Channel;
   TBranch        *b_StartTick_Waveform;
   TBranch        *b_Waveform;
   TBranch        *b_StartTick_Derivative;
   TBranch        *b_Derivative;

   PMTWaveformTree(TTree *tree=0);
   virtual ~PMTWaveformTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   // diffuser_filter: "ALL"    -> all diffusers (default)
   //                  "West_3","East_1",... -> single diffuser only
   virtual void     Loop(const std::string& diffuser_filter = "ALL");
   virtual Bool_t   Notify();
   virtual int      bin_of_halfmin(std::vector<float> der);
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PMTWaveformTree_cxx

PMTWaveformTree::PMTWaveformTree(TTree *tree) : fChain(0)
{
   if (tree == 0) {
      const char* pattern =
         "/pnfs/sbnd/persistent/users/castalyf/LED/PMT_outfiles/"
         "outfile_ana-data_EventBuilder*.root";

      TChain* ch = new TChain("PMTWaveformTree");
      int nfiles = ch->Add(pattern);

      if (nfiles == 0) {
         std::cerr << "[PMTWaveformTree] WARNING: no files matched: " << pattern << "\n"
                   << "  Pass a TChain/TTree explicitly if needed.\n";
      } else {
         std::cout << "[PMTWaveformTree] Auto-loaded " << nfiles << " file(s).\n";
      }
      tree = ch;
   }
   Init(tree);
}

PMTWaveformTree::~PMTWaveformTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PMTWaveformTree::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t PMTWaveformTree::LoadTree(Long64_t entry)
{
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PMTWaveformTree::Init(TTree *tree)
{
   Waveform   = 0;
   Derivative = 0;
   if (!tree) return;
   fChain   = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run",                  &Run,                  &b_Run);
   fChain->SetBranchAddress("Event",                &Event,                &b_Event);
   fChain->SetBranchAddress("Channel",              &Channel,              &b_Channel);
   fChain->SetBranchAddress("StartTick_Waveform",   &StartTick_Waveform,   &b_StartTick_Waveform);
   fChain->SetBranchAddress("Waveform",             &Waveform,             &b_Waveform);
   fChain->SetBranchAddress("StartTick_Derivative", &StartTick_Derivative, &b_StartTick_Derivative);
   fChain->SetBranchAddress("Derivative",           &Derivative,           &b_Derivative);
   Notify();
}

Bool_t PMTWaveformTree::Notify() { return kTRUE; }

int PMTWaveformTree::bin_of_halfmin(std::vector<float> der)
{
   double mean = std::reduce(der.begin(), der.end()) / double(der.size());
   double mn   = *std::min_element(der.begin(), der.end());
   double half_min = mean - (mean - mn) / 2.0;
   auto it = std::find_if(der.begin(), der.end(),
                          [half_min](double val){ return val < half_min; });
   return std::distance(der.begin(), it);
}

void PMTWaveformTree::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t PMTWaveformTree::Cut(Long64_t entry) { return 1; }

#endif // PMTWaveformTree_cxx