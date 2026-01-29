//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 29 16:07:01 2024 by ROOT version 6.28/12
// from TTree PMTWaveformTree/PMTWaveformTree
// found on file: ../outfile.root
//////////////////////////////////////////////////////////

#ifndef PMTWaveformTree_h
#define PMTWaveformTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class PMTWaveformTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           Channel;
   Int_t           StartTick_Waveform;
   vector<int>     *Waveform;
   Int_t           StartTick_Derivative;
   vector<float>   *Derivative;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Channel;   //!
   TBranch        *b_StartTick_Waveform;   //!
   TBranch        *b_Waveform;   //!
   TBranch        *b_StartTick_Derivative;   //!
   TBranch        *b_Derivative;   //!

   PMTWaveformTree(TTree *tree=0);
   virtual ~PMTWaveformTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual int bin_of_halfmin(std::vector<float> der);
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef PMTWaveformTree_cxx
PMTWaveformTree::PMTWaveformTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/exp/sbnd/data/users/castalyf/PMT_timing/outfile_ana-data_EventBuilder6_art1_run19406_1_20250924T161002_decoded-filtered.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/exp/sbnd/data/users/castalyf/PMT_timing/outfile_ana-data_EventBuilder6_art1_run19406_1_20250924T161002_decoded-filtered.root");
      }
      f->GetObject("PMTWaveformTree",tree);

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
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PMTWaveformTree::LoadTree(Long64_t entry)
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

void PMTWaveformTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Waveform = 0;
   Derivative = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Channel", &Channel, &b_Channel);
   fChain->SetBranchAddress("StartTick_Waveform", &StartTick_Waveform, &b_StartTick_Waveform);
   fChain->SetBranchAddress("Waveform", &Waveform, &b_Waveform);
   fChain->SetBranchAddress("StartTick_Derivative", &StartTick_Derivative, &b_StartTick_Derivative);
   fChain->SetBranchAddress("Derivative", &Derivative, &b_Derivative);
   Notify();
}

Bool_t PMTWaveformTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

int PMTWaveformTree::bin_of_halfmin(std::vector<float> der){

  double mean = std::reduce(der.begin(),der.end())/double(der.size());
  double min = *std::min_element(der.begin(), der.end());

  double half_min =  mean - (mean - min)/double(2);

  auto it = std::find_if(der.begin(), der.end(), [half_min](double val) {
      return val < half_min;
    });

  return std::distance(der.begin(), it);

}

void PMTWaveformTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PMTWaveformTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PMTWaveformTree_cxx
