// ------------------------------------------------------------
//  AnalyzerMakinTree.C – UPDATED VERSION
// ------------------------------------------------------------

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <filesystem>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "lardataobj/RawData/OpDetWaveform.h"

using namespace art;
using namespace std;
namespace fs = std::filesystem;

// ------------------------------------------------------------
//  Delay table 
// ------------------------------------------------------------
std::map<int, float> chanID_to_delayNS = {
    {250, 5.44},   {304, 1.22},   {274, 4.1},    {248, 0.25},   {302, 5.54},   {246, 3.18},
    {300, 3.63},   {272, 2.75},   {244, 3.22},   {298, 3.9},    {242, 3.8},    {296, 0},
    {270, 0.76},   {240, 1.8},    {294, 0.1},    {172, 1.95},   {226, 1.58},   {196, 0.2},
    {170, 1.9},    {224, 3.77},   {168, -0.64},  {222, -0.9},   {194, 2.09},   {166, 2.3},
    {220, 3.85},   {164, 1.42},   {218, -0.35},  {192, 0.16},   {162, -0.6},   {216, 5.61},
    {94, 1.75},    {148, -0.93},  {118, -1.9},   {92, 1.2},     {146, 0.63},   {90, 5.15},
    {144, 3.25},   {116, -2.67},  {88, 2.13},    {142, -2},     {86, 0},       {140, -0.05},
    {114, 0.52},   {84, -1.45},   {138, 1.63},   {16, 1.97},    {70, 5.67},    {40, 5.45},
    {14, 3.7},     {68, 5.62},    {12, 1},       {66, 0.65},    {38, 1.93},    {10, -2.25},
    {64, 5.63},    {8, -10.22},   {62, -4.05},   {36, 0.45},    {6, -2.99},    {60, -3.85},
    {71, 0.5},     {17, -3.32},   {41, -1.74},   {69, 2.38},    {15, 0.3},     {67, 2},
    {13, -2.74},   {39, 2.2},     {65, -1.15},   {11, 5.62},    {63, 0.43},    {9, -4.26},
    {37, -3.83},   {61, -1.53},   {7, -3.07},    {149, 2.4},    {95, -1},      {119, 4.2},
    {147, 1.4},    {93, -0.7},    {145, 4.05},   {91, 4.07},    {117, 2.1},    {143, 1.6},
    {89, 3.81},    {141, -3.05},  {87, 0.08},    {115, -0.68},  {139, -1.76},  {85, -11.41},
    {227, 4.16},   {173, 1.71},   {197, 3.28},   {225, -2.26},  {171, -2.45},  {223, 1.91},
    {169, 0.42},   {195, 4.04},   {221, 1.25},   {167, 1.07},   {219, -1.25},  {165, 2},
    {193, -1.7},   {217, 1.27},   {163, 1.4},    {305, 1.27},   {251, 4.27},   {275, 1.84},
    {303, 5.18},   {249, -1.77},  {301, 1.19},   {247, 3.09},   {273, -0.21},  {299, 3.63},
    {245, 1.75},   {297, 5.25},   {243, 5.28},   {271, 0.75},   {295, -2.98},  {241, 1.24}
};

// ------------------------------------------------------------
//  Helper: derivative only around the pulse
// ------------------------------------------------------------
std::vector<std::pair<int,float>> Derivative_at_tick(const raw::OpDetWaveform& opt)
{
  auto itMin = std::min_element(opt.Waveform().begin(), opt.Waveform().end());
  int iMin = std::distance(opt.Waveform().begin(), itMin);

  std::vector<float> num, den;
  std::vector<std::pair<int,float>> out;

  for (int tick = iMin-50; tick < iMin+30; ++tick) {
    if (tick < 0 || tick >= (int)opt.Waveform().size()) continue;

    num.push_back(opt.Waveform()[tick]);
    den.push_back(opt.Waveform()[tick]);

    if (num.size() > 5)  num.erase(num.begin());
    if (den.size() > 10) den.erase(den.begin());

    if (num.size()==5 && den.size()==10) {
      float sumDen = std::accumulate(den.begin(), den.end(), 0.f);
      if (sumDen != 0.f) {
        float deriv = std::accumulate(num.begin(), num.end(), 0.f) / sumDen;
        out.emplace_back(tick-5, deriv);
      }
    }
  }
  return out;
}

// ------------------------------------------------------------
//  Main function
// ------------------------------------------------------------
void AnalyzerMakinTree(const std::string& input)
{
  std::vector<std::string> filenames{input};
  std::cout << "Input file: " << input << "\n" << std::flush;

  // ----- Extract base file name -----
  size_t lastSlash = input.find_last_of('/');
  std::string basename = (lastSlash == std::string::npos) ? input : input.substr(lastSlash + 1);

  size_t dotPos = basename.rfind(".root");
  if (dotPos != std::string::npos) basename.erase(dotPos);

  // ----- Auto plot directory -----
  std::string plotDir = "Plots-" + basename;
  gSystem->mkdir(plotDir.c_str(), true);

  // ----- Output ROOT file -----
  std::string outName = "outfile_ana-" + basename + ".root";
  std::cout << "Output file: " << outName << "\n" << std::flush;

  TFile* outFile = new TFile(outName.c_str(), "RECREATE");

  // ----- Tree -----
  TTree* tree = new TTree("PMTWaveformTree", "PMT waveform + derivative");
  Int_t Run, Event, Channel, StartTick_Waveform, StartTick_Derivative;
  std::vector<int>   *Waveform   = new std::vector<int>;
  std::vector<float> *Derivative = new std::vector<float>;

  tree->Branch("Run",                 &Run);
  tree->Branch("Event",               &Event);
  tree->Branch("Channel",             &Channel);
  tree->Branch("StartTick_Waveform",  &StartTick_Waveform);
  tree->Branch("StartTick_Derivative",&StartTick_Derivative);
  tree->Branch("Waveform",            &Waveform);
  tree->Branch("Derivative",          &Derivative);

  // ----- Histograms -----
  const int NMAX = 2000000;
  TH1F* hWave = new TH1F("hWave", "Raw waveform (all channels)", NMAX, 0, NMAX);
  TH1F* hDer  = new TH1F("hDer",  "Derivative (full)",           NMAX, 0, NMAX);
  TH1F* hDerPulse = new TH1F("hDerPulse", "Derivative near pulse", 401, -200.5, 200.5);
  TH2F* hHeat = new TH2F("hHeat",
                         "Waveform heatmap;Tick;ADC",
                         1000, 0, 5000,
                         800, 12000, 16000);

  hWave->SetDirectory(0);
  hDer->SetDirectory(0);
  hDerPulse->SetDirectory(0);
  hHeat->SetDirectory(0);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111);

  // ---------- Event loop ----------
  int evtCnt = 0;
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {

    std::cout << "=== Event " << evtCnt << " ===\n" << std::flush;

    auto const& waveforms = *ev.getValidHandle<std::vector<raw::OpDetWaveform>>("pmtdecoder:PMTChannels");
    Run   = ev.eventAuxiliary().run();
    Event = ev.eventAuxiliary().event();

    for (const auto& wf : waveforms) {
      Channel = wf.ChannelNumber();
      const auto& adc = wf.Waveform();

      Waveform->clear();
      Derivative->clear();
      StartTick_Waveform = 0;
      StartTick_Derivative = 10;

      int pulseTick = std::distance(adc.begin(), std::min_element(adc.begin(), adc.end()));

      for (auto [t,val] : Derivative_at_tick(wf))
        hDerPulse->Fill(t - pulseTick, val);

      std::vector<float> num, den;
      for (size_t t = 0; t < adc.size(); ++t) {

        int val = adc[t];
        hWave->Fill(t, val);
        if (t < 5000) hHeat->Fill(t, val);

        num.push_back(val);
        den.push_back(val);
        if (num.size() > 5)  num.erase(num.begin());
        if (den.size() > 10) den.erase(den.begin());

        if (t >= 10 && num.size()==5 && den.size()==10) {
          float sumDen = std::accumulate(den.begin(), den.end(), 0.f);
          if (sumDen != 0.f) {
            float deriv = std::accumulate(num.begin(), num.end(), 0.f) / sumDen;
            hDer->Fill(t-5, deriv);
            Derivative->push_back(deriv);
          }
        }
        Waveform->push_back(val);
      }
      tree->Fill();
    }
    ++evtCnt;
  }

  // ---------- Save plots ----------
  auto save = [&](TH1* h, const char* name, const char* opt = "") {
    TCanvas c("c","c", 1400, 700);   // *** wider canvas ***
    c.SetLeftMargin(0.12);           // *** avoids clipping "ADC" on heatmap ***
    c.SetRightMargin(0.12);          // *** room for colorbar ***
    c.cd();
    if (strlen(opt)) h->Draw(opt);
    else h->Draw();
    c.SaveAs(Form("%s/%s.png", plotDir.c_str(), name));
  };

  hWave->GetXaxis()->SetRangeUser(0, 6000);
  hWave->GetYaxis()->SetRangeUser(12000, 16000);
  save(hWave, "waveform", "hist");

  hDer->GetXaxis()->SetRangeUser(0, 6000);
  save(hDer, "derivative", "hist");

  save(hDerPulse, "derivative_pulse", "hist");

  hHeat->SetStats(0);
  save(hHeat, "heatmap", "colz");

  // ----- Write and close -----
  outFile->cd();
  tree->Write();
  outFile->Close();

  std::cout << "\n=== ALL DONE ===\n";
  std::cout << "Tree saved as: " << outName << "\n";
  std::cout << "Plots saved in: " << plotDir << "/\n" << std::flush;
}
