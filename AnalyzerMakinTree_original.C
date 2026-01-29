// /pnfs/uboone/data/uboone/reconstructed/prod_v08_00_00_18/data_extbnb_mcc9.1_v08_00_00_18/run1_reco2_C1/00/00/69/58/PhysicsRun-2016_7_24_11_54_40-0006958-00040_20160724T214243_ext_bnb_20160724T234647_merged_20181107T121731_optfilter_20181224T075125_reco1_postwcct_postdl_20181224T082906_20190723T182808_reco2.root


/// Psuedo code:
/*

  Find all the T0 from acpttrigtagger
  
  then grab the tracks that have an association through acpttrigtagger
  
  using these tracks get the associated caloritmetry data product

  then got through all the calo's trajectory points, mark it's x, then find the associated hit and store it's attributes 
  

 */


// Standard things to include
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h> 
// These are the includes to use "Root" things 
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

// These are the larsoft includes that let you
// have access to data-products and the event 
// details
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
//#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/fwd.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

//I'll need, calo, tracks, hits, anab::T0

#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "lardataobj/RawData/OpDetWaveform.h"


//This way you can be lazy
using namespace art;
using namespace std;

std::vector< std::pair<int, float> > Derivative_at_tick(raw::OpDetWaveform opt);

std::map<int, float> chanID_to_delayNS = {
  {17,-3.75},
  {305,0.8},
  {63,0},
  {297,4.8},
  {16,1.5},
  {36,0},
  {294,-0.35},
  {274,3.65},
  {15,-0.15},
  {65,-1.5},
  {9,-4.75},
  {147,1.05},
  {117,1.75},
  {139,1.55},
  {197,2.75},
  {167,0.6},
  {219,-1.6},
  {303,4.75},
  {299,3.2},
  {295,-3.45},
  {250,5.15},
  {272,2.4},
  {296,-0.35},
  {170,1.55},
  {194,1.6},
  {218,-0.7},
  {94,1.3},
  {116,-3.1},
  {140,-0.4},
  {40,5},
  {12,0.65},
  {60,-4.2} 
};

void AnalyzerMakinTree(string input){
  
  // create a vector of files we want to process
  std::vector<std::string> filenames{Form("%s",input.c_str())};
  
  // read in a file list that we get from a sam-def but remember it 
  // is very long so if we want to run over it all it'll take a while
  // and we'll probably want to break it up on the grid

  //ifstream myfile("elec_standard.list");
  /*
  ifstream myfile(Form("%s.list",input.c_str()));
  copy(istream_iterator<string>(myfile),
       istream_iterator<string>(),
       back_inserter(filenames));
  */
  //We'll just check the first N files for now 
  //  filenames.erase(filenames.begin()+250,filenames.end());
  
  std::cout << "File Count: " << filenames.size() << std::endl; 
  
  
  // Here we will create all of our histograms 
  // I did this crazy inefficiently but I don't really care
  // This is currently only set up for single dimensional 
  // projections but extenting it to 2D will be straight forward
  
  TFile* out = new TFile("outfile_testTree.root","RECREATE");  
/*  TTree* fTree = new TTree("waveformTree","tree");
  int EventNum;
  fTree->Branch("EventNum",&EventNum);
  int OpDetChannel;
  fTree->Branch("OpDetChannel",&OpDetChannel);
  int Tick;
  fTree->Branch("Tick",&Tick);
//  int WaveformValue;
//  fTree->Branch("WaveformValue",&WaveformValue);
  int Derivative;
  fTree->Branch("Derivative",&Derivative);
  int ns;
  fTree->Branch("ns",&ns);
  int calib_ns;
  fTree->Branch("calib_ns",&calib_ns);*/
/*  int Time;
  fTree->Branch("Time",&Time);
  int N_ADC_more_200;
  fTree->Branch("N_ADC_more_200",&N_ADC_more_200);
  int N_ADC_more_1000;
  fTree->Branch("N_ADC_more_1000",&N_ADC_more_1000);
  float RMS;
  fTree->Branch("RMS",&RMS);
  int Plane;
  fTree->Branch("Plane",&Plane);
  int FEMBOnWIB;       // 0:3
  fTree->Branch("FEMBOnWIB",&FEMBOnWIB);
  int FEMBCh;          // channel on FEMB -- 0:127
  fTree->Branch("FEMBCh",&FEMBCh);
  int asic;            // 0:7
  fTree->Branch("asic",&asic);
  int asicchan;        // ASIC channel:  0 to 15
  fTree->Branch("asicchan",&asicchan);
  int WIBCrate;        // 1:4
  fTree->Branch("WIBCrate",&WIBCrate);
  int WIB;             // 1:6
  fTree->Branch("WIB",&WIB);
  int WIBCh;           // 0:895   (7 FEMBs)
  fTree->Branch("WIBCh",&WIBCh);
  int WIBQFSP;         // 1:2
  fTree->Branch("WIBQFSP",&WIBQFSP);
  int QFSPFiber;       // 1:4
  fTree->Branch("QFSPFiber",&QFSPFiber);
  int FEMCrate;        // 1:11
  fTree->Branch("FEMCrate",&FEMCrate);
  int FEM;             // 1:16
  fTree->Branch("FEM",&FEM);
  int FEMCh;           // 0:63 channel in a FEM
  fTree->Branch("FEMCh",&FEMCh);
  int offlchan;        // in gdml and channel sorting convention
  fTree->Branch("offlchan",&offlchan);  
  */

  int e = 0; 
  
  TH1F* wf = new TH1F("wf","wf",1000000,0,1000000);
  TH1F* dev = new TH1F("dev","dev",1000000,0,1000000);
  TH1F* devfunc = new TH1F("devfunc","devfunc",1000000,0,1000000);
      TH2F* h_waveform_heatmap = new TH2F(Form("h_waveform_heatmap_ch%d",0), Form("h_waveform_heatmap_ch%d",0), 400+1, 48000-0.5, 52000+0.5, 1000+1, 7000-0.5, 17000+0.5); //Make a histogram of the waveform for each channel

  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    std::cout << "Event number: " << e << std::endl;
    if(e != 1){
	e++;
	continue;
    }

    auto const &opdets_handle =
      ev.getValidHandle< std::vector<raw::OpDetWaveform> >("pmtdecoder:PMTChannels");

    auto opdets(*opdets_handle);
        
    for(auto opdet : opdets){    
 //     if(opdet.ChannelNumber() != 305) continue;
//	std::cout << "Channel Number: " << opdet.ChannelNumber() << endl;

      std::vector< std::pair<int, float> > thing = Derivative_at_tick(opdet);
      
      for(auto pr : thing){
	devfunc->Fill(pr.first,pr.second);
      }
      
      //std::cout << opdet.ChannelNumber() << " had delay " << chanID_to_delayNS[opdet.ChannelNumber()] << std::endl;
      
      //for(auto tick : opdet.Waveform()){
      std::cout << " N ticks " << opdet.Waveform().size() << std::endl;

      std::vector<float> num;
      std::vector<float> den;

      for(int tick = 0; tick < opdet.Waveform().size(); tick ++){
	
	wf->Fill(tick, opdet.Waveform().at(tick));
	h_waveform_heatmap->Fill(tick,opdet.Waveform().at(tick));
			 
	num.push_back(opdet.Waveform().at(tick));
	den.push_back(opdet.Waveform().at(tick));

	if(num.size() > 5)
	  num.erase(num.begin());

	if(den.size() > 10)
	  den.erase(den.begin());
	
	if(num.size() == 5 && den.size() == 10 && std::accumulate(den.begin(), den.end(), 0) != 0)
	  dev->Fill(tick-5, float(std::accumulate(num.begin(), num.end(), 0))/float(std::accumulate(den.begin(), den.end(), 0)));
	
	//std::cout << opdet.Waveform().at(tick) << std::endl;
	//Assigning the values for the tree:
/*	EventNum=e;
	OpDetChannel=opdet.ChannelNumber();
	Tick=tick;
//	WaveformValue=opdet.Waveform().at(tick);
	Derivative=float(std::accumulate(num.begin(), num.end(), 0))/float(std::accumulate(den.begin(), den.end(), 0));
	ns= 2*tick; //2*tick
	calib_ns= ns + chanID_to_delayNS[opdet.ChannelNumber()]; //ns (see above) + map from above
	fTree->Fill();  */   
      }
    }//Raws	  	  
    e++;
  }//events
  
  TCanvas* c1 = new TCanvas("c1");
  c1->cd();
  wf->Draw("hist");


  TCanvas* c2 = new TCanvas("c2");
  c2->cd();
  dev->Draw("hist");

  TCanvas* c3 = new TCanvas("c3");
  c3->cd();
  devfunc->Draw("hist");

  TCanvas* c4 = new TCanvas("c4");
  c4->cd();
  h_waveform_heatmap->SetStats(0);
  h_waveform_heatmap->Draw("colz");


    out->cd();  
 // fTree->Write();      
}//Close Program


std::vector< std::pair<int, float> > Derivative_at_tick(raw::OpDetWaveform opt){
  
  auto smallestIt = std::min_element(opt.Waveform().begin(), opt.Waveform().end());
  auto smallestIndex = std::distance(opt.Waveform().begin(), smallestIt);

  std::vector<float> num;
  std::vector<float> den;
  
  std::vector< std::pair <int, float> > out; 
  
  for(int tick = smallestIndex-50; tick < smallestIndex+30; tick ++){
    
    num.push_back(opt.Waveform().at(tick));
    den.push_back(opt.Waveform().at(tick));
    
    if(num.size() > 5)
      num.erase(num.begin());
    
    if(den.size() > 10)
      den.erase(den.begin());
    
    if(num.size() == 5 && den.size() == 10 && std::accumulate(den.begin(), den.end(), 0) != 0){
      std::pair<int, float> pr = std::make_pair(tick, float(std::accumulate(num.begin(), num.end(), 0))/float(std::accumulate(den.begin(), den.end(), 0)));
      out.push_back(pr);
    }

  }
  
  return out;

}
