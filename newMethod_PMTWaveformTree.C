#define PMTWaveformTree_cxx
#include "PMTWaveformTree.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TGraph.h>
#include <TGraphErrors.h>


/*std::map<int, float> chanID_to_delayNS = {
//  {17,-3.75},
//  {305,0.8},
//  {63,0},
//  {297,4.8},
//  {16,1.5},
//  {36,0},
//  {294,-0.35},
//  {274,3.65},
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
};*/

std::map<int,float> chanID_to_delayNS = 
{{250,5.01},
{304,0.79},
{274,3.67},
{248,-0.18},
{302,5.11},
{246,2.75},
{300,3.2},
{272,2.32},
{244,2.79},
{298,3.47},
{242,3.37},
{296,-0.43},
{270,0.33},
{240,1.37},
{294,-0.33},
{172,1.52},
{226,1.15},
{196,-0.23},
{170,1.47},
{224,3.34},
{168,-1.07},
//{222,-1.33},
{194,1.66},
{166,1.87},
{220,3.42},
{164,0.99},
{218,-0.78},
{192,-0.27},
{162,-1.03},
{216,5.18},
{94,1.32},
{148,-1.36},
{118,-2.33},
{92,0.77},
{146,0.2},
{90,4.72},
{144,2.82},
{116,-3.1},
{88,1.7},
{142,-2.43},
//{86,-0.43},
{140,-0.48},
{114,0.09},
{84,-1.88},
//{138,1.2},
{16,1.54},
{70,5.24},
{40,5.02},
{14,3.27},
{68,5.19},
{12,0.57},
//{66,0.22},
{38,1.5},
{10,-2.68},
{64,5.2},
{8,-10.65},
{62,-4.48},
{36,0.02},
{6,-3.42},
{60,-4.28},
//{71,0.07},
{17,-3.75},
{41,-2.17},
{69,1.95},
{15,-0.13},
{67,1.57},
{13,-3.17},
{39,1.77},
{65,-1.58},
{11,5.19},
{63,0},
{9,-4.69},
{37,-4.26},
{61,-1.96},
{7,-3.5},
{149,1.97},
{95,-1.43},
{119,3.77},
{147,0.97},
{93,-1.13},
{145,3.62},
{91,3.64},
{117,1.67},
{143,1.17},
{89,3.38},
//{141,-3.48},
//{87,-0.35},
//{115,-1.11},
{139,-1.61},
//{85,-11.84},
{227,3.73},
{173,1.28},
{197,2.85},
{225,-2.69},
{171,-2.88},
{223,1.48},
{169,-0.01},
{195,3.61},
//{221,0.82},
//{167,0.64},
{219,-1.68},
{165,1.57},
{193,-2.13},
{217,0.84},
{163,0.97},
{305,0.84},
{251,4.46}, //checking this one
{275,1.41},
{303,4.75},
{249,-2.2},
{301,0.76},
{247,2.66},
{273,-0.64},
{299,3.2},
{245,1.32},
{297,4.82},
{243,4.85},
{271,0.32},
{295,-3.41},
{241,0.81}};

struct Point3D {
    double x, y, z;
};

std::map<int, Point3D> chanID_to_Point3D = {
 {6 , { -213.4, -175, 27.8742}},
 {7 , { 213.4, -175, 27.8742}},
 {8 , { -213.4, -95, 27.8742}},
 {9 , { 213.4, -95, 27.8742}},
 {10 , { -213.4, -40, 27.8742}},
 {11 , { 213.4, -40, 27.8742}},
 {12 , { -213.4, 40, 27.8742}},
 {13 , { 213.4, 40, 27.8742}},
 {14 , { -213.4, 95, 27.8742}},
 {15 , { 213.4, 95, 27.8742}},
 {16 , { -213.4, 175, 27.8742}},
 {17 , { 213.4, 175, 27.8742}},
 {36 , { -213.4, -135, 57.8742}},
 {37 , { 213.4, -135, 57.8742}},
 {38 , { -213.4, 0, 57.8742}},
 {39 , { 213.4, 0, 57.8742}},
 {40 , { -213.4, 135, 57.8742}},
 {41 , { 213.4, 135, 57.8742}},
 {60 , { -213.4, -175, 87.8742}},
 {61 , { 213.4, -175, 87.8742}},
 {62 , { -213.4, -95, 87.8742}},
 {63 , { 213.4, -95, 87.8742}},
 {64 , { -213.4, -40, 87.8742}},
 {65 , { 213.4, -40, 87.8742}},
 {66 , { -213.4, 40, 87.8742}},
 {67 , { 213.4, 40, 87.8742}},
 {68 , { -213.4, 95, 87.8742}},
 {69 , { 213.4, 95, 87.8742}},
 {70 , { -213.4, 175, 87.8742}},
 {71 , { 213.4, 175, 87.8742}},
 {84 , { -213.4, -175, 161.158}},
 {85 , { 213.4, -175, 161.158}},
 {86 , { -213.4, -95, 161.158}},
 {87 , { 213.4, -95, 161.158}},
 {88 , { -213.4, -40, 161.158}},
 {89 , { 213.4, -40, 161.158}},
 {90 , { -213.4, 40, 161.158}},
 {91 , { 213.4, 40, 161.158}},
 {92 , { -213.4, 95, 161.158}},
 {93 , { 213.4, 95, 161.158}},
 {94 , { -213.4, 175, 161.158}},
 {95 , { 213.4, 175, 161.158}},
 {114 , { -213.4, -135, 191.158}},
 {115 , { 213.4, -135, 191.158}},
 {116 , { -213.4, 0, 191.158}},
 {117 , { 213.4, 0, 191.158}},
 {118 , { -213.4, 135, 191.158}},
 {119 , { 213.4, 135, 191.158}},
 {138 , { -213.4, -175, 221.158}},
 {139 , { 213.4, -175, 221.158}},
 {140 , { -213.4, -95, 221.158}},
 {141 , { 213.4, -95, 221.158}},
 {142 , { -213.4, -40, 221.158}},
 {143 , { 213.4, -40, 221.158}},
 {144 , { -213.4, 40, 221.158}},
 {145 , { 213.4, 40, 221.158}},
 {146 , { -213.4, 95, 221.158}},
 {147 , { 213.4, 95, 221.158}},
 {148 , { -213.4, 175, 221.158}},
 {149 , { 213.4, 175, 221.158}},
 {162 , { -213.4, -175, 288.242}},
 {163 , { 213.4, -175, 288.242}},
 {164 , { -213.4, -95, 288.242}},
 {165 , { 213.4, -95, 288.242}},
 {166 , { -213.4, -40, 288.242}},
 {167 , { 213.4, -40, 288.242}},
 {168 , { -213.4, 40, 288.242}},
 {169 , { 213.4, 40, 288.242}},
 {170 , { -213.4, 95, 288.242}},
 {171 , { 213.4, 95, 288.242}},
 {172 , { -213.4, 175, 288.242}},
 {173 , { 213.4, 175, 288.242}},
 {192 , { -213.4, -135, 318.242}},
 {193 , { 213.4, -135, 318.242}},
 {194 , { -213.4, 0, 318.242}},
 {195 , { 213.4, 0, 318.242}},
 {196 , { -213.4, 135, 318.242}},
 {197 , { 213.4, 135, 318.242}},
 {216 , { -213.4, -175, 348.242}},
 {217 , { 213.4, -175, 348.242}},
 {218 , { -213.4, -95, 348.242}},
 {219 , { 213.4, -95, 348.242}},
 {220 , { -213.4, -40, 348.242}},
 {221 , { 213.4, -40, 348.242}},
 {222 , { -213.4, 40, 348.242}},
 {223 , { 213.4, 40, 348.242}},
 {224 , { -213.4, 95, 348.242}},
 {225 , { 213.4, 95, 348.242}},
 {226 , { -213.4, 175, 348.242}},
 {227 , { 213.4, 175, 348.242}},
 {240 , { -213.4, -175, 421.526}},
 {241 , { 213.4, -175, 421.526}},
 {242 , { -213.4, -95, 421.526}},
 {243 , { 213.4, -95, 421.526}},
 {244 , { -213.4, -40, 421.526}},
 {245 , { 213.4, -40, 421.526}},
 {246 , { -213.4, 40, 421.526}},
 {247 , { 213.4, 40, 421.526}},
 {248 , { -213.4, 95, 421.526}},
 {249 , { 213.4, 95, 421.526}},
 {250 , { -213.4, 175, 421.526}},
 {251 , { 213.4, 175, 421.526}},
 {270 , { -213.4, -135, 451.526}},
 {271 , { 213.4, -135, 451.526}},
 {272 , { -213.4, 0, 451.526}},
 {273 , { 213.4, 0, 451.526}},
 {274 , { -213.4, 135, 451.526}},
 {275 , { 213.4, 135, 451.526}},
 {294 , { -213.4, -175, 481.526}},
 {295 , { 213.4, -175, 481.526}},
 {296 , { -213.4, -95, 481.526}},
 {297 , { 213.4, -95, 481.526}},
 {298 , { -213.4, -40, 481.526}},
 {299 , { 213.4, -40, 481.526}},
 {300 , { -213.4, 40, 481.526}},
 {301 , { 213.4, 40, 481.526}},
 {302 , { -213.4, 95, 481.526}},
 {303 , { 213.4, 95, 481.526}},
 {304 , { -213.4, 175, 481.526}},
 {305 , { 213.4, 175, 481.526}}};
 
// const int thisRun=12813;
 const int refChannel=9;
 



void PMTWaveformTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L PMTWaveformTree.C
//      root> PMTWaveformTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

	const int derivStart_low = 3500;
	const int derivStart_high = 4500;


	const int counter_index=0;
	const int mean_index=1;
	const int rms_index=2;

	  //Make the 1D histograms here too
	std::map<int, TH1F*> map_hists_time;
	std::map<int, TH1F*> map_hists_timeSinceAvg;
	std::map<int, TH2F*> map_hists_timeSinceAvg_VS_avgTime;
	std::map<int, TH1F*> map_hists_timeSinceRef;
	std::map<int, TH1F*> map_hists_timeHalfVsMin;
	std::map<int, std::array<double, 3>> map_channel_mean_RMS;
	std::map<int, std::array<double, 3>> map_channel_mean_RMS_truncated;
	for(const std::pair<int, float> & x : chanID_to_delayNS){
//	   cout << "Testing the loop! " << x.first << "\t" << x.second << endl;
		map_hists_time[x.first] = new TH1F(Form("h_time_Ch%d",x.first),Form("h_time_Ch%d",x.first),2000, 6800,8800);
		map_hists_timeSinceAvg[x.first] = new TH1F(Form("h_timeSinceAvg_Ch%d",x.first),Form("h_timeSinceAvg_Ch%d",x.first), 200, -50,50);
		map_hists_timeSinceAvg_VS_avgTime[x.first] = new TH2F(Form("h_timeSinceAvg_VS_avgTime_Ch%d",x.first),Form("h_timeSinceAvg_VS_avgTime_Ch%d",x.first),2000, 6800,8800, 200, -50,50);
//		map_hists_time[x.first] = new TH1F(Form("h_time_Ch%d",x.first),Form("h_time_Ch%d",x.first),100, 5000,15000);
		map_hists_timeSinceRef[x.first] = new TH1F(Form("h_timeSinceRef_Ch%d",x.first),Form("h_timeSinceRef_Ch%d",x.first),50, -50,50);
		map_hists_timeHalfVsMin[x.first] = new TH1F(Form("h_timeHalfVsMin_Ch%d",x.first),Form("h_timeHalfVsMin_Ch%d",x.first),25, 0,10);
		map_channel_mean_RMS[x.first] = {0,0,0};
		map_channel_mean_RMS_truncated[x.first] = {0,0,0};
	}
	
	TH1F* h_avgTime = new TH1F("h_avgTime","h_avgTime",2000, 6800,8800);


   Long64_t nentries = fChain->GetEntriesFast();

//Internal Variables
	 int wf_length = 100; //starting value -- it actually looks for the proper length later
	 int der_length = 80; //starting value -- it actually looks for the proper length later
	int numPMTs_z=15;
	int numPMTs_y=9;
	
	double der_threshold = 0.48;
	


   Long64_t nbytes = 0, nb = 0;
   int num_events=0;
    int num_channels=0;
     	int channel_index=0;


	 //Map and components for the map
	  	int minBin=0;
		int der_time_halfMin=0;
		std::vector<std::pair< int, int >> EventNums;
	  	std::pair< int, float > ChannelNum_and_Time;
	  	std::vector<std::pair<int, float>> All_ChannelNum_and_Time;
	  	std::map<std::pair< int, int >, std::vector<std::pair<int, float>> > Event_PMTtimes; //Map of event: <channel, time>'s
	  	
	  	std::pair< int, float > ChannelNum_and_PMTsignal;
	  	std::vector<std::pair<int, float>> All_ChannelNum_and_PMTsignal;
	  	std::map<std::pair< int, int >, std::vector<std::pair<int, float>> > Event_PMTsignal; //Map of event: <channel, time>'s	


//Output file and the stuff for it:
	 TFile* outfile = new TFile("output_waveforms_newMethod.root","RECREATE"); //output file

	//The mean and RMS plots of the PMT times
 	TH2F* h_PMTtimes_negX_mean = new TH2F("h_PMTtimes_negX_mean", "Mean Signal Delay - East PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_mean = new TH2F("h_PMTtimes_posX_mean", "Mean Signal Delay - West PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_negX_RMS = new TH2F("h_PMTtimes_negX_RMS", "Signal Delay RMS - East PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_RMS = new TH2F("h_PMTtimes_posX_RMS", "Signal Delay RMS - West PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_negX_counts = new TH2F("h_PMTtimes_negX_counts", "h_PMTtimes_negX_counts", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_counts = new TH2F("h_PMTtimes_posX_counts", "h_PMTtimes_posX_counts", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
	
	//The mean, sigma and chi2 of the guassian fit
 	TH2F* h_PMTtimes_negX_gaus_mean = new TH2F("h_PMTtimes_negX_gaus_mean", "h_PMTtimes_negX_gaus_mean", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_gaus_mean = new TH2F("h_PMTtimes_posX_gaus_mean", "h_PMTtimes_posX_gaus_mean", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 		TH2F* h_PMTtimes_negX_gaus_mean_uncert = new TH2F("h_PMTtimes_negX_gaus_mean_uncert", "h_PMTtimes_negX_gaus_mean_uncert", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 		TH2F* h_PMTtimes_posX_gaus_mean_uncert = new TH2F("h_PMTtimes_posX_gaus_mean_uncert", "h_PMTtimes_posX_gaus_mean_uncert", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_negX_gaus_sigma = new TH2F("h_PMTtimes_negX_gaus_sigma", "h_PMTtimes_negX_gaus_sigma", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_gaus_sigma = new TH2F("h_PMTtimes_posX_gaus_sigma", "h_PMTtimes_posX_gaus_sigma", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
	 	TH2F* h_PMTtimes_negX_gaus_sigma_uncert = new TH2F("h_PMTtimes_negX_gaus_sigma_uncert", "h_PMTtimes_negX_gaus_sigma_uncert", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 		TH2F* h_PMTtimes_posX_gaus_sigma_uncert = new TH2F("h_PMTtimes_posX_gaus_sigma_uncert", "h_PMTtimes_posX_gaus_sigma_uncert", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_negX_gaus_chi2 = new TH2F("h_PMTtimes_negX_gaus_chi2", "h_PMTtimes_negX_gaus_chi2", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_gaus_chi2 = new TH2F("h_PMTtimes_posX_gaus_chi2", "h_PMTtimes_posX_gaus_chi2", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
	
 	TH2F* h_PMTtimes_negX_gaus_mean_minusCh9 = new TH2F("h_PMTtimes_negX_gaus_mean_minusCh9", "h_PMTtimes_negX_gaus_mean_minusCh9", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
 	TH2F* h_PMTtimes_posX_gaus_mean_minusCh9 = new TH2F("h_PMTtimes_posX_gaus_mean_minusCh9", "h_PMTtimes_posX_gaus_mean_minusCh9", numPMTs_z, 20, 490, numPMTs_y, -200, 200);

	
 	int prevEntry_runNum=-1; 	
 	int prevEntry_eventNum=-1;
 	int run_counter=1;
   
	 
	 //go through the events and determine the number of events and channels
		//Loop through the events   
		   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		      Long64_t ientry = LoadTree(jentry);
		      if (ientry < 0) break;
		      nb = fChain->GetEntry(jentry);   nbytes += nb;
			wf_length=Waveform->size();
			der_length=Derivative->size();
			
			std::pair< int, int > Run_Event;

		//Moved to a new event, need to save the completed event map and reset the vector for the next event
		if(prevEntry_eventNum!=Event && prevEntry_eventNum!=-1){
			  Run_Event = std::make_pair(prevEntry_runNum, prevEntry_eventNum);
//			cout << "Run: " << Run << "\tEvent Num: " << Event << endl;
			EventNums.push_back(Run_Event);
		  	Event_PMTtimes.insert({Run_Event, All_ChannelNum_and_Time});
		 	Event_PMTsignal.insert({Run_Event, All_ChannelNum_and_PMTsignal});
			 All_ChannelNum_and_Time.clear();
			 All_ChannelNum_and_PMTsignal.clear();
		}

	

   //Collecting and reorganizing the channels and times of each of the events:
            //Calculate the minimum tick for the time calculation
              int der_time_min = std::distance(Derivative->begin(), (std::min_element(Derivative->begin(), Derivative->end())));
              der_time_halfMin = bin_of_halfmin(*Derivative);

		
		//Make the new pair and add it to the vector (new or existing)
		  ChannelNum_and_Time = std::make_pair(Channel, 2.*(StartTick_Derivative+der_time_halfMin)+chanID_to_delayNS[Channel]);
//		  if(StartTick_Derivative < 2000 && StartTick_Derivative > 800) ChannelNum_and_PMTsignal = std::make_pair(Channel, Derivative->at(der_time_min));
		  if(StartTick_Derivative < derivStart_high && StartTick_Derivative > derivStart_low){
		  	ChannelNum_and_PMTsignal = std::make_pair(Channel, Derivative->at(der_time_min));
				if((Derivative->at(der_time_min))<der_threshold){
					 auto it = map_hists_timeHalfVsMin.find(Channel); //doin' the thing here!
					if(it != map_hists_timeHalfVsMin.end()){
						 map_hists_time[Channel]->Fill(2.*(StartTick_Derivative+der_time_halfMin)+chanID_to_delayNS[Channel]);

				 	}	              
				}
		  }
		  else  ChannelNum_and_PMTsignal = std::make_pair(Channel, 2);
//		  ChannelNum_and_PMTsignal = std::make_pair(Channel, Derivative->at(der_time_min));
		  All_ChannelNum_and_Time.push_back(ChannelNum_and_Time);
		  All_ChannelNum_and_PMTsignal.push_back(ChannelNum_and_PMTsignal);

//		cout << "2" << endl;
           

		 if(jentry==nentries-1){  //Add the last event info when you're at the last entry
			Run_Event = std::make_pair(prevEntry_runNum, prevEntry_eventNum);
		 	Event_PMTtimes.insert({Run_Event, All_ChannelNum_and_Time});
		 	Event_PMTsignal.insert({Run_Event, All_ChannelNum_and_PMTsignal});
			EventNums.push_back(Run_Event);
		 }
		 
		 auto it = map_hists_timeHalfVsMin.find(Channel); //doin' the thing here!
			if(it != map_hists_timeHalfVsMin.end()){
				 map_hists_timeHalfVsMin[Channel]->Fill(der_time_min-der_time_halfMin);
		 	}
		 
		 prevEntry_eventNum=Event;
		 prevEntry_runNum=Run;
		 
		 
//		 		cout << "3" << endl;
		 
			     
		   }
	
	cout << "Done with pass 1!!" << endl;

	
	
	//Map of event to average time

	std::map<std::pair< int, int >,double> event_avgTime;
	
	for(int i_event=0; i_event<EventNums.size(); i_event++){
		double avgTime=0;
		int Run = EventNums.at(i_event).first;
		int Event = EventNums.at(i_event).second;
		std::pair< int, int > Run_Event = std::make_pair(Run, Event);
		num_channels=Event_PMTtimes[EventNums.at(i_event)].size(); 
		int num_channelsWSignals=0;
		
		
		for(int i_channel=0; i_channel<num_channels; i_channel++){
			int chanNum=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).first;
			float time=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).second;
			if((Event_PMTsignal[EventNums.at(i_event)]).at(i_channel).second < der_threshold){ //passed the threshold

//				if((chanID_to_Point3D[chanNum].z > 250 || chanID_to_Point3D[chanNum].y>0) && chanID_to_Point3D[chanNum].x<0) continue; //To select only the Matt PMTs
//				if(chanID_to_Point3D[chanNum].x < 0) continue; //To only use the positive PMTs in the average time

				avgTime+=double(time);
				num_channelsWSignals+=1;
				}

		}
		
		if(num_channelsWSignals==0){
			event_avgTime.insert({Run_Event, 0});
			continue;
		}
		avgTime = avgTime/num_channelsWSignals;
		event_avgTime.insert({Run_Event, avgTime});
		h_avgTime->Fill(avgTime);
		
//		cout << "Run: " << Run << "\tEvent: " << Event << "\tNum Channels: " << num_channelsWSignals << "\tAvg Time: " << avgTime << endl;
	
	}
	
	outfile->cd();
	h_avgTime->GetXaxis()->SetTitle("Average event time [ns]");
	h_avgTime->GetYaxis()->SetTitle("Number of Events");
	h_avgTime->Write();
//	TCanvas* c00 = new TCanvas("c00", "c00");
//	c00->cd();
//	h_avgTime->Draw();
	
	cout << "Done with averaging the events' times, time to compare to the average for each channel!!" << endl;	

//TIME TO COMPARE THE INDIVIDUAL CHANNELS TO THE AVERAGE FOR THAT EVENT
	for(int i_event=0; i_event<EventNums.size(); i_event++){
		int Run = EventNums.at(i_event).first;
		int Event = EventNums.at(i_event).second;
		
		num_channels=Event_PMTtimes[EventNums.at(i_event)].size();
		
		for(int i_channel=0; i_channel<num_channels; i_channel++){
			int chanNum=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).first;
			float time=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).second;
			std::pair< int, int > Run_Event = std::make_pair(Run, Event);
			double runEvent_avgTime = event_avgTime[Run_Event];
			
			 auto it = map_hists_timeSinceAvg.find(chanNum); //checking that the channel exists
			if(it == map_hists_timeSinceAvg.end()) continue;
			
			if(runEvent_avgTime==0) continue;
			if(runEvent_avgTime>=7760) continue; //Here is the additional cut to select the blobs
			
			if((Event_PMTsignal[EventNums.at(i_event)]).at(i_channel).second < der_threshold){//passed the threshold

//				if((chanID_to_Point3D[chanNum].z > 250 || chanID_to_Point3D[chanNum].y>0) && chanID_to_Point3D[chanNum].x<0) continue; //To select only the Matt PMTs
				
				map_hists_timeSinceAvg[chanNum]->Fill(time-runEvent_avgTime);
				map_hists_timeSinceAvg_VS_avgTime[chanNum]->Fill(runEvent_avgTime, time-runEvent_avgTime);

					 
			 	map_channel_mean_RMS[chanNum][counter_index] +=1; //counting the number of signals for this PMT
			 	map_channel_mean_RMS[chanNum][mean_index] += time - runEvent_avgTime; //summing up the PMT signal times relative to the reference PMT (for the mean)
			 	map_channel_mean_RMS[chanNum][rms_index] += pow((time - runEvent_avgTime),2); //summing up the square of the PMT signal times relative to the reference PMT (for the RMS)
			}
		}
	}
	
	for(const std::pair<int, float> & x : chanID_to_delayNS){
		map_channel_mean_RMS[x.first][mean_index] = map_channel_mean_RMS[x.first][mean_index]/map_channel_mean_RMS[x.first][counter_index];
		map_channel_mean_RMS[x.first][rms_index] = sqrt(map_channel_mean_RMS[x.first][rms_index]/map_channel_mean_RMS[x.first][counter_index]);
		
		cout << "Channel  " << x.first << "\tMean  " << map_channel_mean_RMS[x.first][mean_index] << "\tRMS  " << map_channel_mean_RMS[x.first][rms_index] << endl;
	}
	
//TIME TO LOOP THROUGH THE EVENTS AGAIN FOR THE TRUNCATED MEAN AND RMS CALCULATION
	for(int i_event=0; i_event<EventNums.size(); i_event++){
		int Run = EventNums.at(i_event).first;
		int Event = EventNums.at(i_event).second;
		
		num_channels=Event_PMTtimes[EventNums.at(i_event)].size();
		
		for(int i_channel=0; i_channel<num_channels; i_channel++){
			int chanNum=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).first;
			float time=(Event_PMTtimes[EventNums.at(i_event)]).at(i_channel).second;
			std::pair< int, int > Run_Event = std::make_pair(Run, Event);
			double runEvent_avgTime = event_avgTime[Run_Event];
			
			 auto it = map_hists_timeSinceAvg.find(chanNum); //checking that the channel exists
			if(it == map_hists_timeSinceAvg.end()) continue;
			
			if(runEvent_avgTime==0) continue;
			if(runEvent_avgTime>=7760) continue; //Here is the additional cut to select the blobs
			
			if((Event_PMTsignal[EventNums.at(i_event)]).at(i_channel).second < der_threshold){ //passed the threshold
			
//				if((chanID_to_Point3D[chanNum].z > 250 || chanID_to_Point3D[chanNum].y>0) && chanID_to_Point3D[chanNum].x<0) continue; //To select only the Matt PMTs
			
				//if the time relative to the average is less than the mean - RMS continue
				if((time - runEvent_avgTime)<map_channel_mean_RMS[chanNum][mean_index] - map_channel_mean_RMS[chanNum][rms_index]) continue;
				//if the time relative to the average is more than the mean + RMS continue
				if((time - runEvent_avgTime)>map_channel_mean_RMS[chanNum][mean_index] + map_channel_mean_RMS[chanNum][rms_index]) continue;

				//if not, cacluate the the truncated mean and RMS
			 	map_channel_mean_RMS_truncated[chanNum][counter_index] +=1; //counting the number of signals for this PMT
			 	map_channel_mean_RMS_truncated[chanNum][mean_index] += time - runEvent_avgTime; //summing up the PMT signal times relative to the reference PMT (for the mean)
			 	map_channel_mean_RMS_truncated[chanNum][rms_index] += pow((time - runEvent_avgTime),2); //summing up the square of the PMT signal times relative to the reference PMT (for the RMS)
				
			}
		}
	}

	for(const std::pair<int, float> & x : chanID_to_delayNS){
		map_channel_mean_RMS_truncated[x.first][mean_index] = map_channel_mean_RMS_truncated[x.first][mean_index]/map_channel_mean_RMS_truncated[x.first][counter_index];
		map_channel_mean_RMS_truncated[x.first][rms_index] = sqrt(map_channel_mean_RMS_truncated[x.first][rms_index]/map_channel_mean_RMS_truncated[x.first][counter_index]);
		
//		cout << "TRUNCATED Channel  " << x.first << "\tMean  " << map_channel_mean_RMS_truncated[x.first][mean_index] << "\tRMS  " << map_channel_mean_RMS_truncated[x.first][rms_index] << endl;
	}


	//location of the origin of the spark to determine the radial distance for each PMT
/*	const float origin_x = 213;
	const float origin_y = -213;
	const float origin_z = 0;*/
	
	//Best fit point
	const float origin_x = 162;
	const float origin_y = -213;
	const float origin_z = 22;
	
	//vectors for the r vs t_avg plots
       std::vector<Double_t> v_radius_posx;
       std::vector<Double_t> v_radius_err_posx;
       std::vector<Double_t> v_t_bar_posx;
       std::vector<Double_t> v_t_bar_err_posx;
       std::vector<int> v_channel_posx;

       std::vector<Double_t> v_radius_negx;
       std::vector<Double_t> v_radius_err_negx;
       std::vector<Double_t> v_t_bar_negx;
       std::vector<Double_t> v_t_bar_err_negx;
       std::vector<int> v_channel_negx;
       
	
	for(const std::pair<int, TH1F*> & x : map_hists_timeSinceAvg){
		//Fit them
		  TF1 *f_gaus = new TF1("f_gaus", "gaus", -30, 30); // Gaussian function
		 map_hists_timeSinceAvg[x.first]->Fit("f_gaus", "Q");

		double mean = f_gaus->GetParameter(1);
		double mean_uncert = f_gaus->GetParError(1);
		double sigma = f_gaus->GetParameter(2);
		double sigma_uncert = f_gaus->GetParError(2);
		double chiSquare = f_gaus->GetChisquare();

		int bin_z = h_PMTtimes_posX_mean->GetXaxis()->FindBin(chanID_to_Point3D[x.first].z);
		int bin_y = h_PMTtimes_posX_mean->GetYaxis()->FindBin(chanID_to_Point3D[x.first].y);

		double radius = sqrt(pow(chanID_to_Point3D[x.first].x-origin_x,2)+pow(chanID_to_Point3D[x.first].y-origin_y,2)+pow(chanID_to_Point3D[x.first].z-origin_z,2)); //Where we calculate the radius from the "origin" of our choosing (spark for example)
		double radius_err = 10.16;

		
		if(chanID_to_Point3D[x.first].x > 0){
//				if((chanID_to_Point3D[x.first].z > 250 || chanID_to_Point3D[x.first].y>0) && chanID_to_Point3D[x.first].x<0) continue; //To select only the Matt PMTs
			h_PMTtimes_posX_gaus_mean->SetBinContent(bin_z, bin_y, mean);
			h_PMTtimes_posX_gaus_mean_uncert->SetBinContent(bin_z, bin_y, mean_uncert);
			h_PMTtimes_posX_gaus_sigma->SetBinContent(bin_z, bin_y, sigma);
			h_PMTtimes_posX_gaus_sigma_uncert->SetBinContent(bin_z, bin_y, sigma_uncert);
			h_PMTtimes_posX_gaus_chi2->SetBinContent(bin_z, bin_y, chiSquare);

			h_PMTtimes_posX_counts->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][counter_index]);
			h_PMTtimes_posX_mean->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][mean_index]);
			h_PMTtimes_posX_RMS->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][rms_index]);
			

			//add to the vector for the r vs t_avg plots
			v_channel_posx.push_back(x.first);
			v_t_bar_posx.push_back(map_channel_mean_RMS_truncated[x.first][mean_index]);
			v_t_bar_err_posx.push_back(sigma);
			v_radius_posx.push_back(radius);
			v_radius_err_posx.push_back(radius_err);

			
			if(map_channel_mean_RMS_truncated[x.first][rms_index] !=0) cout << "{" << x.first << ", {" << chanID_to_Point3D[x.first].x << "," << chanID_to_Point3D[x.first].y << "," << chanID_to_Point3D[x.first].z << "}}" << endl;
			
		}
		if(chanID_to_Point3D[x.first].x < 0){
//				if((chanID_to_Point3D[x.first].z > 250 || chanID_to_Point3D[x.first].y>0) && chanID_to_Point3D[x.first].x<0) continue; //To select only the Matt PMTs
			h_PMTtimes_negX_gaus_mean->SetBinContent(bin_z, bin_y, mean);
			h_PMTtimes_negX_gaus_mean_uncert->SetBinContent(bin_z, bin_y, mean_uncert);
			h_PMTtimes_negX_gaus_sigma->SetBinContent(bin_z, bin_y, sigma);
			h_PMTtimes_negX_gaus_sigma_uncert->SetBinContent(bin_z, bin_y, sigma_uncert);
			h_PMTtimes_negX_gaus_chi2->SetBinContent(bin_z, bin_y, chiSquare);

			h_PMTtimes_negX_counts->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][counter_index]);
			h_PMTtimes_negX_mean->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][mean_index]);
			h_PMTtimes_negX_RMS->SetBinContent(bin_z, bin_y, map_channel_mean_RMS_truncated[x.first][rms_index]);

			if(map_channel_mean_RMS_truncated[x.first][mean_index]>18) cout << "THIS IS THE CHANNEL: " << x.first << " !!!!!!!!!!!" << " Radius: " << radius << " time delay: " << map_channel_mean_RMS_truncated[x.first][mean_index] << endl;

			//add to the vector for the r vs t_avg plots
			v_channel_negx.push_back(x.first);
			v_t_bar_negx.push_back(map_channel_mean_RMS_truncated[x.first][mean_index]);
			v_t_bar_err_negx.push_back(sigma);
			v_radius_negx.push_back(radius);
			v_radius_err_negx.push_back(radius_err);

		}
		
		//Write them
		outfile->cd();
		map_hists_time[x.first]->GetXaxis()->SetTitle("Time [ns]");
		map_hists_time[x.first]->GetYaxis()->SetTitle("Counts");
		map_hists_time[x.first]->Write();

		map_hists_timeSinceAvg[x.first]->GetXaxis()->SetTitle("Time relative to average [ns]");
		map_hists_timeSinceAvg[x.first]->GetYaxis()->SetTitle("Counts");
		map_hists_timeSinceAvg[x.first]->Write();
		
		map_hists_timeSinceAvg_VS_avgTime[x.first]->GetXaxis()->SetTitle("Average Time [ns]");
		map_hists_timeSinceAvg_VS_avgTime[x.first]->GetYaxis()->SetTitle("Time relative to average [ns]");
		map_hists_timeSinceAvg_VS_avgTime[x.first]->GetZaxis()->SetTitle("Counts");
		map_hists_timeSinceAvg_VS_avgTime[x.first]->Write();
	}
	
      outfile->cd();
       h_PMTtimes_negX_gaus_mean->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_negX_gaus_mean->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_negX_gaus_mean->GetZaxis()->SetTitle("ns");
       h_PMTtimes_negX_gaus_mean->SetStats(0);
       h_PMTtimes_negX_gaus_mean->Write();
       h_PMTtimes_posX_gaus_mean->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_posX_gaus_mean->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_posX_gaus_mean->GetZaxis()->SetTitle("ns");
       h_PMTtimes_posX_gaus_mean->SetStats(0);
       h_PMTtimes_posX_gaus_mean->Write();
	       h_PMTtimes_negX_gaus_mean_uncert->GetXaxis()->SetTitle("Z[cm]");
	       h_PMTtimes_negX_gaus_mean_uncert->GetYaxis()->SetTitle("Y[cm]");
	       h_PMTtimes_negX_gaus_mean_uncert->GetZaxis()->SetTitle("Uncertainty in ns after Average");
	       h_PMTtimes_negX_gaus_mean_uncert->SetStats(0);
	       h_PMTtimes_negX_gaus_mean_uncert->Write();
	       h_PMTtimes_posX_gaus_mean_uncert->GetXaxis()->SetTitle("Z[cm]");
	       h_PMTtimes_posX_gaus_mean_uncert->GetYaxis()->SetTitle("Y[cm]");
	       h_PMTtimes_posX_gaus_mean_uncert->GetZaxis()->SetTitle("Uncertainty in ns after Average");
	       h_PMTtimes_posX_gaus_mean_uncert->SetStats(0);
	       h_PMTtimes_posX_gaus_mean_uncert->Write();
	       
/*		TCanvas* c8 = new TCanvas("c8","c8");
		c8->cd();
		h_PMTtimes_posX_gaus_mean->Draw("colz");
*/

       h_PMTtimes_negX_gaus_sigma->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_negX_gaus_sigma->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_negX_gaus_sigma->GetZaxis()->SetTitle("Sigma");
       h_PMTtimes_negX_gaus_sigma->SetStats(0);
       h_PMTtimes_negX_gaus_sigma->Write();
       h_PMTtimes_posX_gaus_sigma->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_posX_gaus_sigma->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_posX_gaus_sigma->GetZaxis()->SetTitle("Sigma");
       h_PMTtimes_posX_gaus_sigma->SetStats(0);
       h_PMTtimes_posX_gaus_sigma->Write();
	       h_PMTtimes_negX_gaus_sigma_uncert->GetXaxis()->SetTitle("Z[cm]");
	       h_PMTtimes_negX_gaus_sigma_uncert->GetYaxis()->SetTitle("Y[cm]");
	       h_PMTtimes_negX_gaus_sigma_uncert->GetZaxis()->SetTitle("Sigma Uncertainty");
	       h_PMTtimes_negX_gaus_sigma_uncert->SetStats(0);
	       h_PMTtimes_negX_gaus_sigma_uncert->Write();
	       h_PMTtimes_posX_gaus_sigma_uncert->GetXaxis()->SetTitle("Z[cm]");
	       h_PMTtimes_posX_gaus_sigma_uncert->GetYaxis()->SetTitle("Y[cm]");
	       h_PMTtimes_posX_gaus_sigma_uncert->GetZaxis()->SetTitle("Sigma Uncertainty");
	       h_PMTtimes_posX_gaus_sigma_uncert->SetStats(0);
	       h_PMTtimes_posX_gaus_sigma_uncert->Write();

       h_PMTtimes_negX_gaus_chi2->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_negX_gaus_chi2->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_negX_gaus_chi2->GetZaxis()->SetTitle("#Chi^{2}");
       h_PMTtimes_negX_gaus_chi2->SetStats(0);
       h_PMTtimes_negX_gaus_chi2->Write();
       h_PMTtimes_posX_gaus_chi2->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_posX_gaus_chi2->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_posX_gaus_chi2->GetZaxis()->SetTitle("#Chi^{2}");
       h_PMTtimes_posX_gaus_chi2->SetStats(0);
       h_PMTtimes_posX_gaus_chi2->Write();

       outfile->cd();
       h_PMTtimes_negX_mean->GetXaxis()->SetTitle("Z [cm]");
       h_PMTtimes_negX_mean->GetYaxis()->SetTitle("Y [cm]");
//       h_PMTtimes_negX_mean->GetZaxis()->SetRangeUser(-5,30);
       h_PMTtimes_negX_mean->GetZaxis()->SetTitle("Delay from Average [ns]");
       h_PMTtimes_negX_mean->SetStats(0);
       h_PMTtimes_negX_mean->Write();
       h_PMTtimes_posX_mean->GetXaxis()->SetTitle("Z [cm]");
//       h_PMTtimes_posX_mean->GetZaxis()->SetRangeUser(-5,20);
       h_PMTtimes_posX_mean->GetYaxis()->SetTitle("Y [cm]");
       h_PMTtimes_posX_mean->GetZaxis()->SetTitle("Delay from Average [ns]");
       h_PMTtimes_posX_mean->SetStats(0);
       h_PMTtimes_posX_mean->Write();

/*		TText *text = new TText(0, 3, "SBND Preliminary");  // (x, y, text)
//		text->SetTextSize(0.05);
		text->SetTextColor(kBlack);

	   TCanvas *anotherCanvas = new TCanvas("anotherCanvas", "anotherCanvas", 800, 600);
		h_PMTtimes_posX_mean->Draw("COLZ");
		text->Draw("SAME");

		TCanvas *anotherCanvas1 = new TCanvas("anotherCanvas1", "anotherCanvas1", 800, 600);
		h_PMTtimes_negX_mean->Draw("COLZ");
		text->Draw("SAME");*/

       h_PMTtimes_negX_RMS->GetXaxis()->SetTitle("Z [cm]");
       h_PMTtimes_negX_RMS->GetYaxis()->SetTitle("Y [cm]");
       h_PMTtimes_negX_RMS->GetZaxis()->SetTitle("RMS");
       h_PMTtimes_negX_RMS->SetStats(0);
       h_PMTtimes_negX_RMS->Write();
       h_PMTtimes_posX_RMS->GetXaxis()->SetTitle("Z [cm]");
       h_PMTtimes_posX_RMS->GetYaxis()->SetTitle("Y [cm]");
       h_PMTtimes_posX_RMS->GetZaxis()->SetTitle("RMS");
       h_PMTtimes_posX_RMS->SetStats(0);
       h_PMTtimes_posX_RMS->Write();

//		TText *text = new TText(0, 3, "SBND Preliminary");  // (x, y, text)
//		text->SetTextSize(0.05);
//		text->SetTextColor(kBlack);

	   TCanvas *anotherCanvas = new TCanvas("anotherCanvas", "anotherCanvas", 800, 600);
		h_PMTtimes_posX_RMS->Draw("COLZ");
//		text->Draw("SAME");

		TCanvas *anotherCanvas1 = new TCanvas("anotherCanvas1", "anotherCanvas1", 800, 600);
		h_PMTtimes_negX_RMS->Draw("COLZ");
//		text->Draw("SAME");

       h_PMTtimes_negX_counts->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_negX_counts->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_negX_counts->GetZaxis()->SetTitle("Counts");
       h_PMTtimes_negX_counts->SetStats(0);
       h_PMTtimes_negX_counts->Write();
       h_PMTtimes_posX_counts->GetXaxis()->SetTitle("Z[cm]");
       h_PMTtimes_posX_counts->GetYaxis()->SetTitle("Y[cm]");
       h_PMTtimes_posX_counts->GetZaxis()->SetTitle("Counts");
       h_PMTtimes_posX_counts->SetStats(0);
       h_PMTtimes_posX_counts->Write();
       
       
       // Change the vectors into arrays and make the TGraphs
              Double_t *arr_radius_posx = v_radius_posx.data();
              Double_t *arr_radius_err_posx = v_radius_err_posx.data();
              Double_t *arr_t_bar_posx = v_t_bar_posx.data();
              Double_t *arr_t_bar_err_posx = v_t_bar_err_posx.data();

              Double_t *arr_radius_negx = v_radius_negx.data();
              Double_t *arr_radius_err_negx = v_radius_err_negx.data();
              Double_t *arr_t_bar_negx = v_t_bar_negx.data();
              Double_t *arr_t_bar_err_negx = v_t_bar_err_negx.data();
              
	TF1 *f_line_24cmperns_posx = new TF1("f_line_24cmperns_posx", "[0]+x/24", 0,750); // Linear function
	f_line_24cmperns_posx->SetLineColor(kBlue);
              
	TF1 *f_line_posx = new TF1("f_line_posx", "[0]+[1]*x", 0,750); // Linear function
       TGraphErrors *gr_tbar_vs_radius_posx = new TGraphErrors(v_radius_posx.size(), arr_radius_posx, arr_t_bar_posx, arr_radius_err_posx, arr_t_bar_err_posx);
       gr_tbar_vs_radius_posx->Fit(f_line_posx, "R");
//       gr_tbar_vs_radius_posx->Fit(f_line_24cmperns_posx, "R");
//		cout << "velocity according to the positive PMTs = " << 1/(f_line_posx->GetParameter(1)) << "cm/ns" << endl;
       gr_tbar_vs_radius_posx->SetTitle("PMT Mean Delay Time - West PMTs");
       gr_tbar_vs_radius_posx->GetXaxis()->SetTitle("Radius w.r.t. Bottom Southwest Corner [cm]");
       gr_tbar_vs_radius_posx->GetYaxis()->SetTitle("PMT Mean Delay Time [ns]");
	gr_tbar_vs_radius_posx->SetMarkerStyle(20);
	gr_tbar_vs_radius_posx->SetMarkerColor(kBlack);
//	gr_tbar_vs_radius_posx->SetLineWidth(0);
	gr_tbar_vs_radius_posx->Write();
		TCanvas* c2 = new TCanvas("c2","c2");
		c2->cd();
		gr_tbar_vs_radius_posx->Draw("AP");
		f_line_posx->Draw("same");
//		f_line_24cmperns_posx->Draw("same");
//		text->Draw("SAME");



	TF1 *f_line_24cmperns_negx = new TF1("f_line_24cmperns_negx", "[0]+x/24", 0,750); // Linear function
	f_line_24cmperns_negx->SetLineColor(kBlue);

	TF1 *f_line_negx = new TF1("f_line_negx", "[0]+[1]*x", 0,750); // Linear function
       TGraphErrors *gr_tbar_vs_radius_negx = new TGraphErrors(v_radius_negx.size(), arr_radius_negx, arr_t_bar_negx, arr_radius_err_negx, arr_t_bar_err_negx);
//       gr_tbar_vs_radius_negx->Fit(f_line_negx, "R");
//       gr_tbar_vs_radius_negx->Fit(f_line_24cmperns_negx, "R");
//		cout << "velocity according to the negative PMTs = " << 1/(f_line_negx->GetParameter(1)) << "cm/ns" << endl;
       gr_tbar_vs_radius_negx->SetTitle("PMT Mean Delay Time - East PMTs");
       gr_tbar_vs_radius_negx->SetName("PMT Mean Delay Time - East PMTs");
       gr_tbar_vs_radius_negx->GetXaxis()->SetTitle("Radius w.r.t. Bottom Southwest Corner [cm]");
       gr_tbar_vs_radius_negx->GetYaxis()->SetTitle("PMT Mean Delay Time [ns]");
       gr_tbar_vs_radius_negx->GetXaxis()->SetTitleSize(0.07);
       gr_tbar_vs_radius_negx->GetYaxis()->SetTitleSize(0.07);
       gr_tbar_vs_radius_negx->GetXaxis()->SetLabelSize(0.07);
       gr_tbar_vs_radius_negx->GetYaxis()->SetLabelSize(0.07);
	gr_tbar_vs_radius_negx->SetMarkerStyle(20);
	gr_tbar_vs_radius_negx->SetMarkerColor(kBlack);
//	gr_tbar_vs_radius_negx->SetLineWidth(0);
//		TText *text2 = new TText(500, 0, "SBND Preliminary");  // (x, y, text)
//		text->SetTextSize(0.05);
//		text2->SetTextColor(kBlack);
		TCanvas* c3 = new TCanvas("c3","c3");
		c3->cd();
		gr_tbar_vs_radius_negx->Draw("AP");
//		text2->Draw("SAME");

	gr_tbar_vs_radius_negx->Write();

		
//		f_line_negx->Draw("same");
//		f_line_24cmperns_negx->Draw("same");
/*
		TMultiGraph *mg = new TMultiGraph();
		mg->Add(gr_tbar_vs_radius_negx);
		mg->Add(gr_tbar_vs_radius_posx);
    	 	  mg->GetXaxis()->SetTitle("radius [cm]");
   		    mg->GetYaxis()->SetTitle("PMT mean delay time [ns]");

		TCanvas* c4 = new TCanvas("c4","c4");
		c4->cd();
		mg->Draw("AP");
		f_line_negx->Draw("same");
		f_line_posx->Draw("same");
		f_line_24cmperns_negx->Draw("same");

    Double_t chiSquare_posx = f_line_24cmperns_posx->GetChisquare();
    Int_t degreesOfFreedom_posx = v_radius_posx.size() - f_line_24cmperns_posx->GetNpar();
    Double_t reducedChiSquare_posx = chiSquare_posx / degreesOfFreedom_posx;

    Double_t chiSquare_negx = f_line_24cmperns_negx->GetChisquare();
    Int_t degreesOfFreedom_negx = v_radius_negx.size() - f_line_24cmperns_negx->GetNpar();
    Double_t reducedChiSquare_negx = chiSquare_negx / degreesOfFreedom_negx;

	cout << "Positive: Chi^2/NDF: " << chiSquare_posx << " / " << degreesOfFreedom_posx << " = " << reducedChiSquare_posx << endl;
	cout << "Negative: Chi^2/NDF: " << chiSquare_negx << " / " << degreesOfFreedom_negx << " = " << reducedChiSquare_negx << endl;
		
		
//Make the plots looking at the comparison between the data and the velocity lines
//Do something like this: double x = 3.14; double y; y = f->Eval(x); // y = f(3.14)

	//arrays for the r vs t_avg plots
       std::vector<Double_t> v_diff_tbar_fit_posx;
       std::vector<Double_t> v_diff_tbar_24_posx;
//       Double_t *arr_diff_tbar_err_posx[v_radius_posx.size()];

       std::vector<Double_t> v_diff_tbar_fit_negx;
       std::vector<Double_t> v_diff_tbar_24_negx;
 //      Double_t *arr_diff_tbar_err_negx[v_radius_negx.size()];



	TH1F* h_diffs_tbar_fit_posx = new TH1F("h_diffs_tbar_fit_posx","h_diffs_tbar_fit_posx",50, -20, 20);
	TH1F* h_diffs_tbar_24_posx = new TH1F("h_diffs_tbar_24_posx","h_diffs_tbar_24_posx",50, -20, 20);
	TH1F* h_diffs_tbar_fit_negx = new TH1F("h_diffs_tbar_fit_negx","h_diffs_tbar_fit_negx",50, -20, 20);
	TH1F* h_diffs_tbar_24_negx = new TH1F("h_diffs_tbar_24_negx","h_diffs_tbar_24_negx",50, -20, 20);

       
	//looping through the negative PMT channels
	for(int iEntry=0; iEntry<v_radius_negx.size(); iEntry++){
		//for each point, get the time delay and the line value and take the diff
		v_diff_tbar_fit_negx.push_back(arr_t_bar_negx[iEntry]-(f_line_negx->Eval(arr_radius_negx[iEntry])));
		v_diff_tbar_24_negx.push_back(arr_t_bar_negx[iEntry]-(f_line_24cmperns_negx->Eval(arr_radius_negx[iEntry])));

		h_diffs_tbar_fit_negx->Fill(arr_t_bar_negx[iEntry]-(f_line_negx->Eval(arr_radius_negx[iEntry])));
		h_diffs_tbar_24_negx->Fill(arr_t_bar_negx[iEntry]-(f_line_24cmperns_negx->Eval(arr_radius_negx[iEntry])));
	}

	for(int iEntry=0; iEntry<v_radius_posx.size(); iEntry++){
		//for each point, get the time delay and the line value and take the diff
		v_diff_tbar_fit_posx.push_back(arr_t_bar_posx[iEntry]-(f_line_posx->Eval(arr_radius_posx[iEntry])));
		v_diff_tbar_24_posx.push_back(arr_t_bar_posx[iEntry]-(f_line_24cmperns_posx->Eval(arr_radius_posx[iEntry])));

		h_diffs_tbar_fit_posx->Fill(arr_t_bar_posx[iEntry]-(f_line_posx->Eval(arr_radius_posx[iEntry])));
		h_diffs_tbar_24_posx->Fill(arr_t_bar_posx[iEntry]-(f_line_24cmperns_posx->Eval(arr_radius_posx[iEntry])));
	}
	
	//convert to arrays
        Double_t *arr_diff_tbar_fit_negx = v_diff_tbar_fit_negx.data();
        Double_t *arr_diff_tbar_24_negx = v_diff_tbar_24_negx.data();

        Double_t *arr_diff_tbar_fit_posx = v_diff_tbar_fit_posx.data();
        Double_t *arr_diff_tbar_24_posx = v_diff_tbar_24_posx.data();
	
	//Write/Draw them!!!
       TGraph *gr_diff_tbar_fit_negx = new TGraph(v_radius_negx.size(), arr_radius_negx, arr_diff_tbar_fit_negx);
       gr_diff_tbar_fit_negx->SetTitle("Negative TPC PMTs wrt Fit line");
       gr_diff_tbar_fit_negx->GetXaxis()->SetTitle("radius [cm]");
       gr_diff_tbar_fit_negx->GetYaxis()->SetTitle("DIFFERENCE: PMT mean delay time-Fit [ns}");
	gr_diff_tbar_fit_negx->SetMarkerStyle(20);
	gr_diff_tbar_fit_negx->SetMarkerColor(kBlack);
	gr_diff_tbar_fit_negx->SetLineWidth(0);
	gr_diff_tbar_fit_negx->Write();
		TCanvas* c5 = new TCanvas("c5","c5");
		c5->cd();
		gr_diff_tbar_fit_negx->Draw("AP");

       TGraph *gr_diff_tbar_fit_posx = new TGraph(v_radius_posx.size(), arr_radius_posx, arr_diff_tbar_fit_posx);
       gr_diff_tbar_fit_posx->SetTitle("Positive TPC PMTs wrt Fit line");
       gr_diff_tbar_fit_posx->GetXaxis()->SetTitle("radius [cm]");
       gr_diff_tbar_fit_posx->GetYaxis()->SetTitle("DIFFERENCE: PMT mean delay time-Fit [ns}");
	gr_diff_tbar_fit_posx->SetMarkerStyle(20);
	gr_diff_tbar_fit_posx->SetMarkerColor(kBlack);
	gr_diff_tbar_fit_posx->SetLineWidth(0);
	gr_diff_tbar_fit_posx->Write();
		TCanvas* c6 = new TCanvas("c6","c6");
		c6->cd();
		gr_diff_tbar_fit_posx->Draw("AP");



       TGraph *gr_diff_tbar_24_negx = new TGraph(v_radius_negx.size(), arr_radius_negx, arr_diff_tbar_24_negx);
       gr_diff_tbar_24_negx->SetTitle("Negative TPC PMTs wrt 24 cm/ns line");
       gr_diff_tbar_24_negx->GetXaxis()->SetTitle("radius [cm]");
       gr_diff_tbar_24_negx->GetYaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
	gr_diff_tbar_24_negx->SetMarkerStyle(20);
	gr_diff_tbar_24_negx->SetMarkerColor(kBlack);
	gr_diff_tbar_24_negx->SetLineWidth(0);
	gr_diff_tbar_24_negx->Write();
		TCanvas* c7 = new TCanvas("c7","c7");
		c7->cd();
		gr_diff_tbar_24_negx->Draw("AP");

       TGraph *gr_diff_tbar_24_posx = new TGraph(v_radius_posx.size(), arr_radius_posx, arr_diff_tbar_24_posx);
       gr_diff_tbar_24_posx->SetTitle("Positive TPC PMT wrt 24 cm/ns lines");
       gr_diff_tbar_24_posx->GetXaxis()->SetTitle("radius [cm]");
       gr_diff_tbar_24_posx->GetYaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
	gr_diff_tbar_24_posx->SetMarkerStyle(20);
	gr_diff_tbar_24_posx->SetMarkerColor(kBlack);
	gr_diff_tbar_24_posx->SetLineWidth(0);
	gr_diff_tbar_24_posx->Write();
		TCanvas* c8 = new TCanvas("c8","c8");
		c8->cd();
		gr_diff_tbar_24_posx->Draw("AP");

	  TF1 *f_gaus_diff = new TF1("f_gaus_diff", "gaus", -30, 30); // Gaussian function

//Plotting the histograms
		TCanvas* c9 = new TCanvas("c9","c9");
		c9->cd();
		h_diffs_tbar_fit_posx->GetXaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
		h_diffs_tbar_fit_posx->Draw("hist");

		TCanvas* c10 = new TCanvas("c10","c10");
		c10->cd();
		h_diffs_tbar_fit_negx->GetXaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
		h_diffs_tbar_fit_negx->Draw("hist");

		TCanvas* c11 = new TCanvas("c11","c11");
		c11->cd();
		h_diffs_tbar_24_posx->GetXaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
		 h_diffs_tbar_24_posx->Fit("f_gaus_diff", "R");
		h_diffs_tbar_24_posx->Draw("hist");
		f_gaus_diff->Draw("same");
		
		double diffs_tbar_24_posx_mean = f_gaus_diff->GetParameter(1);
		double diffs_tbar_24_posx_mean_uncert = f_gaus_diff->GetParError(1);
		double diffs_tbar_24_posx_sigma = f_gaus_diff->GetParameter(2);
		double diffs_tbar_24_posx_sigma_uncert = f_gaus_diff->GetParError(2);
		double diffs_tbar_24_posx_chi2 = f_gaus_diff->GetChisquare();
		
		TCanvas* c12 = new TCanvas("c12","c12");
		c12->cd();
		h_diffs_tbar_24_negx->GetXaxis()->SetTitle("DIFFERENCE: PMT mean delay time - line [ns}");
		 h_diffs_tbar_24_negx->Fit("f_gaus_diff", "R");
		h_diffs_tbar_24_negx->Draw("hist");
		f_gaus_diff->Draw("same");

		double diffs_tbar_24_negx_mean = f_gaus_diff->GetParameter(1);
		double diffs_tbar_24_negx_mean_uncert = f_gaus_diff->GetParError(1);
		double diffs_tbar_24_negx_sigma = f_gaus_diff->GetParameter(2);
		double diffs_tbar_24_negx_sigma_uncert = f_gaus_diff->GetParError(2);
		double diffs_tbar_24_negx_chi2 = f_gaus_diff->GetChisquare();
*/
	
	cout << "Have a nice day!!" << endl;
 //  outfile->Close();



}

