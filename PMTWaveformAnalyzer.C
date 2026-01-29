#define PMTWaveformTree_cxx
#include "PMTWaveformTree.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TSystem.h>
#include <iostream>
#include <numeric>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

// === HALF-MINIMUM FINDER ===
int bin_of_halfmin(const std::vector<float>& der) {
    if (der.empty()) return -1;
    auto min_it = std::min_element(der.begin(), der.end());
    float min_val = *min_it;
    if (min_val >= 1.0) return -1;
    float half_val = min_val + 0.5 * (1.0 - min_val);
    for (size_t i = 0; i < der.size(); ++i) {
        if (der[i] <= half_val) return i;
    }
    return -1;
}

// === Run-to-Diffuser Location Mapping ===
std::map<int, std::pair<int, std::string>> run_to_diffuser = {
    {19403, {1, "West"}}, {19406, {1, "West"}}, {19407, {1, "West"}}, {19408, {1, "West"}}, {19409, {1, "West"}},
    {19410, {2, "West"}}, {19411, {2, "West"}}, {19413, {2, "West"}}, {19414, {2, "West"}}, {19418, {2, "West"}},
    {19417, {3, "West"}}, {19419, {3, "West"}}, {19420, {3, "West"}}, {19421, {3, "West"}}, {19423, {3, "West"}},
    {19424, {4, "West"}}, {19425, {4, "West"}}, {19426, {4, "West"}}, {19427, {4, "West"}}, {19428, {4, "West"}},
    {19429, {5, "West"}}, {19430, {5, "West"}}, {19431, {5, "West"}}, {19432, {5, "West"}}, {19433, {5, "West"}},
    {19434, {6, "East"}}, {19435, {6, "East"}}, {19436, {6, "East"}}, {19437, {6, "East"}}, {19438, {6, "East"}},
    {19439, {7, "East"}}, {19440, {7, "East"}}, {19441, {7, "East"}}, {19442, {7, "East"}}, {19443, {7, "East"}},
    {19449, {8, "East"}}, {19450, {8, "East"}}, {19451, {8, "East"}}, {19452, {8, "East"}}, {19453, {8, "East"}},
    {19454, {9, "East"}}, {19455, {9, "East"}}, {19457, {9, "East"}}, {19458, {9, "East"}}, {19460, {9, "East"}},
    {19463, {10, "East"}}, {19464, {10, "East"}}, {19465, {10, "East"}}, {19466, {10, "East"}}, {19467, {10, "East"}},
};

// === DELAY MAP ===
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

// === 3D POSITIONS ===
struct Point3D { double x, y, z; };
std::map<int, Point3D> chanID_to_Point3D = {
{6 ,{-213.4,-175,27.8742}},{7,{213.4,-175,27.8742}},{8,{-213.4,-95,27.8742}},
{9,{213.4,-95,27.8742}},{10,{-213.4,-40,27.8742}},{11,{213.4,-40,27.8742}},
{12,{-213.4,40,27.8742}},{13,{213.4,40,27.8742}},{14,{-213.4,95,27.8742}},
{18,{213.4,95,27.8742}},{16,{-213.4,175,27.8742}},{17,{213.4,175,27.8742}},
{36,{-213.4,-135,57.8742}},{37,{213.4,-135,57.8742}},{38,{-213.4,0,57.8742}},
{39,{213.4,0,57.8742}},{40,{-213.4,135,57.8742}},{41,{213.4,135,57.8742}},
{60,{-213.4,-175,87.8742}},{61,{213.4,-175,87.8742}},{62,{-213.4,-95,87.8742}},
{63,{213.4,-95,87.8742}},{64,{-213.4,-40,87.8742}},{65,{213.4,-40,87.8742}},
{66,{-213.4,40,87.8742}},{67,{213.4,40,87.8742}},{68,{-213.4,95,87.8742}},
{69,{213.4,95,87.8742}},{70,{-213.4,175,87.8742}},{71,{213.4,175,87.8742}},
{84,{-213.4,-175,161.188}},{85,{213.4,-175,161.188}},{86,{-213.4,-95,161.188}},
{87,{213.4,-95,161.188}},{88,{-213.4,-40,161.188}},{89,{213.4,-40,161.188}},
{90,{-213.4,40,161.188}},{91,{213.4,40,161.188}},{92,{-213.4,95,161.188}},
{93,{213.4,95,161.188}},{94,{-213.4,175,161.188}},{95,{213.4,175,161.188}},
{114,{-213.4,-135,191.188}},{118,{213.4,-135,191.188}},{116,{-213.4,0,191.188}},
{117,{213.4,0,191.188}},{118,{-213.4,135,191.188}},{119,{213.4,135,191.188}},
{138,{-213.4,-175,221.188}},{139,{213.4,-175,221.188}},{140,{-213.4,-95,221.188}},
{141,{213.4,-95,221.188}},{142,{-213.4,-40,221.188}},{143,{213.4,-40,221.188}},
{144,{-213.4,40,221.188}},{145,{213.4,40,221.188}},{146,{-213.4,95,221.188}},
{147,{213.4,95,221.188}},{148,{-213.4,175,221.188}},{149,{213.4,175,221.188}},
{162,{-213.4,-175,288.242}},{163,{213.4,-175,288.242}},{164,{-213.4,-95,288.242}},
{165,{213.4,-95,288.242}},{166,{-213.4,-40,288.242}},{167,{213.4,-40,288.242}},
{168,{-213.4,40,288.242}},{169,{213.4,40,288.242}},{170,{-213.4,95,288.242}},
{171,{213.4,95,288.242}},{172,{-213.4,175,288.242}},{173,{213.4,175,288.242}},
{192,{-213.4,-135,318.242}},{193,{213.4,-135,318.242}},{194,{-213.4,0,318.242}},
{195,{213.4,0,318.242}},{196,{-213.4,135,318.242}},{197,{213.4,135,318.242}},
{216,{-213.4,-175,348.242}},{217,{213.4,-175,348.242}},{218,{-213.4,-95,348.242}},
{219,{213.4,-95,348.242}},{220,{-213.4,-40,348.242}},{221,{213.4,-40,348.242}},
{222,{-213.4,40,348.242}},{223,{213.4,40,348.242}},{224,{-213.4,95,348.242}},
{225,{213.4,95,348.242}},{226,{-213.4,175,348.242}},{227,{213.4,175,348.242}},
{240,{-213.4,-175,421.526}},{241,{213.4,-175,421.526}},{242,{-213.4,-95,421.526}},
{243,{213.4,-95,421.526}},{244,{-213.4,-40,421.526}},{245,{213.4,-40,421.526}},
{246,{-213.4,40,421.526}},{247,{213.4,40,421.526}},{248,{-213.4,95,421.526}},
{249,{213.4,95,421.526}},{250,{-213.4,175,421.526}},{251,{213.4,175,421.526}},
{270,{-213.4,-135,451.526}},{271,{213.4,-135,451.526}},{272,{-213.4,0,451.526}},
{273,{213.4,0,451.526}},{274,{-213.4,135,451.526}},{275,{213.4,135,451.526}},
{294,{-213.4,-175,481.526}},{295,{213.4,-175,481.526}},{296,{-213.4,-95,481.526}},
{297,{213.4,-95,481.526}},{298,{-213.4,-40,481.526}},{299,{213.4,-40,481.526}},
{300,{-213.4,40,481.526}},{301,{213.4,40,481.526}},{302,{-213.4,95,481.526}},
{303,{213.4,95,481.526}},{304,{-213.4,175,481.526}},{305,{213.4,175,481.526}}};

void PMTWaveformTree::Loop()
{
    if (fChain == 0) return;

    const int counter_index = 0, mean_index = 1, rms_index = 2;
    double der_threshold = 0.95;

    // === INITIALIZE HISTOGRAMS ===
    std::map<int, TH1F*> h_time, h_timeSinceAvg, h_timeHalfVsMin;
    std::map<int, TH2F*> h_timeSinceAvg_VS_avgTime;
    std::map<int, std::array<double, 3>> mean_rms, mean_rms_trunc;

    for (const auto& p : chanID_to_delayNS) {
        int ch = p.first;
        h_time[ch] = new TH1F(Form("h_time_Ch%d", ch), Form("h_time_Ch%d", ch), 2000, 0, 300000);
        h_timeSinceAvg[ch] = new TH1F(Form("h_timeSinceAvg_Ch%d", ch), Form("h_timeSinceAvg_Ch%d", ch), 1000, -500000, 500000);
        h_timeSinceAvg_VS_avgTime[ch] = new TH2F(Form("h_timeSinceAvg_VS_avgTime_Ch%d", ch), Form("h_timeSinceAvg_VS_avgTime_Ch%d", ch), 2000, 0, 300000, 1000, -500000, 500000);
        h_timeHalfVsMin[ch] = new TH1F(Form("h_timeHalfVsMin_Ch%d", ch), Form("h_timeHalfVsMin_Ch%d", ch), 25, 0, 10);
        mean_rms[ch] = {0,0,0};
        mean_rms_trunc[ch] = {0,0,0};
    }

    TH1F* h_avgTime = new TH1F("h_avgTime", "h_avgTime", 2000, -300000, 300000);

    int numPMTs_z = 18, numPMTs_y = 9;
    TH2F* h_PMTtimes_posX_RMS = new TH2F("h_PMTtimes_posX_RMS", "Signal Delay RMS - West PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_negX_RMS = new TH2F("h_PMTtimes_negX_RMS", "Signal Delay RMS - East PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_posX_mean = new TH2F("h_PMTtimes_posX_mean", "Mean Signal Delay - West PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_negX_mean = new TH2F("h_PMTtimes_negX_mean", "Mean Signal Delay - East PMTs", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_posX_counts = new TH2F("h_PMTtimes_posX_counts", "h_PMTtimes_posX_counts", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_negX_counts = new TH2F("h_PMTtimes_negX_counts", "h_PMTtimes_negX_counts", numPMTs_z, 20, 490, numPMTs_y, -200, 200);

    // Gaussian fit result histograms
    TH2F* h_PMTtimes_posX_gaus_mean = new TH2F("h_PMTtimes_posX_gaus_mean", "h_PMTtimes_posX_gaus_mean", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_negX_gaus_mean = new TH2F("h_PMTtimes_negX_gaus_mean", "h_PMTtimes_negX_gaus_mean", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_posX_gaus_sigma = new TH2F("h_PMTtimes_posX_gaus_sigma", "h_PMTtimes_posX_gaus_sigma", numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtimes_negX_gaus_sigma = new TH2F("h_PMTtimes_negX_gaus_sigma", "h_PMTtimes_negX_gaus_sigma", numPMTs_z, 20, 490, numPMTs_y, -200, 200);

    // NEW: Amplitude and Total PE histograms
    TH2F* h_PMTamplitude_posX = new TH2F("h_PMTamplitude_posX", "Signal Amplitude - West PMTs", 
                                          numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTamplitude_negX = new TH2F("h_PMTamplitude_negX", "Signal Amplitude - East PMTs", 
                                          numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtotalPE_posX = new TH2F("h_PMTtotalPE_posX", "Total PE - West PMTs", 
                                        numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTtotalPE_negX = new TH2F("h_PMTtotalPE_negX", "Total PE - East PMTs", 
                                        numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    
    // NEW: Signal count histograms (counts signals, not timing measurements)
    TH2F* h_PMTsignals_posX = new TH2F("h_PMTsignals_posX", "Number of Signals - West PMTs", 
                                        numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    TH2F* h_PMTsignals_negX = new TH2F("h_PMTsignals_negX", "Number of Signals - East PMTs", 
                                        numPMTs_z, 20, 490, numPMTs_y, -200, 200);
    
    // === PASS 1: COLLECT DATA ===
    Long64_t nentries = fChain->GetEntriesFast();
    std::vector<std::pair<int,int>> events;
    std::map<std::pair<int,int>, std::vector<std::pair<int,float>>> times, signals;
    int prev_run = -1, prev_evt = -1;

    cout << "=== PASS 1: Collecting timing data from " << nentries << " entries ===" << endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        LoadTree(i); fChain->GetEntry(i);

        if (!Derivative || Derivative->empty()) continue;

        if (prev_evt != Event && prev_evt != -1) {
            events.push_back({prev_run, prev_evt});
        }

        int t_min = std::distance(Derivative->begin(), std::min_element(Derivative->begin(), Derivative->end()));
        int t_half = bin_of_halfmin(*Derivative);
        if (t_half < 0) continue;

        // ===== COMPUTE AMPLITUDE AND TOTAL PE =====
        // ===== COMPUTE AMPLITUDE AND TOTAL PE =====
        if (Waveform && Waveform->size() > 200 && chanID_to_Point3D.count(Channel)) {

            bool candidate = (t_half >= 0);

            if (candidate) {

                // Calculate baseline, amplitude, totalPE
                float baseline = 0;
                for (int j = 0; j < 100; ++j) {
                    baseline += Waveform->at(j);
                }
                baseline /= 100.0f;
                
                int min_tick = std::distance(Waveform->begin(), 
                                            std::min_element(Waveform->begin(), Waveform->end()));
                float min_val = Waveform->at(min_tick);
                
                float amplitude = baseline - min_val;
                
                float totalPE = 0;
                int pulse_start = std::max(0, min_tick - 50);
                int pulse_end = std::min((int)Waveform->size(), min_tick + 200);
                
                for (int j = pulse_start; j < pulse_end; ++j) {
                    totalPE += (baseline - Waveform->at(j));
                }

                // Decide if it's a real signal — tune these thresholds!
                bool has_valid_signal = 
                    (amplitude > 4.0) &&               // ← start low, increase if too many
                    (totalPE > 20.0) &&                // ← adjust based on your data
                    (Derivative->at(t_min) < 0.90);    // ← mild derivative cut

                if (has_valid_signal) {

                    // Debug print (shows typical values of accepted signals)
                    static int count = 0;
                    if (++count % 50 == 0) {  // print every 50 accepted one
                        cout << "Accepted signal: ch=" << Channel 
                            << "  deriv_min=" << Derivative->at(t_min)
                            << "  amp=" << amplitude 
                            << "  pe=" << totalPE << endl;
                    }

                    // Fill only real signals
                    double z_pos = chanID_to_Point3D[Channel].z;
                    double y_pos = chanID_to_Point3D[Channel].y;
                    double x_pos = chanID_to_Point3D[Channel].x;
                    
                    if (x_pos > 0) {
                        h_PMTamplitude_posX->Fill(z_pos, y_pos, amplitude);
                        h_PMTtotalPE_posX->Fill(z_pos, y_pos, totalPE);
                        h_PMTsignals_posX->Fill(z_pos, y_pos);
                    } else {
                        h_PMTamplitude_negX->Fill(z_pos, y_pos, amplitude);
                        h_PMTtotalPE_negX->Fill(z_pos, y_pos, totalPE);
                        h_PMTsignals_negX->Fill(z_pos, y_pos);
                    }
                }
            }
        }
        // ===== END AMPLITUDE/PE BLOCK =====
        
        float val = Derivative->at(t_min);

        if (chanID_to_delayNS.count(Channel) == 0) continue;

        float t_calc = 2.0 * (StartTick_Derivative + t_half) + chanID_to_delayNS[Channel];
        times[{Run, Event}].emplace_back(Channel, t_calc);
        signals[{Run, Event}].emplace_back(Channel, val);
        
        if (val < der_threshold) {
            h_time[Channel]->Fill(t_calc);
        }
        
        h_timeHalfVsMin[Channel]->Fill(t_min - t_half);

        if (i == nentries-1) events.push_back({Run, Event});
        prev_run = Run; prev_evt = Event;
    }

    cout << "Pass 1 complete: " << events.size() << " events" << endl;

    // === COMPUTE AVERAGE TIME PER EVENT ===
    std::map<std::pair<int,int>, double> avgTime;
    for (auto& e : events) {
        double sum = 0; int n = 0;
        for (size_t j = 0; j < times[e].size(); ++j) {
            if (signals[e][j].second < der_threshold) {
                sum += times[e][j].second;
                ++n;
            }
        }
        avgTime[e] = n > 0 ? sum / n : 0;
        if (avgTime[e] > 0) h_avgTime->Fill(avgTime[e]);
    }

    cout << "\n=== PASS 2: Computing relative delays ===" << endl;

    // === PASS 2: COMPUTE DELAYS RELATIVE TO AVERAGE ===
    for (auto& e : events) {
        double avg = avgTime[e];
        if (avg == 0) continue;
        
        for (size_t j = 0; j < times[e].size(); ++j) {
            int ch = times[e][j].first;
            float t = times[e][j].second;
            if (signals[e][j].second >= der_threshold) continue;

            double dt = t - avg;
            h_timeSinceAvg[ch]->Fill(dt);
            h_timeSinceAvg_VS_avgTime[ch]->Fill(avg, dt);

            mean_rms[ch][counter_index]++;
            mean_rms[ch][mean_index] += dt;
            mean_rms[ch][rms_index] += dt*dt;
        }
    }

    // Normalize mean and RMS
    for (auto& p : mean_rms) {
        if (p.second[counter_index] > 0) {
            p.second[mean_index] /= p.second[counter_index];
            p.second[rms_index] = sqrt(p.second[rms_index] / p.second[counter_index]);
        }
    }

    // === PASS 3: TRUNCATED STATISTICS ===
    for (auto& e : events) {
        double avg = avgTime[e];
        if (avg == 0) continue;
        
        for (size_t j = 0; j < times[e].size(); ++j) {
            int ch = times[e][j].first;
            float t = times[e][j].second;
            if (signals[e][j].second >= der_threshold) continue;
            double dt = t - avg;
            double m = mean_rms[ch][mean_index], r = mean_rms[ch][rms_index];
            if (dt < m - r || dt > m + r) continue;
            mean_rms_trunc[ch][counter_index]++;
            mean_rms_trunc[ch][mean_index] += dt;
            mean_rms_trunc[ch][rms_index] += dt*dt;
        }
    }
    
    for (auto& p : mean_rms_trunc) {
        if (p.second[counter_index] > 0) {
            p.second[mean_index] /= p.second[counter_index];
            p.second[rms_index] = sqrt(p.second[rms_index] / p.second[counter_index]);
        }
    }

    cout << "Statistics computed for all channels" << endl;

    // === FIT & BUILD 2D MAPS ===
    const float origin_x = 162, origin_y = -213, origin_z = 22;
    std::vector<double> r_pos, r_err_pos, t_pos, t_err_pos;
    std::vector<double> r_neg, r_err_neg, t_neg, t_err_neg;

    cout << "\n=== Fitting histograms and building output ===" << endl;

    for (auto& p : h_timeSinceAvg) {
        int ch = p.first;
        TH1F* h = p.second;
        
        if (h->GetEntries() < 1) continue;

        TF1* f = new TF1("f_gaus", "gaus", -30, 30);
        h->Fit(f, "RQ");
        double gmean  = f->GetParameter(1);
        double gsigma = f->GetParameter(2);

        int bz = h_PMTtimes_posX_RMS->GetXaxis()->FindBin(chanID_to_Point3D[ch].z);
        int by = h_PMTtimes_posX_RMS->GetYaxis()->FindBin(chanID_to_Point3D[ch].y);
        double radius = sqrt(pow(chanID_to_Point3D[ch].x - origin_x, 2) + 
                           pow(chanID_to_Point3D[ch].y - origin_y, 2) + 
                           pow(chanID_to_Point3D[ch].z - origin_z, 2));

        if (chanID_to_Point3D[ch].x > 0) {
            h_PMTtimes_posX_mean->SetBinContent(bz, by, mean_rms_trunc[ch][mean_index]);
            h_PMTtimes_posX_RMS->SetBinContent(bz, by, mean_rms_trunc[ch][rms_index]);
            h_PMTtimes_posX_counts->SetBinContent(bz, by, mean_rms_trunc[ch][counter_index]);
            h_PMTtimes_posX_gaus_mean->SetBinContent(bz, by, gmean);
            h_PMTtimes_posX_gaus_sigma->SetBinContent(bz, by, gsigma);
            
            if (mean_rms_trunc[ch][counter_index] > 0) {
                r_pos.push_back(radius);
                r_err_pos.push_back(10.16);
                t_pos.push_back(mean_rms_trunc[ch][mean_index]);
                t_err_pos.push_back(gsigma);
            }
        } else {
            h_PMTtimes_negX_mean->SetBinContent(bz, by, mean_rms_trunc[ch][mean_index]);
            h_PMTtimes_negX_RMS->SetBinContent(bz, by, mean_rms_trunc[ch][rms_index]);
            h_PMTtimes_negX_counts->SetBinContent(bz, by, mean_rms_trunc[ch][counter_index]);
            h_PMTtimes_negX_gaus_mean->SetBinContent(bz, by, gmean);
            h_PMTtimes_negX_gaus_sigma->SetBinContent(bz, by, gsigma);
            
            if (mean_rms_trunc[ch][counter_index] > 0) {
                r_neg.push_back(radius);
                r_err_neg.push_back(10.16);
                t_neg.push_back(mean_rms_trunc[ch][mean_index]);
                t_err_neg.push_back(gsigma);
            }
        }
        delete f;
    }

    cout << "West PMTs: " << r_pos.size() << " points for radius plot" << endl;
    cout << "East PMTs: " << r_neg.size() << " points for radius plot" << endl;

    // === SAVE ALL OUTPUTS ===
    cout << "\n=== Saving outputs ===" << endl;
    TFile* outfile = new TFile("output_waveforms_combined.root", "RECREATE");
    
    // Save individual channel histograms
    for (auto& p : h_time) {
        p.second->GetXaxis()->SetTitle("Time [ns]");
        p.second->GetYaxis()->SetTitle("Counts");
        p.second->Write();
    }
    for (auto& p : h_timeSinceAvg) {
        if (p.second->GetEntries() > 0) {
            p.second->GetXaxis()->SetTitle("Time relative to average [ns]");
            p.second->GetYaxis()->SetTitle("Counts");
            p.second->Write();
        }
    }
    for (auto& p : h_timeSinceAvg_VS_avgTime) {
        if (p.second->GetEntries() > 0) {
            p.second->GetXaxis()->SetTitle("Average Time [ns]");
            p.second->GetYaxis()->SetTitle("Time relative to average [ns]");
            p.second->GetZaxis()->SetTitle("Counts");
            p.second->Write();
        }
    }
    
    h_avgTime->GetXaxis()->SetTitle("Average event time [ns]");
    h_avgTime->GetYaxis()->SetTitle("Number of Events");
    h_avgTime->Write();
    
    // Write 2D heatmaps (timing)
    h_PMTtimes_negX_mean->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_negX_mean->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_negX_mean->GetZaxis()->SetTitle("Delay from Average [ns]");
    h_PMTtimes_negX_mean->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_negX_mean->SetStats(0);
    h_PMTtimes_negX_mean->Write();
    
    h_PMTtimes_posX_mean->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_posX_mean->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_posX_mean->GetZaxis()->SetTitle("Delay from Average [ns]");
    h_PMTtimes_posX_mean->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_posX_mean->SetStats(0);
    h_PMTtimes_posX_mean->Write();
    
    h_PMTtimes_negX_RMS->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_negX_RMS->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_negX_RMS->GetZaxis()->SetTitle("RMS [ns]");
    h_PMTtimes_negX_RMS->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_negX_RMS->SetStats(0);
    h_PMTtimes_negX_RMS->Write();
    
    h_PMTtimes_posX_RMS->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_posX_RMS->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_posX_RMS->GetZaxis()->SetTitle("RMS [ns]");
    h_PMTtimes_posX_RMS->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_posX_RMS->SetStats(0);
    h_PMTtimes_posX_RMS->Write();
    
    h_PMTtimes_negX_gaus_mean->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_negX_gaus_mean->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_negX_gaus_mean->GetZaxis()->SetTitle("Gaussian Mean [ns]");
    h_PMTtimes_negX_gaus_mean->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_negX_gaus_mean->SetStats(0);
    h_PMTtimes_negX_gaus_mean->Write();
    
    h_PMTtimes_posX_gaus_mean->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_posX_gaus_mean->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_posX_gaus_mean->GetZaxis()->SetTitle("Gaussian Mean [ns]");
    h_PMTtimes_posX_gaus_mean->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_posX_gaus_mean->SetStats(0);
    h_PMTtimes_posX_gaus_mean->Write();
    
    h_PMTtimes_negX_gaus_sigma->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_negX_gaus_sigma->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_negX_gaus_sigma->GetZaxis()->SetTitle("Gaussian Sigma [ns]");
    h_PMTtimes_negX_gaus_sigma->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_negX_gaus_sigma->SetStats(0);
    h_PMTtimes_negX_gaus_sigma->Write();
    
    h_PMTtimes_posX_gaus_sigma->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_posX_gaus_sigma->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_posX_gaus_sigma->GetZaxis()->SetTitle("Gaussian Sigma [ns]");
    h_PMTtimes_posX_gaus_sigma->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_posX_gaus_sigma->SetStats(0);
    h_PMTtimes_posX_gaus_sigma->Write();
    
    h_PMTtimes_negX_counts->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_negX_counts->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_negX_counts->GetZaxis()->SetTitle("Counts");
    h_PMTtimes_negX_counts->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_negX_counts->SetStats(0);
    h_PMTtimes_negX_counts->Write();
    
    h_PMTtimes_posX_counts->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtimes_posX_counts->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtimes_posX_counts->GetZaxis()->SetTitle("Counts");
    h_PMTtimes_posX_counts->GetZaxis()->SetTitleOffset(1.4);
    h_PMTtimes_posX_counts->SetStats(0);
    h_PMTtimes_posX_counts->Write();

    // Write amplitude and total PE heatmaps
    h_PMTamplitude_posX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTamplitude_posX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTamplitude_posX->GetZaxis()->SetTitle("Average Amplitude [ADC]");
    h_PMTamplitude_posX->GetZaxis()->SetTitleOffset(1.6);
    h_PMTamplitude_posX->SetStats(0);
    h_PMTamplitude_posX->Write();
    
    h_PMTamplitude_negX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTamplitude_negX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTamplitude_negX->GetZaxis()->SetTitle("Average Amplitude [ADC]");
    h_PMTamplitude_negX->GetZaxis()->SetTitleOffset(1.6);
    h_PMTamplitude_negX->SetStats(0);
    h_PMTamplitude_negX->Write();
    
    h_PMTtotalPE_posX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtotalPE_posX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtotalPE_posX->GetZaxis()->SetTitle("Total PE [ADC*ticks]");
    h_PMTtotalPE_posX->GetZaxis()->SetTitleOffset(1.5);
    h_PMTtotalPE_posX->SetStats(0);
    h_PMTtotalPE_posX->Write();
    
    h_PMTtotalPE_negX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTtotalPE_negX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTtotalPE_negX->GetZaxis()->SetTitle("Total PE [ADC*ticks]");
    h_PMTtotalPE_negX->GetZaxis()->SetTitleOffset(1.5);
    h_PMTtotalPE_negX->SetStats(0);
    h_PMTtotalPE_negX->Write();

    h_PMTtotalPE_negX->Write();

    // NEW: signal count heatmaps
    h_PMTsignals_posX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTsignals_posX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTsignals_posX->GetZaxis()->SetTitle("Number of Signals");
    h_PMTsignals_posX->GetZaxis()->SetTitleOffset(1.4);
    h_PMTsignals_posX->SetStats(0);
    h_PMTsignals_posX->Write();

    h_PMTsignals_negX->GetXaxis()->SetTitle("Z [cm]");
    h_PMTsignals_negX->GetYaxis()->SetTitle("Y [cm]");
    h_PMTsignals_negX->GetZaxis()->SetTitle("Number of Signals");
    h_PMTsignals_negX->GetZaxis()->SetTitleOffset(1.4);
    h_PMTsignals_negX->SetStats(0);
    h_PMTsignals_negX->Write();

    // === CREATE AND SAVE PLOTS TO Plots/ DIRECTORY ===
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("Plots", kTRUE);
    gStyle->SetOptStat(0);
    
    // 1. Average time histogram
    TCanvas* c_avgTime = new TCanvas("c_avgTime", "c_avgTime", 800, 600);
    h_avgTime->Draw();
    c_avgTime->SaveAs("Plots/h_avgTime.png");
    c_avgTime->Write();
    delete c_avgTime;
    
    // 2. RMS heatmaps
    TCanvas* canvas_posX_RMS = new TCanvas("canvas_posX_RMS", "canvas_posX_RMS", 1000, 800);
    canvas_posX_RMS->SetRightMargin(0.18);
    canvas_posX_RMS->cd();
    h_PMTtimes_posX_RMS->Draw("COLZ");
    canvas_posX_RMS->SaveAs("Plots/h_PMTtimes_posX_RMS.png");
    canvas_posX_RMS->Write();
    delete canvas_posX_RMS;

    TCanvas* canvas_negX_RMS = new TCanvas("canvas_negX_RMS", "canvas_negX_RMS", 1000, 800);
    canvas_negX_RMS->SetRightMargin(0.18);
    canvas_negX_RMS->cd();
    h_PMTtimes_negX_RMS->Draw("COLZ");
    canvas_negX_RMS->SaveAs("Plots/h_PMTtimes_negX_RMS.png");
    canvas_negX_RMS->Write();
    delete canvas_negX_RMS;

    // 3. Mean delay heatmaps
    TCanvas* c_posX_mean = new TCanvas("c_posX_mean", "c_posX_mean", 1000, 800);
    c_posX_mean->SetRightMargin(0.18);
    h_PMTtimes_posX_mean->Draw("COLZ");
    c_posX_mean->SaveAs("Plots/h_PMTtimes_posX_mean.png");
    c_posX_mean->Write();
    delete c_posX_mean;

    TCanvas* c_negX_mean = new TCanvas("c_negX_mean", "c_negX_mean", 1000, 800);
    c_negX_mean->SetRightMargin(0.18);
    h_PMTtimes_negX_mean->Draw("COLZ");
    c_negX_mean->SaveAs("Plots/h_PMTtimes_negX_mean.png");
    c_negX_mean->Write();
    delete c_negX_mean;

    // 4. Gaussian fit mean heatmaps
    TCanvas* c_posX_gaus_mean = new TCanvas("c_posX_gaus_mean", "c_posX_gaus_mean", 1000, 800);
    c_posX_gaus_mean->SetRightMargin(0.18);
    h_PMTtimes_posX_gaus_mean->Draw("COLZ");
    c_posX_gaus_mean->SaveAs("Plots/h_PMTtimes_posX_gaus_mean.png");
    c_posX_gaus_mean->Write();
    delete c_posX_gaus_mean;

    TCanvas* c_negX_gaus_mean = new TCanvas("c_negX_gaus_mean", "c_negX_gaus_mean", 1000, 800);
    c_negX_gaus_mean->SetRightMargin(0.18);
    h_PMTtimes_negX_gaus_mean->Draw("COLZ");
    c_negX_gaus_mean->SaveAs("Plots/h_PMTtimes_negX_gaus_mean.png");
    c_negX_gaus_mean->Write();
    delete c_negX_gaus_mean;

    // 5. Gaussian fit sigma heatmaps
    TCanvas* c_posX_gaus_sigma = new TCanvas("c_posX_gaus_sigma", "c_posX_gaus_sigma", 1000, 800);
    c_posX_gaus_sigma->SetRightMargin(0.18);
    h_PMTtimes_posX_gaus_sigma->Draw("COLZ");
    c_posX_gaus_sigma->SaveAs("Plots/h_PMTtimes_posX_gaus_sigma.png");
    c_posX_gaus_sigma->Write();
    delete c_posX_gaus_sigma;

    TCanvas* c_negX_gaus_sigma = new TCanvas("c_negX_gaus_sigma", "c_negX_gaus_sigma", 1000, 800);
    c_negX_gaus_sigma->SetRightMargin(0.18);
    h_PMTtimes_negX_gaus_sigma->Draw("COLZ");
    c_negX_gaus_sigma->SaveAs("Plots/h_PMTtimes_negX_gaus_sigma.png");
    c_negX_gaus_sigma->Write();
    delete c_negX_gaus_sigma;

    // 6. NEW: Amplitude heatmaps
    TCanvas* c_posX_amp = new TCanvas("c_posX_amp", "c_posX_amp", 1000, 800);
    c_posX_amp->SetRightMargin(0.18);
    h_PMTamplitude_posX->Draw("COLZ");
    c_posX_amp->SaveAs("Plots/h_PMTamplitude_posX.png");
    c_posX_amp->Write();
    delete c_posX_amp;

    TCanvas* c_negX_amp = new TCanvas("c_negX_amp", "c_negX_amp", 1000, 800);
    c_negX_amp->SetRightMargin(0.18);
    h_PMTamplitude_negX->Draw("COLZ");
    c_negX_amp->SaveAs("Plots/h_PMTamplitude_negX.png");
    c_negX_amp->Write();
    delete c_negX_amp;

    // 7. NEW: Total PE heatmaps
    TCanvas* c_posX_pe = new TCanvas("c_posX_pe", "c_posX_pe", 1000, 800);
    c_posX_pe->SetRightMargin(0.18);
    h_PMTtotalPE_posX->Draw("COLZ");
    c_posX_pe->SaveAs("Plots/h_PMTtotalPE_posX.png");
    c_posX_pe->Write();
    delete c_posX_pe;

    TCanvas* c_negX_pe = new TCanvas("c_negX_pe", "c_negX_pe", 1000, 800);
    c_negX_pe->SetRightMargin(0.18);
    h_PMTtotalPE_negX->Draw("COLZ");
    c_negX_pe->SaveAs("Plots/h_PMTtotalPE_negX.png");
    c_negX_pe->Write();
    delete c_negX_pe;

    // NEW: Signal count heatmaps
    TCanvas* c_posX_signals = new TCanvas("c_posX_signals", "c_posX_signals", 1000, 800);
    c_posX_signals->SetRightMargin(0.18);
    h_PMTsignals_posX->Draw("COLZ");
    c_posX_signals->SaveAs("Plots/h_PMTsignals_posX.png");
    c_posX_signals->Write();
    delete c_posX_signals;

    TCanvas* c_negX_signals = new TCanvas("c_negX_signals", "c_negX_signals", 1000, 800);
    c_negX_signals->SetRightMargin(0.18);
    h_PMTsignals_negX->Draw("COLZ");
    c_negX_signals->SaveAs("Plots/h_PMTsignals_negX.png");
    c_negX_signals->Write();
    delete c_negX_signals;

    // 8. Radius vs Delay plots
    if (r_pos.size() >= 2) {
        TF1* f_line_posx = new TF1("f_line_posx", "[0]+[1]*x", 0, 750);
        TGraphErrors* gr_tbar_vs_radius_posx = new TGraphErrors(r_pos.size(), r_pos.data(), t_pos.data(),
                                                                r_err_pos.data(), t_err_pos.data());
        if (r_pos.size() >= 3) gr_tbar_vs_radius_posx->Fit(f_line_posx, "R");
        gr_tbar_vs_radius_posx->SetTitle("PMT Mean Delay Time - West PMTs");
        gr_tbar_vs_radius_posx->GetXaxis()->SetTitle("Radius w.r.t. Bottom Southwest Corner [cm]");
        gr_tbar_vs_radius_posx->GetYaxis()->SetTitle("PMT Mean Delay Time [ns]");
        gr_tbar_vs_radius_posx->SetMarkerStyle(20);
        gr_tbar_vs_radius_posx->SetMarkerColor(kBlack);
        gr_tbar_vs_radius_posx->Write();

        TCanvas* c2 = new TCanvas("c2", "c2", 900, 700);
        c2->SetLeftMargin(0.18); c2->SetBottomMargin(0.12);
        gr_tbar_vs_radius_posx->Draw("AP");
        if (r_pos.size() >= 3) f_line_posx->Draw("same");
        c2->SaveAs("Plots/gr_tbar_vs_radius_posx.png");
        c2->Write();
        delete c2;
    }

    if (r_neg.size() >= 2) {
        TF1* f_line_negx = new TF1("f_line_negx", "[0]+[1]*x", 0, 750);
        TGraphErrors* gr_tbar_vs_radius_negx = new TGraphErrors(r_neg.size(), r_neg.data(), t_neg.data(),
                                                                r_err_neg.data(), t_err_neg.data());
        if (r_neg.size() >= 3) gr_tbar_vs_radius_negx->Fit(f_line_negx, "R");
        gr_tbar_vs_radius_negx->SetTitle("PMT Mean Delay Time - East PMTs");
        gr_tbar_vs_radius_negx->GetXaxis()->SetTitle("Radius w.r.t. Bottom Southwest Corner [cm]");
        gr_tbar_vs_radius_negx->GetYaxis()->SetTitle("PMT Mean Delay Time [ns]");
        gr_tbar_vs_radius_negx->SetMarkerStyle(20);
        gr_tbar_vs_radius_negx->SetMarkerColor(kBlack);
        gr_tbar_vs_radius_negx->Write();

        TCanvas* c3 = new TCanvas("c3", "c3", 900, 700);
        c3->SetLeftMargin(0.18); c3->SetBottomMargin(0.12);
        gr_tbar_vs_radius_negx->Draw("AP");
        if (r_neg.size() >= 3) f_line_negx->Draw("same");
        c3->SaveAs("Plots/gr_tbar_vs_radius_negx.png");
        c3->Write();
        delete c3;
    }

    // 9. Combined plot
    if (r_pos.size() >= 2 || r_neg.size() >= 2) {
        TH1F* hframe = new TH1F("hframe","PMT Mean Delay vs Distance from Bottom SW Corner",100,0,700);
        hframe->SetMinimum(-2000);
        hframe->SetMaximum(2000);
        hframe->GetXaxis()->SetTitle("Radius [cm]");
        hframe->GetYaxis()->SetTitle("Mean delay relative to event average [ns]");
        hframe->SetStats(0);

        TGraphErrors* gr_pos = new TGraphErrors(r_pos.size(), r_pos.data(), t_pos.data(),
                                                r_err_pos.data(), t_err_pos.data());
        TGraphErrors* gr_neg = new TGraphErrors(r_neg.size(), r_neg.data(), t_neg.data(),
                                                r_err_neg.data(), t_err_neg.data());
        gr_pos->SetMarkerColor(kRed);   gr_pos->SetMarkerStyle(20); gr_pos->SetMarkerSize(1.3);
        gr_neg->SetMarkerColor(kBlue);  gr_neg->SetMarkerStyle(20); gr_neg->SetMarkerSize(1.3);

        TGraphErrors* gfit_pos = nullptr; TGraphErrors* gfit_neg = nullptr;
        TF1* fit_pos = nullptr; TF1* fit_neg = nullptr;

        { // West
            std::vector<double> x,y,ex,ey;
            for (size_t i=0; i<r_pos.size(); ++i)
                if (t_pos[i] > -2000 && t_pos[i] < 2000) {
                    x.push_back(r_pos[i]); y.push_back(t_pos[i]);
                    ex.push_back(r_err_pos[i]); ey.push_back(t_err_pos[i]);
                }
            if (x.size() >= 3) {
                gfit_pos = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
                fit_pos = new TF1("fit_pos", "[0]+[1]*x", 0, 700);
                fit_pos->SetLineColor(kRed);
                gfit_pos->Fit(fit_pos, "RQ");
            }
        }
        { // East
            std::vector<double> x,y,ex,ey;
            for (size_t i=0; i<r_neg.size(); ++i)
                if (t_neg[i] > -2000 && t_neg[i] < 2000) {
                    x.push_back(r_neg[i]); y.push_back(t_neg[i]);
                    ex.push_back(r_err_neg[i]); ey.push_back(t_err_neg[i]);
                }
            if (x.size() >= 3) {
                gfit_neg = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
                fit_neg = new TF1("fit_neg", "[0]+[1]*x", 0, 700);
                fit_neg->SetLineColor(kBlue);
                gfit_neg->Fit(fit_neg, "RQ");
            }
        }

        TCanvas* c_combined = new TCanvas("c_combined", "Combined Radius vs Delay", 950, 700);
        c_combined->SetLeftMargin(0.14); c_combined->SetBottomMargin(0.12);
        hframe->Draw();
        gr_pos->Draw("P same");
        gr_neg->Draw("P same");
        if (fit_pos) fit_pos->Draw("same");
        if (fit_neg) fit_neg->Draw("same");

        TLegend* leg = new TLegend(0.18, 0.80, 0.35, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
        leg->AddEntry(gr_pos, "West PMTs (x > 0)", "p");
        leg->AddEntry(gr_neg, "East PMTs (x < 0)", "p");
        leg->Draw();

        c_combined->SaveAs("Plots/gr_tbar_vs_radius_combined.png");
        c_combined->Write();

        gr_pos->Write("gr_tbar_vs_radius_posx_combined");
        gr_neg->Write("gr_tbar_vs_radius_negx_combined");

        delete hframe; delete c_combined; delete leg;
        delete gr_pos; delete gr_neg;
        delete gfit_pos; delete gfit_neg;
        if (fit_pos) delete fit_pos;
        if (fit_neg) delete fit_neg;
    }

    // 10. Save all individual channel plots
    for (auto& p : h_timeSinceAvg) {
        if (p.second->GetEntries() > 0) {
            TCanvas* c_ch = new TCanvas(Form("c_ch%d", p.first), Form("c_ch%d", p.first), 800, 600);
            p.second->Draw();
            c_ch->SaveAs(Form("Plots/h_timeSinceAvg_Ch%d.png", p.first));
            delete c_ch;
        }
    }

    outfile->Close();
    
    cout << "\n=== ANALYSIS COMPLETE ===" << endl;
    cout << "Output file: output_waveforms_combined.root" << endl;
    cout << "\nPlots saved in Plots/ directory:" << endl;
    cout << "  - h_avgTime.png" << endl;
    cout << "  - h_PMTtimes_posX/negX_RMS.png" << endl;
    cout << "  - h_PMTtimes_posX/negX_mean.png" << endl;
    cout << "  - h_PMTtimes_posX/negX_gaus_mean.png" << endl;
    cout << "  - h_PMTtimes_posX/negX_gaus_sigma.png" << endl;
    cout << "  - h_PMTamplitude_posX/negX.png (NEW)" << endl;
    cout << "  - h_PMTtotalPE_posX/negX.png (NEW)" << endl;
    cout << "  - h_PMTsignals_posX/negX.png (NEW)" << endl;
    if (r_pos.size() >= 2) cout << "  - gr_tbar_vs_radius_posx.png" << endl;
    if (r_neg.size() >= 2) cout << "  - gr_tbar_vs_radius_negx.png" << endl;
    cout << "  - gr_tbar_vs_radius_combined.png" << endl;
    cout << "  - Individual channel plots (h_timeSinceAvg_Ch*.png)" << endl;
    cout << "\nHave a nice day!!" << endl;
}