// PMTTimingAnalyzer.C
// ============================================================
// Timing calibration analysis for SBND PMTs using LED diffuser runs.
//
// REVISION 2026-11-12  (per-pulse timing, replaces the 2026-11-11
// averaged-waveform CF15% method -- that one is still available with
// kTimingMode = 1 as a cross-check)
// ------------------------------------------------------------------
//   * Every flash x channel is timed individually (half-max of the first
//     pulse, baseline mean/RMS from the pre-LED ticks).
//   * Signal size Q (pulse integral) is stored per flash.  Because the LED
//     pulse is ~200 ns wide, the first detected photon comes earlier for
//     bright channels; that intensity bias is removed with a walk function
//     W(lnQ) fitted from the flash-to-flash Q variation WITHIN each channel
//     (geometry is fixed within a channel, so W cannot eat the r/v slope).
//   * A per-flash common offset o_evt (trigger/LED/DAQ jitter) is fitted at
//     the same time:  t(ch,evt) = T_ch + o_evt + W(lnQ) + noise.
//   * New plots: signal size vs distance and vs time delay (see below).
//
// REVISION 2026-11-12 (b)
// ------------------------------------------------------------------
//   * Size correction at channel level now comes from a GLOBAL fit over all
//     diffusers of one side (ALL mode):
//        R(ch,d) - r/v = c_ch + D_d + Wg(lnQ_ch,d)
//     R = channel mean time, c_ch = channel offset (the calibration
//     constant), D_d = diffuser offset.  The within-channel W above is
//     kept only to reduce the flash-to-flash jitter: on the 2026-11-12 run
//     it was too shallow (common LED-brightness jitter is absorbed by
//     o_evt), leaving ~-28 ns per decade of Q in the corrected times.
//     In single-diffuser mode the within-channel W is used (warning printed).
//   * Channel signal size = median lnQ of its flashes.
//   * Radius plots: plain weighted line fit, NO points removed; fixed
//     y-range [kRadiusYmin, kRadiusYmax]; frameless legend.
//   * size_vs_distance_<dk>  : per-channel points (no 2-D binning gaps).
//     size_vs_delay_<dk>     : size on y, delay on x; left = before,
//                              right = after size correction.
//   * Plot-B quantities (outlier / map / summary plots) are now the channel
//     offset from the r/v expectation, not "deviation from fit B".
//   * Channel type comparison plot + calibration table timing_offsets_<tag>.txt
//
// KEY FIX vs previous version
// ----------------------------
// GetEntriesFast() on a TChain returns LLONG_MAX before the chain
// is fully materialised, causing an infinite loop.
// Use GetEntries() instead -- it forces a proper count across all
// files before the event loop begins.
//
// ANALYSIS OVERVIEW
// -----------------
// For each event:
//   1. Theoretical propagation time: t_theo = dist(diffuser,PMT) / c
//   2. Experimental arrival time from waveform derivative half-minimum,
//      corrected by cable delays.
//   3. Both expressed RELATIVE to their per-event average -> t0 cancels.
//
// Two radius plot types per diffuser:
//   (A) radius_vs_dtexp   -- truncated mean of dt_exp_rel vs distance
//   (B) radius_vs_delay   -- mean of (dt_exp_rel - dt_theo_rel) vs distance
//
// MATCHING CONVENTION (updated):
//   West diffusers  →  West (positive-X) APA only
//   East diffusers  →  East (negative-X) APA only
//   Cross-TPC panels are skipped (no LED light → unreliable).
//
// Outlier analysis:
//   For each diffuser + matching APA half:
//     * Fit a line to the data points in plots A and B.
//     * Compute per-channel deviation from that fit.
//     * Flag channels outside ±OUTLIER_NSIGMA·σ, annotate on plot.
//     * Save per-channel deviation scatter (Ch ID vs deviation).
//   In ALL mode: summary scatter (all diffusers overlaid) + 1-D spread histo.
//
// USAGE
// -----
//   root [0] .L PMTTimingAnalyzer.C
//   root [1] PMTWaveformTree t           // auto-loads all EventBuilder*.root
//   root [2] t.Loop()                    // all 10 diffusers
//   root [2] t.Loop("West_3")            // single diffuser
//   root [2] t.Loop("East_1")
//
// Or with an explicit single file:
//   root [0] .L PMTTimingAnalyzer.C
//   root [1] TFile* f = TFile::Open("outfile.root")
//   root [2] TTree* tr; f->GetObject("PMTWaveformTree",tr)
//   root [3] PMTWaveformTree t(tr)
//   root [4] t.Loop()
// ============================================================

#define PMTWaveformTree_cxx
#include "PMTWaveformTree.h"

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TProfile.h>
#include <TPad.h>
#include <functional>
#include <fstream>
#include <TColor.h>

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>

using namespace std;

// ============================================================
// TUNABLE PARAMETERS (edit here)
// ============================================================
static const double kLightSpeed_cm_ns = 22.0;   // cm/ns
static const int    kMinGoodCh        = 5;        // min channels per event
static const double kOutlierNsigma    = 2.0;      // flag beyond ±N·σ

// ------------------------------------------------------------
// TIMING METHOD
// ------------------------------------------------------------
//  kTimingMode = 0 : PER-PULSE (default) -- see REVISION 2026-11-12 above
//  kTimingMode = 1 : AVERAGED WAVEFORM CF (2026-11-11), cross-check only
static const int    kTimingMode       = 0;
static const double kNsPerTick        = 2.0;     // PMT digitizer (V1730, 500 MHz)
static const double kPulseSign        = -1.0;    // PMT pulses go negative
static const int    kNBaselineTicks   = 800;     // baseline mean/RMS from ticks [0,N)

// ---- per-pulse hit finding (kTimingMode = 0) ----
static const double kMaxBaselineRMS   = 8.0;     // ADC (smoothed); drop flashes with noisy/pulsed baseline
static const int    kHitLo            = 1100;     // search window [ticks]
static const int    kHitHi            = 1500;     //   (LED light onset ~1190, prompt light ends ~1450)
static const int    kHitSmooth        = 1;       // +/- ticks moving average
static const double kHitNsigma        = 5.0;     // threshold = max(N * baseline RMS, kHitMinADC)
static const double kHitMinADC        = 5.0;
static const int    kHitMinAbove      = 2;       // consecutive ticks above threshold
static const int    kHitPeakTicks     = 8;       // first-pulse peak searched within N ticks of the crossing
static const double kHitCF            = 0.5;     // half-max of the first pulse (leading edge)

// ---- signal-size (walk) correction + per-flash offsets ----
static const bool   kUseWalkCorr      = true;
static const int    kWalkOrder        = 3;       // polynomial order in ln(Q)
static const bool   kUseEventOffset   = true;    // fit a common per-flash time offset
static const int    kNIter            = 5;       // T_ch / o_evt / W iterations
static const double kTruncNsigma      = 3.0;     // robust (median/MAD-clipped) per-channel mean
static const int    kMinHitsPerCh     = 20;
static const double kMinHitFrac       = 0.2;     // drop channels that see light in <20% of flashes

// ---- global (all-diffuser) size correction, ALL mode only ----
//   R(ch,d) - r/v = c_ch + D_d + Wg(lnQ)   fitted on channel means of one side.
//   A channel is seen at very different brightness by the 5 diffusers of its
//   side, which separates the size effect Wg from the channel offset c_ch.
static const bool   kUseGlobalSizeCorr = true;
static const int    kGlobalWalkOrder   = 3;      // polynomial order in ln(Q)

// ---- radius plots ----
// Plain weighted straight-line fit, every channel used (no clipping).
static const double kRadiusYmin       = -100.0;  // fixed y-range [ns] for radius_vs_dtexp and radius_vs_delay
static const double kRadiusYmax       =  150.0;

// ---- plot selection ----
// false: only the plots used in the technote are written as PNG.  The less
// informative ones (res_mean per diffuser, res_sigma_All, devhist_*,
// devA_diffuser_vs_channel, summary_scatter_devFromFitA_vs_ch,
// summary_spread_histos, residual_mean_vs_channel) are still written to the
// ROOT file but not saved as PNG.  Set true to get them back.
static const bool   kSaveExtraPlots   = false;
// res_sigma maps: colour scale capped at kSigmaCapFactor x median sigma of the
// map, so one broken channel (e.g. PMT 118, sigma ~300 ns) cannot wash out the
// map.  Channels above the cap are drawn at the cap and listed on the plot.
static const double kSigmaCapFactor   = 3.0;

// ---- averaged waveform (kTimingMode = 1 only) ----
static const int    kWfLo             = 1000;
static const int    kWfHi             = 1600;
static const int    kWfSmooth         = 2;
static const double kCFfrac           = 0.5;     // was 0.15 on 2026-11-11
static const double kMinAvgPeakADC    = 3.0;

// ============================================================
// RUN -> DIFFUSER MAPPING
// ============================================================
struct DiffuserInfo { int index; string side; };
static map<int, DiffuserInfo> run_to_diffuser = {
    {19403,{3,"West"}},{19406,{3,"West"}},{19407,{3,"West"}},{19408,{3,"West"}},{19409,{3,"West"}},
    {19410,{1,"West"}},{19411,{1,"West"}},{19413,{1,"West"}},{19414,{1,"West"}},{19418,{1,"West"}},
    {19417,{5,"West"}},{19419,{5,"West"}},{19420,{5,"West"}},{19421,{5,"West"}},{19423,{5,"West"}},
    {19424,{2,"West"}},{19425,{2,"West"}},{19426,{2,"West"}},{19427,{2,"West"}},{19428,{2,"West"}},
    {19429,{4,"West"}},{19430,{4,"West"}},{19431,{4,"West"}},{19432,{4,"West"}},{19433,{4,"West"}},
    {19434,{1,"East"}},{19435,{1,"East"}},{19436,{1,"East"}},{19437,{1,"East"}},{19438,{1,"East"}},
    {19439,{3,"East"}},{19440,{3,"East"}},{19441,{3,"East"}},{19442,{3,"East"}},{19443,{3,"East"}},
    {19449,{5,"East"}},{19450,{5,"East"}},{19451,{5,"East"}},{19452,{5,"East"}},{19453,{5,"East"}},
    {19454,{4,"East"}},{19455,{4,"East"}},{19457,{4,"East"}},{19458,{4,"East"}},{19460,{4,"East"}},
    {19463,{2,"East"}},{19464,{2,"East"}},{19465,{2,"East"}},{19466,{2,"East"}},{19467,{2,"East"}},
};

// ============================================================
// DIFFUSER POSITIONS  (cm)
// ============================================================
struct Point3D { double x, y, z; };
static map<int, Point3D> west_diffuser_pos = {
    {1,{0.0, 86.36,392.48}},{2,{0.0,-118.97,392.48}},{3,{0.0,-12.95,246.23}},
    {4,{0.0, 86.60,134.95}},{5,{0.0,-117.46,135.10}},
};
static map<int, Point3D> east_diffuser_pos = {
    {1,{0.0, 86.82,131.62}},{2,{0.0,-116.81,131.62}},{3,{0.0,-13.33,261.65}},
    {4,{0.0, 85.40,388.95}},{5,{0.0,-119.98,389.05}},
};

// ============================================================
// CABLE DELAYS  (ns)
// ============================================================
static map<int, double> chanID_to_delayNS = {
    {250,5.44},{304,1.22},{274,4.10},{248,0.25},{302,5.54},{246,3.18},
    {300,3.63},{272,2.75},{244,3.22},{298,3.90},{242,3.80},{296,0.00},
    {270,0.76},{240,1.80},{294,0.10},{172,1.95},{226,1.58},{196,0.20},
    {170,1.90},{224,3.77},{168,-0.64},{222,-0.90},{194,2.09},{166,2.30},
    {220,3.85},{164,1.42},{218,-0.35},{192,0.16},{162,-0.60},{216,5.61},
    {94,1.75},{148,-0.93},{118,-1.90},{92,1.20},{146,0.63},{90,5.15},
    {144,3.25},{116,-2.67},{88,2.13},{142,-2.00},{86,0.00},{140,-0.05},
    {114,0.52},{84,-1.45},{138,1.63},{16,1.97},{70,5.67},{40,5.45},
    {14,3.70},{68,5.62},{12,1.00},{66,0.65},{38,1.93},{10,-2.25},
    {64,5.63},{8,-10.22},{62,-4.05},{36,0.45},{6,-2.99},{60,-3.85},
    {71,0.50},{17,-3.32},{41,-1.74},{69,2.38},{15,0.30},{67,2.00},
    {13,-2.74},{39,2.20},{65,-1.15},{11,5.62},{63,0.43},{9,-4.26},
    {37,-3.83},{61,-1.53},{7,-3.07},{149,2.40},{95,-1.00},{119,4.20},
    {147,1.40},{93,-0.70},{145,4.05},{91,4.07},{117,2.10},{143,1.60},
    {89,3.81},{141,-3.05},{87,0.08},{115,-0.68},{139,-1.76},{85,-11.41},
    {227,4.16},{173,1.71},{197,3.28},{225,-2.26},{171,-2.45},{223,1.91},
    {169,0.42},{195,4.04},{221,1.25},{167,1.07},{219,-1.25},{165,2.00},
    {193,-1.70},{217,1.27},{163,1.40},{305,1.27},{251,4.27},{275,1.84},
    {303,5.18},{249,-1.77},{301,1.19},{247,3.09},{273,-0.21},{299,3.63},
    {245,1.75},{297,5.25},{243,5.28},{271,0.75},{295,-2.98},{241,1.24},
};

// ============================================================
// PMT POSITIONS  (cm)
// ============================================================
static map<int, Point3D> chanID_to_pos = {
    {6,{-213.4,-175,27.8742}},{7,{213.4,-175,27.8742}},
    {8,{-213.4,-95,27.8742}},{9,{213.4,-95,27.8742}},
    {10,{-213.4,-40,27.8742}},{11,{213.4,-40,27.8742}},
    {12,{-213.4,40,27.8742}},{13,{213.4,40,27.8742}},
    {14,{-213.4,95,27.8742}},{18,{213.4,95,27.8742}},
    {16,{-213.4,175,27.8742}},{17,{213.4,175,27.8742}},
    {36,{-213.4,-135,57.8742}},{37,{213.4,-135,57.8742}},
    {38,{-213.4,0,57.8742}},{39,{213.4,0,57.8742}},
    {40,{-213.4,135,57.8742}},{41,{213.4,135,57.8742}},
    {60,{-213.4,-175,87.8742}},{61,{213.4,-175,87.8742}},
    {62,{-213.4,-95,87.8742}},{63,{213.4,-95,87.8742}},
    {64,{-213.4,-40,87.8742}},{65,{213.4,-40,87.8742}},
    {66,{-213.4,40,87.8742}},{67,{213.4,40,87.8742}},
    {68,{-213.4,95,87.8742}},{69,{213.4,95,87.8742}},
    {70,{-213.4,175,87.8742}},{71,{213.4,175,87.8742}},
    {84,{-213.4,-175,161.188}},{85,{213.4,-175,161.188}},
    {86,{-213.4,-95,161.188}},{87,{213.4,-95,161.188}},
    {88,{-213.4,-40,161.188}},{89,{213.4,-40,161.188}},
    {90,{-213.4,40,161.188}},{91,{213.4,40,161.188}},
    {92,{-213.4,95,161.188}},{93,{213.4,95,161.188}},
    {94,{-213.4,175,161.188}},{95,{213.4,175,161.188}},
    {114,{-213.4,-135,191.188}},{118,{213.4,-135,191.188}},
    {116,{-213.4,0,191.188}},{117,{213.4,0,191.188}},
    {119,{213.4,135,191.188}},
    {138,{-213.4,-175,221.188}},{139,{213.4,-175,221.188}},
    {140,{-213.4,-95,221.188}},{141,{213.4,-95,221.188}},
    {142,{-213.4,-40,221.188}},{143,{213.4,-40,221.188}},
    {144,{-213.4,40,221.188}},{145,{213.4,40,221.188}},
    {146,{-213.4,95,221.188}},{147,{213.4,95,221.188}},
    {148,{-213.4,175,221.188}},{149,{213.4,175,221.188}},
    {162,{-213.4,-175,288.242}},{163,{213.4,-175,288.242}},
    {164,{-213.4,-95,288.242}},{165,{213.4,-95,288.242}},
    {166,{-213.4,-40,288.242}},{167,{213.4,-40,288.242}},
    {168,{-213.4,40,288.242}},{169,{213.4,40,288.242}},
    {170,{-213.4,95,288.242}},{171,{213.4,95,288.242}},
    {172,{-213.4,175,288.242}},{173,{213.4,175,288.242}},
    {192,{-213.4,-135,318.242}},{193,{213.4,-135,318.242}},
    {194,{-213.4,0,318.242}},{195,{213.4,0,318.242}},
    {196,{-213.4,135,318.242}},{197,{213.4,135,318.242}},
    {216,{-213.4,-175,348.242}},{217,{213.4,-175,348.242}},
    {218,{-213.4,-95,348.242}},{219,{213.4,-95,348.242}},
    {220,{-213.4,-40,348.242}},{221,{213.4,-40,348.242}},
    {222,{-213.4,40,348.242}},{223,{213.4,40,348.242}},
    {224,{-213.4,95,348.242}},{225,{213.4,95,348.242}},
    {226,{-213.4,175,348.242}},{227,{213.4,175,348.242}},
    {240,{-213.4,-175,421.526}},{241,{213.4,-175,421.526}},
    {242,{-213.4,-95,421.526}},{243,{213.4,-95,421.526}},
    {244,{-213.4,-40,421.526}},{245,{213.4,-40,421.526}},
    {246,{-213.4,40,421.526}},{247,{213.4,40,421.526}},
    {248,{-213.4,95,421.526}},{249,{213.4,95,421.526}},
    {250,{-213.4,175,421.526}},{251,{213.4,175,421.526}},
    {270,{-213.4,-135,451.526}},{271,{213.4,-135,451.526}},
    {272,{-213.4,0,451.526}},{273,{213.4,0,451.526}},
    {274,{-213.4,135,451.526}},{275,{213.4,135,451.526}},
    {294,{-213.4,-175,481.526}},{295,{213.4,-175,481.526}},
    {296,{-213.4,-95,481.526}},{297,{213.4,-95,481.526}},
    {298,{-213.4,-40,481.526}},{299,{213.4,-40,481.526}},
    {300,{-213.4,40,481.526}},{301,{213.4,40,481.526}},
    {302,{-213.4,95,481.526}},{303,{213.4,95,481.526}},
    {304,{-213.4,175,481.526}},{305,{213.4,175,481.526}},
};

static double Dist3(const Point3D& a, const Point3D& b) {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

// Simple stats on a vector: returns {mean, sigma}
// Constant-fraction time (in ticks, relative to the start of `s`) of a
// positive-going, baseline-subtracted average waveform.  Returns -1 if
// the (smoothed) peak is below kMinAvgPeakADC.
static double CFTimeTicks(const vector<double>& s_raw, double frac, double* peakOut=nullptr) {
    int n = (int)s_raw.size();
    if (n < 3) return -1;
    vector<double> s(n, 0.0);
    for (int i = 0; i < n; ++i) {
        int lo = max(0, i-kWfSmooth), hi = min(n-1, i+kWfSmooth);
        double sum = 0; for (int k = lo; k <= hi; ++k) sum += s_raw[k];
        s[i] = sum / (hi-lo+1);
    }
    int ipk = (int)(max_element(s.begin(), s.end()) - s.begin());
    double pk = s[ipk];
    if (peakOut) *peakOut = pk;
    if (pk < kMinAvgPeakADC) return -1;
    double h = frac * pk;
    int j = ipk;
    while (j > 0 && s[j] >= h) --j;                 // walk back down the rising edge
    if (j >= n-1 || s[j+1] == s[j]) return j;
    return j + (h - s[j]) / (s[j+1] - s[j]);        // linear interpolation
}

static pair<double,double> VecStats(const vector<double>& v);
// Write the channel ID (red) next to every point beyond +-kOutlierNsigma*sigma
// of the mean of y (same rule as the outlier_* plots). ylo/yhi = frame range.
static void LabelOutlierChannels(const vector<double>& ch, const vector<double>& y,
                                 double ylo, double yhi, double tsize = 0.030) {
    if (ch.size() < 3 || ch.size() != y.size()) return;
    pair<double,double> ms = VecStats(y);
    double mu = ms.first, sig = ms.second;
    if (sig <= 0) return;
    double off = 0.04*(yhi-ylo);
    TLatex lat; lat.SetTextSize(tsize); lat.SetTextAlign(22); lat.SetTextColor(kRed+1);
    for (size_t k = 0; k < y.size(); ++k) {
        if (fabs(y[k]-mu) <= kOutlierNsigma*sig) continue;
        double ly = (y[k] > mu) ? y[k]+off : y[k]-off;
        ly = min(max(ly, ylo+off), yhi-off);          // keep label inside the frame
        lat.DrawLatex(ch[k], ly, Form("%d", (int)ch[k]));
    }
}

static pair<double,double> VecStats(const vector<double>& v) {
    if (v.empty()) return {0.0, 0.0};
    int n = (int)v.size();
    double sum=0, sum2=0;
    for (double x : v){ sum+=x; sum2+=x*x; }
    double mean  = sum/n;
    double sigma = (n>1) ? sqrt((sum2 - n*mean*mean)/(n-1)) : 0.0;
    return {mean, sigma};
}


// ============================================================
// PER-PULSE HELPERS
// ============================================================
struct PulseHit { bool ok; double tTick; double amp; double charge; double rms; };

// Time the FIRST pulse of one flash on the raw waveform `w`.
//   baseline = mean of [0,kNBaselineTicks); RMS of the smoothed baseline
//   hit      = first kHitMinAbove consecutive ticks above
//              max(kHitNsigma*RMS, kHitMinADC) inside [kHitLo,kHitHi)
//   time     = kHitCF crossing on the leading edge of that pulse (peak
//              searched within kHitPeakTicks), linearly interpolated
//   charge   = integral from the crossing to kHitHi  [ADC*ns]
static PulseHit FindFirstPulse(const vector<int>& w)
{
    PulseHit h{false, -1.0, 0.0, 0.0, 0.0};
    const int n = (int)w.size();
    if (n < kHitHi || kNBaselineTicks > kHitLo) return h;

    double base = 0;
    for (int k = 0; k < kNBaselineTicks; ++k) base += w[k];
    base /= kNBaselineTicks;

    auto sm = [&](int i) {
        int lo = max(0, i-kHitSmooth), hi = min(n-1, i+kHitSmooth);
        double s = 0; for (int k = lo; k <= hi; ++k) s += w[k];
        return kPulseSign * (s/(hi-lo+1) - base);
    };
    double s2 = 0; int nb = 0;
    for (int k = kHitSmooth; k < kNBaselineTicks-kHitSmooth; ++k) { double v = sm(k); s2 += v*v; ++nb; }
    h.rms = nb ? sqrt(s2/nb) : 0.0;
    if (h.rms <= 0 || h.rms > kMaxBaselineRMS) return h;

    const int m = kHitHi - kHitLo;
    vector<double> s(m);
    for (int k = 0; k < m; ++k) s[k] = sm(kHitLo + k);
    const double thr = max(kHitNsigma*h.rms, kHitMinADC);

    int i0 = -1;
    for (int k = 0; k + kHitMinAbove <= m; ++k) {
        bool above = true;
        for (int j = 0; j < kHitMinAbove; ++j) if (s[k+j] <= thr) { above = false; break; }
        if (above) { i0 = k; break; }
    }
    if (i0 < 1) return h;   // no pulse, or already above threshold at the window start

    int iEnd = min(m-1, i0 + kHitPeakTicks), ipk = i0;
    for (int k = i0; k <= iEnd; ++k) if (s[k] > s[ipk]) ipk = k;
    const double lvl = kHitCF * s[ipk];
    int j = ipk;
    while (j > 0 && s[j] >= lvl) --j;
    if (j <= 0) return h;
    double tt = j + ((s[j+1] != s[j]) ? (lvl - s[j]) / (s[j+1] - s[j]) : 0.0);

    double q = 0;
    for (int k = max(0, i0-2); k < m; ++k) q += s[k];

    h.ok     = true;
    h.tTick  = kHitLo + tt;
    h.amp    = *max_element(s.begin(), s.end());
    h.charge = q * kNsPerTick;
    return h;
}

// ============================================================
// PMT TYPE  (sbnd_pds_mapping.json: pd_type = pmt_uncoated)
// ============================================================
static const set<int> kUncoatedPMTs = {
    36,37,38,39,40,41, 114,115,116,117,118,119,
    192,193,194,195,196,197, 270,271,272,273,274,275,
};
static const char* kTypeNameA   = "TPB-coated PMT";
static const char* kTypeNameB   = "uncoated PMT";
static const char* kTypeFileTag = "coated_vs_uncoated";
static bool IsTypeB(int ch) { return kUncoatedPMTs.count(ch) > 0; }

struct LEDHit { int ch; int evt; float t; float u; };   // t [ns] incl. cable delay, u = ln(Q)

struct ChanTime { double mean, err, sigma, meanU, rawMean, rawErr, sdU; int n; double hitFrac; };
//   mean/err : robust mean of t - o_evt - W(lnQ)  (within-channel W)
//   rawMean  : robust mean of t - o_evt           (no size correction)
//   meanU    : MEDIAN lnQ of the channel, sdU = half of its 16-84% spread

struct DiffSolution {
    map<int,ChanTime> ch;
    vector<double>    walk;            // a_1..a_P :  W(u) = sum a_k (u-u0)^k
    double            u0 = 0, uLo = 0, uHi = 0;
    vector<double>    evtOff;
    vector<char>      evtOK;
};

static double WalkEval(const DiffSolution& S, double u)
{
    if (S.walk.empty()) return 0.0;
    double x = min(max(u, S.uLo), S.uHi) - S.u0, p = 1.0, w = 0.0;
    for (double a : S.walk) { p *= x; w += a*p; }
    return w;
}

static double Median(vector<double> v)
{
    if (v.empty()) return 0.0;
    size_t k = v.size()/2;
    nth_element(v.begin(), v.begin()+k, v.end());
    return v[k];
}

// Median/MAD-clipped mean.  Returns false if nothing survives.
static bool RobustMean(const vector<double>& v, double nsig,
                       double& mean, double& err, double& sigma, int& nUsed)
{
    mean = err = sigma = 0; nUsed = 0;
    if (v.empty()) return false;
    double med = Median(v);
    vector<double> d(v.size());
    for (size_t i = 0; i < v.size(); ++i) d[i] = fabs(v[i]-med);
    double cut = nsig * max(1.4826*Median(d), 1e-6);
    double s = 0, s2 = 0; int n = 0;
    for (double x : v) if (fabs(x-med) <= cut) { s += x; s2 += x*x; ++n; }
    if (n == 0) { mean = med; return false; }
    mean  = s/n;
    sigma = (n > 1) ? sqrt(max((s2 - n*mean*mean)/(n-1), 0.0)) : 0.0;
    err   = (n > 1) ? sigma/sqrt((double)n) : sigma;
    nUsed = n;
    return true;
}

// Small dense solver (Gaussian elimination, partial pivoting).
static bool SolveLinear(vector<vector<double>> A, vector<double> b, vector<double>& x)
{
    const int n = (int)b.size();
    for (int c = 0; c < n; ++c) {
        int p = c;
        for (int r = c+1; r < n; ++r) if (fabs(A[r][c]) > fabs(A[p][c])) p = r;
        if (fabs(A[p][c]) < 1e-12) return false;
        swap(A[p], A[c]); swap(b[p], b[c]);
        for (int r = c+1; r < n; ++r) {
            double f = A[r][c]/A[c][c];
            for (int k = c; k < n; ++k) A[r][k] -= f*A[c][k];
            b[r] -= f*b[c];
        }
    }
    x.assign(n, 0.0);
    for (int r = n-1; r >= 0; --r) {
        double s = b[r];
        for (int k = r+1; k < n; ++k) s -= A[r][k]*x[k];
        x[r] = s/A[r][r];
    }
    return true;
}

// Solve  t(ch,evt) = T_ch + o_evt + W(lnQ)  for one diffuser.
//   T_ch  : robust mean of the corrected times
//   o_evt : median over the channels of that flash
//   W     : polynomial in (lnQ-u0), fitted on per-channel-demeaned data
//           (within-channel regression -> geometry cannot leak into W)
static DiffSolution SolveDiffuser(const vector<LEDHit>& hits, int nEvt)
{
    DiffSolution S;
    S.evtOff.assign(nEvt, 0.0);
    S.evtOK.assign(nEvt, 0);
    if (hits.empty()) return S;

    vector<double> us; us.reserve(hits.size());
    for (const auto& h : hits) us.push_back(h.u);
    sort(us.begin(), us.end());
    S.uLo = us[(size_t)(0.02*(us.size()-1))];
    S.uHi = us[(size_t)(0.98*(us.size()-1))];
    S.u0  = us[us.size()/2];

    vector<int> nPerEvt(nEvt, 0);
    for (const auto& h : hits) nPerEvt[h.evt]++;
    for (int e = 0; e < nEvt; ++e) S.evtOK[e] = (nPerEvt[e] >= kMinGoodCh);
    int nEvtOK = 0; for (int e = 0; e < nEvt; ++e) nEvtOK += S.evtOK[e];

    map<int, vector<int>> chHits;
    for (int i = 0; i < (int)hits.size(); ++i)
        if (S.evtOK[hits[i].evt]) chHits[hits[i].ch].push_back(i);
    for (auto it = chHits.begin(); it != chHits.end(); ) {
        int n = (int)it->second.size();
        if (n < kMinHitsPerCh || n < kMinHitFrac*nEvtOK) it = chHits.erase(it);
        else ++it;
    }
    if (chHits.empty()) return S;

    const int P = kUseWalkCorr ? kWalkOrder : 0;
    map<int,double> T, Sig;
    auto corr = [&](int i) {
        return hits[i].t - S.evtOff[hits[i].evt] - WalkEval(S, hits[i].u);
    };

    for (int it = 0; it < kNIter; ++it) {
        // (a) channel times
        for (const auto& kv : chHits) {
            vector<double> v; v.reserve(kv.second.size());
            for (int i : kv.second) v.push_back(corr(i));
            double m, e, sg; int nu;
            RobustMean(v, kTruncNsigma, m, e, sg, nu);
            T[kv.first] = m; Sig[kv.first] = sg;
        }
        // (b) per-flash offsets
        if (kUseEventOffset) {
            vector<vector<double>> per(nEvt);
            for (const auto& kv : chHits)
                for (int i : kv.second)
                    per[hits[i].evt].push_back(hits[i].t - WalkEval(S, hits[i].u) - T.at(kv.first));
            double so = 0; int no = 0;
            for (int e = 0; e < nEvt; ++e) {
                if ((int)per[e].size() >= kMinGoodCh) { S.evtOff[e] = Median(per[e]); so += S.evtOff[e]; ++no; }
                else S.evtOff[e] = 0.0;
            }
            if (no) for (int e = 0; e < nEvt; ++e) if ((int)per[e].size() >= kMinGoodCh) S.evtOff[e] -= so/no;
        }
        // (c) walk: within-channel least squares on clipped hits
        if (P > 0) {
            vector<vector<double>> A(P, vector<double>(P, 0.0));
            vector<double> B(P, 0.0);
            for (const auto& kv : chHits) {
                double Tc = T.at(kv.first), cut = kTruncNsigma*max(Sig.at(kv.first), 1e-3);
                vector<int> in;
                for (int i : kv.second) if (fabs(corr(i) - Tc) <= cut) in.push_back(i);
                if ((int)in.size() < 5) continue;
                vector<vector<double>> bs(in.size(), vector<double>(P));
                vector<double> y(in.size()), bm(P, 0.0);
                double ym = 0;
                for (size_t k = 0; k < in.size(); ++k) {
                    const LEDHit& h = hits[in[k]];
                    double x = min(max((double)h.u, S.uLo), S.uHi) - S.u0, p = 1.0;
                    for (int q = 0; q < P; ++q) { p *= x; bs[k][q] = p; bm[q] += p; }
                    y[k] = h.t - S.evtOff[h.evt]; ym += y[k];
                }
                for (int q = 0; q < P; ++q) bm[q] /= in.size();
                ym /= in.size();
                for (size_t k = 0; k < in.size(); ++k)
                    for (int q = 0; q < P; ++q) {
                        double bq = bs[k][q]-bm[q];
                        B[q] += bq*(y[k]-ym);
                        for (int r = 0; r < P; ++r) A[q][r] += bq*(bs[k][r]-bm[r]);
                    }
            }
            for (int q = 0; q < P; ++q) A[q][q] *= (1.0 + 1e-9);
            vector<double> a;
            if (SolveLinear(A, B, a)) S.walk = a;
        }
    }

    // final per-channel numbers
    for (const auto& kv : chHits) {
        vector<double> v, vr, vu;
        for (int i : kv.second) {
            v.push_back(corr(i));
            vr.push_back(hits[i].t - S.evtOff[hits[i].evt]);
            vu.push_back(hits[i].u);
        }
        ChanTime c{};
        int nu;
        RobustMean(v, kTruncNsigma, c.mean, c.err, c.sigma, nu);
        double s2; int n2;
        RobustMean(vr, kTruncNsigma, c.rawMean, c.rawErr, s2, n2);
        // channel signal size: median lnQ (robust against the low-Q tail of
        // flashes where only a small first pulse was found)
        sort(vu.begin(), vu.end());
        auto qtl = [&](double f){ return vu[(size_t)(f*(vu.size()-1))]; };
        c.meanU = qtl(0.5);
        c.sdU   = 0.5*(qtl(0.84) - qtl(0.16));
        c.n = nu;
        c.hitFrac = nEvtOK ? (double)kv.second.size()/nEvtOK : 0.0;
        S.ch[kv.first] = c;
    }
    return S;
}

// Weighted straight-line fit y = p0 + p1*x (weights 1/ey^2), ALL points used.
struct LineFit { bool ok=false; double p0=0, p1=0, e0=0, e1=0, chi2=0; int ndf=0; };

static LineFit WeightedLineFit(const vector<double>& x, const vector<double>& y,
                               const vector<double>& ey)
{
    LineFit f;
    const int n = (int)x.size();
    double S=0,Sx=0,Sy=0,Sxx=0,Sxy=0;
    for (int i = 0; i < n; ++i) {
        double w = 1.0/max(ey[i]*ey[i], 1e-4);
        S+=w; Sx+=w*x[i]; Sy+=w*y[i]; Sxx+=w*x[i]*x[i]; Sxy+=w*x[i]*y[i];
    }
    double D = S*Sxx - Sx*Sx;
    if (n < 3 || D <= 0) return f;
    f.p1 = (S*Sxy - Sx*Sy)/D;  f.p0 = (Sxx*Sy - Sx*Sxy)/D;
    f.e1 = sqrt(S/D);          f.e0 = sqrt(Sxx/D);
    for (int i = 0; i < n; ++i) {
        double r = y[i] - f.p0 - f.p1*x[i];
        f.chi2 += r*r/max(ey[i]*ey[i], 1e-4);
    }
    f.ndf = n - 2;
    f.ok  = true;
    return f;
}

// ------------------------------------------------------------
// GLOBAL SIZE CORRECTION (one side, all its diffusers)
//   y(ch,d) = R(ch,d) - r/v = c_ch + D_d + sum_p g_p (u - u0)^p ,  u = median lnQ
// Unweighted least squares on the channel means (their statistical errors,
// ~0.5 ns, are far below the real channel-to-channel spread), gauge
// sum_ch c_ch = 0.
// ------------------------------------------------------------
struct GlobalRow  { string dk; int ch; double u, y; };
struct GlobalSize {
    bool ok = false;
    vector<double> g;  double u0 = 0, uLo = 0, uHi = 0;
    map<string,double> D;  map<int,double> c;
    double rms = 0;  int nRows = 0;
};

static double GlobalW(const GlobalSize& G, double u)
{
    if (!G.ok) return 0.0;
    double x = min(max(u, G.uLo), G.uHi) - G.u0, p = 1.0, w = 0.0;
    for (double a : G.g) { p *= x; w += a*p; }
    return w;
}

static GlobalSize FitGlobalSize(const vector<GlobalRow>& rows, int P)
{
    GlobalSize G;
    map<int,int> ci;  map<string,int> di;
    for (const auto& r : rows) { ci.emplace(r.ch, 0); di.emplace(r.dk, 0); }
    if (di.size() < 2 || P < 1) return G;
    int k = 0;
    for (auto& kv : ci) kv.second = k++;
    for (auto& kv : di) kv.second = k++;
    const int nc = (int)ci.size(), nd = (int)di.size(), np = nc + nd + P;
    if ((int)rows.size() < np + 3) return G;

    vector<double> us;
    for (const auto& r : rows) us.push_back(r.u);
    sort(us.begin(), us.end());
    G.uLo = us.front();  G.uHi = us.back();  G.u0 = us[us.size()/2];

    vector<vector<double>> A(np, vector<double>(np, 0.0));
    vector<double> B(np, 0.0);
    auto addRow = [&](const vector<pair<int,double>>& a, double y) {
        for (const auto& p : a) {
            B[p.first] += p.second*y;
            for (const auto& q : a) A[p.first][q.first] += p.second*q.second;
        }
    };
    auto design = [&](const GlobalRow& r) {
        vector<pair<int,double>> a = {{ci.at(r.ch),1.0},{di.at(r.dk),1.0}};
        double x = r.u - G.u0, p = 1.0;
        for (int q = 0; q < P; ++q) { p *= x; a.push_back({nc+nd+q, p}); }
        return a;
    };
    for (const auto& r : rows) addRow(design(r), r.y);
    vector<pair<int,double>> gauge;
    for (int i = 0; i < nc; ++i) gauge.push_back({i, 1.0});
    addRow(gauge, 0.0);

    vector<double> x;
    if (!SolveLinear(A, B, x)) return G;
    for (const auto& kv : ci) G.c[kv.first] = x[kv.second];
    for (const auto& kv : di) G.D[kv.first] = x[kv.second];
    G.g.assign(x.begin()+nc+nd, x.end());
    double s2 = 0;
    for (const auto& r : rows) {
        double m = 0;
        for (const auto& p : design(r)) m += p.second*x[p.first];
        s2 += (r.y - m)*(r.y - m);
    }
    G.nRows = (int)rows.size();
    G.rms   = sqrt(s2/max(1, G.nRows - np + 1));
    G.ok    = true;
    return G;
}

// ============================================================
void PMTWaveformTree::Loop(const std::string& diffuser_filter)
// ============================================================
{
    if (fChain == 0) return;

    // Timing uses the raw Waveform only; skip the unused derivative branch.
    fChain->SetBranchStatus("*", 1);
    fChain->SetBranchStatus("Derivative", 0);

    const bool filterAll = (diffuser_filter == "ALL");

    cout << "\n=================================================" << endl;
    cout << "=== PMTTimingAnalyzer ===" << endl;
    cout << "  diffuser_filter : " << diffuser_filter << endl;

    // Output directory / ROOT file namespaced by filter
    string outTag   = filterAll ? "ALL" : diffuser_filter;
    string plotDir  = "Plots_" + outTag;
    string rootName = "output_" + outTag + ".root";

    // -------------------------------------------------------
    // Pre-compute theoretical dist/c for relevant diffusers
    // -------------------------------------------------------
    map<string, map<int,double>> theo_abs;
    for (int i = 1; i <= 5; ++i) for (int s = 0; s < 2; ++s) {
        const char* side = s ? "East" : "West";
        string dk = Form("%s_%d", side, i);
        if (!filterAll && dk != diffuser_filter) continue;
        const Point3D& diff = s ? east_diffuser_pos.at(i) : west_diffuser_pos.at(i);
        for (const auto& kv : chanID_to_pos)
            if (chanID_to_delayNS.count(kv.first))
                theo_abs[dk][kv.first] = Dist3(diff, kv.second) / kLightSpeed_cm_ns;
    }

    // -------------------------------------------------------
    // Histogram booking  (residual per channel)
    // -------------------------------------------------------
    map<int,TH1F*> h_res;
    for (const auto& kv : chanID_to_delayNS) {
        int ch = kv.first;
        h_res[ch] = new TH1F(Form("h_res_Ch%d",ch),
            Form("Residual Ch%d;#Deltat_{exp,rel}-#Deltat_{theo,rel} [ns];Counts",ch),
            400,-100,100);
    }
    map<string,map<int,TH1F*>> h_res_diff;
    for (const auto& dkv : theo_abs) {
        const string& dk = dkv.first;
        for (const auto& kv : chanID_to_delayNS)
            h_res_diff[dk][kv.first] = new TH1F(
                Form("h_res_%s_Ch%d",dk.c_str(),kv.first),
                Form("Residual %s Ch%d;#Deltat [ns];Counts",dk.c_str(),kv.first),
                400,-100,100);
    }

    // 2-D Y-Z maps
    const int nZ=18; const double Zlo=20,Zhi=490;
    const int nY=9;  const double Ylo=-200,Yhi=200;
    auto make2D = [&](const char* nm, const char* tit) {
        TH2F* h = new TH2F(nm,tit,nZ,Zlo,Zhi,nY,Ylo,Yhi);
        h->GetXaxis()->SetTitle("Z [cm]"); h->GetYaxis()->SetTitle("Y [cm]");
        h->SetStats(0); return h;
    };
    TH2F* h2_mean_posX  = make2D("h2_res_mean_posX",  "Residual mean - West APA");
    TH2F* h2_mean_negX  = make2D("h2_res_mean_negX",  "Residual mean - East APA");
    TH2F* h2_sigma_posX = make2D("h2_res_sigma_posX", "Residual #sigma - West APA");
    TH2F* h2_sigma_negX = make2D("h2_res_sigma_negX", "Residual #sigma - East APA");
    map<string,TH2F*> h2_mp_d,h2_mn_d,h2_sp_d,h2_sn_d;
    for (const auto& dkv : theo_abs) {
        const string& dk = dkv.first;
        h2_mp_d[dk]=make2D(Form("h2_res_mean_posX_%s",dk.c_str()),Form("Res mean West - %s",dk.c_str()));
        h2_mn_d[dk]=make2D(Form("h2_res_mean_negX_%s",dk.c_str()),Form("Res mean East - %s",dk.c_str()));
        h2_sp_d[dk]=make2D(Form("h2_res_sigma_posX_%s",dk.c_str()),Form("Res #sigma West - %s",dk.c_str()));
        h2_sn_d[dk]=make2D(Form("h2_res_sigma_negX_%s",dk.c_str()),Form("Res #sigma East - %s",dk.c_str()));
    }

    // -------------------------------------------------------
    // PASS 1: per-pulse hits (kTimingMode 0) or averaged waveforms (1)
    // GetEntries() not GetEntriesFast(): the latter returns LLONG_MAX on a
    // TChain that has not been scanned yet.
    // Only the diffuser's own side is used (cross-TPC channels see no
    // direct LED light).
    // -------------------------------------------------------
    Long64_t nentries = fChain->GetEntries();
    cout << "=== PASS 1: reading " << nentries << " entries ===" << endl;

    map<string, vector<LEDHit>>            hits;       // [dk]
    map<string, map<pair<int,int>,int>>    evtIndex;   // [dk][(run,event)] -> index
    long nSeen = 0, nNoisyBase = 0, nNoHit = 0;

    const int nWin = kWfHi - kWfLo;
    map<string,map<int,vector<double>>> wfSum, wfSumHalf[2];
    map<string,map<int,long>>           wfN,   wfNHalf[2];

    for (Long64_t i = 0; i < nentries; ++i) {
        if (i % 50000 == 0)
            cout << "  Entry " << i << " / " << nentries << "\r" << flush;

        // Read Channel/Run first and skip unused channels before unpacking waveforms.
        Long64_t local = LoadTree(i);
        if (local < 0) break;
        b_Channel->GetEntry(local);
        b_Run->GetEntry(local);
        if (!chanID_to_delayNS.count(Channel) || !chanID_to_pos.count(Channel)) continue;
        if (!run_to_diffuser.count(Run)) continue;
        const DiffuserInfo& diR = run_to_diffuser.at(Run);
        string dkR = Form("%s_%d", diR.side.c_str(), diR.index);
        if (!filterAll && dkR != diffuser_filter) continue;
        if (!theo_abs.count(dkR)) continue;
        bool chWest = chanID_to_pos.at(Channel).x > 0;
        if (chWest != (diR.side == "West")) continue;

        fChain->GetEntry(i);
        if (!Waveform || (int)Waveform->size() < max(kHitHi, kWfHi)) continue;
        ++nSeen;

        if (kTimingMode == 1) {
            double base = 0;
            for (int k = 0; k < kNBaselineTicks; ++k) base += (*Waveform)[k];
            base /= kNBaselineTicks;
            vector<double>& acc  = wfSum[dkR][Channel];
            vector<double>& accH = wfSumHalf[Event%2][dkR][Channel];
            if (acc.empty())  acc.assign(nWin, 0.0);
            if (accH.empty()) accH.assign(nWin, 0.0);
            for (int k = 0; k < nWin; ++k) {
                double v = kPulseSign*((*Waveform)[kWfLo+k] - base);
                acc[k] += v; accH[k] += v;
            }
            wfN[dkR][Channel]++; wfNHalf[Event%2][dkR][Channel]++;
            continue;
        }

        PulseHit ph = FindFirstPulse(*Waveform);
        if (!ph.ok) { if (ph.rms > kMaxBaselineRMS) ++nNoisyBase; else ++nNoHit; continue; }

        map<pair<int,int>,int>& eidx = evtIndex[dkR];
        pair<int,int> key = {Run, Event};
        auto eit = eidx.find(key);
        int ie;
        if (eit == eidx.end()) { ie = (int)eidx.size(); eidx[key] = ie; }
        else ie = eit->second;

        double t = kNsPerTick*(StartTick_Waveform + ph.tTick) + chanID_to_delayNS.at(Channel);
        hits[dkR].push_back({Channel, ie, (float)t, (float)log(max(ph.charge, 1.0))});
    }
    cout << "\n  waveforms used: " << nSeen;
    if (kTimingMode == 0)
        cout << "   no hit: " << nNoHit << "   noisy baseline: " << nNoisyBase;
    cout << endl;

    struct CPoint { double mean, err; };
    map<string,map<int,CPoint>> chTime;      // [dk][ch] -> time used by the radius plots [ns]
    map<string,DiffSolution>    solution;    // [dk] per-pulse solution (mode 0)
    map<string,GlobalSize>      gsize;       // ["West"/"East"] global size correction

    // Size-correction helpers (mode 0).  With the global fit:
    //   channel time  = rawMean - Wg(median lnQ)
    //   flash time    = t - o_evt - W(lnQ) + W(median lnQ) - Wg(median lnQ)
    //                   (within-channel W only removes flash-to-flash jitter)
    // Without it (single diffuser): within-channel W only.
    auto sideOf    = [](const string& dk) { return dk.substr(0,4); };
    auto useGlobal = [&](const string& dk) {
        auto it = gsize.find(sideOf(dk)); return it != gsize.end() && it->second.ok; };
    auto Wcurve    = [&](const string& dk, double u) {
        return useGlobal(dk) ? GlobalW(gsize.at(sideOf(dk)), u) : WalkEval(solution.at(dk), u); };
    auto chanCorr  = [&](const string& dk, int ch) {
        const ChanTime& c = solution.at(dk).ch.at(ch);
        return useGlobal(dk) ? c.rawMean - Wcurve(dk, c.meanU) : c.mean; };
    auto flashCorr = [&](const string& dk, const LEDHit& h) {
        const DiffSolution& S = solution.at(dk);
        double t = h.t - S.evtOff[h.evt] - WalkEval(S, h.u);
        if (useGlobal(dk)) {
            double U = S.ch.at(h.ch).meanU;
            t += WalkEval(S, U) - Wcurve(dk, U);
        }
        return t;
    };

    // -------------------------------------------------------
    // PER-PULSE: solve T_ch + o_evt + W(lnQ) per diffuser, then fill
    // the per-flash residual histograms with the corrected times.
    // -------------------------------------------------------
    if (kTimingMode == 0) {
        cout << "=== Per-pulse timing: T_ch + o_evt + W(lnQ) per diffuser ===" << endl;
        for (const auto& dkv : hits) {
            const string& dk = dkv.first;
            int nEvt = (int)evtIndex[dk].size();
            DiffSolution S = SolveDiffuser(dkv.second, nEvt);
            int nOK = 0; for (char c : S.evtOK) nOK += c;
            cout << "  " << dk << ": " << dkv.second.size() << " hits, "
                 << nOK << "/" << nEvt << " flashes with >=" << kMinGoodCh
                 << " hits, " << S.ch.size() << " channels" << endl;
            if (!S.walk.empty()) {
                cout << "    walk W(lnQ) coeffs (u0=" << Form("%.2f",S.u0) << "):";
                for (double a : S.walk) cout << " " << Form("%.3g",a);
                cout << Form("   W(lnQ_lo=%.1f)=%.1f ns, W(lnQ_hi=%.1f)=%.1f ns",
                             S.uLo, WalkEval(S,S.uLo), S.uHi, WalkEval(S,S.uHi)) << endl;
            }
            for (const auto& c : S.ch) {
                cout << Form("    Ch%-4d T=%9.2f +- %5.2f ns  sigma=%6.1f  raw=%9.2f  med lnQ=%5.2f  hitFrac=%.2f  N=%d",
                             c.first, c.second.mean, c.second.err, c.second.sigma,
                             c.second.rawMean, c.second.meanU, c.second.hitFrac, c.second.n) << endl;
            }
            if (S.ch.empty()) continue;
            solution[dk] = S;
        }

        // ---------------------------------------------------
        // STAGE 2: global size correction, per side (ALL mode)
        // ---------------------------------------------------
        if (kUseGlobalSizeCorr) {
            for (const char* sd : {"West","East"}) {
                vector<GlobalRow> rows;
                for (const auto& kv : solution) {
                    if (sideOf(kv.first) != sd) continue;
                    for (const auto& c : kv.second.ch)
                        rows.push_back({kv.first, c.first, c.second.meanU,
                                        c.second.rawMean - theo_abs.at(kv.first).at(c.first)});
                }
                if (rows.empty()) continue;
                GlobalSize G = FitGlobalSize(rows, kGlobalWalkOrder);
                if (!G.ok) {
                    cout << "  [size] " << sd << ": global size fit not possible (needs >=2 diffusers of this side,"
                         << " i.e. diffuser_filter = \"ALL\") -> using the within-channel W only" << endl;
                    continue;
                }
                gsize[sd] = G;
                cout << "=== Global size correction " << sd << ": " << G.nRows << " (channel,diffuser) points, "
                     << G.c.size() << " channels, " << G.D.size() << " diffusers ===" << endl;
                cout << "    Wg coeffs (u0=" << Form("%.2f",G.u0) << "):";
                for (double a : G.g) cout << " " << Form("%.4g",a);
                cout << Form("\n    Wg(log10Q=%.2f) = %.1f ns ... Wg(log10Q=%.2f) = %.1f ns,  fit rms = %.1f ns",
                             G.uLo/log(10.), GlobalW(G,G.uLo), G.uHi/log(10.), GlobalW(G,G.uHi), G.rms) << endl;
                for (const auto& d : G.D) cout << Form("    D[%s] = %.1f ns", d.first.c_str(), d.second) << endl;
            }
        }

        // channel times used by the radius plots
        for (const auto& kv : solution)
            for (const auto& c : kv.second.ch)
                chTime[kv.first][c.first] = {chanCorr(kv.first, c.first),
                                             useGlobal(kv.first) ? c.second.rawErr : c.second.err};

        // per-flash residuals (same correction as the channel times)
        for (const auto& kv : solution) {
            const string& dk = kv.first;
            const DiffSolution& S = kv.second;
            double sT = 0, sM = 0;
            for (const auto& c : S.ch) { sT += theo_abs.at(dk).at(c.first); sM += chTime[dk][c.first].mean; }
            double avgT = sT/S.ch.size(), avgM = sM/S.ch.size();
            for (const auto& h : hits.at(dk)) {
                if (!S.evtOK[h.evt] || !S.ch.count(h.ch)) continue;
                double res = (flashCorr(dk, h) - avgM) - (theo_abs.at(dk).at(h.ch) - avgT);
                h_res[h.ch]->Fill(res);
                if (h_res_diff.count(dk)) h_res_diff[dk][h.ch]->Fill(res);
            }
        }
    }

    // -------------------------------------------------------
    // AVERAGED WAVEFORM (kTimingMode 1, cross-check only)
    // t = kNsPerTick*(kWfLo + CF tick) + cable delay   [ns, absolute]
    // -------------------------------------------------------
    if (kTimingMode == 1) {
        cout << "=== Average-waveform CF timing (frac=" << kCFfrac << ") ===" << endl;
        for (const auto& dkv : wfSum) {
            const string& dk = dkv.first;
            for (const auto& chv : dkv.second) {
                int ch = chv.first; long n = wfN[dk][ch];
                if (n < 10) continue;
                vector<double> avg(chv.second);
                for (double& x : avg) x /= n;
                double pk = 0;
                double tk = CFTimeTicks(avg, kCFfrac, &pk);
                if (tk < 0) continue;
                double t = kNsPerTick*(kWfLo + tk) + chanID_to_delayNS.at(ch);
                double th[2]; bool okH = true;
                for (int h = 0; h < 2; ++h) {
                    long nh = wfNHalf[h][dk][ch];
                    if (nh < 5) { okH = false; break; }
                    vector<double> ah(wfSumHalf[h][dk][ch]);
                    for (double& x : ah) x /= nh;
                    double tkh = CFTimeTicks(ah, kCFfrac);
                    if (tkh < 0) { okH = false; break; }
                    th[h] = kNsPerTick*(kWfLo + tkh);
                }
                double err = okH ? max(0.5*fabs(th[0]-th[1]), 0.5) : 50.0;
                chTime[dk][ch] = {t, err};
            }
        }
    }

    // -------------------------------------------------------
    // Gaussian fits -> 2-D Y-Z maps
    // -------------------------------------------------------
    cout << "=== Fitting 2-D maps ===" << endl;
    // Channels that share a map bin are AVERAGED (TH2::Fill would add them up;
    // XA boxes hold 2 channels at the same y,z).  mapCh keeps (channel,value)
    // per bin for the labels on the res_sigma maps.
    map<TH2F*, map<int,pair<double,int>>>            mapAcc;
    map<TH2F*, map<int,vector<pair<int,double>>>>    mapCh;
    auto setAvg = [&](TH2F* hh, double z, double y, double v, int ch) {
        int b = hh->FindBin(z,y);
        auto& a = mapAcc[hh][b]; a.first += v; a.second++;
        hh->SetBinContent(b, a.first/a.second);
        mapCh[hh][b].push_back({ch, v});
    };
    auto fitAndFill = [&](TH1F* h, int ch,
                          TH2F* hMP, TH2F* hMN, TH2F* hSP, TH2F* hSN) {
        if (!h||h->GetEntries()<10) return;
        TF1 fg("fg","gaus",-50,50); h->Fit(&fg,"RQ0");
        double mu=fg.GetParameter(1), sig=fabs(fg.GetParameter(2));
        if (!sig) return;
        double z=chanID_to_pos[ch].z, y=chanID_to_pos[ch].y, x=chanID_to_pos[ch].x;
        if (x>0){setAvg(hMP,z,y,mu,ch);setAvg(hSP,z,y,sig,ch);}
        else    {setAvg(hMN,z,y,mu,ch);setAvg(hSN,z,y,sig,ch);}
    };
    for (const auto& kv : h_res)
        fitAndFill(kv.second,kv.first,h2_mean_posX,h2_mean_negX,h2_sigma_posX,h2_sigma_negX);
    for (const auto& dkv : h_res_diff)
        for (const auto& kv : dkv.second)
            fitAndFill(kv.second,kv.first,
                       h2_mp_d.at(dkv.first),h2_mn_d.at(dkv.first),
                       h2_sp_d.at(dkv.first),h2_sn_d.at(dkv.first));

    // -------------------------------------------------------
    // Open output file, set style
    // -------------------------------------------------------
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(plotDir.c_str(), kTRUE);
    TFile* outfile = new TFile(rootName.c_str(), "RECREATE");

    auto saveColz = [&](TH2F* h, const char* nm, const char* zt, bool savePng = true, bool isSigma = false) {
        if (!h) return;
        h->GetZaxis()->SetTitle(zt); h->GetZaxis()->SetTitleOffset(1.5);
        h->Write();                                   // unmodified values in the ROOT file
        if (!savePng || mapAcc[h].empty()) return;   // nothing to show (empty map)
        TCanvas c("ctmp","",1000,800); c.SetRightMargin(0.18);
        string note;
        if (isSigma) {
            vector<double> v; for (const auto& kv : mapAcc[h]) v.push_back(kv.second.first/kv.second.second);
            if (!v.empty()) {
                sort(v.begin(), v.end());
                double cap = kSigmaCapFactor * v[v.size()/2];
                vector<int> over;
                for (const auto& kv : mapCh[h]) for (const auto& cv : kv.second) if (cv.second > cap) over.push_back(cv.first);
                for (const auto& kv : mapAcc[h])
                    if (h->GetBinContent(kv.first) > cap) h->SetBinContent(kv.first, cap);
                h->SetMinimum(1e-3); h->SetMaximum(cap);
                note = Form("colour scale capped at %.0f ns (%.0f#times median)", cap, kSigmaCapFactor);
                if (!over.empty()) {
                    note += "; above cap: ch";
                    for (int ch : over) note += Form(" %d", ch);
                }
            }
        }
        h->Draw("COLZ");
        if (!note.empty()) {
            TLatex t; t.SetNDC(); t.SetTextSize(0.026); t.SetTextColor(kGray+3);
            t.DrawLatex(0.11, 0.035, note.c_str());
        }
        c.SaveAs(Form("%s/%s.png",plotDir.c_str(),nm));
    };
    saveColz(h2_mean_posX, "res_mean_posX_All",  "Residual mean [ns]");
    saveColz(h2_mean_negX, "res_mean_negX_All",  "Residual mean [ns]");
    saveColz(h2_sigma_posX,"res_sigma_posX_All", "Residual #sigma [ns]", kSaveExtraPlots, true);
    saveColz(h2_sigma_negX,"res_sigma_negX_All", "Residual #sigma [ns]", kSaveExtraPlots, true);
    for (const auto& dkv : theo_abs) {
        const string& dk=dkv.first;
        // same-side matching: West diffusers only see posX channels, East only
        // negX, so the opposite-side map would always be empty -> not drawn
        if (dk.substr(0,4)=="West") {
            saveColz(h2_mp_d[dk],Form("res_mean_posX_%s", dk.c_str()),"Residual mean [ns]", kSaveExtraPlots);
            saveColz(h2_sp_d[dk],Form("res_sigma_posX_%s",dk.c_str()),"Residual #sigma [ns]", true, true);
        } else {
            saveColz(h2_mn_d[dk],Form("res_mean_negX_%s", dk.c_str()),"Residual mean [ns]", kSaveExtraPlots);
            saveColz(h2_sn_d[dk],Form("res_sigma_negX_%s",dk.c_str()),"Residual #sigma [ns]", true, true);
        }
    }
    for (const auto& kv : h_res)
        if (kv.second->GetEntries()>0) kv.second->Write();
    for (const auto& dkv : h_res_diff)
        for (const auto& kv : dkv.second)
            if (kv.second->GetEntries()>0) kv.second->Write();

    // Residual mean/sigma vs channel
    {
        vector<double> ch_v,mu_v,mu_e,sig_v;
        for (const auto& kv : h_res) {
            int ch=kv.first; TH1F* h=kv.second;
            if (!h||h->GetEntries()<10) continue;
            TF1 fg("fg","gaus",-50,50); h->Fit(&fg,"RQ0");
            double mu=fg.GetParameter(1),sig=fabs(fg.GetParameter(2)),n=h->GetEntries();
            ch_v.push_back(ch); mu_v.push_back(mu);
            mu_e.push_back(sig/sqrt(n)); sig_v.push_back(sig);
        }
        if ((int)ch_v.size()>2) {
            TGraphErrors* grMu=new TGraphErrors((int)ch_v.size(),ch_v.data(),
                                mu_v.data(),nullptr,mu_e.data());
            grMu->SetName("gr_res_mean_vs_ch");
            grMu->SetTitle("Residual mean vs channel;Channel ID;"
                           "Mean(#Deltat_{exp,rel}-#Deltat_{theo,rel}) [ns]");
            grMu->SetMarkerStyle(20);grMu->SetMarkerColor(kBlack);grMu->SetMarkerSize(0.8);
            TGraph* grSig=new TGraph((int)ch_v.size(),ch_v.data(),sig_v.data());
            grSig->SetName("gr_res_sigma_vs_ch");
            grSig->SetTitle("Residual #sigma vs channel;Channel ID;#sigma [ns]");
            grSig->SetMarkerStyle(20);grSig->SetMarkerColor(kBlack);grSig->SetMarkerSize(0.8);
            TF1 fz("fz","0",0,400);fz.SetLineColor(kRed);fz.SetLineStyle(2);
            TCanvas c1("c_mu","",1200,600);
            c1.SetLeftMargin(0.10);c1.SetBottomMargin(0.12);
            grMu->Draw("AP"); fz.Draw("same");
            LabelOutlierChannels(ch_v, mu_v, grMu->GetHistogram()->GetMinimum(), grMu->GetHistogram()->GetMaximum());
            if (kSaveExtraPlots) c1.SaveAs(Form("%s/residual_mean_vs_channel.png",plotDir.c_str()));
            c1.Write("c_residual_mean_vs_channel"); grMu->Write();
            TCanvas c2("c_sig","",1200,600);
            c2.SetLeftMargin(0.10);c2.SetBottomMargin(0.12);
            grSig->Draw("AP");
            LabelOutlierChannels(ch_v, sig_v, grSig->GetHistogram()->GetMinimum(), grSig->GetHistogram()->GetMaximum());
            c2.SaveAs(Form("%s/residual_sigma_vs_channel.png",plotDir.c_str()));
            c2.Write("c_residual_sigma_vs_channel"); grSig->Write();
        }
    }

    // ============================================================
    // SIGNAL-SIZE DIAGNOSTICS  (per-pulse mode only)
    //  size_vs_distance_<dk> : channel signal size (median log10 Q, bar =
    //                          16-84% flash-to-flash spread) vs distance
    //  size_vs_delay_<dk>    : signal size (y) vs delay (x)
    //                          left  = BEFORE size correction, with the
    //                                  size-correction curve W
    //                          right = AFTER size correction; the channel
    //                                  points are exactly the plot-B values
    //  *_ALL                 : channel points of every diffuser
    //  delay = channel time - r/v - diffuser average  (light travel removed)
    //  All diffusers share the same axis ranges.
    // ============================================================
    if (kTimingMode == 0 && !solution.empty()) {
        cout << "=== Signal-size diagnostics ===" << endl;
        const double LN10 = log(10.0);
        int dcol[] = {kRed+1,kBlue+1,kGreen+2,kOrange+7,kMagenta+1,
                      kCyan+2,kYellow+3,kViolet+2,kTeal+3,kPink+6};
        struct ChSum { string dk; int ch; double r, lq, slq, before, after, err; };
        vector<ChSum> allCh;
        map<string,double> Dd;   // per-diffuser average of (corrected time - r/v)

        for (const auto& kv : solution) {
            const string& dk = kv.first;
            if (kv.second.ch.empty()) continue;
            const Point3D& diff = (dk.substr(0,4)=="West")
                ? west_diffuser_pos.at(stoi(dk.substr(5)))
                : east_diffuser_pos.at(stoi(dk.substr(5)));
            double s = 0;
            for (const auto& c : kv.second.ch) s += chanCorr(dk, c.first) - theo_abs.at(dk).at(c.first);
            Dd[dk] = s/kv.second.ch.size();
            for (const auto& c : kv.second.ch) {
                double after = chanCorr(dk, c.first) - theo_abs.at(dk).at(c.first) - Dd[dk];
                allCh.push_back({dk, c.first, Dist3(diff, chanID_to_pos.at(c.first)),
                                 c.second.meanU/LN10, c.second.sdU/LN10,
                                 after + Wcurve(dk, c.second.meanU), after, chTime[dk][c.first].err});
            }
        }
        if (!allCh.empty()) {
        // ---- common axis ranges ----
        double rLo=1e9,rHi=-1e9,qLo=1e9,qHi=-1e9,dLo=1e9,dHi=-1e9;
        for (const auto& s : allCh) {
            rLo=min(rLo,s.r); rHi=max(rHi,s.r);
            qLo=min(qLo,s.lq-s.slq); qHi=max(qHi,s.lq+s.slq);
            dLo=min({dLo,s.before,s.after}); dHi=max({dHi,s.before,s.after});
        }
        rLo -= 15; rHi += 15;
        qLo = floor(qLo*2)/2 - 0.5;  qHi = ceil(qHi*2)/2 + 0.5;
        dLo = 20*floor(dLo/20) - 40;  dHi = 20*ceil(dHi/20) + 40;
        const char* kSizeT  = "Signal size  log_{10}(Q / ADC#upointns)";
        const char* kBeforeT = "Delay: t - r/v - #LTt - r/v#GT [ns]";
        const char* kAfterT  = "Delay: #Deltat_{exp,rel} - #Deltat_{theo,rel} [ns]";

        auto typeGraph = [&](const vector<const ChSum*>& v, bool typeB,
                             function<double(const ChSum&)> X, function<double(const ChSum&)> Y,
                             function<double(const ChSum&)> EX, function<double(const ChSum&)> EY,
                             const string& name) -> TGraphErrors* {
            vector<double> x,y,ex,ey;
            for (const ChSum* s : v) if (IsTypeB(s->ch) == typeB) {
                x.push_back(X(*s)); y.push_back(Y(*s)); ex.push_back(EX(*s)); ey.push_back(EY(*s)); }
            if (x.empty()) return nullptr;
            TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
            g->SetName(name.c_str());
            g->SetMarkerStyle(typeB ? 25 : 20); g->SetMarkerSize(typeB ? 1.1 : 0.9);
            int col = typeB ? kRed+1 : kBlack;
            g->SetMarkerColor(col); g->SetLineColor(col);
            return g;
        };
        auto zero = [](const ChSum&) { return 0.0; };

        for (const auto& kv : solution) {
            const string& dk = kv.first;
            if (!Dd.count(dk)) continue;
            const DiffSolution& S = kv.second;
            vector<const ChSum*> mine;
            for (const auto& s : allCh) if (s.dk == dk) mine.push_back(&s);

            // ---------- size vs distance ----------
            TCanvas* c1 = new TCanvas(Form("c_size_vs_distance_%s",dk.c_str()),"",900,650);
            c1->SetLeftMargin(0.12); c1->SetRightMargin(0.04);
            TH2F* f1 = new TH2F(Form("hfr_size_vs_distance_%s",dk.c_str()),
                Form("%s: signal size vs distance;Distance from diffuser [cm];%s",dk.c_str(),kSizeT),
                1,rLo,rHi, 1,qLo,qHi);
            f1->SetStats(0); f1->Draw();
            TLegend* l1 = new TLegend(0.55,0.78,0.94,0.89);
            l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.035);
            for (bool tb : {false,true}) {
                TGraphErrors* g = typeGraph(mine, tb, [](const ChSum& s){return s.r;},
                                            [](const ChSum& s){return s.lq;}, zero,
                                            [](const ChSum& s){return s.slq;},
                                            Form("gr_size_vs_distance_%s_%s",dk.c_str(), tb?"B":"A"));
                if (!g) continue;
                g->Draw("P same"); g->Write();
                l1->AddEntry(g, tb ? kTypeNameB : kTypeNameA, "p");
            }
            l1->Draw();
            TLatex n1; n1.SetNDC(); n1.SetTextSize(0.028); n1.SetTextColor(kGray+2);
            n1.DrawLatex(0.14,0.15,"point = channel median, bar = 16-84% of its flashes");
            c1->SaveAs(Form("%s/size_vs_distance_%s.png",plotDir.c_str(),dk.c_str()));
            c1->Write(); delete c1;

            // ---------- size vs delay: before | after ----------
            TH2F* hB = new TH2F(Form("h2_size_vs_delay_before_%s",dk.c_str()),
                Form("%s: before size correction;%s;%s;Flashes",dk.c_str(),kBeforeT,kSizeT),
                90,dLo,dHi, 60,qLo,qHi);
            TH2F* hA = new TH2F(Form("h2_size_vs_delay_after_%s",dk.c_str()),
                Form("%s: after size correction;%s;%s;Flashes",dk.c_str(),kAfterT,kSizeT),
                90,dLo,dHi, 60,qLo,qHi);
            for (TH2F* hh : {hB,hA}) { hh->SetStats(0); hh->GetZaxis()->SetTitleOffset(1.3); }
            for (const auto& h : hits.at(dk)) {
                if (!S.evtOK[h.evt] || !S.ch.count(h.ch)) continue;
                double tr = theo_abs.at(dk).at(h.ch) + Dd.at(dk);
                hB->Fill(h.t - S.evtOff[h.evt] - tr, h.u/LN10);
                hA->Fill(flashCorr(dk, h) - tr, h.u/LN10);
            }
            TCanvas* c2 = new TCanvas(Form("c_size_vs_delay_%s",dk.c_str()),"",1500,620);
            c2->Divide(2,1);
            for (int pad = 1; pad <= 2; ++pad) {
                c2->cd(pad);
                gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.14); gPad->SetBottomMargin(0.12); gPad->SetLogz();
                TH2F* hh = (pad==1) ? hB : hA;
                hh->GetXaxis()->SetTitleSize(0.040); hh->GetYaxis()->SetTitleSize(0.042);
                hh->Draw("COLZ");
                TLegend* lg = new TLegend(0.37,0.76,0.78,0.89);   // top right, clear of the colour bar
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.034);
                for (bool tb : {false,true}) {
                    TGraphErrors* g = typeGraph(mine, tb,
                        pad==1 ? function<double(const ChSum&)>([](const ChSum& s){return s.before;})
                               : function<double(const ChSum&)>([](const ChSum& s){return s.after;}),
                        [](const ChSum& s){return s.lq;}, [](const ChSum& s){return s.err;}, zero,
                        Form("gr_size_vs_delay_%s_%s_%s", pad==1?"before":"after", dk.c_str(), tb?"B":"A"));
                    if (!g) continue;
                    g->Draw("P same"); g->Write();
                    lg->AddEntry(g, Form("channel average, %s", tb ? kTypeNameB : kTypeNameA), "p");
                }
                if (pad == 1) {
                    vector<double> wx, wy;
                    for (int k = 0; k <= 200; ++k) {
                        double q = qLo + (qHi-qLo)*k/200.0;
                        wx.push_back(Wcurve(dk, q*LN10)); wy.push_back(q);
                    }
                    TGraph* gW = new TGraph((int)wx.size(), wx.data(), wy.data());
                    gW->SetName(Form("gr_size_corr_curve_%s",dk.c_str()));
                    gW->SetLineColor(kMagenta+1); gW->SetLineWidth(3);
                    gW->Draw("L same"); gW->Write();
                    lg->AddEntry(gW, "size correction", "l");
                } else {
                    TLine* l0 = new TLine(0,qLo,0,qHi); l0->SetLineStyle(2); l0->SetLineWidth(2); l0->Draw();
                    lg->AddEntry(l0, "expected: 0", "l");
                }
                lg->Draw();
            }
            c2->SaveAs(Form("%s/size_vs_delay_%s.png",plotDir.c_str(),dk.c_str()));
            c2->Write(); hB->Write(); hA->Write();
            delete c2;

            // ---------- size map: median signal size of every channel in Z-Y ----------
            // The brightest channels sit next to the diffuser -> checks the
            // run -> diffuser mapping.  Channels with no usable light: open circles.
            {
                const bool isW = (sideOf(dk) == "West");
                const map<int,Point3D>& dmap = isW ? west_diffuser_pos : east_diffuser_pos;
                const int dIdx = stoi(dk.substr(dk.find('_')+1));
                TH2F* hm = new TH2F(Form("h2_size_map_%s",dk.c_str()),
                    Form("%s: signal size per channel (star = this diffuser);Z [cm];Y [cm]",dk.c_str()),
                    34,0,510, 23,-230,230);
                hm->SetStats(0);
                hm->GetZaxis()->SetTitle(Form("median %s", kSizeT)); hm->GetZaxis()->SetTitleOffset(1.4);
                map<int,pair<double,int>> acc;
                set<int> seen;
                for (const ChSum* s : mine) {
                    const Point3D& p = chanID_to_pos.at(s->ch);
                    int b = hm->FindBin(p.z, p.y);
                    acc[b].first += s->lq; acc[b].second++;
                    hm->SetBinContent(b, acc[b].first/acc[b].second);
                    seen.insert(s->ch);
                }
                hm->SetMinimum(qLo); hm->SetMaximum(qHi);
                TCanvas* cm = new TCanvas(Form("c_size_map_%s",dk.c_str()),"",1000,800);
                cm->SetLeftMargin(0.11); cm->SetRightMargin(0.17);
                hm->Draw("COLZ");
                vector<double> uz, uy;
                for (const auto& kv : chanID_to_pos)
                    if ((kv.second.x > 0) == isW && !seen.count(kv.first)) { uz.push_back(kv.second.z); uy.push_back(kv.second.y); }
                if (!uz.empty()) {
                    TGraph* gu = new TGraph((int)uz.size(), uz.data(), uy.data());
                    gu->SetMarkerStyle(24); gu->SetMarkerColor(kGray+1); gu->SetMarkerSize(1.6); gu->Draw("P same");
                }
                TLatex lt; lt.SetTextSize(0.035); lt.SetTextAlign(12); lt.SetTextFont(62);
                for (const auto& d : dmap) {
                    const Point3D& dp = d.second;
                    bool me = (d.first == dIdx);
                    TGraph* gd = new TGraph(1, &dp.z, &dp.y);
                    gd->SetMarkerStyle(me ? 29 : 30); gd->SetMarkerSize(me ? 3.5 : 2.5);
                    gd->SetMarkerColor(me ? kRed+1 : kGray+2); gd->Draw("P same");
                    lt.SetTextColor(me ? kRed+1 : kGray+2);
                    lt.DrawLatex(dp.z+12, dp.y, Form("%d", d.first));
                }
                cm->SaveAs(Form("%s/size_map_%s.png",plotDir.c_str(),dk.c_str()));
                cm->Write(); hm->Write();
                delete cm;
            }
        }

        // ---------- all diffusers ----------
        auto drawAll = [&](const char* fname, const char* xt, double xlo, double xhi, int mode) {
            // mode 0: size vs distance, 1: size vs delay before, 2: size vs delay after
            TCanvas* c = new TCanvas(Form("c_%s",fname),"",1000,700);
            c->SetLeftMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.16);
            TH2F* fr = new TH2F(Form("hfr_%s",fname), Form(";%s;%s",xt,kSizeT), 1,xlo,xhi, 1,qLo,qHi);
            fr->SetStats(0); fr->Draw();
            TLegend* lg = new TLegend(0.10,0.86,0.98,0.99); lg->SetNColumns(5);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.03);
            int ic = 0;
            for (const auto& kv : solution) {
                vector<double> xs, ys, ex, ey;
                for (const auto& s : allCh) if (s.dk == kv.first) {
                    xs.push_back(mode==0 ? s.r : (mode==1 ? s.before : s.after));
                    ys.push_back(s.lq);
                    ex.push_back(mode==0 ? 0.0 : s.err);
                    ey.push_back(mode==0 ? s.slq : 0.0);
                }
                if (xs.empty()) { ++ic; continue; }
                TGraphErrors* g = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), ex.data(), ey.data());
                g->SetName(Form("gr_%s_%s",fname,kv.first.c_str()));
                g->SetMarkerStyle(20+(ic%5)); g->SetMarkerColor(dcol[ic%10]); g->SetLineColor(dcol[ic%10]);
                g->SetMarkerSize(0.9); g->Draw("P same"); g->Write();
                lg->AddEntry(g, kv.first.c_str(), "p");
                ++ic;
            }
            if (mode == 1) {
                int ls = 1;
                for (const char* sd : {"West","East"}) {
                    string dk0;
                    for (const auto& kv : solution) if (sideOf(kv.first) == sd && Dd.count(kv.first)) { dk0 = kv.first; break; }
                    if (dk0.empty() || !useGlobal(dk0)) continue;
                    vector<double> wx, wy;
                    for (int k = 0; k <= 200; ++k) {
                        double q = qLo + (qHi-qLo)*k/200.0;
                        wx.push_back(Wcurve(dk0, q*LN10)); wy.push_back(q);
                    }
                    TGraph* gW = new TGraph((int)wx.size(), wx.data(), wy.data());
                    gW->SetName(Form("gr_size_corr_curve_%s",sd));
                    gW->SetLineColor(kGray+2); gW->SetLineWidth(3); gW->SetLineStyle(ls++);
                    gW->Draw("L same"); gW->Write();
                    lg->AddEntry(gW, Form("size corr. %s", sd), "l");
                }
            }
            if (mode == 2) { TLine* l0 = new TLine(0,qLo,0,qHi); l0->SetLineStyle(2); l0->Draw(); }
            lg->Draw();
            c->SaveAs(Form("%s/%s.png",plotDir.c_str(),fname));
            c->Write();
            delete c;
        };
        drawAll("size_vs_distance_ALL", "Distance from diffuser [cm]", rLo, rHi, 0);
        drawAll("size_vs_delay_ALL", kBeforeT, dLo, dHi, 1);
        drawAll("size_vs_delay_corrected_ALL", kAfterT, dLo, dHi, 2);
        }
    }

    // ============================================================
    // PMT POINT STRUCTS
    // ============================================================
    struct PMTPoint {
        int    ch;
        double r;          // 3-D distance from diffuser [cm]
        double theo_rel;   // dist/c - mean(dist/c) [ns]
        double meas_mean;  // truncated mean (A) or residual mean (B) [ns]
        double meas_err;
    };

    // ============================================================
    // BUILD DATA FOR RADIUS PLOTS
    // ============================================================
    cout << "=== Building radius plots ===" << endl;

    // Fill PMTPoint vectors for one APA half + one diffuser
    auto fillHalf = [&](
        const string& dk,
        bool wantPosX,
        bool usePlotA,   // true->Plot A (colleague), false->Plot B (residual)
        vector<PMTPoint>& out)
    {
        out.clear();
        const Point3D& diff = [&]() -> const Point3D& {
            // parse side from dk
            if (dk.substr(0,4)=="West") return west_diffuser_pos.at(stoi(dk.substr(5)));
            else                        return east_diffuser_pos.at(stoi(dk.substr(5)));
        }();

        if (!chTime.count(dk)) return;
        const map<int,CPoint>& ct = chTime.at(dk);
        vector<int> chsW;
        for (const auto& kv : ct)
            if ((chanID_to_pos.at(kv.first).x>0)==wantPosX && theo_abs.at(dk).count(kv.first))
                chsW.push_back(kv.first);
        if (chsW.empty()) return;
        double sT=0, sM=0;
        for (int ch:chsW){ sT+=theo_abs.at(dk).at(ch); sM+=ct.at(ch).mean; }
        double avgT=sT/chsW.size(), avgM=sM/chsW.size();
        for (int ch:chsW) {
            double r        = Dist3(diff, chanID_to_pos.at(ch));
            double theo_rel = theo_abs.at(dk).at(ch) - avgT;
            double meas_rel = ct.at(ch).mean - avgM;
            double mean     = usePlotA ? meas_rel : meas_rel - theo_rel;
            out.push_back({ch, r, theo_rel, mean, ct.at(ch).err});
        }
        sort(out.begin(),out.end(),[](const PMTPoint& a,const PMTPoint& b){return a.r<b.r;});
    };

    // ============================================================
    // DRAW: radius canvas (one diffuser, A or B)
    // Returns {p0,p1} of the best-fit line to each APA half.
    //
    // Updated: draws a single-panel canvas when only one APA half
    // has data (the matching side).  The fixedYRange flag forces
    // the y-axis to [kRadiusYmin, kRadiusYmax].
    // ============================================================
    struct FitResult { double p0, p1; };   // meas = p0 + p1*r

    auto drawRadiusCanvas = [&](
        const vector<PMTPoint>& pts_neg,
        const vector<PMTPoint>& pts_pos,
        const string& dk, int idx, const char* side,
        const char* filePrefix, const char* canvPrefix,
        const char* yAxisTitle,
        bool fixedYRange,
        bool residualPlot,      // plot B: expectation is 0
        FitResult& fitPos, FitResult& fitNeg)
    {
        fitPos={0,0}; fitNeg={0,0};
        if (pts_neg.empty()&&pts_pos.empty()) return;

        // Determine which panels to draw
        bool hasPos = !pts_pos.empty();
        bool hasNeg = !pts_neg.empty();
        int  nPads  = (hasPos && hasNeg) ? 2 : 1;

        double rMin=1e9,rMax=-1e9,yMin=1e9,yMax=-1e9;
        for (const auto* pp:{&pts_neg,&pts_pos})
            for (const auto& p:*pp){
                rMin=min(rMin,p.r); rMax=max(rMax,p.r);
                yMin=min({yMin,p.theo_rel,p.meas_mean-p.meas_err});
                yMax=max({yMax,p.theo_rel,p.meas_mean+p.meas_err});
            }
        if (rMax<=rMin){rMin-=10;rMax+=10;}
        if (yMax<=yMin){yMin-=5;yMax+=5;}
        double rP=0.05*(rMax-rMin), yP=0.12*(yMax-yMin);
        rMin-=rP;rMax+=rP;yMin-=yP;yMax+=yP;

        // Override y-range for dtexp plots
        if (fixedYRange) { yMin = kRadiusYmin; yMax = kRadiusYmax; }

        int canvW = (nPads==2) ? 1400 : 800;
        TCanvas* cRad=new TCanvas(Form("%s_%s",canvPrefix,dk.c_str()),
            Form("%s  %s %d",filePrefix,side,idx), canvW, 620);
        if (nPads==2) cRad->Divide(2,1,0.005,0.01);

        // Build ordered list of panels to draw
        struct PanelJob {
            int padIdx;
            const vector<PMTPoint>* pts;
            FitResult* fit;
            bool isPosX;
        };
        vector<PanelJob> panels;
        if (nPads==2) {
            panels.push_back({1, &pts_pos, &fitPos, true});
            panels.push_back({2, &pts_neg, &fitNeg, false});
        } else if (hasPos) {
            panels.push_back({1, &pts_pos, &fitPos, true});
        } else {
            panels.push_back({1, &pts_neg, &fitNeg, false});
        }

        for (const auto& panel : panels) {
            const vector<PMTPoint>& pts = *panel.pts;
            FitResult& fitOut           = *panel.fit;
            const char* apaLabel = panel.isPosX
                ? "Positive-X (West)" : "Negative-X (East)";
            const char* padTitle = Form("%s APA  |  Diffuser %s %d",
                                        apaLabel, side, idx);

            if (nPads==2) cRad->cd(panel.padIdx);
            gPad->SetLeftMargin(0.14);gPad->SetBottomMargin(0.13);
            gPad->SetRightMargin(0.04);gPad->SetTopMargin(0.10);

            if (pts.empty()) {
                TH2F* hfr=new TH2F(Form("hfr_empty_%s_p%d",dk.c_str(),panel.padIdx),
                    Form("%s;Distance [cm];%s",padTitle,yAxisTitle),1,0,700,1,-20,20);
                hfr->SetStats(0);hfr->Draw();
                TLatex tl;tl.SetNDC();tl.SetTextSize(0.07);tl.DrawLatex(0.3,0.5,"No data");
                continue;
            }
            int np=(int)pts.size();
            vector<double> r_v(np),theo_v(np),meas_v(np),meas_e(np),zero_v(np,0.0);
            for(int k=0;k<np;++k){
                r_v[k]=pts[k].r;theo_v[k]=pts[k].theo_rel;
                meas_v[k]=pts[k].meas_mean;meas_e[k]=pts[k].meas_err;
            }

            TH2F* hfr=new TH2F(Form("hfr_%s_%s_p%d",filePrefix,dk.c_str(),panel.padIdx),
                Form("%s;Distance from diffuser [cm];%s",padTitle,yAxisTitle),
                1,rMin,rMax,1,yMin,yMax);
            hfr->SetStats(0);
            hfr->GetXaxis()->SetTitleSize(0.050);hfr->GetYaxis()->SetTitleSize(0.050);
            hfr->GetXaxis()->SetLabelSize(0.042);hfr->GetYaxis()->SetLabelSize(0.042);
            hfr->GetYaxis()->SetTitleOffset(1.20);hfr->Draw();

            // Expectation: plot A -> r/v - <r/v> ; plot B (residual) -> 0
            vector<pair<double,double>> tsorted;
            for(int k=0;k<np;++k) tsorted.push_back({r_v[k], residualPlot ? 0.0 : theo_v[k]});
            sort(tsorted.begin(),tsorted.end());
            vector<double> rt(np),tt(np);
            for(int k=0;k<np;++k){rt[k]=tsorted[k].first;tt[k]=tsorted[k].second;}
            TGraph* grTheo=new TGraph(np,rt.data(),tt.data());
            grTheo->SetName(Form("gr_theo_%s_%s_p%d",filePrefix,dk.c_str(),panel.padIdx));
            grTheo->SetLineColor(kBlue+1);grTheo->SetLineWidth(2);grTheo->SetMarkerStyle(1);
            grTheo->Draw("L same");

            // Weighted straight-line fit, every channel used
            LineFit lf = WeightedLineFit(r_v, meas_v, meas_e);
            TGraphErrors* grMeas=new TGraphErrors(np,r_v.data(),meas_v.data(),zero_v.data(),meas_e.data());
            grMeas->SetName(Form("gr_meas_%s_%s_p%d",filePrefix,dk.c_str(),panel.padIdx));
            grMeas->SetMarkerStyle(20);grMeas->SetMarkerColor(kBlack);
            grMeas->SetMarkerSize(0.9);grMeas->SetLineColor(kBlack);grMeas->SetLineWidth(1);
            grMeas->Draw("P same");
            int nOut = 0;
            for (int k=0;k<np;++k) if (meas_v[k]<yMin || meas_v[k]>yMax) ++nOut;

            TF1* fMeas=new TF1(Form("fMeas_%s_%s_p%d",filePrefix,dk.c_str(),panel.padIdx),
                               "[0]+[1]*x",rMin,rMax);
            fMeas->SetParameters(lf.p0, lf.p1);
            fMeas->SetParError(0, lf.e0); fMeas->SetParError(1, lf.e1);
            fMeas->SetLineColor(kRed);fMeas->SetLineStyle(2);fMeas->SetLineWidth(2);
            if (lf.ok) fMeas->Draw("same");
            fitOut = {lf.p0, lf.p1};
            // slope error inflated by sqrt(chi2/ndf) when chi2/ndf > 1
            double sclE = (lf.ndf>0 && lf.chi2>lf.ndf) ? sqrt(lf.chi2/lf.ndf) : 1.0;
            cout << Form("  [fit] %-16s %-8s slope = %8.4f +- %.4f ns/cm  chi2/ndf = %.1f/%d  (%d points outside y-range)",
                         filePrefix, dk.c_str(), lf.p1, lf.e1*sclE, lf.chi2, lf.ndf, nOut) << endl;

            TLegend* leg=new TLegend(0.15,0.66,0.80,0.89);
            leg->SetBorderSize(0);leg->SetFillStyle(0);
            leg->SetTextSize(0.034);
            leg->AddEntry(grTheo, residualPlot
                ? Form("Expected: 0 (r/v removed, 1/%g = %.4f ns/cm)", kLightSpeed_cm_ns, 1.0/kLightSpeed_cm_ns)
                : Form("Theory: 1/%g = %.4f ns/cm", kLightSpeed_cm_ns, 1.0/kLightSpeed_cm_ns),"l");
            leg->AddEntry(grMeas, kTimingMode==0 ? "Per-channel data"
                                                 : Form("Per-channel data (averaged waveform, CF %.0f%%)", 100*kCFfrac),"p");
            leg->AddEntry(fMeas, lf.ok ? Form("Fit: %.4f #pm %.4f ns/cm", lf.p1, lf.e1*sclE)
                                       : "Fit: fewer than 3 points", "l");
            if (lf.ok) leg->AddEntry((TObject*)nullptr, Form("(#chi^{2}/ndf = %.1f/%d)", lf.chi2, lf.ndf), "");
            leg->Draw();
            if (nOut) {
                TLatex to; to.SetNDC(); to.SetTextSize(0.030); to.SetTextColor(kGray+2);
                to.DrawLatex(0.16,0.16,Form("%d point(s) outside the y-range (still in the fit)", nOut));
            }

            grTheo->Write();grMeas->Write();fMeas->Write();
        }
        cRad->SaveAs(Form("%s/%s_%s.png",plotDir.c_str(),filePrefix,dk.c_str()));
        cRad->Write(Form("%s_%s",canvPrefix,dk.c_str()));
        delete cRad;
    };

    // ============================================================
    // DRAW: outlier plot for one APA half
    // Shows: (top) deviation from fit A per channel
    //        (bottom) deviation from fit B per channel
    // Channels beyond ±kOutlierNsigma·σ are annotated by ID.
    // ============================================================
    auto drawOutlierHalf = [&](
        const vector<PMTPoint>& ptsA,
        const vector<PMTPoint>& ptsB,
        const FitResult& fitA,
        const FitResult& fitB,
        const string& dk,
        const char* apaSide)
    {
        if (ptsA.empty()&&ptsB.empty()) return;

        // compute deviations
        auto makeDevs = [](const vector<PMTPoint>& pts,
                           const FitResult& fit,
                           vector<double>& ch_v,
                           vector<double>& dv,
                           vector<double>& ev)
        {
            for (const auto& p : pts) {
                ch_v.push_back((double)p.ch);
                dv.push_back(p.meas_mean - (fit.p0 + fit.p1 * p.r));
                ev.push_back(p.meas_err);
            }
        };

        vector<double> chA,dA,eA, chB,dB,eB;
        makeDevs(ptsA,fitA,chA,dA,eA);
        makeDevs(ptsB,fitB,chB,dB,eB);

        auto [muA,sigA] = VecStats(dA);
        auto [muB,sigB] = VecStats(dB);

        // console summary
        cout << "  [outlier] " << dk << " " << apaSide
             << "  sigA=" << Form("%.2f",sigA) << " ns"
             << "  sigB=" << Form("%.2f",sigB) << " ns"
             << "  | flagged (plot A, >" << kOutlierNsigma << "sig): ";
        int nflag=0;
        for (size_t k=0;k<dA.size();++k)
            if (fabs(dA[k]-muA)>kOutlierNsigma*sigA){cout<<(int)chA[k]<<" ";++nflag;}
        if (!nflag) cout<<"(none)";
        cout<<endl;

        // canvas: two pads stacked vertically
        string cname=Form("c_outlier_%s_%s",dk.c_str(),apaSide);
        TCanvas* c=new TCanvas(cname.c_str(),
            Form("Outlier  %s  APA=%s",dk.c_str(),apaSide), 1100, 820);
        c->Divide(1,2,0,0.005);

        auto drawOutPad = [&](
            int padIdx,
            const vector<double>& ch_v,
            const vector<double>& d_v,
            const vector<double>& e_v,
            double mu, double sig,
            const char* ytitle, const char* label)
        {
            c->cd(padIdx);
            gPad->SetLeftMargin(0.10);gPad->SetRightMargin(0.04);
            gPad->SetBottomMargin(0.14);gPad->SetTopMargin(0.10);

            int n=(int)ch_v.size();
            if (n==0) return;

            double chMin=*min_element(ch_v.begin(),ch_v.end())-8;
            double chMax=*max_element(ch_v.begin(),ch_v.end())+8;
            double dMin=*min_element(d_v.begin(),d_v.end());
            double dMax=*max_element(d_v.begin(),d_v.end());
            double pad=max(0.15*(dMax-dMin),3.0);
            dMin-=pad; dMax+=pad;

            TH2F* hfr=new TH2F(
                Form("hfr_out_%s_%s_p%d",dk.c_str(),apaSide,padIdx),
                Form("%s  APA=%s  |  %s;Channel ID;%s",
                     dk.c_str(),apaSide,label,ytitle),
                1,chMin,chMax,1,dMin,dMax);
            hfr->SetStats(0);
            hfr->GetXaxis()->SetTitleSize(0.048);hfr->GetYaxis()->SetTitleSize(0.048);
            hfr->GetXaxis()->SetLabelSize(0.040);hfr->GetYaxis()->SetLabelSize(0.040);
            hfr->Draw();

            vector<double> zero_v(n,0.0);
            TGraphErrors* gr=new TGraphErrors(n,ch_v.data(),d_v.data(),
                                              zero_v.data(),e_v.data());
            gr->SetName(Form("gr_out_%s_%s_p%d",dk.c_str(),apaSide,padIdx));
            gr->SetMarkerStyle(20);gr->SetMarkerColor(kBlack);gr->SetMarkerSize(0.9);
            gr->Draw("P same");
            gr->Write();

            // Zero line
            TLine* l0=new TLine(chMin,0,chMax,0);
            l0->SetLineColor(kGray+2);l0->SetLineStyle(1);l0->SetLineWidth(1);l0->Draw();

            // ±nsigma band
            double lo=mu-kOutlierNsigma*sig, hi=mu+kOutlierNsigma*sig;
            for (double yy : {lo, hi}) {
                TLine* ll=new TLine(chMin,yy,chMax,yy);
                ll->SetLineColor(kOrange+7);ll->SetLineStyle(2);ll->SetLineWidth(2);
                ll->Draw();
            }

            // Annotate outlier channels
            TLatex lat; lat.SetTextSize(0.040); lat.SetTextAlign(22);
            lat.SetTextColor(kRed+1);
            for (int k=0;k<n;++k) {
                if (fabs(d_v[k]-mu)>kOutlierNsigma*sig) {
                    double ly=(d_v[k]>mu) ? d_v[k]+0.05*(dMax-dMin)
                                          : d_v[k]-0.05*(dMax-dMin);
                    lat.DrawLatex(ch_v[k], ly, Form("%d",(int)ch_v[k]));
                }
            }

            // Stats annotation
            TLatex stats; stats.SetNDC(); stats.SetTextSize(0.040);
            stats.DrawLatex(0.12,0.85,
                Form("#mu=%.2f ns   #sigma=%.2f ns   "
                     "#pm%.0f#sigma threshold shown",mu,sig,kOutlierNsigma));
        };

        drawOutPad(1, chA, dA, eA, muA, sigA,
            "Δt_{exp} deviation from fit [ns]",
            Form("Outliers of #Delta t_{exp} plot (slope = %.4f ns/cm)", fitA.p1));

        drawOutPad(2, chB, dB, eB, muB, sigB,
            "Offset from r/v expectation [ns]",
            "Channel offsets = plot-B values (expectation 0)");
        c->SaveAs(Form("%s/outlier_%s_%s.png",plotDir.c_str(),dk.c_str(),apaSide));
        c->Write(cname.c_str());
        delete c;

        // Also draw a 1-D histogram of the deviations for this half
        {
            TCanvas* ch=new TCanvas(Form("c_devhist_%s_%s",dk.c_str(),apaSide),
                Form("Deviation hist  %s %s",dk.c_str(),apaSide),900,400);
            ch->Divide(2,1);
            auto makeH=[&](int pad, const vector<double>& dv,
                           double mu, double sig, const char* title){
                ch->cd(pad);
                double rng=max(3.0*sig,15.0);
                TH1F* h=new TH1F(
                    Form("h_devhist_%s_%s_p%d",dk.c_str(),apaSide,pad),
                    Form("%s  APA=%s;Deviation from fit A / offset (B) [ns];Channels",
                         title,apaSide),
                    40,-rng,rng);
                for (double v:dv) h->Fill(v);
                h->SetFillColor(pad==1?kAzure-4:kOrange-3);
                h->Draw();
                // Gaussian fit overlay
                if (h->GetEntries()>=4) {
                    TF1* fg=new TF1(Form("fg_dh_%d",pad),"gaus",-rng,rng);
                    h->Fit(fg,"RQ0");
                    fg->SetLineColor(kRed);fg->Draw("same");
                }
                TLine* lmu=new TLine(mu,0,mu,h->GetMaximum());
                lmu->SetLineColor(kGreen+2);lmu->SetLineWidth(2);lmu->Draw();
                h->Write();
            };
            makeH(1,dA,muA,sigA,Form("Plot A  %s",dk.c_str()));
            makeH(2,dB,muB,sigB,Form("Offsets (plot B)  %s",dk.c_str()));
            if (kSaveExtraPlots) ch->SaveAs(Form("%s/devhist_%s_%s.png",plotDir.c_str(),dk.c_str(),apaSide));
            ch->Write(Form("c_devhist_%s_%s",dk.c_str(),apaSide));
            delete ch;
        }
    };

    // ============================================================
    // MAIN LOOP OVER DIFFUSERS
    // ============================================================

    // For summary scatter (ALL mode)
    struct SumPt { string dk; string apa; int ch; double dA; double dB; double eB; };
    const FitResult kZeroFit = {0,0};   // plot B: expectation is 0
    vector<SumPt> summaryPts;

    for (int i=1;i<=5;++i) for (int s=0;s<2;++s) {
        const char* side = s?"East":"West";
        string dk = Form("%s_%d",side,i);
        if (!filterAll && dk!=diffuser_filter) continue;

        // ----------------------------------------------------------
        // Only fill matching APA half:
        //   West diffusers (s==0) → posX (West APA)
        //   East diffusers (s==1) → negX (East APA)
        // Cross-TPC side is intentionally left empty.
        // ----------------------------------------------------------
        bool isDiffuserWest = (s == 0);

        vector<PMTPoint> pA_pos,pA_neg,pB_pos,pB_neg;
        if (isDiffuserWest) {
            fillHalf(dk, true,  true,  pA_pos);   // posX, Plot A
            fillHalf(dk, true,  false, pB_pos);   // posX, Plot B
            // pA_neg, pB_neg left empty (cross-TPC, no LED light)
        } else {
            fillHalf(dk, false, true,  pA_neg);   // negX, Plot A
            fillHalf(dk, false, false, pB_neg);   // negX, Plot B
            // pA_pos, pB_pos left empty (cross-TPC, no LED light)
        }

        // ---- Plot A (dtexp) ---- fixed y-range [kRadiusYmin, kRadiusYmax]
        FitResult fA_pos={0,0}, fA_neg={0,0};
        drawRadiusCanvas(pA_neg,pA_pos,dk,i,side,
            "radius_vs_dtexp","c_radius_vs_dtexp",
            "T_{ch} - #bar{T} [ns]",
            true, false,           // fixedYRange, residualPlot
            fA_pos, fA_neg);

        // ---- Plot B (delay) ---- same fixed y-range
        FitResult fB_pos={0,0}, fB_neg={0,0};
        drawRadiusCanvas(pB_neg,pB_pos,dk,i,side,
            "radius_vs_delay","c_radius_vs_delay",
            "#Deltat_{exp,rel} - #Deltat_{theo,rel} [ns]",
            true, true,            // fixedYRange, residualPlot
            fB_pos, fB_neg);

        // ---- Outlier plots (only matching APA half has data) ----
        drawOutlierHalf(pA_pos,pB_pos,fA_pos,kZeroFit,dk,"posX");
        drawOutlierHalf(pA_neg,pB_neg,fA_neg,kZeroFit,dk,"negX");

        // ---- Collect for summary ----
        auto collect=[&](const vector<PMTPoint>& pA,
                         const vector<PMTPoint>& pB,
                         const FitResult& fitA,
                         const FitResult& fitB,
                         const char* apaLabel)
        {
            // fitB is kZeroFit -> dB = channel offset from the r/v expectation
            map<int,pair<double,double>> dB_map;
            for (const auto& p:pB) dB_map[p.ch]={p.meas_mean-(fitB.p0+fitB.p1*p.r), p.meas_err};
            for (const auto& p:pA) {
                if (!dB_map.count(p.ch)) continue;
                double dA_val=p.meas_mean-(fitA.p0+fitA.p1*p.r);
                summaryPts.push_back({dk,apaLabel,p.ch,dA_val,dB_map.at(p.ch).first,dB_map.at(p.ch).second});
            }
        };
        collect(pA_pos,pB_pos,fA_pos,kZeroFit,"posX");
        collect(pA_neg,pB_neg,fA_neg,kZeroFit,"negX");
    }

    // ============================================================
    // 2-D: Diffuser (x) vs Channel ID (y)
    //      z = measured dt_exp minus fit value at that channel's radius
    //          (deviation from the linear fit in radius_vs_dtexp, Plot A)
    // ============================================================
    {
        cout << "=== Building 2D diffuser-vs-channel deviation map ===" << endl;

        // ---- Build ordered diffuser list (only those with data) ----
        vector<string> dkOrder;
        {
            set<string> seen;
            for (int s = 0; s < 2; ++s)
                for (int i = 1; i <= 5; ++i) {
                    string dk = Form("%s_%d", s ? "East" : "West", i);
                    if (seen.count(dk)) continue;
                    bool hasData = false;
                    for (const auto& sp : summaryPts)
                        if (sp.dk == dk) { hasData = true; break; }
                    if (hasData) { dkOrder.push_back(dk); seen.insert(dk); }
                }
        }

        // ---- Collect unique channel IDs, sorted ----
        set<int> chSet;
        for (const auto& sp : summaryPts) chSet.insert(sp.ch);
        vector<int> chVec(chSet.begin(), chSet.end());
        sort(chVec.begin(), chVec.end());

        int nDK = (int)dkOrder.size();
        int nCH = (int)chVec.size();
        if (nDK == 0 || nCH == 0) {
            cout << "  No data for 2D map – skipping." << endl;
            goto end_2d_map;
        }
        {
            map<string, int> dkIdx;
            for (int i = 0; i < nDK; ++i) dkIdx[dkOrder[i]] = i;
            map<int, int> chIdx;
            for (int i = 0; i < nCH; ++i) chIdx[chVec[i]] = i;

            // ---- Create histograms (one for Plot A deviation, one for Plot B) ----
            auto makeH2 = [&](const char* name, const char* title) -> TH2F* {
                TH2F* h = new TH2F(name, title, nDK, 0, nDK, nCH, 0, nCH);
                h->SetStats(0);
                for (int i = 0; i < nDK; ++i)
                    h->GetXaxis()->SetBinLabel(i + 1, dkOrder[i].c_str());
                for (int i = 0; i < nCH; ++i)
                    h->GetYaxis()->SetBinLabel(i + 1, Form("%d", chVec[i]));
                h->GetXaxis()->SetLabelSize(0.042);
                h->GetYaxis()->SetLabelSize(max(0.012, min(0.022, 18.0 / nCH)));
                h->GetXaxis()->LabelsOption("v");   // rotate x labels 90 deg
                h->GetXaxis()->SetTitle("");
                h->GetYaxis()->SetTitle("Channel ID");
                h->GetZaxis()->SetTitleOffset(1.4);
                return h;
            };
            TH2F* hDevA = makeH2("h2_devA_diff_ch",
                "dt_{exp} deviation from fit;"
                "Diffuser;Channel ID;Deviation [ns]");
            TH2F* hDevB = makeH2("h2_devB_diff_ch",
                "Channel offset from r/v expectation;Diffuser;Channel ID;Offset [ns]");

            // ---- Accumulate deviations (average if channel appears in both
            //      posX and negX for the same diffuser, which shouldn't happen,
            //      but handles it gracefully) ----
            map<pair<int,int>, pair<double,int>> accA, accB;
            for (const auto& sp : summaryPts) {
                if (!dkIdx.count(sp.dk) || !chIdx.count(sp.ch)) continue;
                auto key = make_pair(dkIdx.at(sp.dk), chIdx.at(sp.ch));
                accA[key].first += sp.dA;  accA[key].second++;
                accB[key].first += sp.dB;  accB[key].second++;
            }
            for (const auto& kv : accA)
                if (kv.second.second > 0)
                    hDevA->SetBinContent(kv.first.first + 1, kv.first.second + 1,
                                        kv.second.first / kv.second.second);
            for (const auto& kv : accB)
                if (kv.second.second > 0)
                    hDevB->SetBinContent(kv.first.first + 1, kv.first.second + 1,
                                        kv.second.first / kv.second.second);

            // ---- Symmetric z-range (same scale for both maps) ----
            auto getZmax = [](TH2F* h) -> double {
                double zm = 0;
                for (int bx = 1; bx <= h->GetNbinsX(); ++bx)
                    for (int by = 1; by <= h->GetNbinsY(); ++by)
                        zm = max(zm, fabs(h->GetBinContent(bx, by)));
                return max(zm, 5.0);
            };
            double zmA = getZmax(hDevA);
            double zmB = getZmax(hDevB);
            hDevA->GetZaxis()->SetRangeUser(-zmA, zmA);
            hDevB->GetZaxis()->SetRangeUser(-zmB, zmB);
            hDevA->GetZaxis()->SetTitle("dt_{exp} deviation from fit [ns]");
            hDevB->GetZaxis()->SetTitle("Offset from r/v expectation [ns]");

            // ---- Blue-white-red diverging palette ----
            // (saves and restores the current palette afterwards)
            const int kNStops = 5;
            double stops[kNStops] = {0.00, 0.25, 0.50, 0.75, 1.00};
            double reds[kNStops]  = {0.00, 0.40, 1.00, 1.00, 0.60};
            double gres[kNStops]  = {0.20, 0.60, 1.00, 0.40, 0.00};
            double bles[kNStops]  = {0.60, 1.00, 1.00, 0.40, 0.00};
            int fgi = TColor::CreateGradientColorTable(kNStops, stops, reds, gres, bles, 255);
            int devPalette[255];
            for (int ci = 0; ci < 255; ++ci) devPalette[ci] = fgi + ci;
            gStyle->SetPalette(255, devPalette);   // apply the just-created table
            gStyle->SetNumberContours(255);

            // ---- Canvas sizing and label geometry ----
            const int    kPxPerBin  = 28;
            int canvH = max(600, nCH * kPxPerBin + 220);
            int canvW = max(600, nDK * 90 + 280);

            auto drawAndSave = [&](TH2F* h, const char* fname, const char* cname) {
                TCanvas* c = new TCanvas(cname, cname, canvW, canvH);
                c->SetLeftMargin(0.09);
                c->SetRightMargin(0.18);
                double botMargin = max(0.10, min(0.22, (double)(nDK * 14 + 30) / canvH));
                c->SetBottomMargin(botMargin);
                c->SetTopMargin(0.06);

                h->GetXaxis()->SetLabelFont(43);
                h->GetXaxis()->SetLabelSize(20);
                h->GetYaxis()->SetLabelFont(43);
                h->GetYaxis()->SetLabelSize(20);
                h->GetYaxis()->SetTickLength(0.005);
                h->GetXaxis()->LabelsOption("v");
                h->Draw("COLZ");

                gPad->SetGridy(1);
                gPad->SetGridx(1);

                for (int bx = 1; bx < nDK; ++bx) {
                    TLine* vl = new TLine(bx, 0, bx, nCH);
                    vl->SetLineColor(kGray + 1);
                    vl->SetLineWidth(1);
                    vl->SetLineStyle(3);
                    vl->Draw();
                }
                c->SaveAs(Form("%s/%s.png", plotDir.c_str(), fname));
                c->Write(cname);
                h->Write();
                delete c;
            };

            if (kSaveExtraPlots) drawAndSave(hDevA, "devA_diffuser_vs_channel", "c_devA_diffuser_vs_channel");
            else hDevA->Write();
            drawAndSave(hDevB, "offset_diffuser_vs_channel", "c_offset_diffuser_vs_channel");

            // Restore ROOT's default rainbow palette
            gStyle->SetPalette(kBird);
            gStyle->SetNumberContours(99);
        }
        end_2d_map:;
    }

    // ============================================================
    // SUMMARY (ALL mode): scatter of deviation vs channel,
    //                     one color per diffuser
    // ============================================================
    if (filterAll && !summaryPts.empty()) {
        cout << "=== Summary scatter plots (ALL mode) ===" << endl;

        int colors[]={kRed,kBlue+1,kGreen+2,kOrange+7,kMagenta+1,
                      kCyan+1,kYellow+2,kViolet+2,kTeal+3,kPink+6};
        int nCol=sizeof(colors)/sizeof(int);

        map<string,vector<pair<double,double>>> scA,scB;
        for (const auto& sp:summaryPts){
            scA[sp.dk].push_back({(double)sp.ch,sp.dA});
            scB[sp.dk].push_back({(double)sp.ch,sp.dB});
        }

        auto drawSummary=[&](const map<string,vector<pair<double,double>>>& sc,
                             const char* fname, const char* yt)
        {
            TCanvas* c=new TCanvas(fname,fname,1400,600);
            c->SetLeftMargin(0.09);c->SetRightMargin(0.04);
            c->SetBottomMargin(0.12);c->SetTopMargin(0.09);
            TH2F* hfr=new TH2F(Form("hfr_%s",fname),
                Form("%s;Channel ID;%s",fname,yt),1,-10,320,1,-60,60);
            hfr->SetStats(0);hfr->Draw();
            TLine* l0=new TLine(-10,0,320,0);
            l0->SetLineColor(kGray+2);l0->SetLineWidth(1);l0->Draw();
            TLegend* leg=new TLegend(0.10,0.77,0.60,0.91);
            leg->SetNColumns(5);leg->SetTextSize(0.028);
            leg->SetBorderSize(0);leg->SetFillStyle(0);
            int ic=0;
            vector<double> allX, allY;
            for (const auto& kv:sc){
                int n=(int)kv.second.size(); if(!n){++ic;continue;}
                vector<double> xv(n),yv(n);
                for(int k=0;k<n;++k){xv[k]=kv.second[k].first;yv[k]=kv.second[k].second;}
                allX.insert(allX.end(),xv.begin(),xv.end()); allY.insert(allY.end(),yv.begin(),yv.end());
                TGraph* gr=new TGraph(n,xv.data(),yv.data());
                gr->SetName(Form("gr_%s_%s",fname,kv.first.c_str()));
                gr->SetMarkerStyle(20+(ic%10));
                gr->SetMarkerColor(colors[ic%nCol]);
                gr->SetMarkerSize(0.9);
                gr->Draw("P same");
                leg->AddEntry(gr,kv.first.c_str(),"p");
                gr->Write();
                ++ic;
            }
            LabelOutlierChannels(allX, allY, hfr->GetYaxis()->GetXmin(), hfr->GetYaxis()->GetXmax(), 0.028);
            leg->Draw();
            c->SaveAs(Form("%s/%s.png",plotDir.c_str(),fname));
            c->Write(fname);
            delete c;
        };
        if (kSaveExtraPlots) drawSummary(scA,"summary_scatter_devFromFitA_vs_ch",
                    "dtexp deviation from fit A [ns]");
        drawSummary(scB,"summary_scatter_offset_vs_ch",
                    "offset from r/v expectation (plot B) [ns]");

        // 1-D spread histograms (all diffusers combined)
        TH1F* hspA=new TH1F("h_spread_A_ALL",
            "All diffusers: dtexp dev from fit A;#Delta [ns];Entries",80,-40,40);
        TH1F* hspB=new TH1F("h_spread_B_ALL",
            "All diffusers: offset from r/v expectation (plot B);#Delta [ns];Entries",80,-40,40);
        for (const auto& sp:summaryPts){hspA->Fill(sp.dA);hspB->Fill(sp.dB);}
        TCanvas* csp=new TCanvas("c_spread_ALL","Spread all diffs",1200,500);
        csp->Divide(2,1);
        csp->cd(1);hspA->SetFillColor(kAzure-4);hspA->Draw();
        csp->cd(2);hspB->SetFillColor(kOrange-3);hspB->Draw();
        if (kSaveExtraPlots) csp->SaveAs(Form("%s/summary_spread_histos.png",plotDir.c_str()));
        csp->Write("c_spread_ALL");
        hspA->Write();hspB->Write();
        delete csp;
    }

    // ============================================================
    // CHANNEL TIMING OFFSETS  and  TYPE COMPARISON (TPB-coated vs uncoated PMTs)
    //   offset(ch) = mean over diffusers of the plot-B value
    //                #Deltat_exp,rel - #Deltat_theo,rel  (size corrected,
    //                relative to the average of that diffuser's channels)
    //   error      = spread over diffusers / sqrt(N)   (N >= 2)
    //                statistical error of the single value (N = 1)
    //   Table written to timing_offsets_<tag>.txt  (calibration output)
    // ============================================================
    if (!summaryPts.empty()) {
        cout << "=== Channel timing offsets (TPB-coated vs uncoated PMTs) ===" << endl;
        map<int, vector<pair<double,double>>> offs;   // ch -> {offset, stat err} per diffuser
        for (const auto& sp : summaryPts) offs[sp.ch].push_back({sp.dB, sp.eB});
        struct Off { int ch; double m, e; int n; };
        vector<Off> tab;
        for (const auto& kv : offs) {
            int n = (int)kv.second.size();
            double s = 0; for (const auto& p : kv.second) s += p.first;
            double m = s/n, e;
            if (n >= 2) {
                double v = 0; for (const auto& p : kv.second) v += (p.first-m)*(p.first-m);
                e = sqrt(v/(n-1))/sqrt((double)n);
            } else e = kv.second[0].second;
            tab.push_back({kv.first, m, e, n});
        }
        ofstream ftab(Form("timing_offsets_%s.txt", outTag.c_str()));
        ftab << "# channel  type  offset_ns  err_ns  n_diffusers\n"
             << "# offset = mean over diffusers of (dt_exp,rel - dt_theo,rel), size corrected (plot B values)\n";
        for (const auto& o : tab)
            ftab << Form("%5d  %-12s %8.2f  %6.2f  %d\n", o.ch, IsTypeB(o.ch) ? kTypeNameB : kTypeNameA, o.m, o.e, o.n);
        ftab.close();
        cout << "  table written to timing_offsets_" << outTag << ".txt" << endl;

        vector<double> xA,yA,eA,xB,yB,eB;
        double chLo=1e9,chHi=-1e9,yLo=-40,yHi=40;
        for (const auto& o : tab) {
            bool tb = IsTypeB(o.ch);
            (tb?xB:xA).push_back(o.ch); (tb?yB:yA).push_back(o.m); (tb?eB:eA).push_back(o.e);
            chLo=min(chLo,(double)o.ch); chHi=max(chHi,(double)o.ch);
            yLo=min(yLo,o.m-o.e-10); yHi=max(yHi,o.m+o.e+10);
        }
        double chPad = max(2.0, 0.03*(chHi-chLo));
        auto meanSE = [](const vector<double>& v, double& m, double& se) {
            m = se = 0; if (v.empty()) return;
            for (double x : v) m += x;
            m /= v.size();
            if (v.size() < 2) return;
            double s = 0; for (double x : v) s += (x-m)*(x-m);
            se = sqrt(s/(v.size()-1))/sqrt((double)v.size());
        };
        double mA,sA,mB,sB; meanSE(yA,mA,sA); meanSE(yB,mB,sB);
        cout << Form("  %s: %d ch, mean offset %+.2f +- %.2f ns", kTypeNameA, (int)yA.size(), mA, sA) << endl;
        cout << Form("  %s: %d ch, mean offset %+.2f +- %.2f ns", kTypeNameB, (int)yB.size(), mB, sB) << endl;

        TCanvas* c = new TCanvas("c_type_offsets","",1500,560);
        TPad* p1 = new TPad("p_type_1","",0.00,0,0.70,1); p1->Draw();
        TPad* p2 = new TPad("p_type_2","",0.70,0,1.00,1); p2->Draw();
        p1->cd(); p1->SetLeftMargin(0.09); p1->SetRightMargin(0.02); p1->SetBottomMargin(0.12);
        TH2F* fr = new TH2F("hfr_type_offsets",
            "Channel timing offset (mean over diffusers);Channel ID;Offset: #Deltat_{exp,rel} - #Deltat_{theo,rel} [ns]",
            1,chLo-chPad,chHi+chPad, 1,yLo,yHi);
        fr->SetStats(0); fr->GetYaxis()->SetTitleOffset(1.1); fr->Draw();
        TLine* l0 = new TLine(chLo-chPad,0,chHi+chPad,0); l0->SetLineColor(kGray+1); l0->Draw();
        TLegend* lg = new TLegend(0.11,0.80,0.80,0.89);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035); lg->SetNColumns(2);
        TH1F* hA = new TH1F("h_type_offsets_A",";Offset [ns];Channels",30,-60,60);
        TH1F* hB = new TH1F("h_type_offsets_B",";Offset [ns];Channels",30,-60,60);
        for (int t = 0; t < 2; ++t) {
            const vector<double>& X = t ? xB : xA;
            if (X.empty()) continue;
            const vector<double>& Y = t ? yB : yA;
            const vector<double>& E = t ? eB : eA;
            vector<double> z(X.size(), 0.0);
            TGraphErrors* g = new TGraphErrors((int)X.size(), X.data(), Y.data(), z.data(), E.data());
            g->SetName(t ? "gr_type_offsets_B" : "gr_type_offsets_A");
            int col = t ? kRed+1 : kBlack;
            g->SetMarkerStyle(t ? 25 : 20); g->SetMarkerSize(t ? 1.1 : 0.9);
            g->SetMarkerColor(col); g->SetLineColor(col);
            g->Draw("P same"); g->Write();
            lg->AddEntry(g, Form("%s (%d ch): %+.1f #pm %.1f ns", t ? kTypeNameB : kTypeNameA,
                                 (int)X.size(), t ? mB : mA, t ? sB : sA), "p");
            TH1F* h = t ? hB : hA;
            for (double y : Y) h->Fill(y);
            h->SetLineColor(col); h->SetLineWidth(2);
        }
        {
            vector<double> allX(xA), allY(yA);
            allX.insert(allX.end(), xB.begin(), xB.end()); allY.insert(allY.end(), yB.begin(), yB.end());
            LabelOutlierChannels(allX, allY, yLo, yHi, 0.032);
        }
        lg->Draw();
        p2->cd(); p2->SetLeftMargin(0.17); p2->SetRightMargin(0.04); p2->SetBottomMargin(0.12);
        hA->SetStats(0); hB->SetStats(0);
        hA->SetMaximum(1.25*max(hA->GetMaximum(), hB->GetMaximum()));
        hA->SetTitle("All channels (under/overflow not drawn)");
        hA->Draw("HIST"); hB->Draw("HIST same");
        c->SaveAs(Form("%s/%s_offsets.png", plotDir.c_str(), kTypeFileTag));
        c->Write(); hA->Write(); hB->Write();
        delete c;
    }

    // ============================================================
    // PMT + DIFFUSER POSITION MAP
    // ============================================================
    {
        for (int s=0;s<2;++s) {
            const char* side=s?"East":"West";
            if (!filterAll) {
                string fside=diffuser_filter.substr(0,diffuser_filter.find('_'));
                if (fside!=side) continue;
            }
            const map<int,Point3D>& dpos=s?east_diffuser_pos:west_diffuser_pos;
            vector<double> py,pz; vector<int> pchs;
            for (const auto& kv:chanID_to_pos)
                if ((kv.second.x<0)==(!s)){
                    py.push_back(kv.second.y);pz.push_back(kv.second.z);
                    pchs.push_back(kv.first);
                }
            TCanvas* c=new TCanvas(Form("c_map_%s",side),
                Form("PMT+diffuser map - %s APA",side),1060,850);
            c->SetLeftMargin(0.11);c->SetBottomMargin(0.11);
            c->SetRightMargin(0.04);c->SetTopMargin(0.09);
            TH2F* hf=new TH2F(Form("hframe_map_%s",side),
                Form("%s APA: PMT+diffuser positions;Z [cm];Y [cm]",side),
                1,0,510,1,-230,230);
            hf->SetStats(0);
            hf->GetXaxis()->SetTitleSize(0.048);hf->GetYaxis()->SetTitleSize(0.048);
            hf->GetXaxis()->SetLabelSize(0.040);hf->GetYaxis()->SetLabelSize(0.040);
            hf->GetYaxis()->SetTitleOffset(1.1);hf->Draw();
            TLatex latCh;latCh.SetTextSize(0.018);latCh.SetTextAlign(22);
            latCh.SetTextColor(kWhite);latCh.SetTextFont(62);
            for (size_t j=0;j<py.size();++j){
                TGraph* gp=new TGraph(1,&pz[j],&py[j]);
                gp->SetMarkerStyle(20);gp->SetMarkerColor(kAzure+1);gp->SetMarkerSize(4.5);
                gp->Draw("P same");
                latCh.DrawLatex(pz[j],py[j],Form("%d",pchs[j]));
            }
            TLatex latD;latD.SetTextSize(0.050);latD.SetTextAlign(12);latD.SetTextFont(62);
            for (int i=1;i<=5;++i){
                const Point3D& dp=dpos.at(i);
                TGraph* gd=new TGraph(1,&dp.z,&dp.y);
                gd->SetMarkerStyle(29);gd->SetMarkerColor(kRed+1);gd->SetMarkerSize(3.5);
                gd->Draw("P same");
                latD.SetTextColor(kRed+1);
                latD.DrawLatex(dp.z+10,dp.y,Form("%d",i));
            }
            TLegend* leg=new TLegend(0.68,0.80,0.96,0.91);
            leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(0.034);
            TGraph* gl1=new TGraph(1);gl1->SetMarkerStyle(20);
            gl1->SetMarkerColor(kAzure+1);gl1->SetMarkerSize(4.5);
            leg->AddEntry(gl1,"PMT (channel ID)","p");
            TGraph* gl2=new TGraph(1);gl2->SetMarkerStyle(29);
            gl2->SetMarkerColor(kRed+1);gl2->SetMarkerSize(3.5);
            leg->AddEntry(gl2,"Diffuser","p");leg->Draw();
            c->SaveAs(Form("%s/PMT_diffuser_map_%s.png",plotDir.c_str(),side));
            c->Write(Form("c_PMT_diffuser_map_%s",side));
            delete hf;delete leg;delete gl1;delete gl2;delete c;
        }
    }

    outfile->Close();
    cout << "\n=== Analysis complete ===" << endl;
    cout << "Output ROOT : " << rootName  << endl;
    cout << "Plots dir   : " << plotDir   << "/" << endl;
    cout << "Have a nice day!!" << endl;
}