# SBND PMT timing calibration with the LED diffusers — full walk-through

Written 2026-11-13. Matches `PMTTimingAnalyzer.C` as of the same date
(XA version: `../XA_timing/XATimingAnalyzer.C`, same logic with 16 ns ticks).

Runs 19403–19467: 5 West diffusers + 5 East diffusers, mounted on the cathode.
West diffusers light the +X PMTs, East diffusers light the −X PMTs (same-side only).

All plot paths below are relative to `/exp/sbnd/data/users/castalyf/PMT_timing/`.

---------------------------------------------------------------------------

## 0. The goal in one picture

The LED fires; light travels from the diffuser to every PMT; every PMT records
a pulse. The time each PMT records is

    t(ch, e) = T0(e)  +  r(ch,d)/v  +  c_ch  +  W(Q)  +  noise

| piece      | what it is                                                   |
|------------|--------------------------------------------------------------|
| T0(e)      | moment the LED fired on flash e (unknown, jitters)           |
| r(ch,d)/v  | light travel time from diffuser d to PMT ch (from geometry)  |
| **c_ch**   | **fixed delay of the channel (cable, PMT, electronics) = what we want** |
| W(Q)       | extra delay because the pulse is small (few photons)         |

We can only measure differences between channels (T0 is unknown), so the
answer is c_ch relative to the average channel: Σ c_ch = 0.

**Plots**
- `Plots_ALL/PMT_diffuser_map_West.png`, `Plots_ALL/PMT_diffuser_map_East.png`
  — where the PMTs and diffusers are (this is where r comes from).
- `Plots_ALL/size_map_West_1.png` … `size_map_East_5.png`
  — how bright each PMT was for one diffuser. The brightest region should sit
  next to the red star (= this diffuser). Checks the run → diffuser mapping.

---------------------------------------------------------------------------

## 1. Symbols (used consistently below)

| symbol        | meaning                                                         | unit   |
|---------------|-----------------------------------------------------------------|--------|
| ch            | channel (one PMT)                                               | –      |
| d             | diffuser, e.g. West_3                                           | –      |
| e             | flash number (one LED flash = one event)                        | –      |
| tick          | one digitizer sample. **PMT: 1 tick = 2 ns** (500 MHz). XA: 16 ns | –    |
| t(ch,e)       | time channel ch sees light from flash e                         | ns     |
| A             | peak height of the pulse                                        | ADC    |
| Q(ch,e)       | size (area) of the pulse                                        | ADC·ns |
| u(ch,e)       | ln Q — natural log of the size                                  | –      |
| r(ch,d)       | distance diffuser d → PMT ch                                    | cm     |
| v             | speed of light in liquid argon, we use 22 cm/ns                 | cm/ns  |
| T_ch          | channel's average time within one diffuser run                  | ns     |
| o_e           | how early/late the LED fired on flash e vs. its average         | ns     |
| W(u)          | extra delay because the pulse is small (per-diffuser version)   | ns     |
| u0            | reference size where W = 0 (median of u)                        | –      |
| a1,a2,a3      | coefficients of W                                               | ns     |
| R(ch,d)       | channel's average time for diffuser d, LED jitter removed       | ns     |
| ū(ch,d)       | channel's typical size for diffuser d (median of u)             | –      |
| c_ch          | channel's fixed delay = **calibration result**                  | ns     |
| D_d           | one constant per diffuser run (each run has its own clock zero) | ns     |
| W_g(u), g1..g3| size correction fitted on all diffusers of a side together      | ns     |
| t̂(ch,d)       | fully corrected channel time                                    | ns     |
| δ(ch,d)       | what's left after removing distance: the channel's delay        | ns     |

---------------------------------------------------------------------------

## 2. Step 1 — measure every pulse  (`FindFirstPulse`)

Each waveform is a list of ADC values, one every 2 ns (1 tick).
PMT pulses go **down**, so we flip them: s(k) = baseline − waveform(k),
so the pulse becomes positive.

1. **Baseline.** b = average of ticks 0–799. Also its RMS (noise level).
   If RMS > 8 ADC the flash is noisy → dropped for this channel.
2. **Find the pulse.** In ticks 1100–1500, find the **first** place where s
   stays above max(5 × RMS, 5 ADC) for 2 ticks in a row. This is the
   **first** pulse, not the biggest one.
3. **Peak.** A = highest s within 8 ticks after that point.
4. **Time = half-height on the rising edge.** Walk back from the peak to the
   tick k where s(k) < A/2 ≤ s(k+1). Draw a straight line between those two
   ticks:

       t_tick = k + (A/2 − s(k)) / (s(k+1) − s(k))
       t      = 2 ns × t_tick + cable delay of ch

   Example: s(k) = 10, s(k+1) = 30, A/2 = 20 → t_tick = k + 0.5.
   That is why we get times finer than 2 ns.
5. **Size.** Q = area under the pulse = 2 ns × Σ s(k), summed from 2 ticks
   before the pulse start to tick 1500.

Big Q = many photons = PMT close to / facing the diffuser.
Small Q = few photons = far away.

**Plots**
- `Plots/Plots-data_EventBuilder5_art1_run<RUN>_<N>_..._decoded_filtered/heatmap.png`
  — **the useful one.** All channels' raw waveforms of that file overlaid
  (x = tick, y = ADC, colour = how many samples). You can see:
  the baseline (~14250 ADC), the LED pulses going **down** starting near
  tick ~1200, and the analysis window 1100–1500. Good to show what a "pulse"
  and "baseline" are.
  Note: in the run 19403 example the dip continues to ~tick 1650 and
  overshoots upward after that; Q is only integrated up to tick 1500.
  There is also a narrow spike near tick ~900 in some channels (before our
  window, after the baseline region).
- `.../waveform.png`, `.../derivative.png`, `.../derivative_pulse.png` in the
  same folders — **not usable**: they show flat boxes (they were filled with
  the tick index instead of the ADC value). Don't show them.
- (No per-pulse plot is made by the current analyzer. If you want one for the
  technote, a single-channel waveform with the baseline, threshold, A/2 line
  and Q area marked would be the ideal Step-1 figure — ask and I'll add it.)

---------------------------------------------------------------------------

## 3. Why the pulse size matters, and why ln Q

**Time walk.** A dim pulse is made of a few photons arriving one by one, so
it rises slowly and reaches half-height **later**. A bright pulse rises fast.
So a small Q looks late even when the light arrived at the same moment.
Far PMTs are both far *and* dim — if we don't remove the size effect, they
look extra late for the wrong reason.

**Why ln Q.** Q spans ~100× between near and far PMTs. The delay depends on
*how many times* bigger/smaller a pulse is: 100 → 200 gives about the same
change as 1000 → 2000 (both "×2"). A logarithm turns "×2" into "+ a fixed
amount":

    ln(2Q) = ln Q + ln 2 = ln Q + 0.69

So we use **u = ln Q** as the size variable — data spread evenly, a smooth
curve fits. Plots show log10 Q instead (easier to read; log10 Q = 3 means
Q = 1000). Conversion: log10 Q = u / 2.303.

**Plots**
- `Plots_ALL/size_vs_distance_West_1.png` … `_East_5.png`, `size_vs_distance_ALL.png`
  — typical size of each channel vs its distance: size drops with distance,
  i.e. distance and size are mixed together. This is *why* a size correction
  is needed. Error bars = middle 68% of that channel's flashes.

---------------------------------------------------------------------------

## 4. Step 2 — per-diffuser fit  (`SolveDiffuser`, one diffuser at a time)

### Model

    t(ch,e) = T_ch + o_e + W(u(ch,e))

    W(u) = a1·(u − u0) + a2·(u − u0)² + a3·(u − u0)³

- T_ch: where this channel sits on average (contains r/v and c_ch, which
  cannot be separated within one diffuser).
- o_e: shift shared by all channels on flash e (LED fires a bit early/late).
- W: extra delay vs. a typical pulse. u0 = median u of this diffuser, so
  W(u0) = 0. u outside the 2%–98% range is clamped to the edge.

### Solving: fix two pieces, solve the third, repeat 5 times
Start with o_e = 0 and W = 0.

**(a) Channel average T_ch.** For each channel remove the current flash shift
and size effect:

    y(ch,e) = t(ch,e) − o_e − W(u(ch,e))

T_ch = **robust average** of that channel's y values:
1. m = median of y
2. MAD = median of |y − m|; spread σ ≈ 1.4826 × MAD
   (for a bell curve this equals the usual standard deviation)
3. drop any y with |y − m| > 3σ
4. T_ch = ordinary mean of what remains

**(b) Flash shift o_e.** For each flash, for every channel that saw it:

    z(ch,e) = t(ch,e) − W(u(ch,e)) − T_ch
    o_e     = median over channels of z(ch,e)

then subtract the average of all o_e so they average to 0.
Same job as Olivia's "event average", but the median means one bad channel
can't move it.

**(c) Size curve W.** Learn a1, a2, a3 only from changes *inside one channel*
(same PMT, slightly dimmer on some flashes → slightly later?).
Use pulses that passed the 3σ cut in (a). For each channel subtract its own
average from both the time and the size terms:

    ỹ   = (t − o_e) − [channel average of (t − o_e)]
    x̃_k = (u − u0)^k − [channel average of (u − u0)^k],   k = 1,2,3

and choose a1, a2, a3 to minimise (least squares)

    Σ over all pulses ( ỹ − a1·x̃1 − a2·x̃2 − a3·x̃3 )²

Removing each channel's own average stops W from absorbing
"far PMTs are late AND dim" (that's distance, not size).

### Output per channel, per diffuser
- R(ch,d) = robust average over flashes of [ t(ch,e) − o_e ]
  (LED jitter removed, size **not yet** corrected — Step 3 does it)
- ū(ch,d) = median of u(ch,e) over flashes

Quality cuts: channel needs ≥ 20 pulses and light in ≥ 20% of flashes;
a flash needs ≥ 5 channels.

**Plots:** none directly; Step 2 feeds Steps 3–7.

---------------------------------------------------------------------------

## 5. Step 3 — global size correction  (`FitGlobalSize`, all 5 diffusers of one side)

### Why Step 2's W isn't enough
Within one run a PMT's size changes only a little between flashes → W is
learned from a narrow range and comes out too flat. Across diffusers the
**same** PMT is very bright (diffuser next to it) or very dim (diffuser far
away). That wide range pins down the size effect.

### Model (one equation per channel–diffuser pair)

    R(ch,d) − r(ch,d)/v = c_ch + D_d + W_g(ū(ch,d))

    W_g(u) = g1·(u − u0) + g2·(u − u0)² + g3·(u − u0)³,
    u0 = median of all ū(ch,d)

Unknowns: ~150 c_ch, 5 D_d, g1..g3. Solve all at once by least squares:

    minimise  Σ over (ch,d) [ R − r/v − c_ch − D_d − W_g(ū) ]²
    with the rule Σ c_ch = 0

The rule is needed because adding +5 ns to every c and −5 ns to every D
fits equally well; only differences between channels are physical.

**Toy example.** PMT 42:
- West_1: ū = ln(3000) = 8.0, R − r/v = 10 ns
- West_3: ū = ln(150)  = 5.0, R − r/v = 16 ns

c_42 is the same in both rows. If the two D's are similar, the 6 ns must be
W_g(5.0) − W_g(8.0). Hundreds of channels doing this fix the shape of W_g.

### Corrected channel time

    t̂(ch,d) = R(ch,d) − W_g(ū(ch,d))

⚠ **Caveat:** r/v is subtracted *before* fitting W_g, so the speed of light
v is an **assumption** here, not a measurement. If v were also free it could
not be told apart from W_g.

**Plots**
- `Plots_ALL/size_vs_delay_West_1.png` … `_East_5.png`
  x = delay (R − r/v, minus its average), y = log10 Q.
  Left panel = before correction; magenta curve = W_g (points should follow it).
  Right panel = after subtracting W_g; points should line up at 0.
- `Plots_ALL/size_vs_delay_ALL.png`, `Plots_ALL/size_vs_delay_corrected_ALL.png`
  — same, all diffusers together.

---------------------------------------------------------------------------

## 6. Step 4 — does time grow with distance at 1/v?  (plot A)

Within one diffuser, compare to that diffuser's average (D_d cancels):

    Δt_exp,rel(ch)  = t̂(ch,d) − [average of t̂ over channels of d]
    Δt_theo,rel(ch) = r(ch,d)/v − [average of r/v over channels of d]

Fit a straight line Δt_exp,rel = p0 + p1·r, each point weighted by 1/σ²:

    minimise  Σ (Δt_exp,rel − p0 − p1·r)² / σ²

Expected slope p1 = 1/v = 1/22 = **0.0455 ns/cm**.
χ²/ndf ≈ 1 means a line fits well; if χ²/ndf > 1 the slope error is
multiplied by √(χ²/ndf).
Current PMT results: +0.010 to +0.062 ns/cm (±0.002).

**Plots**
- `Plots_ALL/radius_vs_dtexp_West_1.png` … `_East_5.png`
  — main validation plot (Olivia's "radius vs average time delay").

---------------------------------------------------------------------------

## 7. Step 5 — remove distance: what's left is the channel delay  (plot B)

    δ(ch,d) = Δt_exp,rel − Δt_theo,rel

Should be flat in r, scattered around 0; the scatter is the channels' fixed
delays c_ch.

**Plots**
- `Plots_ALL/radius_vs_delay_West_1.png` … `_East_5.png` — δ vs r (expect flat at 0).
- `Plots_ALL/outlier_West_1_posX.png` … `outlier_East_5_negX.png`
  — δ vs channel ID; channels more than 2σ from the mean are labelled red.

---------------------------------------------------------------------------

## 8. Step 6 — final calibration constants  (`timing_offsets_ALL.txt`)

A channel seen by N diffusers (N ≥ 2):

    c_ch  = ( δ(ch,d1) + δ(ch,d2) + … + δ(ch,dN) ) / N
    error = s / √N,   s = standard deviation of the N values of δ

**Plots**
- `Plots_ALL/offset_diffuser_vs_channel.png`, `Plots_ALL/summary_scatter_offset_vs_ch.png`
  — do different diffusers give the same δ for the same channel?
  **Strongest test**: each diffuser sees that channel at a different distance
  and brightness.
- `Plots_ALL/coated_vs_uncoated_offsets.png` — final c_ch per channel,
  TPB-coated vs uncoated PMTs.
- `Plots_ALL/res_mean_posX_All.png`, `Plots_ALL/res_mean_negX_All.png`
  — offsets drawn on the Z–Y map (look for patterns: one cable bundle / board off).

**Output file:** `timing_offsets_ALL.txt`

---------------------------------------------------------------------------

## 9. Step 7 — single-flash timing resolution

For each single pulse:

    res(ch,e) = [ t(ch,e) − o_e − size correction − average ] − [ r/v − average r/v ]

Histogram res per channel, fit a Gaussian (bell curve). Its width **σ** is
that channel's single-flash resolution. Dim (far) channels → fewer photons →
larger σ. So the σ maps are brightest **far** from the diffuser, not on it.

**Plots**
- `Plots_ALL/residual_sigma_vs_channel.png` — σ for every channel.
- `Plots_ALL/res_sigma_posX_West_1.png` … `res_sigma_negX_East_5.png`
  — σ on the Z–Y map. Colour scale capped at 3× median; channels above the cap
  (e.g. ch 118) are listed at the bottom of the plot.

---------------------------------------------------------------------------

## 10. One line per step

1. Measure time (half-height) and size (area) of each pulse.
2. Per diffuser: remove the LED's flash-to-flash jitter.
3. All diffusers of a side together: learn and remove the small-pulse delay.
4. Check that time grows with distance at 1/v.
5. Subtract distance → channel delay δ.
6. Average δ over diffusers → calibration constant c_ch.
7. Single-flash scatter → resolution σ.

---------------------------------------------------------------------------

## 11. Suggested figure order for the technote

1. `Plots/Plots-data_..._decoded_filtered/heatmap.png` (raw data) + diffuser map + one `size_map`
2. `size_vs_distance_ALL.png` (why size correction is needed)
3. `size_vs_delay_<d>.png` before/after (it works)
4. `radius_vs_dtexp_<d>.png` + table of all 10 slopes (validation)
5. `radius_vs_delay_<d>.png` and `outlier_<d>.png` for one diffuser
6. `offset_diffuser_vs_channel.png` (consistency across diffusers)
7. `coated_vs_uncoated_offsets.png` + offsets table (result)
8. `residual_sigma_vs_channel.png` (resolution)

## 12. Caveats / open items

- v is assumed (22 cm/ns = c/n, n ≈ 1.36 at 128 nm). Group velocity at 128 nm
  is ~13.4 cm/ns; visible light ~24 cm/ns. LED wavelength still to be
  confirmed with Olivia.
- The size fit uses r/v as input → slope check in Step 4 is not fully
  independent. Independent check: consistency across diffusers (Step 6).
- ch 118: failed Gaussian fit, σ ≈ 300 ns.
- XA (separate section): 16 ns ticks; most XAs see no LED light; the y = +135
  boxes are 60–85 ns late, which sets the sign of the XA slopes.
- Dropped plots (still in the ROOT file; set `kSaveExtraPlots = true` to get
  PNGs back): old copies in `Plots_ALL/unused_old/`.
