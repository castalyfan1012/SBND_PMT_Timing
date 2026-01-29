# PMT Timing Calibration

> **Note:** This analysis script is designed for internal use within the SBND experiment.

## Prerequisites

Before beginning the analysis, ensure you have properly set up the SL7 container with `sbndcode`:

0. Enter the SBND virtual machine / terminal prompt.
1. Set up the SL7 container.
2. Set up `sbndcode` and `ROOT` with the proper versions.
3. Navigate to your working directory: `cd {path_to_/larsoft}`
4. Source the setup: `source localProducts*/setup`
5. Initialize MRB: `mrbslp`

## Getting Started

## File Preparation

Navigate to your working directory e.g.
```
cd {path_to_your_repo}/SBND_PMT_timing/LED_data
```

### 0. Creating a File List

Generate a list of input files using the provided script:
```
# Edit get_file_paths.sh to set your RUN_RANGE (e.g., "19403-19467")
# then run:
./get_file_paths.sh
```

This will create a `filelist.txt` containing the paths to your raw data files.


### 1. Decode Raw Files

Decode the raw data file (`/pnfs/.../data_EventBuilder....root`):
```
lar -c run_decoders_job.fcl -s {raw_file}
```

**Output:** `data_EventBuilder_..._decoded-filtered.root`

### 2. Retrieve Waveform Information (Optional)

Extract waveform data from the decoded file:
```
lar -c /exp/sbnd/data/users/acastill/WaveformCalibration/run_wvfmana_undeco.fcl -s {decoded_file}
```

**Output:** `wvfm_ana_undeco.root` (rename as needed)

### 3. Extract Channel Information (Optional)

Use the Jupyter notebook `PMT_waveform_test.ipynb` to identify channels based on PMT responses.

## Analysis Workflow

### Step 1: Create TTree from Decoded Data

Navigate to the analysis directory:
```
cd /exp/sbnd/data/users/castalyf/PMT_timing
```

Generate a TTree containing the following branches:
- Run number (data collection session ID)
- Event number (LED flash number)
- Channel number (PMT ID: 6-305)
- Waveform start tick
- Waveform vector (100 raw ADC samples)
- Derivative start tick
- Derivative vector (80 normalized values)

Execute the analyzer:
```
root 'AnalyzerMakinTree.C("{decoded_root}")'
```

**Input:** `data_EventBuilder..._decoded-filtered.root`

**Outputs:**
- `outfile_ana-..._decoded-filtered.root`
- Waveform/pulse plots in `plots-..._decoded-filtered/` directory

### Step 2: Configure Analysis Tree

Edit `PMTWaveformTree.h` and update the input filename to match the output from Step 1 (`outfile_ana-..._decoded-filtered.root`).

### Step 3: Run Timing Analysis

The main analysis script `PMTWaveformAnalyzer.C` performs:
- PMT waveform timing extraction via derivative half-minimum method
- Per-channel time delay corrections
- Event-wise averaging of PMT times
- Per-PMT delay statistics with truncated (±1 RMS) distributions
- Gaussian fitting of delay distributions
- Geometric mapping of timing statistics (East/West detector separation)
- Radial delay analysis from detector origin

Execute the analysis in a ROOT session (you may just type `root` to enter the terminal):

```cpp
.L PMTWaveformAnalyzer.C
PMTWaveformTree t;
t.Loop();
```

**Outputs:**
- ROOT file: `output_waveforms_combined_LED1.root`
- Per-channel absolute and relative time histograms
- Average event time histogram
- 2D heatmaps: mean/RMS PMT delays (East/West)
- 2D heatmaps: Gaussian fit parameters
- PMT delay vs. radius scatter plots (East/West)
- Combined radius vs. delay plot
- PNG plots saved to `Plots/` directory

Summary statistics are printed to stdout during execution.
