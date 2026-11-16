# LED Timing Calibration

Relative timing offsets of the SBND PMTs from the LED diffuser runs 19403-19467 (5 West + 5 East diffusers; each side lights only its own PMTs). A detailed description of the method is in `docs/PMT_timing_calibration_guide.md`.

## Get files prepared
0. Regular setup: `setup_container` --> `sbnd_login` --> `cd {path_to_/larsoft}` --> `source localProducts*/setup` --> `mrbslp`
1. *[optional]* List the raw LED files with samweb (already saved in `LED_data/file_list.txt`):
    ```
    ./get_file_paths.sh
    ```
2. Decode the raw file (`/pnfs/.../data_EventBuilder....root`):
    ```
    lar -c run_decoders_job.fcl -s {raw_file}
    ```
    The output to use is `data_EventBuilder_..._decoded-filtered.root`.
    * To decode many files, use `decode_LED_data.sh` (outputs go to `/pnfs/sbnd/persistent/users/castalyf/LED/LED_data`, set `DECODED_DST` to change):
    ```
    DRYRUN=1 ./decode_LED_data.sh          # only print what would run
    MAXFRAG=3 ./decode_LED_data.sh         # test first: 3 fragments/run
    ./decode_LED_data.sh                   # full set
    ```
3. *[optional]* To look at the waveforms:
    ```
    lar -c /exp/sbnd/data/users/acastill/WaveformCalibration/run_wvfmana_undeco.fcl -s {decoded_file}
    ```
    This outputs `wvfm_ana_undeco.root`, be sure to rename it properly.

## Analysis
1. Make a TTree for each decoded file with `AnalyzerMakinTree.C`. The tree holds the Run #, Event # (LED flash), Channel # (PMT ID 6-305), the full raw waveform (2 ns/tick) and a normalized derivative (≈0.5 on baseline, drops below 0.5 on a pulse's falling edge).
    ```
    MACRO=$PWD/AnalyzerMakinTree.C PLOTS_DIR=$PWD/Plots ./make_trees_PMT.sh
    # single file: root -l -b -q 'AnalyzerMakinTree.C("{decoded_file}")'
    ```
    This gives `outfile_ana-..._decoded_filtered.root` (written to `/pnfs/sbnd/persistent/users/castalyf/LED/PMT_outfiles`) and a quick-look folder `Plots/Plots-..._decoded_filtered/`. The useful plot there is `heatmap.png` (all waveforms of all channels overlaid); `heatmap_zoom.py` makes a cleaner zoomed version.
2. If your trees are somewhere else, change the file path in the constructor of `PMTWaveformTree.h`.
3. Run the main analysis `PMTTimingAnalyzer.C`. In short, it:
    * finds the first LED pulse in each waveform and times it at half of its peak height, and also records the pulse size;
    * removes the flash-to-flash LED jitter and the cable delays;
    * corrects for time walk (small pulses from dim, far-away PMTs come out later) with a fit of time vs pulse size, done over all diffusers of one side;
    * subtracts the expected light travel time from the diffuser (22 cm/ns) and checks the time vs distance slope;
    * averages what is left over the 5 diffusers to get one offset per channel.

    Open `root` and run:
    ```
    .L PMTTimingAnalyzer.C
    PMTWaveformTree t;
    t.Loop();            // all 10 diffusers
    t.Loop("West_3");    // or a single diffuser (quick check only, no global size correction)
    ```
    Outputs:
    * `output_ALL.root`: all histograms and graphs
    * `timing_offsets_ALL.txt`: offset per channel (the calibration result)
    * `Plots_ALL/`: radius vs delay per diffuser, size vs distance/delay (before/after correction), size and residual maps, outliers, offsets per channel and coated vs uncoated comparison

    Settings (timing method, pulse window, thresholds, size-correction order, extra plots) are at the top of `PMTTimingAnalyzer.C`.

## Other files
* `newMethod_PMTWaveformTree.C`: original derivative half-minimum method, kept for reference
* `PMTTimingAnalyzer.py`: PMT position maps and cable-delay maps
