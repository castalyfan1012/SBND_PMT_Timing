import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

fin, outbase = sys.argv[1], sys.argv[2]

# Same binning as AnalyzerMakinTree.C hHeat: 1000 bins in [0,5000], 800 bins in [12000,16000]
xedges = np.linspace(0, 5000, 1001)
yedges = np.linspace(12000, 16000, 801)
H = np.zeros((1000, 800))

t = uproot.open(fin)["PMTWaveformTree"]
nwf = 0
for arr in t.iterate(["Waveform"], step_size="200 MB", library="ak"):
    for wf in arr["Waveform"]:
        w = np.asarray(wf)[:5000]
        h, _, _ = np.histogram2d(np.arange(len(w)), w, bins=[xedges, yedges])
        H += h
        nwf += 1
print("waveforms:", nwf)

fig, ax = plt.subplots(figsize=(8, 4))
Hm = np.ma.masked_where(H.T == 0, H.T)
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad("white")
pc = ax.pcolormesh(xedges, yedges, Hm, cmap=cmap, norm=LogNorm(vmin=1, vmax=H.max()))
ax.set_xlim(0, 2500)
ax.set_ylim(12500, 15300)
ax.minorticks_on()
ax.set_xlabel("Time [tick]")
ax.set_ylabel("ADC")
ax.set_title("PMT waveforms, all channels (run 19410)")
fig.colorbar(pc, ax=ax, label="Entries")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{outbase}.{ext}", dpi=200)
print("saved", outbase)
