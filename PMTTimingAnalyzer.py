#!/usr/bin/env python3
"""
PMT Cable Delay Heatmap + PMT Position Map Generator
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os

os.makedirs('Plots', exist_ok=True)

# ============================================================================
# DATA
# ============================================================================

pmt_positions = {
    6: {'x': -213.4, 'y': -175, 'z': 27.8742},
    7: {'x': 213.4, 'y': -175, 'z': 27.8742},
    8: {'x': -213.4, 'y': -95, 'z': 27.8742},
    9: {'x': 213.4, 'y': -95, 'z': 27.8742},
    10: {'x': -213.4, 'y': -40, 'z': 27.8742},
    11: {'x': 213.4, 'y': -40, 'z': 27.8742},
    12: {'x': -213.4, 'y': 40, 'z': 27.8742},
    13: {'x': 213.4, 'y': 40, 'z': 27.8742},
    14: {'x': -213.4, 'y': 95, 'z': 27.8742},
    15: {'x': 213.4, 'y': 95, 'z': 27.8742},
    16: {'x': -213.4, 'y': 175, 'z': 27.8742},
    17: {'x': 213.4, 'y': 175, 'z': 27.8742},
    36: {'x': -213.4, 'y': -135, 'z': 57.8742},
    37: {'x': 213.4, 'y': -135, 'z': 57.8742},
    38: {'x': -213.4, 'y': 0, 'z': 57.8742},
    39: {'x': 213.4, 'y': 0, 'z': 57.8742},
    40: {'x': -213.4, 'y': 135, 'z': 57.8742},
    41: {'x': 213.4, 'y': 135, 'z': 57.8742},
    60: {'x': -213.4, 'y': -175, 'z': 87.8742},
    61: {'x': 213.4, 'y': -175, 'z': 87.8742},
    62: {'x': -213.4, 'y': -95, 'z': 87.8742},
    63: {'x': 213.4, 'y': -95, 'z': 87.8742},
    64: {'x': -213.4, 'y': -40, 'z': 87.8742},
    65: {'x': 213.4, 'y': -40, 'z': 87.8742},
    66: {'x': -213.4, 'y': 40, 'z': 87.8742},
    67: {'x': 213.4, 'y': 40, 'z': 87.8742},
    68: {'x': -213.4, 'y': 95, 'z': 87.8742},
    69: {'x': 213.4, 'y': 95, 'z': 87.8742},
    70: {'x': -213.4, 'y': 175, 'z': 87.8742},
    71: {'x': 213.4, 'y': 175, 'z': 87.8742},
    84: {'x': -213.4, 'y': -175, 'z': 161.158},
    85: {'x': 213.4, 'y': -175, 'z': 161.158},
    86: {'x': -213.4, 'y': -95, 'z': 161.158},
    87: {'x': 213.4, 'y': -95, 'z': 161.158},
    88: {'x': -213.4, 'y': -40, 'z': 161.158},
    89: {'x': 213.4, 'y': -40, 'z': 161.158},
    90: {'x': -213.4, 'y': 40, 'z': 161.158},
    91: {'x': 213.4, 'y': 40, 'z': 161.158},
    92: {'x': -213.4, 'y': 95, 'z': 161.158},
    93: {'x': 213.4, 'y': 95, 'z': 161.158},
    94: {'x': -213.4, 'y': 175, 'z': 161.158},
    95: {'x': 213.4, 'y': 175, 'z': 161.158},
    114: {'x': -213.4, 'y': -135, 'z': 191.158},
    115: {'x': 213.4, 'y': -135, 'z': 191.158},
    116: {'x': -213.4, 'y': 0, 'z': 191.158},
    117: {'x': 213.4, 'y': 0, 'z': 191.158},
    118: {'x': -213.4, 'y': 135, 'z': 191.158},
    119: {'x': 213.4, 'y': 135, 'z': 191.158},
    138: {'x': -213.4, 'y': -175, 'z': 221.158},
    139: {'x': 213.4, 'y': -175, 'z': 221.158},
    140: {'x': -213.4, 'y': -95, 'z': 221.158},
    141: {'x': 213.4, 'y': -95, 'z': 221.158},
    142: {'x': -213.4, 'y': -40, 'z': 221.158},
    143: {'x': 213.4, 'y': -40, 'z': 221.158},
    144: {'x': -213.4, 'y': 40, 'z': 221.158},
    145: {'x': 213.4, 'y': 40, 'z': 221.158},
    146: {'x': -213.4, 'y': 95, 'z': 221.158},
    147: {'x': 213.4, 'y': 95, 'z': 221.158},
    148: {'x': -213.4, 'y': 175, 'z': 221.158},
    149: {'x': 213.4, 'y': 175, 'z': 221.158},
    162: {'x': -213.4, 'y': -175, 'z': 288.242},
    163: {'x': 213.4, 'y': -175, 'z': 288.242},
    164: {'x': -213.4, 'y': -95, 'z': 288.242},
    165: {'x': 213.4, 'y': -95, 'z': 288.242},
    166: {'x': -213.4, 'y': -40, 'z': 288.242},
    167: {'x': 213.4, 'y': -40, 'z': 288.242},
    168: {'x': -213.4, 'y': 40, 'z': 288.242},
    169: {'x': 213.4, 'y': 40, 'z': 288.242},
    170: {'x': -213.4, 'y': 95, 'z': 288.242},
    171: {'x': 213.4, 'y': 95, 'z': 288.242},
    172: {'x': -213.4, 'y': 175, 'z': 288.242},
    173: {'x': 213.4, 'y': 175, 'z': 288.242},
    192: {'x': -213.4, 'y': -135, 'z': 318.242},
    193: {'x': 213.4, 'y': -135, 'z': 318.242},
    194: {'x': -213.4, 'y': 0, 'z': 318.242},
    195: {'x': 213.4, 'y': 0, 'z': 318.242},
    196: {'x': -213.4, 'y': 135, 'z': 318.242},
    197: {'x': 213.4, 'y': 135, 'z': 318.242},
    216: {'x': -213.4, 'y': -175, 'z': 348.242},
    217: {'x': 213.4, 'y': -175, 'z': 348.242},
    218: {'x': -213.4, 'y': -95, 'z': 348.242},
    219: {'x': 213.4, 'y': -95, 'z': 348.242},
    220: {'x': -213.4, 'y': -40, 'z': 348.242},
    221: {'x': 213.4, 'y': -40, 'z': 348.242},
    222: {'x': -213.4, 'y': 40, 'z': 348.242},
    223: {'x': 213.4, 'y': 40, 'z': 348.242},
    224: {'x': -213.4, 'y': 95, 'z': 348.242},
    225: {'x': 213.4, 'y': 95, 'z': 348.242},
    226: {'x': -213.4, 'y': 175, 'z': 348.242},
    227: {'x': 213.4, 'y': 175, 'z': 348.242},
    240: {'x': -213.4, 'y': -175, 'z': 421.526},
    241: {'x': 213.4, 'y': -175, 'z': 421.526},
    242: {'x': -213.4, 'y': -95, 'z': 421.526},
    243: {'x': 213.4, 'y': -95, 'z': 421.526},
    244: {'x': -213.4, 'y': -40, 'z': 421.526},
    245: {'x': 213.4, 'y': -40, 'z': 421.526},
    246: {'x': -213.4, 'y': 40, 'z': 421.526},
    247: {'x': 213.4, 'y': 40, 'z': 421.526},
    248: {'x': -213.4, 'y': 95, 'z': 421.526},
    249: {'x': 213.4, 'y': 95, 'z': 421.526},
    250: {'x': -213.4, 'y': 175, 'z': 421.526},
    251: {'x': 213.4, 'y': 175, 'z': 421.526},
    270: {'x': -213.4, 'y': -135, 'z': 451.526},
    271: {'x': 213.4, 'y': -135, 'z': 451.526},
    272: {'x': -213.4, 'y': 0, 'z': 451.526},
    273: {'x': 213.4, 'y': 0, 'z': 451.526},
    274: {'x': -213.4, 'y': 135, 'z': 451.526},
    275: {'x': 213.4, 'y': 135, 'z': 451.526},
    294: {'x': -213.4, 'y': -175, 'z': 481.526},
    295: {'x': 213.4, 'y': -175, 'z': 481.526},
    296: {'x': -213.4, 'y': -95, 'z': 481.526},
    297: {'x': 213.4, 'y': -95, 'z': 481.526},
    298: {'x': -213.4, 'y': -40, 'z': 481.526},
    299: {'x': 213.4, 'y': -40, 'z': 481.526},
    300: {'x': -213.4, 'y': 40, 'z': 481.526},
    301: {'x': 213.4, 'y': 40, 'z': 481.526},
    302: {'x': -213.4, 'y': 95, 'z': 481.526},
    303: {'x': 213.4, 'y': 95, 'z': 481.526},
    304: {'x': -213.4, 'y': 175, 'z': 481.526},
    305: {'x': 213.4, 'y': 175, 'z': 481.526},
}

cable_delays = {
    250: 5.44, 304: 1.22, 274: 4.1, 248: 0.25, 302: 5.54, 246: 3.18,
    300: 3.63, 272: 2.75, 244: 3.22, 298: 3.9, 242: 3.8, 296: 0,
    270: 0.76, 240: 1.8, 294: 0.1, 172: 1.95, 226: 1.58, 196: 0.2,
    170: 1.9, 224: 3.77, 168: -0.64, 222: -0.9, 194: 2.09, 166: 2.3,
    220: 3.85, 164: 1.42, 218: -0.35, 192: 0.16, 162: -0.6, 216: 5.61,
    94: 1.75, 148: -0.93, 118: -1.9, 92: 1.2, 146: 0.63, 90: 5.15,
    144: 3.25, 116: -2.67, 88: 2.13, 142: -2, 86: 0, 140: -0.05,
    114: 0.52, 84: -1.45, 138: 1.63, 16: 1.97, 70: 5.67, 40: 5.45,
    14: 3.7, 68: 5.62, 12: 1, 66: 0.65, 38: 1.93, 10: -2.25,
    64: 5.63, 8: -10.22, 62: -4.05, 36: 0.45, 6: -2.99, 60: -3.85,
    71: 0.5, 17: -3.32, 41: -1.74, 69: 2.38, 15: 0.3, 67: 2,
    13: -2.74, 39: 2.2, 65: -1.15, 11: 5.62, 63: 0.43, 9: -4.26,
    37: -3.83, 61: -1.53, 7: -3.07, 149: 2.4, 95: -1, 119: 4.2,
    147: 1.4, 93: -0.7, 145: 4.05, 91: 4.07, 117: 2.1, 143: 1.6,
    89: 3.81, 141: -3.05, 87: 0.08, 115: -0.68, 139: -1.76, 85: -11.41,
    227: 4.16, 173: 1.71, 197: 3.28, 225: -2.26, 171: -2.45, 223: 1.91,
    169: 0.42, 195: 4.04, 221: 1.25, 167: 1.07, 219: -1.25, 165: 2,
    193: -1.7, 217: 1.27, 163: 1.4, 305: 1.27, 251: 4.27, 275: 1.84,
    303: 5.18, 249: -1.77, 301: 1.19, 247: 3.09, 273: -0.21, 299: 3.63,
    245: 1.75, 297: 5.25, 243: 5.28, 271: 0.75, 295: -2.98, 241: 1.24
}

# Diffuser positions (Y, Z) — X=0 on CPA
# West diffusers face the West APA (x = -213.4)
west_diffusers = {
    1: {'y':   86.36, 'z': 123.42},
    2: {'y': -118.97, 'z': 123.42},
    3: {'y':  -12.95, 'z': 246.23},
    4: {'y':   86.60, 'z': 380.95},
    5: {'y': -117.46, 'z': 380.80},
}
# East diffusers face the East APA (x = +213.4)
east_diffusers = {
    1: {'y':   86.82, 'z': 131.62},
    2: {'y': -116.81, 'z': 131.62},
    3: {'y':  -13.33, 'z': 261.65},
    4: {'y':   85.40, 'z': 388.95},
    5: {'y': -119.98, 'z': 389.05},
}

# Diffuser label colours — one per index 1-5
diff_colors = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3']

# ============================================================================
# HELPER: draw diffuser overlays on an axis
# ============================================================================
def draw_diffusers(ax, diffusers, label_side='right'):
    """
    Plot diffuser positions as solid coloured circles with index labels.
    label_side: 'right' or 'left' for number label placement.
    """
    ha = 'left' if label_side == 'right' else 'right'
    x_offset = 14 if label_side == 'right' else -14

    for idx, pos in diffusers.items():
        color = diff_colors[idx - 1]
        # Solid circle, slightly larger than a PMT dot, with white edge
        ax.scatter(pos['z'], pos['y'],
                   s=350, marker='o',
                   facecolors=color, edgecolors='white',
                   linewidths=1.8, zorder=5)
        # Bold number label offset to the side
        ax.annotate(
            str(idx),
            xy=(pos['z'], pos['y']),
            xytext=(pos['z'] + x_offset, pos['y']),
            fontsize=14, fontweight='bold', color=color,
            ha=ha, va='center', zorder=6,
            path_effects=[pe.withStroke(linewidth=2, foreground='white')]
        )

# ============================================================================
# HELPER: draw PMT positions with channel-number labels on an axis
# ============================================================================
def draw_pmts_labelled(ax, pmts_dict, scatter_kwargs=None):
    """
    Scatter all PMTs in pmts_dict and annotate each with its channel number.
    Returns the scatter object (for colorbar if needed).
    """
    if scatter_kwargs is None:
        scatter_kwargs = {}

    zs = [p['z'] for p in pmts_dict.values()]
    ys = [p['y'] for p in pmts_dict.values()]
    chs = list(pmts_dict.keys())

    sc = ax.scatter(zs, ys, zorder=3, **scatter_kwargs)

    for ch, z, y in zip(chs, zs, ys):
        ax.annotate(
            str(ch),
            xy=(z, y), xytext=(0, 0), textcoords='offset points',
            fontsize=6.5, ha='center', va='center',
            color='black', fontweight='bold', zorder=4,
            path_effects=[pe.withStroke(linewidth=1.5, foreground='white')]
        )
    return sc

# ============================================================================
# SEPARATE PMTs BY SIDE
# ============================================================================
west_pmts = {ch: p for ch, p in pmt_positions.items() if p['x'] < 0}
east_pmts = {ch: p for ch, p in pmt_positions.items() if p['x'] > 0}

west_z_d = [pmt_positions[ch]['z'] for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] < 0]
west_y_d = [pmt_positions[ch]['y'] for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] < 0]
west_d   = [cable_delays[ch]       for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] < 0]

east_z_d = [pmt_positions[ch]['z'] for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] > 0]
east_y_d = [pmt_positions[ch]['y'] for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] > 0]
east_d   = [cable_delays[ch]       for ch in cable_delays if ch in pmt_positions and pmt_positions[ch]['x'] > 0]

vmin, vmax = -12, 6

# ============================================================================
# PLOT 1 & 2 — CABLE DELAY HEATMAPS (original, unchanged)
# ============================================================================
print("\nGenerating cable delay heatmaps...")

for side_label, zs, ys, ds, title, fname in [
    ('West', west_z_d, west_y_d, west_d,
     'Cable Delays - West PMTs (x = -213.4 cm)', 'Plots/cable_delays_west.png'),
    ('East', east_z_d, east_y_d, east_d,
     'Cable Delays - East PMTs (x = +213.4 cm)', 'Plots/cable_delays_east.png'),
]:
    fig, ax = plt.subplots(figsize=(12, 8))
    sc = ax.scatter(zs, ys, c=ds, s=700, cmap='RdYlBu_r',
                    edgecolors='black', linewidths=1.0, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Cable Delay (ns)', fontsize=18)
    ax.set_xlabel('Z Position (cm)', fontsize=18, fontweight='bold')
    ax.set_ylabel('Y Position (cm)', fontsize=18, fontweight='bold')
    ax.set_title(title, fontsize=24, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=1.2)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    print(f'✓ Saved: {fname}')
    plt.close()

# ============================================================================
# PLOT 3 — WEST APA PMT MAP + WEST DIFFUSERS
# West APA is at x = -213.4.
# West diffusers (x=0, CPA) face toward the West APA.
# Both share the same Y-Z plane, so diffuser Y/Z maps directly onto the PMT map.
# ============================================================================
print("\nGenerating PMT position maps with diffuser overlays...")

fig, ax = plt.subplots(figsize=(14, 8))

# PMTs: filled blue circles, large enough to fit channel number inside
sc = draw_pmts_labelled(
    ax, west_pmts,
    scatter_kwargs=dict(s=800, color='tomato', edgecolors='darkred',
                        linewidths=0.8, alpha=0.85)
)

# Diffusers: open stars, coloured by index
draw_diffusers(ax, west_diffusers, label_side='right')

ax.set_xlim(-20, 520)
ax.set_ylim(-220, 240)
ax.set_xlabel('Z Position (cm)', fontsize=18, fontweight='bold')
ax.set_ylabel('Y Position (cm)', fontsize=18, fontweight='bold')
ax.set_title('PMT Layout – West APA (x = −213.4 cm)\nwith West Diffuser Positions (CPA, x = 0)',
             fontsize=18, fontweight='bold', pad=12)
ax.grid(True, alpha=0.25, linestyle='--', linewidth=1.0)
ax.tick_params(axis='both', which='major', labelsize=13)

# Legend: one entry per diffuser index
from matplotlib.lines import Line2D
pmt_handle = Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='tomato', markeredgecolor='darkred',
                    markersize=15, label='PMT (channel # inside)')
diff_handle = Line2D([0], [0], marker='o', color='w',
                     markerfacecolor='none', markeredgecolor='#444',
                     markeredgewidth=2.0, markersize=14, label='Diffuser (# labelled)')
ax.legend(handles=[pmt_handle, diff_handle],
          fontsize=12, loc='upper right', framealpha=0.9)

plt.tight_layout()
plt.savefig('Plots/pmt_map_west.png', dpi=300, bbox_inches='tight')
print('✓ Saved: Plots/pmt_map_west.png')
plt.close()

# ============================================================================
# PLOT 4 — EAST APA PMT MAP + EAST DIFFUSERS
# East APA is at x = +213.4.
# East diffusers face toward the East APA.
# ============================================================================
fig, ax = plt.subplots(figsize=(14, 8))

sc = draw_pmts_labelled(
    ax, east_pmts,
    scatter_kwargs=dict(s=800, color='steelblue', edgecolors='navy',
                        linewidths=0.8, alpha=0.85)
)

draw_diffusers(ax, east_diffusers, label_side='right')

ax.set_xlim(-20, 520)
ax.set_ylim(-220, 240)
ax.set_xlabel('Z Position (cm)', fontsize=18, fontweight='bold')
ax.set_ylabel('Y Position (cm)', fontsize=18, fontweight='bold')
ax.set_title('PMT Layout – East APA (x = +213.4 cm)\nwith East Diffuser Positions (CPA, x = 0)',
             fontsize=18, fontweight='bold', pad=12)
ax.grid(True, alpha=0.25, linestyle='--', linewidth=1.0)
ax.tick_params(axis='both', which='major', labelsize=13)

pmt_handle = Line2D([0], [0], marker='o', color='w',
                    markerfacecolor='steelblue', markeredgecolor='navy',
                    markersize=15, label='PMT (channel # inside)')
diff_handle = Line2D([0], [0], marker='o', color='w',
                     markerfacecolor='none', markeredgecolor='#444',
                     markeredgewidth=2.0, markersize=14, label='Diffuser (# labelled)')
ax.legend(handles=[pmt_handle, diff_handle],
          fontsize=12, loc='upper right', framealpha=0.9)

plt.tight_layout()
plt.savefig('Plots/pmt_map_east.png', dpi=300, bbox_inches='tight')
print('✓ Saved: Plots/pmt_map_east.png')
plt.close()

# ============================================================================
# SUMMARY
# ============================================================================
delays = list(cable_delays.values())
print(f"\n{'='*70}")
print("SUMMARY")
print(f"{'='*70}")
print(f"Total PMTs with delays: {len(cable_delays)}")
print(f"  West PMTs: {len(west_d)}   East PMTs: {len(east_d)}")
print(f"\nDelay statistics:")
print(f"  Overall: {np.mean(delays):.2f} ± {np.std(delays):.2f} ns  "
      f"(range: {min(delays):.2f} to {max(delays):.2f} ns)")
print(f"  West:    {np.mean(west_d):.2f} ± {np.std(west_d):.2f} ns")
print(f"  East:    {np.mean(east_d):.2f} ± {np.std(east_d):.2f} ns")
print(f"\n{'='*70}")
print("COMPLETE — 4 plots saved to Plots/")
print(f"{'='*70}\n")