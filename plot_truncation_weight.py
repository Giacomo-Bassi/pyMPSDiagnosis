import re
import sys
import argparse
import numpy as np
import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

# Show tick values as-is: no offset (e.g. "+1.727e3"), no automatic rescaling
matplotlib.rcParams['axes.formatter.useoffset'] = False
matplotlib.rcParams['axes.formatter.use_mathtext'] = False
matplotlib.rcParams['axes.formatter.limits'] = (-99, 99)  # never switch to sci notation on linear axes



def load_from_hdf5(filepath):
    with h5py.File(filepath, 'r') as f:
        data        = f['truncated_weights'][:]
        site_labels = [s.decode() for s in f['site_labels'][:]]
        sweep_nums  = f['sweep_numbers'][:].tolist()
        energies    = f['sweep_energies'][:].tolist() if 'sweep_energies' in f else [None] * len(sweep_nums)
    return data, site_labels, sweep_nums, energies


def site_label_to_index(label):
    """Extract the left site index from a label like -3-4-> or <-3-4-."""
    nums = re.findall(r'\d+', label)
    return int(nums[0]) if nums else None


def plot_truncated_weights(h5_file, savepath, plot_name='truncated_weights_plot.png', threshold=0):
    data, site_labels, sweep_nums, energies = load_from_hdf5(h5_file)
    # data shape: (n_sites, n_sweeps)
    zero_level = np.min(data[np.isfinite(data) & (data > threshold)]) * 0.5
    print(f"Zero level for invalid points: {zero_level:.2e}")
    # Split into forward and backward labels (preserve original row index)
    fwd_rows  = [(i, lbl) for i, lbl in enumerate(site_labels) if '->' in lbl]
    bwd_rows  = [(i, lbl) for i, lbl in enumerate(site_labels) if '<-' in lbl]

    # Sort forward ascending by bond index, backward descending (right-to-left order)
    fwd_rows.sort(key=lambda t: site_label_to_index(t[1]))
    bwd_rows.sort(key=lambda t: site_label_to_index(t[1]), reverse=True)

    n_fwd = len(fwd_rows)
    n_bwd = len(bwd_rows)
    n_total = n_fwd + n_bwd          # x axis: 1 .. n_total

    # Build x positions and tick labels
    # Forward:  positions 1 .. n_fwd
    # Backward: positions n_fwd+1 .. n_total
    x_fwd   = np.arange(1, n_fwd + 1)
    x_bwd   = np.arange(n_fwd + 1, n_total + 1)

    tick_positions = np.concatenate([x_fwd, x_bwd])
    tick_labels    = [lbl for _, lbl in fwd_rows] + [lbl for _, lbl in bwd_rows]

    fwd_data_idx  = [i for i, _ in fwd_rows]
    bwd_data_idx  = [i for i, _ in bwd_rows]

    n_sweeps = len(sweep_nums)
    colors   = cm.viridis(np.linspace(0, 1, n_sweeps))

    fig, ax = plt.subplots(figsize=(max(12, n_total * 0.5), 5))

    for col_idx, (sweep, color, energy) in enumerate(zip(sweep_nums, colors, energies)):
        # forward half
        lbl = f'Sweep {sweep}' if energy is None else f'Sweep {sweep}  (E = {energy:.6f})'
        y_fwd = data[fwd_data_idx, col_idx]
        valid = np.isfinite(y_fwd) & (y_fwd > threshold)
        ax.semilogy(x_fwd[valid], y_fwd[valid],
                    marker='o', markersize=3, linewidth=1.2,
                    color=color, label=lbl)
        ax.plot(x_fwd[~valid], zero_level * np.ones_like(x_fwd[~valid]),  # plot invalid points at half the smallest valid value
                marker='x', markersize=5, linestyle='None',
                color='red')
        # backward half – same color, dashed, no extra legend entry
        y_bwd = data[bwd_data_idx, col_idx]
        valid = np.isfinite(y_bwd) & (y_bwd > threshold)
        ax.semilogy(x_bwd[valid], y_bwd[valid],
                    marker='o', markersize=3, linewidth=1.2,
                    color=color)
        ax.plot(x_bwd[~valid], zero_level * np.ones_like(x_bwd[~valid]),  # plot invalid points at half the smallest valid value
                marker='x', markersize=5, linestyle='None',
                color='red')

    # Vertical separator between forward and backward half
    if n_fwd > 0 and n_bwd > 0:
        ax.axvline(n_fwd + 0.5, color='gray', linestyle=':', linewidth=1.0)
        ax.text(n_fwd * 0.5 + 0.5, ax.get_ylim()[0], '→ forward',
                ha='center', va='bottom', fontsize=12, color='gray')
        ax.text(n_fwd + n_bwd * 0.5 + 0.5, ax.get_ylim()[0], '← backward',
                ha='center', va='bottom', fontsize=12, color='gray')

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=60, fontsize=12, 
                       ha='right', va='top', fontweight='bold')
    ax.set_xlim(0.5, n_total + 0.5)
    ax.set_xlabel('Bond')
    ax.set_ylabel('Truncated weight')
    ax.set_title('DMRG truncated weights per sweep')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)

    fig.tight_layout()

    output = os.path.join(savepath, plot_name)
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {output}")


def plot_max_truncation_per_sweep(h5_file, savepath, plot_name='max_truncation_per_sweep.png'):
    """Plot the maximum truncated weight across all sites for each sweep."""
    data, site_labels, sweep_nums, energies = load_from_hdf5(h5_file)
    # data shape: (n_sites, n_sweeps)


    def safe_max(rows, col):
        vals = data[rows, col]
        valid = vals[np.isfinite(vals) & (vals > 0)]
        return valid.max() if len(valid) > 0 else np.nan

    max_all = np.array([safe_max(list(range(data.shape[0])), c) for c in range(len(sweep_nums))])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax2 = ax.twinx()
    x = np.array(sweep_nums)

    ax.semilogy(x, max_all, marker='D', linewidth=1.5, markersize=6,
                color='black', label='Max truncated weight')
    ax2.plot(x, energies, marker='o', linestyle='-', color='steelblue', label='Sweep energy')

    ax.set_xlabel('Sweep')
    ax.set_ylabel('Max truncated weight')
    ax2.set_ylabel('Sweep energy')
    ax.set_title('Maximum truncated weight per sweep')
    ax.set_xticks(x)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper right')
    ax.grid(True, which='both', linestyle='--', alpha=0.4)

    fig.tight_layout()

    output = os.path.join(savepath, plot_name)
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {output}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot DMRG truncated weights from an HDF5 file.')
    parser.add_argument('--h5_file', required=True,
                        help='Path to the HDF5 file produced by parse_truncation_weight.py')
    parser.add_argument('--output_root', '-o', default='truncated_weights_plot.png', required=False,
                        help='Root name for output files (e.g. plot). If omitted, show interactively.')
    parser.add_argument('--savepath', required=True,
                        help='Directory to save the plots.')
            
    args = parser.parse_args()

    plot_truncated_weights(args.h5_file, args.savepath, plot_name=args.output_root+"_truncated_weights_plot.png")
    # plot_max_truncation_per_sweep(args.h5_file, args.savepath, plot_name=args.output_root+"_max_truncation_per_sweep.png")

