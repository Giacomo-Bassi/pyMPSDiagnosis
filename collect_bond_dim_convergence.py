"""
collect_bond_dim_convergence.py

Reads multiple HDF5 files produced by parse_truncation_weight.py (one per bond
dimension), extracts the sweep energies and maximum truncated weight per sweep,
and saves a new HDF5 file with 2-D arrays indexed by (sweep, bond_dim).

Usage
-----
    python collect_bond_dim_convergence.py \
        --h5_files  path/to/calc_m64_hf.QCMaquis_truncated_weights.h5 \
                    path/to/calc_m128_hf.QCMaquis_truncated_weights.h5 \
                    path/to/calc_m256_hf.QCMaquis_truncated_weights.h5 \
        --output    convergence.h5

Bond dimensions are extracted automatically from the filename via the '_m{N}_'
pattern (e.g. '..._m1024_hf.QCMaquis_truncated_weights.h5' → 1024).
You can also override them explicitly with --bond_dims.
"""

import re
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


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_h5(filepath):
    """Return (data, site_labels, sweep_nums, energies) from an h5 file."""
    with h5py.File(filepath, 'r') as f:
        data        = f['truncated_weights'][:]
        site_labels = [s.decode() for s in f['site_labels'][:]]
        sweep_nums  = f['sweep_numbers'][:].tolist()
        energies    = (f['sweep_energies'][:].tolist()
                       if 'sweep_energies' in f
                       else [float('nan')] * len(sweep_nums))
    return data, site_labels, sweep_nums, energies


def max_truncation_per_sweep(data, n_sweeps):
    """Return array of shape (n_sweeps,) with the max truncated weight per sweep."""
    result = np.full(n_sweeps, np.nan)
    for col in range(n_sweeps):
        vals = data[:, col]
        valid = vals[np.isfinite(vals)]
        if len(valid) > 0:
            result[col] = valid.max()
    return result

def sum_truncation_per_sweep(data, n_sweeps):
    """Return array of shape (n_sweeps,) with the sum of truncated weights per sweep."""
    result = np.full(n_sweeps, np.nan)
    for col in range(n_sweeps):
        vals = data[:, col]
        valid = vals[np.isfinite(vals)]
        if len(valid) > 0:
            result[col] = valid.sum()
    return result

def extract_bond_dim(filepath):
    """Extract bond dimension from a path containing '_m{N}_'."""
    m = re.search(r'_m(\d+)_', filepath)
    if m:
        return int(m.group(1))
    raise ValueError(
        f"Cannot extract bond dimension from '{filepath}'. "
        "Expected a '_m<N>_' substring. Use --bond_dims to set them explicitly."
    )


# ---------------------------------------------------------------------------
# main logic
# ---------------------------------------------------------------------------

def collect(h5_files, bond_dims, output):
    records = []   # list of (bond_dim, sweep_nums, energies_arr, max_trunc_arr, sum_trunc_arr)

    for path, bdim in zip(h5_files, bond_dims):
        data, _, sweep_nums, energies = load_h5(path)
        n_sweeps = len(sweep_nums)
        max_trunc = max_truncation_per_sweep(data, n_sweeps)
        sum_trunc = sum_truncation_per_sweep(data, n_sweeps)
        energies_arr = np.array(energies, dtype=float)
        records.append((bdim, sweep_nums, energies_arr, max_trunc, sum_trunc))

    # Union of all sweep numbers, sorted
    all_sweeps = sorted({s for _, sw, _, _, _ in records for s in sw})
    all_bdims  = sorted({bd for bd, *_ in records})

    n_sweeps = len(all_sweeps)
    n_bdims  = len(all_bdims)

    sweep_idx = {s: i for i, s in enumerate(all_sweeps)}
    bdim_idx  = {b: i for i, b in enumerate(all_bdims)}

    energy_matrix    = np.full((n_sweeps, n_bdims), np.nan)
    max_trunc_matrix = np.full((n_sweeps, n_bdims), np.nan)
    sum_trunc_matrix = np.full((n_sweeps, n_bdims), np.nan)

    for bdim, sweep_nums, energies_arr, max_trunc, sum_trunc in records:
        bi = bdim_idx[bdim]
        for k, s in enumerate(sweep_nums):
            si = sweep_idx[s]
            energy_matrix[si, bi]    = energies_arr[k]
            max_trunc_matrix[si, bi] = max_trunc[k]
            sum_trunc_matrix[si, bi] = sum_trunc[k]

    with h5py.File(output, 'w') as f:
        f.create_dataset('sweep_numbers',          data=np.array(all_sweeps, dtype=int))
        f.create_dataset('bond_dimensions',        data=np.array(all_bdims,  dtype=int))
        f.create_dataset('sweep_energies',         data=energy_matrix)
        f.create_dataset('max_truncated_weight',   data=max_trunc_matrix)
        f.create_dataset('sum_truncated_weight',   data=sum_trunc_matrix)
        f.attrs['description'] = (
            'DMRG convergence data: rows=sweeps, cols=bond_dimensions. '
            'Datasets: sweep_energies, max_truncated_weight, sum_truncated_weight.'
        )
        f.attrs['sweep_numbers_info'] = 'shape (n_sweeps,)'
        f.attrs['bond_dimensions_info'] = 'shape (n_bdims,)'
        f.attrs['matrix_shape_info'] = '(n_sweeps, n_bdims)'

    print(f"\nSaved to {output}")
    print(f"  sweep_numbers        : {all_sweeps}")
    print(f"  bond_dimensions      : {all_bdims}")


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def load_convergence_h5(filepath):
    """Load the convergence HDF5 produced by collect()."""
    with h5py.File(filepath, 'r') as f:
        sweep_nums = f['sweep_numbers'][:].tolist()
        bond_dims  = f['bond_dimensions'][:].tolist()
        energies   = f['sweep_energies'][:]
        max_trunc  = f['max_truncated_weight'][:]
        sum_trunc  = f['sum_truncated_weight'][:]
    return sweep_nums, bond_dims, energies, max_trunc, sum_trunc


def plot_energy_convergence(conv_h5, savepath, plot_name='energy_convergence.png'):
    """x = sweep number, y = energy, one line per bond dimension."""
    sweep_nums, bond_dims, energies, _, _ = load_convergence_h5(conv_h5)
    x = np.array(sweep_nums)
    colors = cm.viridis(np.linspace(0, 1, len(bond_dims)))

    fig, ax = plt.subplots(figsize=(7, 4))
    for bi, (bdim, color) in enumerate(zip(bond_dims, colors)):
        y = energies[:, bi]
        valid = np.isfinite(y)
        ax.plot(x[valid], y[valid], marker='o', markersize=5, linewidth=1.4,
                color=color, label=f'm = {bdim}')

    ax.set_xlabel('Sweep')
    ax.set_ylabel('Energy')
    ax.set_title('Sweep energy convergence vs bond dimension')
    ax.set_xticks(x)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.4)
    fig.tight_layout()

    out = os.path.join(savepath, plot_name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {out}")


def plot_max_trunc_convergence(conv_h5, savepath, plot_name='max_trunc_convergence.png'):
    """x = sweep number, y = max truncated weight (log scale), one line per bond dimension."""
    sweep_nums, bond_dims, _, max_trunc, _ = load_convergence_h5(conv_h5)
    x = np.array(sweep_nums)
    colors = cm.viridis(np.linspace(0, 1, len(bond_dims)))

    fig, ax = plt.subplots(figsize=(7, 4))
    for bi, (bdim, color) in enumerate(zip(bond_dims, colors)):
        y = max_trunc[:, bi]
        with np.errstate(invalid='ignore'):
            valid = np.isfinite(y) & (y > 0)
        ax.semilogy(x[valid], y[valid], marker='o', markersize=5, linewidth=1.4,
                    color=color, label=f'm = {bdim}')

    ax.set_xlabel('Sweep')
    ax.set_ylabel('Max truncated weight')
    ax.set_title('Max truncated weight convergence vs bond dimension')
    ax.set_xticks(x)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    fig.tight_layout()

    out = os.path.join(savepath, plot_name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {out}")

def plot_sum_trunc_convergence(conv_h5, savepath, plot_name='sum_trunc_convergence.png'):
    """x = sweep number, y = sum of truncated weights (log scale), one line per bond dimension."""
    sweep_nums, bond_dims, _, _, sum_trunc = load_convergence_h5(conv_h5)
    x = np.array(sweep_nums)
    colors = cm.viridis(np.linspace(0, 1, len(bond_dims)))

    fig, ax = plt.subplots(figsize=(7, 4))
    for bi, (bdim, color) in enumerate(zip(bond_dims, colors)):
        y = sum_trunc[:, bi]
        with np.errstate(invalid='ignore'):
            valid = np.isfinite(y) & (y > 0)
        ax.semilogy(x[valid], y[valid], marker='o', markersize=5, linewidth=1.4,
                    color=color, label=f'm = {bdim}')

    ax.set_xlabel('Sweep')
    ax.set_ylabel('Sum of truncated weights')
    ax.set_title('Sum of truncated weights convergence vs bond dimension')
    ax.set_xticks(x)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    fig.tight_layout()

    out = os.path.join(savepath, plot_name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {out}")

def plot_sum_trunc_vs_energy(conv_h5, savepath, plot_name='sum_trunc_vs_energy.png'):
    """x = sum of truncated weights (log scale), y = |E - E_min| (log scale), one point per sweep colored by bond dimension."""
    sweep_nums, bond_dims, energies, _, sum_trunc = load_convergence_h5(conv_h5)
    colors = cm.viridis(np.linspace(0, 1, len(bond_dims)))
    e_min = np.nanmin(energies)

    fig, ax = plt.subplots(figsize=(7, 4))
    for bi, (bdim, color) in enumerate(zip(bond_dims, colors)):
        x = sum_trunc[:, bi]
        y = np.abs(energies[:, bi] - e_min)
        with np.errstate(invalid='ignore'):
            valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        ax.loglog(x[valid], y[valid], marker='o', markersize=5,
                  linestyle='None', color=color, label=f'm = {bdim}')
        valid_idx = np.where(valid)[0]
        if len(valid_idx) > 0:
            last = valid_idx[-1]
            ax.loglog(x[last], y[last], marker='o', markersize=11,
                      linestyle='None', markerfacecolor='none',
                      markeredgecolor='red', markeredgewidth=1.5)

    ax.set_xlabel('Sum of truncated weights')
    ax.set_ylabel('|E - E$_{min}$|')
    ax.set_title('|E - E$_{min}$| vs sum of truncated weights colored by bond dimension')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    fig.tight_layout()

    out = os.path.join(savepath, plot_name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {out}")

def plot_max_trunc_vs_energy(conv_h5, savepath, plot_name='max_trunc_vs_energy.png'):
    """x = max truncated weight (log scale), y = |E - E_min| (log scale), one point per sweep colored by bond dimension."""
    sweep_nums, bond_dims, energies, max_trunc, _ = load_convergence_h5(conv_h5)
    colors = cm.viridis(np.linspace(0, 1, len(bond_dims)))
    e_min = np.nanmin(energies)

    fig, ax = plt.subplots(figsize=(7, 4))
    for bi, (bdim, color) in enumerate(zip(bond_dims, colors)):
        x = max_trunc[:, bi]
        y = np.abs(energies[:, bi] - e_min)
        with np.errstate(invalid='ignore'):
            valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        ax.loglog(x[valid], y[valid], marker='o', markersize=5,
                  linestyle='None', color=color, label=f'm = {bdim}')
        valid_idx = np.where(valid)[0]
        if len(valid_idx) > 0:
            last = valid_idx[-1]
            ax.loglog(x[last], y[last], marker='o', markersize=11,
                      linestyle='None', markerfacecolor='none',
                      markeredgecolor='red', markeredgewidth=1.5)

    ax.set_xlabel('Max truncated weight')
    ax.set_ylabel('|E - E$_{min}$|')
    ax.set_title('|E - E$_{min}$| vs max truncated weight colored by bond dimension')
    ax.legend(fontsize=9)
    ax.grid(True, which='both', linestyle='--', alpha=0.4)
    fig.tight_layout()

    out = os.path.join(savepath, plot_name)
    fig.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Saved figure to {out}")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collect truncation-weight HDF5 files across bond dimensions '
                    'into a single convergence HDF5.')
    parser.add_argument('--h5_files', nargs='+', required=True,
                        help='List of HDF5 files, one per bond dimension.')
    parser.add_argument('--bond_dims', nargs='+', type=int, default=None,
                        help='Bond dimensions, in the same order as --h5_files. '
                             'If omitted, extracted from filenames via _m<N>_ pattern.')
    parser.add_argument('--output', '-o', required=True,
                        help='Path for the output HDF5 file.')
    parser.add_argument('--savepath', '-s', default=None,
                        help='Directory for the plots. Defaults to the directory of --output.')
    parser.add_argument('--plot_root', default='convergence',
                        help='Prefix for output plot filenames (default: convergence).')
    args = parser.parse_args()

    if args.bond_dims is not None:
        if len(args.bond_dims) != len(args.h5_files):
            parser.error('--bond_dims must have the same length as --h5_files')
        bond_dims = args.bond_dims
    else:
        bond_dims = [extract_bond_dim(p) for p in args.h5_files]

    collect(args.h5_files, bond_dims, args.output)

    savepath = args.savepath if args.savepath else os.path.dirname(os.path.abspath(args.output))
    os.makedirs(savepath, exist_ok=True)

    plot_energy_convergence(
        args.output, savepath,
        plot_name=f'{args.plot_root}_energy.png')
    # plot_max_trunc_convergence(
    #     args.output, savepath,
    #     plot_name=f'{args.plot_root}_max_trunc.png')
    plot_sum_trunc_convergence(
        args.output, savepath,
        plot_name=f'{args.plot_root}_sum_trunc.png')
    plot_sum_trunc_vs_energy(
        args.output, savepath,
        plot_name=f'{args.plot_root}_sum_trunc_vs_energy.png')
    # plot_max_trunc_vs_energy(
    #     args.output, savepath,
    #     plot_name=f'{args.plot_root}_max_trunc_vs_energy.png')

