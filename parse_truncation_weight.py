import re
import sys
import numpy as np
import pandas as pd
import h5py
import argparse


def parse_dmrg_log(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Split into sweeps
    sweep_blocks = re.split(r'SWEEP (\d+);', content)
    # sweep_blocks: [preamble, sweep_num, block, sweep_num, block, ...]

    sweeps = {}
    for i in range(1, len(sweep_blocks), 2):
        sweep_num = int(sweep_blocks[i])
        block = sweep_blocks[i + 1]

        # Match site labels like -0-1-> or <-11-12- followed (anywhere after)
        # by "Truncated weight = <number> [eV]"
        # Works for both old format (concatenated lines) and new format (clean lines + " eV")
        pattern = (
            r'(-\d+-\d+->|<-\d+-\d+-)'   # site label (forward or backward)
            r'.*?'                         # anything in between (non-greedy)
            r'Truncated weight\s*=\s*'     # key
            r'([0-9eE+\-\.]+)'             # numeric value
            r'(?:\s*eV)?'                  # optional unit
        )
        matches = re.findall(pattern, block, re.DOTALL)

        site_weights = {}
        for site_label, weight in matches:
            site_label = site_label.strip()
            site_weights[site_label] = float(weight)

        # Extract sweep-level energy
        energy_match = re.search(
            r'SWEEP ENERGY\s*=\s*([0-9eE+\-\.]+)', block)
        sweep_energy = float(energy_match.group(1)) if energy_match else float('nan')

        sweeps[sweep_num] = (site_weights, sweep_energy)

    # Build DataFrame: rows = site labels, columns = sweep numbers
    all_sites = []
    for site_weights, _ in sweeps.values():
        for site in site_weights:
            if site not in all_sites:
                all_sites.append(site)

    sorted_sweeps = sorted(sweeps.keys())
    df = pd.DataFrame(index=all_sites, columns=sorted_sweeps, dtype=float)
    df.index.name = 'site'

    sweep_energies = {}
    for sweep_num, (site_weights, sweep_energy) in sweeps.items():
        for site, weight in site_weights.items():
            df.loc[site, sweep_num] = weight
        sweep_energies[sweep_num] = sweep_energy

    df.columns.name = 'sweep'
    energies = np.array([sweep_energies[s] for s in sorted_sweeps])
    return df, energies


def save_to_hdf5(df, energies, filepath):
    with h5py.File(filepath, 'w') as f:
        f.create_dataset('truncated_weights',
                         data=df.values.astype(float))
        f.create_dataset('site_labels',
                         data=np.array(df.index.tolist(), dtype='S'))
        f.create_dataset('sweep_numbers',
                         data=np.array(df.columns.tolist(), dtype=int))
        f.create_dataset('sweep_energies',
                         data=energies.astype(float))
        f.attrs['description'] = 'DMRG truncated weights and sweep energies: rows=sites, cols=sweeps'

    print(f"Saved to {filepath}")
    print(f"  truncated_weights : {df.values.shape}")
    print(f"  site_labels       : {df.index.tolist()}")
    print(f"  sweep_numbers     : {df.columns.tolist()}")
    print(f"  sweep_energies    : {energies.tolist()}")


def load_from_hdf5(filepath):
    """Reload the DataFrame and sweep energies from the HDF5 file."""
    with h5py.File(filepath, 'r') as f:
        data        = f['truncated_weights'][:]
        site_labels = [s.decode() for s in f['site_labels'][:]]
        sweep_nums  = f['sweep_numbers'][:].tolist()
        energies    = f['sweep_energies'][:] if 'sweep_energies' in f else None
    df = pd.DataFrame(data, index=site_labels, columns=sweep_nums)
    df.index.name   = 'site'
    df.columns.name = 'sweep'
    return df, energies


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse DMRG log file for truncated weights.')
    parser.add_argument('--log_file', help='Path to the DMRG log file', required=True)
    args = parser.parse_args()

    filepath = args.log_file
    df, energies = parse_dmrg_log(filepath)

    print(df.to_string())
    print(f"\nShape: {df.shape}  (sites x sweeps)")
    print(f"\nSweep energies: {energies.tolist()}")

    out = filepath.replace('.log', '_truncated_weights.h5')
    save_to_hdf5(df, energies, out)