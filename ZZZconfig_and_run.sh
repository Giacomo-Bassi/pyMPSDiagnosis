#!/bin/bash
set -euo pipefail


source /dss/dsshome1/03/go67sec2/code/pyMPSDiagnosis/run_qcmaquis_dev.sh
# ---------------------------------------------------------------------------
# Main — edit only this section for different runs
# ---------------------------------------------------------------------------

declare -A CFG=(
    [dmrg_exe]=/dss/dsshome1/03/go67sec2/code/src_qcmaquis/maquis-dmrg-core_valence_separated/build/qcmaquis
    [pyscript_folder]=/dss/dsshome1/03/go67sec2/code/pyMPSDiagnosis
    [folder]=/dss/dsshome1/03/go67sec2/Projects/DMRG_calcs/NiC10H10/qcmaquis_CIonly_1616
    [fcidump_name]=NiC10H10_1616
    [kind]=NiC10H10_dmrgscf_1616
    [initial_guess]=random
    [do_dmrg]=false
    [bond_dim_list]="64 128 256 512"
    # [chkp_name]=chkp/FeH12O6.noCVS.checkpoint_state.0.h5   # only needed for initial_guess=chkp
)

run_wrapper_qcmaquis CFG