#!/bin/bash

DMRG_EXE=/home/giacomo/code/maquis-dmrg-core_valence_separated/build/qcmaquis

folder=/home/giacomo/dmrg_calcs/DMRG_DIAGNOSIS/FeH12O6
fcidump_name=FCIDUMP_2s2p3p_groundORBS

molecule_name=FeH12O6

prj_name=${molecule_name}_16e22o.000_GS
bond_dim=1024

# initial_guess can be 'hf' or 'random' or 'chkp'
initial_guess=hf
# chkp_name=chkp/FeH12O6.noCVS.checkpoint_state.0.h5

prj_name=${prj_name}_m${bond_dim}_${initial_guess}

calc_dir=${folder}/${prj_name}
save_folder=${folder}/AA_plots

mkdir -p $calc_dir
mkdir -p $save_folder

cd $calc_dir

python /home/giacomo/code/pyMPSDiagnosis/plot_truncation_weight.py \
    --h5_file ${calc_dir}/${prj_name}.QCMaquis_truncated_weights.h5 \
    --savepath ${save_folder} \
    --output ${prj_name}_truncated_weights_plot.png