#!/bin/bash

run_dmrg() {

local DMRG_EXE=/home/giacomo/code/maquis-dmrg-core_valence_separated/build/qcmaquis
local FOLDER=$1
local KIND=$2
local PRJ_NAME=$3
local BOND_DIM=$4
local INITIAL_GUESS=$5
local FCIDUMP_NAME=$6
local DO_DMRG=$7
local CHKP_NAME=$8


PRJ_NAME=${PRJ_NAME}_m${BOND_DIM}_${INITIAL_GUESS}

CALC_DIR=${FOLDER}/${PRJ_NAME}
mkdir -p $CALC_DIR
cd $CALC_DIR || { echo "Error changing directory to ${CALC_DIR}, check if the path is correct."; exit 1; }

LOG_FILE=${CALC_DIR}/${PRJ_NAME}.QCMaquis.log
echo "Running DMRG for project ${PRJ_NAME} with bond dimension ${BOND_DIM} and initial guess ${INITIAL_GUESS}..."
echo "Parameters: FOLDER=${FOLDER}, FCIDUMP_NAME=${FCIDUMP_NAME}, PRJ_NAME=${PRJ_NAME}, BOND_DIM=${BOND_DIM}, INITIAL_GUESS=${INITIAL_GUESS}, CHKP_NAME=${CHKP_NAME}, DO_DMRG=${DO_DMRG}"
echo "Am I running DMRG? $DO_DMRG"
if [ "$DO_DMRG" = true ]; then
        cp ${FOLDER}/${FCIDUMP_NAME}.FCIDUMP ${CALC_DIR}/current.FCIDUMP || { echo "Error copying FCIDUMP file, check if the file exists and the path is correct."; exit 1; }
        cp ${FOLDER}/${KIND}_template.inp ${CALC_DIR}/${PRJ_NAME}.inp || { echo "Error copying template input file, check if the file exists and the path is correct."; exit 1; }
        sed -i "s|max_bond_dimension = XXX|max_bond_dimension = ${BOND_DIM}|" ${CALC_DIR}/${PRJ_NAME}.inp
        if [ "$INITIAL_GUESS" = "chkp" ]; then
                if [ -z "$CHKP_NAME" ]; then
                        echo "Error: CHKP_NAME variable is not set for checkpoint initial guess."
                        exit 1
                fi
                rm -rf ${CALC_DIR}/current.checkpoint
                cp -r ${FOLDER}/${CHKP_NAME} ${CALC_DIR}/current.checkpoint || { echo "Error copying checkpoint file, check if the file exists and the path is correct."; exit 1; }
        elif [ "$INITIAL_GUESS" = "hf" ]; then
                sed -i "s|// init_type = 'hf'|init_type = 'hf'|" ${CALC_DIR}/${PRJ_NAME}.inp
                sed -i "s|// hf_occ = |hf_occ = |" ${CALC_DIR}/${PRJ_NAME}.inp
                sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_hf.checkpoint'|" ${CALC_DIR}/${PRJ_NAME}.inp
                rm -rf ${CALC_DIR}/current_hf.checkpoint
        elif [ "$INITIAL_GUESS" = "random" ]; then
                sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_random.checkpoint'|" ${CALC_DIR}/${PRJ_NAME}.inp   
                rm -rf ${CALC_DIR}/current_random.checkpoint
        elif [ "$INITIAL_GUESS" = "coreX" ]; then
                sed -i "s|// init_type = 'hf'|init_type = 'coreX'|" ${CALC_DIR}/${PRJ_NAME}.inp
                sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_coreX.checkpoint'|" ${CALC_DIR}/${PRJ_NAME}.inp
        fi
        $DMRG_EXE ${PRJ_NAME}.inp &> $LOG_FILE || { echo "Error in DMRG run, check the log file for details"; exit 1; }
fi

python /home/giacomo/code/pyMPSDiagnosis/parse_truncation_weight.py \
        --log_file $LOG_FILE > ${CALC_DIR}/${PRJ_NAME}_truncation_weights.parse.log \
         || { echo "Error in parsing truncation weights, check the log file for details"; exit 1; }

mkdir -p ${FOLDER}/AA_plots

python /home/giacomo/code/pyMPSDiagnosis/plot_truncation_weight.py \
    --h5_file ${CALC_DIR}/${PRJ_NAME}.QCMaquis_truncated_weights.h5 \
    --savepath ${FOLDER}/AA_plots \
    --output_root ${PRJ_NAME} || { echo "Error in plotting truncation weights, check the log file for details"; exit 1; }
}




folder=/home/giacomo/dmrg_calcs/DMRG_DIAGNOSIS/FeH12O6
fcidump_name=FCIDUMP_2s2p3p_groundORBS

do_dmrg=true

kind=FeH12O6_2pcvsXXS_reverse
prj_name=${kind}


# initial_guess can be 'hf' or 'random' or 'chkp'
# chkp_name=chkp/FeH12O6.noCVS.checkpoint_state.0.h5


bond_dim_list=(64)

initial_guess=chkp
chkp_name=FeH12O6_2pcvsXXS_m1024_coreX.checkpoint


#####################################################


for bond_dim in "${bond_dim_list[@]}"; do
    run_dmrg $folder $kind $prj_name $bond_dim \
            $initial_guess $fcidump_name $do_dmrg $chkp_name || { echo "Error in DMRG run for bond dimension ${bond_dim}, check the log file for details"; exit 1; }
done
h5_files=""
for bond_dim in "${bond_dim_list[@]}"; do
    h5_files="${h5_files} ${folder}/${prj_name}_m${bond_dim}_${initial_guess}/${prj_name}_m${bond_dim}_${initial_guess}.QCMaquis_truncated_weights.h5"
done
python /home/giacomo/code/pyMPSDiagnosis/collect_bond_dim_convergence.py \
    --h5_files $h5_files \
    --output ${folder}/${prj_name}_${initial_guess}_convergence.h5 \
    --savepath  ${folder}/AA_plots \
    --plot_root ${prj_name}_${initial_guess}

