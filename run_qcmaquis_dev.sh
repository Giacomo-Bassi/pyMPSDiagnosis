#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# QCMaquis DMRG runner
# ---------------------------------------------------------------------------
# Config is passed between functions as an associative array (nameref).
# Call pattern:  func_name cfg  (where cfg is declared with declare -A)
# ---------------------------------------------------------------------------

setup_initial_guess() {          # args: array-name  prj_name  calc_dir
    local cfg_name=$1
    local -n _cfg=$cfg_name
    local prj_name=$2
    local calc_dir=$3

    case "${_cfg[initial_guess]}" in
        chkp)
            [[ -z "${_cfg[chkp_name]:-}" ]] && { echo "Error: chkp_name must be set for checkpoint initial guess."; exit 1; }
            rm -rf "${calc_dir}/current.checkpoint"
            cp -r  "${_cfg[folder]}/${_cfg[chkp_name]}" "${calc_dir}/current.checkpoint"
            ;;
        hf)
            sed -i "s|// init_type = 'hf'|init_type = 'hf'|"   "${calc_dir}/${prj_name}.inp"
            sed -i "s|// hf_occ = |hf_occ = |"                 "${calc_dir}/${prj_name}.inp"
            sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_hf.checkpoint'|" \
                                                                 "${calc_dir}/${prj_name}.inp"
            rm -rf "${calc_dir}/current_hf.checkpoint"
            ;;
        random)
            sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_random.checkpoint'|" \
                                                                 "${calc_dir}/${prj_name}.inp"
            rm -rf "${calc_dir}/current_random.checkpoint"
            ;;
        coreX)
            sed -i "s|// init_type = 'hf'|init_type = 'coreX'|" "${calc_dir}/${prj_name}.inp"
            sed -i "s|chkpfile = 'current.checkpoint'|chkpfile = 'current_coreX.checkpoint'|" \
                                                                  "${calc_dir}/${prj_name}.inp"
            ;;
        *)
            echo "Error: unknown initial_guess '${_cfg[initial_guess]}'"; exit 1 ;;
    esac
}

run_single_dmrg() {              # args: array-name  bond_dim
    local cfg_name=$1
    local -n _cfg=$cfg_name
    local bond_dim=$2

    local prj="${_cfg[kind]}_m${bond_dim}_${_cfg[initial_guess]}"
    local calc_dir="${_cfg[folder]}/${prj}"
    local log_file="${calc_dir}/${prj}.QCMaquis.log"

    mkdir -p "$calc_dir"
    cd       "$calc_dir"

    echo "Running DMRG: project=${prj}  m=${bond_dim}  guess=${_cfg[initial_guess]}"
    echo "  folder=${_cfg[folder]}  fcidump=${_cfg[fcidump_name]}  do_dmrg=${_cfg[do_dmrg]}"

    if [[ "${_cfg[do_dmrg]}" == true ]]; then
        cp "${_cfg[folder]}/${_cfg[fcidump_name]}.FCIDUMP" "${calc_dir}/current.FCIDUMP"
        cp "${_cfg[folder]}/${_cfg[kind]}_template.inp"    "${calc_dir}/${prj}.inp"
        sed -i "s|max_bond_dimension = XXX|max_bond_dimension = ${bond_dim}|" "${calc_dir}/${prj}.inp"

        setup_initial_guess "$cfg_name" "$prj" "$calc_dir"

        "${_cfg[dmrg_exe]}" "${prj}.inp" &> "$log_file" \
            || { echo "Error in DMRG run — see ${log_file}"; exit 1; }
    fi

    mkdir -p "${_cfg[folder]}/AA_plots"

    python "${_cfg[pyscript_folder]}/parse_truncation_weight.py" \
        --log_file "$log_file" \
        > "${calc_dir}/${prj}_truncation_weights.parse.log" \
        || { echo "Error parsing truncation weights"; exit 1; }

    python "${_cfg[pyscript_folder]}/plot_truncation_weight.py" \
        --h5_file     "${calc_dir}/${prj}.QCMaquis_truncated_weights.h5" \
        --savepath    "${_cfg[folder]}/AA_plots" \
        --output_root "$prj" \
        || { echo "Error plotting truncation weights"; exit 1; }
}

run_wrapper_qcmaquis() {         # args: array-name
    local cfg_name=$1
    local -n _cfg=$cfg_name

    for bond_dim in ${_cfg[bond_dim_list]}; do
        run_single_dmrg "$cfg_name" "$bond_dim" \
            || { echo "Error in DMRG run for m=${bond_dim}"; exit 1; }
    done

    local h5_files=""
    for bond_dim in ${_cfg[bond_dim_list]}; do
        h5_files+=" ${_cfg[folder]}/${_cfg[kind]}_m${bond_dim}_${_cfg[initial_guess]}/${_cfg[kind]}_m${bond_dim}_${_cfg[initial_guess]}.QCMaquis_truncated_weights.h5"
    done

    python "${_cfg[pyscript_folder]}/collect_bond_dim_convergence.py" \
        --h5_files  $h5_files \
        --output    "${_cfg[folder]}/${_cfg[kind]}_${_cfg[initial_guess]}_convergence.h5" \
        --savepath  "${_cfg[folder]}/AA_plots" \
        --plot_root "${_cfg[kind]}_${_cfg[initial_guess]}" \
        || { echo "Error collecting bond-dim convergence data"; exit 1; }
}

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