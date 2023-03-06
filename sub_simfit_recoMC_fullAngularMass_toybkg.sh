#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub/xgbv8 ]; then mkdir -p logs_parSub/xgbv8; fi

nbins=$(wc -l ../confSF/KDE_SF_all.list | cut -d " " -f1)

while read -a line; do
    bin=${line[0]}
    # Creation of the submit HTCondor file
    cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
Executable  = run_simfit_recoMC_fullAngularMass_toybkg.sh
nsamp       = \$(ProcId)
Arguments   = \$INT(nsamp) ${bin} 
Log         = logs_parSub/xgbv8/sub_\$(ClusterId).log
Output      = logs_parSub/xgbv8/simfit_recoMC_fullAngularMass_toybkg_\$INT(nsamp)_${bin}.out
Error       = logs_parSub/xgbv8/simfit_recoMC_fullAngularMass_toybkg_\$INT(nsamp)_${bin}.err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

    if [ "${USER}" == "fiorendi" ]; then
        echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin${bin}.sub
    fi
    echo "Queue $((100))">>temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin${bin}.sub


    # Compilation, submission and file removal
    if make simfit_recoMC_fullAngularMass_toybkg
    then condor_submit temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
    fi
    rm temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
done < ../confSF/KDE_SF_all.list
    
