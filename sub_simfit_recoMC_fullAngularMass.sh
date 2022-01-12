#!/bin/bash

sigstat=${1}
totsamp=${2}

# Create directory for log files
if [ ! -d logs_parSub_fracTest ]; then mkdir logs_parSub_fracTest; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_oneBin.sub
Executable  = run_simfit_recoMC_fullAngularMass.sh
nsamp       = ( \$(ProcId) / 6 ) + 1
bin         = \$(ProcId) % 6
Arguments   = \$INT(nsamp) \$INT(bin) ${sigstat}
Log         = logs_parSub_fracTest/sub_\$(ClusterId).log
Output      = logs_parSub_fracTest/simfit_recoMC_fullAngularMass_\$INT(nsamp)_\$INT(bin)_${sigstat}.out
Error       = logs_parSub_fracTest/simfit_recoMC_fullAngularMass_\$INT(nsamp)_\$INT(bin)_${sigstat}.err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngularMass_oneBin.sub
fi
echo "Queue $((6*${totsamp}))">>temp_sub_simfit_recoMC_fullAngularMass_oneBin.sub

# Compilation, submission and file removal
if make simfit_recoMC_fullAngularMass
then condor_submit temp_sub_simfit_recoMC_fullAngularMass_oneBin.sub
fi
rm temp_sub_simfit_recoMC_fullAngularMass_oneBin.sub
