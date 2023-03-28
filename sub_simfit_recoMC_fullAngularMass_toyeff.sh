#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_toyeff_oneBin.sub
Executable  = run_simfit_recoMC_fullAngularMass_toyeff.sh
bin         = \$(ProcId)
Arguments   = \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_recoMC_fullAngularMass_toyeff_\$INT(bin).out
Error       = logs_parSub/simfit_recoMC_fullAngularMass_toyeff_\$INT(bin).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngularMass_toyeff_oneBin.sub
fi
echo 'Queue 6'>>temp_sub_simfit_recoMC_fullAngularMass_toyeff_oneBin.sub

# Compilation, submission and file removal
if make simfit_recoMC_fullAngularMass_toyeff
then condor_submit temp_sub_simfit_recoMC_fullAngularMass_toyeff_oneBin.sub
fi
rm temp_sub_simfit_recoMC_fullAngularMass_toyeff_oneBin.sub
