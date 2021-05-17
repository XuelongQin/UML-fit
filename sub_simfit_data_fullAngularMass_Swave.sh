#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
Executable  = run_simfit_data_fullAngularMass_Swave.sh
nsamp       = 0
bin         = \$(ProcId) / 4
yearConf    = \$(ProcId) % 4
Arguments   = \$INT(nsamp) \$INT(bin) \$INT(yearConf)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf)_noMeanCon.out
Error       = logs_parSub/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf)_noMeanCon.err
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
echo 'Queue 8'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# Compilation, submission and file removal
if make simfit_data_fullAngularMass_Swave
then condor_submit temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
rm temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
