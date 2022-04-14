#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
Executable  = run_simfit_data_fullAngularMass_Swave.sh
nsamp       = 0
bin         = 4
yearConf    = 3
q2stat      = 5
save	    = 2
Arguments   = \$INT(nsamp) \$INT(bin) \$INT(yearConf) \$INT(q2stat) \$INT(save)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf).out
Error       = logs_parSub/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
echo 'Queue 1'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# Compilation, submission and file removal
if make simfit_data_fullAngularMass_Swave
then condor_submit temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
rm temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# bin         = ( \$(ProcId) / 4 ) * 2 + 4
# yearConf    = \$(ProcId) % 4
