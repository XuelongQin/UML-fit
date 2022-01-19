#!/bin/bash

# Create directory for log files
logdir="logs_parSub_freeSigmas"
if [ ! -d ${logdir} ]; then mkdir ${logdir}; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
Executable  = run_simfit_data_fullAngularMass_Swave.sh
nsamp       = 0
bin         = \$(ProcId) * 2 + 4
yearConf    = 0
q2stat      = 5
save	    = 2
Arguments   = \$INT(nsamp) \$INT(bin) \$INT(yearConf) \$INT(q2stat) \$INT(save)
Log         = ${logdir}/sub_\$(ClusterId).log
Output      = ${logdir}/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf).out
Error       = ${logdir}/simfit_data_fullAngularMass_Swave_\$INT(nsamp)_\$INT(bin)_\$INT(yearConf).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
echo 'Queue 2'>>temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub

# Compilation, submission and file removal
if make simfit_data_fullAngularMass_Swave
then condor_submit temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
fi
rm temp_sub_simfit_data_fullAngularMass_Swave_oneBin.sub
