#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

nbins=$(wc -l ../confSF/KDE_SF.list | cut -d " " -f1)

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
Executable  = run_simfit_recoMC_fullAngularMass_toybkg.sh
nsamp       = ( \$(ProcId) / ${nbins} ) + 1
bin         = \$(ProcId) % ${nbins}
Arguments   = \$INT(nsamp) \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_\$INT(nsamp)_\$INT(bin).out
Error       = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_\$INT(nsamp)_\$INT(bin).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
fi
echo "Queue $((100*${nbins}))">>temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub

# Compilation, submission and file removal
if make simfit_recoMC_fullAngularMass_toybkg
then condor_submit temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
fi
rm temp_sub_simfit_recoMC_fullAngularMass_toybkg_oneBin.sub
