#!/bin/bash

# Create directory for log files
if [ ! -d logs_parSub ]; then mkdir logs_parSub; fi

# Creation of the submit HTCondor file
cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub
Executable  = run_simfit_recoMC_fullAngularMass_toybkg_toyeff.sh
neff        = 3
bin         = 4
Arguments   = \$INT(neff) \$INT(bin)
Log         = logs_parSub/sub_\$(ClusterId).log
Output      = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_toyeff_\$INT(neff)_\$INT(bin).out
Error       = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_toyeff_\$INT(neff)_\$INT(bin).err
transfer_output_files = ""
+JobFlavour = "testmatch"
EOF

if [ "${USER}" == "fiorendi" ]; then
    echo '+AccountingGroup = "group_u_CMST3.all"'>>temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub
fi
echo 'Queue 1'>>temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub
# echo 'Queue 60'>>temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub

# Compilation, submission and file removal
if make simfit_recoMC_fullAngularMass_toybkg_toyeff
then condor_submit temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub
fi
rm temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub

# # Creation of the submit HTCondor file
# cat << EOF > temp_sub_simfit_recoMC_fullAngularMass_toybkg_toyeff_oneBin.sub
# Executable  = run_simfit_recoMC_fullAngularMass_toybkg_toyeff.sh
# neff        = ( \$(ProcId) / 6 )
# bin         = \$(ProcId) % 6
# Arguments   = \$INT(neff) \$INT(bin)
# Log         = logs_parSub/sub_\$(ClusterId).log
# Output      = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_toyeff_\$INT(neff)_\$INT(bin).out
# Error       = logs_parSub/simfit_recoMC_fullAngularMass_toybkg_toyeff_\$INT(neff)_\$INT(bin).err
# transfer_output_files = ""
# +JobFlavour = "testmatch"
# EOF
