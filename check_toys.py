import pdb
import pandas as pd
#from tabulate import tabulate
import os
import glob
#import numpy as np
#from pathlib import Path

# bins = [2,5]
# ntoys = 52

bins = [0,1,2,3,5,7]
# bins = [3]
ntoys = 100

for ibin in bins:
    print ('############ Processing bin', ibin, '################' )
    count_bad = 0
    count_converged = 0
    count_missing_all = 0
    count_missing_root = 0
    count_missing_root_converged = 0
    count_error = 0
    count_else = 0
    for itoy in range(ntoys):
        
        path_root = 'simFitResults4d/xgbv8/simFitResult_recoMC_fullAngularMass201620172018_dataStat-%s_b%s_XGBv8.root'%(itoy,ibin)
        #logs_parSub/xgbv8/
        path_out = 'logs_parSub/xgbv8/simfit_recoMC_fullAngularMass_*_*_%s_%s.*'%(itoy,ibin)
        out_files = glob.glob(path_out)
        root_files = glob.glob(path_root)
        
#         print (len(root_files), len(out_files))

        if len(root_files) == 0 and len(out_files)>0:
            print ('missing root file for toy %s bin %s'%(itoy,ibin))
            count_missing_root += 1
            out_file = [f for f in out_files if 'out' in f]
            with open(out_file[0]) as fout:
                out_lines = fout.readlines()
                for iline in out_lines:
                    if 'Converged' in iline:
                        count_missing_root_converged += 1
            
        elif len(root_files) == 0 and len(out_files) == 0:
#             print ('missing everything for toy %s bin %s'%(itoy,ibin))
            count_missing_all +=1

        elif len(root_files) == 1 and len(out_files) == 2:
#             print ('have everything for toy %s bin %s'%(itoy,ibin))
            out_file = [f for f in out_files if 'out' in f]
            found = False
            with open(out_file[0]) as fout:
                out_lines = fout.readlines()
                for iline in out_lines:
                    if 'Not converged' in iline:
#                         print ('not --->', iline)
                        count_bad = count_bad + 1
                        found = True
                        break
                    elif 'Converged' in iline:
#                         print ('yes --->', iline)
                        count_converged = count_converged + 1
                        found = True
                        break
                if not found:
                    print ('Error in toy %s bin %s:'%(itoy,ibin)) 
                    print (out_lines[-1])      
                    count_error +=1  

        else:
            print ('have %s root file and %s log files for toy %s bin %s'%(len(root_files), len(out_files), itoy,ibin))
            count_else +=1
                            
    print ('######## Summary for bin', ibin, '########' )
    print ('Converged:',  count_converged )
    print ('Failed:',  count_bad )
    print ('Missing root file:',  count_missing_root )
    print ('\t of which converged:',  count_missing_root_converged )
    print ('Missing all:',  count_missing_all )
    print ('Error/incomplete:',  count_error )
    print ('Not classified:',  count_else )
    print ('Total:',  count_converged+count_bad+count_missing_root+count_missing_all+count_error+count_else )
    
    
#     print ('######################################' )

