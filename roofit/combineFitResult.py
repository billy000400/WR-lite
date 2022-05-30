# @Author: Billy Li <billyli>
# @Date:   05-23-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 05-30-2022



from pathlib import Path
import numpy as np
import pandas as pd

dists = ["ExpmCB", "DSCB"]
recoMethods = ['ee', 'mumu']
for dist in dists:
    print(f"[INFO]: Combining fit results for {dist}")
    for recoMethod in recoMethods:
        data = {} # dict for making pd df
        catchVarName = False
        varNames = []

        # list all result files
        result_dir = Path.cwd().joinpath("results_"+dist+"_"+recoMethod)
        files = result_dir.glob('*.txt')
        for file in files:
            f = open(file, mode='r')
            lines = f.readlines()
            f.close()

            # catch all variable names at the first run
            if not catchVarName:
                data["WR"] = []
                data["N"] = []

                current_line_num=9 # RooFit record starts
                while (len(lines[current_line_num]) > 1):
                    current_line = lines[current_line_num]
                    varName = current_line.split(None,1)[0]
                    varNames.append(varName)
                    current_line_num += 1
                    # print(current_line_num, lines[current_line_num], len(lines[current_line_num]))

                for varName in varNames:
                    data[varName] = []
                    data[f"{varName}_err"] = []

                data["chi2"] = []
                catchVarName = True

            # extract WR mass
            WR_line = lines[0]
            WR_text = WR_line.split(":")[1]
            data['WR'].append(float(WR_text))

            # extract N mass
            N_line = lines[1]
            N_text = N_line.split(":")[1]
            data['N'].append(float(N_text))

            # loop through variable names to extract final values and errs
            current_line_num=9 # RooFit record starts
            for varName in varNames:
                current_line = lines[current_line_num]
                _,_,final,_,err,_=current_line.split()
                data[varName].append(final)
                data[f"{varName}_err"].append(err)

            # extract chi2
            line_chi2 = lines[-1]
            _, chi2 = line_chi2.split()
            data["chi2"].append(float(chi2))

        df = pd.DataFrame.from_dict(data)
        print(df)
        df.to_csv("csv/"+dist+"_"+recoMethod+"_fitResult.csv")
