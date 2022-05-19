from pathlib import Path
import numpy as np
import pandas as pd

dist = "ExpmCB"
recoMethods = ['ee', 'mumu']
for recoMethod in recoMethods:
    # dict for df
    data = {"WR":[], "N":[],\
            "alpha":[], "alpha_err":[],\
            "n":[], "n_err":[],\
            "beta":[], "beta_err":[],\
            "m":[], "m_err":[],\
            "mean":[], "mean_err":[],\
            "sigma":[], "sigma_err":[]
            }

    result_dir = Path.cwd().joinpath("results_"+dist+"_"+recoMethod)
    files = result_dir.glob('*.txt')
    for file in files:
        f = open(file, mode='r')
        lines = f.readlines()
        f.close()

        # extract WR mass
        if "WR:" in lines[0]:
            WRMass = int(lines[0][3:])
        elif "WR" in lines[0]:
            WRMass = int(lines[0][2:])
        data['WR'].append(WRMass)

        # extract N mass
        if "N:" in line[1]:
            NMass = int(lines[1][2:])
        elif "N" in line[1]:
            NMass = int(lines[1][1:])
        data['N'].append(NMass)

        # extract alpha
        line_alpha = lines[9]
        name_alpha, alpha_init, alpha_result = line_alpha.split("    ")
        alpha_final, alpha_errCorr = alpha_result.split(" +/- ")
        alpha_error, _ = alpha_errCorr.split("  ")
        data['alpha'].append(float(alpha_final))
        data['alpha_err'].append(float(alpha_error))

        # extract n
        line_n = lines[13]
        name_n, n_init, n_result = line_n.split("    ")
        n_final, n_errCorr = n_result.split(" +/- ")
        n_error, _ = n_errCorr.split("  ")
        data['n'].append(float(n_final))
        data['n_err'].append(float(n_error))

        # extract beta
        line_beta = lines[10]
        name_beta, beta_init, beta_result = line_beta.split("    ")
        beta_final, beta_errCorr = beta_result.split(" +/- ")
        beta_error, _ = beta_errCorr.split("  ")
        data['beta'].append(float(beta_final))
        data['beta_err'].append(float(beta_error))

        # extract m
        line_m = lines[11]
        name_m, m_init, m_result = line_m.split("    ")
        m_final, m_errCorr = m_result.split(" +/- ")
        m_error, _ = m_errCorr.split("  ")
        data['m'].append(float(m_final))
        data['m_err'].append(float(m_error))

        # extract mean
        line_mean = lines[12]
        name_mean, mean_init, mean_result = line_mean.split("    ")
        mean_final, mean_errCorr = mean_result.split(" +/- ")
        mean_error, _ = mean_errCorr.split("  ")
        data['mean'].append(float(mean_final))
        data['mean_err'].append(float(mean_error))

        # extract sigma
        line_sigma = lines[14]
        name_sigma, sigma_init, sigma_result = line_sigma.split("    ")
        sigma_final, sigma_errCorr = sigma_result.split(" +/- ")
        sigma_error, _ = sigma_errCorr.split("  ")
        data['sigma'].append(float(sigma_final))
        data['sigma_err'].append(float(sigma_error))

    df = pd.DataFrame.from_dict(data)
    print(df)
