# @Author: Billy Li <billyli>
# @Date:   06-24-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 06-24-2022


# this script is for calcualting the effective fsig at (1600, 800)
# from the 95 CL Xsec in the old paper
xSec = 0.8 #/fb
lumi = 137 #/fb
selEff_mumu = 0.5127393291316945 # select effiency
selEff_ee = 0.40538845246273786

N_sig_mumu = xSec*lumi*selEff_mumu
N_sig_ee = xSec*lumi*selEff_ee

N_bg_mumu = 30098
N_bg_ee = 21416

fsig_mumu = N_sig_mumu/N_bg_mumu
fsig_ee = N_sig_ee/N_bg_ee

print(fsig_mumu)
print(fsig_ee)
