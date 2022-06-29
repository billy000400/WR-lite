# @Author: Billy Li <billyli>
# @Date:   06-24-2022
# @Email:  li000400@umn.edu
# @Last modified by:   billyli
# @Last modified time: 06-29-2022


# this script is for calcualting the effective fsig at (1600, 800)
# from the published Xsec expectation (see arxiv 2112.03949)
xSec_mumu = 0.8 #/fb
xSec_ee = 1.2 #/fb
lumi = 137 #/fb
selEff_mumu = 0.5127393291316945 # select effiency
selEff_ee = 0.40538845246273786

N_sig_mumu = xSec_mumu*lumi*selEff_mumu
N_sig_ee = xSec_ee*lumi*selEff_ee

N_bg_mumu = 30098
N_bg_ee = 21416

fsig_mumu = N_sig_mumu/N_bg_mumu
fsig_ee = N_sig_ee/N_bg_ee

print(fsig_mumu)
print(fsig_ee)
